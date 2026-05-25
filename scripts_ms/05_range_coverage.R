# Script to calculate the range coverage of a species before and after BTG

# load libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(terra)
library(gpkg)
library(BIEN)

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# load range coverage function
source("scripts_ms/functions/FUN_calc_range_coverage_change.R")

# list of taxa groups
taxagroups = list.files("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/")
inatgroups = c("Amphibia", "Aves", "Insecta", "Mammalia", "Plantae", "Reptilia")

# let's not do birds:
taxagroups = taxagroups[-2]
inatgroups = inatgroups[-2]

# Load the databases ===========================================================

canada_poly = terra::vect("data/base-layers/canada-polygon/canada.outline.shp")
canada_poly = terra::project(canada_poly, "EPSG:6933")

# load canada base grid
canada = terra::rast("data/base-layers/canada-basegrids/canada.base.5k.tiff")
canada = terra::aggregate(canada, factor = 2) # 10 km^2 cells
canada = terra::project(canada, "EPSG:6933")

# load parquet
inat.pre = arrow::read_parquet("data/heavy/BTG-data/inaturalist-canada-dec2025_PRE_JUN12025.parquet")
inat.post = arrow::read_parquet("data/heavy/BTG-data/inaturalist-canada-dec2025_JUN1OCT1.parquet")

# list of species to calculate range coverage for
df = lapply(as.list(paste0("outputs/range-coverage/range_coverage_", taxagroups, ".rds")), readRDS)
names(df) = taxagroups
df = bind_rows(df, .id = "Taxa")

summaries = vector("list", length = length(taxagroups))
names(summaries) = taxagroups
for(t in c(1:3,5)){ # we do plants separately below
  
  # Load species names that have previously been included
  species = dplyr::filter(df, Taxa == taxagroups[t]) |>
    select(species) |>
    as.matrix() |> as.vector()
  
  # Prepare to load IUCN ranges
  filepath = paste0("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/", taxagroups[t], "/")
  
  summaries[[t]] = list()
  for(i in 1:length(species)){
    tryCatch({  
      temp = calc_range_coverage_change(species_name = species[i], 
                                        range = terra::rast(paste0(filepath, species[i], ".tif")),
                                        inat.pre = inat.pre,
                                        inat.post = inat.post,
                                        canada_poly = canada_poly
      )
      
      terra::writeRaster(temp[[1]], 
                         paste0("~/McGill University/Laura's Lab_Group - range-coverage/",taxagroups[t],"/", species[i], "_JUN1OCT12025.tif"), overwrite=T)
      summaries[[t]][[i]] = temp[[2]]
    }, error = function(e) {
      message("An error occurred: ", e$message)
      return(NA) # Return NA if an error occurs
    })
  }
}

# bind results into a data frame
res_df = summaries |> lapply(bind_rows) |> bind_rows(.id = "taxa")
write.csv(res_df, "outputs/temporary_rangecoverage.csv")

## Plants version

for(t in c(4)){
  
  # Load species names that have previously been included
  species = dplyr::filter(df, Taxa == taxagroups[t]) |>
    select(species) |>
    as.matrix() |> as.vector()
  
  for(i in 1:length(species)){
    tryCatch({  
      
      range_sf <- BIEN::BIEN_ranges_load_species(species[i])
      
      # go to next species if there isn't a range polygon
      if(nrow(range_sf) == 0){ next
      } else {
        range_rast = range_sf |>
          terra::vect() |> 
          terra::project(crs(canada)) |>
          terra::rasterize(canada)
      }
      
      temp = calc_range_coverage_change(species_name = species[i], 
                                        range = range_rast,
                                        inat.pre = inat.pre,
                                        inat.post = inat.post, 
                                        canada_poly = canada_poly
      )
      
      terra::writeRaster(temp[[1]], 
                         paste0("~/McGill University/Laura's Lab_Group - range-coverage/",taxagroups[t],"/", species[i], "_JUN1OCT1.tif"), overwrite=T)
      summaries[[t]][[i]] = temp[[2]]
    }, error = function(e) {
      message("An error occurred: ", e$message)
      return(NA) # Return NA if an error occurs
    })
  }
}

# bind results into a data frame
res_df = summaries |> lapply(bind_rows) |> bind_rows(.id = "taxa")
res_df$prop_range_can = res_df$range_size_newcells_km/res_df$range_size_canada_km
res_df$prop_range_total = res_df$range_size_newcells_km/res_df$range_size_full_km
write.csv(res_df, "outputs/range-coverage/after/range_coverage_ALL_JUN1OCT12025.rds")

## Plot ------------------------------------------------------------------------

# cleaning
res_df = filter(res_df, sp != "Corispermum hookeri") # BIEN range seems to small and this makes the proportion gain huge
res_df = filter(res_df, sp != "Rosa glauca") # canadian range also too small. introduced sp 

# summaries for labels and sorting
res_summary = res_df |>
  filter(prop_range_can>0) |>
  group_by(taxa) |>
  summarise(
    "mu" = mean(prop_range_can),
    "max_sp" = sp[which.max(prop_range_can)]
  )
res_df$Taxa = factor(res_df$taxa, 
                     levels = res_summary$taxa[order(res_summary$mu)])

# Panel: jitter plot
(p.perc = ggplot(data = res_df) +
    geom_jitter(aes(x = 100*prop_range_can,
                    y = Taxa,
                    group = sp,
                    col = 100*prop_range_can),
                size = 3, alpha = .7, height = .2) +
    colorspace::scale_color_continuous_sequential("Batlow", 
                                                  rev = F, 
                                                  begin = .1, end = .9,
                                                  trans = "sqrt") +
    labs(x = "Gained range coverage (%)", 
         y = "") +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14) +
    theme(legend.position = "none",
          axis.title = element_text(size = 16)))
ggsave("figures/gained_coverage_percentage_JUN1OCT1.png", width = 6.19, height = 3.29)


# Map --------------------------------------------------------------------------

# map the gains in species richness per cell (i.e., new species-level range coverage)

sr = terra::rast("outputs/richnessmaps/richnessmap_JUN1OCT12025_alltaxa.tif")
