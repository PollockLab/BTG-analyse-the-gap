# Script to calculate the range coverage of a species as the number of
# 10km2 cells that have at least 1, 3, or 10 observations in them

# load libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(terra)
library(gpkg)
library(taxize)
library(BIEN)

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# load range coverage function
source("scripts/range-coverage/calc_range_coverage.R")

# Load the databases ===========================================================

canada_poly = terra::vect("~/Documents/GitHub/ciee/blitz-the-gap/00_rawdata/base-layers/canada-polygon/canada.outline.shp")

# load canada base grid
canada = terra::rast("~/Documents/GitHub/ciee/blitz-the-gap/00_rawdata/base-layers/canada.base.5k.tiff")
canada = terra::aggregate(canada, factor = 2) # 10 km^2 cells

# load parquet
inat_pq <- arrow::open_dataset("~/McGill University/Laura's Lab_Group - BioBlitz/data/raw/biodiversity-data/inat-canada/iNat_non_sensitive_data_Jan2025.parquet")
inat_pq$schema

# set projections
canada_poly = terra::project(canada_poly, crs(mammals))
canada = terra::project(canada, crs(canada_poly))

# list of taxa groups
taxagroups = list.files("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/")
inatgroups = c("Amphibia", "Aves", "Insecta", "Mammalia", "Plantae", "Reptilia")

# Calculate range coverage
for(t in c(2:4,6)){

# get list of species in iNat Canada with at least 3 observations
query <- inat_pq |> 
  # filter by species group:
  filter(iconic_taxon_name == inatgroups[t],
         captive_cultivated == FALSE) |>
  # Summarise number of obs per species, per species group
  group_by(iconic_taxon_name, scientific_name) |>
  summarize(total_obs = n()) |> 
  filter(total_obs >= 3) |>
  # load the query into our R session
  collect()

# count words in each string to find the species-level ones and make clean sp list
species = query$scientific_name[which(lengths(strsplit(query$scientific_name, "\\W+"))==2)]

# list iucn range polygons that have a match
filepath = paste0("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/", taxagroups[t], "/")
ranges_list = gsub(".tif", "", list.files(filepath))
ranges_inat = ranges_list[which(ranges_list %in% species)]

range_coverage = data.frame("species" = ranges_inat,
                            "coverage_min1" = NA,
                            "coverage_min3" = NA,
                            "coverage_min10" = NA,
                            "range_size_canada" = NA,
                            "range_size_full" = NA)

for(i in 1:length(ranges_inat)){
  temp = calc_range_coverage(species_name = ranges_inat[i], 
                             range = terra::rast(paste0(filepath, ranges_inat[i], ".tif")),
                             inat = query)
  range_coverage$coverage_min1[i] = temp$coverage[1]
  range_coverage$coverage_min3[i] = temp$coverage[2]
  range_coverage$coverage_min10[i] = temp$coverage[3]
  range_coverage$range_size_canada[i] = temp$n_cells_canada
  range_coverage$range_size_full[i] = temp$n_cells_fullrange
  
}

saveRDS(range_coverage, paste0("outputs/range-coverage/range_coverage_", taxagroups[t], ".rds"))

}

## Plot! =======================================================================

df = lapply(as.list(paste0("outputs/range-coverage/range_coverage_", taxagroups[-5], ".rds")), readRDS)
names(df) = taxagroups[-5]
df = bind_rows(df, .id = "Taxa")

ggplot(data = df) +
  geom_jitter(aes(col = Taxa, 
                 x = 100*coverage_min1/range_size_canada,
                 y = Taxa)) +
  theme(legend.position = "none") +
  scale_x_continuous() +
  labs(x = "Range coverage", y = "")


## Plants - BIEN ranges --------------------------------------------------------

# can_list = BIEN::BIEN_list_country("Canada")

# get list of species in iNat Canada with at least 3 observations
query <- inat_pq |> 
  # filter by species group:
  filter(iconic_taxon_name == "Plantae",
         captive_cultivated == FALSE) |>
  # Summarise number of obs per species, per species group
  group_by(iconic_taxon_name, scientific_name) |>
  summarize(total_obs = n()) |> 
  filter(total_obs >= 3) |>
  # load the query into our R session
  collect()

# count words in each string to find the species-level ones and make clean sp list
species = query$scientific_name[which(lengths(strsplit(query$scientific_name, "\\W+"))==2)]

# Run for all species
for(sp in 1:length(species)){
  
  # Load the range map from BIEN
  
  range_sf <- BIEN_ranges_load_species(species[sp])
  
  # go to next species if there isn't a range polygon
  if(nrow(range_sf) == 0){ next
  } else {
    range_rast = range_sf |>
      terra::vect() |> 
      terra::project(crs(canada)) |>
      terra::rasterize(canada)
  }
  
  # Calculate range coverage
  temp = calc_range_coverage(species_name = species[sp], 
                             range = range_rast,
                             inat = query)
  # save as a new row to append
  new_row = c(species[sp], temp$coverage[1], temp$coverage[2], temp$coverage[3], 
              temp$n_cells_canada)
  
  # start data frame OR append to existing data frame
  if(sp == 1){
    range_coverage = data.frame(t(new_row))
    colnames(range_coverage) = c("species",
                                "coverage_min1",
                                "coverage_min3",
                                "coverage_min10",
                                "range_size_canada")
  } else {
  range_coverage = rbind(range_coverage, new_row)
  }
  
  # save every 500 species to avoid losing progress
  if(sp %in% seq(500, 7000, 500)){
    saveRDS(range_coverage, paste0("outputs/range-coverage/range_coverage_PLANTS.rds"))
  }
}
# final save
saveRDS(range_coverage, paste0("outputs/range-coverage/range_coverage_", taxagroups[t], ".rds"))





query$scientific_name[1] %in% species

(range_bien <- BIEN_ranges_load_species(query$scientific_name[3])) 
range_bien |> nrow()
plot(range_bien[1])
