# Calculate change in range coverage -------------------------------------------

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
library(mapview)

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# load range coverage function
source("scripts/range-coverage/calc_range_coverage_change.R")

# list of taxa groups
taxagroups = list.files("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/")
inatgroups = c("Amphibia", "Aves", "Insecta", "Mammalia", "Plantae", "Reptilia")


# Load the databases ===========================================================

canada_poly = terra::vect("~/Documents/GitHub/ciee/blitz-the-gap/00_rawdata/base-layers/canada-polygon/canada.outline.shp")
canada_poly = terra::project(canada_poly, "EPSG:6933")

# load canada base grid
canada = terra::rast("~/Documents/GitHub/ciee/blitz-the-gap/00_rawdata/base-layers/canada.base.5k.tiff")
canada = terra::aggregate(canada, factor = 2) # 10 km^2 cells
canada = terra::project(canada, "EPSG:6933")

# load parquet
inat_pq <- arrow::open_dataset("~/McGill University/Laura's Lab_Group - BioBlitz/data/raw/biodiversity-data/inat-canada/iNat_non_sensitive_data_Jan2025.parquet")

# list of species to calculate range coverage for
df = lapply(as.list(paste0("outputs/range-coverage/range_coverage_", taxagroups, ".rds")), readRDS)
names(df) = taxagroups
df = bind_rows(df, .id = "Taxa")

for(t in c(1:4,6)){

  # Load species names that have previously been included
  species = dplyr::filter(df, Taxa == taxagroups[t]) |>
    select(species) |>
    as.matrix() |> as.vector()
  
  # Prepare to load IUCN ranges
  filepath = paste0("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/", taxagroups[t], "/")
  
  for(i in 1:length(species)){
    tryCatch({  
    temp = calc_range_coverage_change(species_name = species[i], 
                               range = terra::rast(paste0(filepath, species[i], ".tif")),
                               inat = inat_pq#, 
                               #new_data = newdata
                               )
    
    terra::writeRaster(temp, 
                 paste0("~/McGill University/Laura's Lab_Group - range-coverage/",taxagroups[t],"/", species[i], ".tif"), overwrite=T)
    }, error = function(e) {
      message("An error occurred: ", e$message)
      return(NA) # Return NA if an error occurs
    })
  }
}
# 
# 
# 
# rc_before = lapply(paste0("outputs/range-coverage/range_coverage_", taxagroups[c(1:4,6)], ".rds"), readRDS)
# rc_after = lapply(paste0("outputs/range-coverage/after/range_coverage_", taxagroups[c(1:4,6)], ".rds"), readRDS)
# names(rc_before) = taxagroups[c(1:4,6)]
# names(rc_after) = taxagroups[c(1:4,6)]
# rc_before = rc_before |> bind_rows(.id = "Taxa")
# rc_after = rc_after |> bind_rows(.id = "Taxa")
# rc_compare = left_join(rc_before, rc_after, 
#                        by = c("Taxa", "species", "range_size_canada", "range_size_full"), 
#                        suffix = c(".before", ".after"))
# 
# 
# ggplot(data = rc_compare) +
#   geom_point(aes(col = Taxa, 
#                 y = (coverage_min1.after-coverage_min1.before)/range_size_canada,
#                 x = Taxa))

## reminder: code the plants version