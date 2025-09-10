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

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# load range coverage function
source("scripts/range-coverage/calc_range_coverage_change.R")

# list of taxa groups
taxagroups = list.files("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/")
inatgroups = c("Amphibia", "Aves", "Insecta", "Mammalia", "Plantae", "Reptilia")


# connect to sharepoint ========================================================

# sign in to authenticate
odb <- get_business_onedrive()
# create folder (only run once!)
# odb$create_folder("BTG_analyses/range-coverage-rasters/")
# odb$create_folder("BTG_analyses/range-coverage-rasters/AMPHIBIANS/")


# Load the databases ===========================================================

canada_poly = terra::vect("~/Documents/GitHub/ciee/blitz-the-gap/00_rawdata/base-layers/canada-polygon/canada.outline.shp")

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


for(t in c(2:4,6)){ #c(2:4,6)){

  # Load species names that have previously been included
  species = dplyr::filter(df, Taxa == taxagroups[t]) |>
    select(species) |>
    as.matrix() |> as.vector()
  
  # Prepare to load IUCN ranges
  filepath = paste0("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/", taxagroups[t], "/")
  
  # Prepare data frame to store results
  range_coverage = data.frame("species" = species,
                              "coverage_min1" = NA,
                              "coverage_min3" = NA,
                              "coverage_min10" = NA,
                              "range_size_canada" = NA,
                              "range_size_full" = NA)
  
  for(i in 1:length(species)){
    tryCatch({  
    temp = calc_range_coverage_change(species_name = species[i], 
                               range = terra::rast(paste0(filepath, species[i], ".tif")),
                               inat = inat_pq)
    range_coverage$coverage_min1[i] = temp$coverage[1]
    range_coverage$coverage_min3[i] = temp$coverage[2]
    range_coverage$coverage_min10[i] = temp$coverage[3]
    range_coverage$range_size_canada[i] = temp$n_cells_canada
    range_coverage$range_size_full[i] = temp$n_cells_fullrange
    
    terra::writeRaster(temp$rasters, 
                 paste0("~/McGill University/Laura's Lab_Group - range-coverage/",taxagroups[t],"/", species[i], ".tif"), overwrite=T)
    }, error = function(e) {
      message("An error occurred: ", e$message)
      return(NA) # Return NA if an error occurs
    })
  }
  saveRDS(range_coverage, 
          paste0("outputs/range-coverage/after/range_coverage_", taxagroups[t], ".rds"))
}


# function to map the change in cells that meet a threshold
cvg_raster = function(before, after, threshold = 3){
  
  after_min = after
  after_min[after >= threshold] <- 1
  after_min[after < threshold] <- 0
  
  before_min = before
  before_min[before >= threshold] <- 1
  before_min[before < threshold] <- 0
  
  temp = c(before_min, after_min, after_min-before_min)
  names(temp) = c("before", "after", "change")
  return(temp)
}

temp = cvg_raster(temp2$before, temp2$after)
plot(temp)
sum(values(temp$after), na.rm = T)
sum(values(temp$before), na.rm = T)
sum(values(temp$change), na.rm = T)
