# Script to compare richness maps (to see where more species were recorded in 2025)

# load libraries
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(arrow)

# Input data -------------------------------------------------------------------

# load base grid
base1k = terra::rast("data/base-layers/canada-basegrids/canada.base.1k.tiff")
# read the Canada polygon
canada = sf::read_sf("data/base-layers/canada-polygon/canada.outline.shp")

# load parquet file
df = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# Richness 2008-2024

taxa = unique(count.year.province.group$iconic_taxon_name)
years = 2008:2025 # 2008 is the year iNaturalist launched

# make observations layer
temp = df |>
  dplyr::filter(year < 2025) |>
  select(c(longitude, latitude, year, iconic_taxon_name, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
# reproject
temp = project(temp, crs(base1k))
  
# get cells
temp$cell = extract(base1k, temp, cell=T)[,3]
  
# thin to one species occ per cell (i.e. convert density to just presence/absence)
temp.thin <- temp[-which(duplicated(paste(temp$scientific_name,temp$cell)))]
richmap.pre2025 <- rasterize(temp.thin, base1k, fun="length")
  
# Richness 2008-2025 -----------------------------------------------------------

# make observations layer
temp = df |>
  # dplyr::filter(year < 2025) |> # all years!
  select(c(longitude, latitude, year, iconic_taxon_name, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
# reproject
temp = project(temp, crs(base1k))

# get cells
temp$cell = extract(base1k, temp, cell=T)[,3]

# thin to one species occ per cell (i.e. convert density to just presence/absence)
temp.thin <- temp[-which(duplicated(paste(temp$scientific_name,temp$cell)))]
richmap.at2025 <- rasterize(temp.thin, base1k, fun="length")

# calculate difference
richmap.diff = richmap.at2025 - richmap.pre2025

richmap_stack = rast(list("sr.pre2025" = richmap.pre2025, "sr.at2025" = richmap.at2025, "sr.difference" = richmap.diff))
plot(richmap_stack)
writeRaster(richmap_stack, "outputs/comparison-maps/richness/srchange_alltaxa_2025_vs_20082024.tif", overwrite=TRUE)

