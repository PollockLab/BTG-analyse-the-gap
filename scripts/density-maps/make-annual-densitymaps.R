## Script to make annual density maps from the iNaturalist data

# load libraries
library(tidyverse)
library(sf)
library(cartogram)
library(terra)
library(tidyterra)
library(ggplot2)
library(arrow)

# set theme for all plots
theme_set(hrbrthemes::theme_ipsum_rc())

# Input data -------------------------------------------------------------------

# load base grid
base1k = terra::rast("data/base-layers/canada-basegrids/canada.base.1k.tiff") 
# read the Canada polygon
canada = sf::read_sf("data/base-layers/canada-polygon/canada.outline.shp")
# load ecoregions
ecoreg = vect("data/base-layers/ecoregions.shp")

# load parquet file
df = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# Count data -------------------------------------------------------------------

# count.year = df |>
#   group_by(year) |>
#   summarise("n" = n()) |> collect() |> na.omit()
# count.year$year = as.numeric(count.year$year)
# count.year = filter(count.year, year >= 2008)
# 
# count.year.province = df |>
#   group_by(year, place_state_name) |>
#   summarise("n" = n()) |> collect() |> na.omit()
# count.year.province$year = as.numeric(count.year.province$year)
# count.year.province = filter(count.year.province, year >= 2008)
# 
# count.year.province.group = df |>
#   group_by(year, place_state_name, iconic_taxon_name) |>
#   summarise("n" = n()) |> collect() |> na.omit()
# count.year.province.group$year = as.numeric(count.year.province.group$year)
# count.year.province.group = filter(count.year.province.group, year >= 2008)

# ggplot(data = count.year) +
#   geom_line(aes(x = year, y = n))
# ggplot(data = count.year.province) +
#   geom_line(aes(x = year, y = n, col = place_state_name)) +
#   scale_y_sqrt()
# ggplot(data = count.year.province.group) +
#   geom_line(aes(x = year, y = n, col = place_state_name)) +
#   scale_y_sqrt() +
#   facet_wrap(~iconic_taxon_name)


# make some maps ---------------------------------------------------------------

taxa = unique(count.year.province.group$iconic_taxon_name)
years = 2008:2025 # 2008 is the year iNaturalist launched

densmap = list()
for(y in 1:length(years)){
  
  # make observations layer
  temp = df |>
    dplyr::filter(year == years[y]) |>
    select(c(longitude, latitude, year, iconic_taxon_name, scientific_name)) |>
    collect() |>
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
  
  temp = project(temp, crs(base1k))
  densmap[[y]] = rasterize(temp, base1k, fun = "length")
}
# stack the maps and save
map_stack = densmap
map_stack = rast(densmap)
names(map_stack) = years
terra::writeRaster(map_stack, "outputs/densitymaps/densitymap_yearly_alltaxa.tif", overwrite = TRUE)

# map per taxa group
densmap_taxa = list()
for(t in 1:length(taxa)){
  
  # make observations layer
  temp = df |>
    dplyr::filter(iconic_taxon_name == taxa[t]) |>
    select(c(longitude, latitude, year, iconic_taxon_name, scientific_name)) |>
    collect() |>
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
  # project the observations layer
  temp = project(temp, crs(base1k))
  
  # run yearly
  map_temp = list()
  for(y in 1:length(years)){
    temp2 = temp |>
      filter(year == years[y])
    map_temp[[y]] = rasterize(temp2, base1k, fun = "length")
  }
  names(map_temp) = years
  densmap_taxa[[y]] = rast(map_temp)
  terra::writeRaster(densmap_taxa[[y]], paste0("outputs/densitymaps/densitymap_yearly_",taxa[t],".tif"), overwrite = TRUE)
}

# make total density map
temp = df |>
  select(c(longitude, latitude, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
temp = project(temp, crs(base1k))
densmap_total = rasterize(temp, base1k, fun = "length")
terra::writeRaster(densmap_total, paste0("outputs/densitymaps/densitymap_20082025_alltaxa.tif"), overwrite = TRUE)
