# Script to make annual species richness maps 

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

# load parquet file
df = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# Count data -------------------------------------------------------------------

count.year = df |>
  group_by(year) |>
  select(year, scientific_name) |>
  distinct() |>
  summarise("n" = n()) |> collect() |> na.omit()
count.year$year = as.numeric(count.year$year)
count.year = filter(count.year, year >= 2008)

count.year.province = df |>
  group_by(year, place_state_name) |>
  select(year, place_state_name, scientific_name) |>
  distinct() |>
  summarise("n" = n()) |> collect() |> na.omit()
count.year.province$year = as.numeric(count.year.province$year)
count.year.province = filter(count.year.province, year >= 2008)

count.year.province.group = df |>
  group_by(year, place_state_name, iconic_taxon_name) |>
  select(year, place_state_name, , iconic_taxon_name, scientific_name) |>
  distinct() |>
  summarise("n" = n()) |> collect() |> na.omit()
count.year.province.group$year = as.numeric(count.year.province.group$year)
count.year.province.group = filter(count.year.province.group, year >= 2008)

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

richmap = list()
for(y in 1:length(years)){
  
  # make observations layer
  temp = df |>
    dplyr::filter(year == years[y]) |>
    select(c(longitude, latitude, year, iconic_taxon_name, scientific_name)) |>
    collect() |>
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
  
  temp = project(temp, crs(base1k))

  # get cells
  temp$cell = extract(base1k, temp, cell=T)[,3]
  
  # thin to one species occ per cell (i.e. convert density to just presence/absence)
  temp.thin <- temp[-which(duplicated(paste(temp$scientific_name,temp$cell)))]
  richmap[[y]] <- rasterize(temp.thin, base1k, fun="length")

}
# stack the maps and save
map_stack = richmap
map_stack = rast(richmap)
names(map_stack) = years
terra::writeRaster(map_stack, "outputs/richnessmaps/richnessmap_yearly_alltaxa.tif", overwrite = TRUE)

# map per taxa group
richmap_taxa = list()
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
    
    # get cells
    temp2$cell = extract(base1k, temp2, cell=T)[,3]
    
    # thin to one species occ per cell (i.e. convert density to just presence/absence)
    temp2.thin <- temp2[-which(duplicated(paste(temp2$scientific_name,temp2$cell)))]
    map_temp[[y]] = rasterize(temp2.thin, base1k, fun="length")
  }
  
  names(map_temp) = years
  richmap_taxa[[y]] = rast(map_temp)
  terra::writeRaster(richmap_taxa[[y]], paste0("outputs/richnessmaps/richnessmap_yearly_",taxa[t],".tif"), overwrite = TRUE)
}