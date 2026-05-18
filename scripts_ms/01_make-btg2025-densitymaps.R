## Script to make BTG 2025 (JUN1 to OCT1) density maps from the iNaturalist data

# load libraries
library(tidyverse)
library(sf)
library(cartogram)
library(terra)
library(tidyterra)
library(ggplot2)
library(arrow)

# set theme for all plots
# set ggplot theme
theme_set( hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                                      axis = FALSE, 
                                      axis_text_size = 1,
                                      ticks = FALSE,
                                      base_size = 14) +
             theme(legend.position = "top", legend.key.width = unit(2.5,"cm")))

# Input data -------------------------------------------------------------------

# load base grid
base1k = terra::rast("data/base-layers/canada-basegrids/canada.base.1k.tiff") 
# read the Canada polygon
canada = sf::read_sf("data/base-layers/canada-polygon/canada.outline.shp")
# load ecoregions
ecoreg = vect("data/base-layers/ecoregions.shp")
# load canada polygon
canada_poly = terra::vect("data/base-layers/canada-polygon/canada.outline.shp")

# load parquet file
df = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_JUN1OCT1.parquet")

# make maps --------------------------------------------------------------------

# make observations points layer layer
temp = df |>
  select(c(longitude, latitude, year, iconic_taxon_name, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
temp = project(temp, crs(base1k))
# make density raster
densmap = rasterize(temp, base1k, fun = "length")
# write file
terra::writeRaster(densmap, "outputs/densitymaps/densitymap_JUN1OCT12025_alltaxa.tif", overwrite = TRUE)

ggplot() +
  geom_sf(data = canada, col = "grey90") +
  geom_spatraster(data = densmap, aes(fill = V1_length)) +
  colorspace::scale_fill_continuous_sequential("Batlow", 
                                           trans = "log10", 
                                           na.value = "transparent", 
                                           rev = F) +
  labs(fill = "Observations") +
  theme_void() +
  theme(legend.position = "top", 
        legend.key.width = unit(2, "cm"),
        text = element_text(family = "Roboto Condensed"))
ggsave("figures/density_map_JUN1OCT12025_DEC1DATA.png", width = 6.41, height = 5.77)


# density maps by taxa ---------------------------------------------------------

# get list of taxa
taxa.df = df |> group_by(iconic_taxon_name) |> summarise("n" = n()) |> collect()
taxa = na.omit(taxa.df$iconic_taxon_name)
# map per taxa group
densmap_taxa = list()
for(t in 1:length(taxa)){
  
  # make observations layer
  temp = df |>
    dplyr::filter(iconic_taxon_name == taxa[t]) |>
    select(c(longitude, latitude, iconic_taxon_name, scientific_name)) |>
    collect() |>
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
  # project the observations layer
  temp = project(temp, crs(base1k))
  
  map_temp = rasterize(temp, base1k, fun = "length")
  
  terra::writeRaster(map_temp, paste0("outputs/densitymaps/densitymap_JUN1OCT12025_",taxa[t],".tif"), overwrite = TRUE)
}


## make density map of pre-BTG for comparison (everything up until April 30 2025)

# load parquet file
pq = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_PRE_JUN12025.parquet")
# make observations points layer layer
temp = pq |>
  select(c(longitude, latitude, iconic_taxon_name, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
temp = project(temp, crs(base1k))
# make density raster
densmap = rasterize(temp, base1k, fun = "length")
# write file
terra::writeRaster(densmap, "outputs/densitymaps/densitymap_PREJUN12025_alltaxa.tif", overwrite = TRUE)



## plot added density to the map -----------------------------------------------

densmap.btg = rast("outputs/densitymaps/densitymap_JUN1OCT12025_alltaxa.tif")
densmap.prebtg = rast("outputs/densitymaps/densitymap_PREJUN12025_alltaxa.tif")

ggplot() +
  geom_sf(data = canada, col = "grey90") +
  geom_spatraster(data = densmap.btg, aes(fill = V1_length)) +
  colorspace::scale_fill_continuous_sequential("Batlow", 
                                               trans = "log10", 
                                               na.value = "transparent", 
                                               rev = F) +
  labs(fill = "Observations") +
  theme_void() +
  theme(legend.position = "top", 
        legend.key.width = unit(2, "cm"),
        text = element_text(family = "Roboto Condensed"))


# new cells --------------------------------------------------------------------

sampled.btg = densmap.btg
sampled.btg[densmap.btg > 1] <- 1
sampled.btg[is.na(densmap.btg)] <- 0

sampled.pre = densmap.prebtg
sampled.pre[densmap.prebtg > 1] <- 1
sampled.pre[is.na(sampled.pre)] <- 0

sampled.change = sampled.btg - sampled.pre
sampled.new = sampled.change
sampled.new[sampled.change<1] <- NA

# assign densities in BTG to these new cells -----------------------------------
densmap.new = densmap.btg
densmap.new[is.na(sampled.new)] <- NA

# project canada polygon
canada_poly = project(canada_poly, "EPSG:3347")
# make points layer
obsdens.pts = as.points(densmap.new)
obsdens.pts = project(obsdens.pts, "EPSG:3347")
# project densmap
densmap.new = project(densmap.new, "EPSG:3347")

# map the new cells with number of observations --------------------------------

# trickery to make the palette prettier

# first, make a layer of everything < 100 obs
densmap.under100 = filter(densmap.new, V1_length < 100)
densmap.under100[densmap.under100<100] <- 100

ggplot() +
  geom_sf(data = sf::st_as_sf(canada_poly), col = "grey90", fill = "grey90") +
  geom_spatraster(data = densmap.new, 
                  aes(fill = V1_length),
                  interpolate = F) +
  geom_sf(data = filter(obsdens.pts, V1_length >= 99), 
          aes(fill = V1_length, size = V1_length), 
          pch = 21, alpha = .65, linewidth = .05) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                               name = "Observations", 
                                               trans = "sqrt",
                                               na.value = "transparent", 
                                               rev = F, 
                                               breaks = c(100, 500, 1000, 1500)) +
  scale_size_area(name = "", 
                  max_size = 11, 
                  breaks = c(100, 500, 1000),
                  guide = guide_legend(reverse = FALSE, order = 2)) +
  theme_void() +
  theme(text = element_text(family = "Roboto Condensed", size = 14),
    legend.box = "horizontal",
        legend.position = "top",
        legend.key.width = unit(1.3, 'cm')) +
  guides(
    fill = guide_colorbar(order = 1),
    area = guide_legend(order = 2)
  )
ggsave("figures/gain_cells_obsdens_alltaxa_map_JUN1OCT1gains.png", width = 7.45, height = 6.85)



