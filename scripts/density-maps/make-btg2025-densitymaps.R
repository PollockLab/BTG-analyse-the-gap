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
theme_set(hrbrthemes::theme_ipsum_rc())

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
canada_poly = project(canada_poly, crs(densmap.new))

obsdens.pts = as.points(densmap.new)

# set ggplot theme
theme_set( hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                                      axis = FALSE, 
                                      axis_text_size = 1,
                                      ticks = FALSE,
                                      base_size = 14) +
             theme(legend.position = "top", legend.key.width = unit(2.5,"cm")))

ggplot() +
  geom_sf(data = sf::st_as_sf(canada_poly), col = "grey90", fill = "grey90") +
  geom_spatraster(data = filter(densmap.new, V1_length < 100),
                  aes(fill = V1_length),
                  interpolate = F) +
  geom_sf(data = filter(obsdens.pts, V1_length >=50), 
          aes(fill = V1_length, size = V1_length), 
          pch = 21, alpha = .8, linewidth = .1) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                               name = "Observations", 
                                               trans = "log10",
                                               na.value = "transparent", rev = F) +
  scale_size_area(name = "", 
                  max_size = 11, 
                  breaks = c(100, 250, 500, 750, 1000),
                  guide = guide_legend(reverse = TRUE)) +
  theme(legend.box = "horizontal",
        legend.position = "right",
        legend.key.width = unit(0.5, 'cm'),
        legend.key.height = unit(1.25, 'cm'),
        legend.spacing.x = unit(0.2, 'cm'))
ggsave("figures/gain_cells_obsdens_alltaxa_map.png", width = 7.57, height = 6.35)



