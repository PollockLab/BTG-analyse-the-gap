# Script to compare richness maps (to see where more species were recorded in 2025)

# load libraries
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(arrow)

theme_set(hrbrthemes::theme_ipsum_rc())

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

richmap_stack = rast(list("sr.pre2025" = richmap.pre2025, 
                          "sr.at2025" = richmap.at2025, 
                          "sr.difference" = richmap.diff))
plot(richmap_stack)
writeRaster(richmap_stack, "outputs/comparison-maps/richness/srchange_alltaxa_2025_vs_20082024.tif", overwrite=TRUE)

# summarise the richness gains
richvals = values(richmap_stack)
richvals = richvals |> na.omit() |> data.frame()
richvals$cellID = terra::cells(richmap.diff)
richvals$x = richxy[,1]
richvals$y = richxy[,2]

richxy = terra::xyFromCell(object = richmap.diff, cell = richvals$cellID)

# make point layer of cells where at least 50% of species were found in 2025
richmap.biggains = richmap.diff/richmap.at2025
richmap.onlygains = richmap.biggains
richmap.onlygains[richmap.onlygains == 0] <- NA
richmap.onlygains.pts = as.points(richmap.onlygains)

# map the gains ----------------------------------------------------------------

ggplot() +
  tidyterra::geom_spatvector(data = canada, linewidth = 0) +
  geom_spatvector(data = richmap.onlygains.pts,
                  aes(col = V1_length)) +
  colorspace::scale_color_continuous_sequential("Batlow", rev = FALSE) +
  labs(col = "Relative\ngain (%)")
ggsave("figures/gained_richness/map_relativegains_allgains.png", 
       width = 8.54, height = 7.33)

ggplot() +
  tidyterra::geom_spatvector(data = canada, linewidth = 0) +
  geom_spatvector(data = filter(richmap.onlygains.pts, V1_length >= 0.75),
                  aes(col = V1_length)) +
  colorspace::scale_color_continuous_sequential("Batlow", 
                                                rev = FALSE, begin = 0.75) +
  labs(col = "Relative\ngain (%)")
ggsave("figures/gained_richness/map_relativegains_75pgains.png", 
       width = 8.54, height = 7.33)

ggplot() +
  tidyterra::geom_spatvector(data = canada, linewidth = 0) +
  geom_spatvector(data = filter(richmap.onlygains.pts, V1_length > 0.9),
                  aes(col = V1_length)) +
  colorspace::scale_color_continuous_sequential("Batlow", 
                                                rev = FALSE, begin = 0.9) +
  labs(col = "Relative\ngain (%)") 
ggsave("figures/gained_richness/map_relativegains_90pgains.png", 
       width = 8.54, height = 7.33)

# plot the gains ---------------------------------------------------------------

richgains = richvals |> filter(sr.difference > 0)

ggplot(data = richgains) +
  geom_histogram(aes(x = 100*sr.difference/sr.at2025), 
               fill = "green4", alpha = .4, linewidth = .1) +
  labs(x = "Proportion of total observed species that were found in 2025 (%)", 
       y = "Density (sites)")
ggsave("figures/gained_richness/sr_diffVStotal_histogram.png", 
       width = 7.3, height = 4.6)

ggplot(data = richgains) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, col = "grey") +
  geom_point(aes(y = sr.difference, x = sr.at2025, 
                 col = 100*sr.difference/sr.at2025), 
             alpha = .5, size = 2) +
  scale_color_viridis_c(option = "turbo") +
  labs(y = "New species found in 2025", 
       x = "Observed species richness",
       col = "Relative\ngain (%)") +
  theme(legend.position = "bottom")
ggsave("figures/gained_richness/sr_diffVStotal_progressrainbow.png", 
       width = 6.3, height = 6.5)

ggplot(data = filter(richgains, sr.difference >=100)) +
  geom_histogram(aes(x = sr.difference), bins = 20) +
  labs(x = "Gained species richness", y = "Number of cells") +
  scale_y_sqrt()
