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
base1k = terra::rast("data/base-layers/canada-basegrids/canada.base.100k.tiff")
# read the Canada polygon
canada = sf::read_sf("data/base-layers/canada-polygon/canada.outline.shp")

# load parquet file
inat.pre = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_PRE_JUN12025.parquet")
inat.post = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_JUN1OCT1.parquet")
inat.post = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_")

# Richness 2008-2024 ###########################################################

taxa = unique(count.year.province.group$iconic_taxon_name)

# make observations layer
temp = inat.pre |>
  select(c(longitude, latitude, iconic_taxon_name, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
# reproject
temp = project(temp, crs(base1k))

# get cells
temp$cell = extract(base1k, temp, cell=T)[,3]

# thin to one species occ per cell (i.e. convert density to just presence/absence)
temp.thin <- temp[-which(duplicated(paste(temp$scientific_name,temp$cell)))]
richmap.prebtg <- rasterize(temp.thin, base1k, fun="length")

# Richness 2008-2025 -----------------------------------------------------------

# make observations layer
temp2 = inat.post |>
  # dplyr::filter(year < 2025) |> # all years!
  select(c(longitude, latitude, iconic_taxon_name, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
# reproject
temp2 = project(temp2, crs(base1k))

# get cells
temp2$cell = extract(base1k, temp2, cell=T)[,3]

# thin to one species occ per cell (i.e. convert density to just presence/absence)
temp2.thin <- temp2[-which(duplicated(paste(temp2$scientific_name,temp2$cell)))]
richmap.btg <- rasterize(temp2.thin, base1k, fun="length")

## filter out repeated species to only show new taxa added per cell ------------
temp2.thin <- temp2[-which(duplicated(paste(temp2$scientific_name,temp2$cell)))]
index1 = paste(temp2.thin$scientific_name, temp2.thin$cell)
index2 = paste(temp.thin$scientific_name, temp.thin$cell)
index = which(index1 %in% index2)
new.thin <- temp2.thin[-index]
richmap.new <- rasterize(new.thin, base1k, fun="length")
terra::writeRaster(richmap.new, "outputs/richnessmaps/richnessmap_JUN1OCT12025_allNEWtaxa.tif")

# summarise the richness gains
richvals = values(richmap_stack)
richvals = richvals |> na.omit() |> data.frame()
richvals$cellID = terra::cells(richmap.diff)
richxy = terra::xyFromCell(object = richmap.diff, cell = richvals$cellID)
richvals$x = richxy[,1]
richvals$y = richxy[,2]

# make point layer of cells where at least 50% of species were found in 2025
richmap.biggains = richmap_stack$sr.difference/richmap_stack$sr.at2025
richmap.onlygains = richmap.biggains
richmap.onlygains[richmap.onlygains == 0] <- NA
richmap.onlygains.pts = as.points(richmap.onlygains)


# map the gains ----------------------------------------------------------------
theme_set( hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                                      axis = FALSE, axis_text_size = 1,
                                      ticks = FALSE,
                                      base_size = 14) +
             theme(legend.position = "top", legend.key.width = unit(2.5,"cm")))
north = filter(canada, NAME %in% c("Nunavut", 
                                   "Yukon Territory / Territoire du Yukon", 
                                   "Northwest Territories / Territoires du Nord-Ouest",
                                   "Newfoundland and Labrador / Terre-Neuve-et-Labrador"))
north_additions = crop(richmap.new, north, mask = T)
north_additions.pts = as.points(north_additions)

(A = ggplot() +
    tidyterra::geom_spatvector(data = canada, linewidth = 0, fill = "grey90") +
    tidyterra::geom_spatraster(
      data = richmap.new,
      aes(fill = V1_length)
    ) +
    geom_sf(data = filter(north_additions.pts, V1_length >= 100),
            aes(fill = V1_length), pch = 21, size = 1.5) +
    colorspace::scale_fill_continuous_sequential("Batlow", 
                                                 name = "",
                                                 rev = FALSE,
                                                 na.value = "transparent",
                                                 trans = "log10"))
ggsave("figures/gains_speciesrichness_JUN1OCT12025_northernhighlights.png", width = 7.81, height = 7.37)