# Script to map new cells in 2025

## minimum year map

# load libraries
library(ggplot2)
library(dplyr)
library(terra)
library(gpkg)
library(plotly)
library(hrbrthemes)
library(tidyterra)

# set common theme for all maps
theme_set( hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                                      axis = FALSE, axis_text_size = 1,
                                      ticks = FALSE,
                                      base_size = 14) +
             theme(legend.position = "top", legend.key.width = unit(2.5,"cm")))

# read the Canada polygon
canada = sf::read_sf("data/base-layers/canada-polygon/canada.outline.shp")

# load density maps
PATH = "outputs/densitymaps/"
FILES = list.files(PATH)[grepl("yearly", list.files(PATH))]
maps = lapply(paste0(PATH, FILES), rast)
names(maps) = gsub("densitymap_yearly_", "", list.files(PATH)[-1])
names(maps) = gsub(".tif", "", names(maps))

# project canada
canada = st_transform(canada, crs(maps[[1]]))

# map
ggplot() +
  geom_sf(data = canada, fill = "grey10", color = "grey10") +
  geom_spatraster(data = maps[[1]], 
                  aes(fill = 2025)) +
  scale_fill_viridis_c(na.value = "transparent", begin = .5, end = 1) 

# extract observations in 2025
maps.2025 = list()
for(i in 1:length(maps)){
  maps.2025[[i]] = maps[[i]]$`2025`
}
names(maps.2025) = names(maps)
maps.2025 = rast(maps.2025)

# aggregate to make them a bit bigger!
maps.2025 = aggregate(maps.2025, factor = 5, FUN = "sum", na.rm = T)
maps.2025 = trim(maps.2025)

# map
ggplot() +
  geom_sf(data = canada, fill = "grey10", color = "grey10") +
  geom_spatraster(data = maps.2025, 
                  aes(fill = alltaxa)) +
  scale_fill_viridis_c(option = "turbo", 
                       na.value = "transparent", 
                       begin = .2, end = .8,
                       trans = "log10",
                       name = "Observations (log10)") 
ggsave("figures/map_observations_2025_alltaxa.png", width = 12.5, height = 7.9)


## species richness version ----------------------------------------------------

# load density maps
PATH = "outputs/richnessmaps/"
FILES = list.files(PATH)[grepl("yearly", list.files(PATH))]
maps = lapply(paste0(PATH, FILES), rast)
names(maps) = gsub("richnessmap_yearly_", "", list.files(PATH))
names(maps) = gsub(".tif", "", names(maps))

# extract observations in 2025
maps.2025 = list()
for(i in 1:length(maps)){
  maps.2025[[i]] = maps[[i]]$`2025`
}
names(maps.2025) = names(maps)
maps.2025 = rast(maps.2025)

# aggregate to make them a bit bigger!
maps.2025 = aggregate(maps.2025, factor = 5, FUN = "sum", na.rm = T)
maps.2025 = trim(maps.2025)

# map
ggplot() +
  geom_sf(data = canada, fill = "grey10", color = "grey10") +
  geom_spatraster(data = maps.2025, 
                  aes(fill = alltaxa)) +
  scale_fill_viridis_c(option = "turbo", 
                       na.value = "transparent", 
                       begin = .2, end = .8,
                       name = "Species observed") 
ggsave("figures/map_observations_2025_alltaxa.png", width = 12.5, height = 7.9)
