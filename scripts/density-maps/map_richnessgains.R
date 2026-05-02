# Script to show cells that were newly sampled this year

# load libraries
library(ggplot2)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(dplyr)
library(sf)

# Canada polygon
canada <- ne_countries(scale = "medium", country = "Canada", returnclass = "sf")
canada <- st_transform(canada, crs = "EPSG:3347")

# load new cell file
newcells = terra::rast("outputs/comparison-maps/density/newcells_alltaxa_2025_vs_20082024.tif")
newSR = terra::rast("outputs/comparison-maps/richness/srchange_alltaxa_2025_vs_20082024.tif")

# project
newcells = project(newcells, crs(canada))
newSR = project(newSR, crs(canada))
newSR = trim(newSR)
Sr = newSR
newSR = filter(newSR, sr.difference > 0)

# convert newnewSR_pts# convert newcells to points
newcells_pts = as.points(newcells$alltaxa)
newcells_pts = st_as_sf(newcells_pts)

# convert newSR to points
newSR_pts = as.points(newSR$sr.difference)
newSR_pts = st_as_sf(newSR_pts)
newSR_pts = newSR_pts |> filter(sr.difference>0)

# resample for nicer visualisation
newcells_big = aggregate(newcells, factor = 25, fun = "max", na.rm = T)
newSR_big = aggregate(newcells, factor = 50, fun = "max", na.rm = T)

# make plotting function
plot_fun = function(dataset, type = "points"){
  
  if(type == "raster"){
  
  p = ggplot() +
    geom_sf(data = canada, fill = "grey10", col = "grey10") +
    geom_spatraster(data = dataset, aes(fill = sr.difference)) +
    scale_fill_viridis_c(option = "turbo",
                         name = "Gain in\nobserved\nrichness",
                         trans = "sqrt",
                         limits = c(1,1000),
                         end = .8, begin = .2,
                         na.value = "transparent") +
    hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                               axis = FALSE, axis_text_size = 1,
                               ticks = FALSE,
                               base_size = 14) +
    theme(legend.position = "top", legend.key.width = unit(2,"cm"))
  } else {
    
    p = ggplot() +
      geom_sf(data = canada, fill = "grey10", col = "grey10") +
      geom_sf(data = dataset, aes(fill = sr.difference), 
              size = 1.5, pch = 21, col = "grey10", linewidth = .1) +
      scale_fill_viridis_c(option = "turbo",
                           name = "Gain in\nobserved\nrichness",
                           trans = "sqrt",
                           limits = c(1,1000),
                           end = .8, begin = .2,
                           na.value = "transparent") +
      hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                                 axis = FALSE, axis_text_size = 1,
                                 ticks = FALSE,
                                 base_size = 14) +
      theme(legend.position = "top", legend.key.width = unit(2,"cm"))
      
    }
  return(p)
  
}

plot_fun(newSR, type = "raster")
ggsave("figures/gain_observedrichness_map.png", width = 12.5, height = 7.9)

plot_fun(filter(newSR_pts, sr.difference >= 50), type = "points")
ggsave("figures/gain_observedrichness_map_atleast50new.png", width = 12.5, height = 7.9)

plot_fun(filter(newSR_pts, sr.difference >= 100), type = "points")
ggsave("figures/gain_observedrichness_map_atleast100new.png", width = 12.5, height = 7.9)

plot_fun(filter(newSR_pts, sr.difference >= 300), type = "points")
ggsave("figures/gain_observedrichness_map_atleast300new.png", width = 12.5, height = 7.9)




## inspect these high-gain cells: are they where people went for bioblitzes?
biggains = filter(newSR, sr.difference >= 100)

mapview::mapview(biggains)

## geog distance and clim distance of these high-gain cells from the previous nearest neighbours

