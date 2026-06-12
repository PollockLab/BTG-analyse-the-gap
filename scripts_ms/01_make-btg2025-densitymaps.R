## Script to make BTG 2025 (JUN1 to OCT1) density maps from the iNaturalist data

# load libraries
library(tidyverse)
library(sf)
library(cartogram)
library(terra)
library(tidyterra)
library(ggplot2)
library(arrow)

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

# map per taxa group
densmap_taxa = list()
for(t in 1:length(taxa)){
  
  # make observations layer
  temp = pq |>
    dplyr::filter(iconic_taxon_name == taxa[t]) |>
    select(c(longitude, latitude, iconic_taxon_name, scientific_name)) |>
    collect() |>
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
  # project the observations layer
  temp = project(temp, crs(base1k))
  
  map_temp = rasterize(temp, base1k, fun = "length")
  
  terra::writeRaster(map_temp, paste0("outputs/densitymaps/densitymap_PREJUN12025_",taxa[t],".tif"), overwrite = TRUE)
}

## ALL TAXA --------------------------------------------------------------------

# change in density ------------------------------------------------------------

densmap.btg = rast("outputs/densitymaps/densitymap_JUN1OCT12025_alltaxa.tif")
densmap.prebtg = rast("outputs/densitymaps/densitymap_PREJUN12025_alltaxa.tif")

densmap.change = densmap.btg - densmap.prebtg
densmap.change[densmap.change<0] <- 0 # change all cells that did not gain anything to 0
writeRaster(densmap.change, "outputs/comparison-maps/density/changecells_alltaxa_JUN1OCT12025.tif")

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
writeRaster(densmap.new, "outputs/comparison-maps/density/newcells_alltaxa_JUN1OCT12025.tif")

# prepare spatial layers for mapping -------------------------------------------

# project 
canada_poly = project(canada_poly, "EPSG:3347")
densmap.new = project(densmap.new, "EPSG:3347")
# make points layer
obsdens.pts = as.points(densmap.new)
obsdens.pts = project(obsdens.pts, "EPSG:3347")


# map the new cells with number of observations --------------------------------

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


## Calculate new spatial coverage (area of new cells) --------------------------

densmap.new = project(densmap.new, base1k)
# get cell area
cell_area = cellSize(densmap.new, unit = "km", mask = TRUE) |>
  values() |> na.omit() |> as.vector()
# summarise the new cells
res_df = data.frame(
  "area_km_newcells" = sum(cell_area),
  "n_newcells" = length(cell_area)
)

# calculate area of Canada base map
base1k_area = cellSize(base1k, unit = "km")
base1k_area[is.na(base1k)] <- NA
can_area = values(base1k_area) |> sum(na.rm = T)
res_df$area_km_newcells*100/can_area

# proportion of gained coverage
newcells_dens = values(densmap.new) |> na.omit()
length(newcells_dens)
100*length(which(newcells_dens==1))/length(newcells_dens) 
100*length(which(newcells_dens>1))/length(newcells_dens) 
100*length(which(newcells_dens>=10))/length(newcells_dens)  
100*length(which(newcells_dens>=100))/length(newcells_dens)  
100*length(which(newcells_dens>=500))/length(newcells_dens)  

## PER TAXA --------------------------------------------------------------------

# load all the rasters ---------------------------------------------------------

# pre Blitz the Gap (name the layers by taxa)
densmaps.prebtg = lapply(taxa,
                         function(x) { rast(paste0("outputs/densitymaps/densitymap_PREJUN12025_",x,".tif"))})
names(densmaps.prebtg) = taxa
# stack the rasters
densmaps.prebtg = rast(densmaps.prebtg)

# during Blitz the Gap (name the layers by taxa)
densmaps.btg = lapply(taxa,
                         function(x) { rast(paste0("outputs/densitymaps/densitymap_JUN1OCT12025_",x,".tif"))})
names(densmaps.btg) = taxa
# stack the rasters
densmaps.btg = rast(densmaps.btg)

# change in density ============================================================

densmaps.change = densmaps.btg - densmaps.prebtg
densmaps.change[densmaps.change<0] <- 0 # change all cells that did not gain anything to 0
writeRaster(densmap.change, "outputs/comparison-maps/density/changecells_pertaxa_JUN1OCT12025.tif", overwrite=TRUE)

# new cells --------------------------------------------------------------------

sampled.btg = densmaps.btg
sampled.btg[densmaps.btg > 1] <- 1
sampled.btg[is.na(densmaps.btg)] <- 0

sampled.pre = densmaps.prebtg
sampled.pre[densmaps.prebtg > 1] <- 1
sampled.pre[is.na(densmaps.prebtg)] <- 0

sampled.change = sampled.btg - sampled.pre
sampled.new = sampled.change
sampled.new[sampled.change<1] <- NA

# assign densities in BTG to these new cells -----------------------------------

densmaps.new = densmaps.btg
densmaps.new[is.na(sampled.new)] <- NA
writeRaster(densmaps.new, "outputs/comparison-maps/density/newcells_pertaxa_JUN1OCT12025.tif", overwrite = T)

## Calculate new spatial coverage (area of new cells) --------------------------

densmaps.new = project(densmaps.new, base1k)
# get cell area
cell_area = lapply(densmaps.new, cellSize, unit = "km", mask = TRUE) |>
  lapply(values) |> 
  lapply(na.omit) |> 
  lapply(as.vector)
names(cell_area) = names(densmaps.new)
area_df = cell_area |>
  lapply(as.data.frame) |>
  bind_rows(.id = "taxa")
colnames(area_df)[2] = "cell_area"

# summarise the new cells
res_df = area_df |>
  group_by(taxa) |>
  summarise(
    "area_km_newcells" = sum(cell_area),
    "n_newcells" = n()
  )
res_df$prop_canadaarea = 100*res_df$area_km_newcells/can_area
write.csv(res_df, "outputs/spatial-coverage-per-taxa.csv")

# plot the spatial coverage gains per group ------------------------------------

res_df$taxa = factor(res_df$taxa, levels = res_df$taxa[order(res_df$area_km_newcells)])
ggplot(data = filter(res_df, taxa != "alltaxa")) +
  geom_bar(aes(x = area_km_newcells, y = taxa, fill = area_km_newcells), 
           stat = "identity") +
  scale_fill_viridis_c(option = "viridis", trans = "sqrt", end = .95) +
  labs(x = "Coverage gain (km²)", y = "Iconic Taxon Groups") +
  hrbrthemes::theme_ipsum_rc() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "none") 
ggsave("figures/gain_totalrangecoverage_area_barplot.png", width = 6.39, height = 3.05)


## Density map of 2024 (BTG window: June 1 to Oct 1) for comparisons ===========

pq = arrow::read_parquet("data/heavy/BTG-data/inaturalist-canada-dec2025_PRE_JUN12025.parquet")
v24 = pq |>
  dplyr::filter(year == 2024, month %in% 6:10)
# make observations points layer layer
temp = v24 |>
  select(c(longitude, latitude, year, iconic_taxon_name, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
temp = project(temp, crs(base1k))
# make density raster
densmap = rasterize(temp, base1k, fun = "length")
# write file
terra::writeRaster(densmap, "outputs/densitymaps/densitymap_JUN1OCT12024_alltaxa.tif", overwrite = TRUE)

