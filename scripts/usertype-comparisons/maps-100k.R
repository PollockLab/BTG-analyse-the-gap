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
df = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")


# Richness 2008-2024 ###########################################################

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
writeRaster(richmap_stack, "outputs/comparison-maps/richness/srchange_alltaxa_2025_vs_20082024_100k.tif", overwrite=TRUE)
richmap_stack = terra::rast("outputs/comparison-maps/richness/srchange_alltaxa_2025_vs_20082024_100k.tif")
richmap.diff = richmap_stack$sr.difference

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
north_additions = crop(richmap_stack$sr.difference, north, mask = T)
north_additions.pts = as.points(north_additions)

(A = ggplot() +
  tidyterra::geom_spatvector(data = canada, linewidth = 0, fill = "grey90") +
  tidyterra::geom_spatraster(
    data = filter(richmap_stack$sr.difference, sr.difference > 0),
    aes(fill = sr.difference)
  ) +
    geom_sf(data = filter(north_additions.pts, sr.difference >= 100),
            aes(fill = sr.difference), pch = 21, size = 1.5) +
  colorspace::scale_fill_continuous_sequential("Batlow", 
                                               name = "",
                                               rev = FALSE,
                                               na.value = "transparent",
                                               trans = "log10"))
ggsave("figures/gains_speciesrichness_northernhighlights.png", width = 7.81, height = 7.37)

(B = ggplot() +
    tidyterra::geom_spatvector(data = canada, linewidth = 0, fill = "grey90") +
    geom_sf(
      data = north_additions.pts,
      aes(col = sr.difference), col = "transparent"
    ) +
    geom_sf(data = filter(north_additions.pts, sr.difference > 100),
            aes(col = sr.difference)) +
    colorspace::scale_color_continuous_sequential("Batlow", 
                                                 name = "Increased observed\nspecies richness",
                                                 rev = FALSE,
                                                 na.value = "transparent",
                                                 trans = "log10") #+
    #coord_sf(xlim = c(-100000, 350000), ylim = c(2000000,4275000))
  )


## DENSITY #####################################################################

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
terra::writeRaster(map_stack, "outputs/densitymaps/densitymap_yearly_alltaxa_100k.tif", overwrite = TRUE)

# separate 
dens.2025 = map_stack$`2025`
dens.pre2025 = map_stack[[-18]]
dens.diff = dens.2025
