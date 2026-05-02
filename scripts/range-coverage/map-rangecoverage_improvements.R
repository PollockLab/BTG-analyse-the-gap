# Script to map the cells where range coverage improved for X species

# load libraries
library(ggplot2)
library(dplyr)
library(terra)
library(gpkg)
library(plotly)
library(hrbrthemes)
library(tidyterra)

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# load canada base grid
canada = terra::rast("data/base-layers/canada-basegrids/canada.base.5k.tiff")
canada = terra::aggregate(canada, factor = 2) # 10 km^2 cells
canada = project(canada, "EPSG:3347")

# load canada polygon
canada_poly = terra::vect("data/base-layers/canada-polygon/canada.outline.shp")
canada_poly = project(canada_poly, "EPSG:3347")

# list of taxa groups
taxagroups = list.files("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/")
inatgroups = c("Amphibia", "Aves", "Insecta", "Mammalia", "Plantae", "Reptilia")

# let's not do birds:
taxagroups = taxagroups[-2]
inatgroups = inatgroups[-2]

# List range coverage rasters
spp = lapply(as.list(paste0("~/McGill University/Laura's Lab_Group - range-coverage/",taxagroups,"/")), list.files)
names(spp) = taxagroups
df2 = spp |> lapply(as.data.frame) |> bind_rows(.id = "Taxa") #|> na.omit()
colnames(df2)[2] = "species"

for(t in 1:length(taxagroups)){
  
  # Load species names 
  species = dplyr::filter(df2, Taxa == taxagroups[t]) |>
    select(species) |>
    as.matrix() |> as.vector()
  
  # import all rasters in this group
  rast.list = lapply(species, function(x) terra::rast(paste0("~/McGill University/Laura's Lab_Group - range-coverage/",taxagroups[t],"/", x))$change)
  
  # reproject the base grid for resampling so they all match
  canada = terra::project(canada, rast.list[[1]])
  rast.resamp = lapply(rast.list, resample, canada, method = "sum")
  for(i in 1:length(rast.resamp)){
    rast.resamp[[i]][rast.resamp[[i]]>1] = 1
  }
  # stack 'em
  rast.stack = terra::rast(rast.resamp)
  
  # sum stacked rasters (giving: new cell coverage for X number of species)
  rast.sum = sum(rast.stack, na.rm = T)
  
  # crop to canada
  rast.sum = crop(rast.sum, canada, mask = TRUE)
  
  # save the raster
  terra::writeRaster(rast.sum, paste0("outputs/range-coverage/summary-results/map_upgradedcells_sumspecies_",taxagroups[t],".tif"), overwrite = TRUE)
}

# read them back in to make plots and a full summary stack
all = lapply(taxagroups, function(x) terra::rast(paste0("outputs/range-coverage/summary-results/map_upgradedcells_sumspecies_",x,".tif")))
names(all) = taxagroups
all = lapply(all, project, "EPSG:3347")
canada = project(canada, "EPSG:3347")
all.resamp = lapply(all, resample, canada, method = "sum")
# stack 'em
all.stack = terra::rast(all.resamp)

# map

# ggplot() +
#   #geom_spatraster(data = canada, fill = "grey90",) +
#   geom_spatraster(data = all.stack) + 
#   facet_wrap(~lyr) +
#   colorspace::scale_fill_continuous_sequential("Batlow",
#                                                name = "",
#                                                #name = "Climate sampling balance",
#                                                na.value = "transparent", rev = F) +
#   hrbrthemes::theme_ipsum_rc(grid = FALSE, 
#                              axis = FALSE, axis_text_size = 1,
#                              ticks = FALSE,
#                              base_size = 14) +
#   theme(legend.position = "top", legend.key.width = unit(1.5,"cm"))


# sum stacked rasters (giving: new cell coverage for X number of species)
all.sum = sum(all.stack, na.rm = T)
all.sum = round(all.sum)


# plot to check it out
plot(all.sum)
mapview::mapview(all.sum)

# remove zeros
all.sum.nozeros = all.sum
all.sum.nozeros[all.sum.nozeros==0] <- NA
# save
terra::writeRaster(all.sum.nozeros, paste0("outputs/range-coverage/summary-results/map_upgradedcells_sumspecies_all.tif"), overwrite = TRUE)

# read in the raster
all.sum.nozeros = terra::rast(paste0("outputs/range-coverage/summary-results/map_upgradedcells_sumspecies_all.tif"))

# interactive map
mapview::mapviewOptions(basemaps = "OpenStreetMap",
               na.color = "transparent")

# plot it!
pal = viridis::turbo(5)
(m = mapview::mapview(all.sum.nozeros, 
        col.regions = pal,
        layer.name = "Species with new sightings", 
        na.color = "transparent", alpha.regions = .8) )
htmlwidgets::saveWidget(m@map, 
                        file = "outputs/range-coverage/interactive/allgroups_newsightings.html",
                        selfcontained = TRUE)
# published to RPubs
# https://rpubs.com/blitzthegap/newsightings


df = project(all.sum.nozeros, "EPSG:3347")
canada_poly = project(canada_poly, "EPSG:3347")


theme_set( hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                                      axis = FALSE, axis_text_size = 1,
                                      ticks = FALSE,
                                      base_size = 14) +
             theme(legend.position = "top", legend.key.width = unit(2.5,"cm")))

ggplot() +
  geom_sf(data = sf::st_as_sf(canada_poly), col = "grey90", fill = "grey90") +
  geom_spatraster(data = all.sum.nozeros, aes(fill = sum), interpolate = F) +
    colorspace::scale_fill_continuous_sequential("Batlow",
                                                 name = "Newly observed\nspecies per cell", trans = "sqrt",
                                                 na.value = "transparent", rev = F) 
ggsave("figures/gain_rangecoverage_alltaxa_map.png", width = 12.5, height = 7.9)
