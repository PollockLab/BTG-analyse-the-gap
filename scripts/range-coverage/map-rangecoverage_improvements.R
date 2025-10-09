# Script to map the cells where range coverage improved for X species

# load libraries
library(ggplot2)
library(dplyr)
library(terra)
library(gpkg)
library(plotly)
library(hrbrthemes)

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# load canada base grid
canada = terra::rast("~/Documents/GitHub/ciee/blitz-the-gap/00_rawdata/base-layers/canada.base.5k.tiff")
canada = terra::aggregate(canada, factor = 2) # 10 km^2 cells

# load canada polygon
canada_poly = terra::vect("~/Documents/GitHub/ciee/blitz-the-gap/00_rawdata/base-layers/canada-polygon/canada.outline.shp")

# list of taxa groups
taxagroups = list.files("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/")
inatgroups = c("Amphibia", "Aves", "Insecta", "Mammalia", "Plantae", "Reptilia")

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
  rast.list = lapply(species, function(x) terra::rast(paste0("~/McGill University/Laura's Lab_Group - range-coverage/upgraded_cells/",taxagroups[t],"/min_1obs/", x))$upgraded)
  
  # reproject the base grid for resampling so they all match
  canada = terra::project(canada, rast.list[[1]])
  rast.resamp = lapply(rast.list, resample, canada, method = "sum")
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
temp = terra::rast(paste0("outputs/range-coverage/summary-results/map_upgradedcells_sumspecies_",taxagroups[t],".tif"))
plot(rast.sum)
mapview::mapview(temp)
