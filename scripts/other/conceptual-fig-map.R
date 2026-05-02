# Maps for a concept figure (sampling curves project)

library(tidyverse)
library(sf)
library(terra)

# load ecoregions and project
ecoreg = st_read("data/base-layers/ecoregions.shp")

ecoreg |> mapview::mapview()

reg = ecoreg |> filter(REGION_NAM == "Southern Laurentians")

mapview::mapview(reg, col.regions = "orange")


# load raster

basegrid = terra::rast("data/base-layers/canada-basegrids/canada.base.100k.tiff")
crs(basegrid)
basegrid = terra::project(basegrid, crs(reg))
basegrid = crop(basegrid, reg)
plot(basegrid)
lines(reg)


reg_raster = rasterize(reg, basegrid)
plot(reg_raster)

cellvals = values(reg_raster)

reg_raster[!is.na(cellvals)] = rnorm(length(cellvals[!is.na(cellvals)]), 10, 5)

plot(reg_raster)
mapview::mapview(reg_raster, na.color = "transparent")
