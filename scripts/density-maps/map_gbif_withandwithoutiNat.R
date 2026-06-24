# Map of Canada GBIF with and without iNaturalist Canada

library(rgbif)
library(terra)
library(tidyterra)
library(sf)
library(rnaturalearth)
library(dplyr)
library(ggplot2)
library(gbifdb)
library(duckdb)

# set common theme for all maps
theme_set( hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                                       axis = FALSE, axis_text_size = 1,
                                       ticks = FALSE,
                                       base_size = 14) +
              theme(legend.position = "top", legend.key.width = unit(2.5,"cm")))

# Canada polygon
canada <- st_read("data/base-layers/canada-polygon/canada.outline.shp")
canada <- st_transform(canada, crs = "EPSG:3347")

# base grid
base100k = rast("data/base-layers/canada-basegrids/canada.base.100k.tiff")
base100k = project(base100k, crs(canada))

## GBIF - this takes a while but eventually works. Run once.

# # establish connection
# gbif <- gbif_remote()
# 
# # filter to canada
# gbifcan = gbif |>
#   filter(countrycode == "CA") |>
#   # keep only the useful columns for mapping
#   select(c(decimallatitude, decimallongitude, datasetkey, class)) |>
#   # round lat long to make some "cells"
#   mutate(latitude = round(decimallatitude,2),
#          longitude = round(decimallongitude,2)) |>
#   # count!
#   group_by(latitude, longitude, datasetkey, class) |>
#   summarise("n_obs" = n())
# #   
# # # collect it into R
# gbifcan = gbifcan |> collect()
# arrow::write_parquet(gbifcan, "data/heavy/BTG-data/gbif_all_latlong.parquet")

# read parquet file
gbifcan = arrow::open_dataset("data/heavy/BTG-data/gbif_all_latlong.parquet") |> 
  collect()

# map it
inatkey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7"

# all of GBIF (without birds) --------------------------------------------------

r_gbif_nobirds.pts = gbifcan |>
  filter(class != "Aves") |>
  group_by(longitude, latitude) |>
  summarise("n" = sum(n_obs, na.rm = T)) |>
  na.omit() |>
  vect(crs = "EPSG:4326",
       geom = c("longitude", "latitude"))
r_gbif_nobirds.pts = project(r_gbif_nobirds.pts, "EPSG:3347")
base100k = project(base100k, "EPSG:3347")

r_gbif_nobirds.rast = rasterize(r_gbif_nobirds.pts, 
                                base100k, 
                                field = "n",
                                fun = "sum", 
                                na.rm = T)
r_gbif_nobirds.rast = crop(r_gbif_nobirds.rast, canada, mask = TRUE)


## gbif without birds and without inat

n_inat = gbifcan |>
  filter(class != "Aves") |>
  filter(datasetkey != inatkey) |>
  summarise("n" = n())
n_gbif = gbifcan |>
  filter(class != "Aves") |>
  summarise("n" = n())


r_gbif_nobirds_noinat.pts = gbifcan |>
  filter(class != "Aves") |>
  filter(datasetkey != inatkey) |>
  group_by(longitude, latitude) |>
  summarise("n" = sum(n_obs, na.rm = T)) |>
  na.omit() |>
  vect(crs = "epsg:4326",
       geom = c("longitude", "latitude"))
r_gbif_nobirds_noinat.pts = project(r_gbif_nobirds_noinat.pts, crs(base100k))

r_gbif_nobirds_noinat.rast = rasterize(r_gbif_nobirds_noinat.pts, 
                                base100k, 
                                field = "n",
                                fun = "sum", 
                                na.rm = T)
r_gbif_nobirds_noinat.rast = crop(r_gbif_nobirds_noinat.rast, canada, mask = TRUE)

# only inat
r_gbif_nobirds_inat = gbifcan |>
  filter(class != "Aves") |>
  filter(datasetkey == inatkey) |>
  group_by(longitude, latitude) |>
  summarise("n" = sum(n_obs, na.rm = T)) |> 
  na.omit() |>
  vect(crs = "epsg:4326", 
       geom = c("longitude", "latitude"))
r_gbif_nobirds_inat = project(r_gbif_nobirds_inat, crs(base100k))

r_gbif_nobirds_inat.rast = rasterize(r_gbif_nobirds_inat, 
                                       base100k, 
                                       field = "n",
                                       fun = "sum", 
                                       na.rm = T)
r_gbif_nobirds_inat.rast = crop(r_gbif_nobirds_inat.rast, canada, mask = TRUE)

## make a map of the contributions of inat to each cell

# proportion of observations represented by inat
diffmap = r_gbif_nobirds_inat.rast/r_gbif_nobirds.rast
diffmap2 = resample(diffmap, base100k, method = "bilinear")
terra::writeRaster(diffmap2, "outputs/summaries/map_gbif_inaturalistcontribution.tif")