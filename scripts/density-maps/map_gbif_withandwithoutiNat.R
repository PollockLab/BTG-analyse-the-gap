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

# Canada polygon
canada <- st_read("data/base-layers/canada-polygon/canada.outline.shp")
canada <- st_transform(canada, crs = "EPSG:3347")

# base grid
base100k = rast("data/base-layers/canada-basegrids/canada.base.100k.tiff")
base100k = project(base100k, crs(canada))

# establish connection
gbif <- gbif_remote()

# filter to canada
gbifcan = gbif |>
  filter(countrycode == "CA") |>
  # keep only the useful columns for mapping
  select(c(decimallatitude, decimallongitude, datasetkey, class)) |>
  # round lat long to make some "cells"
  mutate(latitude = round(decimallatitude,1),
         longitude = round(decimallongitude,1)) |>
  # count!
  group_by(latitude, longitude, datasetkey, class) |>
  summarise("n_obs" = n())
  
# collect it into R
gbifcan = gbifcan |> collect()
arrow::write_parquet(gbifcan, "data/heavy/BTG-data/gbif_all_latlong.parquet")

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
  vect(crs = "epsg:4326",
       geom = c("longitude", "latitude"))
r_gbif_nobirds.pts = project(r_gbif_nobirds.pts, crs(base100k))
r_gbif_nobirds.rast = rasterize(r_gbif_nobirds.pts, 
                                base100k, 
                                field = "n",
                                fun = "sum", 
                                na.rm = T)
r_gbif_nobirds.rast = crop(r_gbif_nobirds.rast, canada, mask = TRUE)


## gbif without birds and without inat

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
plot(diffmap)

# aggregate a little bit to make it cleaner
diffmap_lg = aggregate(diffmap, factor = 1, fun = "mean", na.rm = T)

# summary table - base100k
cellvals = values(diffmap)
totalcells =  nrow(cellvals)
empties = length(which(is.na(cellvals)))
df = data.frame(
  "summary" = c("inat100p",
                "inat90p",
                "inat75p",
                "inat50p",
                "inatunder50p"),
  "n_cells" = c(
    length(which(cellvals==1)),
    length(which(cellvals>=.9)),
    length(which(cellvals>=.75)),
    length(which(cellvals>=.5)),
    length(which(cellvals<.5))
  )
)
df$perc = 100*df$n_cells/(totalcells-empties)
write.table(df, "outputs/map_gbif_inaturalistcontribution.csv")

# Plotting
ggplot() +
  geom_sf(data = canada, fill = "grey10", color = "grey10") +
  geom_spatraster(data = 100*diffmap_lg,
                  interpolate = FALSE) +
  scale_fill_viridis_c(option = "turbo",
                       name = "iNaturalist\ncontribution (%)",
                       end = .8, begin = .2,
                       limits = c(0,100),
                       na.value = "transparent") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2.5,"cm"))
ggsave("figures/map_gbif_inaturalistcontribution.png", width = 12.5, height = 7.9)

ggplot() +
  geom_sf(data = canada, fill = "grey10", color = "grey10") +
  geom_spatraster(data = 100*filter(diffmap_lg, sum  == 1),
                  interpolate = FALSE) +
  scale_fill_viridis_c(option = "turbo",
                       name = "iNaturalist\ncontribution (%)",
                       end = .8, begin = .2,
                       limits = c(0,100),
                       na.value = "transparent") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2.5,"cm"))
ggsave("figures/map_gbif_only100pinaturalist.png", width = 12.5, height = 7.9)


ggplot() +
  geom_sf(data = canada, fill = "grey10", color = "grey10") +
  geom_spatraster(data = 100*filter(diffmap_lg, sum  > .9),
                  interpolate = FALSE) +
  scale_fill_viridis_c(option = "turbo",
                       name = "iNaturalist\ncontribution (%)",
                       end = .8, begin = .2,
                       limits = c(0,100),
                       na.value = "transparent") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2.5,"cm"))
ggsave("figures/map_gbif_morethan90pinaturalist.png", width = 12.5, height = 7.9)

ggplot() +
  geom_sf(data = canada, fill = "grey10", color = "grey10") +
  geom_spatraster(data = 100*filter(diffmap_lg, sum  >= .75),
                  interpolate = FALSE) +
  scale_fill_viridis_c(option = "turbo",
                       name = "iNaturalist\ncontribution (%)",
                       end = .8, begin = .2,
                       limits = c(0,100),
                       na.value = "transparent") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2.5,"cm"))
ggsave("figures/map_gbif_morethan75pinaturalist.png", width = 12.5, height = 7.9)


ggplot() +
  geom_sf(data = canada, fill = "grey10", color = "grey10") +
  geom_spatraster(data = 100*filter(diffmap_lg, sum  >= .5),
                  interpolate = FALSE) +
  scale_fill_viridis_c(option = "turbo",
                       name = "iNaturalist\ncontribution (%)",
                       end = .8, begin = .2,
                       limits = c(0,100),
                       na.value = "transparent") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2.5,"cm"))
ggsave("figures/map_gbif_morethan50pinaturalist.png", width = 12.5, height = 7.9)
