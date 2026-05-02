# summary stats

library(sf)
library(terra)
library(tidyterra)

# Canada polygon
canada <- st_read("data/base-layers/canada-polygon/canada.outline.shp")
canada <- st_transform(canada, crs = "EPSG:3347")

# base grid
base100k = rast("data/base-layers/canada-basegrids/canada.base.100k.tiff")
base100k = project(base100k, crs(canada))

# load gbif raster
dat = terra::rast("outputs/summaries/map_gbif_inaturalistcontribution.tif")
plot(dat)
plot(base100k)

# full canada area
area_can = cellSize(base100k, unit = "km")
can_area_km = sum(values(area_can), na.rm = T)

# area with data
area_dat = cellSize(dat, unit = "km")
area_dat[is.na(dat)] <- NA 
plot(area_dat)
dat_area_km = sum(values(area_dat), na.rm = T)

area_dat[dat<.5] <- NA
dat_area_km_overhalf = sum(values(area_dat), na.rm = T)

area_dat[dat<1] <- NA
dat_area_km_allinat = sum(values(area_dat), na.rm = T)

# proportion of area with data
dat_area_km/can_area_km
dat_area_km_overhalf/can_area_km
dat_area_km_allinat/can_area_km

# map
ggplot() +
  geom_sf(data = canada, fill = "grey95", col = "grey95") +
  geom_spatraster(data = dat) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                               na.value = "transparent",
                                               end = 1, begin = .1,
                                               name = "Proportion of GBIF data\ncontributed via iNaturalist") +
  theme_void() +
  theme(legend.position = "bottom")
ggsave("figures/map_inatcontribution_datadistribution_SI.png", width = 7.25, height = 7)
