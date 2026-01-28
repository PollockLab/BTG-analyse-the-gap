# Cartogram of Canadian observation densities in 2025

library(tidyverse)
library(sf)
library(cartogram)
library(gdalcubes)
library(rstac)
library(knitr)
library(stars)
library(terra)
library(ggplot2)
library(mapview)
library(arrow)

# load base grid
base100k = terra::rast("data/base-layers/canada-basegrids/canada.base.100k.tiff") 
# read the Canada polygon
canada = sf::read_sf("data/base-layers/canada-polygon/canada.outline.shp")
# load ecoregions
ecoreg = vect("data/base-layers/ecoregions.shp")

## user convex hull geographically - bigger than last years?
## new users' behaviour??

# load parquet file
df_total = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025.parquet")
df_2025 = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_subset2025.parquet")

# remove 2025 data from the total dataset
df_no2025 = df_total |> dplyr::filter(!grepl("2025", observed_on_string))
df_2024 = df_total |> dplyr::filter(grepl("2024", observed_on_string))

# load 2025 data
obs = df_2025 |>
  select(c(latitude, longitude, scientific_name, iconic_taxon_name)) |>
  collect()
obs.2024 = df_2024 |>
  select(c(latitude, longitude, scientific_name, iconic_taxon_name)) |>
  collect()

# spatialise the observations
inat = vect(obs, geom = c("longitude", "latitude"), crs = "epsg:4326")
inat = project(inat, crs(base100k))
inat_rast = rasterize(inat, base100k, fun = "length")
# spatialise the observations (2024)
inat.2024 = vect(obs.2024, geom = c("longitude", "latitude"), crs = "epsg:4326")
inat.2024 = project(inat.2024, crs(base100k))
inat.2024 = crop(inat.2024, base100k)
inat_rast.2024 = terra::rasterize(inat.2024, base100k, fun = "length")

# project to match ecoregions
inat_rast = project(inat_rast, crs(ecoreg))
inat_rast.2024 = project(inat_rast.2024, crs(ecoreg))

# sum number of obs per polygon
sum_by_region = terra::zonal(x = inat_rast, z = ecoreg, 
                             fun = "sum", na.rm = TRUE, 
                             as.polygons = TRUE, touches = TRUE)
# sum number of obs per polygon
sum_by_region.2024 = terra::zonal(x = inat_rast.2024, z = ecoreg, 
                             fun = "sum", na.rm = TRUE, 
                             as.polygons = TRUE, touches = TRUE)
plot(sqrt(sum_by_region.2024$V1_length), sqrt(sum_by_region$V1_length))

# convert to sf and prepare for plotting
sum_df = st_as_sf(sum_by_region)
sum_df$n_obs = sum_df$V1_length


# make a cartogram
sum_df = st_transform(sum_df, crs(canada))
saveRDS(sum_df, "outputs/cartogram/sum_df_2025.rds")

# calculate cartogram 
df_cartogram_raw <- cartogram_cont(sum_df, "n_obs", itermax = 50)
sum_df$n_obs_sqrt = sqrt(sum_df$n_obs)
df_cartogram_sqrt <- cartogram_cont(sum_df, "n_obs_sqrt", itermax = 50)

# plot
ggplot() +
  geom_sf(data = st_as_sf(sum_df), 
          aes(fill = n_obs),
          linewidth = .1) +
  labs(fill = "Observations") +
  colorspace::scale_fill_continuous_sequential("BuGn") +
  theme_void() +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(family = "Open Sans"),
        legend.title = element_text(family = "Open Sans", face = "bold"))
ggsave("figures/cartograms/map_ecoregions_original_2025.png", width = 8, height = 6)

ggplot() +
  geom_sf(data = st_as_sf(sum_df), 
          aes(fill = n_obs_sqrt),
          linewidth = .1) +
  labs(fill = "Observations\n(sqrt)") +
  colorspace::scale_fill_continuous_sequential("BuGn") +
  theme_void() +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(family = "Open Sans"),
        legend.title = element_text(family = "Open Sans", face = "bold"))
ggsave("figures/cartograms/map_ecoregions_originalsqrt_2025.png", width = 8, height = 6)


ggplot(df_cartogram_raw) +
  geom_sf(aes(fill = sqrt(n_obs)), linewidth = .1) +
  colorspace::scale_fill_continuous_sequential("BuGn") +
  labs(fill = "Observations\n(sqrt)") +
  theme_void() +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(family = "Open Sans"),
        legend.title = element_text(family = "Open Sans", face = "bold"))
ggsave("figures/cartograms/map_ecoregions_nobs_2025.png", width = 8, height = 6)


ggplot(df_cartogram_sqrt) +
  geom_sf(aes(fill = n_obs_sqrt), linewidth = .1) +
  colorspace::scale_fill_continuous_sequential("BuGn") +
  labs(fill = "Observations") +
  theme_void() +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(family = "Open Sans"),
        legend.title = element_text(family = "Open Sans", face = "bold"))
ggsave("figures/cartograms/map_ecoregions_nobssqrt_2025.png", width = 8, height = 6)
