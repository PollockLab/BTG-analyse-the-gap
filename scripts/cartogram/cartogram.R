# Cartogram of Canadian observation densities

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

# set up
stac_obj <- stac("https://io.biodiversite-quebec.ca/stac/")

# get iNaturalist heatmaps for Canada
it_obj <- stac_obj |>
  stac_search(collections = "inat_canada_heatmaps") |>
  post_request() |> items_fetch()
it_obj

# See layers in the object
df <- data.frame(description=character())
for (f in it_obj[['features']]){
  df <- rbind(df,data.frame(description=f$properties$description)) 
}
df$feature_number = 1:nrow(df)

# 25 is the number to change if you want another layer from the STAC
inat <- read_stars(paste0('/vsicurl/', it_obj[['features']][[25]]$assets[[1]]$href), 
                   proxy = FALSE) 
inat = terra::rast(inat)
# aggregate into bigger cells for faster computation
inat = aggregate(x = inat, factor = 100, fun = "sum", na.rm = TRUE)

# read the Canada polygon
canada = sf::read_sf("data/base-layers/canada-polygon/canada.outline.shp")

# crop the raster to canada (to cut off the ocean, USA, etc.)
canada = st_transform(canada, crs = terra::crs(inat))
canada = vect(canada)

# crop density map to Canada and project
inat = terra::crop(inat, canada, mask = TRUE)
inat = project(inat, "epsg:4326")

# load ecoregions and project
ecoreg = st_read("~/Desktop/Ecoregions/ecoregions.shp")
ecoreg = st_transform(ecoreg, crs = "epsg:4326")
st_write(ecoreg, "data/base-layers/ecoregions.shp", append = FALSE)
ecoreg = vect("data/base-layers/ecoregions.shp")

# sum number of obs per polygon
sum_by_region = terra::zonal(x = inat, z = ecoreg, 
                             fun = "sum", na.rm = TRUE, 
                             as.polygons = TRUE, touches = TRUE)

# convert to sf and prepare for plotting
sum_df = st_as_sf(sum_by_region)
sum_df$n_obs = sum_df$All_density_inat_1km.tif


# make a cartogram
sum_df = st_transform(sum_df, crs(canada))
df_cartogram_raw <- cartogram_cont(sum_df, "n_obs", itermax = 50)

# plot
ggplot() +
  geom_sf(data = st_as_sf(ecoreg), 
          aes(fill = n_obs),
          linewidth = .1) +
  labs(fill = "Observations") +
  colorspace::scale_fill_continuous_sequential("BuGn") +
  theme_void() +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(family = "Open Sans"),
        legend.title = element_text(family = "Open Sans", face = "bold"))
ggsave("figures/cartograms/map_ecoregions_original.png", width = 8, height = 6)

ggplot(df_cartogram_raw) +
  geom_sf(aes(fill = n_obs), linewidth = .1) +
  colorspace::scale_fill_continuous_sequential("BuGn") +
  labs(fill = "Observations") +
  theme_void() +
  theme(legend.position = "bottom", 
        legend.key.width = unit(2, "cm"),
        legend.text = element_text(family = "Open Sans"),
        legend.title = element_text(family = "Open Sans", face = "bold"))
ggsave("figures/cartograms/map_ecoregions_nobs.png", width = 8, height = 6)