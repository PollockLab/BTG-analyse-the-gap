# Sampling priority index derived from applying a general curve from 
# filled cells in an ecoregion 

library(tidyverse)
library(sf)
library(terra)
library(iNEXT)
library(arrow)

theme_set(ggpubr::theme_pubr())

# load inat data
inat_pq <- arrow::open_dataset("~/McGill University/Laura's Lab_Group - BioBlitz/data/raw/biodiversity-data/inat-canada/iNat_non_sensitive_data_Jan2025.parquet")

# load ecoregions and project
ecoreg = st_read("data/base-layers/ecoregions.shp")

# load base grid
base5k = terra::rast("data/base-layers/canada-basegrids/canada.base.5k.tiff")


## Select an ecoregion ---------------------------------------------------------

# pick a region
reg = ecoreg |> filter(REGION_NAM == "Central Laurentians")
plot(reg)

# get the bounding box
bbox = st_bbox(reg)

# convert to terra and adjust projections
reg = terra::vect(reg)
base5k = terra::project(base5k, crs(reg)) |> aggregate(factor = 10)

# rasterize
reg_grid = rasterize(reg, base5k, field = "ECOZONE") |> trim()

## Subset iNat to the ecoregion ------------------------------------------------

reg_obs <- inat_pq |> 
  dplyr::filter(captive_cultivated == FALSE,
                between(latitude, bbox[2], bbox[4]),
                between(longitude, bbox[1], bbox[3])) |> 
  select(iconic_taxon_name, scientific_name,
         latitude, longitude,
         quality_grade, user_login, observed_on_string) |>
  # load the query into our R session
  collect()

# spatialise
obs = terra::vect(reg_obs, geom = c("longitude", "latitude"), crs = "epsg:4326")
obs = project(obs, crs = crs(reg_grid))

# rasterize
obs_grid = rasterize(obs, reg_grid, fun = "length")

# remove points outside the region
obs_grid[is.na(reg_grid)] <- NA

# plot
plot(obs_grid)


## Select well-sampled cells ===================================================

hist(obs_grid)
quantile(values(obs_grid), na.rm = T)

# threshold cells (to define more thoughtfully later)
sampled_well = obs_grid
sampled_well[obs_grid < 4000] <- NA
plot(sampled_well)

# extract observations in these cells (NA = not selected as well sampled)
obs_sel = terra::extract(sampled_well, obs) 

# convert to sf
obs_sf = st_as_sf(obs)
obs_sf$selected = obs_sel$V1_length

# get unique cell values
cells = unique(obs_sf$selected)

# Rarefy in one cell -----------------------------------------------------------

# subset data
df = obs_sf |> 
  filter(selected == na.omit(cells)) |>
  mutate(group = as.character(selected)) |>
  select(c(group, scientific_name)) |>
  group_by(group, scientific_name) |>
  summarise(n = n()) |>
  sf::st_drop_geometry()

# prep data for iNEXT
groups = unique(df$group)
# there are nicer ways to do this but whatever! for now :)
occ = list()
for(i in 1:length(groups)){
  df2 = df |> filter(group == groups[i])
  occ[[i]] = df2$n |> sort(decreasing = TRUE)
}
names(occ) = groups
# rarefaction
out_cell = iNEXT(occ, q = 0, datatype = "abundance")
plot(out_cell) 
ggiNEXT(out_cell)
ggiNEXT(out_cell, type = 2)
