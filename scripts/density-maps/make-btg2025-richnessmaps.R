# Script to make annual species richness maps 

# load libraries
library(tidyverse)
library(sf)
library(cartogram)
library(terra)
library(tidyterra)
library(ggplot2)
library(arrow)

# set theme for all plots
theme_set(hrbrthemes::theme_ipsum_rc())

# Input data -------------------------------------------------------------------

# load base grid
base1k = terra::rast("data/base-layers/canada-basegrids/canada.base.1k.tiff")
# read the Canada polygon
canada = sf::read_sf("data/base-layers/canada-polygon/canada.outline.shp")

# load parquet file
df = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_JUN1OCT1.parquet")

# make some maps ---------------------------------------------------------------

# make observations layer
temp = df |>
  select(c(longitude, latitude, iconic_taxon_name, scientific_name)) |>
  collect() |>
  terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
temp = project(temp, crs(base1k))

# get list of species
sciname = temp$scientific_name |> unique()
# check for names with two strings to narrow down to species only
gbifname = taxize::get_gbifid(sciname, ask = FALSE)
ranks = taxize::tax_rank(sciname, db = "gbif", ask = FALSE)

# get list of taxa
taxa = temp$iconic_taxon_name |> unique()

# get cells
temp$cell = extract(base1k, temp, cell=T)[,3]

# thin to one species occ per cell (i.e. convert density to just presence/absence)
temp.thin <- temp[-which(duplicated(paste(temp$scientific_name,temp$cell)))]
richmap <- rasterize(temp.thin, base1k, fun="length")
# save
terra::writeRaster(richmap, "outputs/richnessmaps/richnessmap_JUN1OCT12025_alltaxa.tif", overwrite = TRUE)

ggplot() +
  geom_sf(data = canada, col = "grey90") +
  geom_spatraster(data = richmap, aes(fill = V1_length)) +
  colorspace::scale_fill_continuous_sequential("Batlow", 
                                               trans = "log10", 
                                               na.value = "transparent", 
                                               rev = F) +
  labs(fill = "Taxa") +
  theme_void() +
  theme(legend.position = "top", 
        legend.key.width = unit(2, "cm"),
        text = element_text(family = "Roboto Condensed"))
ggsave("figures/richness_map_JUN1OCT12025_DEC1DATA.png", width = 6.41, height = 5.77)


# get list of taxa
taxa.df = df |> group_by(iconic_taxon_name) |> summarise("n" = n()) |> collect()
taxa = na.omit(taxa.df$iconic_taxon_name)

# map per taxa group
richmap_taxa = list()
for(t in 1:length(taxa)){
  
  # make observations layer
  temp = df |>
    dplyr::filter(iconic_taxon_name == taxa[t]) |>
    select(c(longitude, latitude, iconic_taxon_name, scientific_name)) |>
    collect() |>
    terra::vect(geom = c("longitude", "latitude"), crs = "epsg:4326")
  # project the observations layer
  temp = project(temp, crs(base1k))
  
    # get cells
    temp$cell = extract(base1k, temp, cell=T)[,3]
    
    # thin to one species occ per cell (i.e. convert density to just presence/absence)
    temp.thin <- temp[-which(duplicated(paste(temp$scientific_name, temp$cell)))]
    map_temp[[t]] = rasterize(temp.thin, base1k, fun="length")
  }
  names(map_temp) = taxa
  terra::writeRaster(map_temp, paste0("outputs/richnessmaps/richnessmap_yearly_",taxa[t],".tif"), overwrite = TRUE)
