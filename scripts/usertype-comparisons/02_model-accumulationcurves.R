# Script to model the accumulation curve for users per cell 

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr)
library(duckdb)
library(terra)
library(iNEXT)

# set seed for randomisations
set.seed(22)

# set ggplot theme
theme_set(ggpubr::theme_pubr())

# load parquet
inat_pq <- arrow::open_dataset("data/heavy/iNat_non_sensitive_data_Jan2025.parquet")

# load user categories
users = readRDS("outputs/users/inat_users_categories.rds")

# Canada polygon
canada_poly = terra::vect("data/base-layers/canada-polygon/canada.outline.shp")

# load a modelled sample completeness map
sc = terra::rast("data/richness-maps/plant.richness.stacks.tif")
plot(sc)


# Subsample the users ----------------------------------------------------------

user.samp = users |>
  group_by(category) |>
  slice_sample(n = 10) # select 10 random users per group

# filter the parquet file to these users

user.inat <- inat_pq |> 
  # filter by species group:
  dplyr::filter(user_login %in% user.samp$user_login, 
                captive_cultivated == "false",
                iconic_taxon_name == "Plantae") |> 
  # select needed columns
  dplyr::select(c(latitude, longitude,
                  user_login,
                  scientific_name, 
                  iconic_taxon_name, 
                  observed_on_string, geometry)) |> 
  collect() 

# assign categories to each user
user.inat = left_join(user.inat, select(users, c(user_login, category)))

# Find some grid cells to train the rarefaction curve on -----------------------

# Make some rounded up grid cells from lat/long coordinates (lazy way)
user.inat = user.inat |>
  mutate(lat_round = round(latitude),
         long_round = round(longitude)) 
user.inat$cell = paste(user.inat$long_round, user.inat$lat_round, sep = "_")

# Count obs per cell
obs_cells = user.inat |>
  group_by(cell, user_login, scientific_name, category) |>
  summarise(total_obs = n())

# number of cells
cellids = obs_cells$cell |> unique()

# number of users per cell
cellObs = obs_cells |> group_by(cell) |> summarise("n_obs" = sum(total_obs))

# number of users per cell
cellUsers = obs_cells |> group_by(cell) |> summarise("n_users" = length(unique(user_login)))

# number of user groups per cell
cellCats = obs_cells |> group_by(cell) |> summarise("n_cat" = length(unique(category)))

# number of species per cell
cellSR = obs_cells |> group_by(cell) |> summarise("n_sp" = length(unique(scientific_name)))

# choose cells that have at least 4 categories of user types, for better comparison
cells.multicat = cellCats[which(cellCats$n_cat >= 4),]

# choose top 10 cells by number of observations
cells = cellObs |>
  filter(cell %in% cells.multicat$cell) |>
  slice_max(order_by = n_obs, n = 10)


## Rarefy observations in each cell --------------------------------------------

# function for a rarefaction curve to apply per cell, per user group

rarefy_groups = function(user.inat = user.inat, user_category, cell_index){
  
  # choose a cell
  temp = user.inat |> 
    dplyr::filter(cell == cells$cell[cell_index]) |> # cell index goes here
    group_by(category, user_login, scientific_name) |>
    summarise("total_obs" = n())
  
  # prepare data for rarefaction
  temp2 = temp |>
    dplyr::filter(category == user_category) |> # this is where the user category goes
    group_by(user_login) |>
    group_split() |>
    as.list() |>
    lapply(select, c(total_obs)) |>
    lapply(as.vector) |>
    lapply(unlist)
  
  # remove sites with less than 5 species
  index = unlist(lapply(temp2, length))
  temp2 = temp2[which(index>=5)]
  
  # rarefaction
  out_cell = iNEXT(temp2, q = 0, datatype = "abundance")
  
  return(out_cell)
  
}

ac_superuser = vector("list", 10)
for(i in 1:10){
  ac_superuser[[i]] = rarefy_groups("superuser", i)
  ggiNEXT(ac_superuser[[i]], type = 1,) + hrbrthemes::theme_ipsum_pub()
}


ac_expert = vector("list", 10)
for(i in 1:10){
  ac_expert[[i]] = rarefy_groups(user_category = "expert", cell_index =  i)
  ggiNEXT(ac_expert[[i]], type = 1,) + hrbrthemes::theme_ipsum_pub()
  ggsave("figures/usertype_comparisons/rarefaction-curves/by-user-group/expert_",i,".png")
}
