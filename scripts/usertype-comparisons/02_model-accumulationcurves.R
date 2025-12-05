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
  slice_sample(n = 100) # select 10 random users per group

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
cells.multicat = cellCats[which(cellCats$n_cat >= 5),]

# choose top 10 cells by number of observations
cells = cellObs |>
  filter(cell %in% cells.multicat$cell) |>
  slice_max(order_by = n_obs, n = 10)


## Rarefy observations in each cell --------------------------------------------

# function for a rarefaction curve to apply per cell, per user group -----------

rarefy_groups = function(data = user.inat, data_cells = cells, user_category, cell_index){
  
  # choose a cell
  temp = data |> 
    dplyr::filter(cell == data_cells$cell[cell_index]) |> # cell index goes here
    group_by(category, user_login, scientific_name) |>
    summarise("total_obs" = n())
  
  # prepare data for rarefaction
  temp2 = temp |>
    dplyr::filter(category == user_category) |> # this is where the user category goes
    group_by(user_login) |>
    group_split() |>
    as.list() 
  names(temp2) = lapply(temp2, function(x) return(unique(x[,colnames(x) == "user_login"]))) |> unlist()
  
  temp2 = temp2 |>
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

# function to run this for each user category ----------------------------------

run_rarefaction = function(CATEGORY){
  ac = vector("list", length(cells$cell))
  names(ac) = cells$cell
  for(i in 1:10){
    tryCatch({
      ac[[i]] = rarefy_groups(user_category = CATEGORY, 
                                         cell_index =  i)
      ggiNEXT(ac[[i]], type = 1,) + 
        hrbrthemes::theme_ipsum_pub()
      ggsave(paste0("figures/usertype_comparisons/rarefaction-curves/by-user-group/",CATEGORY, "_", i,".png"),
             width = 6, height = 6)
    }, error = function(e) {ac[[i]] = NA} )
  }
  return(ac)
}

# run!!
ac_superuser = run_rarefaction("superuser")
ac_expert = run_rarefaction("expert")
ac_enthusiast = run_rarefaction("enthusiast")
ac_dabbler = run_rarefaction("dabbler")
ac_casual = run_rarefaction("casual")

ac = list(ac_superuser, ac_expert, ac_enthusiast, ac_dabbler, ac_casual)
saveRDS(ac, "outputs/users/usergroups_rarefactions.rds")

## Fit a model to these curves -------------------------------------------------

# function to extract asymptotes of the curves -----

cell_Asy = function(ac_usergroup){
  
  asy = vector("list", length(ac_usergroup))
  names(asy) = names(ac_usergroup)

  for(i in 1:length(ac_usergroup)){
    
    tryCatch({
      
    asy[[i]] = ac_usergroup[[i]]$AsyEst |> 
      filter(Diversity == "Species richness") |>
      arrange(-Observed)
    }, error = function(e) {asy[[i]] = NA})
  }
    
  asy = asy |> bind_rows(.id = "cell")
  
  return(asy)
}

# apply to the iNEXT results
asy = lapply(ac, cell_Asy)
names(asy) = c("superuser", "expert", "enthusiast", "dabbler", "casual")
asy = bind_rows(asy, .id = "category")


# Model and predict the rarefaction curves -------------------------------------

# function to extract sample coverage from the iNEXT results

cell_SC = function(ac_usergroup){
  
  sc = vector("list", length(ac_usergroup))
  names(sc) = names(ac_usergroup)
  
  for(i in 1:length(ac_usergroup)){
    
    tryCatch({
      
      sc[[i]] = ac_usergroup[[i]]$iNextEst$coverage_based |> 
        filter(SC >= 0.5) |>
        group_by(Assemblage) |>
        slice_head() |>
        mutate("coverage" = "50%")
    }, error = function(e) {sc[[i]] = NA})
  }
  
  sc = sc |> bind_rows(.id = "cell")
  
  return(sc)
}

# apply to the results
sc = lapply(ac, cell_SC)
names(sc) = c("superuser", "expert", "enthusiast", "dabbler", "casual")
sc = bind_rows(sc, .id = "category")

# plot the lines!
ggplot(data = sc) +
  geom_point(aes(y = qD, x = m, col = category, fill = category)) +
  geom_smooth(aes(y = qD, x = m, col = category, fill = category), 
              method = "loess", se = T) +
  labs(y = "Estimated total SR", 
       x = "Number of observations")




# Fit a GAM to the accumulation curves ========================================

## this is where i need to edit 
# model the curves
m50 = mgcv::gam(m ~ s(qD), data = SC_obsneeded_cov50)

rz_m = list(m50, m75, m90, m100)

# function to predict the number of observations needed to find the expected SR
# given the rarefaction curves modelled above
rz_effort = function(SR){
  x = lapply(rz_m, predict, newdata = data.frame("qD" = SR)) |>
    rbind() |>
    as.data.frame() 
  colnames(x) = c("p50", "p75", "p90", "p100")
  return(x)
}