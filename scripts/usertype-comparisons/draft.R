# Script to analyse the contributions of different users to Blitz the Gap

# libraries
library(rinat)
library(rgbif)
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(taxize)
library(rinat)


# List iNaturalist projects linked to Blitz the Gap ----------------------------

btg.projects = c(
  "blitz-the-gap",
  "blitz-the-gap-revisiting-the-past",
  "blitz-the-gap-missing-canadian-animals",
  "blitz-the-gap-missing-canadian-fungi",
  "blitz-the-gap-missing-canadian-plants",
  "blitz-the-gap-missing-canadian-species-insects",
  "blitz-the-gap-made-in-canada",
  "blitz-the-gap-seaweeds-of-canada",
  "blitz-the-gap-conservation-priorities-in-canada-maybas",
  "blitz-the-gap-too-hot-to-handle",
  "blitz-the-gap-canada-s-most-wanted",
  "blitz-the-gap-hey-birders-don-t-look-up",
  "blitz-the-gap-make-a-splash",
  "blitz-the-gap-trailblazers",
  "blitz-the-gap-more-than-the-monarch",
  "blitz-the-gap-closing-the-climate-gap",
  "blitz-the-gap-the-other-99",
  "blitz-the-gap-getting-even"
)

btg.bioblitzes = c(
  "expedition-fiord-arctic-bioblitz",
  "cypress-hills-bioblitz-2025",
  "bioblitz-the-gap-cba-2025",
  "bioblitz-the-gap-gault-nature-reserve",
  "bioblitz-the-gap-station-de-biologie-des-laurentides",
  "irving-nature-park-ocean-week-blitz-the-gap",
  "fundy-trail-ocean-week-fundy-blitz-the-gap",
  "foret-ouareau-bioblitz-2025",
  "bioblitz-the-gap-csee-ile-du-marais",
  "blitz-the-gap-quebec-soil-fauna",
  "bioblitz-the-gap-csee-johnville-bog-forest-park",
  "bioblitz-the-gap-csee-parc-national-du-mont-megantic",
  "ram-mountain-biodiversity-blitz-2025",
  "ubc-farm-bioblitz-blitz-the-gap"
)


# get ID of these projects
info.btg = lapply(btg.projects, get_inat_obs_project, type = "info", raw = FALSE)
info.bioblitz = lapply(btg.bioblitzes, get_inat_obs_project, type = "info", raw = FALSE)
info.qcbs = get_inat_obs_project("blitz-the-gap-qcbs-champions", type = "info", raw = FALSE)

# get top observers of the project
users.btg = get_inat_user_stats(project = info.btg[[1]]$id, 
                                date_range = c("2025-06-01", paste0(Sys.Date())))
saveRDS(users.btg, "outputs/users/btg_users_allinfo.rds")                                                        
users.btg = readRDS("outputs/users/btg_users_allinfo.rds")                                                        
View(btg_users$most_species)    
View(btg_users$most_observations)                            
btg_users$most_species |> nrow()
btg_users$most_species$count |> hist()


# get top observers of the bioblitzes
users.bioblitzes = lapply(info.bioblitz, function(x) get_inat_user_stats(project = x$id, 
                                date_range = c("2025-06-01", paste0(Sys.Date()))))
names(users.bioblitzes) = btg.bioblitzes

ids = users.bioblitzes$`expedition-fiord-arctic-bioblitz`$most_observations$user$id
logins = users.bioblitzes$`expedition-fiord-arctic-bioblitz`$most_observations$user$login

taxa_id = get_inat_taxon_stats(uid = "katherinehebert")

date.ranges = list(
  c("2025-06-01", "2025-09-01"),
  c("2024-06-01", "2024-09-01"),
  c("2023-06-01", "2023-09-01"),
  c("2022-06-01", "2022-09-01"),
  c("2021-06-01", "2021-09-01")
)

temp = get_inat_taxon_stats(date_range = date.ranges[[1]], 
                            uid = ids[[1]])
temp2 = get_inat_ (date_range = date.ranges[[1]], 
                            uid = ids[[1]])
user.taxons = vector("list", length = length(ids))
names(user.taxons) = logins

steps = seq(1, length(ids), by = 25)
for(i in 21:31){
  
  # make a data frame
  user.taxons[[i]] = data.frame(
    "year" = c(2021:2025),
    "n_sp" = NA
  )
  
  # get yearly count of species
  yearly = list()
  for(y in 1:length(date.ranges)){
    user.taxons[[i]][y,2] = get_inat_taxon_stats(date_range = date.ranges[[y]], 
                                       uid = ids[[i]])$total
  }
}
user.taxons[1:35]
## load inat parquet file


# get list of users in:
# all BTG
# qcbs (incentivized)
# inat parquet

# Parquet:
# identify top 100 users (superusers) - or quantile it?

# BTG newbies:
# check for new users this year. more or less than usual?
