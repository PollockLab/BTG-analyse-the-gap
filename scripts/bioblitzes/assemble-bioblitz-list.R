# Script to retrieve Bioblitzes

library(rinat)
library(tidyverse)

# info_wrap = function to reformat multi-project information lists -------------

# category: label for the bioblitz grouping
# multi: are there multiple projects in the ID list ?
info_wrap = function(x, multi = FALSE, category = NA){
  
  if(multi == TRUE){
  res = lapply(x, 
                   function(y) {unlist(y) |> 
                       as.matrix() |> t() |>
                       as.data.frame()}) |> bind_rows()
  } else {
    res = x |> unlist() |> as.matrix() |> t() |> as.data.frame()
  }
  res$category = category
  return(res)
}


## CAMPUS NATURE CHALLENGE -----------------------------------------------------

#Umbrella: https://www.inaturalist.org/projects/campus-nature-challenge-2025-defi-nature-campus
cnc_urls = c(
  "mount-allison-university-campus-nature-challenge-2025",
  "collaboration-heritage-college-cegep-de-l-outaouais-campus-nature-challenge-2025",
  "bioblitz-du-cegep-de-l-outaouais-septembre-2025",
  "fanshawe-college-campus-nature-challenge-2025",
  "vanier-college-campus-nature-challenge-2025",
  "mcgill-university-campus-nature-challenge-2025",
  "universite-laval-defi-nature-campus-2025",
  "2025-campus-nature-challenge-niagara-college",
  "inrs-defi-nature-campus-2025",
  "uqam-defi-nature-campus-2025",
  "bioblitz-armand-frappier-inrs-canopee-pour-le-defi-nature-campus-et-grand-bioblitz-canadien-2025",
  "uqo-defi-nature-campus-2025",
  "bioblitz-avec-le-cegep-garneau-2025",
  "universite-de-sherbrooke-defi-nature-campus-2025",
  "uqar-defi-nature-campus-2025",
  "bioblitz-inrs-teluq-uq-defi-nature-campus-et-grand-bioblitz-canadien-2025",
  "ets-ecole-de-technologie-superieure-defi-nature-campus-2025",
  "cegep-du-vieux-montreal-defi-nature-campus-2025",
  "rrc-polytech-campus-nature-challenge-2025",
  "st-thomas-university-nb-campus-nature-challenge-2025"
)

proj1 = lapply(cnc_urls, 
               function(x) rinat::get_inat_obs_project(grpid = x, type = "info"))
proj1b = info_wrap(proj1, multi = TRUE, category = "Campus Nature Challenge")


## CWF & BTG BIOBLITZES --------------------------------------------------------------

# ID list from James
project_ids = c(248464,227409,240788,241287,243158,243161,244779,241300,248412,241298,241296,241314,249523,246017,255676,256903,256906,256909,256899,257184,257310,257307)
proj2 = lapply(project_ids, function(x) rinat::get_inat_obs_project(grpid = x, type = "info"))
proj2b = info_wrap(proj2, multi = TRUE)
proj2b$category[1:14] <- "Blitz the Gap: Local bioblitzes"
proj2b$category[15:nrow(proj2b)] <- "Great Canadian Bioblitz"
proj2b = proj2b[-which(proj2b$title == "BioBlitz the Gap: Heartbeet Community Farm"),] # this bioblitz did not happen


## CWF GRANTEES ----------------------------------------------------------------

cwf.funded = readxl::read_xlsx("data/heavy/bioblitzes/project-info/GCB funded blitz obs - James.xlsx",sheet = 1)
cwf.funded2 = readxl::read_xlsx("data/heavy/bioblitzes/project-info/GCB funded blitz obs - James.xlsx",sheet = 2)
cwf.experts = readxl::read_xlsx("data/heavy/bioblitzes/project-info/GCB expert obs - from James.xlsx")

# get unique urls
cwf.urls = c(unique(cwf.funded$`link to GCB obs`), 
             unique(cwf.funded2$`iNat proj`[grep("www", cwf.funded2$`iNat proj`)]),
             unique(cwf.experts$`link to GCB obs`)) |> 
  unique()
# break up a multi-url string
index = grep("\r\n", cwf.urls)
cwf.url.split = cwf.urls[index] |> strsplit(split = "\r\n") |> unlist()
cwf.urls = cwf.urls[-index]
cwf.urls = c(cwf.urls, cwf.url.split)
# missing urls (manually by searching the bioblitz name)
missing.urls = c(
  "duc-acadia-university-campus-club-great-canadian-bioblitz",
  "unb-wetlands-conservation-society-great-canadian-bioblitz",
  "univert-laval-le-grand-bioblitz-canadien",
  "wetlanders-duc-campus-club-great-canadian-bioblitz-2025",
  "university-of-guelph-great-canadian-bioblitz",
  "mst-2025-fall-walks-and-forays"
)
cwf.urls = c(cwf.urls, missing.urls)
cwf.urls = gsub("https://www.inaturalist.org/projects/", "", cwf.urls)
cwf.urls = gsub("https://www.inaturalist.ca/projects/", "", cwf.urls)
cwf.urls = gsub("\\?tab=observers", "", cwf.urls)
# remove bug battle umbrella project
cwf.urls = cwf.urls[-which(cwf.urls == "bug-battle-2025-kings-and-annapolis-results")]

# set aside the observer observations urls
cwf.urls.observers = cwf.urls[grep("observations", cwf.urls)]
cwf.urls = cwf.urls[-grep("observations", cwf.urls)]
cwf.urls = unique(cwf.urls)

# extract project information
proj3 = lapply(cwf.urls, 
               function(x) rinat::get_inat_obs_project(grpid = x, type = "info"))
proj3b = info_wrap(proj3, multi = TRUE, "Great Canadian Bioblitz")

proj2b = proj2b[-which(proj2b$id %in% proj3b$id),]

## add Limestone Barrens Bioblitz ----------------------------------------------

proj4 = rinat::get_inat_obs_project(grpid = "limestone-barrens-of-newfoundland", type = "info")
proj4b = info_wrap(proj4, multi = FALSE, "Blitz the Gap")


## assemble into a table -------------------------------------------------------

df = bind_rows(list(proj1b, proj2b, proj3b, proj4b))
df = relocate(df, category)
df$id = as.numeric(df$id)
df$bioblitz = TRUE
# this isn't a concentrated effort bioblitz, it's an ongoing project
df$bioblitz[which(df$slug == "foray-nl-mushroom-lichen-diversity")] = FALSE 
df$bioblitz[which(df$slug == "swan-lake-christmas-hill-nature-sanctuary")] = FALSE 
df$bioblitz[which(df$slug == "limestone-barrens-of-newfoundland")] = FALSE 
df$bioblitz[which(df$slug == "marais-de-la-riviere-aux-cerises")] = FALSE 
df$bioblitz[which(df$slug == "mst-2025-fall-walks-and-forays")] = FALSE 

# add dates for the ongoing projects (dates for bioblitz efforts included in them)
df$date_start_if_ongoing = NA
df$date_end_if_ongoing = NA
# Foray
df$date_start_if_ongoing[which(df$slug == "foray-nl-mushroom-lichen-diversity")] = "22/09/2025"
df$date_end_if_ongoing[which(df$slug   == "foray-nl-mushroom-lichen-diversity")] = "26/09/2025"
# Swan Lake
df$date_start_if_ongoing[which(df$slug == "swan-lake-christmas-hill-nature-sanctuary")] = "21/09/2025"
df$date_end_if_ongoing[which(df$slug   == "swan-lake-christmas-hill-nature-sanctuary")] = "28/09/2025"
# Limestone barrens
df$date_start_if_ongoing[which(df$slug == "limestone-barrens-of-newfoundland")] = "07/07/2025"
df$date_end_if_ongoing[which(df$slug   == "limestone-barrens-of-newfoundland")] = "14/07/2025"
# Marais
df$date_start_if_ongoing[which(df$slug == "marais-de-la-riviere-aux-cerises")] = "21/09/2025"
df$date_end_if_ongoing[which(df$slug   == "marais-de-la-riviere-aux-cerises")] = "28/09/2025"
# Marais
df$date_start_if_ongoing[which(df$slug == "mst-2025-fall-walks-and-forays")] = "21/09/2025"
df$date_end_if_ongoing[which(df$slug   == "mst-2025-fall-walks-and-forays")] = "28/09/2025"

# rename columns
colnames(df)[which(colnames(df) == "taxa_count")] <- "obs_count_rinat"
# remove empty column from rinat...
df = df[,-which(colnames(df) == "taxa_number")]
# write to file!
write.csv(df, "outputs/bioblitzes/master-bioblitz-list.csv", row.names = FALSE)


## extract observations that are part of these bioblitzes ----------------------

# can download whole projects of the non-ongoing projects
# will have to manually download the date ranges for the ongoing projects

df_download = df |> filter(bioblitz == TRUE)

obs.list = vector("list", length = nrow(df_download))
names(obs.list) = df_download$slug
for(i in 1:nrow(df_download)) {
  obs.list[[i]] = get_inat_obs_project(grpid = df_download$id[i], 
                                       type = "observations")
}
obs.list = lapply(obs.list, function(x){
                  select(x, c("id", "uuid", 
                         "observed_on", "user_login", 
                         "longitude", "latitude", 
                         "quality_grade", "captive",
                         "iconic_taxon_name", 
                         "taxon.id", 
                         "taxon.name", "taxon.rank"))})
saveRDS(obs.list, "outputs/bioblitzes/master-bioblitz-obs-list.rds")


# flag the ones stopped by the API (they tap out at 10000 obs)
api_block = which(unlist(lapply(obs.list, nrow)) == 10000)
# "expedition-fiord-arctic-bioblitz"
# manually exported from iNaturalist:
# "data/heavy/bioblitzes/expedition-fiord-2025.zip"

# add these data to the list
bb1 = read.csv("data/heavy/bioblitzes/expedition-fiord-2025/observations-735804.csv")
colnames(bb1)

bb1 = rename(bb1, c("captive" = "captive_cultivated",
               "taxon.id" = "taxon_id",
               "taxon.name" = "taxon_species_name"))
bb1 = select(bb1, c("id", "uuid", 
                         "observed_on", "user_login", 
                         "longitude", "latitude", 
                         "quality_grade", "captive",
                         "iconic_taxon_name", 
                         "taxon.id", 
                         "taxon.name"))

obs.list = obs.list |> lapply(select, -c("taxon.rank"))

obs.list[[length(obs.list) + 1]] <- bb1
names(obs.list)[length(obs.list)] <- "expedition-fiord-arctic-bioblitz"


# Manually download specific dates for these ongoing projects ==================

# "foray-nl-mushroom-lichen-diversity"       
# "swan-lake-christmas-hill-nature-sanctuary"
# "limestone-barrens-of-newfoundland"    



# "data/heavy/bioblitzes/limestone-barrens-2025-07-0714.zip"