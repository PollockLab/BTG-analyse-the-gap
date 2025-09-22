# Script to plot the count number of species per taxa with at least 1, 10, 30, or 100 observations

# libraries
library(rgbif)
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(taxize)
library(rinat)

# load parquet
inat_pq <- arrow::open_dataset("~/McGill University/Laura's Lab_Group - BioBlitz/data/raw/biodiversity-data/inat-canada/iNat_non_sensitive_data_Jan2025.parquet")

# get iNat project IDs
btg_info = get_inat_obs_project("blitz-the-gap", type = "info", raw = FALSE)


## Before BTG ==================================================================

# Count number of observations per species, per taxa group
n_obs <- inat_pq |> 
  dplyr::filter(captive_cultivated == FALSE,
                quality_grade == "research") |> 
  group_by(iconic_taxon_name, scientific_name) |>
  summarize(total_obs = n()) |> 
  # load the query into our R session
  collect()

which(n_obs$total_obs >= 1) |> length()
which(n_obs$total_obs >= 10) |> length()
which(n_obs$total_obs >= 30) |> length()
which(n_obs$total_obs >= 100) |> length()

## After BTG ===================================================================

# get taxa information from the Blitz the Gap project
# btg_taxa = get_inat_taxon_stats(project = btg_info$id,
#                                 date_range = c("2025-04-01", 
#                                                # today's date
#                                                paste0(Sys.Date())))


# so this is a conservative estimate of new data... the rinat package isn't
# well suited to download obs per project, so I will have to go with 
# research-grade observations summarised using rgbif

new_obs = n_obs
new_obs$btg_obs = NA

for(i in 1:nrow(new_obs)){
  new_obs$btg_obs[i] = rgbif::occ_count(scientificName = n_obs$scientific_name[i],
                                        country = "CA",
                                        year = 2025,
                                        month = "4;5;6;7;8;9",
                                        datasetKey = "50c9509d-22c7-4a22-a47d-8c48425ef4a7")
}
saveRDS(new_obs, "outputs/n_observations/n_obs_btg_15092025.rds")


# Compare before and after BTG =================================================

# Before

pre.btg = n_obs
pre.btg$min10 = pre.btg$total_obs >= 10
pre.btg$min30 = pre.btg$total_obs >= 30
pre.btg$min100 = pre.btg$total_obs >= 100

pre.btg = pre.btg |>
  group_by(iconic_taxon_name) |>
  summarise("n_min10" = sum(min10, na.rm = T),
            "n_min30" = sum(min30, na.rm = T),
            "n_min100" = sum(min100, na.rm = T),
            "n_total" = n())

# After

post.btg = new_obs
post.btg$new_total = post.btg$total_obs + post.btg$btg_obs
post.btg$min10 = post.btg$new_total >= 10
post.btg$min30 = post.btg$new_total >= 30
post.btg$min100 = post.btg$new_total >= 100

post.btg = post.btg |>
  group_by(iconic_taxon_name) |>
  summarise("n_min10" = sum(min10, na.rm = T),
            "n_min30" = sum(min30, na.rm = T),
            "n_min100" = sum(min100, na.rm = T),
            "n_total" = n())

btg = left_join(pre.btg, post.btg, 
                by = "iconic_taxon_name",
                suffix = c(".pre", ".post"))
saveRDS(btg, "outputs/n_observations/n_obs_btg_15092025_summary.rds")

# Plotting is done in plot_n_observations.R