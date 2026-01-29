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
df <- arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025.parquet")

# make a dataset without 2025
df_no2025 = df |> dplyr::filter(!grepl("2025", observed_on_string))


## Before BTG ==================================================================

# Count number of observations per species, per taxa group
n_obs <- df_no2025 |> 
  dplyr::filter(captive_cultivated == FALSE,
                quality_grade == "research",
                place_country_name == "Canada") |> 
  group_by(iconic_taxon_name, scientific_name) |>
  summarize(total_obs = n()) |> 
  # load the query into our R session
  collect()
# retain taxa that have two words in the name, i.e. that are species or lower level
n_obs$species_level <- strsplit(n_obs$scientific_name,split = " ") |> 
  lapply(length) |> 
  unlist()
# convert to true/false
n_obs$species_level = n_obs$species_level>1
# filter to minimum species-level obs
n_obs = filter(n_obs, species_level == TRUE)
which(n_obs$total_obs >= 1) |> length()
which(n_obs$total_obs >= 10) |> length()
which(n_obs$total_obs >= 30) |> length()
which(n_obs$total_obs >= 100) |> length()

## After BTG ===================================================================

# Count number of observations per species, per taxa group
n_obs.2025 <- df |> 
  dplyr::filter(captive_cultivated == FALSE,
                quality_grade == "research",
                place_country_name == "Canada") |> # blitz the gap months
  group_by(iconic_taxon_name, scientific_name) |>
  summarize(total_obs = n()) |> 
  # load the query into our R session
  collect()
# retain taxa that have two words in the name, i.e. that are species or lower level
n_obs.2025$species_level <- strsplit(n_obs.2025$scientific_name,split = " ") |> 
  lapply(length) |> 
  unlist()
# convert to true/false
n_obs.2025$species_level = n_obs.2025$species_level>1
# filter to minimum species-level obs
n_obs.2025 = filter(n_obs.2025, species_level == TRUE)

which(n_obs.2025$total_obs >= 1) |> length()
which(n_obs.2025$total_obs >= 10) |> length()
which(n_obs.2025$total_obs >= 30) |> length()
which(n_obs.2025$total_obs >= 100) |> length()


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

post.btg = n_obs.2025
post.btg$min10 = post.btg$total_obs >= 10
post.btg$min30 = post.btg$total_obs >= 30
post.btg$min100 = post.btg$total_obs >= 100

post.btg = post.btg |>
  group_by(iconic_taxon_name) |>
  summarise("n_min10" = sum(min10, na.rm = T),
            "n_min30" = sum(min30, na.rm = T),
            "n_min100" = sum(min100, na.rm = T),
            "n_total" = n())

btg = left_join(pre.btg, post.btg, 
                by = "iconic_taxon_name",
                suffix = c(".pre", ".post"))
saveRDS(btg, "outputs/n_observations/n_obs_btg_summary.rds")

# Plotting is done in plot_n_observations.R

