# Script to plot the effect of Blitz the Gap on observation rates

library(tidyverse)
library(arrow)

# load BTG member list (any user in any BTG related project) from btg-user-contributions.R
btg = read.csv("outputs/users/btg_users_groups.csv", row.names = 1)
length(unique(btg$user_login)) # 1123 people

# load parquet file 
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smAller.parquet") |>
  filter(year %in% c("2024", "2025"), month %in% c("06", "07", "08", "09"))

# load observer categories
user_categories = readRDS("outputs/users/observer_categories.rds")
user_new = read.csv("outputs/users/btg_users_newin2025.csv")

## Observation rate during BTG =================================================

# filter to users that participated in btg

inat.notbtg = inat |> collect() |>
  filter(!c(user_login %in% unique(btg$user_login)))

# get days active per user per year 
days.active = inat.notbtg |>  
  select(c(user_login, year, month, day)) |> 
  distinct() |>
  group_by(user_login, year) |>
  summarise(
    "days_active" = n()
  ) |> collect()


# get observations per year

obs.daily = inat.notbtg |>
  group_by(user_login, year, iconic_taxon_name, scientific_name) |>
  summarise(
    "n_obs" = n()) |> 
  ungroup() |>
  group_by(user_login, year) |>
  summarise(
    "n_obs" = sum(n_obs),
    "n_sp" = n()
  ) |>
  left_join(days.active) |>
  collect()

# calculate daily rate
obs.daily$rate = obs.daily$n_obs/obs.daily$days_active


# average daily rate
mean(obs.daily$rate) # 2.870469

# calculate baseline expectation of number of observations per person during BTG
# based on days active in 2025 x their average rate in 2024

# get average rate in 2024
rate.baseline = obs.daily |> 
  filter(year == "2024") |>
  select(c(user_login, rate, year)) |>
  rename("rate_2024" = "rate")

# multiply by active days in 2025
obs.compare = obs.daily |>
  filter(year == "2025") |>
  select(-c(year)) |>
  left_join(rate.baseline, by = "user_login") |>
  select(-c(year)) |>
  mutate("n_obs_expected" = rate_2024*days_active)

# add user classifications from 00_categorise_observers.R
obs.compare = obs.compare |> 
  left_join(select(user_categories, - c(n_obs)), by = "user_login")

# categorise new users
obs.compare$groups = obs.compare$user_category
obs.compare$groups[which(obs.compare$user_login %in% user_new$user_login)] = "new"

# TOTAL: summary
btg_total = sum(obs.compare$n_obs, na.rm = T)
exp_total = sum(obs.compare$n_obs_expected, na.rm = T)
boost_total = btg_total-exp_total
ratio_total = btg_total/exp_total
# PER GROUP: summary
group_obs = obs.compare |>
  summarise(
    "BTG" = sum(n_obs, na.rm = T),
    "Baseline" = sum(n_obs_expected, na.rm = T)
  ) |>
  mutate("BTG Effect" = BTG-Baseline,
         "ratio" = (BTG/Baseline)) |>
  mutate("New" = is.infinite(ratio))
group_obs$ratio[which(group_obs$New == TRUE)] <- NA

# long version for plotting
group_obs_l = group_obs  |>
  group_by(New) |>
  summarise("BTG" = sum(BTG, na.rm = T),
            "Baseline" = sum(Baseline, na.rm = T),
            "BTG Effect" = sum(`BTG Effect`, na.rm = T))
  select(-c(BTG, ratio)) |>
  pivot_longer(cols = c(`BTG Effect`, Baseline),
               names_to = "period",
               values_to = "n_obs") 
  group_obs_l$ratio = group_obs_l$BTG/group_obs_l$Baseline
saveRDS(group_obs_l, "outputs/btg-effect-figs/non-btg-members.rds")