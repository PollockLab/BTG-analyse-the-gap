# Script to plot the effect of Blitz the Gap on observation rates

library(tidyverse)
library(arrow)

# load BTG member list (any user in any BTG related project) from btg-user-contributions.R
btg = read.csv("outputs/users/btg_users_groups.csv", row.names = 1)
length(unique(btg$user_login)) # 1123 people

# load parquet file 
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet") |>
  filter(year %in% c("2024", "2025"), month %in% c("06", "07", "08", "09"))

# load observer categories
user_categories = readRDS("outputs/users/observer_categories.rds")
user_new = read.csv("outputs/users/btg_users_newin2025.csv")


## Observation rate during BTG =================================================

# filter to users that participated in btg
inat.btg = inat |> 
  filter(user_login %in% unique(btg$user_login))

# get days active per user per year 
days.active = inat.btg |>  
  select(c(user_login, year, month, day)) |> 
  distinct() |>
  group_by(user_login, year) |>
  summarise(
    "days_active" = n()
  ) |> collect()

# get observations per year

obs.daily = inat.btg |>
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

# PER GROUP: summary
group_obs = obs.compare |>
  group_by(groups) |>
  summarise(
    "BTG" = sum(n_obs, na.rm = T),
    "Baseline" = sum(n_obs_expected, na.rm = T)
  ) |>
  mutate("BTG Effect" = BTG-Baseline,
         "ratio" = (100*BTG/Baseline)-100)
group_obs$ratio[which(group_obs$groups == "new")] <- NA

# long version for plotting
group_obs_l = group_obs |>
  select(-c(BTG, ratio)) |>
  pivot_longer(cols = c(`BTG Effect`, Baseline),
              names_to = "period",
              values_to = "n_obs")
# add row for the totals
total_tobind = data.frame(
  "groups" = c("all", "all"),
  "period" = c("BTG Effect", "Baseline"),
  "n_obs" = c(boost_total, exp_total)
)
group_obs_l = rbind(group_obs_l, total_tobind)
group_obs_l$period = factor(group_obs_l$period, levels = c("BTG Effect", "Baseline"))
group_obs_l$groups = factor(group_obs_l$groups, levels = c("new", "casual", "dabbler", 
                                                       "enthusiast", "superuser", "all"))
(A = ggplot(data = group_obs_l) +
    geom_bar(aes(
      x = groups,
      y = n_obs,
      fill = period
    ), stat = "identity", position = "stack") +
    geom_text(data = filter(group_obs, groups != "new"),
              aes(
                y = BTG,
                label = paste0("+",round(ratio),"%"),
                x = groups), 
              size = 5,  vjust = -1, fontface = "bold") +
    geom_text( y = btg_total,
               label = paste0("+",round((100*btg_total/exp_total)-100),"%"),
               x = "all", 
               size = 5,  vjust = -1, fontface = "bold") +
    scale_fill_manual(values = c("goldenrod1", "grey10"), name = "") +
    coord_cartesian(ylim = c(0,6e5)) +
    labs(y = "Observations", 
         x = "",
         title = "Blitz the Gap observers") +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14,
                               axis_title_face = "bold") +
    theme(legend.position = "top",
          panel.grid.major.y = element_blank()) ) 
ggsave("figures/btg-effect/barplot-btgeffect.png", width = 6.6, height = 6.03)





## NOTES =======================================================================
# baseline: expected number of observations given average efficiency
# baseline = rate * active_days

# boost: difference from baseline
# boost = btg - baseline

# totaled across everyone!

## BTG members: compare users to BTG 2025 and same window but 2024
