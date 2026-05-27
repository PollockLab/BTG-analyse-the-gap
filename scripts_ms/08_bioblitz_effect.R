# Script to plot the effect of bioblitzes on observation rates

library(tidyverse)
library(arrow)

# load bioblitz data
obs = read.csv("outputs/bioblitzes/master-bioblitz-obs.csv") |>
  filter(quality_grade != "casual",
         captive == FALSE)

# unique users
bioblitzers = unique(obs$user_login) # 1465
# user and date combos
id = paste(obs$user_login, obs$observed_on, sep = "_") |> unique()

# load parquet file for BTG 2025
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_JUN1OCT1.parquet")

# load observer categories
user_categories = readRDS("outputs/users/observer_categories.rds")
user_new = read.csv("outputs/users/btg_users_newin2025.csv")

## Observation rate during BTG =================================================

# filter to users that participated in btg
inat.btg = inat |> 
  filter(user_login %in% bioblitzers) |>
  mutate(ID = paste0(user_login, "_", year, "-", month, "-", day)) |>
  collect()

# subset to data included / not included in the bioblitzes
inat.btg = inat.btg |>
  mutate("bioblitz" = if_else(ID %in% id, TRUE, FALSE))
table(inat.btg$bioblitz)

# get days active per user per year 
days.active = inat.btg |>  
  select(c(user_login, year, month, day, bioblitz)) |> 
  distinct() |>
  group_by(user_login, bioblitz) |>
  summarise(
    "days_active" = n()
  ) |> collect()

# get observations per year

obs.daily = inat.btg |>
  group_by(user_login, bioblitz, iconic_taxon_name, scientific_name) |>
  summarise(
    "n_obs" = n()) |> 
  ungroup() |>
  group_by(user_login, bioblitz) |>
  summarise(
    "n_obs" = sum(n_obs),
    "n_sp" = n()
  ) |>
  left_join(days.active) |>
  collect()

# calculate daily rate
obs.daily$rate = obs.daily$n_obs/obs.daily$days_active

# calculate baseline expectation of number of observations per person during bioblitz
# based on days active in the bioblitz x their average rate for the rest of BTG

# get average rate in 2024
rate.baseline = obs.daily |> 
  filter(bioblitz == FALSE) |>
  select(c(user_login, rate)) |>
  rename("rate_btg" = "rate")

# multiply by active days in 2025
obs.compare = obs.daily |>
  filter(bioblitz == TRUE) |>
  left_join(rate.baseline) |>
  mutate("n_obs_expected" = rate_btg*days_active)

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
    "Bioblitz" = sum(n_obs, na.rm = T),
    "Baseline" = sum(n_obs_expected, na.rm = T)
  ) |>
  mutate("Bioblitz Effect" = Bioblitz-Baseline,
         "ratio" = (Bioblitz/Baseline))
group_obs$groups = str_to_title(group_obs$groups)

# long version for plotting
group_obs_l = group_obs |>
  select(-c(Bioblitz, ratio)) |>
  pivot_longer(cols = c(`Bioblitz Effect`, Baseline),
               names_to = "period",
               values_to = "n_obs")
# add row for the totals
total_tobind = data.frame(
  "groups" = c("All", "All"),
  "period" = c("Bioblitz Effect", "Baseline"),
  "n_obs" = c(boost_total, exp_total)
)
group_obs_l = rbind(group_obs_l, total_tobind)
group_obs_l$period = factor(group_obs_l$period, levels = c("Bioblitz Effect", "Baseline"))
group_obs_l$groups = str_to_title(group_obs_l$groups)
group_obs_l$groups = factor(group_obs_l$groups, 
                            levels = c("New", "Casual", "Dabbler", 
                                       "Enthusiast", "Superuser", "All"))
(A = ggplot(data = group_obs_l) +
    geom_bar(aes(
      x = groups,
      y = n_obs,
      fill = period
    ), stat = "identity", position = "stack") +
    geom_text(data = group_obs,
              aes(
                y = Bioblitz,
                label = paste0(signif(ratio, 2),"x"),
                x = groups), 
              size = 5,  vjust = -1, fontface = "bold") +
    geom_text( y = btg_total,
               label = paste0(signif((btg_total/exp_total), digits = 2),"x"),
               x = "All", 
               size = 5,  vjust = -1, fontface = "bold") +
    scale_fill_manual(values = c("#8FBCBB", "grey10"), name = "") +
    scale_y_continuous(labels = scales::label_number(scale = 0.001, suffix = "k")) +
    coord_cartesian(ylim = c(0,7e4)) +
    labs(y = "Observations", 
         x = "Bioblitz members") +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14,
                               axis_title_face = "bold") +
    theme(legend.position = "top",
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank()) ) 
ggsave("figures/btg-effect/barplot-bioblitzeffect.png", width = 6.6, height = 6.03)
saveRDS(A, "outputs/btg-effect-figs/barplot-bioblitzeffect.rds")



