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

# number of observers during BTG ===============================================

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
         "ratio" = (BTG/Baseline))
group_obs$ratio[which(group_obs$groups == "new")] <- NA
group_obs$groups = str_to_title(group_obs$groups)

# long version for plotting
group_obs_l = group_obs |>
  select(-c(BTG, ratio)) |>
  pivot_longer(cols = c(`BTG Effect`, Baseline),
               names_to = "period",
               values_to = "n_obs")
# add row for the totals
total_tobind = data.frame(
  "groups" = c("All", "All"),
  "period" = c("BTG Effect", "Baseline"),
  "n_obs" = c(boost_total, exp_total)
)
group_obs_l = rbind(group_obs_l, total_tobind)
group_obs_l$period = factor(group_obs_l$period, levels = c("BTG Effect", "Baseline"))
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
    geom_text(data = filter(group_obs, groups != "New"),
              aes(
                y = BTG,
                label = paste0(signif(ratio, 2),"x"),
                x = groups), 
              size = 5,  vjust = -1, fontface = "bold") +
    geom_text( y = btg_total,
               label = paste0(signif((btg_total/exp_total), digits = 2),"x"),
               x = "All", 
               size = 5,  vjust = -1, fontface = "bold") +
    scale_fill_manual(values = c("#A3BE8C", "grey10"), name = "") +
    scale_y_continuous(labels = scales::label_number(scale = 0.001, suffix = "k")) +
    coord_cartesian(ylim = c(0,6e5)) +
    labs(y = "Observations", 
         x = "BTG members") +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14,
                               axis_title_face = "bold") +
    theme(legend.position = "top",
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()))  
ggsave("figures/btg-effect/barplot-btgeffect-nonmembers.png", width = 6.6, height = 6.03)


saveRDS(A, "outputs/btg-effect-figs/barplot-btgeffect-nonmembers.rds")
