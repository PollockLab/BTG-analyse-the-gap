# Script to plot the effect of grants on observation rates

library(tidyverse)
library(arrow)

# load parquet file 
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smAller.parquet") |>
  filter(year %in% c("2024", "2025"), month %in% c("06", "07", "08", "09"))

# load funded observers
funded = read.csv("outputs/users/btg_users_participants.csv", row.names = 1) |>
  filter(primary_group %in% c("QCBS Grantees", "BC Biodiversity 2025 Team", "GCB Grantees")) |>
  select(-group)

# unique users
fundees = unique(funded$user_login) # 131

## Observation rate during BTG =================================================

# filter to users that participated in btg
inat.btg = inat |> 
  filter(user_login %in% fundees)

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
  collect() |>
  distinct()

# calculate daily rate
obs.daily$rate = obs.daily$n_obs/obs.daily$days_active

# attach funding categories 
obs.daily = left_join(obs.daily, funded, multiple = "any")

# calculate baseline expectation of number of observations per person during BTG
# based on days active in 2025 x their average rate in 2024

# get average rate in 2024
rate.baseline = obs.daily |> 
  filter(year == "2024") |>
  select(c(user_login, rate, year, primary_group)) |>
  rename("rate_2024" = "rate")

# multiply by active days in 2025
obs.compare = obs.daily |>
  filter(year == "2025") |>
  select(-c(year)) |>
  left_join(rate.baseline, by = c("user_login", "primary_group")) |>
  select(-c(year)) |>
  mutate("n_obs_expected" = rate_2024*days_active)
obs.compare$rate_2024[which(is.na(obs.compare$rate_2024))] <- 0
obs.compare$n_obs_expected[which(is.na(obs.compare$rate_2024))] <- 0

# TOTAL: summary
btg_total = sum(obs.compare$n_obs, na.rm = T)
exp_total = sum(obs.compare$n_obs_expected, na.rm = T)
boost_total = btg_total-exp_total

# PER GROUP: summary
group_obs = obs.compare |>
  group_by(primary_group) |>
  summarise(
    "Bioblitz" = sum(n_obs, na.rm = T),
    "Baseline" = sum(n_obs_expected, na.rm = T)
  ) |>
  mutate("Funding Effect" = Bioblitz-Baseline,
         "ratio" = (Bioblitz/Baseline))

# long version for plotting
group_obs_l = group_obs |>
  select(-c(Bioblitz, ratio)) |>
  pivot_longer(cols = c(`Funding Effect`, Baseline),
               names_to = "period",
               values_to = "n_obs")
# add row for the totals
total_tobind = data.frame(
  "primary_group" = c("All", "All"),
  "period" = c("Funding Effect", "Baseline"),
  "n_obs" = c(boost_total, exp_total)
)
group_obs_l = rbind(group_obs_l, total_tobind)
group_obs_l$period = factor(group_obs_l$period, levels = c("Funding Effect", "Baseline"))
group_obs_l$primary_group = factor(group_obs_l$primary_group, 
                                   levels = c("QCBS Grantees", "GCB Grantees", "BC Biodiversity 2025 Team", "All"))
(A = ggplot(data = group_obs_l) +
    geom_bar(aes(
      x = primary_group,
      y = n_obs,
      fill = period
    ), stat = "identity", position = "stack") +
    geom_text(data = group_obs,
              aes(
                y = Bioblitz,
                label = paste0(signif(ratio, 2),"x"),
                x = primary_group), 
              size = 5,  vjust = -1, fontface = "bold") +
    geom_text( y = btg_total,
               label = paste0(signif((btg_total/exp_total), digits = 2),"x"),
               x = "All", 
               size = 5,  vjust = -1, fontface = "bold") +
    scale_fill_manual(values = c("#B48EAD", "grey10"), name = "") +
    coord_cartesian(ylim = c(0,2e5)) +
    # wrap long text in the x axis
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    labs(y = "Observations", 
         x = "Funding") +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14,
                               axis_title_face = "bold") +
    theme(legend.position = "top",
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),) ) 
ggsave("figures/btg-effect/barplot-fundingeffect.png", width = 6.6, height = 6.03)

saveRDS(A, "outputs/btg-effect-figs/barplot-fundingeffect.rds")
