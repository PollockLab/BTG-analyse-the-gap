# Script to plot the effect of grants on observation rates

library(tidyverse)
library(arrow)

# load parquet file 
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet") |>
  filter(year %in% c("2018", "2024", "2025"), month %in% c("06", "07", "08", "09"))

# load funded observers
funded = read.csv("outputs/users/btg_users_participants.csv", row.names = 1) |>
  filter(primary_group %in% c("QCBS Grantees", "BC Biodiversity 2025 Team", "GCB Grantees")) |>
  select(-group) |> distinct()

# unique users
fundees = unique(funded$user_login) # 131

# bc biodiversity program users who had accounts in 2018
bc = funded |> filter(primary_group == "BC Biodiversity 2025 Team")

pre.bc = inat |> 
  filter(year == "2018", 
         user_login %in% bc$user_login) |> 
  collect()

bc$user_login %in% unique(pre.bc$user_login) # 3 people did not have accounts

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
rate.BCbaseline = obs.daily |> 
  filter(year == "2018",
         primary_group == "BC Biodiversity 2025 Team") |>
  select(c(user_login, rate, year, primary_group)) |>
  rename("rate_2018" = "rate")

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

# add 2018 baseline for the BC Team
obs.compare = obs.compare |> 
  left_join(rate.BCbaseline, 
            by = c("user_login", "primary_group")) |>
  mutate("n_obs_expected_2018baseline" = rate_2018*days_active)
index = which(obs.compare$primary_group == "BC Biodiversity 2025 Team")
obs.compare$n_obs_expected[index] <- obs.compare$n_obs_expected_2018baseline[index]

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
# add row for the totals
total_tobind = data.frame(
  "groups" = c("All"),
  "Bioblitz" = btg_total,
  "Baseline" = exp_total,
  "`Bioblitz Effect`" = btg_total-exp_total,
  "ratio" = btg_total/exp_total
)
colnames(total_tobind) = colnames(group_obs)
group_obs = rbind(group_obs, total_tobind)
group_obs$primary_group = gsub(" Grantees", "", group_obs$primary_group) 
group_obs$primary_group[which(group_obs$primary_group == "BC Biodiversity 2025 Team")] <- "BCBP"

# long version for plotting
group_obs_l = group_obs |>
  select(-c(Bioblitz, ratio)) |>
  pivot_longer(cols = c(`Funding Effect`, Baseline),
               names_to = "period",
               values_to = "n_obs")
group_obs_l$period = factor(group_obs_l$period, levels = c("Funding Effect", "Baseline"))
group_obs_l$primary_group = factor(group_obs_l$primary_group, 
                                   levels = c("QCBS", 
                                              "GCB", 
                                              "BCBP", 
                                              "All"))
group_obs$primary_group = factor(group_obs$primary_group, 
                                   levels = c("QCBS", 
                                              "GCB", 
                                              "BCBP", 
                                              "All"))
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
    scale_y_continuous(labels = scales::label_number(scale = 0.001, suffix = "k")) +
    labs(y = "Observations", x = "") +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14,
                               axis_title_face = "bold") +
    theme(legend.position = "top",
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()))
ggsave("figures/btg-effect/barplot-fundingeffect.png", width = 6.6, height = 6.03)
saveRDS(A, "outputs/btg-effect-figs/barplot-fundingeffect.rds")

(A_nolabels = ggplot(data = group_obs_l) +
    geom_bar(aes(
      x = primary_group,
      y = n_obs,
      fill = period
    ), stat = "identity", position = "stack") +
    scale_fill_manual(values = c("#B48EAD", "grey10"), name = "") +
    coord_cartesian(ylim = c(0,2e5)) +
    # wrap long text in the x axis
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    scale_y_continuous(labels = scales::label_number(scale = 0.001, suffix = "k")) +
    labs(y = "Observations") +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14,
                               axis_title_face = "bold") +
    theme(legend.position = "top",
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()))
ggsave("figures/btg-effect/barplot-fundingeffect-nolabels.png", width = 6.6, height = 6.03)
saveRDS(A_nolabels, "outputs/btg-effect-figs/barplot-fundingeffect-nolabels.rds")

(B = ggplot(data = group_obs) +
    geom_bar(aes(
      x = primary_group,
      y = ratio), 
      stat = "identity", 
      position = "stack",
      fill = "#B48EAD") +
    geom_hline(yintercept =  2, lty = 2) +
    labs(y = "Observations", 
         x = "") +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14,
                               axis_title_face = "bold") +
    theme(legend.position = "none",
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank()))  
ggsave("figures/btg-effect/barplot-fundingeffect-ratio.png", width = 6.6, height = 6.03)
saveRDS(B, "outputs/btg-effect-figs/barplot-fundingeffect-ratio.rds")
