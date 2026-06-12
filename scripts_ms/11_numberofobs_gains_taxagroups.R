# Script to plot the count number of species per taxa with at least 1, 10, 30, or 100 observations

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)


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


# Plotting =====================================================================

# load data
btg = readRDS("outputs/n_observations/n_obs_btg_summary.rds")
btg = filter(btg, !is.na(iconic_taxon_name))

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# Prepare data for plotting
btg_temp = filter(btg, !is.na(iconic_taxon_name))
btg_temp = btg_temp |>
  mutate("diff" = n_min100.post-n_min100.pre,
         "diffperc" = 100*(n_min100.post-n_min100.pre)/n_min100.post) |>
  mutate("Percentage" = paste0(round(diffperc, digits = 0), "%"),
         "Species" = paste0(diff, " sp."))
btg_temp$iconic_taxon_name = gsub("Animalia", "Other Animalia", btg_temp$iconic_taxon_name)
btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$n_min100.post/btg_temp$n_total.post)])

# manually adding total species counts from inat.ca filtering all observations 
# on june 12, 2026 (as not all taxa in the parquet are species level)
# taking the Other Animalia count from the inat parquet since it is a catch-all group
btg_temp$total_sp = c(574, 54, 1659, 1398, 759, 891, 9064, 19460, 215, 950, 9086, 519, 65)

# Relative gain with n species labels ------------------------------------------

ggplot(data = btg_temp, 
       aes(fill = iconic_taxon_name)) +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100),
           stat = "identity", fill = "black") +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100*n_min100.pre/n_total.post),
           stat = "identity", alpha = 1) +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100*n_min100.post/n_total.post),
           stat = "identity", alpha = .8) +
  geom_text(aes(label = paste0("+", (n_min100.post-n_min100.pre)), 
                x = 100*n_min100.post/n_total.post, 
                y = iconic_taxon_name), 
            size = 5,  hjust = -.2, color = "white") +
  geom_text(aes(label = paste0(n_total.post, " "), 
                x = 100, 
                y = iconic_taxon_name), 
            size = 4,  hjust = +1.1, color = "grey40") +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Species (%)", 
       y = "",
       title = "Reached 100 observations") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min100_proportionoftotalspecies_labelledn.png", width = 7.94, height = 5.73)


# Basic number of species with 100 obs plot ------------------------------------

btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$diff)])


## count number of species that have at least 100 obs now!
btg_temp$diff |> sum() # 423

# Minimum 10 observations -----------------------------------------------------

# Prepare data for plotting
btg_temp = filter(btg, !is.na(iconic_taxon_name))
btg_temp$total_sp = c(574, 54, 1659, 1398, 759, 891, 9064, 19460, 215, 950, 9086, 519, 65)
btg_temp = btg_temp |>
  mutate("diff" = n_min10.post-n_min10.pre,
         "diffperc" = 100*(n_min10.post-n_min10.pre)/n_total.post) |>
  mutate("Percentage" = paste0(round(diffperc, digits = 0), "%"),
         "Species" = paste0(diff, " sp."))
btg_temp$iconic_taxon_name = gsub("Animalia", "Other Animalia", btg_temp$iconic_taxon_name)
btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$n_min10.post/btg_temp$n_total.post)])


# Relative gain with n species labels ------------------------------------------

ggplot(data = btg_temp, 
       aes(fill = iconic_taxon_name)) +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100),
           stat = "identity", fill = "black") +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100*n_min10.pre/n_total.post),
           stat = "identity", alpha = 1) +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100*n_min10.post/n_total.post),
           stat = "identity", alpha = .8) +
  geom_text(aes(label = paste0("+", (n_min10.post-n_min10.pre)), 
                x = 100*n_min10.post/n_total.post, 
                y = iconic_taxon_name), 
            size = 5,  hjust = -.2, color = "white") +
  geom_text(aes(label = paste0(n_total.post, ""), 
                x = 100, 
                y = iconic_taxon_name), 
            size = 4,  hjust = +1.1, color = "grey40") +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Species (%)", 
       y = "",
       title = "Reached 10 observations") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min10_proportionoftotalspecies_labelledn.png", width = 7, height = 6.21)


# Basic number of species with 100 obs plot ------------------------------------

btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$diff)])
## count number of species that have at least 100 obs now!
btg_temp$diff |> sum() # 601


# Minimum 30 observations ------------------------------------------------------

# Prepare data for plotting
btg_temp = filter(btg, !is.na(iconic_taxon_name))
btg_temp = btg_temp |>
  mutate("diff" = n_min30.post-n_min30.pre,
         "diffperc" = 100*(n_min30.post-n_min30.pre)/n_min30.post) |>
  mutate("Percentage" = paste0(round(diffperc, digits = 0), "%"),
         "Species" = paste0(diff, " sp."))
btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$n_min30.post/btg_temp$n_total.post)])


# Relative gain as a percentage ------------------------------------------------

ggplot(data = btg_temp, 
       aes(fill = iconic_taxon_name)) +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100),
           stat = "identity", fill = "black") +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100*n_min30.pre/n_total.post),
           stat = "identity", alpha = 1) +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100*n_min30.post/n_total.post),
           stat = "identity", alpha = .8) +
  geom_text(aes(label = paste0("+", round(100*(n_min30.post-n_min30.pre)/n_total.post, digits = 1), "%"), 
                x = 100*n_min30.post/n_total.post, 
                y = iconic_taxon_name), 
            size = 5,  hjust = -.1, color = "white") +
  geom_text(aes(label = paste0(n_total.post, " sp"), 
                x = 100, 
                y = iconic_taxon_name), 
            size = 4,  hjust = +1.1, color = "grey40") +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Taxa (%)", 
       y = "",
       title = "Reached 30 observations") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min30_proportionoftotalspecies.png", width = 7, height = 6.21)


# Relative gain with n species labels ------------------------------------------

ggplot(data = btg_temp, 
       aes(fill = iconic_taxon_name)) +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100),
           stat = "identity", fill = "black") +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100*n_min30.pre/n_total.post),
           stat = "identity", alpha = 1) +
  geom_bar(aes(y = iconic_taxon_name, 
               x = 100*n_min30.post/n_total.post),
           stat = "identity", alpha = .8) +
  geom_text(aes(label = paste0("+", (n_min30.post-n_min30.pre)), 
                x = 100*n_min30.post/n_total.post, 
                y = iconic_taxon_name), 
            size = 5,  hjust = -.2, color = "white") +
  geom_text(aes(label = paste0(n_total.post, " sp"), 
                x = 100, 
                y = iconic_taxon_name), 
            size = 4,  hjust = +1.1, color = "grey40") +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Taxa (%)", 
       y = "",
       title = "Reached 30 observations") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min30_proportionoftotalspecies_labelledn.png", width = 7, height = 6.21)