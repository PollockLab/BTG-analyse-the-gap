# Script to plot the gains in number of species per taxa with at least 1, 10, 30, or 100 observations

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)

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
btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$n_min100.post/btg_temp$n_total.post)])


# Relative gain as a percentage ------------------------------------------------

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
  geom_text(aes(label = paste0("+", round(100*(n_min100.post-n_min100.pre)/n_total.post, digits = 1), "%"), 
                x = 100*n_min100.post/n_total.post, 
                y = iconic_taxon_name), 
            size = 5,  hjust = -.1, color = "white") +
  geom_text(aes(label = paste0(n_total.post, " sp"), 
                x = 100, 
                y = iconic_taxon_name), 
            size = 4,  hjust = +1.1, color = "grey40") +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Species (%)", 
       y = "",
       title = "Getting species to >100 observations!",
       subtitle = "Species that reached 100 observations during Blitz the Gap.") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min100_proportionoftotalspecies.png", width = 7, height = 6.21)


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
  geom_text(aes(label = paste0(n_total.post, " sp"), 
                x = 100, 
                y = iconic_taxon_name), 
            size = 4,  hjust = +1.1, color = "grey40") +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Species (%)", 
       y = "",
       title = "Getting species to >100 observations!",
       subtitle = "Species that reached 100 observations during Blitz the Gap.") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min100_proportionoftotalspecies_labelledn.png", width = 7, height = 6.21)


# Basic number of species with 100 obs plot ------------------------------------

btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$diff)])
(pp = ggplot(data = btg_temp) +
    geom_bar(aes(y = iconic_taxon_name, 
                 x = diff,
                 fill = iconic_taxon_name,
                 group = Species),
             stat = "identity", 
             alpha = 1) +
    geom_text(aes(label = diff, x = diff, y = iconic_taxon_name), 
              size = 5,  hjust = 1.2, color = "white") +
    colorspace::scale_color_discrete_qualitative() +
    scale_x_sqrt() +
    labs(x = "Species", 
         y = "",
         title = "Getting species to >100 observations!",
         subtitle = "Species that reached 100 observations during Blitz the Gap.") +
    theme(legend.position = "none",
          panel.grid.major.y = element_blank()))
ggsave("figures/n_observations/min100_numberofspecies.png", width = 7, height = 6.21)

## count number of species that have at least 100 obs now!
btg_temp$diff |> sum() # 589



# Minimum 10 observations -----------------------------------------------------

# Prepare data for plotting
btg_temp = filter(btg, !is.na(iconic_taxon_name))
btg_temp = btg_temp |>
  mutate("diff" = n_min10.post-n_min10.pre,
         "diffperc" = 100*(n_min10.post-n_min10.pre)/n_min10.post) |>
  mutate("Percentage" = paste0(round(diffperc, digits = 0), "%"),
         "Species" = paste0(diff, " sp."))
btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$n_min10.post/btg_temp$n_total.post)])


# Relative gain as a percentage ------------------------------------------------

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
  geom_text(aes(label = paste0("+", round(100*(n_min10.post-n_min10.pre)/n_total.post, digits = 1), "%"), 
                x = 100*n_min10.post/n_total.post, 
                y = iconic_taxon_name), 
            size = 5,  hjust = -.1, color = "white") +
  geom_text(aes(label = paste0(n_total.post, " sp"), 
                x = 100, 
                y = iconic_taxon_name), 
            size = 4,  hjust = +1.1, color = "grey40") +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Species (%)", 
       y = "",
       title = "Getting species to >30 observations!",
       subtitle = "Species that reached 30 observations during Blitz the Gap.") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min10_proportionoftotalspecies.png", width = 7, height = 6.21)


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
  geom_text(aes(label = paste0(n_total.post, " sp"), 
                x = 100, 
                y = iconic_taxon_name), 
            size = 4,  hjust = +1.1, color = "grey40") +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Species (%)", 
       y = "",
       title = "Getting species to >30 observations!",
       subtitle = "Species that reached 30 observations during Blitz the Gap.") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min10_proportionoftotalspecies_labelledn.png", width = 7, height = 6.21)


# Basic number of species with 100 obs plot ------------------------------------

btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$diff)])
(pp = ggplot(data = btg_temp) +
    geom_bar(aes(y = iconic_taxon_name, 
                 x = diff,
                 fill = iconic_taxon_name,
                 group = Species),
             stat = "identity", 
             alpha = 1) +
    geom_text(aes(label = diff, x = diff, y = iconic_taxon_name), 
              size = 5,  hjust = 1.2, color = "white") +
    colorspace::scale_color_discrete_qualitative() +
    scale_x_sqrt() +
    labs(x = "Species", 
         y = "",
         title = "Getting species to >30 observations!",
         subtitle = "Species that reached 30 observations during Blitz the Gap.") +
    theme(legend.position = "none",
          panel.grid.major.y = element_blank()))
ggsave("figures/n_observations/min10_numberofspecies.png", width = 7, height = 6.21)

## count number of species that have at least 100 obs now!
btg_temp$diff |> sum() # 842


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
  labs(x = "Species (%)", 
       y = "",
       title = "Getting species to >30 observations!",
       subtitle = "Species that reached 30 observations during Blitz the Gap.") +
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
  labs(x = "Species (%)", 
       y = "",
       title = "Getting species to >30 observations!",
       subtitle = "Species that reached 30 observations during Blitz the Gap.") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())
ggsave("figures/n_observations/min30_proportionoftotalspecies_labelledn.png", width = 7, height = 6.21)


# Basic number of species with 100 obs plot ------------------------------------

btg_temp$iconic_taxon_name = factor(btg_temp$iconic_taxon_name,
                                    levels = btg_temp$iconic_taxon_name[order(btg_temp$diff)])
(pp = ggplot(data = btg_temp) +
    geom_bar(aes(y = iconic_taxon_name, 
                 x = diff,
                 fill = iconic_taxon_name,
                 group = Species),
             stat = "identity", 
             alpha = 1) +
    geom_text(aes(label = diff, x = diff, y = iconic_taxon_name), 
              size = 5,  hjust = 1.2, color = "white") +
    colorspace::scale_color_discrete_qualitative() +
    scale_x_sqrt() +
    labs(x = "Species", 
         y = "",
         title = "Getting species to >30 observations!",
         subtitle = "Species that reached 30 observations during Blitz the Gap.") +
    theme(legend.position = "none",
          panel.grid.major.y = element_blank()))
ggsave("figures/n_observations/min30_numberofspecies.png", width = 7, height = 6.21)

## count number of species that have at least 100 obs now!
btg_temp$diff |> sum() # 709
