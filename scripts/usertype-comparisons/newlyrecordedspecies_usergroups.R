# Who found the new species?

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(sf)
library(patchwork)

# load parquet file
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# read the species list
sp_2025 = read.csv("outputs/found_species_in_2025.csv")

# tally findings by user logins
df = inat |>
  filter(scientific_name %in% sp_2025$scientific_name,
         iconic_taxon_name != "NA",
         quality_grade == "research",
         year == "2025") |>
  group_by(scientific_name, user_login) |>
  summarise("n" = n()) |> 
  collect()

# unique users
users = unique(df$user_login)

# number of observations and species, years active, etc.
users.info = inat |> 
  filter(user_login %in% users) |>
  group_by(user_login, scientific_name, iconic_taxon_name, year) |>
  summarise("obs" = n()) |>
  collect()

n.obs = users.info |>
  group_by(user_login) |>
  summarise("obs_total" = sum(obs, na.rm = TRUE))

n.sp = users.info |>
  group_by(user_login) |>
  distinct(scientific_name) |>
  summarise("sp_total" = n())

n.year = users.info |>
  group_by(user_login) |>
  distinct(year) |>
  mutate(year.numeric = as.numeric(year)) |>
  summarise("min.year" = min(year.numeric, na.rm = T),
            "years.active" = n(),
            "duration.active" = 2025-min(year.numeric, na.rm = T))

users.df = left_join(n.obs, n.sp) |> left_join(n.year)

new.sp = df |>
  group_by(user_login) |>
  summarise("new.sp" = n())
users.df = left_join(users.df, new.sp)

# calculate yearly rate of obs and species
users.df$yearly.obs = users.df$obs_total/users.df$years.active

# BC Biodiversity Program teams
bc = read.csv("data/user_information/bigteam_userlogins.csv", row.names = 1)
bc = bc$user_login

# QCBS users
# import qcbs logins
qcbs <- readRDS("~/Documents/GitHub/storymap-qcbs/data/championteams_userlogins.rds")
qcbs = c(qcbs, "katherinehebert") # i'm missing!
qcbs = gsub("and abrunet27", "abrunet27", qcbs)
write.csv(qcbs, "data/user_information/championteams_userlogins.csv")

# assign categories
users.df$category = "Other observers"
users.df$category[which(users.df$obs_total >= 1000)] <- "Enthusiast"
users.df$category[which(users.df$obs_total >= 5000)] <- "Superuser"
users.df$category[which(users.df$user_login %in% qcbs)] <- "QCBS"
users.df$category[which(users.df$user_login %in% bc)] <- "BC Biodiversity"
users.df$category = factor(users.df$category,
                           levels = (c("Other observers", "Enthusiast", 
                                       "Superuser", "QCBS", "BC Biodiversity")))
write.csv(users.df, "outputs/users/btg_users_newspecies.csv")



A = ggplot(data = users.df) +
  geom_jitter(aes(y = category, 
                  x = sp_total, 
                  colour = category),
              size = 3, alpha = .6, height = .1) +
  geom_boxplot(aes(x = sp_total, 
                   y = category,
                   fill = category), alpha = .6, 
               outliers = FALSE, width = .5) +
  colorspace::scale_fill_discrete_qualitative() + 
  colorspace::scale_color_discrete_qualitative() + 
  labs(x = "Total species observed (all-time) in iNaturalist Canada", 
      y = "",
       fill = "") + 
  hrbrthemes::theme_ipsum_rc(axis_text_size = 18,
                             axis_title_size = 20,
                             axis = F, 
                             base_size = 18) +
  theme(legend.position = "none") 


B = ggplot(data = users.df) +
  geom_jitter(aes(y = category, 
                  x = new.sp, 
                  colour = category),
              size = 3, alpha = .6, height = .1) +
  geom_boxplot(aes(x = new.sp, 
                   y = category,
                   fill = category), alpha = .6, 
               outliers = FALSE, width = .5) +
  colorspace::scale_fill_discrete_qualitative() + 
  colorspace::scale_color_discrete_qualitative() + 
  labs(x = "Newly-recorded species", 
       y = "",
       fill = "") +
  hrbrthemes::theme_ipsum_rc(axis_text_size = 18,
                             axis_title_size = 20,
                             axis = T, 
                             base_size = 12) +
  theme(legend.position = "none",
        axis.text.y = element_blank()) 

A + B + plot_annotation(tag_levels = "a") 
ggsave("figures/usertype_comparisons/newlyrecordedspecies_usergroups.png", width = 6.92, height = 7.42)


users.df |>
  group_by(category) |>
  summarise("mean" = mean(new.sp),
            "sd" = sd(new.sp),
            "max" = max(new.sp),
            "n" = n())
