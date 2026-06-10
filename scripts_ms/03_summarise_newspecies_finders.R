# Script to summarise the users who found new species records during BTG

library(tidyverse)
library(sf)

# load new species records list
newsp = readRDS("outputs/missing-species/newsp_finders_sf.rds")
newsp$user_login |> unique() |> length()

# first: were any of these found during bioblitzes?
bioblitz = read.csv("outputs/bioblitzes/master-bioblitz-obs.csv", row.names = 1)
newsp.blitz = filter(bioblitz, taxon.name %in% newsp$scientific_name)
write.csv(newsp.blitz, "outputs/missing-species/newsp_finders_bioblitzes.csv")
unique(newsp.blitz$taxon.name)
# load user groups 
users = read.csv("outputs/users/btg_users_groups.csv", row.names = 1)

# join
users.simple = select(users, c(primary_group, user_login)) |> distinct()
df = left_join(newsp, users.simple)
df$primary_group[which(is.na(df$primary_group))] <- "iNaturalist"

# categorize the iNaturalist users ---------------------------------------------

users.cat = readRDS("outputs/users/observer_categories.rds") |> select(-c(n_obs))
df = left_join(df, users.cat)
df = sf::st_drop_geometry(df)
# save table
write.csv(df, "outputs/missing-species/newsp_finders_categorized.csv")

# summarise how many species were found per group ------------------------------

df.groups = df |> 
  group_by(primary_group, user_login, scientific_name) |>
  distinct() |>
  ungroup()

nsp.groups = df.groups |>
  group_by(primary_group) |>
  distinct(scientific_name) |>
  summarise("n_sp" = n())

df.groups2 = df.groups |>
  group_by(user_login, user_category) |>
  summarise("n_sp" = n(),
            "primary_group" = unique(primary_group))


## average new species for BC Biodiversity Program

bc.avg = df.groups2 |>
  filter(primary_group %in% c("BC Biodiversity", "BC Biodiversity 2025 Team")) 
data.frame("mean_sp" = mean(bc.avg$n_sp),
            "sd_sp" = sd(bc.avg$n_sp),
            "max_sp" = max(bc.avg$n_sp),)
# avg for other groups 
group.avg = df.groups2 |>
  group_by(primary_group) |>
  summarise("n_people" = n(),
            "mean_sp" = mean(n_sp),
            "sd_sp" = sd(n_sp),
            "max_sp" = max(n_sp))

# avg for inat user groups
inat.avg = df.groups2 |>
  filter(primary_group %in% c("iNaturalist")) |>
  group_by(user_category) |>
  summarise("n_people" = n(),
            "mean_sp" = mean(n_sp),
            "sd_sp" = sd(n_sp),
            "max_sp" = max(n_sp))

# how many were found during bioblitzes?

newsp.names = newsp.blitz |>
  filter(quality_grade == "research") |>
  distinct(taxon.name)
nrow(newsp.names)
