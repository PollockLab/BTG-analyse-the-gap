## Script to quantile users by n_obs and n_sp they observe

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr)
library(duckdb)

# set ggplot theme
theme_set(ggpubr::theme_pubr())

# load parquet
inat_pq <- arrow::open_dataset("data/heavy/iNat_non_sensitive_data_Jan2025.parquet")


# exploring the data to figure out some useful cutoffs for categories ----------

# cutoff: how many observations does a user need to be included in the quantile?
cutoff = 50

# Count number of observations per user
user_summ <- inat_pq |> 
  dplyr::filter(captive_cultivated == FALSE) |> # should we filter out certain taxa here?
  group_by(user_login) |>
  summarize(n_obs = n()) |> 
  collect()

# quantiles
quants = quantile(user_summ$n_obs[which(user_summ$n_obs >= cutoff)],
                  probs = seq(0,1,0.1))
hist(log(user_summ$n_obs), freq = FALSE, col = "navyblue", border = "white")
abline(v = log(quants), col = "coral")


# count number of observations per species per user
user_summ_sp <- inat_pq |> 
  dplyr::filter(captive_cultivated == FALSE) |> 
  group_by(user_login, scientific_name) |>
  summarize(total_sp = n()) |> 
  collect()

# count number of species per user
user_summ_nsp = user_summ_sp |>
  group_by(user_login) |>
  summarise(n_sp = n()) |> 
  collect()
# join with n_obs table
users = left_join(user_summ, user_summ_nsp)
# quantiles
quants = quantile(users$n_sp[which(users$n_sp>=100)],
                  probs = c(seq(0,1,0.2)))
hist(log10(users$n_sp), freq = FALSE, col = "navyblue", border = "white")
abline(v = log10(quants), col = "coral")


## Categorise users ------------------------------------------------------------

## min number of species -----

sp_superuser = 5000
sp_expert = 1000
sp_enthusiast = 500
sp_dabbler = 100
# under: casual

# min number of observations -----

obs_superuser = 50000
obs_expert = 10000
obs_enthusiast = 1000
obs_dabbler = 100
# under: casual

# categorise by both criteria
users$category = NA
users$category[which(users$n_sp <  sp_dabbler | users$n_obs < obs_dabbler)] <- "casual"
users$category[which(users$n_sp >= sp_dabbler | users$n_obs >= obs_dabbler)] <- "dabbler"
users$category[which(users$n_sp >= sp_enthusiast | users$n_obs >= obs_enthusiast)] <- "enthusiast"
users$category[which(users$n_sp >= sp_expert | users$n_obs >= obs_expert)] <- "expert"
users$category[which(users$n_sp >= sp_superuser | users$n_obs >= obs_superuser)] <- "superuser"
users$category = factor(users$category,
                        levels = c("casual", "dabbler", "enthusiast", "expert", "superuser"))

# check out the categorisation
ggplot(data = users) +
  geom_boxplot(aes(x = category, y = n_sp))
ggplot(data = users) +
  geom_boxplot(aes(x = category, y = n_obs)) +
  scale_y_sqrt()

# plot the blocked user groups in terms of number of species and observations
ggplot(data = users) +
  geom_point(aes(x = n_obs, y = n_sp, col = category), alpha = .5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Observations", y = "Species", col = "User\ncategories") +
  theme(legend.position = "right")
ggsave("figures/usertype_comparisons/user_categories.png", width = 10, height = 10)

# save the user dataframe
saveRDS(users, "outputs/users/inat_users_categories.rds")