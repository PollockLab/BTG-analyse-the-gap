# Script to classify users based on their number of observations

# load parquet file
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# count number of observations per observer
user.obs = inat |>
  group_by(user_login) |>
  summarise(
    "n_obs" = n()
  ) |> 
  collect()

# categorize
user.obs$user_category = NA
user.obs$user_category[which(user.obs$n_obs < 100)] = "casual"
user.obs$user_category[which(user.obs$n_obs >= 100)] = "enthusiast"
user.obs$user_category[which(user.obs$n_obs >= 5000)] = "superuser"
# check
table(user.obs$user_category)

# save
saveRDS(user.obs, "outputs/users/observer_categories.rds")
