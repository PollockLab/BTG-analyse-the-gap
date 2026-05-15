# Script to make parquet file from the inat dataset

library(arrow)
library(tidyverse)

df = open_csv_dataset("data/heavy/inaturalist-canada-5-6Dec2025-nousergeopriv/inaturalist-canada-5-observations-no-usergeopriv.csv",
                      parse_options = list("newlines_in_values" = TRUE))
df$observed_on |> max(na.rm = T)
df$observed_on |> min(na.rm = T)
df$observed_on |> range(na.rm = T)

df_test = df |>
  group_by(iconic_taxon_name) |>
  summarise(n = n()) |> 
  collect()

arrow::write_parquet(df, "data/heavy/BTG-data/inaturalist-canada-dec2025.parquet")
df = arrow::read_parquet("data/heavy/BTG-data/inaturalist-canada-dec2025.parquet")

# format to be lighter and to get a "year" column
df = df |>  
  dplyr::filter(captive_cultivated == FALSE,
                place_country_name == "Canada") |>
  select(c(observed_on,
           user_login, 
           quality_grade, 
           longitude, latitude, 
           coordinates_obscured, 
           place_state_name, 
           scientific_name, 
           common_name, iconic_taxon_name,
           taxon_id)) |>
  collect()
df = df |>
  separate(observed_on, into = c("year", "month", "day"), sep = "-")
# load parquet file
arrow::write_parquet(df, "data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# read in pq file
df = arrow::read_parquet("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# filter to 2025
df_2025 = df |>
  dplyr::filter(year == "2025") |>
  dplyr::filter(quality_grade != "casual") # remove obs that are not verifiable
arrow::write_parquet(df_2025, "data/heavy/BTG-data/inaturalist-canada-dec2025_subset2025.parquet")

# filter to APR 1 TO OCT 1
df_APROCT = df_2025 |>
  dplyr::filter(month %in% c("04", "05", "06", "07", "08", "09", "10")) 
# remove days in oct after oct 1
df_APROCT = df_APROCT[-which(df_APROCT$month == "10" & df_APROCT$day != "01"),]
arrow::write_parquet(df_APROCT, "data/heavy/BTG-data/inaturalist-canada-dec2025_APR1OCT1.parquet")

# filter to JUNE 1 TO OCT 1
df_BTG = df_APROCT |> dplyr::filter(month != "04") 
arrow::write_parquet(df_BTG, "data/heavy/BTG-data/inaturalist-canada-dec2025_JUN1OCT1.parquet")
