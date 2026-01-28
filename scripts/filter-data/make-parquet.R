# Script to make parquet file from the inat dataset

library(arrow)
library(tidyverse)

df = open_csv_dataset("data/heavy/inaturalist-canada-5-6Dec2025-nousergeopriv/inaturalist-canada-5-observations-no-usergeopriv.csv",
                      parse_options = list("newlines_in_values" = TRUE))
df

df_test = df |>
  group_by(iconic_taxon_name) |>
  summarise(n = n()) |> 
  collect()

arrow::write_parquet(df, "data/heavy/BTG-data/inaturalist-canada-dec2025.parquet")

# format to be lighter and to get a "year" column
df = df |>  
  dplyr::filter(captive_cultivated == FALSE,
                place_country_name == "Canada") |>
  select(c(observed_on_string, observed_on,
           user_login, quality_grade, 
           longitude, latitude, coordinates_obscured, 
           place_state_name, scientific_name, 
           common_name, iconic_taxon_name,
           taxon_id)) |>
  collect()
df = df |>
  separate(observed_on, into = c("year", "month", "day"), sep = "-")
# load parquet file
arrow::write_parquet(df, "data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")


# filter to 2025
df_2025 = df |>  
  dplyr::filter(captive_cultivated == FALSE,
                place_country_name == "Canada") |>
  dplyr::filter(grepl("2025", observed_on_string)) |>
  select(c(observed_on_string, observed_on,
           user_login, quality_grade, 
           longitude, latitude, coordinates_obscured, 
           place_state_name, scientific_name, 
           common_name, iconic_taxon_name,
           taxon_id)) |>
  collect()
arrow::write_parquet(df_2025, "data/heavy/BTG-data/inaturalist-canada-dec2025_subset2025.parquet")
