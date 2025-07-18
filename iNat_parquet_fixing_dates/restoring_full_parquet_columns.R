# Goal: Take the iNat parquet with standardized dates and limited columns
# (the product of updating_dates_full_iNat_parquet.R), and add back in the columns
# to reform our original full parquet

# Ryan Hull, Quantitative Biodiversity Lab, McGill University
# July 2025

library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(stringr)


# Loading the parquets
rm(list=ls())
inat_pq_limited_cols <- arrow::open_dataset("C:/Users/Dell/OneDrive - McGill University/Laura's Lab_Group - Blitz the Gap/iNaturalist Canada parquet/iNat_non_sensitive_data_Jan2025_dates_standardized.parquet")
original_inat_pq <- arrow::open_dataset("C:/Users/Dell/OneDrive - McGill University/Laura's Lab_Group - Blitz the Gap/iNaturalist Canada parquet/iNat_non_sensitive_data_Jan2025.parquet")

limited_cols_df <- inat_pq_limited_cols |>
  select(c(id, observed_on_string, standardized_ddmmyyyy, is_ambiguous)) |>
  collect()

original_pq_df <- original_inat_pq |>
  collect()

names(limited_cols_df)[names(limited_cols_df) == "is_ambiguous"] <- "standardized_date_is_ambiguous"

# left joining. left join will complain about duplicates, so let's see them first.
obs_duplicates <- original_pq_df[duplicated(original_pq_df) | duplicated(original_pq_df, fromLast="TRUE")]

full_cols <- dplyr::left_join(limited_cols_df, original_pq_df, by="id")

