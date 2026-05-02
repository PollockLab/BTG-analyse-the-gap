# Script to filter the GBIF April-Oct 2025 data download
# (research-grade only)

# load libraries
library(data.table)
library(tidyverse)

# fread the big dataset
dat = data.table::fread("data/heavy/gbif/0019803-260108223611665/occurrence.txt")

# select columns we need
dat_filter = dat |>
  select(c(year, month, day, 
           stateProvince,
           decimalLatitude, decimalLongitude, 
           coordinateUncertaintyInMeters,
           taxonID, scientificName, kingdom:genus, taxonRank, species, 
           speciesKey, iucnRedListCategory))
saveRDS(dat_filter, "data/heavy/gbif/iNatCAN_AprOct2025_selcols.rds")
write.csv(dat_filter, "data/heavy/gbif/iNatCAN_AprOct2025_selcols.csv", row.names = FALSE)