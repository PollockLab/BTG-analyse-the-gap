# Script to check for 2025 data in the iNat pull
library(data.table)

dat = data.table::fread("data/heavy/inaturalist-canada-5-6Dec2025-nousergeopriv/inaturalist-canada-5-observations-no-usergeopriv.csv",
                        select = c("observed_on", "observed_on_string"))
# get latest date
dat$observed_on |> max(na.rm = T)

dat_2025 = dat[grep("2025", dat$observed_on_string),]
row_to_keep = grep("2025", dat$observed_on_string)

data.table::fread("data/heavy/inaturalist-canada-5-6Dec2025-nousergeopriv/inaturalist-canada-5-observations-no-usergeopriv.csv",
                        nrow = 1) |> colnames()


df = data.table::fread("data/heavy/inaturalist-canada-5-6Dec2025-nousergeopriv/inaturalist-canada-5-observations-no-usergeopriv.csv",
                        select = c("observed_on", "observed_on_string", 
                                   "latitude", "longitude", "geoprivacy", 
                                   "private_latitude", "private_longitude",
                                   "scientific_name", "common_name",
                                   "iconic_taxon_name", "taxon_id", "id"))[row_to_keep,]
write.csv(df, "data/heavy/BTG-data/inaturalist-canada-2025.csv")