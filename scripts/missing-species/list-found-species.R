# Script to count number of species per taxa that were previously missing in iNat Canada

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(sf)

# load manually downloaded project data
# animals = read.csv("data/heavy/missing-species/missing-animals/observations-619233.csv")
# insects = read.csv("data/heavy/missing-species/missing-insects.csv/observations-619217.csv")
# plants = read.csv("data/heavy/missing-species/missing-plants/observations-619220.csv")
# fungi = read.csv("data/heavy/missing-species/missing-fungi/observations-619222.csv")

# load parquet file
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# find species that were not found before 2025 ---------------------------------

missing_sp = inat |>
  filter(quality_grade == "research") |>
  group_by(iconic_taxon_name, scientific_name, year) |>
  distinct(year) |>
  na.omit() |>
  collect()
missing_sp$year = as.numeric(missing_sp$year)

sp_minyear = missing_sp |>
  group_by(scientific_name) |>
  summarise("min_year" = min(year),
            "max_year" = max(year))

# keep species whose minimum year is in 2025
sp_2025 = sp_minyear |>
  filter(min_year == 2025)

# filter to species names (genus + species)
temp = stringr::str_count(sp_2025$scientific_name, "\\w+")
# find species with nothing in the species column
to_remove = which(temp != 2)
to_remove = c(to_remove, grep("×", sp_2025$scientific_name)) # note: this is not a normal "x"
to_remove = c(to_remove, grep("Mustela furo", sp_2025$scientific_name)) # note: this is not a normal "x"

# remove species that have > 2 words, or that include the hybrid symbol
sp_2025 <- sp_2025[-to_remove,]
rm(temp)

# once again, check that the species were not seen before 2025
sp_2025 = sp_2025 |>
  group_by(scientific_name) |>
  summarise("min_year" = min(min_year),
            "max_year" = max(max_year))

# keep species whose minimum year is in 2025
sp_2025 = sp_2025 |>
  filter(min_year == 2025)

# save the species list
write.csv(sp_2025, "outputs/found_species_in_2025.csv", row.names = FALSE)

# Plots and summaries ----------------------------------------------------------
  
# how many missing species were found? (with corrected species list)

# double check by counting obs per species
correct.sp = inat |>
  filter(scientific_name %in% sp_2025$scientific_name) |>
  group_by(scientific_name, year) |>
  summarise("n" = n()) |>
  summarise("min_year" = min(year, na.rm = T)) |>
  filter(min_year == 2025) |> collect()

# record total number of new-to-inat species
n.sp = unique(correct.sp$scientific_name) |> length()

# summarise
new.groups = inat |>
  filter(scientific_name %in% correct.sp$scientific_name) |>
  group_by(iconic_taxon_name) |>
  distinct(scientific_name) |>
  summarise("n" = n()) |>
  collect()
new.groups = filter(new.groups, iconic_taxon_name != "NA")
# reorder for plotting
new.groups$iconic_taxon_name = factor(new.groups$iconic_taxon_name,
                                      levels = new.groups$iconic_taxon_name[order(new.groups$n)])
total.sp = sum(new.groups$n)
ggplot(data = new.groups) +
  geom_bar(aes(
    y = iconic_taxon_name,
    x = n,
    fill = iconic_taxon_name
  ), stat = "identity", position = "dodge") +
  geom_text(aes(
    y = iconic_taxon_name,
    label = n,
    x = n,
    col = iconic_taxon_name
  ), size = 5,  hjust = -.2,fontface = "bold") +
  colorspace::scale_fill_discrete_qualitative() +
  colorspace::scale_color_discrete_qualitative() +
  labs(y = "", 
       x = "Species",
       title = paste0(total.sp, " new species were recorded in 2025")) +
  hrbrthemes::theme_ipsum_rc(base_size = 14,
                             axis_title_size = 14,
                             axis_title_face = "bold") +
  coord_cartesian(xlim = c(0, max(new.groups$n+12))) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank()) 
ggsave("figures/missing-species-pergroup_researchgrade.png", width = 7.61, height = 4.71)



# which ones are the most observed? --------------------------------------------

# summarise
obs.newsp = inat |>
  filter(scientific_name %in% correct.sp$scientific_name,
         iconic_taxon_name != "NA") |>
  group_by(iconic_taxon_name, scientific_name) |>
  summarise("n" = n()) |>
  collect()
write.csv(obs.newsp, "outputs/missing-species/newsp_nobs_2025.csv")


# map the new observations -----------------------------------------------------

newsp.pts = inat |>
  filter(scientific_name %in% correct.sp$scientific_name,
         iconic_taxon_name != "NA",
         quality_grade == "research") |>
  select(c(longitude, latitude, scientific_name, iconic_taxon_name, user_login)) |>
  collect() |>
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = "epsg:4326")
newsp.pts = st_transform(newsp.pts, "EPSG:3347")
# reorder for plotting
newsp.pts$iconic_taxon_name = factor(newsp.pts$iconic_taxon_name,
                                      levels = new.groups$iconic_taxon_name[order(new.groups$n)])
# canada polygon
canada = st_read("data/base-layers/canada-polygon/canada.outline.shp")
canada = st_transform(canada, "EPSG:3347")

# plot!
ggplot() +
  geom_sf(data = canada, col = "grey20", fill = "grey20") +
  geom_sf(data = newsp.pts,
          aes(col = iconic_taxon_name)) +
  colorspace::scale_color_discrete_qualitative(name = "") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) 
ggsave("figures/map_newspecies2025.png", width = 10, height = 8)

mapview::mapview(newsp.pts, zcol = "iconic_taxon_name")



# who saw these first (research-grade) ? ---------------------------------------

finders.spnames = newsp.pts |>
  group_by(user_login, scientific_name) |>
  summarise("n" = n())

finders = newsp.pts |>
  group_by(user_login) |>
  distinct(scientific_name) |>
  summarise("n_sp" = n())

finders$name = NA
names.finders = list()
for(n in 1:nrow(finders)){
  names.finders[[n]] = rinat::get_inat_user_stats(uid = finders$user_login[n])
}
finders$name = names.finders

# import qcbs logins
qcbs <- readRDS("~/Documents/GitHub/storymap-qcbs/data/championteams_userlogins.rds")
qcbs.finders = finders$user_login[which(finders$user_login %in% qcbs)]

qc = st_read("data/base-layers/prov_territo/prov_territo.shp") |>
  filter(PROV_TERRI == "QC")

# plot!
ggplot() +
  geom_sf(data = qc, col = "grey20", fill = "grey20") +
  geom_sf(data = filter(newsp.pts, user_login %in% qcbs.finders),
          aes(col = iconic_taxon_name)) +
  colorspace::scale_color_discrete_qualitative(name = "") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +theme_bw()

mapview::mapview(filter(newsp.pts, user_login %in% qcbs.finders), 
                 zcol = "scientific_name", 
                 layer = "Species", 
                 cex = 10,
                 map.types = "Esri.WorldImagery",
                 alpha.regions = 1,
                 col.regions = PNWColors::pnw_palette("Bay", 10)[3:10])

new.groups.qcbs = filter(newsp.pts, user_login %in% qcbs.finders) |>
  group_by(iconic_taxon_name, scientific_name) |>
  summarise("n" = n())

total.sp.qcbs = sum(new.groups.qcbs$n)
new.groups.qcbs 
new.groups.qcbs$iconic_taxon_name = factor(new.groups.qcbs$iconic_taxon_name,
                                      levels = new.groups$iconic_taxon_name[order(new.groups$n)])

ggplot(data = new.groups.qcbs) +
  geom_bar(aes(
    y = iconic_taxon_name,
    x = n,
    fill = iconic_taxon_name
  ), stat = "identity", position = "dodge") +
  geom_text(aes(
    y = iconic_taxon_name,
    label = n,
    x = n,
    col = iconic_taxon_name
  ), size = 5,  hjust = -.2,fontface = "bold") +
  colorspace::scale_fill_discrete_qualitative() +
  colorspace::scale_color_discrete_qualitative() +
  labs(y = "", 
       x = "Species",
       title = paste0(total.sp, " species")) +
  hrbrthemes::theme_ipsum_rc(base_size = 14,
                             axis_title_size = 14,
                             axis_title_face = "bold") +
  coord_cartesian(xlim = c(0, max(new.groups.qcbs$n+5))) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank()) 
g


# # make into one big df
# df = list(animals, insects, plants, fungi) |>
#   lapply(select, c(user_login, 
#                    quality_grade, 
#                    latitude, longitude,
#                    place_state_name,
#                    iconic_taxon_name,
#                    scientific_name, 
#                    common_name,
#                    taxon_subspecies_name)) |>
#   bind_rows()

# Some species in these original projects do not belong (they had many observations)
# before BTG... might have been included because of taxonomic glitches

# Species to remove ------------------------------------------------------------

# load parquet file of pre-BTG to double check
inat_pq <- arrow::open_dataset("~/McGill University/Laura's Lab_Group - BioBlitz/data/raw/biodiversity-data/inat-canada/iNat_non_sensitive_data_Jan2025.parquet")

# Count number of observations per species, per taxa group
n_obs <- inat_pq |> 
  dplyr::filter(scientific_name %in% df$scientific_name) |> 
  group_by(iconic_taxon_name, scientific_name) |>
  summarize(total_obs = n()) |> 
  # load the query into our R session
  collect()

# filter to species with research-grade obs
n_obs.research = inat_pq |> 
  dplyr::filter(scientific_name %in% df$scientific_name,
                quality_grade == "research") |> 
  group_by(iconic_taxon_name, scientific_name) |>
  summarize(total_obs = n()) |> 
  # load the query into our R session
  collect()

# species that have previously been observed in iNat
to.remove = n_obs$scientific_name

# filter out those species
df = df[-which(df$scientific_name %in% to.remove),]

# how many missing species were found?
new.sp = df |>
  filter(quality_grade == "research") |>
  group_by(iconic_taxon_name, scientific_name) |>
  summarise("n" = n())

# manually remove ones that had obs in inat (suspiciously high ones)
df = df[-which(df$scientific_name== "Boechera polyantha"),]





# check how many were not in GBIF before 2025
library(rgbif)
gbif = lapply(correct.sp$scientific_name, 
              function(x) tryCatch(occ_count(scientificName = x, country = "CA", facet="year"), error=function(e) NULL))
names(gbif) = correct.sp$scientific_name
gbif_df = bind_rows(gbif, .id = "scientific_name")

# which ones are from 2025 only?
# first: these are ones that are not even in gbif
not_in_gbif = lapply(gbif, length) |> unlist() 
not_in_gbif = correct.sp$scientific_name[which(not_in_gbif == 0)]
