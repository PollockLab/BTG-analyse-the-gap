# Script to count number of species per taxa that were previously missing in iNat Canada

# libraries
library(rgbif)
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(taxize)
library(rinat)

# load manually downloaded project data
animals = read.csv("data/heavy/missing-species/missing-animals/observations-619233.csv")
insects = read.csv("data/heavy/missing-species/missing-insects.csv/observations-619217.csv")
plants = read.csv("data/heavy/missing-species/missing-plants/observations-619220.csv")
fungi = read.csv("data/heavy/missing-species/missing-fungi/observations-619222.csv")

# make into one big df
df = list(animals, insects, plants, fungi) |>
  lapply(select, c(user_login, 
                   quality_grade, 
                   latitude, longitude,
                   place_state_name,
                   iconic_taxon_name,
                   scientific_name, 
                   common_name,
                   taxon_subspecies_name)) |>
  bind_rows()

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

# how many missing species were found? (with corrected species list)
new.sp = df |>
  filter(quality_grade == "research") |>
  group_by(iconic_taxon_name, scientific_name) |>
  summarise("n" = n())

# how many species per group?
new.groups = new.sp |>
  group_by(iconic_taxon_name) |>
  summarise("n" = n())

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
       x = "Species found",
       title = paste0("Missing species found this summer: ", total.sp, " ")) +
  hrbrthemes::theme_ipsum_rc(base_size = 14,
                             axis_title_size = 14,
                             axis_title_face = "bold") +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank()) +
  coord_cartesian(xlim = c(0,100))
ggsave("figures/missing-species-pergroup_researchgrade.png", width = 7.61, height = 4.71)


# who saw these first (research-grade) ? ---------------------------------------

finders = df |>
  filter(quality_grade == "research") |>
  group_by(user_login) |>
  distinct(scientific_name) |>
  summarise("n_sp" = n())

finders$name = NA
names.finders = list()
for(n in 1:nrow(finders)){
  names.finders[[n]] = rinat::get_inat_user_stats(uid = finders$user_login[n])
}
finders$name = names.finders

finders.rg = df |>
  filter(quality_grade == "research") |>
  group_by(user_login) |>
  distinct(scientific_name) |>
  summarise("n_sp" = n())

finders2 = df |>
  group_by(user_login, scientific_name, quality_grade) |>
  summarise("n_obs" = n())

temp = rinat::get_inat_user_stats(uid = "kalvinchan")
temp$most_observations$user$name

