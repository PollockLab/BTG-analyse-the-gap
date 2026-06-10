# Script to count number of species per taxa that were previously missing in iNat Canada

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(sf)

# load parquet file
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# custom list by Nathan and Brian of species that had obs prior to 2025 despite this filtering
sp_to_remove_checked = c("Drepanocladus polygamus",
                         "Castilleja minor",
                         "Cicindela fulgida",
                         "Cicindela parowana",
                         "Crataegus laurentiana",
                         "Crataegus chrysocarpa",
                         "Eubosmina longispina",
                         "Bosmina coregoni",
                         "Hyssopus thymus",
                         "Hyssopus officinalis",
                         "Numenius phaeopus",
                         "Poa leptocoma",
                         "Poa paucispicula",
                         "Podosphaera clandestina",
                         "Podosphaera aucupariae",
                         "Veronica agrestis",
                         "Veronica hederifolia",
                         "Veronica alpina",
                         "Veronica wormskjoldii")

# find species that were not found before 2025 ---------------------------------

missing_sp = inat |>
  filter(quality_grade == "research") |>
  group_by(iconic_taxon_name, scientific_name, year) |>
  distinct(year) |>
  na.omit() |>
  collect()
missing_sp$year = as.numeric(missing_sp$year)
# change Animalia to "Other Animalia"
missing_sp$iconic_taxon_name = gsub("Animalia", "Other Animalia", missing_sp$iconic_taxon_name)

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
to_remove = c(to_remove, grep("Mustela furo", sp_2025$scientific_name)) 
to_remove = c(to_remove, which(sp_2025$scientific_name %in% sp_to_remove_checked))

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

# read it back in
sp_2025 = read.csv("outputs/found_species_in_2025.csv")

## species found during Blitz the Gap (April to Oct 1) -------------------------

# load parquet for missing species challenge timeframe
btg.pq = read_parquet("data/heavy/BTG-data/inaturalist-canada-dec2025_APR1OCT1.parquet")
btg = btg.pq |>
  group_by(iconic_taxon_name, scientific_name) |>
  summarise("n" = n())
btg = collect(btg)
# change Animalia to "Other Animalia"
btg$iconic_taxon_name = gsub("Animalia", "Other Animalia", btg$iconic_taxon_name)
# filter for species in the Apr-Oct timeframe parquet
sp_btg = sp_2025[which(sp_2025$scientific_name %in% btg$scientific_name),]

# Plots and summaries ----------------------------------------------------------

# how many missing species were found? (with corrected species list)

# double check by counting obs per species
correct.sp = inat |>
  filter(scientific_name %in% sp_btg$scientific_name) |>
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
new.groups$iconic_taxon_name = gsub("Animalia", "Other Animalia", new.groups$iconic_taxon_name)
# reorder for plotting
new.groups$iconic_taxon_name = factor(new.groups$iconic_taxon_name,
                                      levels = new.groups$iconic_taxon_name[order(new.groups$n)])
total.sp = sum(new.groups$n)
(A = ggplot(data = new.groups) +
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
          panel.grid.major.y = element_blank()) )
ggsave("figures/missing-species-pergroup_researchgrade_APR1OCT1.png", width = 7.61, height = 4.71)



# which ones are the most observed? --------------------------------------------

# summarise
obs.newsp = btg.pq |>
  filter(scientific_name %in% correct.sp$scientific_name,
         iconic_taxon_name != "NA") |>
  group_by(iconic_taxon_name, scientific_name) |>
  summarise("n" = n()) |>
  collect()
obs.newsp$iconic_taxon_name = gsub("Animalia", "Other Animalia", obs.newsp$iconic_taxon_name)

write.csv(obs.newsp, "outputs/missing-species/newsp_nobs_btg_APR1OCT1.csv")

# average obs per species
mean(obs.newsp$n)
sd(obs.newsp$n)

# this accounts for how much of the BTG data?
obs.newsp.sum = sum(obs.newsp$n)
totalobs = btg.pq |>
  filter(iconic_taxon_name != "NA") |>
  group_by(iconic_taxon_name, scientific_name) |>
  summarise("n" = n()) |>
  summarise("nn" = sum(n)) |>
  collect()
totalobs$iconic_taxon_name = gsub("Animalia", "Other Animalia", totalobs$iconic_taxon_name)
obs.newsp.sum/sum(totalobs$nn)

# map the new observations -----------------------------------------------------

newsp.pts = btg.pq |>
  filter(scientific_name %in% correct.sp$scientific_name,
         iconic_taxon_name != "NA",
         quality_grade == "research") |>
  select(c(longitude, latitude, scientific_name, iconic_taxon_name, user_login)) |>
  collect() |>
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = "epsg:4326")
newsp.pts = st_transform(newsp.pts, "EPSG:3347")
# reorder for plotting
newsp.pts$iconic_taxon_name = gsub("Animalia", "Other Animalia", newsp.pts$iconic_taxon_name)
newsp.pts$iconic_taxon_name = factor(newsp.pts$iconic_taxon_name,
                                     levels = new.groups$iconic_taxon_name[order(new.groups$n)])
# canada polygon
canada = st_read("data/base-layers/canada-polygon/canada.outline.shp")
canada = st_transform(canada, "EPSG:3347")

# plot!
(B = ggplot() +
    geom_sf(data = canada, col = "grey90", fill = "grey90") +
    geom_sf(data = newsp.pts,
            aes(fill = iconic_taxon_name), 
            size = 3, alpha = .9, pch = 21) +
    colorspace::scale_fill_discrete_qualitative(name = "") +
    hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                               axis = FALSE, axis_text_size = 1,
                               ticks = FALSE,
                               base_size = 14)) 
ggsave("figures/map_newspeciesAPR1OCT12025.png", width = 10, height = 8)


# who saw these first (research-grade) ? ---------------------------------------

finders.spnames = newsp.pts |>
  group_by(user_login, scientific_name) |>
  summarise("n" = n())
saveRDS(finders.spnames, "outputs/missing-species/newsp_finders_sf.rds")

finders = newsp.pts |>
  group_by(user_login) |>
  distinct(scientific_name) |>
  summarise("n_sp" = n())
write.csv(finders, "outputs/missing-species/newsp_finders_nsp.csv")

## which ones are new Globally? ------------------------------------------------
# According to the CWF's list here: https://www.inaturalist.ca/projects/1-global-observation-in-canada-1-observation-globale-au-canada

obs.newsp = read.csv("outputs/missing-species/newsp_nobs_btg_APR1OCT1.csv")
global = read.csv("data/heavy/missing-species/1 global observation in Canada | 1 observation globale au Canada/observations-733976.csv")
global.2025 = global[grepl("2025", global$observed_on),]

new.globally = obs.newsp |> filter(scientific_name %in% global.2025$scientific_name)
new.globally = new.globally[order(new.globally$n, decreasing = FALSE),]
write.csv(new.globally, "outputs/missing-species/newsp_globally_nobs_APR1OCT1.csv")

# summarise
new.groups = btg.pq |>
  filter(scientific_name %in% new.globally$scientific_name) |>
  group_by(iconic_taxon_name) |>
  distinct(scientific_name) |>
  summarise("n" = n()) |>
  collect()
new.groups = filter(new.groups, iconic_taxon_name != "NA")
# reorder for plotting
new.groups$iconic_taxon_name = gsub("Animalia", "Other Animalia", new.groups$iconic_taxon_name)
new.groups$iconic_taxon_name = factor(new.groups$iconic_taxon_name,
                                      levels = new.groups$iconic_taxon_name[order(new.groups$n)])
total.sp = sum(new.groups$n)
(A = ggplot(data = new.groups) +
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
          panel.grid.major.y = element_blank()) )
ggsave("figures/missing-species-pergroup_GLOBAL_researchgrade_APR1OCT12025.png", width = 7.61, height = 4.71)
