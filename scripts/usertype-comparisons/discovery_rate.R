# Script to count number of new species per year to see how many are added

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(sf)
library(taxize)
library(ape)

# set common theme for all maps
theme_set( hrbrthemes::theme_ipsum_rc(base_size = 14, axis_title_size = 16))

# load parquet file
inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")

# find species that were not found before 2025 ---------------------------------

missing_sp = inat |>
  filter(quality_grade == "research") |>
  select(iconic_taxon_name, scientific_name, year) |>
  group_by(scientific_name) |>
  distinct(year) |>
  na.omit() |>
  collect()
missing_sp$year = as.numeric(missing_sp$year)

sp_ids = inat |>
  select(taxon_id, iconic_taxon_name, scientific_name) |>
  distinct(taxon_id, .keep_all = T) |>
  collect()

# observations per year
obs_year = inat |>
  filter(quality_grade == "research") |>
  group_by(year) |>
  summarise("obsperyear" = n()) |> collect()
obs_year$year = as.numeric(obs_year$year)
obs_year = na.omit(obs_year)

# observations per year
obs_year.groups = inat |>
  filter(quality_grade == "research") |>
  group_by(year, iconic_taxon_name) |>
  summarise("obsperyear" = n()) |> 
  collect()
obs_year.groups$year = as.numeric(obs_year.groups$year)
obs_year.groups = na.omit(obs_year.groups)

# when was each species first and last found?
sp_minyear = missing_sp |>
  group_by(scientific_name) |>
  summarise("min_year" = min(year),
            "max_year" = max(year))

# filter to species names (genus + species)
temp = stringr::str_count(sp_minyear$scientific_name, "\\w+")
# find species with nothing in the species column
to_remove = which(temp != 2)
to_remove = c(to_remove, grep("×", sp_minyear$scientific_name)) # note: this is not a normal "x"

# remove species that have > 2 words, or that include the hybrid symbol
sp_minyear <- sp_minyear[-to_remove,]
rm(temp)
sp_minyear = filter(sp_minyear, min_year != 203) |> filter(!is.na(min_year))

# match with species ids to make the tree visualisation later
sp_fortree = left_join(sp_minyear, sp_ids)


# count number of species discovered each year
sp_yearly = sp_minyear |>
  group_by(min_year) |>
  summarise("n" = n())
sp_yearly$cumsum = cumsum(sp_yearly$n)

# plot
ggplot() +
  geom_line(data = sp_yearly,
            aes(x = min_year, y = cumsum)) +
  labs(y = "Species recorded", x = "Year") +
  coord_cartesian(xlim = c(1960,2025))
ggsave("figures/timeline_recordedspecies_all_alltime.png", width = 9, height = 6)

ggplot() +
  geom_line(data = sp_yearly,
            aes(x = min_year, y = cumsum)) +
  labs(y = "Species recorded", x = "Year") +
  coord_cartesian(xlim = c(2008,2025))
ggsave("figures/timeline_recordedspecies_all.png", width = 9, height = 6)

ggplot() +
  geom_line(data = sp_yearly,
            aes(x = min_year, y = n)) +
  labs(y = "New species recorded", x = "Year") +
  coord_cartesian(xlim = c(2008,2025))
ggsave("figures/timeline_recordednewspecies_all.png", width = 9, height = 6)

## per group

# count number of species discovered each year
sp_yearly.group = sp_fortree |>
  group_by(min_year, iconic_taxon_name) |>
  summarise("n" = n())
sp_yearly.group$cumsum = cumsum(sp_yearly.group$n)

# plot
ggplot() +
  geom_line(data = na.omit(sp_yearly.group),
            aes(x = min_year, y = cumsum, col = iconic_taxon_name)) +
  labs(x = "Year", y = "Recorded species", col = "")+
  coord_cartesian(xlim = c(2008,2025))
ggsave("figures/timeline_recordedspecies_groups.png", width = 9, height = 6)


## discovery rate - make this relative to the number of observations!

discovery = left_join(sp_yearly, obs_year, by = c("min_year" = "year"))
discovery$rate = discovery$n/discovery$obsperyear
discovery$obsneeded = discovery$obsperyear/discovery$n

# plot
ggplot() +
  geom_line(data = filter(discovery, min_year > 2009),
            aes(x = min_year, y = rate)) +
  labs(x = "Year", y = "Discovery rate (species/observation)", col = "") 
ggplot() +
  geom_line(data = filter(discovery, min_year > 2009),
            aes(x = min_year, y = obsneeded)) +
  labs(x = "Year", y = "Observations per new species", title = "It's becoming harder to record new species", col = "") 
ggsave("figures/timeline_effort_for_discovery.png", width = 9, height = 6)

# group-wise
discovery.groups = left_join(sp_yearly.group, obs_year.groups, 
                             by = c("min_year" = "year", "iconic_taxon_name"))
discovery.groups$rate = discovery.groups$n/discovery.groups$obsperyear
discovery.groups$obsneeded = discovery.groups$obsperyear/discovery.groups$n
discovery.groups = na.omit(discovery.groups)

# plot
ggplot() +
  geom_line(data = filter(discovery.groups),
            aes(x = min_year, y = obsneeded, col = iconic_taxon_name),
            linewidth = .5) +
  labs(x = "Year", y = "Observations / new species record", col = "") +
  geom_line(data = filter(discovery, min_year > 2007),
            aes(x = min_year, y = obsneeded), linewidth = 1) +
  scale_y_sqrt(breaks = c(100, 500, 1000, 2500, 5000, 10000, 20000, 30000,50000)) +
  ggrepel::geom_text_repel(data = filter(discovery.groups, min_year == 2025),
                             aes(x = 2025, y = obsneeded, 
                                 label = iconic_taxon_name, 
                                 col = iconic_taxon_name),
                             size = 3.5, hjust = -.25, min.segment.length = 1,
                           direction = "y", family = "Roboto Condensed") +
  colorspace::scale_color_discrete_qualitative() +
  coord_cartesian(xlim = c(2010, 2027)) +
  theme(legend.position = "none")
ggsave("figures/discoveryrate_obspernewsp_canada.png", width = 8.52, height = 7.77)

## smoothed version ----
library(broom)
label_df <- discovery.groups |> 
  nest(data = -iconic_taxon_name) |> 
  mutate(model = map(data, ~loess(obsneeded ~ min_year, .x)),
         augmented = map(model, broom::augment),
         fitted = map(augmented, ".fitted") |> map_dbl(last),
         year = map(data, "min_year") |> map_int(last)) |> 
  select(iconic_taxon_name, fitted, year)



# plot
ggplot() +
  geom_smooth(data = discovery.groups,
            aes(x = min_year, y = obsneeded, 
                col = iconic_taxon_name, 
                fill = iconic_taxon_name),
            linewidth = .5, alpha = .1) +
  labs(x = "Year", y = "Observations / new species record", col = "") +
  scale_y_sqrt(breaks = c(100, 500, 1000, 2500, 5000, 10000, 20000, 30000)) +
  ggrepel::geom_text_repel(data = label_df,
                           aes(x = year, y = fitted, 
                               label = iconic_taxon_name, 
                               col = iconic_taxon_name),
                           size = 4, min.segment.length = .5,
                           direction = "y", hjust = -2, family = "Roboto Condensed") +
  colorspace::scale_color_discrete_qualitative() +
  coord_cartesian(xlim = c(2010, 2030), ylim = c(0,30000)) +
  theme(legend.position = "none")
ggsave("figures/discoveryrate_obspernewsp_smoothed_canada.png", width = 8.8, height = 7.75)

# add zero for new amphibians in 2025
complete_df = discovery.groups[which(discovery.groups$iconic_taxon_name == "Amphibia")[1],]
complete_df$min_year = 2025
complete_df$n = 0 
complete_df[,5:7] = NA
discovery.groups = bind_rows(discovery.groups, complete_df)
ggplot(data = discovery.groups,
       aes(x = min_year, y = n, 
           col = iconic_taxon_name, 
           fill = iconic_taxon_name)) +
  # moving window average: 5 years
  tidyquant::geom_ma(n = 5, lty = 1, size = .7) +
  ggrepel::geom_text_repel(
    data = filter(discovery.groups, min_year %in% 2020:2025) |>
      group_by(iconic_taxon_name) |>
      summarise("label_y" = mean(n, na.rm = T)),
    aes(x = 2025, y = label_y, 
                               label = iconic_taxon_name, 
                               col = iconic_taxon_name),
                           size = 4, min.segment.length = 1,
                           direction = "y", hjust = -.5, family = "Roboto Condensed") +
  #labs(x = "Year", y = "Observations / new species record", col = "") +
  scale_y_sqrt(breaks = c(100, 250, 500, 1000)) +
  colorspace::scale_color_discrete_qualitative() +
  coord_cartesian(xlim = c(2008, 2030)) +
  labs(y = "Newly-recorded species", x = "") +
  theme(legend.position = "none")
ggsave("figures/discoveryrate_obspernewsp_mw5avg_canada.png", width = 7.61, height = 6.6)


discovery.groups$year = NA
discovery.groups$year[discovery.groups$min_year < 2015] <- "2009-2014"
discovery.groups$year[discovery.groups$min_year >= 2015] <- "2015-2019"
discovery.groups$year[discovery.groups$min_year >= 2020] <- "2020-2025"
write.csv(discovery.groups, "outputs/n_observations/discovery-rate.csv")

toplot = discovery.groups |>
  group_by(year, iconic_taxon_name) |>
  summarise("mean_obs" = mean(obsneeded, na.rm = T),
            "max_obs" = max(obsneeded, na.rm = T),
            "mean_rate" = mean(rate, na.rm = T),
            "max_rate" = max(rate, na.rm = T))
temp = toplot |> filter(year == "2020-2025") 


toplot$iconic_taxon_name = factor(toplot$iconic_taxon_name,
                                  levels = temp$iconic_taxon_name[order(temp$mean_obs)])
ggplot(data = toplot) +
  geom_line(aes(y = iconic_taxon_name, x = mean_obs),
             size = .3) +
  geom_point(aes(y = iconic_taxon_name, x = mean_obs, fill = year),
             size = 3, pch = 21) +
  scale_x_sqrt(breaks = c(100, 500, 1000, 2500, 5000, 10000, 20000)) +
  colorspace::scale_fill_discrete_divergingx(palette = "Zissou 1",
                                              name = "Period") +
  theme(legend.position = "top") +
  labs(x = "Observations before a new species record", y = "")
ggsave("figures/discoveryrate_obspernewsp_canada_pergroup_periods.png", width = 9.57, height = 6)

ggplot() +
  geom_line(data = filter(discovery.groups, min_year > 2007),
            aes(x = min_year, y = obsneeded, col = iconic_taxon_name),
            linewidth = .5) +
  labs(x = "Year", y = "Observations / new species record", col = "") +
  geom_line(data = filter(discovery, min_year > 2007),
            aes(x = min_year, y = obsneeded), linewidth = 1) +
  scale_y_sqrt(breaks = c(100, 500, 1000, 2500, 5000, 10000, 20000, 30000)) +
  ggrepel::geom_text_repel(data = filter(discovery.groups, min_year == 2025),
                           aes(x = 2025, y = obsneeded, 
                               label = iconic_taxon_name, 
                               col = iconic_taxon_name),
                           size = 3.5, hjust = -.25, min.segment.length = 1,
                           direction = "y", family = "Roboto Condensed") +
  colorspace::scale_color_discrete_qualitative() +
  coord_cartesian(xlim = c(2008, 2027)) +
  theme(legend.position = "none", panel.grid.major = element_line(linewidth = .1))
ggsave("figures/discoveryrate_obspernewsp_canada.png", width = 8.52, height = 7.77)


ggplot() +
  geom_line(data = filter(discovery.groups, min_year > 2007),
            aes(x = min_year, y = rate, col = iconic_taxon_name),
            linewidth = .5) +
  labs(x = "Year", y = "Discovery rate", col = "") +
  geom_line(data = filter(discovery, min_year > 2007),
            aes(x = min_year, y = rate), linewidth = 1) +
  scale_y_sqrt() +
  ggrepel::geom_text_repel(data = filter(discovery.groups, min_year == 2008),
                           aes(x = 2008, y = rate, 
                               label = iconic_taxon_name, 
                               col = iconic_taxon_name),
                           size = 3.5, min.segment.length = 1,
                           direction = "y", hjust = 1.5, family = "Roboto Condensed") +
  colorspace::scale_color_discrete_qualitative() +
  coord_cartesian(xlim = c(2005, 2025), ylim = c(0, 0.2)) +
  theme(legend.position = "none")
ggsave("figures/discoveryrate_canada.png", width = 8.52, height = 7.77)



## research to needs id ratio through time -------------------------------------

df = inat |>
  group_by(quality_grade, scientific_name, year, iconic_taxon_name) |>
  summarise("n" = n()) |>
  collect() |> na.omit()
df$year = as.numeric(df$year)

df.grade = df |>
  group_by(quality_grade, year) |>
  summarise("n_obs" = sum(n),
            "n_sp" = n())

df.research = df |>
  filter(quality_grade == "research") |>
  group_by(year) |>
  summarise("n_obs" = sum(n, na.rm = T),
            "n_sp" = n())
df.needsid = df |>
  filter(quality_grade == "needs_id") |>
  group_by(year) |>
  summarise("n_obs" = sum(n, na.rm = T),
            "n_sp" = n())

df.ratio.all = left_join(df.research, df.needsid, 
                         by = "year", 
                         suffix = c(".research", ".needsid"))
df.ratio.all$prop = df.ratio.all$n_obs.needsid/(df.ratio.all$n_obs.needsid+df.ratio.all$n_obs.research)
ggplot(data = df.ratio.all) +
  geom_line(aes(x = year, y = prop)) +
  coord_cartesian(xlim = c(2009, 2025), ylim = c(0,1)) +
  labs(y = "Proportion of species not yet IDed", y = "Year")
ggsave("figures/identification-gap/timeline_needsIDratio_allinone.png", width = 9, height = 6)


# with taxon groups

df.grade = df |>
  filter(quality_grade != "casual") |>
  group_by(quality_grade, year, iconic_taxon_name) |>
  summarise("n_obs" = sum(n),
            "n_sp" = n())

df.research = df.grade |>
  filter(quality_grade == "research") |>
  group_by(year, iconic_taxon_name) |>
  select(-c(quality_grade))
df.needsid = df.grade |>
  filter(quality_grade == "needs_id") |>
  group_by(year, iconic_taxon_name)|>
  select(-c(quality_grade))

df.ratio = left_join(df.research, 
                     df.needsid, 
                     by = c("year", "iconic_taxon_name"), 
                     suffix = c(".research", ".needsid"))
df.ratio$prop = df.ratio$n_obs.needsid/(df.ratio$n_obs.needsid+df.ratio$n_obs.research)
df.ratio$year = as.numeric(df.ratio$year)

# plot the proportion of needs id data through time ----------------------------

ggplot() +
  geom_line(data = df.ratio.all,
            aes(x = year, y = prop), linewidth = 1) +
  geom_line(data = df.ratio,
            aes(x = year, y = prop, col = iconic_taxon_name)) +
  geom_text(data = filter(df.ratio, year == 2025),
            aes(x = 2025, y = prop, 
                label = iconic_taxon_name, 
                col = iconic_taxon_name),
            hjust = -.1) +
  colorspace::scale_color_discrete_qualitative() +
  coord_cartesian(xlim = c(2009, 2030), ylim = c(0,1)) +
  labs(y = "Proportion of species not yet IDed", y = "Year",
       col = "") +
  theme(legend.position = "none")
ggsave("figures/identification-gap/timeline_needsIDratio_groups.png",
       width = 7, height = 10)


