# Script to calculate which cells now have enough observations to be past a given
# threshold

# load libraries
library(ggplot2)
library(dplyr)
library(terra)
library(gpkg)
library(plotly)

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# load range coverage function
source("scripts/range-coverage/calc_range_coverage_change.R")

# list of taxa groups
taxagroups = list.files("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/")
inatgroups = c("Amphibia", "Aves", "Insecta", "Mammalia", "Plantae", "Reptilia")


# function to calculate and map the change in cells that meet a threshold ======
coverage_threshold = function(before, after, threshold = 1, aggregation_factor = 1){
  
  after_min = after
  after_min[after >= threshold] <- 1
  after_min[after < threshold] <- 0
  
  before_min = before
  before_min[before >= threshold] <- 1
  before_min[before < threshold] <- 0
  
  change = after_min-before_min
  
  rasters = c(before_min, after_min, change)
  names(rasters) = c("before", "after", "upgraded")
  return(rasters)
}


## Calculate coverage change with a thresholded # of obs per cell ==============

# List range coverage rasters
spp = lapply(as.list(paste0("~/McGill University/Laura's Lab_Group - range-coverage/",taxagroups,"/")), list.files)
names(spp) = taxagroups
df2 = spp |> lapply(as.data.frame) |> bind_rows(.id = "Taxa") |> na.omit()
colnames(df2)[2] = "species"
df2$n_upgraded_min01 = NA
df2$n_upgraded_min05 = NA
df2$n_upgraded_min10 = NA
df2$before_km = NA
df2$after_km = NA
df2$change_km = NA
df2$range_km = NA
df2$n_cells = NA

# Run change in range coverage for all species!
for(t in c(1:4,6)){
  
  # Load species names that have previously been included
  species = dplyr::filter(df2, Taxa == taxagroups[t]) |>
    select(species) |>
    as.matrix() |> as.vector()
  
  # Prepare to load range coverage results
  filepath = paste0("~/McGill University/Laura's Lab_Group - range-coverage/",taxagroups[t],"/")
  
  for(i in 1:length(species)){
    
    cvg = terra::rast(paste0(filepath, species[i]))
    
    tryCatch({  
      index = which(df2$Taxa == taxagroups[t] & df2$species == species[i])
      
      # Minimum 1 observation
      temp = coverage_threshold(cvg$before,
                                cvg$after,
                                threshold = 1)
      df2$n_upgraded_min01[index] = sum(values(temp), na.rm = T)
      terra::writeRaster(temp, 
                         paste0("~/McGill University/Laura's Lab_Group - range-coverage/upgraded_cells/",taxagroups[t],"/min_1obs/", species[i]), overwrite=T)
      
      # Minimum 5 observations
      temp = coverage_threshold(cvg$before,
                                cvg$after,
                                threshold = 5)
      df2$n_upgraded_min05[index] = sum(values(temp), na.rm = T)
      terra::writeRaster(temp, 
                         paste0("~/McGill University/Laura's Lab_Group - range-coverage/upgraded_cells/",taxagroups[t],"/min_5obs/", species[i]), overwrite=T)
      
      # Minimum 10 observations
      temp = coverage_threshold(cvg$before,
                                cvg$after,
                                threshold = 10)
      df2$n_upgraded_min10[index] = sum(values(temp), na.rm = T)
      terra::writeRaster(temp, 
                         paste0("~/McGill University/Laura's Lab_Group - range-coverage/upgraded_cells/",taxagroups[t],"/min_10obs/", species[i]), overwrite=T)
      
      ## save number of cells in the range -------------------------------------
      
      n_cells = terra::values(cvg$before, na.rm = TRUE) |> length()
      
      ## save coverage as an area ----------------------------------------------
      
      # change zeros to NAs
      cvg_mod = cvg
      cvg_mod$before[cvg_mod$before < 1] <- NA
      cvg_mod$after[cvg_mod$after < 1] <- NA
      cvg_mod$change[cvg_mod$change < 1] <- NA
      
      # calculate summed area of non-NA cells
      cvg_area = terra::expanse(cvg_mod, unit = "km")[,2]
      range_area = terra::expanse(cvg$before, unit = "km")[2]
      
      # save in the data frame
      df2$before_km[index] = cvg_area[1]
      df2$after_km[index] = cvg_area[2]
      df2$change_km[index] = cvg_area[3]
      df2$range_km[index] = range_area
      df2$n_cells[index] = n_cells
      
    }, error = function(e) {
      message("An error occurred: ", e$message)
      return(NA) # Return NA if an error occurs
    })
  }
} 
df2$change_km = df2$after_km - df2$before_km
df2$range_km = unlist(df2$range_km)
# save the results
saveRDS(df2, "outputs/range-coverage/summary-results/upgraded_coverage_thresholds.rds")

df2_toplot = df2
df2_toplot$species = gsub(".tif", "", df2_toplot$species)
(p = ggplot(data = df2_toplot, 
            aes(x = Taxa, group = species,
                col = Taxa,
                y = n_upgraded_min05)) +
  geom_jitter(width = .4, size = 2, alpha = .7) +
  coord_flip() +
  labs(y = "Number of upgraded cells", 
       title = "Gains in range coverage",
       ) +
  theme(legend.position = "none")
  )
ggsave("figures/upgraded_cells_min05.png", width = 6.07, height = 4.18)


ggplotly(p, tooltip = "species")
plotly::gg2list(p, tooltip = "species") |> 
  plotly::as_widget() |>
  htmlwidgets::saveWidget(file = "figures/upgraded_cells_min05.html")


# Plot range coverage gains by PROPORTION (area) ===============================

df3 = df2 |> 
  dplyr::filter(range_km > 0) |> # Delete rows with 0 range area
  select(c(Taxa, species, before_km:range_km)) 

# calculate a percentage of range covered in km2
df3$before_perc = 100*df3$before_km/df3$range_km
df3$after_perc = 100*df3$after_km/df3$range_km

# order species names by the change in range for nicer plotting later
df3 = df3 |>
  group_by(Taxa) |>
  mutate("species_fct" = factor(species, 
                                level = species[order(change_km, decreasing = TRUE)]))
df3$species = gsub(".tif", "", df3$species)

# reset theme for the plotly part
theme_set(theme_classic())

(p.perc = ggplot(data = df3) +
  geom_jitter(aes(x = 100*change_km/range_km,
                  y = Taxa,
                  group = species,
                 col = Taxa),
             size = 3, alpha = .7) +
  colorspace::scale_color_discrete_qualitative() +
  labs(x = "Gained range coverage (% of the total range)", 
       y = "",
       title = "Gains in range coverage"
  ) +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14) +
    theme(legend.position = "none",
          axis.title = element_text(size = 16)))

plotly::gg2list(p.perc, tooltip = "species") |> 
  plotly::as_widget() |>
  htmlwidgets::saveWidget(file = "figures/gained_coverage_percentage.html")


# Plot range coverage gains by AREA ============================================

(p.area = ggplot(data = df3) +
    geom_jitter(aes(x = change_km,
                    y = Taxa,
                    group = species,
                    col = Taxa),
                size = 3, alpha = .7) +
    colorspace::scale_color_discrete_qualitative() +
    labs(x = "Gained range coverage (km²)", 
         y = "",
         title = "Gains in range coverage"
    ) +
    hrbrthemes::theme_ipsum_rc(base_size = 14,
                               axis_title_size = 14) +
    theme(legend.position = "none",
          axis.title = element_text(size = 16)))

plotly::gg2list(p.area, tooltip = "species") |> 
  plotly::as_widget() |>
  htmlwidgets::saveWidget(file = "figures/gained_coverage_area.html")


# Top 10 range coverage gains in each taxa =====================================

# select the top 10 species in each taxa
top10 = df3 |>
  group_by(Taxa) |>
  group_split(.keep = TRUE)
indices = lapply(top10, function(x) order(x$change_km, decreasing = TRUE)[1:10])
for(i in 1:length(top10)){
  top10[[i]] = top10[[i]][indices[[i]],]
}

# bind rows and order species names by final percentage of coverage
top10 = top10 |> 
  bind_rows() |>
  group_by(Taxa) |>
  mutate("species_fct" = factor(species, levels = species[order(after_perc)]))

# convert to long format to make groupings based on before or after BTG
top10 = top10 |> tidyr::pivot_longer(cols = c(before_perc, after_perc),
                      names_to = "when",
                      values_to = "coverage_perc")
top10$when = factor(top10$when,
                    levels = c("before_perc", "after_perc"))

# make the plot
(p.anim = ggplot(data = top10) +
    geom_bar(aes(y = species_fct, 
                    x = coverage_perc,
                 alpha = when,
                    fill = Taxa),
             stat = "identity") +
    colorspace::scale_fill_discrete_qualitative() +
    scale_alpha_discrete(range = c(1, .7)) +
    labs(x = "Range coverage (%)", 
         y = "",
         title = "Top 10: Gained range coverage",
         subtitle = "Proportion of species' ranges covered by at least 1 observation per 10km cell before (pale) and after (dark) Blitz the Gap."
    ) +
    hrbrthemes::theme_ipsum_es(base_size = 14,
                               axis_title_size = 14) +
    theme(legend.position = "none",
          axis.text.y = element_text(face = "italic"))) +
  facet_wrap(~Taxa, scales = "free_y", ncol = 3)
ggsave("figures/gained_coverage_barplot.png", width = 12, height = 10)