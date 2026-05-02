# Script to compare 2025 observation density to prior density maps

# load libraries
library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(arrow)

# function to count cells that were not previously sampled but are now sampled
new_cells = function(map.stack, n_years = 18){

  # get last year of density map
  pa.2025 = map.stack[[n_years]]
  pa.2025[pa.2025 >= 1] <- 1
  pa.2025[is.na(pa.2025)] <- 0
  
  # sum all previous years 
  pa.pre2025 = sum(map.stack[[-n_years]], na.rm=T)
  pa.pre2025[pa.pre2025 >= 1] <- 1
  pa.pre2025[is.na(pa.pre2025)] <- 0
  
  # difference between 2025 filled cells and pre 2025 filled cells
  pa.diff = pa.2025 - pa.pre2025
  fullmap = sum(map.stack, na.rm=T)
  pa.diff[is.na(fullmap)] <- NA
  return(pa.diff)
}

# run!
files = list.files("outputs/densitymaps/")[grepl("yearly", list.files("outputs/densitymaps/"))]
maps = lapply(files, function(x) {rast(paste0("outputs/densitymaps/", x))})
cell_change = lapply(maps, new_cells)
# stack the rasters
cell_change = rast(cell_change)
names(cell_change) = gsub("densitymap_yearly_", "", files)
names(cell_change) = gsub(".tif", "", names(cell_change))
# save!
writeRaster(cell_change, "outputs/comparison-maps/density/changecells_alltaxa_2025_vs_20082024.tif", overwrite = TRUE)


## save a version that is only the new cells

# stack the rasters
cell_onlynew = list()
res_df = list()
for(i in 1:dim(cell_change)[3]){
  temp = cell_change[[i]]
  temp[temp!=1] <- NA   # change all non-1 values to NA
  cell_onlynew[[i]] <- temp
  # get cell area
  cell_area = cellSize(temp, unit = "km", mask = TRUE) |>
    values() |> na.omit() |> as.vector()
  # summarise the new cells
  res_df[[i]] = data.frame(
    "area_km_newcells" = sum(cell_area),
    "n_newcells" = length(cell_area)
  )
}
# stack the rasters
cell_onlynew = rast(cell_onlynew)
# save!
writeRaster(cell_onlynew, "outputs/comparison-maps/density/newcells_alltaxa_2025_vs_20082024.tif", overwrite=TRUE)

# save the summary results
names(res_df) = names(cell_change)
res_df = res_df |> bind_rows(.id = "taxa")
res_df$first_year = 2008
res_df$last_year = 2025
write.csv(res_df, "outputs/summaries/newcells_alltaxa_2025_vs_20082024.csv")

# map the new cells
cell_onlynew = rast("outputs/comparison-maps/density/newcells_alltaxa_2025_vs_20082024.tif")
cell_onlynew$alltaxa |> plot(col = "black")
cell_change = rast("outputs/comparison-maps/density/changecells_alltaxa_2025_vs_20082024.tif")
cell_change$alltaxa |> plot()

cell_change_pts = as.points(cell_change)
cell_new_pts = as.points(cell_onlynew$alltaxa)

# load canada polygon
canada_poly = terra::vect("data/base-layers/canada-polygon/canada.outline.shp")
canada_poly = project(canada_poly, crs(cell_change$alltaxa))

ggplot() +
  geom_sf(data = sf::st_as_sf(canada_poly), col = "black", fill = "black") +
  geom_spatraster(data = cell_change, aes(fill = alltaxa), interpolate = F) +
  scale_fill_viridis_c(option = "viridis",
                       name = "",
                       end = 1, begin = .1,
                       na.value = "transparent")  +
  theme_void()
ggsave("figures/cellchange_alltaxa_map.png", width = 12.5, height = 7.9)

# how many observations are there in new cells?
obsdens = maps[[2]] |> select(`2025`)
new = cell_onlynew |> select(alltaxa)
newarea = cell_onlynew |> select(alltaxa) |> cellSize(unit = "km") 
newarea[is.na(new)] <- NA
newarea = newarea |> values() |> sum(na.rm = T)

# observation density in new cells only
obsdens[is.na(new)] <- NA
newcells_dens = obsdens |> values() |> na.omit()
hist(log10(newcells_dens))
quantile(newcells_dens[newcells_dens>10])

# check out the maximum obs density cells to see the biggest gains
max.cell = obsdens
max.cell[obsdens<1000] <- NA
max.cell = as.points(max.cell)
mapview::mapview(max.cell)

# proportion of gained coverage
length(newcells_dens)
which(newcells_dens>1) |> length() # 20097/68605
which(newcells_dens>=10) |> length() # 2805/68605
which(newcells_dens>=30) |> length() # 720/68605
which(newcells_dens>=50) |> length() # 382/68605
which(newcells_dens>=100) |> length() # 143/68605
which(newcells_dens>=500) |> length() # 16/68605

# calculate area of Canada base map
base1k = terra::rast("data/base-layers/canada-basegrids/canada.base.1k.tiff") 
base1k_area = cellSize(base1k, unit = "km")
base1k_area[is.na(base1k)] <- NA
can_area = values(base1k_area) |> sum(na.rm = T)
newarea*100/can_area

# MAP NUMBER OF OBS AS WELL

obsdens.pts = as.points(obsdens)
ggplot() +
  geom_sf(data = sf::st_as_sf(canada_poly), col = "grey10", fill = "grey10") +
  geom_spatraster(data = filter(obsdens, `2025`< 100),
                  aes(fill = `2025`),
                  interpolate = F) +
  geom_sf(data = filter(obsdens.pts, `2025`>=100), aes(fill = `2025`, size = `2025`), 
          pch = 21, alpha = .8, linewidth = .1) +
  scale_fill_viridis_c(option = "plasma",
                       name = "Observations",
                       end = .3, begin = 1,
                       limits = c(1,max(values(obsdens$`2025`, na.rm = T))),
                       na.value = "transparent")  +
  scale_size_area(name = "", max_size = 11) +
  hrbrthemes::theme_ipsum_rc() +
  theme(legend.position = "bottom")
ggsave("figures/gain_cells_obsdens_alltaxa_map.png", width = 7.57, height = 6.35)

# Zoom to British Columbia

ggplot() +
  geom_sf(data = sf::st_as_sf(canada_poly), col = "grey10", fill = "grey10") +
  geom_spatraster(data = filter(obsdens, `2025`< 100),
                  aes(fill = `2025`),
                  interpolate = F) +
  geom_sf(data = filter(obsdens.pts, `2025`>=100), aes(fill = `2025`, size = `2025`),
          pch = 21, alpha = .8, linewidth = .1) +
  scale_fill_viridis_c(option = "plasma",
                       name = "Observations",
                       end = .3, begin = 1,
                       limits = c(1,max(values(obsdens$`2025`, na.rm = T))),
                       na.value = "transparent")  +
  coord_sf(xlim = c(-2146000, -1300000), ylim = c(500000, 1300000)) +
  scale_size_area(name = "", max_size = 11) +
  theme_void()
ggsave("figures/gain_cells_obsdens_alltaxa_map_BCinset.png", width = 6.45, height = 5.6)


# Zoom to Maritimes and Fleuve

ggplot() +
  geom_sf(data = sf::st_as_sf(canada_poly), col = "grey10", fill = "grey10") +
  geom_spatraster(data = filter(obsdens, `2025`< 100),
                  aes(fill = `2025`),
                  interpolate = F) +
  geom_sf(data = filter(obsdens.pts, `2025`>=100), aes(fill = `2025`, size = `2025`),
          pch = 21, alpha = .8, linewidth = .1) +
  scale_fill_viridis_c(option = "plasma",
                       name = "Observations",
                       end = .3, begin = 1,
                       limits = c(1,max(values(obsdens$`2025`, na.rm = T))),
                       na.value = "transparent")  +
  coord_sf(xlim = c(1800000, 3335000), ylim = c(-218000, 1700000)) +
  scale_size_area(name = "", max_size = 11) +
  theme_void() +
  theme(legend.position = "none")
ggsave("figures/gain_cells_obsdens_alltaxa_map_Eastinset.png", width = 6.45, height = 5.6)



## annual timeline of new cells (all groups, to put btg in context)

## compare to 2024

cell_change.2024 = lapply(maps, new_cells, n_years = 17)

res_df$taxa = factor(res_df$taxa, levels = res_df$taxa[order(res_df$area_km_newcells)])
ggplot(data = filter(res_df, taxa != "alltaxa")) +
  geom_bar(aes(x = area_km_newcells, y = taxa, fill = area_km_newcells), 
           stat = "identity") +
  scale_fill_viridis_c(option = "viridis", trans = "sqrt", end = .95) +
  labs(x = "Coverage gain (km²)", y = "Iconic Taxon Groups") +
  hrbrthemes::theme_ipsum_rc() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "none") 
ggsave("figures/gain_totalrangecoverage_area_barplot.png", width = 6.39, height = 3.05)
