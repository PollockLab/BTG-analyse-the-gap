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


## annual timeline of new cells (all groups, to put btg in context)