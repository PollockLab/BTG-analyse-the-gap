# function to calculate range coverage -----------------------------------------

# inat = query of the inat parquet file (filtered for the group of interest, etc.)
# three thresholds are set: 1, 3, or 10 obs per cell

calc_range_coverage = function(species_name, range, inat = query){
  
  # choose a species
  # sp = "Panax quinquefolius"
  sp = species_name
  
  # select species range
  range[range>0] = 1
  range = terra::project(range, crs(canada))
  # cut range to Canada only
  range_can = terra::mask(range, canada_poly)

  # get occurrence points of the species
  occ = inat_pq |> 
    filter(scientific_name == sp) |>
    select(c(longitude, latitude, scientific_name)) |>
    collect()

    # make a spatial layer
  occ = terra::vect(occ, 
                    geom = c("longitude", "latitude"), 
                    crs = crs(canada))
  
  # count number of cells in the range
  ncells_range = sum(values(range), na.rm = TRUE)
  ncells_range_can = sum(values(range_can), na.rm = TRUE)
  
  # count number of points per grid cell
  range_occ = rasterize(occ, range_can, fun = "length", background = 0)
  range_occ = mask(range_occ, range_can)
  
  # set ranges to 1 or 0 if enough observations have been found
  range_enoughdata_min1 = range_occ
  range_enoughdata_min1[range_occ >= 1] <- 1
  range_enoughdata_min1[range_occ < 1] <- 0
  
  # set ranges to 1 or 0 if enough observations have been found
  range_enoughdata_min3 = range_occ
  range_enoughdata_min3[range_occ >= 3] <- 1
  range_enoughdata_min3[range_occ < 3] <- 0
  
  # set ranges to 1 or 0 if enough observations have been found
  range_enoughdata_min10 = range_occ
  range_enoughdata_min10[range_occ >= 10] <- 1
  range_enoughdata_min10[range_occ < 10] <- 0
  
  # count number of cells with enough data
  ncells_enoughdata = c(
    sum(values(range_enoughdata_min1), na.rm = TRUE),
    sum(values(range_enoughdata_min3), na.rm = TRUE),
    sum(values(range_enoughdata_min10), na.rm = TRUE)
  )
  
  # make into a percentage
  coverage = ncells_enoughdata
  
  return(list("coverage" = coverage,
              "n_cells_canada" = ncells_range_can,
              "n_cells_fullrange" = ncells_range))
}

