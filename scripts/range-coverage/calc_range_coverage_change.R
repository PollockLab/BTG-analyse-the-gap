# function to calculate range coverage -----------------------------------------

# inat = query of the inat parquet file (filtered for the group of interest, etc.)
# three thresholds are set: 1, 3, or 10 obs per cell

# tester
# species_name = "Anaxyrus americanus"
# range = terra::rast(paste0(filepath, species_name, ".tif"))

calc_range_coverage_change = function(species_name, range, inat = inat_pq){
  
  # choose a species
  # sp = "Panax quinquefolius"
  sp = species_name
  
  # select species range
  range[range>0] = 1
  canada = project(canada, crs(range))
  # cut range to Canada only
  range_can = terra::resample(range, canada)
  
  # get all occurrence points of the species from inat parquet
  occ = inat |> 
    filter(scientific_name == sp) |>
    select(c(longitude, latitude, scientific_name)) |>
    collect()
  
  # Make a "before" raster -----------------------------------------------------
  
  before = occ
  before$longitude = before$longitude*100000
  before$latitude = before$latitude*100000
  
  # make a spatial layer
  before = terra::vect(before, 
                       geom = c("longitude", "latitude"), 
                       crs = crs(range_can))
  
  # count number of points per grid cell
  range_before = rasterize(before, range_can, fun = "length", background = 0)

  # Make a "after" raster ------------------------------------------------------
  
  # get new inat obs from BTG
  btg.occ <- get_inat_obs(taxon_name = species_name, 
                          year = 2025, 
                          geo = TRUE,
                          maxresults = 100,
                          bounds = bounds,
                          quality = "research") |>
    tidyr::separate("observed_on",
                    into =  c("year", "month", "day")) |>
    filter(month %in% c("05", "06", "07", "08", "09"),
           captive_cultivated == "false",
           scientific_name == species_name) |>
    select(c(longitude, latitude, scientific_name))
  
  # join pre-BTG and during BTG datasets
  occ = rbind(occ, btg.occ)
  occ$longitude = occ$longitude*100000
  occ$latitude = occ$latitude*100000
  
  # make a spatial layer
  occ = terra::vect(occ, 
                    geom = c("longitude", "latitude"), 
                    crs = crs(range_can))
  
  # count number of cells in the range
  ncells_range = sum(values(range), na.rm = TRUE)
  ncells_range_can = sum(values(range_can), na.rm = TRUE)
  
  # count number of points per grid cell
  range_after = rasterize(occ, range_can, fun = "length", background = 0)

  # set ranges to 1 or 0 if enough observations have been found
  range_enoughdata_min1 = range_after
  range_enoughdata_min1[range_after >= 1] <- 1
  range_enoughdata_min1[range_after < 1] <- 0
  
  # set ranges to 1 or 0 if enough observations have been found
  range_enoughdata_min3 = range_after
  range_enoughdata_min3[range_after >= 3] <- 1
  range_enoughdata_min3[range_after < 3] <- 0
  
  # set ranges to 1 or 0 if enough observations have been found
  range_enoughdata_min10 = range_after
  range_enoughdata_min10[range_after >= 10] <- 1
  range_enoughdata_min10[range_after < 10] <- 0
  
  # count number of cells with enough data
  ncells_enoughdata = c(
    sum(values(range_enoughdata_min1), na.rm = TRUE),
    sum(values(range_enoughdata_min3), na.rm = TRUE),
    sum(values(range_enoughdata_min10), na.rm = TRUE)
  )
  
  # make into a percentage
  coverage = ncells_enoughdata
  raster_stack = c(range_before, 
                   range_after,
                   range_enoughdata_min1,
                   range_enoughdata_min3,
                   range_enoughdata_min10)
  names(raster_stack) = c("before", "after", "min1", "min3", "min10")
  raster_stack = terra::crop(raster_stack, range_can, mask = TRUE)
  
  return(list("coverage" = coverage,
              "n_cells_canada" = ncells_range_can,
              "n_cells_fullrange" = ncells_range,
              "rasters" = raster_stack))
}
