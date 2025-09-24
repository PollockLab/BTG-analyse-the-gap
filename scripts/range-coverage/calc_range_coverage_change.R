# function to calculate range coverage -----------------------------------------

# inat = query of the inat parquet file (filtered for the group of interest, etc.)
# three thresholds are set: 1, 3, or 10 obs per cell

# tester
# species_name = "Anaxyrus americanus"
# filepath = paste0("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/", taxagroups[1], "/")
# range = terra::rast(paste0(filepath, species_name, ".tif"))

calc_range_coverage_change = function(species_name, range, inat = inat_pq){
  
  # choose a species
  # sp = "Panax quinquefolius"
  sp = species_name
  
  # select species range
  range[range>0] = 1
  
  # cut range to Canada only
  canada_poly = project(canada_poly, crs(range))
  range_can = terra::crop(range, canada_poly, mask = TRUE, touches = FALSE)
  
  # get all occurrence points of the species from inat parquet
  occ = inat |> 
    filter(scientific_name == sp) |>
    select(c(longitude, latitude, scientific_name)) |>
    collect()
  
  # Make a "before" raster -----------------------------------------------------
  
  # make a spatial layer
  before = terra::vect(occ, 
                       geom = c("longitude", "latitude"), 
                       crs = "epsg:4326")
  before = terra::project(before, crs(range_can))
  
  # count number of points per grid cell
  range_before = rasterize(before, range_can, fun = "count", background = 0) |>
    terra::crop(range_can, mask = TRUE, touches = FALSE)

  # Make an "after" raster -----------------------------------------------------
  
  btg.occ <- rinat::get_inat_obs(taxon_name = species_name, 
                          year = 2025, 
                          geo = TRUE,
                          maxresults = 10000,
                          bounds = c(41.310824,-157.500000,83.420215,-49.921875),
                          quality = "research") |>
    tidyr::separate("observed_on",
                    into =  c("year", "month", "day")) |>
    filter(month %in% c("05", "06", "07", "08", "09"),
           captive_cultivated == "false",
           scientific_name == species_name) |>
    select(c(longitude, latitude, scientific_name))
  
  # join pre-BTG and during BTG datasets
  occ = rbind(occ, btg.occ)
  
  # make a spatial layer
  occ = terra::vect(occ, 
                       geom = c("longitude", "latitude"), 
                       crs = "epsg:4326")
  occ = terra::project(occ, crs(range_can))
  
  # count number of cells in the range
  ncells_range = sum(values(range), na.rm = TRUE)
  ncells_range_can = sum(values(range_can), na.rm = TRUE)
  
  # count number of points per grid cell
  range_after = terra::rasterize(occ, range_can, fun = "count", background = 0) |>
  terra::crop(range_can, mask = TRUE, touches = FALSE)
  
  # change per pixel
  range_change = range_after - range_before
  
  # make into a percentage
  raster_stack = c(range_before, 
                   range_after,
                   range_change)
  names(raster_stack) = c("before", "after", "change")
  
  return(raster_stack)
              
}
