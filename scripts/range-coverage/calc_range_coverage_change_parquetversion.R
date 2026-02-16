# function to calculate range coverage -----------------------------------------

# Version for when you have a full parquet file withall required years

# inat = query of the inat parquet file (filtered for the group of interest, etc.)
# three thresholds are set: 1, 3, or 10 obs per cell

# tester
# species_name = "Anaxyrus americanus"
# filepath = paste0("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/", taxagroups[1], "/")
# range = terra::rast(paste0(filepath, species_name, ".tif"))
# before = 
calc_range_coverage_change = function(species_name, range, inat = inat_pq, YEAR = "2025"){
  
  # choose a species
  # sp = "Ambystoma gracile"
  sp = species_name
  
  # select species range
  range[range>0] = 1
  
  # cut range to Canada only
  canada_poly = project(canada_poly, crs(range))
  range_can = terra::crop(range, canada_poly, mask = TRUE, touches = FALSE)
  
  # get all occurrence points of the species from inat parquet
  occ = inat |> 
    filter(scientific_name == sp,
           year != YEAR) |>
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
  
  # make a spatial layer
  after = inat |>
    filter(scientific_name == sp) |>
    select(c(longitude, latitude, scientific_name)) |>
    collect() |>
    terra::vect(
      geom = c("longitude", "latitude"),
      crs = "epsg:4326")
  after = terra::project(after, crs(range_can))
  
  # count number of points per grid cell
  range_after = rasterize(after, range_can, fun = "count", background = 0) |>
    terra::crop(range_can, mask = TRUE, touches = FALSE)

  # count number of cells in the range
  ncells_range = sum(values(range), na.rm = TRUE)
  ncells_range_can = sum(values(range_can), na.rm = TRUE)
  
  # change per pixel
  range_change = range_after - range_before
  
  # make into a percentage
  raster_stack = c(range_before, 
                   range_after,
                   range_change)
  names(raster_stack) = c("before", "after", "change")
  
  return(raster_stack)
  
}
