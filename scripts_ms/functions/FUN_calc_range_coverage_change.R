# function to calculate range coverage -----------------------------------------

# tester
# species_name = "Gyrinophilus porphyriticus"
# filepath = paste0("~/McGill University/Laura's Lab_Group - IUCN Ranges/Unclipped/EASE2.0_12.5km/AMPHIBIANS/")
# range = terra::rast(paste0(filepath, species_name, ".tif"))
# inat.pre = arrow::read_parquet("data/heavy/BTG-data/inaturalist-canada-dec2025_PRE_JUN12025.parquet")
# inat.post = arrow::read_parquet("data/heavy/BTG-data/inaturalist-canada-dec2025_JUN1OCT1.parquet")
# canada_poly = terra::vect("data/base-layers/canada-polygon/canada.outline.shp")
  
calc_range_coverage_change = function(species_name, range, inat.pre, inat.post, canada_poly){
  
  # choose a species
  sp = species_name
  
  # select species range
  range[range>0] = 1
  
  # cut range to Canada only
  canada_poly = project(canada_poly, crs(range))
  range_can = terra::crop(range, canada_poly, mask = TRUE, touches = FALSE)
  
  # get all occurrence points of the species from inat parquet
  occ = inat.pre |> 
    filter(scientific_name == sp) |>
    select(c(longitude, latitude, scientific_name)) |>
    collect()
  
  # Make a "before" raster -----------------------------------------------------
  
  # make a spatial layer
  before = terra::vect(occ, 
                       geom = c("longitude", "latitude"), 
                       crs = "epsg:4326")
  before = terra::project(before, crs(range))
  
  # count number of points per grid cell
  range_before = rasterize(before, range_can, fun = "count", background = 0) |>
    terra::crop(range_can, mask = TRUE, touches = FALSE)
  
  # Make an "after" raster -----------------------------------------------------
  
  # make a spatial layer
  after = inat.post |>
    filter(scientific_name == sp) |>
    select(c(longitude, latitude, scientific_name)) |>
    collect() |>
    terra::vect(
      geom = c("longitude", "latitude"),
      crs = "epsg:4326")
  after = terra::project(after, crs(range))
  
  # count number of points per grid cell
  range_change = rasterize(after, range_can, fun = "count", background = 0) |>
    terra::crop(range_can, mask = TRUE, touches = FALSE)
  range_after = range_before + range_change
  
  # identify and count new cells
  range_new = range_change
  range_new[range_before>0] <- 0
  range_new[range_new != 0] <- 1
  
  # measure area coverage gains
  range_area = cellSize(range, unit = "km")
  range_area[is.na(range)] <- NA

  range_can_area = cellSize(range_can, unit = "km")
  range_can_area[is.na(range_can)] <- NA

  range_new_area = cellSize(range_new, unit = "km")
  range_new_area[is.na(range_new)] <- NA
  range_new_area[range_new == 0] <- NA

  area_range = sum(values(range_area), na.rm = TRUE)
  area_range_can = sum(values(range_can_area), na.rm = TRUE)
  area_new = sum(values(range_new_area), na.rm = TRUE)
  
  res = data.frame(
    "sp" = species_name,
    "range_size_full_km" = area_range,
    "range_size_canada_km" = area_range_can,
    "range_size_newcells_km" = area_new
  )
  raster_stack = c(range_before, 
                   range_after,
                   range_change,
                   range_new)
  names(raster_stack) = c("before", "after", "change", "new")
  
  return(list("raster_stack" = raster_stack,
              "summary" = res))
  
}
