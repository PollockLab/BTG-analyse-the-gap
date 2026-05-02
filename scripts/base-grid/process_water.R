# This script is by Maximiliane Jousse 
# base layer saved in canada-basegrids folder

#' Processing `WaterUrbanBuiltMask.tif` & create base layer.
#' 
#' raster layer for US-CAN with water (value == 1) and urban areas (value == 2). 
#' Values are continuous (?isaac? what does this mean). 
#' 
#' 
#' Process:
#' 1) Make categorical
#' 2) Resample to desired resolution (1k, 5k)
#' 3) Check and Save
#' 4) create a base grid


source("./fn/coarse.grid.R")


x <- rast("./WaterUrbanBuiltMask.tif")

y <- rast("./Can_1km.tiff")

y.5k <- coarse.grid(5, y)


#1. CHANGE TO CATEGORICAL---------------------
## I.e. 1 -> Water, 2 -> cities
unique(values(x))
x <- floor(x)

#2. RESAMPLE to get same res :) -------------
#' Method: `"near"`, `"mode"` or `"min"`/`"max"`.

# Exploration:
plot(resample(x, y, method = "near"))
plot(resample(x, y, method = "max"))
plot(resample(x, y, method = "mode"))

# Resampling
x.1k <- resample(x, y, method = "mode")
x.5k <- resample(x, y.5k, method = "mode")

#3. CROP ----------------------------------
# check same CRS 
crs(x.1k) == crs(y)
crs(x.5k) == crs(y.5k)

# can crop with x using mask :)
x.1k <- mask(x.1k, y)
x.5k <- mask(x.5k, y.5k)

#4. CHECK + SAVE ----------------------
ext(x.1k) == ext(y)
ext(x.5k) == ext(y.5k)

plot(x.1k)
plot(x.5k)

writeRaster(x.1k, "./WaterUrban_1km.tiff", overwrite = TRUE)
writeRaster(x.5k, "./WaterUrban_5km.tiff", overwrite = TRUE)


#5. CREATE BASE LAYER (1k, 5k) with WATER removed only.
x.1k[x.1k == 1] <- NA
x.5k[x.5k == 1] <- NA

y <- mask(y, x.1k)
y[y == 2] <- 1        # the baselayer we have has a value of "2". for pretty we change to 1 :)
writeRaster(y, "./Can_1km_basegrid.tiff")


y.5k <- mask(y.5k, x.5k)
y.5k[y.5k == 2] <- 1
writeRaster(y.5k, "./Can_5km_basegrid.tiff")
