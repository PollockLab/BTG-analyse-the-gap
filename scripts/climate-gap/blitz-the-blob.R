# climate space of canada coverage!!
# also who wins in terms of climate and spatial distance from other obs this year? 

library(sf)
library(tidyverse)
library(terra)
library(tidyterra)
library(geodata)
library(patchwork)
library(biscale)

# set common theme for all maps
theme_set( hrbrthemes::theme_ipsum_rc(base_size = 13) +
             theme(legend.position = "top", legend.key.width = unit(2.5,"cm")))


# climate space of canada

# climate polygon
canada = st_read("data/base-layers/canada-polygon/canada.outline.shp")

# import worldclim
PATH = "~/Documents/GitHub/sampling-scenarios/sampling-scenarios/data-raw/climate/wc2.1_10m/"
FILES = list.files(PATH)[c(1,4)] # annual mean temp and annual precip
clim = lapply(paste0(PATH,FILES), rast)

# mask the climate layers to canada
canada = st_transform(canada, crs(clim[[1]]))
clim = lapply(clim, crop, canada, mask = TRUE)
clim = rast(clim)

# scale and center (instead of PCA)
clim = scale(clim, scale = TRUE, center = TRUE)

# extract values
clim_vals <- values(clim, dataframe = TRUE, na.rm = TRUE)
clim_IDs <- cells(clim)


# make a PCA
#pca <- prcomp(clim_vals, center = TRUE, scale = TRUE)

# This creates a new raster where layers are PC1, PC2, etc.
# pca_raster <- predict(clim, pca)

# map the PCA raster
pc1 = ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = clim, aes(fill = wc2.1_10m_bio_1)) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                               na.value = "transparent", 
                                               rev = F,
                                               name = "Annual Mean\nTemperature") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
pc2 = ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = clim, aes(fill = wc2.1_10m_bio_12)) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                               na.value = "transparent", 
                                               rev = F,
                                               name = "Annual\nPrecipitation") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
pc1 + pc2
# ggsave("figures/climate gap/map_PC1PC2.png", width = 12.3, height = 8)
ggsave("figures/climate gap/map_BIO1BIO12.png", width = 12.3, height = 8)

# raster of sampled cells
newcells = rast("outputs/comparison-maps/density/newcells_alltaxa_2025_vs_20082024.tif")$alltaxa
newcells = project(newcells, crs(clim))
newcells = resample(newcells, clim[[1]], method = "max")
newcells[newcells > 0] <- 1
clim_newcells = clim[[1:2]]*newcells
plot(clim_newcells)
df_new = values(clim_newcells, data.frame = TRUE, na.rm = TRUE)

# raster of cell sampling change
change = rast("outputs/comparison-maps/density/changecells_alltaxa_2025_vs_20082024.tif")$alltaxa
changecells = project(change, crs(clim))
changecells = resample(changecells, clim[[1]], method = "max")
changecells[changecells == 1] <- NA # this is a map of cells that were sampled from previous years in or not in 2025
# i.e. removed cells that were only sampled in 2025
changecells[!is.na(changecells)] <- 1
plot(changecells)

clim_changecells = clim[[1:2]]*changecells
plot(clim_changecells)
df_change = values(clim_changecells, data.frame = TRUE, na.rm = TRUE)


## PLOTTING --------------------------------------------------------------------

# plot Canada's climate space
df = clim_vals

# full climate space
(p = ggplot(data = df) +
  geom_bin_2d(aes(x = wc2.1_10m_bio_1, y = wc2.1_10m_bio_12), 
              bins = 50) +
    colorspace::scale_fill_continuous_sequential("Batlow",
                                                 na.value = "transparent", 
                                                 rev = F,
                                                 name = "Climate frequency", 
                                                 trans = "sqrt") +
    labs(x = "Annual mean temperature", y = "Annual precipitation")) 
bins_canada <- ggplot_build(p)$data[[1]]
ggsave("figures/climate gap/pca_canada.png", width = 8.26, height = 7.75)

# sampled climate space before 2025

(p2 = ggplot() +
  geom_bin_2d(data = df,
              aes(x = wc2.1_10m_bio_1, y = wc2.1_10m_bio_12, fill = NA),
              bins = 50, fill = "grey") +
  geom_bin_2d(data = df_change, 
              aes(x = wc2.1_10m_bio_1, y = wc2.1_10m_bio_12),
              bins = 50) +
    colorspace::scale_fill_continuous_sequential("Batlow",
                                                 na.value = "transparent", 
                                                 rev = F,
                                                 name = "Climate frequency", 
                                                 trans = "sqrt") +
    labs(x = "Annual mean temperature", y = "Annual precipitation"))
bins_pre2025 <- ggplot_build(p2)$data[[2]]
ggsave("figures/climate gap/pca_sampledbefore2025.png", width = 8.26, height = 7.75)


## map it geographically




# sampled climate space in 2025

(p3 = ggplot() +
  geom_bin_2d(data = df, aes(x = wc2.1_10m_bio_1, y = wc2.1_10m_bio_12, fill = NA), 
              bins = 50, fill = "grey") +
  geom_bin_2d(data = df_new, aes(x = wc2.1_10m_bio_1, y = wc2.1_10m_bio_12), 
              bins = 50) +
    colorspace::scale_fill_continuous_sequential("Batlow",
                                                 na.value = "transparent", 
                                                 rev = F,
                                                 name = "Climate frequency", 
                                                 trans = "sqrt") +
    labs(x = "Annual mean temperature", y = "Annual precipitation"))
ggsave("figures/climate gap/pca_sampledin2025.png", width = 8.26, height = 7.75)
bins_2025 <- ggplot_build(p3)$data[[2]]

## calculate difference between sampled
bins1 = bins_canada 
bins2 = bins_pre2025 |> select(c(x,y,count))
bins3 = bins_2025 |> select(c(x,y,count))

bins = left_join(bins1, bins2, by = c("x","y"), suffix = c(".canada", ".pre2025"))
bins = left_join(bins, bins3, by = c("x","y"), suffix = c("", ".2025"))

bins$diff = (bins$count-bins$count.pre2025)/(bins$count+bins$count.pre2025)
ggplot(bins, aes(x = x, y = y, fill = diff)) +
  geom_tile() +
  colorspace::scale_fill_continuous_divergingx("Spectral", 
                                               name = "",
                                               rev = T,
                                               # name = "2025/Past",
                                               na.value = "grey90",
                                               limits = c(-1,1)) +
labs(x = "Annual mean temperature", y = "Annual precipitation") +
       #caption = ">0 = more visited in 2025 than past years\n<0 = less visited in 2025 than past years") +
  hrbrthemes::theme_ipsum_rc(axis_text_size = 18,
                             axis_title_size = 20,
                             axis = F, 
                             base_size = 18) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/pca_relativesamplingcells.png", width = 7.45, height = 6.7)


## map the differences in geographical space

bins_raster = clim[[1]] 
bins_raster[!is.na(bins_raster)] <- 0
plot(bins_raster)

bins$BIN = paste0(bins$xbin, "_", bins$ybin)

for(i in 1:nrow(bins)){
  index = clim$wc2.1_10m_bio_1 >= bins$xmin[i] & clim$wc2.1_10m_bio_1 < bins$xmax[i] & clim$wc2.1_10m_bio_12 >= bins$ymin[i] & clim$wc2.1_10m_bio_12 < bins$ymax[i] 
  bins_raster[index] <- bins$diff[i]
}
plot(bins_raster)
writeRaster(bins_raster, "outputs/climate-gap/map-sampledclimates-bins-diff2025pre2025.tif", overwrite = T)
bins_raster = rast("outputs/climate-gap/map-sampledclimates-bins-diff2025pre2025.tif")

# pretty map
ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = bins_raster) +
  colorspace::scale_fill_continuous_divergingx("Spectral", 
                                               name = "",
                                               rev = T,
                                               # name = "2025/Past",
                                               na.value = "transparent",
                                               limits = c(-1,1)) +
  # scale_fill_distiller(palette = "Spectral",
  #                                              name = "Relative sampling\n of climate space",
  #                                              na.value = "transparent",
  #                                              limits = c(-1,1)) +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 20) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/map_relativesampling_climatespace.png", width = 12.5, height = 7.9)

bins_raster_onlygains = bins_raster
bins_raster_onlygains[bins_raster_onlygains<=0] <- NA
ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = bins_raster_onlygains) +
  colorspace::scale_fill_continuous_divergingx("Zissou 1", 
                                               name = "",
                                               rev = F,
                                               # name = "2025/Past",
                                               na.value = "transparent",
                                               limits = c(0,1)) +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/map_relativesampling_gainsonly_climatespace.png", width = 12.5, height = 7.9)

## map the density in geographical space
bins_raster = clim[[1]] 
bins_raster[!is.na(bins_raster)] <- 0
bins$BIN = paste0(bins$xbin, "_", bins$ybin)
for(i in 1:nrow(bins)){
  index = clim$wc2.1_10m_bio_1 >= bins$xmin[i] & clim$wc2.1_10m_bio_1 < bins$xmax[i] & clim$wc2.1_10m_bio_12 >= bins$ymin[i] & clim$wc2.1_10m_bio_12 < bins$ymax[i] 
  bins_raster[index] <- bins$count.canada[i]
}
plot(bins_raster)
writeRaster(bins_raster, "outputs/climate-gap/map-sampledclimates-bins-countcanada.tif", overwrite = T)
bins_raster = rast("outputs/climate-gap/map-sampledclimates-bins-countcanada.tif")

# pretty map
ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = bins_raster) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                                name = "Climate frequency",
                                                na.value = "transparent", rev = F) +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/map_climate_binfrequency.png", width = 12.5, height = 7.9)


## map the climate sampling density in geographical space
bins_raster = clim[[1]] 
bins_raster[!is.na(bins_raster)] <- 0
bins$BIN = paste0(bins$xbin, "_", bins$ybin)
for(i in 1:nrow(bins)){
  index = clim$wc2.1_10m_bio_1 >= bins$xmin[i] & clim$wc2.1_10m_bio_1 < bins$xmax[i] & clim$wc2.1_10m_bio_12 >= bins$ymin[i] & clim$wc2.1_10m_bio_12 < bins$ymax[i] 
  bins_raster[index] <- bins$count.pre2025[i]
}
plot(bins_raster)
writeRaster(bins_raster, "outputs/climate-gap/map-sampledclimates-bins-countpre2025.tif", overwrite = T)
bins_raster = rast("outputs/climate-gap/map-sampledclimates-bins-countpre2025.tif")

# pretty map
ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = bins_raster) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                               name = "Climate sampling density",
                                               na.value = "transparent", rev = F) +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/map_climate_samplingfrequency_pre2025.png", width = 12.5, height = 7.9)


## map the climate sampling density in geographical space
bins_raster = clim[[1]] 
bins_raster[!is.na(bins_raster)] <- 0
bins$BIN = paste0(bins$xbin, "_", bins$ybin)
for(i in 1:nrow(bins)){
  index = clim$wc2.1_10m_bio_1 >= bins$xmin[i] & clim$wc2.1_10m_bio_1 < bins$xmax[i] & clim$wc2.1_10m_bio_12 >= bins$ymin[i] & clim$wc2.1_10m_bio_12 < bins$ymax[i] 
  bins_raster[index] <- bins$count[i]
}
plot(bins_raster)
writeRaster(bins_raster, "outputs/climate-gap/map-sampledclimates-bins-count2025.tif", overwrite = T)
bins_raster = rast("outputs/climate-gap/map-sampledclimates-bins-count2025.tif")

# pretty map
ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = bins_raster) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                               name = "Climate sampling density",
                                               na.value = "transparent", rev = F) +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/map_climate_samplingfrequency_2025.png", width = 12.5, height = 7.9)


## Residuals -------------------------------------------------------------------

bins$clim_perc = bins$count.canada/sum(bins$count.canada, na.rm = T)
bins$samp_perc = (bins$count.pre2025/sum(bins$count.pre2025, na.rm = T))
bins$resids = bins$samp_perc - bins$clim_perc

ggplot(data=bins) +
  geom_point(aes(x = count.canada/sum(bins$count.canada, na.rm = T),
                 y = count.pre2025/sum(bins$count.pre2025, na.rm = T))) +
  geom_abline(intercept = 0, slope = 1) +
  coord_equal(ylim = c(0, .025), xlim = c(0,.025))


## map the climate sampling density in geographical space
bins_raster = clim[[1]] 
bins_raster[!is.na(bins_raster)] <- 0
bins$BIN = paste0(bins$xbin, "_", bins$ybin)
for(i in 1:nrow(bins)){
  index = clim$wc2.1_10m_bio_1 >= bins$xmin[i] & clim$wc2.1_10m_bio_1 < bins$xmax[i] & clim$wc2.1_10m_bio_12 >= bins$ymin[i] & clim$wc2.1_10m_bio_12 < bins$ymax[i] 
  bins_raster[index] <- bins$resids[i]
}
writeRaster(bins_raster, "outputs/climate-gap/map-sampledclimates-bins-resids.tif", overwrite = T)
bins_raster = rast("outputs/climate-gap/map-sampledclimates-bins-resids.tif")

canada = st_transform(canada, "EPSG:3347")
bins_raster = project(bins_raster, "EPSG:3347") |> trim()

max.lims = bins$resids |> range(na.rm=T) |> abs() |> max()
ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = bins_raster) +
  colorspace::scale_fill_continuous_divergingx("Spectral",
                                               name = "",
                                               #name = "Climate sampling balance",
                                               na.value = "transparent", rev = T, 
                                               limits = c(-max.lims, max.lims)) +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(1.5,"cm"))
ggsave("figures/climate gap/map_climate_samplingbalance.png", width = 12.5, height = 7.9)

# plot in climate space

ggplot(bins, aes(x = x, y = y, fill = resids)) +
  geom_tile() +
  colorspace::scale_fill_continuous_divergingx("Spectral", 
                                               name = "",
                                               rev = T,
                                               na.value = "grey90",
                                               limits = c(-max.lims, max.lims)) +
  labs(x = "Annual mean temperature", y = "Annual precipitation") +
  #caption = ">0 = more visited in 2025 than past years\n<0 = less visited in 2025 than past years") +
  hrbrthemes::theme_ipsum_rc(axis_text_size = 18,
                             axis_title_size = 20,
                             axis = F, 
                             base_size = 18) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/pca_relativesampling_alltime.png", width = 7.45, height = 6.7)



## RESIDS WITH 2025 -----------------------------------------------------------

bins$samp_perc.2025 = ((bins$count + bins$count.pre2025)/sum(bins$count+bins$count.pre2025, na.rm = T))
bins$resids.2025 = bins$samp_perc.2025 - bins$clim_perc
bins$resids.diff = bins$resids-bins$resids.2025

## map the climate sampling density in geographical space
bins_raster = clim[[1]] 
bins_raster[!is.na(bins_raster)] <- 0
bins$BIN = paste0(bins$xbin, "_", bins$ybin)
for(i in 1:nrow(bins)){
  index = clim$wc2.1_10m_bio_1 >= bins$xmin[i] & clim$wc2.1_10m_bio_1 < bins$xmax[i] & clim$wc2.1_10m_bio_12 >= bins$ymin[i] & clim$wc2.1_10m_bio_12 < bins$ymax[i] 
  bins_raster[index] <- bins$resids.diff[i]
}
writeRaster(bins_raster, "outputs/climate-gap/map-sampledclimates-bins-residsdiff.tif", overwrite = T)
bins_raster = rast("outputs/climate-gap/map-sampledclimates-bins-residsdiff.tif")

canada = st_transform(canada, "EPSG:3347")
bins_raster = project(bins_raster, "EPSG:3347") |> trim()

max.lims = bins$resids.diff |> range(na.rm=T) |> abs() |> max()
max.lims = max.lims + .0001
ggplot() +
  geom_sf(data = canada, fill = "grey90", col = "grey90") +
  geom_spatraster(data = bins_raster) +
  colorspace::scale_fill_continuous_divergingx("Spectral",
                                               name = "",
                                               #name = "Climate sampling balance",
                                               na.value = "transparent", rev = T, 
                                               limits = c(-max.lims, max.lims)) +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(1.5,"cm"))
ggsave("figures/climate gap/map_climate_relativesamplingbalance.png", width = 12.5, height = 7.9)


# plot in climate space

ggplot(bins, aes(x = x, y = y, fill = resids.diff)) +
  geom_tile() +
  colorspace::scale_fill_continuous_divergingx("Spectral", 
                                               name = "",
                                               rev = T,
                                               na.value = "grey90",
                                               limits = c(-max.lims, max.lims)) +
  labs(x = "Annual mean temperature", y = "Annual precipitation") +
  #caption = ">0 = more visited in 2025 than past years\n<0 = less visited in 2025 than past years") +
  hrbrthemes::theme_ipsum_rc(axis_text_size = 18,
                             axis_title_size = 20,
                             axis = F, 
                             base_size = 18) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/pca_relativesamplingbalance.png", width = 7.45, height = 6.7)


