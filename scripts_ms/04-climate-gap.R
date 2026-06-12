# climate gap

library(terra)
library(tidyterra)
library(tidyverse)
library(sf)
library(ggnewscale)

# set common theme for all figs
theme_set( hrbrthemes::theme_ipsum_rc(base_size = 13) +
             theme(legend.position = "top", legend.key.width = unit(2.5,"cm")))


# DATA =========================================================================

# density maps from June 1 to Oct 1 2024 or 2025
v24 = rast("outputs/densitymaps/densitymap_PREJUN12025_alltaxa.tif")
v25 = rast("outputs/densitymaps/densitymap_JUN1OCT12025_alltaxa.tif")

# Canada polygon
canada = st_read("data/base-layers/canada-polygon/canada.outline.shp")

# climate data
PATH = "~/Documents/GitHub/sampling-scenarios/sampling-scenarios/data-raw/climate/wc2.1_10m/"
FILES = list.files(PATH)[c(1,4)] # annual mean temp and annual precip
clim = lapply(paste0(PATH,FILES), rast)

# SPATIAL OPERATIONS ===========================================================

# mask the climate layers to Canada
canada = st_transform(canada, crs(clim[[1]]))
clim = lapply(clim, crop, canada, mask = TRUE)
clim = rast(clim)

# scale and center (instead of PCA)
clim = scale(clim, scale = TRUE, center = TRUE)

# extract values
clim_vals <- values(clim, dataframe = TRUE, na.rm = TRUE)
clim_IDs <- cells(clim)

# project layers to match the climate
v24 = project(v24, crs(clim[[1]]))
v25 = project(v25, crs(clim[[1]]))

# resample to match resolution
v24 = resample(v24, y = clim[[1]], method = "sum")
v25 = resample(v25, y = clim[[1]], method = "sum")
vpost = sum(v24, v25, na.rm = T)

## EXTRACT CLIMATE =============================================================

# convert density maps to binary
b24 = v24
b24[v24>0] <- 1
b25 = v25
b25[v25>0] <- 1
bpost = vpost
bpost[vpost>0] <- 1

# extract climate from the sampled cells (=1)
clim_24 = clim[[1:2]]*b24
df_24 = values(clim_24, data.frame = TRUE, na.rm = TRUE) |> as.data.frame()
clim_25 = clim[[1:2]]*b25
df_25 = values(clim_25, data.frame = TRUE, na.rm = TRUE) |> as.data.frame()
clim_post = clim[[1:2]]*bpost
df_post = values(clim_post, data.frame = TRUE, na.rm = TRUE) |> as.data.frame()


## Make bins to define climate zones ===========================================

## Grid Canada's full climate space
(p = ggplot(data = clim_vals) +
   geom_bin_2d(aes(x = wc2.1_10m_bio_1, y = wc2.1_10m_bio_12), 
               bins = 50) +
   colorspace::scale_fill_continuous_sequential("Batlow",
                                                na.value = "transparent", 
                                                rev = F,
                                                name = "Climate frequency", 
                                                trans = "sqrt") +
   labs(x = "Annual mean temperature", y = "Annual precipitation")) 
bins_canada <- ggplot_build(p)$data[[1]] |>
  select(c(xbin, ybin, value, count, density, xmin:ymax))
colnames(bins_canada)[grep("count", colnames(bins_canada))] <- "count.canada"
colnames(bins_canada)[grep("density", colnames(bins_canada))] <- "density.canada"
bins_canada$bin = paste0(bins_canada$xbin, "_", bins_canada$ybin)

# use these limits to grid the pre 2025 and 2025 versions ----------------------

# bin inat data before BTG
df_24$xbin = NA
df_24$ybin = NA
for(i in 1:nrow(df_24)){
  index = df_24$wc2.1_10m_bio_1 >= bins_canada$xmin[i] & df_24$wc2.1_10m_bio_1 < bins_canada$xmax[i] & df_24$wc2.1_10m_bio_12 >= bins_canada$ymin[i] & df_24$wc2.1_10m_bio_12 < bins_canada$ymax[i] 
  df_24$xbin[index] <- bins_canada$xbin[i]
  df_24$ybin[index] <- bins_canada$ybin[i]
}
df_24$bin = paste0(df_24$xbin, "_", df_24$ybin)
bins_24 = df_24 |>
  group_by(bin) |>
  summarise("count.pre" = n(),
            "density.pre" = n()/nrow(df_24))

# bin the BTG samples
df_25$xbin = NA
df_25$ybin = NA
for(i in 1:nrow(df_25)){
  index = df_25$wc2.1_10m_bio_1 >= bins_canada$xmin[i] & df_25$wc2.1_10m_bio_1 < bins_canada$xmax[i] & df_25$wc2.1_10m_bio_12 >= bins_canada$ymin[i] & df_25$wc2.1_10m_bio_12 < bins_canada$ymax[i] 
  df_25$xbin[index] <- bins_canada$xbin[i]
  df_25$ybin[index] <- bins_canada$ybin[i]
}
df_25$bin = paste0(df_25$xbin, "_", df_25$ybin)
bins_25 = df_25 |>
  group_by(bin) |>
  summarise("count.25" = n(),
            "density.25" = n()/nrow(df_25))

# bin the post-BTG data
df_post$xbin = NA
df_post$ybin = NA
for(i in 1:nrow(df_post)){
  index = df_post$wc2.1_10m_bio_1 >= bins_canada$xmin[i] & df_post$wc2.1_10m_bio_1 < bins_canada$xmax[i] & df_post$wc2.1_10m_bio_12 >= bins_canada$ymin[i] & df_post$wc2.1_10m_bio_12 < bins_canada$ymax[i] 
  df_post$xbin[index] <- bins_canada$xbin[i]
  df_post$ybin[index] <- bins_canada$ybin[i]
}
df_post$bin = paste0(df_post$xbin, "_", df_post$ybin)
bins_post = df_post |>
  group_by(bin) |>
  summarise("count.post" = n(),
            "density.post" = n()/nrow(df_post))

# join into one big table
df = left_join(bins_canada, bins_24, by = "bin") |>
  left_join(bins_25, by = "bin") |>
  left_join(bins_post, by = "bin")

# calculate proportion of climate sampled
# 1 = fully sampled at least once
# 0 = none of the zone was sampled
df$sampled_proportion_post = df$count.post/df$count.canada
df$sampled_proportion_btg = df$count.25/df$count.canada
df$sampled_proportion_pre = df$count.pre/df$count.canada
df$sampled_proportion_diff = df$sampled_proportion_post - df$sampled_proportion_pre
df$sampled_proportion_sum = df$sampled_proportion_post + df$sampled_proportion_pre

# change NAs to zeros
df$sampled_proportion_pre[which(is.na(df$sampled_proportion_pre))] = 0
df$sampled_proportion_btg[which(is.na(df$sampled_proportion_btg))] = 0
df$sampled_proportion_post[which(is.na(df$sampled_proportion_post))] = 0
df$sampled_proportion_diff[which(is.na(df$sampled_proportion_diff))] = 0
df$sampled_proportion_sum[which(is.na(df$sampled_proportion_sum))] = 0

# categorise the bins based on sampling progress ===============================

df$progress = NA
df$progress[which(df$sampled_proportion_pre == 0 & df$sampled_proportion_post == 0)] <- "Never"
df$progress[which(df$sampled_proportion_pre > 0 & df$sampled_proportion_post == df$sampled_proportion_pre)] <- "No progress"
df$progress[which(df$sampled_proportion_pre == 0 & df$sampled_proportion_post > 0)] <- "First progress"
df$progress[which(df$sampled_proportion_pre > 0 & df$sampled_proportion_post > df$sampled_proportion_pre)] <- "Progress"
df$progress[which(df$sampled_proportion_pre == 1 & df$sampled_proportion_post == 1)] <- "Already complete"
df$progress[which(df$sampled_proportion_pre < 1 & df$sampled_proportion_post == 1)] <- "Completed"


# BTG progress -----------------------------------------------------

clim.progbtg = df |> 
  filter(progress %in% c("Completed", "First progress", "Progress", "No progress"))
rast_progbtg = bins_raster
for(i in 1:nrow(clim.progbtg)){
  index = clim$wc2.1_10m_bio_1 >= clim.progbtg$xmin[i] & clim$wc2.1_10m_bio_1 < clim.progbtg$xmax[i] & clim$wc2.1_10m_bio_12 >= clim.progbtg$ymin[i] & clim$wc2.1_10m_bio_12 < clim.progbtg$ymax[i] 
  rast_progbtg[index] <- clim.progbtg$sampled_proportion_diff[i]
}
plot(rast_progbtg)
rast_progbtg[rast_progbtg==0] <- NA

# coverage status after BTG
rast_status = bins_raster
for(i in 1:nrow(df)){
  index = clim$wc2.1_10m_bio_1 >= df$xmin[i] & clim$wc2.1_10m_bio_1 < df$xmax[i] & clim$wc2.1_10m_bio_12 >= df$ymin[i] & clim$wc2.1_10m_bio_12 < df$ymax[i] 
  rast_status[index] <- df$sampled_proportion_post[i]
}
plot(rast_status)
rast_status[rast_status==0] <- NA

# coverage status before BTG
rast_status_pre = bins_raster
for(i in 1:nrow(df)){
  index = clim$wc2.1_10m_bio_1 >= df$xmin[i] & clim$wc2.1_10m_bio_1 < df$xmax[i] & clim$wc2.1_10m_bio_12 >= df$ymin[i] & clim$wc2.1_10m_bio_12 < df$ymax[i] 
  rast_status_pre[index] <- df$sampled_proportion_pre[i]
}
plot(rast_status_pre)

# resample to finer res for prettier map
base1k = rast("data/base-layers/canada-basegrids/canada.base.1k.tiff")
base1k = project(base1k, crs(rast_status))
rast_status_1k = resample(rast_status, base1k, method = "bilinear")
rast_status_pre_1k = resample(rast_status_pre, base1k, method = "bilinear")

ggplot() +
  geom_spatraster(data = 100*rast_status_1k, aes(fill = wc2.1_10m_bio_1)) +
  colorspace::scale_fill_continuous_sequential("Batlow",
                                               name = "Coverage (%)",
                                               rev = F,
                                               na.value = "transparent") +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 20) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm"))
ggsave("figures/climate gap/map_climatespace_coveragepostBTG.png", width = 12.5, height = 7.9)



rast_status_diff = rast_status_1k - rast_status_pre_1k
rast_status_sum = rast_status_1k + rast_status_pre_1k
# 2 = totally covered + redundant sampling
# > 1 = well covered
# < 1 = undersampled
# < 0.5 = very undersampled

# clip raster to clean up the edges
canada = st_transform(canada, crs(rast_status_sum))
rast_status_sum_clip = terra::mask(rast_status_sum, canada)
rast_status_sum_clip = trim(rast_status_sum_clip)
rast_status_sum_clip[rast_status_sum_clip==0] <- NA
canada = sf::st_transform(canada, crs = "EPSG:3347")
rast_status_sum_clip = project(rast_status_sum_clip, crs = crs(canada))

(map.climgap = ggplot() +
  geom_sf(data = canada, fill = "grey80", col = "white", linewidth = .1) +
  geom_spatraster(data = rast_status_sum_clip, aes(fill = wc2.1_10m_bio_1)) +
  colorspace::scale_fill_continuous_divergingx("Temps", mid = 1,
                                               name = "Coverage",
                                               rev = T,
                                               na.value = "transparent")  +
  hrbrthemes::theme_ipsum_rc(grid = FALSE, 
                             axis = FALSE, 
                             axis_text_size = 1,
                             ticks = FALSE,
                             base_size = 20) +
  theme(legend.position = "top", legend.key.width = unit(2,"cm")))
ggsave("figures/climate gap/map_climatespace_coveragegainBTG.png", width = 12.5, height = 7.9)

# pca plot of sampling 
rast_status_sum_clip_df = values(rast_status_sum_clip, as.data.frame = TRUE)

df$x_coord = (df$xmin+df$xmax)/2
df$y_coord = (df$ymin+df$ymax)/2

(p = ggplot() +
    # All climate bins
    geom_tile(data = df,
                aes(x = x_coord, y = y_coord,
                    fill = sampled_proportion_sum)) +
    ## Grey out the unsampled climate bins
    geom_tile(data = filter(df, sampled_proportion_sum == 0),
              aes(x = x_coord, y = y_coord, fill = count.canada), 
              fill = "grey85") +
    ## Outline the UNsampled climates during BTG
    geom_tile(data = filter(df, sampled_proportion_btg == 0 & sampled_proportion_sum > 0),
              aes(x = x_coord, y = y_coord, fill = NULL),  
               col = "black", 
              alpha = 0, linewidth = .25, linetype = 1) +
    colorspace::scale_fill_continuous_divergingx("Temps",
                                                 na.value = "transparent", 
                                                 rev = T,
                                                 mid = 1,
                                                 name = "") +
    labs(x = "Annual mean temperature", y = "Annual precipitation"))  +
  hrbrthemes::theme_ipsum_rc(base_size = 20,
                             axis_title_size = 24) +
  theme(legend.position = "top", legend.key.width = unit(3,"cm"))
ggsave("figures/climate gap/pca_coveragegainsum.png", width = 7.83, height = 8.5)