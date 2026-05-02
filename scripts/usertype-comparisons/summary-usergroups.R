# Script to make some summary plots of user group contributions to inat

# SR/obs
# NewSp / year
# geog coverage / year
# clim coverage / year

# user groups: superusers, qcbs, bc biodiversity teams

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr)
library(duckdb)

# set ggplot theme
theme_set(ggpubr::theme_pubr())

# user group logins

# Incentivized: QCBS
qcbs <- readRDS("~/Documents/GitHub/storymap-qcbs/data/championteams_userlogins.rds")
qcbs = c(qcbs, "katherinehebert") # i'm missing!
qcbs = gsub("and abrunet27", "abrunet27", qcbs)
# Employed: BC Biodiversity Program
bcteam = read.csv("~/Documents/GitHub/BCParks/outputs/bigteam_userlogins.csv", row.names = 1)
bcteam =  bcteams$user_login

# load parquet
inat <- arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")