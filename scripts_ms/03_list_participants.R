# Measure and plot the engagement of different user groups compared to 2024
# including BTG members

# libraries
library(ggplot2)
library(dplyr)
library(arrow)
library(dbplyr, warn.conflicts = FALSE)
library(duckdb)
library(sf)
library(patchwork)

# set ggplot theme
theme_set(hrbrthemes::theme_ipsum_rc())

# load parquet file ------------------------------------------------------------

inat = arrow::open_dataset("data/heavy/BTG-data/inaturalist-canada-dec2025_smaller.parquet")
pq = inat

# load user groups -------------------------------------------------------------

# members on the iNat BTG project (umbrella and collections)
btg = read.csv("data/user_information/blitz_the_gap_members.csv")
btg = btg$user_login

# BC Biodiversity Program teams
bc = read.csv("data/user_information/bigteam_userlogins.csv", row.names = 1)
bc = bc$user_login
bc2025 = read.csv("data/user_information/BC Big Summer Teams - Years.csv")
bc2025 = bc2025 |> filter(Year == 2025)
bc2025 = bc2025$User

# QCBS
qcbs = read.csv("data/user_information/championteams_userlogins.csv")
qcbs = qcbs$x

# CWF's Great Canadian Bioblitz grantees
gcb = readxl::read_xlsx("data/user_information/GCB expert obs - fromJames.xlsx")
gcb = gcb$Username

# CWF GCBioblitz participants summary
# these include QCBS bioblitzers. Important for sorting later.
cwf = read.csv("outputs/bioblitzes/bioblitzers-obscounts.csv") |>
  filter(category == "Great Canadian Bioblitz")
cwf = cwf$user_login

# campus bioblitzes
cnc = read.csv("outputs/bioblitzes/bioblitzers-obscounts.csv") |>
  filter(category == "Campus Nature Challenge")
cnc = cnc$user_login

# make user ID data frame
users = data.frame(
  "group" = c(
    rep("Blitz the Gap", length(btg)),
    rep("BC Biodiversity Program", length(bc)),
    rep("BC Biodiversity 2025 Team", length(bc2025)),
    rep("QCBS Grantees", length(qcbs)),
    rep("GCB Grantees", length(gcb)),
    rep("CWF Bioblitzes", length(cwf)),
    rep("Campus Nature Challenge", length(cnc))
  ),
  "user_login" = c(btg, bc, bc2025, qcbs, gcb, cwf, cnc)
)
# cut to just the unique IDs 
users.unique = unique(users$user_login) # these are "BTG/CWF people" 1771 people

# add primary group (some people were in multiple)
# the order is important here:
# first, fill large general groups, then 2025 team, grantees, and fill gaps with BTG members
users$primary_group = NA
users$primary_group[which(users$user_login %in% btg)] = "Blitz the Gap"
users$primary_group[which(users$user_login %in% cnc)] = "Campus Nature Challenge"
users$primary_group[which(users$user_login %in% cwf)] = "CWF Bioblitzes"
users$primary_group[which(users$user_login %in% bc)] = "BC Biodiversity"
users$primary_group[which(users$user_login %in% bc2025)] = "BC Biodiversity 2025 Team"
users$primary_group[which(users$user_login %in% gcb)] = "GCB Grantees"
users$primary_group[which(users$user_login %in% qcbs)] = "QCBS Grantees"
# manually correct users we know
users$primary_group[which(users$user_login == "ajones_mcgill")] <- "QCBS Grantees"
users_all = users
write.csv(users_all, "outputs/users/btg_users_participants.csv")

# distinct users with primary groupings
users_unique = users |>
  select(user_login, primary_group) |> 
  distinct()
write.csv(users_unique, "outputs/users/btg_users_participants_primarygroups.csv")
