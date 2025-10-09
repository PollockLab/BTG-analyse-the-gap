# Script to get user IDs from big summer teams

# load library
library(rinat)
library(dplyr)

# get iNaturalist project id
info = rinat::get_inat_obs_project("bc-big-summer-teams", type = "info")

# get stats per user
# note: this could be a good way to analyse user contributions if we want to, later!
year = 2019:2025
users = list()
for(t in 1:length(year)){
  users[[t]] = rinat::get_inat_user_stats(project = info$id,
                                          date_range = c(paste0(year[t], "-05-01", year[t], "-12-31")))
}

# just keep the logins
logins = list()
for(t in 1:length(year)){
  logins[[t]] = users[[t]]$most_species$user$login
}
names(logins) = year

# convert to data frames
logins = lapply(logins, as.data.frame)

# collapse the list into a data frame with a column "year" that holds each list element's name
logins = dplyr::bind_rows(logins, .id = "year")

# clean up
colnames(logins)[2] = "user_login"
logins$year = as.numeric(logins$year)

# write to csv
write.csv(logins, "data/user_information/bcparks_summerteams_userlogins.csv", row.names = FALSE)