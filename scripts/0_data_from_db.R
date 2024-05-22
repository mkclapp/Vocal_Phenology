# This script takes as inputs the .db files created in Audiodash and writes them to an RData object. Databases include:
# 1) Sampling effort metadata
# 2) BirdNET output
# 3) Manual verifications done for any species within the Olympic project

# To the Sampling Effort data, we add spatial metadata provided by OLYM.

library(DBI)
library(RSQLite)
library(tidyverse)
library(ggpubr)


# Sampling effort data ----------------------------------------------------

my_db_file <- "input/20240429-123-audio_metadata-database.db"
mydb <- dbConnect(RSQLite::SQLite(), my_db_file)

# to see what tables are in the database:
dbListTables(mydb)
# "flatdetectionmodel" = all birdnet detections
# "flatverificationmodel" is the corresponding table from Audiodash
# "projectfilesmodel" is for spatial metadata

# read in the table
sample_data <- dbReadTable(mydb, "projectfilesmodel")

# Once finished interacting with the db always close the connection
dbDisconnect(mydb)

head(sample_data)
# Each row corresponds to a file in the dataset. We need these because they represent the actual inputs, even if birdnet didn't identify anything in them

unique(sample_data$location_name)
# need to fix the weirdly named one (as we did in the birdnet output)

sample_data$location_name[sample_data$location_name=="OLY_41511-1"] <- "41511-1"
unique(sample_data$location_name) # we good
n_distinct(sample_data$location_name) # 191 sites

sample_data %>% dplyr::select(project_name, location_name, duration_seconds, captured_local_date, captured_local_time) %>% group_by(location_name, captured_local_date) %>% summarise(nfile=n(), nsec = sum(duration_seconds))

# BirdNET detections ------------------------------------------------------


my_db_file <- "input/20231127-76-detections-database.db"
mydb <- dbConnect(RSQLite::SQLite(), my_db_file)
det_data <- dbReadTable(mydb, "flatdetectionmodel")
dbDisconnect(mydb)

head(det_data)
unique(det_data$location)
det_data$location[det_data$location=="OLY_41511-1"] <- "41511-1"
n_distinct(det_data$location)



# Manual Verifications ----------------------------------------------------

my_db_file <- "input/20240429-122-verifications-database.db"
mydb <- dbConnect(RSQLite::SQLite(), my_db_file)

# If you want to see what tables are in
# the database you can use this
dbListTables(mydb)
# This "verif" table is where everything
# lives for the verifications from colab;
# "flatverificationmodel" is the corresponding table from Audiodash
#"projectfilesmodel" is for spatial metadata

# To read in the table to R use this function
verif_data <- dbReadTable(mydb, "flatverificationmodel")
# Once finished interacting with the db
# always close the connection
dbDisconnect(mydb)


# Spatial Metadata --------------------------------------------------------


loc_data <- read_csv("input/OLYM_coords_UTM_latlong.csv")

head(loc_data)


# subsample to verified species -------------------------------------------

unique(verif_data$focal_sound_class)

# preliminary list of species we will try to model. 
spp_list <- c("American Robin", "Hermit Thrush", "Hammond's Flycatcher", "Olive-sided Flycatcher", "Varied Thrush", "Pacific Wren", "Red Crossbill", "Evening Grosbeak", "Western Tanager", "Brown Creeper", "Pacific-slope Flycatcher", "Dark-eyed Junco", "Townsend's Warbler", "Chestnut-backed Chickadee")

mod_data <- det_data[det_data$common_name %in% spp_list,]

save(loc_data, sample_data, det_data, verif_data, mod_data, spp_list, file = "output/OLYM_inputs.RData")
