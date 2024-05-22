# 

# requires "mod_data" (BirdNET output for all modeled species);
# requires "thresholds" (calculated thresholds for precision) 

head(mod_data)
library(tidyverse)
library(mgcv)
library(gratia)