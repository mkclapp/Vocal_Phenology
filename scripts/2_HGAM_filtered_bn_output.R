# Script to format BirdNET data into format for HGAM

# 1- threshold the data to a 99% precision rate
# 2- threshold sites to ones w/ detections on 3+ days
# 3- convert counts to hits | trials
# 4- model hits using a Binomial HGAM with elevation as a predictor variable


# requires "mod_data" (BirdNET output for all modeled species);
# requires "thresholds" (calculated thresholds for precision) 

head(mod_data)
library(tidyverse)
library(mgcv)
library(gratia)


elev_data <- read_excel("input/IBP_OLYM_2021_Metadata.xlsx") %>% mutate(location=paste(Hexagon_ID, Station, sep="-"))
# many more locations than the ones we have
# we do not have lat-long for data sharing reasons but do have a shapefile of hexes. each location is 1 replicate station out of 5 within the hex.
unique(elev_data$location)
head(elev_data)

# 1- threshold data -------------------------------------------------------


# single species 

allheth <- mod_data %>% filter(common_name=="Hermit Thrush") %>% dplyr::select(project_name, location, captured_local_date, captured_local_time, duration_seconds, common_name, confidence, logit)

thresholds$heth[thresholds$threshold==0.99]

heth.f <- allheth %>% filter(confidence > thresholds$heth[thresholds$threshold==0.95]) 

heth.f %>% group_by(location) %>% summarise(nday = n_distinct(captured_local_date)) %>% arrange(nday) %>% filter(nday>2)

unique(heth.f$location)
head(heth.f)

allpsfl <- mod_data %>% filter(common_name=="Pacific-slope Flycatcher") %>% dplyr::select(project_name, location, captured_local_date, captured_local_time, duration_seconds, common_name, confidence, logit)

thresholds$psfl[thresholds$threshold==0.99]

psfl.f <- allpsfl %>% filter(confidence > thresholds$psfl[thresholds$threshold==0.95]) 

trials <- psfl.f %>% group_by(location, captured_local_date) %>% summarise(nhits = n())

# need to get 0s for hours with no hits
alldates <- sample_data[sample_data$location_name %in% trials$location,] %>% group_by(location_name, captured_local_date) %>% summarise(ntrials = sum(duration_seconds)/3) #%>% ggplot() + geom_point(aes(x=ymd(captured_local_date), y=ntrials)) + facet_wrap(~location_name)

n_distinct(alldates$location_name)
n_distinct(trials$location)

trials <- full_join(trials, alldates, by=c("location"="location_name", "captured_local_date"="captured_local_date"))
trials$nhits[is.na(trials$nhits)] <- 0
trials$captured_local_date <- ymd(trials$captured_local_date)
ggplot(trials) +
  geom_col(aes(x=captured_local_date, y=ntrials), alpha=0.5) + 
  geom_col(aes(x=captured_local_date, y=nhits), fill="pink") +
  facet_wrap(~location)

trials <- left_join(trials, elev_data)
trials$location <- as.factor(trials$location)

# HGAM --------------------------------------------------------------------

voc_days <- trials %>% group_by(location) %>% summarise(ndays = n_distinct(captured_local_date[nhits>0])) %>% arrange(desc(ndays)) %>% left_join(elev_data)

ggplot(voc_days) + 
  geom_point(aes(x=GIS_Elev_ft, y=ndays))

mod1 <- gam(nhits/ntrials ~ te(yday(captured_local_date), GIS_Elev_ft, bs=c("cc", "tp"), m=2) + 
                s(location, bs="re"),
              data=trials, method="REML", family = "binomial", weights = ntrials/mean(ntrials), knots=list(DateS=c(-91, 251))
)
summary(mod1)
draw(mod1)

# bin elevations
trials$ElevM <- trials$GIS_Elev_ft*3.28
trials$ElevBin <- ifelse(trials$ElevM < 650, "Low", ifelse(trials$ElevM > 1350, "High", "Med")) # from handbook
trials$ElevBin <- factor(trials$ElevBin, levels=c("Low", "Med", "High"))
table(trials$ElevBin)

trials$JDay <- yday(trials$captured_local_date)

# TODO: investigate how to specify weights... below (i think) uses option 2 in the Rhelp file: "As a numerical vector with values between 0 and 1, interpreted as the proportion of successful cases (with the total number of cases given by the weights)" and from mgcv: "If you want to reweight the contributions of each datum without changing the overall magnitude of the log likelihood, then you should normalize the weights (e.g. weights <- weights/mean(weights))."

# this model removes elevation as a predictor in the model (HETH are only detected in the "High" Elevation bins)
mod2 <- gam(nhits/ntrials ~ s(JDay, bs="cc", m=2) + 
                s(location, bs="re"),
              data=trials, method="REML", family = "binomial", weights = ntrials/mean(ntrials), knots=list(JDay=c(-91, 251))
)
summary(mod2)
draw(mod2)
