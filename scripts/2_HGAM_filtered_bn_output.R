# Script to format BirdNET data into format for HGAM

# 1- threshold the data to a 99% precision rate
# 2- threshold sites to ones w/ detections on 3+ days
# 3- convert counts to hits | trials
# 4- model hits using a Binomial HGAM with elevation as a predictor variable


# requires "mod_data" (BirdNET output for all modeled species);
# requires "thresholds" (calculated thresholds for precision) 

library(tidyverse)
library(mgcv)
library(gratia)
library(readxl)
library(ggpubr)

trials.0fill <- readRDS("output/trials.0fill.rds")
thresholds <- read_csv("output/precision_thresholds.csv")
elev_data <- read_excel("input/IBP_OLYM_2021_Metadata.xlsx") %>% mutate(location=paste(Hexagon_ID, Station, sep="-"))
# many more locations than the ones we have
# we do not have lat-long for data sharing reasons but do have a shapefile of hexes. each location is 1 replicate station out of 5 within the hex.
unique(elev_data$location)
head(elev_data)


# bin elevations
elev_data$ElevM <- elev_data$GIS_Elev_ft/3.28
elev_data$ElevBin <- ifelse(elev_data$ElevM < 650, "Low", ifelse(elev_data$ElevM > 1350, "High", "Mid")) # from handbook
elev_data$ElevBin <- factor(elev_data$ElevBin, levels=c("Low", "Mid", "High"))
elev_data$ElevLabs <- factor(ifelse(elev_data$ElevBin=="Low", "Low (<650m)", ifelse(elev_data$ElevBin=="Mid", "Mid (650-1350m)", "High (>1350m)")), levels = c("Low (<650m)", "Mid (650-1350m)", "High (>1350m)"))
table(elev_data$ElevLabs)


# 1- threshold data -------------------------------------------------------

unique(mod_data$common_name)
spp_codes <- data.frame(common_name = unique(mod_data$common_name), aou4 = c("psfl", "brcr", "towa", "cbch", "pawr", "evgr", "vath", "deju", "recr", "amro", "weta", "hafl", "osfl", "heth")) %>% arrange(common_name)

allspp <- mod_data %>% left_join(spp_codes) %>%
  dplyr::select(project_name, location, captured_local_date, captured_local_time, duration_seconds, common_name, aou4, confidence)

thresholds.l <- thresholds %>% pivot_longer(2:ncol(thresholds), names_to = "aou4", values_to = "confidence")

thresh.data <- allspp[NULL, ]
for (i in 1:dim(spp_codes)[1]) {
  tmp <- allspp[allspp$aou4 == spp_codes[i, 2], ]
  tmp <- tmp[tmp$confidence >= thresholds.l$confidence[which(thresholds.l$threshold == 0.99 & thresholds.l$aou4==tmp$aou4[1])], ] # using a 99% threshold
  thresh.data <- rbind(thresh.data, tmp)
}

trials <- thresh.data %>% group_by(location, captured_local_date, common_name, aou4) %>% summarise(nhits = n())

# need to get 0s for hours with no hits
alldates <- sample_data[sample_data$location_name %in% trials$location,] %>% group_by(location_name, captured_local_date) %>% summarise(ntrials = sum(duration_seconds)/3) #%>% ggplot() + geom_point(aes(x=ymd(captured_local_date), y=ntrials)) + facet_wrap(~location_name)

n_distinct(alldates$location_name)
n_distinct(trials$location)

# trials <- full_join(trials, alldates, by=c("location"="location_name", "captured_local_date"="captured_local_date"))
# trials$nhits[is.na(trials$nhits)] <- 0
# trials$captured_local_date <- ymd(trials$captured_local_date)

# trials <- left_join(trials, elev_data)
# trials$location <- as.factor(trials$location)

alldates$loc.date <- paste(alldates$location_name, alldates$captured_local_date, sep = "_")
trials.0fill <- trials
for (i in 1:dim(spp_codes)[1]) {
  tmp <- trials[trials$aou4 == spp_codes[i, 2], ]
  tmp$loc.date <- paste(tmp$location, tmp$captured_local_date, sep = "_")
  need.zero <- alldates[!alldates$loc.date %in% unique(tmp$loc.date), ]
  add.zeros <- data.frame(location =  need.zero$location_name,
                          captured_local_date = need.zero$captured_local_date,
                          common_name = tmp$common_name[1],
                          aou4 = tmp$aou4[1],
                          nhits = 0,
                          # ntrials = need.zero$ntrials,
                          loc.date = need.zero$loc.date)
  trials.0fill <- rbind(trials.0fill, add.zeros)
}
trials.0fill$loc.date <- paste(trials.0fill$location, trials.0fill$captured_local_date, sep = "_")

trials.0fill$ntrials <- round(alldates$ntrials[match(trials.0fill$loc.date, alldates$loc.date)], 0) # rounding number of trials

trials.0fill <- left_join(trials.0fill, elev_data) %>% dplyr::select(location, captured_local_date, common_name, aou4, nhits, ntrials, ElevM, ElevBin, ElevLabs)


trials.0fill$JDay <- yday(trials.0fill$captured_local_date)
summary(trials.0fill$JDay)

ggplot(trials.0fill) +
 # geom_col(aes(x=JDay, y=ntrials, group=common_name), alpha=0.5) + 
  geom_col(aes(x=JDay, y=nhits, fill=common_name)) +
  theme_pubclean() +
  labs(x="Date", y="# trials") +
  facet_wrap(~common_name, scales = "free_y")

trials.0fill$effort <- trials.0fill$nhits/trials.0fill$ntrials

trials.0fill %>% filter(aou4=="psfl") %>%
ggplot() +
  #geom_col(aes(x=JDay, y=ntrials, group=common_name), alpha=0.5) + 
  geom_line(aes(x=as.Date(captured_local_date), y=effort, group=location, color=location), alpha=0.9) +
  theme_pubclean() +
  labs(x="Date", y="prop. successes") +
  scale_x_date(date_breaks = "1 month", date_minor_breaks = "1 week",
               date_labels = "%B") +
  theme(legend.position = "none") +
  facet_wrap(~ElevBin, ncol=1, scales = "free_y")

table(trials.0fill$common_name, (trials.0fill$nhits>0)*1)
trials.0fill %>% group_by(ElevBin) %>% summarise(nsites=n_distinct(location))

saveRDS(trials.0fill, "output/trials.0fill.rds")


# 2- HGAM -----------------------------------------------------------------

# data for binomial family can be specified as a ratio of successes/total trials, providing the total # of trials in weights

glimpse(trials.0fill)
trials.0fill$location <- as.factor(trials.0fill$location)
trials.0fill$aou4 <- as.factor(trials.0fill$aou4)

# no elevation, modeling bins separately
mod1 <- gam(nhits/ntrials ~ s(JDay, bs="tp", k=15) + 
                s(location, bs="re"),
              data=trials.0fill[trials.0fill$common_name=="Varied Thrush" & trials.0fill$ElevBin=="Med (650-1350m)",], 
            method="REML", family = "binomial", weights = ntrials)

summary(mod1)
plot(mod1, trans = plogis)

# interaction between elevation and JDay
mod2 <- gam(nhits/ntrials ~ te(JDay, ElevM, bs=c("cc", "tp")) + 
              s(location, bs="re"),
            data=trials.0fill[trials.0fill$common_name=="Pacific Wren",], 
            method="REML", family = "binomial", weights = ntrials)

voc_days <- trials.0fill %>% group_by(location, ElevM, common_name) %>% summarise(ndays = n_distinct(JDay[nhits>0])) %>% arrange(desc(ndays)) 

ggplot(voc_days) + 
  geom_point(aes(x=ElevM, y=ndays, color=location), alpha=0.5) +
  facet_wrap(~common_name) +
  theme_pubclean() +
  theme(legend.position = "none")

ggplot(trials) +
  geom_point(aes(x=JDay, y=nhits))

# no elevation
mod2 <- gam(nhits/ntrials ~ s(JDay, bs="tp", m=2) + 
                s(location, bs="re"),
              data=trials, method="REML", family = "binomial", weights = ntrials, knots=list(JDay=c(-91, 251))
)
summary(mod2)
draw(mod2)
trials.0fill$aou4 <- as.factor(trials.0fill$aou4)


bird_modGI <- gam(count ∼ species + te(week, latitude, bs=c("cc", "tp"), k=c(10, 10), m=2) + te(week, latitude, by=species, bs= c("cc", "tp"), k=c(10, 10), m=1), data=bird_move, method="REML", family="poisson", knots=list(week=c(0, 52)))


bird_mod <- gam(nhits/ntrials~te(JDay, ElevM, bs=c("cc", "tp"), k=c(10, 10), m=2) + 
                    t2(JDay, ElevM, aou4, bs=c("cc", "tp", "re"), k=c(10, 10, 14), m=2, full=TRUE), 
                  data=trials.0fill, 
                  method="REML", 
                  family="binomial", weights=ntrials,
                  knots=list(JDay=c(90, 251)))


bird_modGS <- gam(nhits/ntrials~te(JDay, ElevM, bs=c("cc", "tp"), k=c(10, 10), m=2) + 
                    t2(JDay, ElevM, aou4, bs=c("cc", "tp", "re"), k=c(10, 10, 14), m=2, full=TRUE), 
                  data=trials.0fill, 
                  method="REML", 
                  family="binomial", weights=ntrials,
                  knots=list(JDay=c(90, 251)))
summary(bird_modGS)
draw(bird_modGS)


# OLD ---------------------------------------------------------------------


# HETH --------------------------------------------------------------------


allheth <- mod_data %>% filter(common_name=="Hermit Thrush") %>% dplyr::select(project_name, location, captured_local_date, captured_local_time, duration_seconds, common_name, confidence, logit)

heth.f <- allheth %>% filter(confidence > thresholds$heth[thresholds$threshold==0.95]) 

heth.f %>% group_by(location) %>% summarise(nday = n_distinct(captured_local_date)) %>% arrange(nday) %>% filter(nday>2)

unique(heth.f$location)
head(heth.f)


# PSFL --------------------------------------------------------------------


allpsfl <- mod_data %>% filter(common_name=="Pacific-slope Flycatcher") %>% dplyr::select(project_name, location, captured_local_date, captured_local_time, duration_seconds, common_name, confidence, logit)

thresholds$psfl[thresholds$threshold==0.99]

psfl.f <- allpsfl %>% filter(confidence > thresholds$psfl[thresholds$threshold==0.95]) 
