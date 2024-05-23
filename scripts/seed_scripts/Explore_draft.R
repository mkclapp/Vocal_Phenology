# birdnet data

library(DBI)
library(RSQLite)
library(tidyverse)
library(mgcv)
library(gratia)

splist <- read_csv("../input/OLYM_specieslist_filtered_proofed.csv")

# Set your directory to where
# you want to create a connection
setwd("C:/Users/Mary/Documents/RProjects/OlympicNP/")

# Specify the name of your database
my_db_file <- "input/20231127-76-detections-database.db"

# Open connection to your database
mydb <- dbConnect(RSQLite::SQLite(), my_db_file)

# Bring in my verification information
det_data <- dbReadTable(mydb, "flatdetectionmodel")

# Close connection
dbDisconnect(mydb)



unique(det_data$location)
det_data$location[det_data$location=="OLY_41511-1"] <- "41511-1"

sitelist2 <- data.frame(SiteID=unique(det_data$location))
sitelist2$SiteID <-gsub('-', '_', sitelist2$SiteID)
sitelist2 <- sitelist2 %>% left_join(all)

sitelist2 %>% group_by(ElevBin) %>% summarise(n())

quantiles <- quantile(sitelist2$ElevM, probs = c(0, 0.25, 0.5, 0.75, 1))

sitelist2$quantile <- cut(sitelist2$ElevM, breaks = quantiles, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)

bnsp <- det_data %>% group_by(common_name) %>% summarise(nsite = n_distinct(location), ncall = n(), meanCS = mean(confidence), maxCS=max(confidence))
bnsp$insplist <- ifelse(bnsp$common_name %in% splist$Common.Name, "Y", "N")

setdiff(bnsp$common_name, splist$Common.Name) # in birdnet, not in OLYM species list
setdiff(splist$Common.Name, bnsp$common_name) # in species list, not in birdnet


hits95 <- det_data %>% filter(common_name=="Hermit Thrush", confidence>0.2)

glimpse(hits95)

hits95



bydayTOWA <- det_data %>% group_by(captured_local_date) %>% filter(common_name=="Townsend's Warbler") %>% summarise(none = n(), `0.95` = sum(confidence>0.87),`0.99` = sum(confidence>0.97))

#long format for plot
byday.l <- byday %>% pivot_longer(cols = none:`0.99`, names_to = "threshold", values_to = "nVoc")


ggplot(bydayTOWA) + 
  geom_point(aes(x=as.Date(captured_local_date), y=`0.95`), alpha=0.8) +
  theme_bw() +
  labs(y="raw counts of BirdNET hits for TOWA", x="") +
  theme(legend.position = "bottom", text = element_text(size=10, angle = 60, hjust=1)) +
  scale_x_date(date_breaks = '1 week', date_labels = "%b %d") 

bydayHETH <- det_data %>% group_by(captured_local_date) %>% filter(common_name=="Hermit Thrush") %>% summarise(none = n(), `0.95` = sum(confidence>0.29),`0.99` = sum(confidence>0.63))

#long format for plot
byday.l <- byday %>% pivot_longer(cols = none:`0.99`, names_to = "threshold", values_to = "nVoc")


ggplot(bydayHETH) + 
  geom_point(aes(x=as.Date(captured_local_date), y=`0.95`), alpha=0.8) +
  theme_bw() +
  labs(y="raw counts of BirdNET hits for HETH", x="") +
  theme(legend.position = "bottom", text = element_text(size=10, angle = 60, hjust=1)) +
  scale_x_date(date_breaks = '1 week', date_labels = "%b %d") 


# TOWA

bystation <- det_data %>%  group_by(captured_local_date, location) %>% filter(common_name=="Townsend's Warbler") %>% summarise(none = n(), nVoc95 = sum(confidence>0.87), nVoc99 = sum(confidence>0.97))

unique(bystation$location) # make sure they all look reasonable

bystation$location<- gsub("-", "_", bystation$location)

bystation <- left_join(bystation, sitelist2, by=c("location"="SiteID"))

class(bystation$ElevBin)
bystation$ElevBin <- factor(bystation$ElevBin, levels=c("Low", "Med", "High"))

bystation %>% group_by(ElevBin) %>% summarise(n_distinct(location))
bystation$Date <- ymd(bystation$captured_local_date)

bystation$lomid <- cut(bystation$ElevM, breaks = c(0,650,1600), labels=c("low", "mid"))

bystation %>% 
  ggplot() + 
  geom_point(aes(x=as.Date(captured_local_date), y=nVoc95), alpha=0.8) +
  theme_bw() +
  labs(y="raw counts of BirdNET hits for HETH") + 
  facet_wrap(~ElevBin) +
  scale_x_date(date_breaks = '1 month', date_labels = "%b %d") 
bystation$location <- as.factor(bystation$location)


low95 <- bystation %>% filter(lomid=="low") %>% mutate(nVoc = nVoc95, JDay = yday(ymd(captured_local_date)))
low95$location <- as.factor(low95$location)

summary(modlo <- gam(nVoc ~ s(JDay, bs="tp") + s(location, bs="re"), 
                    data=low95, 
                    family="poisson",
                    knots = list(JDay=c(91, 252))))
draw(modlo)

med95 <- bystation %>% filter(lomid=="mid") %>% mutate(nVoc = nVoc95, JDay = yday(ymd(captured_local_date)))
med95$location <- as.factor(med95$location) # need for mgcv
med95.no0 <- med95 %>% filter(nVoc > 1) %>% dplyr::select(location, JDay, nVoc)

summary(modmid <- gam(nVoc ~ s(JDay, bs="tp") + s(location, bs="re"), 
                    data=med95, 
                    family="poisson",
                    knots = list(JDay=c(91, 252))))
draw(modmid)

mod2 <- gam(nVoc95 ~ te(yday(captured_local_date), ElevM, bs=c("cc", "tp")) + s(location, bs="re"),
                    data=bystation, 
                    method="REML", 
                    family="poisson",
                    knots=list(JDay=c(91, 252))) 
summary(mod2)
draw(mod2)

mod2.nb <- gam(nVoc95 ~ te(yday(captured_local_date), ElevM, bs=c("cc", "tp")) + s(location, bs="re"),
            data=bystation, 
            method="REML", 
            family="nb",
            knots=list(JDay=c(91, 252))) 
summary(mod2.nb)
draw(mod2.nb)

loc0 <- bystation %>% group_by(location) %>% summarise(nday = n_distinct(captured_local_date[nVoc95>0])) %>% filter(nday>0)
bystation.no0 <- bystation %>% filter(location %in% loc0)

mod3 <- gam(nVoc95 ~ s(yday(captured_local_date), bs="cc") + s(ElevM, bs="tp") + s(location, bs="re"),
                    data=bystation, 
                    method="REML", 
                    family="poisson",
                    knots=list(JDay=c(89, 257))) 
summary(mod3)
draw(mod3)

mod3.nb <- gam(nVoc95 ~ s(yday(captured_local_date), bs="cc") + s(ElevM, bs="tp") + s(location, bs="re"),
            data=bystation, 
            method="REML", 
            family="nb",
            knots=list(JDay=c(89, 257))) 
summary(mod3.nb)
draw(mod3)

###################PREDICT on new data

newdates <- data.frame(JDay=seq(91, 250, 1), SiteID="11111_1")
pred2 <- predict.gam(mod2, newdates[1], exclude="s(SiteID)",newdata.guaranteed=TRUE, se.fit=T)
pred2$JDay <- newdates$JDay
pred2 <- as.data.frame(pred2)
pred3 <- predict.gam(mod3, newdates[1], exclude="s(SiteID)",newdata.guaranteed=TRUE, se.fit=T)
pred3$JDay <- newdates$JDay
pred3 <- as.data.frame(pred3)
ggplot() +
  geom_line(data=pred3, aes(x=JDay, y=fit), color="darkgreen") +
  geom_ribbon(data=pred3, aes(x=JDay, ymax=fit+se.fit, ymin=fit-se.fit),alpha=0.2, fill="darkgreen") +
  geom_line(data=pred2, aes(x=JDay, y=fit), color="blue") +
  geom_ribbon(data=pred2, aes(x=JDay, ymax=fit+se.fit, ymin=fit-se.fit),alpha=0.2, fill="blue") +
  theme_bw()



# FROM OTHER SCRIPTS
det_data

n_distinct(det_data$location)
herthr <- det_data %>% filter(common_name=="Hermit Thrush")
herthr$logit <- cs2logit(herthr$confidence)
# based on calculating false + thresholds from verifications
byday <- herthr %>% group_by(captured_local_date) %>% summarise(none = n(), `0.95` = sum(confidence>0.20), `0.99` = sum(confidence>0.44))

#long format for plot
byday.l <- byday %>% pivot_longer(cols = none:`0.99`, names_to = "threshold", values_to = "nVoc")
# check against below to see sum is working correctly
# herthr %>% group_by(Date, HexStn) %>% filter(logit > 0.476) %>%  summarise(n95 = n())

byday.l$threshold <- factor(byday.l$threshold, levels=c("none", "0.90", "0.95", "0.99"))

hethday.p <- ggplot(byday.l) + 
  geom_point(aes(x=captured_local_date, y=nVoc, color=threshold), size=1.2, alpha=0.8) +
  geom_line(aes(x=captured_local_date, y=nVoc, color=threshold), size=1.2, alpha=0.8) +
  theme_bw() +
  labs(y="raw counts of BirdNET hits for HETH", color="Precision threshold") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d() 

byday.l %>% filter(threshold==0.95) %>%
ggplot() + 
  geom_point(aes(x=captured_local_date, y=nVoc), size=1.2, alpha=0.8) +
  geom_line(aes(x=captured_local_date, y=nVoc), size=1.2, alpha=0.8) +
  theme_bw() +
  labs(y="raw counts of BirdNET hits for HETH", color="Precision threshold") +
  theme(legend.position = "bottom") #+
 # scale_color_viridis_d() 


bystation <- herthr %>% group_by(Date, HexStn) %>% summarise(none = n(), `0.90` = sum(Confidence>0.231), `0.95` = sum(Confidence>0.476), `0.99` = sum(Confidence>0.91))

# first rule: stations must have >3 detections at >95% threshold per day

0.05*0.05*0.05 #0.000125 = probability of all day-detections being entirely false-positive if the number of day-detections = 3 and the threshold is 0.95% "confidence"-- in other words (for this species), the probability of 3 samples with CS > 0.476 ALL being false positive.

0.1*0.1*0.1 #0.001

statf <- bystation %>% group_by(HexStn) %>% filter(`0.95` > 3) %>% summarise(day1=min(Date), ndays = n_distinct(Date)) %>% arrange(day1)

bystation.l <- bystation %>% pivot_longer(cols = none:`0.99`, names_to = "threshold", values_to = "nVoc") %>% 
  group_by(Date, threshold) %>% summarise(nStation = n_distinct(HexStn[which(nVoc>0)]))
bystation.l$threshold <- factor(bystation.l$threshold, levels=c("none", "0.90", "0.95", "0.99"))

ggplot(bystation.l) +
  geom_line(aes(x=Date, y=nStation, color=threshold), linewidth=1.5, alpha=0.6) +
  theme_bw() +
  labs(y="# of stations with BirdNET hits for HETH", color="Precision threshold") +
  theme(legend.position = "bottom") +
  scale_color_viridis_d()

nstat <- bystation %>% pivot_longer(cols = none:`0.99`, names_to = "threshold", values_to = "nVoc") %>% 
  group_by(threshold) %>% summarise(nStation = n_distinct(HexStn[which(nVoc>0)]))
nstat$threshold <- factor(nstat$threshold, levels=c("none", "0.90", "0.95", "0.99"))

bystation.p <- ggplot(nstat) +
  geom_point(aes(x=threshold, y=nStation), size=7) +
  geom_point(aes(x=threshold, y=nStation, color=threshold), size=5) +
  theme_bw() +
  labs(y="# of stations with BirdNET hits for HETH", fill="Precision threshold") +
  #  scale_x_discrete(expand = c(0.4,0.4)) +
  theme(legend.position = "bottom", plot.margin = unit(c(.5,.5,.5,.5), 'cm'), text = element_text( size=24)) +
  scale_color_viridis_d()

heth.bn.p <- ggarrange(hethday.p, bystation.p, 
                       widths=c(2,1),
                       common.legend = T)


# 0.476 is the threshold for 95% precision
summary(herthr %>% group_by(Date, SiteID) %>% summarise(nvoc=sum(Confidence>0.476)))
summary(herthr %>% group_by(Date, SiteID) %>% summarise(nvoc=n()))

str(heth95.byday <- herthr  %>%
      group_by(Date, SiteID) %>% summarise(nVoc=sum(Confidence>0.476), prob01 = round((1-0.95)^nVoc, digits=4)) %>%
      left_join(sitelist) %>% dplyr::select(Date, SiteID, ElevM, nVoc, prob01))

visits %>% filter(SiteID=="42145_3")

dethist_long <- visits %>% left_join(heth95.byday)

# -------------------------------------------------------------------------


hist(heth95.byday$nVoc, breaks=100)

heth95.byday[is.na(heth95.byday$ElevM==T),]
heth95.byday[is.na(heth95.byday$Date==T),]

herthr %>% filter(SiteID=="41511_1")

heth95.byday$ElevS <- as.numeric(scale(heth95.byday$ElevM))
heth95.byday$DateS <- as.numeric(scale(heth95.byday$Date))
heth95.byday$JDay <- yday(heth95.byday$Date)
heth95.byday$ElevBin <- ifelse(heth95.byday$ElevM < 650, "Low", ifelse(heth95.byday$ElevM > 1350, "High", "Med")) # from handbook
heth95.byday$ElevBin <- factor(heth95.byday$ElevBin, levels=c("Low", "Med", "High"))
heth95.byday$SiteID <- as.factor(heth95.byday$SiteID)
med <- heth95.byday %>% filter(ElevBin=="Med")
high <- heth95.byday %>% filter(ElevBin=="High")
low <- heth95.byday %>% filter(ElevBin=="Low")



# without random effect for site
summary(m2 <- gam(nVoc ~ s(JDay, bs="tp", k=15), data=med, method="REML", family="nb"
                  #,knots=list(JDay=c(91, 251)))
))
draw(m2)
est2 <- as.data.frame(evaluate_smooth(m2,smooth = "JDay"))

summary(m3 <- gam(nVoc ~ s(JDay, bs="tp", k=15), data=high, method="REML", family="nb"
                  #,knots=list(JDay=c(91, 251)))
))
draw(m3)
est3 <- as.data.frame(evaluate_smooth(m3,smooth = "JDay"))


est2$JDay[max(est2$est)]
est3$JDay[max(est3$est)]

# without random effects for site
ggplot() + 
  geom_line(data=est2, aes(x=JDay, y=est), linewidth=1.0, color="blue") +
  geom_ribbon(data=est2, aes(x=JDay, ymax=est+se, ymin=est-se), alpha=0.20, fill="blue") +
  geom_line(data=est3, aes(x=JDay, y=est), linewidth=1.0, color="darkgreen") +
  geom_ribbon(data=est3, aes(x=JDay, ymax=est+se, ymin=est-se), alpha=0.20, fill="darkgreen") +
  geom_hline(aes(yintercept=0)) +
  theme_bw() +
  labs(x="Julian Day", y="effect of Date on vocal activity")

# with random effect for SiteID

summary(mod2 <- gam(nVoc ~ s(JDay, bs="tp", k=15) + s(SiteID, bs="re"), data=med, method="REML", family="nb"
                    #,knots=list(JDay=c(91, 251)))
))
draw(mod2)

summary(mod3 <- gam(nVoc ~ s(JDay, bs="tp", k=15) + s(SiteID, bs="re"),
                    data = high, method = "REML", family="nb"))
draw(mod3)

newdates <- data.frame(JDay=seq(91, 250, 1), SiteID="11111_1")
pred2 <- predict.gam(mod2, newdates[1], exclude="s(SiteID)",newdata.guaranteed=TRUE, se.fit=T)
pred2$JDay <- newdates$JDay
pred2 <- as.data.frame(pred2)
pred3 <- predict.gam(mod3, newdates[1], exclude="s(SiteID)",newdata.guaranteed=TRUE, se.fit=T)
pred3$JDay <- newdates$JDay
pred3 <- as.data.frame(pred3)
ggplot() +
  geom_line(data=pred3, aes(x=JDay, y=fit), color="darkgreen") +
  geom_ribbon(data=pred3, aes(x=JDay, ymax=fit+se.fit, ymin=fit-se.fit),alpha=0.2, fill="darkgreen") +
  geom_line(data=pred2, aes(x=JDay, y=fit), color="blue") +
  geom_ribbon(data=pred2, aes(x=JDay, ymax=fit+se.fit, ymin=fit-se.fit),alpha=0.2, fill="blue") +
  theme_bw()

