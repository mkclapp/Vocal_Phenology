# This script digests the verification data available for all species in the Olympic dataset and calculates their precision curves and metrics. 
# It requires "OLYM_inputs.RData"

library(tidyverse)
library(ggplot2)
library(ggpubr)

# some functions to convert logits to confidence scores and vice versa
# confidence score (y) to logit score (x)
# x = log(y/(1-y))

cs2logit <- function(y) {
  x <- log(y/(1-y))
  return(x)
}

# logit score (x) to confidence score (y)
# y = 1/(1+exp(-x))

logit2cs <- function(x) {
  y <- 1/(1+exp(-x))
  return(y)
}

load("input/OLYM_inputs.RData")
rm(det_data) # remove this to save space if you only need the modeled species ("mod_data")

# TODO: missing code to resolve "needs review" samples. 
verif_bind <- left_join(verif_data, mod_data, by=c("detection_id"="id", "project_id"="project_id", "project_name"="project_name")) %>% filter(common_name %in% spp_list)
verif_bind$logit <- cs2logit(verif_bind$confidence)

# calculate precision curve for a single species


# PROOF verifications -----------------------------------------------------


# TODO: parallelize this across species? make its own script?

### AMERICAN ROBIN
amro.v <- verif_bind %>% filter(common_name=="American Robin" & verification_group_name != "OLY AMRO HETH REVIEW")
unique(amro.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
amro.v %>% group_by(detection_id) %>% summarise(nrev = n()) %>% filter(nrev>1) %>% View()
unique(amro.v$verification_group_name)

amro.v %>% filter(needs_review==1) %>% dplyr::select(detection_id, verification_group_name, user, is_species_present, needs_review, focal_sound_class, alternative_species_common_name, alternative_sound_class, common_name, confidence, comments) %>% View()

disagree <- amro.v %>% group_by(detection_id) %>% 
  summarise(nrev = n_distinct(user), disagree=sum(is_species_present)) %>% filter(disagree==1 & nrev==2) # 1 = if one of us said it was there and one of us said it wasn't
amro.v %>% filter(detection_id %in% disagree$detection_id) %>% dplyr::select(detection_id, verification_group_name, user, is_species_present, needs_review, focal_sound_class, alternative_species_common_name, alternative_sound_class, common_name, confidence, comments) %>% arrange(detection_id) %>% View()
# not sure how to code this resolution. for now, going to filter out everything marked "needs review" because there is a corresponding ID by the second reviewer, and for the ones that disagree, we want to keep the second reviewer's ID.
amro <- amro.v %>% filter(needs_review==0)

### VARIED THRUSH
vath.v <- verif_bind %>% filter(common_name=="Varied Thrush", user != "mclapp@birdpop.org")

unique(vath.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
vath.v %>% group_by(detection_id) %>% summarise(nrev = n()) %>% filter(nrev>1) %>% View()
unique(vath.v$verification_group_name)

#vath.v %>% filter(needs_review==1) %>% dplyr::select(detection_id, verification_group_name, user, is_species_present, needs_review, focal_sound_class, alternative_species_common_name, alternative_sound_class, common_name, confidence, comments) %>% View()

disagree <- vath.v %>% group_by(detection_id) %>% 
  summarise(nrev = n_distinct(user), disagree=sum(is_species_present)) %>% filter(disagree==1 & nrev==2) # 1 = if one of us said it was there and one of us said it wasn't
#vath.v %>% filter(detection_id %in% disagree$detection_id) %>% dplyr::select(detection_id, verification_group_name, user, is_species_present, needs_review, focal_sound_class, alternative_species_common_name, alternative_sound_class, common_name, confidence, comments) %>% arrange(detection_id) %>% View()
# not sure how to code this resolution. for now, going to filter out everything marked "needs review" because there is a corresponding ID by the second reviewer, and for the ones that disagree, we want to keep the second reviewer's ID.
vath <- vath.v %>% filter(needs_review==0)

### TOWNSEND'S WARBLER
towa.v <- verif_bind %>% filter(common_name=="Townsend's Warbler")

unique(towa.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
#towa.v %>% group_by(detection_id) %>% summarise(nrev = n()) %>% filter(nrev>1) %>% View()
towa <- towa.v # no needs review

### BROWN CREEPER
brcr.v <- verif_bind %>% filter(common_name=="Brown Creeper")

unique(brcr.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
#brcr.v %>% group_by(detection_id) %>% summarise(nrev = n()) %>% filter(nrev>1) %>% View()
brcr <- brcr.v

### WESTERN FLYCATCHER (psfl)
psfl.v <- verif_bind %>% filter(common_name=="Pacific-slope Flycatcher")

unique(psfl.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
psfl.v %>% filter(needs_review==1)

psfl <- psfl.v 

### HAMMOND'S FLYCATCHER
hafl.v <- verif_bind %>% filter(common_name=="Hammond's Flycatcher")

#unique(hafl.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
#hafl.v %>% group_by(detection_id) %>% summarise(nrev = n()) %>% filter(nrev>1) %>% View()
hafl <- hafl.v

### DARK-EYED JUNCO
deju.v <- verif_bind %>% filter(common_name=="Dark-eyed Junco")

unique(deju.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
deju.v %>% filter(needs_review==1)
deju <- deju.v

### HERMIT THRUSH
heth.v <- verif_bind %>% filter(common_name=="Hermit Thrush" & verification_group_name != "OLY AMRO HETH REVIEW")

unique(heth.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
#heth.v %>% group_by(detection_id) %>% summarise(nrev = n()) %>% filter(nrev>1) %>% View()

disagree <- heth.v %>% group_by(detection_id) %>% 
  summarise(nrev = n_distinct(user), disagree=sum(is_species_present)) %>% filter(disagree==1 & nrev==2) # 1 = if one of us said it was there and one of us said it wasn't
#heth.v %>% filter(detection_id %in% disagree$detection_id) %>% dplyr::select(detection_id, verification_group_name, user, is_species_present, needs_review, focal_sound_class, alternative_species_common_name, alternative_sound_class, common_name, confidence, comments) %>% arrange(detection_id) %>% View()
# not sure how to code this resolution. for now, going to filter out everything marked "needs review" because there is a corresponding ID by the second reviewer, and for the ones that disagree, we want to keep the second reviewer's ID.
heth <- heth.v %>% filter(needs_review==0)

### EVENING GROSBEAK
evgr.v <- verif_bind %>% filter(common_name=="Evening Grosbeak")

unique(evgr.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
#evgr.v %>% group_by(detection_id) %>% summarise(nrev = n()) %>% filter(nrev>1) %>% View()
evgr.v %>% filter(needs_review==1) # i remember verifying these
evgr <- evgr.v %>% filter(needs_review==0)

### OLIVE-SIDED FLYCATCHER

osfl.v <- verif_bind %>% filter(common_name=="Olive-sided Flycatcher", user != "mclapp@birdpop.org")

unique(osfl.v$user) # if more than 1, review discrepancies. rule of thumb should be that second reviewer's ID replaces first reviewer's ID where there are discrepancies
#osfl.v %>% group_by(detection_id) %>% summarise(nrev = n()) %>% filter(nrev>1) %>% View()
osfl.v %>% filter(needs_review==1) 

disagree <- osfl.v %>% group_by(detection_id) %>% 
  summarise(nrev = n_distinct(user), disagree=sum(is_species_present)) %>% filter(disagree==1 & nrev==2) # 1 = if one of us said it was there and one of us said it wasn't
#osfl.v %>% filter(detection_id %in% disagree$detection_id) %>% dplyr::select(detection_id, verification_group_name, user, is_species_present, needs_review, focal_sound_class, alternative_species_common_name, alternative_sound_class, common_name, confidence, comments) %>% arrange(detection_id) %>% View()
# not sure how to code this resolution. for now, going to filter out everything marked "needs review" because there is a corresponding ID by the second reviewer, and for the ones that disagree, we want to keep the second reviewer's ID.
osfl <- osfl.v %>% filter(needs_review!=1)


# MODEL precision ---------------------------------------------------------

# TODO: is there a way to propagate the uncertainty in THIS model into the estimates?

summary(m_amro <- glm(is_species_present ~ logit, data = amro, family=binomial()))
summary(m_vath <- glm(is_species_present ~ logit, data = vath, family=binomial()))
summary(m_towa <- glm(is_species_present ~ logit, data = towa, family=binomial()))
summary(m_brcr <- glm(is_species_present ~ logit, data = brcr, family=binomial()))
summary(m_psfl <- glm(is_species_present ~ logit, data = psfl, family=binomial()))
summary(m_deju <- glm(is_species_present ~ logit, data = deju, family=binomial()))
summary(m_hafl <- glm(is_species_present ~ logit, data = hafl, family=binomial()))
summary(m_heth <- glm(is_species_present ~ logit, data = heth, family=binomial()))
summary(m_osfl <- glm(is_species_present ~ logit, data = osfl, family=binomial()))
summary(m_evgr <- glm(is_species_present ~ logit, data = evgr, family=binomial()))
solvefor <- function(p, mod) {
  logodds = (log(p/(1-p)) - coef(mod)[1])/coef(mod)[2]
  Cscore = logit2cs(logodds)
}

thresholds <- data.frame(threshold=seq(0.500, 0.999, 0.01))

#solvefor(c(0.9, 0.95, 0.975, 0.99))
thresholds$amro <- solvefor(p = thresholds$threshold, mod=m_amro)
thresholds$vath <- solvefor(p = thresholds$threshold, mod=m_vath)
thresholds$towa <- solvefor(p = thresholds$threshold, mod=m_towa)
thresholds$brcr <- solvefor(p = thresholds$threshold, mod=m_brcr)
thresholds$psfl <- solvefor(p = thresholds$threshold, mod=m_psfl)
thresholds$deju <- solvefor(p = thresholds$threshold, mod=m_deju)
thresholds$hafl <- solvefor(p = thresholds$threshold, mod=m_hafl)
thresholds$heth <- solvefor(p = thresholds$threshold, mod=m_heth)
thresholds$osfl <- solvefor(p = thresholds$threshold, mod=m_osfl)
thresholds$evgr <- solvefor(p = thresholds$threshold, mod=m_evgr)

write_csv(thresholds, "./output/precision_thresholds.csv")


# GRAPH precision curves --------------------------------------------------

# TODO: iterate over all species and summarise 
# for now, just HETH as an example

newdat <- with(heth, 
               data.frame(logit = seq(min(heth$logit), max(heth$logit), length=200)))

preds <- predict(m_heth, newdata=newdat, type="link", se.fit = TRUE)

m1 <- m_heth

critval <- 1.96 ## approx 95% CI
upr <- preds$fit + (critval * preds$se.fit)
lwr <- preds$fit - (critval * preds$se.fit)
fit <- preds$fit

fit2 <- m1$family$linkinv(fit)
upr2 <- m1$family$linkinv(upr)
lwr2 <- m1$family$linkinv(lwr)

newdat$Conf <- logit2cs(newdat$logit)
newdat$fit <- fit2
newdat$lwr <- lwr2 
newdat$upr <- upr2 


ggplot(data=heth, mapping=aes(confidence, is_species_present)) + 
  geom_point(aes(color=factor(is_species_present))) +         
  geom_line(data=newdat, aes(x=Conf, y=fit)) +
  #  stat_smooth(method="glm", method.args=list(family=binomial)) + 
  geom_line(data=newdat, mapping=aes(x=Conf, y=upr2)) + 
  geom_line(data=newdat, mapping=aes(x=Conf, y=lwr2)) +
  geom_hline(aes(yintercept=0.9)) +
  scale_color_manual(values = c("darkred", "skyblue")) +
  labs(color="species present?")

# specificiity ~= precision
# sensitivity ~= recall

colnames(mod_data)
mod_data$logit <- qlogis(mod_data$confidence)
allheth <- mod_data %>% filter(common_name=="Hermit Thrush")
hist(allheth$logit,breaks=100)
hist(allheth$logit[allheth$logit>0.277], breaks=100)
# 95% HETH = 0.277

head(thresholds)

towa11 <- towa[towa$is_species_present==1,]
towa01 <- towa[towa$is_species_present==0,]
ggplot() +
  geom_histogram(data= towa11, aes(logit), fill="blue", alpha=0.5) +
  geom_histogram(data=towa01, aes(logit), fill="darkred", alpha=0.8)

truepos <- verif_bind[verif_bind$is_species_present==1,]
falsepos <- verif_bind[verif_bind$is_species_present==0,]

ggplot() +
  geom_histogram(data= truepos, aes(logit), fill="blue", alpha=0.5) +
  geom_histogram(data=falsepos, aes(logit), fill="darkred", alpha=0.8) + 
  facet_wrap(~common_name)
