library(tidyverse)
library(readxl)
library(lubridate)
library(chron)
library(ggrepel)
library(ggpubr)
library(mgcv)
library(gratia)

# GAMs --------------------------------------------------------------------
head(herthr)
head(sitelist) # from OlympicSubsetEDA.R


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


# with an interaction between Elevation and JDay
# not sure if tp or cc is more appropriate for the two predictors
# also not entirely sure what m designates
# needs more data
# TODO: try to predict on a more complete set of values

mod2.2 <- gam(nVoc ~ te(yday(Date), ElevM, bs=c("cc", "tp"), m=2) + 
                s(SiteID, bs="re"),
            data=med, method="REML", family = "nb", knots=list(DateS=c(-91, 251))
)
summary(mod2.2)
draw(mod2.2)

mod3.2 <- gam(nVoc ~ te(yday(Date), ElevM, bs=c("cc", "tp"), m=2) + 
              s(SiteID, bs="re"),
            data=high, method="REML", family="nb")
summary(mod3.2)
draw(mod3.2)


# plot raw data
no0 <- heth95.byday %>% filter(nVoc>0, ElevBin != "Low")

ggplot(no0) +
  geom_point(aes(x=Date, y=log(nVoc), color=ElevM)) + 
  geom_smooth(aes(x=Date, y=log(nVoc))) +
 facet_wrap(~ElevBin) +
  theme_bw()
#  theme(legend.position = "none")


# species will vary in how many elevation bins they occupy, which will influence detection.

# to use ggplot2 to create modifiable output plots from mcgv:
#https://stackoverflow.com/questions/65527152/how-do-i-create-a-ggplot-in-r-from-a-non-linear-model-using-the-mgcv-package





# GLM ---------------------------------------------------------------------

head(heth95.byday)
heth95.byday$JDay <- yday(heth95.byday$Date)
heth95.byday$JDay.s <- as.numeric(scale(heth95.byday$JDay))

library(lme4)
library(sjPlot)
library(emmeans)
library(ggeffects)

foo <- glm(nVoc ~ JDay.s + I(JDay.s^2) + scale(ElevM), data=heth95.byday, family="poisson")
summary(foo)
foo2 <- glmer(nVoc ~ JDay.s + I(JDay.s^2) + I(JDay.s^3) + (1|ElevBin), data=heth95.byday, family="poisson")
summary(foo)
anova(foo, foo2)

plot(ggpredict(foo, terms=c("JDay.s[all]", "ElevBin")))


xDate <- seq(-2.8, 2, 0.01)

#Now we use the predict() function to create the model for all of the values of xweight.

yDate <- predict(foo, xDate)

plot(mtcars$wt, mtcars$vs, pch = 16, xlab = "WEIGHT (g)", ylab = "VS")