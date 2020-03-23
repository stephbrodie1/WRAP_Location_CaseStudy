#This code uses output from SimulatedWorld_ROMS_TrophicInteraction_function to build a GAM
#Species A: distribution and abundance drivn by SST and Chl-a
#Species B: distribution and abundance drivn by SST and Species A
#EM for Species b: only have chl-a and temp as covariates. 

#Note only using GFDL for now

#----Directories----
#Set your working directory
setwd('~/PROJECTS/WRAP Location/')

#----Load Library & Function----
library(mgcv)
library(ggplot2)
library(viridis)
source('WRAP_Location_CaseStudy/SimulatedWorld_ROMS_TrophicInteraction_function.R') #load ROMS simulation function

#Set parameters for functions
abund_enviro <- "lnorm_low" #can be "lnorm_low" (SB); "lnorm_high" (EW); or "poisson" (JS)
PA_shape <- "logistic" #can be "logistic" (SB); "logistic_prev","linear" (JS)
dir <- "~/Dropbox/WRAP Location^3/Rasters_2d_monthly/" #directory where ROMS data is stored (on dropbox, email steph for access)

#Run this function
dat <- SimulateWorld_ROMS_TrophicInteraction(PA_shape = PA_shape, abund_enviro = abund_enviro, dir = dir ) #takes a few mins
colnames(dat)[1:2] <- c("Lon","Lat")
names(dat)[names(dat) == 'sst'] <- 'temp' #matching roms names. Quick temporary fix. 

#Create dataframe with historical/forecast data
dat_hist <- dat[dat$year<=2020,]
dat_fcast <- dat[dat$year>2020,]

#Make some quick plots to explore the data
#All Years
par(mfrow=c(2,2))
plot(aggregate(suitability~year,dat,FUN="mean"),type="l", lwd=2, ylab="Suitability",col="dark grey")
lines(aggregate(suitability~year,dat[dat$year<=2020,],FUN="mean"),col="blue")
plot(aggregate(abundance~year,dat,FUN="sum"),type="l",  lwd=2,ylab="Abundance", col="dark grey")
lines(aggregate(abundance~year,dat[dat$year<=2020,],FUN="sum"),col="blue")
plot(aggregate(temp~year,dat,FUN="min"),type="l",ylab="Temperature",ylim=c(8,30), col="dark grey")
plot(aggregate(chla~year,dat,FUN="mean"),type="l",ylab="Chl",col="dark grey")

#----Build GAM Models----
#Run if lognormal response was simulated
dat_hist$log_abundance <- log(dat_hist$abundance)

gam1.p <- gam(pres ~ s(temp,bs='gp') + s(chla,bs='gp') , data=dat_hist, family=binomial)
gam1.a <- gam(log_abundance ~ s(temp,bs='gp') + s(chla,bs='gp') , data=dat_hist[dat_hist$abundance>0,], family=gaussian)
summary(gam1.p)
summary(gam1.a)
plot(gam1.p)
plot(gam1.a)

#Historical predictions
dat_hist$gam1.p <- predict(gam1.p,dat_hist,type='response')
dat_hist$gam1.a <- predict(gam1.a,dat_hist,type="response")
dat_hist$gam1 <- dat_hist$gam1.p*exp(dat_hist$gam1.a)

#Future predictions
dat_fcast$gam1.p <- predict(gam1.p,dat_fcast,type='response')
dat_fcast$gam1.a <- predict(gam1.a,dat_fcast,type="response")
dat_fcast$gam1 <- dat_fcast$gam1.p*exp(dat_fcast$gam1.a)

#Quick plot: historical
plot(aggregate(abundance~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance")
lines(aggregate(gam1~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="blue")

#Quick plot: future
plot(aggregate(abundance~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance")
lines(aggregate(gam1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='blue')

  
  