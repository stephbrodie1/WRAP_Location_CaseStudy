#This code uses output from SimulatedWorld_ROMS_TrophicInteraction_function to build a GAM
#It is just a simple example to check the output of the OM is sensible

#Simple trophic interaction
#Species A: distribution and abundance drivn by SST and Chl-a
#Species B: distribution and abundance drivn by SST and Species A
#EM for Species b: only have chl-a and temp as covariates. 

#----Directories----
#Set your working directory
setwd('~/PROJECTS/WRAP Location/')


#----Load Library & Function----
library(mgcv)
library(ggplot2)
library(viridis)
source('WRAP_Location_CaseStudy/SimulatedWorld_ROMS_Albacore_Function.R') #load ROMS simulation function
source('WRAP_Location_CaseStudy/SimulatedWorld_ROMS_Anchovy_Function.R') #load ROMS simulation function

#----Run Operating Model---- 
dir <- "~/Dropbox/WRAP Location^3/Rasters_2d_Spring/had/" #local directory where ROMS data is stored
# dat <- SimulateWorld_ROMS_Albacore(dir = dir ) #takes a few mins
# dat <- SimulateWorld_ROMS_Anchovy(dir = dir ) #takes a few mins
# logit.function <- function(x, a, b, c){
#   exp(a*x^2 + b*x + c)/(1+exp(a*x^2 + b*x + c))
# }
dat <- SimulateWorld_ROMS_Groundfish(dir = dir)
names(dat)[names(dat) == 'sst'] <- 'temp' #matching roms names. Quick temporary fix. 
dat_all <- dat
dat <- dat[dat$sample==1,]

#----Explore Data----
#Create dataframe with historical/forecast data
dat$log_abundance <- log(dat$abundance)
dat_hist <- dat[dat$year<=2010,]
dat_fcast <- dat[dat$year>2010,]

# # #Make some quick plots to explore the data
# # #All Years
par(mfrow=c(2,3))
plot(aggregate(suitability~year,dat,FUN="mean"),type="l", lwd=2, main="Mean Suitability",col="dark grey", ylab="")
lines(aggregate(suitability~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(pres~year,dat,FUN="mean"),type="l", lwd=2, main="Prevalence",col="dark grey", ylab="")
lines(aggregate(pres~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(pop_norm~year,dat,FUN="mean"),type="l", lwd=2, main="Prevalence",col="dark grey", ylab="")
lines(aggregate(pop_norm~year,dat[dat$year<=2010,],FUN="mean"),col="blue")

plot(aggregate(abundance~year,dat,FUN="sum"),type="l",  lwd=2,main="Total Biomass", col="dark grey", ylab="")
lines(aggregate(abundance~year,dat[dat$year<=2010,],FUN="sum"),col="blue")
plot(aggregate(abundance~year,dat_all,FUN="sum"),type="l",  lwd=2,main="Total Biomass", col="dark grey", ylab="")
lines(aggregate(abundance~year,dat_all[dat_all$year<=2010,],FUN="sum"),col="blue")

# plot(aggregate(abundance_2~year,dat,FUN="sum"),type="l",  lwd=2,main="Total Biomass", col="dark grey", ylab="")
# lines(aggregate(abundance_2~year,dat[dat$year<=2010,],FUN="sum"),col="blue")
# plot(aggregate(abundance_2~year,dat_all,FUN="sum"),type="l",  lwd=2,main="Total Biomass", col="dark grey", ylab="")
# lines(aggregate(abundance_2~year,dat_all[dat_all$year<=2010,],FUN="sum"),col="blue")


plot(aggregate(temp~year,dat,FUN="mean"),type="l",main="Mean Temperature", col="dark grey", ylab="")
lines(aggregate(temp~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(chla_surface~year,dat,FUN="mean"),type="l",main="Mean Surface Chl",col="dark grey", ylab="")
lines(aggregate(chla_surface~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(mld~year,dat,FUN="mean"),type="l",main="Mean MLD",col="dark grey")
lines(aggregate(mld~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(zoo_200~year,dat,FUN="mean"),type="l",main="Mean Zooplankton",col="dark grey", ylab="")
lines(aggregate(zoo_200~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(btemp~year,dat,FUN="mean"),type="l",main="Mean bottom Temperature", col="dark grey", ylab="")
lines(aggregate(btemp~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(O2~year,dat,FUN="mean"),type="l",main="Mean oxygen", col="dark grey", ylab="")
lines(aggregate(O2~year,dat[dat$year<=2010,],FUN="mean"),col="blue")


#----Build GAM Models----
hist(dat$abundance)
hist(log(dat$abundance_2[dat$abundance_2>0]))

# gam1.p <- gam(pres ~ s(temp,bs='gp') + s(mld, bs='gp') + s(chla_surface, bs="gp") , data=dat_hist, family=binomial)
# gam1.a <- gam(log_abundance ~ s(temp,bs='gp') + s(mld, bs='gp') + s(chla_surface, bs="gp") , data=dat_hist[dat_hist$abundance>0,], family=gaussian)
# # 
# gam2.p <- gam(pres ~ s(temp,bs='gp') + s(chla_surface, bs="gp") , data=dat_hist, family=binomial)
# gam2.a <- gam(abundance ~ s(temp,bs='gp') + s(chla_surface, bs="gp") , data=dat_hist[dat_hist$abundance>0,], family=gaussian)

#Anchovy
# gam1.p <- gam(pres ~ s(temp,bs='gp')+ s(z,bs="gp") + s(chla_surface, bs="gp") , data=dat_hist, family=binomial)
# gam1.a <- gam(log(abundance) ~ s(temp,bs='gp')+ s(z,bs="gp") +s(chla_surface, bs="gp") , data=dat_hist[dat_hist$abundance>0,], family=gaussian)

# gam2.p <- gam(pres ~ s(temp,bs='gp') + s(z,bs="gp")  , data=dat_hist, family=binomial)
# gam2.a <- gam(abundance ~ s(temp,bs='gp') + s(z,bs="gp")  , data=dat_hist[dat_hist$abundance>0,], family=gaussian)

#Groundfish
gam1.p <- gam(pres ~ s(btemp,bs='gp')+ s(z,bs="gp") + s(O2, bs="gp") , data=dat_hist, family=binomial)
gam1.a <- gam(log(abundance) ~ s(btemp,bs='gp')+ s(z,bs="gp") +s(O2, bs="gp") , data=dat_hist[dat_hist$abundance>0,], family=gaussian)




#----Explore GAM outputs-----
summary(gam1.p)
summary(gam1.a)
plot(gam1.p, scale=0, rug=TRUE)
plot(gam1.a, scale=0)

# summary(gam2.p)
# summary(gam2.a)
# # plot(gam2.p, scale=0)
# plot(gam2.a, scale=0)

#----Make predictions-----
#Historical predictions
dat_hist$gam1.p <- predict(gam1.p,dat_hist,type='response')
dat_hist$gam1.a <- predict(gam1.a,dat_hist,type="response")
dat_hist$gam1 <- dat_hist$gam1.p*exp(dat_hist$gam1.a)

#Future predictions
dat_fcast$gam1.p <- predict(gam1.p,dat_fcast,type='response')
dat_fcast$gam1.a <- predict(gam1.a,dat_fcast,type="response")
dat_fcast$gam1 <- dat_fcast$gam1.p*exp(dat_fcast$gam1.a)

# # #Historical predictions
# dat_hist$gam1.p <- predict(gam1.p,dat_hist,type='response')
# dat_hist$gam2.a <- predict(gam2.a,dat_hist,type="response")
# dat_hist$gam2 <- dat_hist$gam1.p*exp(dat_hist$gam2.a)
# # 
# # #Future predictions
# dat_fcast$gam1.p <- predict(gam1.p,dat_fcast,type='response')
# dat_fcast$gam2.a <- predict(gam2.a,dat_fcast,type="response")
# dat_fcast$gam2 <- dat_fcast$gam1.p*exp(dat_fcast$gam2.a)

#-----Performance Metrics----
#Time-series of abundance
#Quick plot: historical
plot(aggregate(abundance~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Biomass (mt)")
lines(aggregate(gam1~year,dat_hist,FUN="sum"),type="l",  lwd=1,ylab="Biomass (mt)",col="blue")

#Quick plot: future
plot(aggregate(abundance~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Biomass (mt)")
lines(aggregate(gam1~year,dat_fcast,FUN="sum"),type="l",  lwd=1,ylab="Biomass (mt)",col='blue')

# #Quick plot: historical
# plot(aggregate(abundance_2~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Biomass (mt)")
# lines(aggregate(gam2~year,dat_hist,FUN="sum"),type="l",  lwd=1,ylab="Biomass (mt)",col="red", ylim=c())
# 
# #Quick plot: future
# plot(aggregate(abundance_2~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Biomass (mt)")
# lines(aggregate(gam2~year,dat_fcast,FUN="sum"),type="l",  lwd=1,ylab="Biomass (mt)",col='red')
# # 
# # # 
# 
# t <- as.data.frame(aggregate(abundance~year,dat,FUN="sum")[2])

#---spatial predictions---

sst <- raster('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/had/sst_monthly/sst_monthly_had_SpringMean_2100.grd')
ild <- raster('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/had/ild_0.5C/ild_0.5C_had_SpringMean_2100.grd')
chla_surface <- raster('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/had/chl_surface/chl_surface_had_SpringMean_2100.grd')
chla_surface <- log(chla_surface)
z <- raster('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/gfdl/bottom_layer_depth.grd')
z <- z*-1
st <- stack(sst, chla_surface, z)
names(st) <- c("temp","chla_surface", "z")

test.p <- raster::predict(st,gam1.p, type='response')
plot(test.p)
test.a <- raster::predict(st,gam1.a, type='response')
plot(test.a)
test <- test.p * (test.a)
plot(test)

