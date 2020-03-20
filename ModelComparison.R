#Code for Location WRAP Workshop
#Compare models using simulated data

#----Directories----
#Set your working directory
setwd('~/PROJECTS/WRAP Location/')
dir.create('Sim1') #data and model outputs will be saved locally here
Sim1 <- ('Sim1/')

#----Load Library & Function----
library(dismo)
library(gbm)
library(mgcv)
library(ggplot2)
library(viridis)
library(BBmisc)
library(neuralnet)
source('WRAP_Location_CaseStudy/SimulatedWorld_Function.R') #load simulation function
source('WRAP_Location_CaseStudy/SimulatedWorld_ROMS_Function.R') #load ROMS simulation function

#-----Simulate data----
#Run a function that simulates species distribution and abundance
#IMPORTANT:
#There are two simulation functions that can be built. One uses ROMS (SimulatedWorld_ROMS_Function.R), one uses simulated environmental data (SimulatedWorld_Function.R):

#Set parameters for functions
abund_enviro <- "lnorm_low" #can be "lnorm_low" (SB); "lnorm_high" (EW); or "poisson" (JS)
PA_shape <- "logistic" #can be "logistic" (SB); "logistic_prev","linear" (JS)
temp_spatial <- "matern" #can be "simple" (SB); or "matern" (EW)
temp_diff <- c(1,4,3,7) #specifies min and max temps at year 1 and year 100 (e.g. temp_diff=c(1,3,5,7) means year 1 varies from 1-3C and year 100 from 5-7C). For non-ROMS data. 
dir <- "~/Dropbox/WRAP Location^3/Rasters_2d_monthly/" #directory where ROMS data is stored (on dropbox, email steph for access)

#Run this function
dat <- SimulateWorld_ROMS(PA_shape = PA_shape, abund_enviro = abund_enviro, dir = dir ) #takes a few mins
#OR this function
dat <- SimulateWorld(temp_diff = temp_diff,  temp_spatial = temp_spatial, PA_shape = PA_shape, abund_enviro = abund_enviro) #takes a few minutes

#make headers consistent (Steph needs to update functions to fix this)
colnames(dat)[1:2] <- c("Lon","Lat")
names(dat)[names(dat) == 'sst'] <- 'temp' #matching roms names. Quick temporary fix. 
names(dat)[names(dat) == 'presabs'] <- 'pres' #matching roms names. Quick temporary fix. 

#Save data
saveRDS(dat, paste0(Sim1,'Sim1.rds')) #save data 
# dat <- readRDS(paste0(Sim1,'Sim1.rds')) #load in data if needed

#Create dataframe with historical/forecast data
dat_hist <- dat[dat$year<=2020,]
dat_fcast <- dat[dat$year>2020,]

#Make some quick plots to explore the data
#All Years
par(mfrow=c(2,2))
plot(aggregate(suitability~year,dat,FUN="mean"),type="l", lwd=2, ylab="Suitability",col="dark grey")
lines(aggregate(suitability~year,dat[dat$year<=2020,],FUN="mean"),col="blue")
plot(aggregate(pres~year,dat,FUN="mean"),type="l", lwd=2,ylab="Presence",col="dark grey")
lines(aggregate(pres~year,dat[dat$year<=2020,],FUN="mean"),col="blue")
plot(aggregate(abundance~year,dat,FUN="sum"),type="l",  lwd=2,ylab="Abundance", col="dark grey")
lines(aggregate(abundance~year,dat[dat$year<=2020,],FUN="sum"),col="blue")
plot(aggregate(temp~year,dat,FUN="min"),type="l",ylab="Temperature",ylim=c(8,30), col="dark grey")
lines(aggregate(temp~year,dat,FUN="max"),type="l",col="dark grey")
lines(aggregate(temp~year,dat,FUN="mean"),type="l")

#----Optional: Sampling Program----
#Optional and up for discussion: subsetting data to simulation imperfect data collection
#JS: randomly subset the data - simulates ~imperfect data collection
# num_obs <- 500
# dat_histc <- dat_hist[sample(1:nrow(dat_hist), num_obs, replace=F),]

#SB: I'm not completely convinced we need to add another level of complexity here
#SB: I see value in testing the effect of different sampling programs, but might be outside the scope. 

#----Build GAM Models----
#Run if lognormal response was simulated
if (abund_enviro == "lnorm_low" | abund_enviro == "lnorm_high"){
  dat_hist$log_abundance <- log(dat_hist$abundance)
  #SB: simple gams with gaussian process smooths
  gam1.p <- gam(pres ~ s(temp,bs='gp') , data=dat_hist, family=binomial)
  gam1.a <- gam(log_abundance ~ s(temp,bs='gp')  , data=dat_hist[dat_hist$abundance>0,], family=gaussian)
  summary(gam1.p)
  summary(gam1.a)
  plot(gam1.p)
  plot(gam1.a)
}

#Run if poisson response was simulated
if (abund_enviro == "poisson"){
  # JS: default smoothness
  gam2.a <- gam(round(abundance) ~ s(temp), data=dat_hist, family=poisson)
  #JS: gam with restricted smoothness
  gam3.a <- gam(round(abundance) ~ s(temp, k=4), data=dat_hist, family=poisson)
  summary(gam2.a)
  summary(gam3.a)
  plot(gam2.a)
  plot(gam3.a)
}


#----Boosted Regression Tree----
#Make sure >1000 trees fitted

#function to extract explained deviance from BRT
dev_eval=function(model_object){
  null <- model_object$self.statistics$mean.null
  res <- model_object$self.statistics$mean.resid
  dev=((null - res)/null)*100 
  return(dev)
}

#Run if lognormal response was simulated
if (abund_enviro == "lnorm_low" | abund_enviro == "lnorm_high"){
  brt1.a <- gbm.step(data=dat_hist[dat_hist$abundance>0,], gbm.x = c("temp"),gbm.y = 'log_abundance',family = "gaussian",tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6)
  brt1.p <- gbm.step(data=dat_hist, gbm.x = c("temp"),gbm.y = 'pres',family = "bernoulli",tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6)
  # saveRDS(brt1.a,paste0(Sim1,'BRT_Sim1_lognorm.rds'))
  # saveRDS(brt1.p,paste0(Sim1,'BRT_Sim1_binom.rds'))
  # brt1.a <- readRDS(paste0(Sim1,'BRT_Sim1_lognorm.rds'))#read in model object if required
  # brt1.p <- readRDS(paste0(Sim1,'BRT_Sim1_binom.rds'))
  dev_eval(brt1.p)
  dev_eval(brt1.a)
  plot(brt1.p)
  plot(brt1.a)
  
}

#Run if poisson response was simulated
if (abund_enviro == "poisson"){
  dat_hist$abundance_round <- round(dat_hist$abundance) #rounded only needed for Poisson
  brt2.a <- gbm.step(data=dat_hist, gbm.x = c("temp"),gbm.y = 'abundance_round',family = "poisson",tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6)
  dev_eval(brt2.a)
  plot(brt2.a)
}


#----Build MLP Models----
#Code from Barb Muhling
# 3 neurons in the hidden layer is good for this particular simulation. If it gets much more complicated, can try 
# different values. But training time increases exponentially with additional neurons...
# Important: variables must be normalized to a common scale before training the MLPs (doesn't matter for current 1 variable SDMs,
# but will if we add more predictors later)

#Run if lognormal response was simulated
if (abund_enviro == "lnorm_low" | abund_enviro == "lnorm_high"){
  dat_norm <- dat
  dat_norm$temp <- BBmisc::normalize(dat$temp) #covariates must be normalised for MLP
  dat_norm$log_abundance <- log(dat_norm$abundance)
  dat_hist_norm  <- subset(dat_norm, dat_norm$year <= 2020)
  dat_fcast_norm <- subset(dat_norm, dat_norm$year >  2020)
  
  # Train the two MLPs
  mlp1.p <- neuralnet(pres ~ temp, data = dat_hist_norm, 
                      hidden = c(3), linear.output = F, algorithm = "rprop+", threshold = 0.2)
  mlp1.a <- neuralnet(log_abundance ~ temp, data = subset(dat_hist_norm, !is.infinite(dat_hist_norm$log_abundance)), 
                      hidden = c(3), linear.output = T, algorithm = "rprop+", threshold = 0.2)
  # mlp1.a sometimes doesn't converge, or predicts essentially a constant value (the two main ways that MLPs can fail)
  # This code below asks it to keep trying to build the MLP until it succeeds
  # If it still won't converge for some datasets, try increasing the threshold 
  # For this simulation, mlp1.p always converges easily, but you could also institute this process for that MLP too if needed:
  # test.mlp1.a <- try(predict(mlp1.a, dat_hist_norm))
  # while(inherits(test.mlp1.a, 'try-error') | diff(range(test.mlp1.a < 0.1))) { 
  #   print("mlp1.a failed")
  #   mlp1.a <- neuralnet(log_abundance ~ temp, data = subset(dat_hist_norm, !is.infinite(dat_hist_norm$log_abundance)), 
  #                       hidden = c(3), linear.output = T, algorithm = "rprop+", threshold = 0.2)
  #   test.mlp1.a <- try(predict(mlp1.a, dat_hist_norm))
  # }
  
  # saveRDS(mlp1.p, paste0(Sim1,'MLP_Sim1_binom.rds'))
  # saveRDS(mlp1.a, paste0(Sim1,'MLP_Sim1_lognorm.rds'))
  # mlp1.p <- readRDS( paste0(Sim1,'MLP_Sim1_binom.rds')) #read in model object if required
  # mlp1.a <- readRDS( paste0(Sim1,'MLP_Sim1_lognorm.rds'))
  
  # Quick code to show partial relationships, the "plot" fn for mlps doesn't give the same plots as it does for GAMs/BRTs
  dat_hist_norm$mlpPres <- predict(mlp1.p, dat_hist_norm)
  dat_hist_norm$mlpAbun <- predict(mlp1.a, dat_hist_norm)
  plot(dat_hist$temp, dat_hist_norm$mlpPres) 
  plot(dat_hist$temp, dat_hist_norm$mlpAbun)
}

#Run if poisson response was simulated
if (abund_enviro == "poisson"){
  dat_norm <- dat
  dat_norm$temp <- BBmisc::normalize(dat$temp) #covariates must be normalised for MLP
  dat_norm$abundance_round <- round(dat_norm$abundance)
  dat_hist_norm  <- subset(dat_norm, dat_norm$year <= 2020)
  dat_fcast_norm <- subset(dat_norm, dat_norm$year >  2020)
  
  # Train the two MLPs
  mlp2.a <- neuralnet(abundance_round ~ temp, data = dat_hist_norm, 
                      hidden = c(3), linear.output = T, algorithm = "rprop+", threshold = 0.2)
  
  # Quick code to show partial relationships, the "plot" fn for mlps doesn't give the same plots as it does for GAMs/BRTs
  dat_hist_norm$mlpAbun <- predict(mlp2.a, dat_hist_norm)
  plot(dat_hist$temp, dat_hist_norm$mlpAbun)
  
}

#----Optional: Exploration of poor estimation of upper thermal limit, and methods to constrain it----
#Written by JS
# gam with restricted smoothness and zeros added at upper thermal limit (~8C)
dat_upper <- dat_hist[1:(nrow(dat_hist)*0.05),]  #add 5% extra rows as zeros  ***need a smart way to calculate penalty here; even very few data points can have big impact
dat_upper[] <- 0
dat_upper$temp <- 8  #estimated upper thermal limit
dat_upper$abundance <- 0  #all zeros
dat_hist2 <- rbind(dat_hist, dat_upper)
gam4 <- gam(round(abundance) ~ s(temp, k=4), data=dat_hist2, family=poisson) gams
#summary(gam4)
#plot(gam4)

##JS: PLOT responses
par(mfrow=c(3,2), mar=c(3,4,4,2))
ylim2 <- 35
new_dat <- data.frame(temp=seq(0,max(dat_hist$temp),length=100))
new_dat2 <- data.frame(temp=seq(0,7,length=100))
#actual TPC
xx <- seq(0, 7, length=100)
yy <- dnorm(xx, mean=4, sd=1)  #Must match function in SimulatedWorld function
plot(xx, yy, type="l", lty=2, main="Actual TPC", col="red", xlim=c(0,8), ylab="suitability", xlab="Temp")
xlim <- round(100*(max(dat_hist$temp)/7))
lines(xx[1:xlim], yy[1:xlim], lwd=2)
#gam 1
plot(new_dat2$temp, predict(gam2.a, newdata=new_dat2, type="response"), type="l",
     main="Poisson GAM", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
points(dat_hist$temp, dat_hist$abundance, col="grey")
lines(new_dat$temp, predict(gam2.a, newdata=new_dat, type="response"), lwd=2)
#gam 2
plot(new_dat2$temp, predict(gam3.a, newdata=new_dat2, type="response"), type="l",
     main="Poisson GAM, k=4", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
points(dat_hist$temp, dat_hist$abundance, col="grey")
lines(new_dat$temp, predict(gam3.a, newdata=new_dat, type="response"), lwd=2)
#gam 3
plot(new_dat2$temp, predict(gam4, newdata=new_dat2, type="response"), type="l",
     main="Poisson GAM, k=4, upperTL", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
points(dat_hist2$temp, dat_hist2$abundance, col="grey")
lines(new_dat$temp, predict(gam4, newdata=new_dat, type="response"), lwd=2)
#BRT
plot(new_dat2$temp, predict(brt2.a, newdata=new_dat2, type="response", n.trees=brt2.a$gbm.call$best.trees), type="l",
     main="Poisson BRT", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
points(dat_hist$temp, dat_hist$abundance, col="grey")
lines(new_dat$temp, predict(brt2.a, newdata=new_dat, type="response", n.trees=brt2.a$gbm.call$best.trees), lwd=2)
#MLP
dat_norm <- dat
dat_norm$temp_C <- dat_norm$temp
dat_norm$temp <- BBmisc::normalize(dat_norm$temp)
temp_mlp <- dat_norm[order(dat_norm$temp),]
temp_mlp_hist <- temp_mlp[temp_mlp$year <= 2020,]
plot(temp_mlp$temp_C, predict(mlp2.a, temp_mlp), type="l",
     main="Poisson MLP", xlim=c(0,8), col="red", lty=2, ylim=c(0,ylim2), ylab="Abundance", xlab="Temp")
points(dat_hist$temp, dat_hist$abundance, col="grey")
lines(temp_mlp_hist$temp_C, predict(mlp2.a, temp_mlp_hist), lwd=2)
par(mfrow=c(1,1))





#----Make Predictions for the future----

#Run if lognormal response was simulated
if (abund_enviro == "lnorm_low" | abund_enviro == "lnorm_high"){
  #GAM Hindcast (aka Fitted values)
  dat_hist$gam1.p <- predict(gam1.p,dat_hist,type='response')
  dat_hist$gam1.a <- predict(gam1.a,dat_hist,type="response")
  dat_hist$gam1 <- dat_hist$gam1.p*exp(dat_hist$gam1.a)
  dat_hist$brt1.p <- predict(brt1.p,dat_hist,n.trees=brt1.p$gbm.call$best.trees,type='response')
  dat_hist$brt1.a <- predict(brt1.a,dat_hist,n.trees=brt1.a$gbm.call$best.trees,type='response')
  dat_hist$brt1 <- dat_hist$brt1.p*exp(dat_hist$brt1.a)
  dat_hist$mlp1.p <- predict(mlp1.p,dat_hist_norm)
  dat_hist$mlp1.a <- predict(mlp1.a,dat_hist_norm)
  dat_hist$mlp1 <- dat_hist$mlp1.p*exp(dat_hist$mlp1.a)
  
  #GAM Forecast
  dat_fcast$gam1.p <- predict(gam1.p,dat_fcast,type='response')
  dat_fcast$gam1.a <- predict(gam1.a,dat_fcast,type="response")
  dat_fcast$gam1 <- dat_fcast$gam1.p*exp(dat_fcast$gam1.a)
  dat_fcast$brt1.p <- predict(brt1.p,dat_fcast,n.trees=brt1.p$gbm.call$best.trees,type='response')
  dat_fcast$brt1.a <- predict(brt1.a,dat_fcast,n.trees=brt1.a$gbm.call$best.trees,type='response')
  dat_fcast$brt1 <- dat_fcast$brt1.p*exp(dat_fcast$brt1.a)
  dat_fcast$mlp1.p <- predict(mlp1.p,dat_fcast_norm)
  dat_fcast$mlp1.a <- predict(mlp1.a,dat_fcast_norm)
  dat_fcast$mlp1 <- dat_fcast$mlp1.p*exp(dat_fcast$mlp1.a)
  
  #Standard errors Historical
  testCI1.a <- predict(gam1.a, dat_hist, type="response",se.fit=TRUE)
  testCI1.p <- predict(gam1.p, dat_hist, type="response",se.fit=TRUE)
  dat_hist$gam1.high <- exp((testCI1.a$fit + (testCI1.a$se.fit))) * (testCI1.p$fit + (testCI1.p$se.fit))
  dat_hist$gam1.low <- exp((testCI1.a$fit - (testCI1.a$se.fit))) * (testCI1.p$fit - (testCI1.p$se.fit))
  
  #Standard errors Forecast
  testCI1.a <- predict(gam1.a, dat_fcast, type="response",se.fit=TRUE)
  testCI1.p <- predict(gam1.p, dat_fcast, type="response",se.fit=TRUE)
  dat_fcast$gam1.high <- exp((testCI1.a$fit + (testCI1.a$se.fit))) * (testCI1.p$fit + (testCI1.p$se.fit))
  dat_fcast$gam1.low <- exp((testCI1.a$fit - (testCI1.a$se.fit))) * (testCI1.p$fit - (testCI1.p$se.fit))
  
  #Errors from BRT can be generated, I just haven't added the code yet (steph)
  #Need to ask Barb if errors can be generated from neural networks. 

}

#Run if poisson response was simulated
if (abund_enviro == "poisson"){

  #GAM Hindcast (aka Fitted values)
  dat_hist$gam2.a <- predict(gam2.a,dat_hist,type="response")
  dat_hist$gam3.a <- predict(gam2.a,dat_hist,type="response")
  dat_hist$brt2.a <- predict(brt2.a,dat_hist,n.trees=brt2.a$gbm.call$best.trees,type='response')
  dat_hist$mlp2.a <- predict(mlp2.a,dat_hist_norm)
  
  #GAM Forecast
  dat_fcast$gam2.a <- predict(gam2.a,dat_fcast,type="response")
  dat_fcast$gam3.a <- predict(gam3.a,dat_fcast,type="response")
  dat_fcast$brt2.a <- predict(brt2.a,dat_fcast,n.trees=brt2.a$gbm.call$best.trees,type='response')
  dat_fcast$mlp2.a <- predict(mlp2.a,dat_fcast_norm)
}


#----Compare abundance predictons----
#Quick and dirty plots (will convert to ggplot at some point)

#Run if lognormal response was simulated
if (abund_enviro == "lnorm_low" | abund_enviro == "lnorm_high"){
  #Historical patterns
  plot(aggregate(abundance~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance")
  lines(aggregate(gam1~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="blue")
  lines(aggregate(gam1.high~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="light blue")
  lines(aggregate(gam1.low~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="light blue")
  lines(aggregate(brt1~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="red")
  lines(aggregate(mlp1~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="green")
  legend("topright",c("Truth","GAM","BRT","MLP"),lty=1,col=c("black","blue","red", "green"),bty="n")
  
  #Future patterns
  plot(aggregate(abundance~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",ylim=c(0,1700))
  lines(aggregate(gam1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='blue')
  lines(aggregate(gam1.high~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='light blue')
  lines(aggregate(gam1.low~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='light blue')
  lines(aggregate(brt1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='red')
  lines(aggregate(mlp1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='green')
  legend("topright",c("Truth","GAM","BRT", "MLP"),lty=1,col=c("black","blue","red", "green"),bty="n")
}

#Run if poisson response was simulated
if (abund_enviro == "poisson"){
  #Historical patterns
  plot(aggregate(round(abundance)~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance")
  lines(aggregate(gam2.a~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="blue")
  lines(aggregate(gam3.a~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="light blue")
  lines(aggregate(brt2.a~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="red")
  lines(aggregate(mlp2.a~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="green")
  legend("topright",c("Truth","GAM_normal","GAM_restricted","BRT","MLP"),lty=1,col=c("black","blue","light blue","red", "green"),bty="n")
  
  #Future patterns
  plot(aggregate(round(abundance)~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",ylim=c(0,1700))
  lines(aggregate(gam2.a~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='blue')
  lines(aggregate(gam3.a~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='light blue')
  lines(aggregate(brt2.a~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='red')
  lines(aggregate(mlp2.a~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='green')
  legend("topright",c("Truth","GAM_normal","GAM_restricted","BRT","MLP"),lty=1,col=c("black","blue","light blue","red", "green"),bty="n")
}


#-----Calculate and plot centre of gravity-----

#Run if lognormal response was simulated
if (abund_enviro == "lnorm_low" | abund_enviro == "lnorm_high"){
  #Historical COG
  cog_hist_lat <- as.data.frame(matrix(NA,nrow=nrow(dat_hist),ncol=8))
  colnames(cog_hist_lat) <- c("year","truth","gam1","gam1.p","brt1","brt1.p","mlp1","mlp1.p")
  counter=1
  for (y in unique(dat_hist$year)){
    cog_hist_lat[counter,1] <- y
    cog_hist_lat[counter,2] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$abundance[dat_hist$year==y])
    cog_hist_lat[counter,3] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam1[dat_hist$year==y])
    cog_hist_lat[counter,4] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam1.p[dat_hist$year==y])
    cog_hist_lat[counter,5] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt1[dat_hist$year==y])
    cog_hist_lat[counter,6] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt1.p[dat_hist$year==y])
    cog_hist_lat[counter,7] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$mlp1[dat_hist$year==y])
    cog_hist_lat[counter,8] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$mlp1.p[dat_hist$year==y])
    counter = counter + 1
  }
  head(cog_hist_lat)
  plot(cog_hist_lat$year,cog_hist_lat$truth, type='b')
  lines(cog_hist_lat$year,cog_hist_lat$gam1, type='b', col="blue")
  lines(cog_hist_lat$year,cog_hist_lat$brt1, type='b', col="red")
  lines(cog_hist_lat$year,cog_hist_lat$mlp1, type='b', col="green")
  
  #Future COG
  cog_fcast_lat <- as.data.frame(matrix(NA,nrow=nrow(dat_fcast),ncol=8))
  colnames(cog_fcast_lat) <- c("year","truth","gam1","gam1.p","brt1","brt1.p","mlp1","mlp1.p")
  counter=1
  for (y in unique(dat_fcast$year)){
    cog_fcast_lat[counter,1] <- y
    cog_fcast_lat[counter,2] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$abundance[dat_fcast$year==y])
    cog_fcast_lat[counter,3] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam1[dat_fcast$year==y])
    cog_fcast_lat[counter,4] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam1.p[dat_fcast$year==y])
    cog_fcast_lat[counter,5] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt1[dat_fcast$year==y])
    cog_fcast_lat[counter,6] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt1.p[dat_fcast$year==y])
    cog_fcast_lat[counter,7] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$mlp1[dat_fcast$year==y])
    cog_fcast_lat[counter,8] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$mlp1.p[dat_fcast$year==y])
    counter = counter + 1
  }
  head(cog_fcast_lat)
  plot(cog_fcast_lat$year,cog_fcast_lat$truth, type='b')
  lines(cog_fcast_lat$year,cog_fcast_lat$gam1, type='b', col="blue")
  lines(cog_fcast_lat$year,cog_fcast_lat$brt1, type='b', col="red")
  lines(cog_fcast_lat$year,cog_fcast_lat$mlp1, type='b', col="green")
}


#Run if poisson response was simulated
if (abund_enviro == "poisson"){
  #Historical COG
  cog_hist_lat <- as.data.frame(matrix(NA,nrow=20,ncol=6))
  colnames(cog_hist_lat) <- c("year","truth","gam2.a","gam3.a","brt2.a","mlp2.a")
  counter=1
  for (y in 2001:2020){
    cog_hist_lat[counter,1] <- y
    cog_hist_lat[counter,2] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$abundance[dat_hist$year==y])
    cog_hist_lat[counter,3] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam2.a[dat_hist$year==y])
    cog_hist_lat[counter,4] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam3.a[dat_hist$year==y])
    cog_hist_lat[counter,5] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt2.a[dat_hist$year==y])
    cog_hist_lat[counter,6] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$mlp2.a[dat_hist$year==y])
    counter = counter + 1
  }
  head(cog_hist_lat)
  plot(cog_hist_lat$year,cog_hist_lat$truth, type='b')
  lines(cog_hist_lat$year,cog_hist_lat$gam2.a, type='b', col="blue")
  lines(cog_hist_lat$year,cog_hist_lat$brt2.a, type='b', col="red")
  lines(cog_hist_lat$year,cog_hist_lat$mlp2.a, type='b', col="green")
  
  #Future COG
  cog_fcast_lat <- as.data.frame(matrix(NA,nrow=80,ncol=6))
  colnames(cog_fcast_lat) <- c("year","truth","gam2.a","gam3.a","brt2.a","mlp2.a")
  counter=1
  for (y in 2021:2100){
    cog_fcast_lat[counter,1] <- y
    cog_fcast_lat[counter,2] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$abundance[dat_fcast$year==y])
    cog_fcast_lat[counter,3] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam2.a[dat_fcast$year==y])
    cog_fcast_lat[counter,4] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam3.a[dat_fcast$year==y])
    cog_fcast_lat[counter,5] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt2.a[dat_fcast$year==y])
    cog_fcast_lat[counter,6] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$mlp2.a[dat_fcast$year==y])
    counter = counter + 1
  }
  head(cog_fcast_lat)
  plot(cog_fcast_lat$year,cog_fcast_lat$truth, type='b')
  lines(cog_fcast_lat$year,cog_fcast_lat$gam2.a, type='b', col="blue")
  lines(cog_fcast_lat$year,cog_fcast_lat$brt3.a, type='b', col="red")
  lines(cog_fcast_lat$year,cog_fcast_lat$mlp2.a, type='b', col="green")
}
#-----Plot Surface Predictions-----
#Note these are point predictions for the ROMS data (not a prediction on the whole surface)
#Future
Y = 2021
#Truth
ggplot(dat_fcast[dat_fcast$year==Y,],aes(Lon,Lat))+
  geom_tile(aes(fill=abundance)) +
  theme_classic() +
  ggtitle("Truth")+
  labs(y="Latitude") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust=0.5),
        # axis.text.x=element_blank(),
        # axis.ticks=element_blank(),
        # axis.title.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_viridis()

#Gam
ggplot(dat_fcast[dat_fcast$year==Y,],aes(Lon,Lat))+
  geom_tile(aes(fill=gam1)) +
  theme_classic() +
  ggtitle("GAM")+
  labs(y="Latitude") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust=0.5),
        # axis.text.x=element_blank(),
        # axis.ticks=element_blank(),
        # axis.title.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_viridis()

#BRT
ggplot(dat_fcast[dat_fcast$year==Y,],aes(Lon,Lat))+
  geom_tile(aes(fill=brt1)) +
  theme_classic() +
  ggtitle("BRT")+
  labs(y="Latitude") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust=0.5),
        # axis.text.x=element_blank(),
        # axis.ticks=element_blank(),
        # axis.title.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_viridis()

#MLP
ggplot(dat_fcast[dat_fcast$year==Y,],aes(Lon,Lat))+
  geom_tile(aes(fill=mlp1)) +
  theme_classic() +
  ggtitle("MLP")+
  labs(y="Latitude") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous( expand = c(0, 0)) +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust=0.5),
        # axis.text.x=element_blank(),
        # axis.ticks=element_blank(),
        # axis.title.x=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_fill_viridis()

