#Code for Location WRAP Workshop
#Compare models using simulated data

#----Directories----
#Set your working directory!
setwd('~/PROJECTS/WRAP Location/')
dir.create('Sim1') #data and model outputs will be saved locally here
Sim1 <- ('Sim1/')

#----Load Library & Function----
library(dismo)
library(gbm)
library(mgcv)
library(ggplot2)
library(viridis)
source('WRAP_Location_CaseStudy/SimulatedWorld_Function.R') #load simulation function

#-----Simulate data----
dat <- SimulateWorld() #takes a few minutes
colnames(dat)[1:2] <- c("Lon","Lat")
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
plot(aggregate(temp~year,dat,FUN="min"),type="l",ylab="Temperature",ylim=c(-2,15), col="dark grey")
lines(aggregate(temp~year,dat,FUN="max"),type="l",col="dark grey")
lines(aggregate(temp~year,dat,FUN="mean"),type="l")

#----Build GAM Models----
dat_hist$log_abundance <- log(dat_hist$abundance)
gam1.p <- gam(pres ~ s(temp,bs='gp') , data=dat_hist, family=binomial)
gam1.a <- gam(log_abundance ~ s(temp,bs='gp')  , data=dat_hist[dat_hist$abundance>0,], family=gaussian)
saveRDS(gam1.p, paste0(Sim1,'GAM_Sim1_binom.rds'))
saveRDS(gam1.a, paste0(Sim1,'GAM_Sim1_lognorm.rds'))
# gam1.p <- readRDS( paste0(Sim1,'GAM_Sim1_binom.rds')) #read in model object if required
# gam1.a <- readRDS( paste0(Sim1,'GAM_Sim1_lognorm.rds'))

summary(gam1.p)
summary(gam1.a)
plot(gam1.p)
plot(gam1.a)

#----Boosted Regression Tree----
#Make sure >1000 trees fitted
brt1.a <- gbm.step(data=dat_hist[dat_hist$abundance>0,], gbm.x = c("temp"),gbm.y = 'log_abundance',family = "gaussian",tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6)
brt1.p <- gbm.step(data=dat_hist, gbm.x = c("temp"),gbm.y = 'pres',family = "bernoulli",tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.6)
saveRDS(brt1.a,paste0(Sim1,'BRT_Sim1_lognorm.rds'))
saveRDS(brt1.p,paste0(Sim1,'BRT_Sim1_binom.rds'))
# brt1.a <- readRDS(paste0(Sim1,'BRT_Sim1_lognorm.rds'))#read in model object if required
# brt1.p <- readRDS(paste0(Sim1,'BRT_Sim1_binom.rds'))

dev_eval=function(model_object){
  null <- model_object$self.statistics$mean.null
  res <- model_object$self.statistics$mean.resid
  dev=((null - res)/null)*100 
  return(dev)
}
dev_eval(brt1.p)
dev_eval(brt1.a)

plot(brt1.p)
plot(brt1.a)

#----Make Predictions for the future----
#GAM Hindcast (aka Fitted values)
dat_hist$gam1.p <- predict(gam1.p,dat_hist,type='response')
dat_hist$gam1.a <- predict(gam1.a,dat_hist,type="response")
dat_hist$gam1 <- dat_hist$gam1.p*exp(dat_hist$gam1.a)
dat_hist$brt1.p <- predict(brt1.p,dat_hist,n.trees=brt1.p$gbm.call$best.trees,type='response')
dat_hist$brt1.a <- predict(brt1.a,dat_hist,n.trees=brt1.a$gbm.call$best.trees,type='response')
dat_hist$brt1 <- dat_hist$brt1.p*exp(dat_hist$brt1.a)

#GAM Forecast
dat_fcast$gam1.p <- predict(gam1.p,dat_fcast,type='response')
dat_fcast$gam1.a <- predict(gam1.a,dat_fcast,type="response")
dat_fcast$gam1 <- dat_fcast$gam1.p*exp(dat_fcast$gam1.a)
dat_fcast$brt1.p <- predict(brt1.p,dat_fcast,n.trees=brt1.p$gbm.call$best.trees,type='response')
dat_fcast$brt1.a <- predict(brt1.a,dat_fcast,n.trees=brt1.a$gbm.call$best.trees,type='response')
dat_fcast$brt1 <- dat_fcast$brt1.p*exp(dat_fcast$brt1.a)

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

#----Compare abundance predictons----
#Quick and dirty plots (will convert to ggplot at some point)
#Historical patterns
plot(aggregate(abundance~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance")
lines(aggregate(gam1~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="blue")
lines(aggregate(gam1.high~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="light blue")
lines(aggregate(gam1.low~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="light blue")
lines(aggregate(brt1~year,dat_hist,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col="red")
legend("topright",c("Truth","GAM","BRT"),lty=1,col=c("black","blue","red"),bty="n")

#Future patterns
plot(aggregate(abundance~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",ylim=c(0,1700))
lines(aggregate(gam1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='blue')
lines(aggregate(gam1.high~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='light blue')
lines(aggregate(gam1.low~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='light blue')
lines(aggregate(brt1~year,dat_fcast,FUN="sum"),type="l",  lwd=2,ylab="Abundance",col='red')
legend("topright",c("Truth","GAM","BRT"),lty=1,col=c("black","blue","red"),bty="n")

#Historical COG
cog_hist_lat <- as.data.frame(matrix(NA,nrow=20,ncol=6))
colnames(cog_hist_lat) <- c("year","truth","gam1","gam1.p","brt1","brt1.p")
counter=1
for (y in 2001:2020){
  cog_hist_lat[counter,1] <- y
  cog_hist_lat[counter,2] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$abundance[dat_hist$year==y])
  cog_hist_lat[counter,3] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam1[dat_hist$year==y])
  cog_hist_lat[counter,4] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$gam1.p[dat_hist$year==y])
  cog_hist_lat[counter,5] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt1[dat_hist$year==y])
  cog_hist_lat[counter,6] <- weighted.mean(dat_hist$Lat[dat_hist$year==y],w=dat_hist$brt1.p[dat_hist$year==y])
  counter = counter + 1
}
head(cog_hist_lat)
plot(cog_hist_lat$year,cog_hist_lat$truth, type='b')
lines(cog_hist_lat$year,cog_hist_lat$gam1, type='b', col="blue")
lines(cog_hist_lat$year,cog_hist_lat$brt1, type='b', col="red")

#Future COG
cog_fcast_lat <- as.data.frame(matrix(NA,nrow=80,ncol=6))
colnames(cog_fcast_lat) <- c("year","truth","gam1","gam1.p","brt1","brt1.p")
counter=1
for (y in 2021:2100){
  cog_fcast_lat[counter,1] <- y
  cog_fcast_lat[counter,2] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$abundance[dat_fcast$year==y])
  cog_fcast_lat[counter,3] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam1[dat_fcast$year==y])
  cog_fcast_lat[counter,4] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$gam1.p[dat_fcast$year==y])
  cog_fcast_lat[counter,5] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt1[dat_fcast$year==y])
  cog_fcast_lat[counter,6] <- weighted.mean(dat_fcast$Lat[dat_fcast$year==y],w=dat_fcast$brt1.p[dat_fcast$year==y])
  counter = counter + 1
}
head(cog_fcast_lat)
plot(cog_fcast_lat$year,cog_fcast_lat$truth, type='b')
lines(cog_fcast_lat$year,cog_fcast_lat$gam1, type='b', col="blue")
lines(cog_fcast_lat$year,cog_fcast_lat$brt1, type='b', col="red")

#-----Plot Surface Predictions-----
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


