#############################################################
######    Location^3 Species Distribution Modelling    ######
######             'Albacore' Archetype                ######
######        April 2020; built with R v3.6.3          ######
#############################################################

# --- Notes ---

# * All delta models are fit to normally distributed abundance
# * Performance is measured using a data subset (500 obs per year), but testing shows this is enough to accurately estimate RMSE etc
# * Take care comparing fitted responses to the specified suitability functions, bc they may not represent the fundamental niche when there are multiple enviro drivers
# * It appears that the rate of degregation of the forecast is less about which model, and more about the (non-linear?) rate of environmental change, and
#   ... whether the fitting period covers a broad range of enviro conditions likely to be encountered by the species
#   ... In the albacore case, it seems the N-ward shift in COG during 1985-2020 continued at same rate into future, making forecasting quite successful
# * The observations are still quite 'knife-edge' in terms of presence, and we could explore other options for the 'convertToPA' function in VirtualSpecies
# * Adding lat, lon or year often has little impact for albacore, bc there is not much spatially- or temporally-structured residual information
#   ... which encouarges comparing models with diff numbers of covariates (e.g. full model vs temp only)

# --- To Do ---

# * Consider comparison of:
#   ... 1) rate of change of warming (and change rate from historical period)
#   ... 2) number of covariates, and adding 'unimportant' ones that vary in collinearity with important ones
#   ... 3) duration of fitting period, or remove warmest temps from sampled observations


#################################
##       Load Libraries        ##
#################################

library(mgcv)
library(gbm)
library(dismo)
library(BBmisc)
library(neuralnet)
# devtools::install_github("pbs-assess/sdmTMB")
# library(sdmTMB)
library(sp)
library(dplyr)
library(foreach)
library(doParallel, quietly = TRUE)


#################################
##       Data Simulation       ##
#################################
# load ROMS simulation function
source('~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/SimulatedWorld_ROMS_Albacore_Function.R')

#Three ESM directories
# had <- "~/Dropbox/WRAP Location^3/Rasters_2d_Spring/had/"
# ipsl <- "~/Dropbox/WRAP Location^3/Rasters_2d_Spring/ipsl/"
gfdl <- "~/Dropbox/WRAP Location^3/Rasters_2d_Spring/gfdl/"

# create abundance data - choose ESM directory
dir_name <- "gfdl"     #of either "had", "ipsl", "gfdl"
dat <- SimulateWorld_ROMS_Albacore(dir = gfdl)  

# name changes
names(dat)[names(dat) == 'sst'] <- 'temp'
names(dat)[names(dat) == 'chla_surface'] <- 'chla'

# add scaled and 'range normalised' variables for GLMs and MLPs
dat$temp_s = scale(dat$temp)
dat$chla_s = scale(dat$chla)
dat$mld_s = scale(dat$mld)
dat$temp_n <- BBmisc::normalize(dat$temp, method = "range", range = c(0, 1))  #covariates must be normalised for MLP
dat$chla_n <- BBmisc::normalize(dat$chla, method = "range", range = c(0, 1))
dat$mld_n <- BBmisc::normalize(dat$mld, method = "range", range = c(0, 1))
dat$lat_n <- BBmisc::normalize(dat$lat, method = "range", range = c(0, 1))
dat$lon_n <- BBmisc::normalize(dat$lon, method = "range", range = c(0, 1))
dat$year_n <- BBmisc::normalize(dat$year, method = "range", range = c(0, 1))

# transform and add factor
dat$log_abundance <- log(dat$abundance) #Must use log abundance to get some EMs to run
dat$fYear <- as.factor(dat$year)

# remove less certain ROMS years
dat <- dat[dat$year >= 1985,]

# create data subset for model fitting and evaluation (500 locations per year [*random, but same every year])
dat_sub <- dat[dat$sample == 1,]

# # save
saveRDS(dat, paste0("~/PROJECTS/WRAP Location/Albacore EMs/",dir_name,"/dat_albacore_all.rds"))
saveRDS(dat_sub, paste0("~/PROJECTS/WRAP Location/Albacore EMs/",dir_name,"/dat_albacore_sub.rds"))


#################################
##        Model Fitting        ##
#################################

rm(list = ls())  #remove all objects and re-load necessary ones; consider restarting R too

# load data subset for model fitting
dir_name <- "ipsl"     #of either "had", "ipsl", "gfdl"
dat <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/",dir_name,"/dat_albacore_sub.rds"))
# year_fcast <- 2010
year_fcast <- 2040
dat_hist <- dat[dat$year <= year_fcast,]  #subset for model fitting
dat_fcast <- dat[dat$year > year_fcast,]  #subset for forecast evaluation

# #Test extrapolation:
# dat_hist$sample <- ifelse(dat_hist$lat<=43 & dat_hist$lat>=35 & dat_hist$lat<=year_fcast,dat_hist$sample==1,dat_hist$sample==0)
# dat_hist <- dat_hist[dat_hist$sample==1,]
# dat <- rbind(dat_hist,dat_fcast)


# --- GAMS ---
# * fitting ECor will take 1-2 hours
type <- "delta"  #delta or tweedie
covs <- c("E","S","ES","EST","ECor")  #covariate combinations, options: c("E","S","ES","EST","ECor")
env_formula <- "s(temp) + s(chla) + s(mld)"  #full: "s(temp) + s(chla) + s(mld)"
# env_formula <- "s(temp)"

source("~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/Fitting_GAMs.R")

# --- GLMs ---
# * fitting Sr and ESr will take 2-6 hours, and if data are changed may not fit w/o changing model structure
# * warnings are ok, but if there are model failures... seek help
# * all GLMs assume at least one quadratic (i.e. specify no intercept and have quadratic_roots=T)

type <- "delta"  #delta or tweedie
covs <- c("E", "ESt", "Sr", "ESr")   #options: c("E", "Sr", "ESt", "ESr") 
env_formula <- "temp_s + I(temp_s^2) + mld_s + I(mld_s^2) + chla_s + I(chla_s^2)"  #for tweedie models and presence part of delta; the '_s' indicates scaled variables
env_formula_deltaN <- "temp_s + I(temp_s^2) + mld_s + I(mld_s^2) + chla_s"  #no quadratic for chl in abundance part of delta
# env_formula <- "temp_s + I(temp_s^2)"  #temp only model
# env_formula_deltaN <- "temp_s + I(temp_s^2)"  #temp only model

t1 <- Sys.time()
source("~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/Fitting_GLMs.R")
t2 <- Sys.time()
t2 - t1

# --- BRTs & MLPs ---
# * check MLP response curves - if some are flat (constant), run them again... then seek help

type <- "delta"  #only delta at this stage
covs <- c("E","ES","EST")  #options: c("E","S","ES","EST")
env_formula_BRT <- c("temp", "chla", "mld")  #full: c("temp", "chla", "mld")
env_formula_MLP <- "temp_n + chla_n + mld_n"  #full: "temp_n + chla_n + mld_n"; the '_n' indicates range normalised variables
# env_formula_BRT <- c("temp")
# env_formula_MLP <- "temp_n"  #full: "temp_n + chla_n + mld_n"; the '_n' indicates range normalised variables


t1 <- Sys.time()
source("~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/Fitting_BRTs.R")
source("~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/Fitting_MLPs.R")
t2 <- Sys.time()
t2 - t1

# --- Save Results and Models ---
# * take care not to delete some models by overwriting with a model subset (change saved names if doing testing)
setwd(paste0('~/PROJECTS/WRAP Location/Albacore EMs 2040train/',dir_name,'/'))
saveRDS(dat_hist, "dat_hist_results_full.rds")  # * 'full' [full models] or 'temp' [temp-only models]
saveRDS(dat_fcast, "dat_fcast_results_full.rds")

all_mods <- c("gam_E","gam_S","gam_ES","gam_EST","gam_ECor",
              "brt_E","brt_S","brt_ES","brt_EST",
              "mlp_E","mlp_S","mlp_ES","mlp_EST",
              "glm_E","glm_Sr","glm_ESt","glm_ESr")

save(list = ls(pattern = paste0(all_mods, collapse="|")),
     file = "saved_models_full.RData")  # * 'full' or 'temp'; this will save all models of these names, even when outdated, if they're in the environment


#hist_save <- dat_hist
#fcast_save <- dat_fcast

dat_fcast <- fcast_save
dat_fcast <- dat_fcast[,c(1:32,34:38)]

dat_hist <- dat_hist[,c(1:22,27:31,23:26,32:37)]



###################################
##       Model Performance       ##
###################################

rm(list = ls())  #remove all objects and re-load necessary ones

# --- load functions, data, and model objects ---
dir_name <- "had"
setwd('~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/')
source("Performance_rmse_fun.R")
source("Performance_cog_fun.R")
source("Performance_eao_fun.R")

# * load these, for full models
dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs TempOnly/",dir_name,"/dat_hist_results_full.rds"))
dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs TempOnly/",dir_name,"/dat_fcast_results_full.rds"))
load(paste0("~/PROJECTS/WRAP Location/Albacore EMs TempOnly/",dir_name,"/saved_models_full.RData"))  #full models
# * note that the glm_ESr full model has a different structure to other glms (abundance part of delta has no quadratics - this was done to speed up fitting)

# * or load these, for temp-only models
# dat_hist <- readRDS("dat_hist_results_temp.rds")  #temp-only models
# dat_fcast <- readRDS("dat_fcast_results_temp.rds")
# load("saved_models_temp.RData")  #temp-only models



# --- RMSE ---
ylim=c(1000,8000)
rmse_results <- sdm_rmse(dat_hist = dat_hist,
                         dat_fcast = dat_fcast)
# * space-only models are flat lines, bc the subset of locations sampled each year are identical

# RMSE of all observations (500 per year)
rmse_results$rmse_abund_all[order(rmse_results$rmse_abund_all$rmse_hist),]
#^ Hist: ordered best -> worst for historical
rmse_results$rmse_abund_all[order(rmse_results$rmse_abund_all$rmse_fcast),]
#^ Fcast: ordered best -> worst for forecast

# RMSE of annual total abundance (sum of 500 cells per year)
rmse_results$rmse_abund_annual[order(rmse_results$rmse_abund_annual$rmse_hist),]
#^ Hist: ordered best -> worst for historical
rmse_results$rmse_abund_annual[order(rmse_results$rmse_abund_annual$rmse_fcast),]
#^ Fcast: ordered best -> worst for forecast

# Annual RMSE for all observations (plotted)
# rmse_results$rmse_abund_all_year


# --- COG Latitude ---

cog_results <- sdm_cog(dat_hist = dat_hist,
                       dat_fcast = dat_fcast)

cog_results$rmse_cog[order(cog_results$rmse_cog$rmse_hist),]
#Hist: glm_Sr best, then glm_ESr, brt_EST, brt_ES, gam_EST ...
cog_results$rmse_cog[order(cog_results$rmse_cog$rmse_fcast),]
#Fcast: lm_ESr best, then glm_ESt, brt_E, glm_E, brt_ES ...

# Difference in cog over time (plotted)
#cog_results$cog_lat_fcast_y  #not often a trend for enviro+ SDMs


# --- Effective Area Occupied ---

eao_results <- sdm_eao(dat_hist = dat_hist,
                       dat_fcast = dat_fcast)

eao_results$eao_hist_y[]
eao_results$eao_fcast_y
# * Most models are positively biased - may be due to lack of zeros in modelled data?


# --- Maps of Spatial Distribution ---
# load full data
dat_all <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/",dir_name,"/dat_albacore_all.rds"))
source("Performance_plotMaps.R")
performance_plotMaps_function(om ="hms") #om options are "hms","cps","grf". Needed for GLMpredict function within maps



#-----summary plots of OM----

# # #Make some quick plots to explore the data
# # #All Years
par(mfrow=c(2,3))
dat_all <- dat
dat_sub <- dat[dat$sample==1,]
plot(aggregate(suitability~year,dat,FUN="mean"),type="l", lwd=2, main="Mean Suitability",col="dark grey", ylab="")
lines(aggregate(suitability~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(pres~year,dat,FUN="mean"),type="l", lwd=2, main="Prevalence",col="dark grey", ylab="")
lines(aggregate(pres~year,dat[dat$year<=2010,],FUN="mean"),col="blue")

plot(aggregate(abundance~year,dat_sub,FUN="sum"),type="l",  lwd=2,main="Total Biomass", col="dark grey", ylab="")
lines(aggregate(abundance~year,dat_sub[dat_sub$year<=2010,],FUN="sum"),col="blue")
plot(aggregate(abundance~year,dat_all,FUN="sum"),type="l",  lwd=2,main="Total Biomass", col="dark grey", ylab="")
lines(aggregate(abundance~year,dat_all[dat_all$year<=2010,],FUN="sum"),col="blue")


plot(aggregate(temp~year,dat,FUN="mean"),type="l",main="Mean Temperature", col="dark grey", ylab="")
lines(aggregate(temp~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(mld~year,dat,FUN="mean"),type="l",main="Mean MLD",col="dark grey")
lines(aggregate(mld~year,dat[dat$year<=2010,],FUN="mean"),col="blue")
plot(aggregate(zoo_200~year,dat,FUN="mean"),type="l",main="Mean Zooplankton",col="dark grey", ylab="")
lines(aggregate(zoo_200~year,dat[dat$year<=2010,],FUN="mean"),col="blue")




