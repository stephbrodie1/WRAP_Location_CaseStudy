#############################################################
######    Location^3 Species Distribution Modelling    ######
######             'Anchovy' Archetype                ######
######        May 2020; built with R v3.6.3          ######
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
library(sdmTMB)
library(sp)
library(dplyr)

#################################
##       Data Simulation       ##
#################################

# load ROMS simulation function
source('~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/SimulatedWorld_ROMS_Anchovy_Function.R')
# * JS: with random cell selection there were too many zeros, so I changed 'convertToPA' options to alpha=-0.1, beta=0.4, among other changes

#Three ESM directories
# had <- "~/Dropbox/WRAP Location^3/Rasters_2d_Spring/had/"
# ipsl <- "~/Dropbox/WRAP Location^3/Rasters_2d_Spring/ipsl/"
gfdl <- "~/Dropbox/WRAP Location^3/Rasters_2d_Spring/gfdl/"

# create abundance data - choose ESM directory
dir_name <- "gfdl"     #of either "had", "ipsl", "gfdl"
dat <- SimulateWorld_ROMS_Anchovy(dir = gfdl)  #simulate data (takes a few mins)

# name changes
names(dat)[names(dat) == 'sst'] <- 'temp'
names(dat)[names(dat) == 'chla_surface'] <- 'chla'

# add scaled and 'range normalised' variables for GLMs and MLPs
dat$temp_s = scale(dat$temp)
dat$chla_s = scale(dat$chla)
dat$z_s = scale(dat$z)

dat$temp_n <- BBmisc::normalize(dat$temp, method = "range", range = c(0, 1))  #covariates must be normalised for MLP
dat$chla_n <- BBmisc::normalize(dat$chla, method = "range", range = c(0, 1))
dat$z_n <- BBmisc::normalize(dat$z, method = "range", range = c(0, 1))
dat$lat_n <- BBmisc::normalize(dat$lat, method = "range", range = c(0, 1))
dat$lon_n <- BBmisc::normalize(dat$lon, method = "range", range = c(0, 1))
dat$year_n <- BBmisc::normalize(dat$year, method = "range", range = c(0, 1))

# remove less certain ROMS years
dat <- dat[dat$year >= 1985,]

# transform and add factor
dat$log_abundance <- log(dat$abundance)
dat$fYear <- as.factor(dat$year)

# create data subset for model fitting and evaluation (500 locations per year [*random, but same every year])
dat_sub <- dat[dat$sample == 1,]

# # save
saveRDS(dat, paste0("~/PROJECTS/WRAP Location/Anchovy EMs/",dir_name,"/dat_anchovy_all.rds"))
saveRDS(dat_sub, paste0("~/PROJECTS/WRAP Location/Anchovy EMs/",dir_name,"/dat_anchovy_sub.rds"))


#################################
##        Model Fitting        ##
#################################

rm(list = ls())  #remove all objects and re-load necessary ones; consider restarting R too

# load data subset for model fitting
dir_name <- "ipsl"     #of either "had", "ipsl", "gfdl"
dat <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/",dir_name,"/dat_anchovy_sub.rds"))
# year_fcast <- 2010
year_fcast <- 2040
dat_hist <- dat[dat$year <= year_fcast,]  #subset for model fitting
dat_fcast <- dat[dat$year > year_fcast,]  #subset for forecast evaluation


# --- GAMS ---
# * fitting ECor will take 1-2 hours
type <- "delta"  #delta or tweedie
covs <- c("E","S","ES","EST", "ECor")  #covariate combinations, options: c("E","S","ES","EST","ECor")
env_formula <- "s(temp) + s(chla) + s(z)"  #full: "s(temp) + s(chla) + s(mld)"
# env_formula <- "s(temp)"
t1 <- Sys.time()
source("~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/Fitting_GAMs.R")
t2 <- Sys.time()
t2-t1

# --- GLMs ---
# * fitting Sr and ESr will take 2-6 hours, and if data are changed may not fit w/o changing model structure
# * warnings are ok, but if there are model failures... seek help
# * all GLMs assume at least one quadratic (i.e. specify no intercept and have quadratic_roots=T)

type <- "delta"  #delta or tweedie
covs <- c("E", "ESt", "Sr", "ESr")   #options: c("E", "Sr", "ESt", "ESr") 
env_formula <- "temp_s + I(temp_s^2) + z_s + I(z_s^2) + chla_s + I(chla_s^2)"  #for tweedie models and presence part of delta; the '_s' indicates scaled variables
env_formula_deltaN <- "temp_s + I(temp_s^2) + z_s + I(z_s^2) + chla_s"  #no quadratic for chl in abundance part of delta
# env_formula <- "temp_s + I(temp_s^2)"  #temp only model
# env_formula_deltaN <- "temp_s + I(temp_s^2)"  #temp only model

t1 <- Sys.time()
source("~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/Fitting_GLMs.R")
t2 <- Sys.time()
t2-t1

# --- BRTs & MLPs ---
# * check MLP response curves - if some are flat (constant), run them again... then seek help

type <- "delta"  #only delta at this stage
covs <- c("E","ES","EST")  #options: c("E","S","ES","EST") #'S' is not really needed, it's run in GAMs as a simple NULL model
env_formula_BRT <- c("temp", "chla", "z")  #full: c("temp", "chla", "mld")
env_formula_MLP <- "temp_n + chla_n + z_n"  #full: "temp_n + chla_n + mld_n"; the '_n' indicates range normalised variables
# env_formula_BRT <- c("temp")  #full: c("temp", "chla", "mld")
# env_formula_MLP <- "temp_n"

source("~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/Fitting_BRTs.R")
source("~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/Fitting_MLPs.R")


# --- Save Results and Models ---
# * take care not to delete some models by overwriting with a model subset (change saved names if doing testing)
setwd(paste0('~/PROJECTS/WRAP Location/Anchovy EMs 2040train/',dir_name,'/'))
saveRDS(dat_hist, "dat_hist_results_full.rds")  # * 'full' [full models] or 'temp' [temp-only models]
saveRDS(dat_fcast, "dat_fcast_results_full.rds")

all_mods <- c("gam_E","gam_S","gam_ES","gam_EST","gam_ECor",
              "brt_E","brt_S","brt_ES","brt_EST",
              "mlp_E","mlp_S","mlp_ES","mlp_EST",
              "glm_E","glm_Sr","glm_ESt","glm_ESr")

save(list = ls(pattern = paste0(all_mods, collapse="|")),
     file = "saved_models_full.RData")  # * 'full' or 'temp'; this will save all models of these names, even when outdated, if they're in the environment


###################################
##       Model Performance       ##
###################################

rm(list = ls())  #remove all objects and re-load necessary ones

# --- load functions, data, and model objects ---
dir_name <- "gfdl"
setwd('~/PROJECTS/WRAP Location/WRAP_Location_CaseStudy/')
source("Performance_rmse_fun.R")
source("Performance_cog_fun.R")
source("Performance_eao_fun.R")

# * load these, for full models
dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/",dir_name,"/dat_hist_results_full.rds"))
dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/",dir_name,"/dat_fcast_results_full.rds"))
load(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/",dir_name,"/saved_models_full.RData"))  #full models
# * note that the glm_ESr full model has a different structure to other glms (abundance part of delta has no quadratics - this was done to speed up fitting)

# * or load these, for temp-only models
# dat_hist <- readRDS("dat_hist_results_temp.rds")  #temp-only models
# dat_fcast <- readRDS("dat_fcast_results_temp.rds")
# load("saved_models_temp.RData")  #temp-only models


# --- RMSE ---
ylim = c(140000,260000)
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
#rmse_results$rmse_abund_all_year  


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

eao_results$eao_hist_y
eao_results$eao_fcast_y
# * Most models are positively biased - may be due to lack of zeros in modelled data?


# --- Maps of Spatial Distribution ---
dir_name <- "gfdl"
dat_all <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/",dir_name,"/dat_anchovy_all.rds"))
source("Performance_plotMaps.R")
performance_plotMaps_function(om ="cps") #om options are "hms","cps","grf". Needed for GLMpredict function within maps







