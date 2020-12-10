#This code is the operating model for a Groundfish archetype, represented by sablefish
#It uses average spring conditions from downscales ROMS projections using Hadley 
#IMPORTANT: Download average spring ROMS data from here: https://www.dropbox.com/sh/pj1vjz4h1n27hnb/AAD9sySDudMKmq3l4yLPvjQra?dl=0

#Code contributions from Kelly Andrews and Steph

dir="~/Dropbox/WRAP Location^3/Rasters_2d_Spring/gfdl/"


SimulateWorld_ROMS_Groundfish <- function(dir){
  #dir is the local directory that points to where ROMS data is stored 
  
  #required libraries
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  library(sf)
  library(dplyr)
  
  #----Create output file----
  #Assuming 500 'samples' are taken each year, from 1980-2100
  #This will be the information passed to the estimation model
  gridcells <- 4012 
  output <- as.data.frame(matrix(NA, nrow=(gridcells*121),ncol=10))
  colnames(output) <- c("lon","lat","year","pres","suitability","btemp","O2", "z",'stemp', "sampled")
  
  #----Identify raster files for adults----
  #These are the average spring conditions from the downscaled Hadley earth system model
  files_stemp <- list.files(paste0(dir,'sst_monthly'), full.names = TRUE, pattern=".grd") 
  files_btemp <- list.files(paste0(dir,'temp_bottom'), full.names = TRUE, pattern=".grd") 
  files_O2 <- list.files(paste0(dir,'oxygen_bottom'), full.names = TRUE, pattern=".grd") 
  years <- seq(1980,2100,1)
  
    #---Load in time series of age-0 recruitment abundance estimates from stock assessment----
  #Recruitment drives population growth
  pop <-read.csv('~/PROJECTS/WRAP Location/Groundfish_2019_AssessmentRecruits.csv', header = TRUE)
  pop <- pop[,c(1,6,10,11)]
  pop$year <- pop$Year-1979
  # plot(pop$Year,pop$Age.0,type='l')
  # means <- as.data.frame(cbind(pop$Year, rep(mean(pop$Age.0),length(pop$Year))))
  # lines(means$V1,means$V2,type='l',col='blue')
  
  #Need to extend time series to 2100 by drawing from either:
  # 1)average-to-good recruitment years or 2) average to bad recruitment years over the next 80 years in 20-year cycles
  #Step 1. Group lower and higher recruitment cycles - both get 'average' years
  lower <- subset(pop,pop$Category=="bad"|pop$Category=="average")
  lower.mean <- mean(lower$Age.0.logged)
  lower.sd <- sd(lower$Age.0.logged)
  
  higher <- subset(pop,pop$Category=="good"|pop$Category=="average")
  higher.mean <- mean(higher$Age.0.logged)
  higher.sd <- sd(higher$Age.0.logged)
  
  #plots of these distributions
  # plot(density(lower$Age.0),ylim=c(0,0.00015))
  # lines(density(rlnorm(10000, lower.mean, lower.sd)),col='blue')
  # lines(density(higher$Age.0),col='red')
  # lines(density(rlnorm(10000, higher.mean, higher.sd)),col='green')

  #2. Draw from these distributions to add to the Age0 time series for the next four time period cycles
  Age.0p <- pop$Age.0
  zero <- rlnorm(21, lower.mean, lower.sd) #1980-2000
  half <-rlnorm(20, higher.mean, higher.sd) #2001-2020
  first <- rlnorm(20, lower.mean, lower.sd) #2021-2040
  second <-rlnorm(20, higher.mean, higher.sd) #2041-2060
  third <-rlnorm(20, lower.mean, lower.sd) #2061-2080
  fourth <-rlnorm(20, higher.mean, higher.sd) #2061-2100
  Age.0p <- c(zero,half,first,second,third,fourth)
  Age.0p <- as.data.frame(cbind(years,Age.0p))
  Age.0p$year <- Age.0p$years-1979
  # plot(Age.0p$years,Age.0p$Age.0p, type='l',xlab="years",ylab="Age-0 abundance")
  # lines(Age.0p$years[1:39],Age.0p$Age.0p[1:39],col='blue')
  Age.0p$Age.0p_norm <- BBmisc::normalize(Age.0p$Age.0p,method = "range",range=c(0.7,1)) #use this range to vary prevalence according to popn. dynamics below. KA: using this range creates prevalence values that are >1 down below which throws an error. I've ifelse'd this to simply use "1" when >1, but this doesn't allow for boom years to be accounted for properly???
  # plot(Age.0p$years,Age.0p$Age.0p_norm, type='l',xlab="years",ylab="Age-0 normalized abundance")
  # lines(Age.0p$years[1:39],Age.0p$Age.0p_norm[1:39],col='blue')

  
  #----Loop through each year----
  for (y in 1:121){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    
    #----Load in environmental rasters for a specific year
    stemp <- raster(files_stemp[y])
    btemp <- raster(files_btemp[y])
    O2 <- raster(files_O2[y])
    z <- raster("~/Dropbox/WRAP Location^3/Rasters_2d_Spring/had/bottom_layer_depth.grd")
    # rug <- focal(z, w=matrix(1,nrow=7,ncol=7), fun=sd, na.rm=TRUE, pad=T)
    
    if(y<=2){
      stemp_t2 <- raster(files_stemp[y])
      btemp_t2 <- raster(files_btemp[y])
      O2_t2 <- raster(files_O2[y])
    } else{
      stemp_t2 <- raster(files_stemp[y-2])
      btemp_t2 <- raster(files_btemp[y-2])
      O2_t2 <- raster(files_O2[y-2])
    }
    
    #----create a mask of a smaller domain-----
    #This is to force OM to operate within a domain appropriate to sablefish.
    #Two options here: 1) Use bottom trawl survey area or 2) Use a generalized constrained geographic area
    #Option 1 code:
    # survey <- sf::st_read("~/Downloads/Groundfish OM/resurveysamplingareashapefile/WCGBTS_Grid_v2008.shp")
    # survey <- st_transform(survey,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    # survey <- st_union(survey)
    # survey <- st_sf(survey)
    # survey <- as(st_geometry(survey),"Spatial")
    
    #Moving forward with option 2:
    #This is to force OM to operate within a domain appropriate to anchovy
    scb_coords=matrix(c(-126,48,
                        -115,48,
                        -115,30,
                        -119,30,
                        -126,40),ncol=2,byrow = T)
    p=Polygon(scb_coords)
    ps=Polygons(list(p),1)
    sps = SpatialPolygons(list(ps))
    # plot(sps)
    proj4string(sps) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    data = data.frame(f=99.9)
    test = SpatialPolygonsDataFrame(sps,data)
    r <- raster(ncol=185,nrow=180)
    extent(r) <- extent(z)
    rt <- rasterize(test,r)
    # plot(rt,xlim=c(-134,-115),ylim=c(30,48),col="light blue")
    # map('worldHires',add=T)
    
    # rt <- z
    # rt[rt>1800] <- NA
    
    stemp <- raster::mask(stemp,rt)
    btemp <- raster::mask(btemp,rt)
    O2 <- raster::mask(O2,rt)
    stemp_t2 <- raster::mask(stemp_t2,rt)
    btemp_t2 <- raster::mask(btemp_t2,rt)
    O2_t2 <- raster::mask(O2_t2,rt)
    z <- raster::mask(z,rt)
    # rug <- raster::mask(rug,rt)
    #Optional: plot environmental layers
    # par(mfrow=c(1,3),mar=c(5, 4, 4, 5))
    # plot(btemp, xlim=c(-126,-116),ylim=c(30,48), main= paste0('Bottom temperature ',years[y]))
    # plot(O2, xlim=c(-126,-116),ylim=c(30,48), main = paste0('Bottom oxygen ',years[y]))
    # plot(z, xlim=c(-126,-116),ylim=c(30,48), main = 'Bathymetry')

    #----assign response curves----
    #Stack rasters
    spA_stack <- stack(z,btemp,O2)
    names(spA_stack) <- c("z","btemp","O2")
    
    spA_stack_t2 <- stack(z,btemp_t2,O2_t2)
    names(spA_stack_t2) <- c("z","btemp_t2","O2_t2")

    #Assign preferences
    spA_parameters <- formatFunctions(z = c(fun = 'dnorm', mean = 900 , sd = 1600),
                                      btemp = c(fun = 'dnorm',mean = 4 , sd = 3),
                                      O2 = c(fun = 'dnorm',mean = 57 , sd = 62))
    spA_parameters_t2 <- formatFunctions(z = c(fun = 'dnorm', mean = 900 , sd = 1600),
                                      btemp_t2 = c(fun = 'dnorm',mean = 4 , sd = 3),
                                      O2_t2 = c(fun = 'dnorm',mean = 57 , sd = 62))
                                     
    spA_suitability <- generateSpFromFun(spA_stack,
                                         parameters=spA_parameters, 
                                         rescale = FALSE,
                                         rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop. 
    spA_suitability_t2 <- generateSpFromFun(spA_stack_t2,
                                         parameters=spA_parameters_t2, 
                                         rescale = FALSE,
                                         rescale.each.response = FALSE)
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability) #plot response curves
    # #manually rescale -
    
    #logit equation for rescaling btemp and O2: exp(at*x^2 + bt*x + ct)/(1+exp(at*x^2 + bt*x + ct))
    # ref_max_z <- 0.9168028 #max value of z will always be the same since it doesn't vary annually
    # ref_max_btemp <- exp(spA_parameters$btemp$args[1]*(btemp@data@min)^2 + spA_parameters$btemp$args[2]*btemp@data@min + spA_parameters$btemp$args[3])/(1+exp(spA_parameters$btemp$args[1]*(btemp@data@min)^2 + spA_parameters$btemp$args[2]*btemp@data@min + spA_parameters$btemp$args[3])) #KA: this only works because probability increases with decreases in temperature so minimum value works here
    # ref_max_O2 <- exp(spA_parameters$O2$args[1]*(O2@data@min)^2 + spA_parameters$O2$args[2]*O2@data@min + spA_parameters$O2$args[3])/(1+exp(spA_parameters$O2$args[1]*(O2@data@min)^2 + spA_parameters$O2$args[2]*O2@data@min + spA_parameters$O2$args[3])) #KA: this only works because probability increases with decreases in O2 so minimum value works here
    ref_max_z <- dnorm(spA_parameters$z$args[1], mean=spA_parameters$z$args[1], sd=spA_parameters$z$args[2])
    ref_max_btemp <- dnorm(spA_parameters$btemp$args[1], mean=spA_parameters$btemp$args[1], sd=spA_parameters$btemp$args[2])
    ref_max_O2 <-dnorm(spA_parameters$O2$args[1], mean=spA_parameters$O2$args[1], sd=spA_parameters$O2$args[2])
    ref_max <-  ref_max_btemp * ref_max_O2  * ref_max_z 
    spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
      spA_suitability_t2$suitab.raster <- (1/ref_max)*spA_suitability_t2$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # print(spA_suitability$suitab.raster)
    # print(spA_suitability_t2$suitab.raster)
    # mean(spA_suitability_t2$suitab.raster@data@values, na.rm=T)
    # plot(spA_suitability$suitab.raster, xlim=c(-126,-116),ylim=c(30,48), main=paste0("Habitat Suitability ",years[y])) #plot habitat suitability
    # plot(spA_suitability_t2$suitab.raster, xlim=c(-126,-116),ylim=c(30,48), main=paste0("Habitat Suitability ",years[y])) #plot habitat suitability
    
    #----Convert suitability to Presence-Absence----
    #Use a specific function to convert suitability (0-1) to presence or absence (1 or 0)
    set.seed(y)
    suitability_PA <- virtualspecies::convertToPA(spA_suitability, 
                                                  PA.method = "probability",
                                                  prob.method = "logistic", 
                                                  beta = 0.4, 
                                                  alpha = -0.07, 
                                                  species.prevalence = NULL, 
                                                  plot = FALSE)
    set.seed(y)
    suitability_PA_t2 <- virtualspecies::convertToPA(spA_suitability_t2, 
                                                  PA.method = "probability",
                                                  prob.method = "logistic", 
                                                  beta = 0.4, 
                                                  alpha = -0.07, 
                                                  species.prevalence = NULL, 
                                                  plot = FALSE)
    # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    # plot(suitability_PA$pa.raster)

    
    #-----Sample Presences and generate Biomass-----  
    # set.seed(y)
    # find_prevalence <- sampleOccurrences(suitability_PA,
    #                                      n = 500,
    #                                      type = "presence-absence",
    #                                      detection.probability = 1,
    #                                      error.probability=0,
    #                                      plot = FALSE)
    #                                      # sampling.area = survey
    #                                      # replacement = TRUE)
    # suit_prev <- as.numeric(find_prevalence$sample.prevalence[1])
    # pop_prev <- ifelse(y==1,Age.0p$Age.0p_norm[Age.0p$year==y],Age.0p$Age.0p_norm[Age.0p$year==(y-1)]) # the prevalence of sampling adult sablefish is based on recruitment estimates from year t-1
    # sp_prev <- round(pop_prev * suit_prev,1)
    
    set.seed(y)
    presence.points <- sampleOccurrences(suitability_PA,
                                         n = 500,
                                         type = "presence-absence",
                                         detection.probability = 1,
                                         error.probability=0,
                                         plot = FALSE)
                                         # sample.prevalence = sp_prev)
                                         # sampling.area = survey,
                                         # replacement = TRUE)
    
    df <- cbind(as.data.frame(presence.points$sample.points$x),as.data.frame(presence.points$sample.points$y))
    colnames(df) <- c("lon","lat")
    df$sampled <- 1
    
    #----Extract data for each year----
    print("Extracting suitability")
    ei <- gridcells*y #end location in output grid to index to
    se <- ei - (gridcells-1) #start location in output grid to index to
    output$lat[se:ei] <- rasterToPoints(btemp)[,2]
    output$lon[se:ei] <- rasterToPoints(btemp)[,1]
    output$year[se:ei] <- rep(years[y],gridcells)
    output$pres[se:ei] <- rasterToPoints(suitability_PA$pa.raster)[,3] 
    output$suitability[se:ei] <- rasterToPoints(suitability_PA$suitab.raster)[,3] 
    output$suitability_t2[se:ei] <- rasterToPoints(suitability_PA_t2$suitab.raster)[,3]
    output$btemp[se:ei] <-  rasterToPoints(btemp)[,3]
    output$O2[se:ei] <- rasterToPoints(O2)[,3]
    output$z[se:ei] <-  rasterToPoints(z)[,3]
    output$stemp[se:ei] <-  rasterToPoints(stemp)[,3]
    #temp file
    temp_sampled <- left_join(output[se:ei,1:8], df, by=c('lon','lat'))
    #write to output
    output$sampled[se:ei] <- temp_sampled$sampled
    output$pop_notnorm[se:ei] <- Age.0p$Age.0p[Age.0p$year==y]
    output$pop_norm[se:ei] <-  Age.0p$Age.0p_norm[Age.0p$year==y]
    if(y<=2){
      output$pop_norm_t2[se:ei] <-  Age.0p$Age.0p_norm[Age.0p$year==y]
    } else {
      output$pop_norm_t2[se:ei] <-  Age.0p$Age.0p_norm[Age.0p$year==y-2]
    }
    
  }
  
  #Convert NAs in 'sampled' columns to zeros
  output$sampled <- ifelse(is.na(output$sampled),0, output$sampled)
  
  #----Create abundance as a function of sablefish total biomass estimates from the stock assessment and the environment----
  
  #For abundance I see two options: 1) Do as Steph did for Anchovy using the total biomass estimates from the stock assessment or 2) Draw from log-normal distribution based on numbers that Eric Ward suggested.
  
  #1: Use total biomass estimates mean and sd: a couple choices: 1) using the entire time series gives a mean of ~228K mt, but the estimates from the whole timeseries (mean=228K mt) and other periods (1890-1979: 271K mt; 1980-1999: 278K mt) are much higher than the most recent 20-year period (2000-2019: 178K mt).
  # ssb <- read.csv('Sablefish abundance time series from stock assessment.csv', header=TRUE)
  # a_presences <- nrow(filter(output,pres==1)) #this identifies the number of grid cells with 'presences' in which to average the biomass across
  
  a_mean <- log((228000 / gridcells )*3) #divide by number of grid cells, but multiply by 2 to offset multiplying by suitablity & pop below
  a_sd <- log((85152 / gridcells)/17) 
  # hist(rlnorm(gridcells,a_mean,a_sd))
  # # a_mean <- 258256 /150
  # a_sd <- 85152 / 150 / 5
  # output$abundance_2 <- ifelse(output$pres==1,rlnorm(nrow(output),a_mean, a_sd)*output$suitability,0)
  
  
  
  # output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),a_mean, a_sd)*output$suitability*output$pop_norm,0)
  output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),a_mean, a_sd)*output$suitability_t2*output$pop_norm_t2,0)
  # sum(output$abundance, na.rm=T)
  # output$abundance_2 <- output$abundance * output$pop_prev

  return(output)
}

