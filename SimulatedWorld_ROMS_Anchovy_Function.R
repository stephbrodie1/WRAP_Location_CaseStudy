#This code is the operating model for a CPS archetype, represented by anchovy
#It uses average spring conditions from downscales ROMS projections using Hadley 
#IMPORTANT: Download average spring ROMS data from here: https://www.dropbox.com/sh/pj1vjz4h1n27hnb/AAD9sySDudMKmq3l4yLPvjQra?dl=0

#We include boom-bust population dynamics from simulations provided by Andre Punt

SimulateWorld_ROMS_Anchovy <- function(dir){
  #dir is the local directory that points to where ROMS data is stored 
  
  #required libraries
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  
  #----Create output file----
  #Assuming 400 'samples' are taken each year, from 1980-2100
  #This will be the information passed to the estimation model
  output <- as.data.frame(matrix(NA, nrow=48400,ncol=10))
  colnames(output) <- c("lon","lat","year","pres","suitability","sst","chla_surface", "z","zoo", "abundance")
  
  #---Load in population simulation ----
  #A MICE model for Anchovy biomass in the CCS. 50 simulations, 2000 years each
  #From Andre Punt: e.g. https://www.sciencedirect.com/science/article/pii/S0304380016302034?via%3Dihub
  #Ignore first 50years (burn in), so just grab last 121 years of first simulation
  #Need some kind of error, so getting standard deviation across the 50 simulations for the last 121 years
  pop <- read.csv('~/PROJECTS/WRAP Location/SUMMARYApreycleanversion.csv', header = TRUE)
  pop <- pop[pop$Year>=1880,]
  pop$year <- pop$Year - 1879
  error <- sd(pop$AnchovySSB) #pretty large
  pop <- pop[pop$Sim==1,]
  # plot(pop$Year,pop$AnchovySSB, type='l')
  
  #----Load in rasters----
  #These are the average spring conditions from the downscaled Hadley earth system model
  files_sst <- list.files(paste0(dir,'sst_monthly'), full.names = TRUE, pattern=".grd") 
  files_chl <- list.files(paste0(dir,'chl_surface'), full.names = TRUE, pattern=".grd") 
  files_zoo200 <- list.files(paste0(dir,'zoo_200m'), full.names = TRUE, pattern=".grd") 
  years <- seq(1980,2100,1)
  
  #----Loop through each year----
  for (y in 1:121){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    
    #Load in environmental rasters for a specific year
    sst <- raster(files_sst[y])
    chla <- raster(files_chl[y])
    chla <- log(chla)
    zoo <- raster(files_zoo200[y])
    z <- raster('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/gfdl/bottom_layer_depth.grd')
    z <- z * -1
    
    #Optional: plot environmental layers
    # plot(sst, main= 'SST')
    # plot(chla, main = 'Surface Chl-a')
    # plot(z, main = 'Bathymetry')
    # plot(zoo)
    
    #----assign response curves----
    #Stack rasters
    spA_stack <- stack(sst, zoo, z)
    names(spA_stack) <- c('sst', 'zoo', 'z')
    
    #Assign preferences
    spA_parameters <- formatFunctions(sst = c(fun="dnorm",mean=15,sd=7),
                                      # chla = c(fun="dnorm",mean=0.5,sd=1.5),
                                      zoo = c(fun="logisticFun",alpha=-6,beta=50),
                                      z = c(fun="logisticFun",alpha=-250,beta=-2000))
    spA_suitability <- generateSpFromFun(spA_stack,parameters=spA_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability) #plot response curves

    #manually rescale
    ref_max_sst <- dnorm(spA_parameters$sst$args[1], mean=spA_parameters$sst$args[1], sd=spA_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spA_parameters$zoo$args[2])/spA_parameters$zoo$args[1]))
    ref_max_z <- 1 / (1 + exp(((z@data@max) - spA_parameters$z$args[2])/spA_parameters$z$args[1]))
    ref_max <- ref_max_sst * ref_max_zoo #* ref_max_z
    spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # print(spA_suitability$suitab.raster)
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    
    #----Convert suitability to Presence-Absence----
    #Use a specific function to convert suitability (0-1) to presence or absence (1 or 0)
    # prev <- pop$AnchovySSB[pop$year==y]
    # prev <- (pop$AnchovySSB - min(pop$AnchovySSB)) / (max(pop$AnchovySSB) - min(pop$AnchovySSB))
    suitability_PA <- virtualspecies::convertToPA(spA_suitability, PA.method = "probability", beta = 0.5,
                                                  alpha = -0.05, species.prevalence = NULL, plot = TRUE)
    # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    # plot(suitability_PA$pa.raster)
    
    #-----Sample Presences and generate Biomass-----
    prev <- pop$AnchovySSB[pop$year==y]
    prev <- (pop$AnchovySSB - min(pop$AnchovySSB)) / (max(pop$AnchovySSB) - min(pop$AnchovySSB))
    presence.points <- sampleOccurrences(suitability_PA,n = 400,type = "presence-absence",
                                         detection.probability = 1,error.probability=0, plot = FALSE,
                                         # sample.prevalence = 0.5)
                                         sample.prevalence = round(prev[y],1)) 
    df <- cbind(as.data.frame(round(presence.points$sample.points$x,1)),as.data.frame(round(presence.points$sample.points$y,1)))
    colnames(df) <- c("x","y")
    
    #----Extract data for each year----
    print("Extracting suitability")
    ei <- 400*y #end location in output grid to index to
    se <- ei - 399 #start location in output grid to index to
    output$lat[se:ei] <- df$y
    output$lon[se:ei] <- df$x
    output$year[se:ei] <- rep(years[y],400)
    output$pres[se:ei] <- presence.points$sample.points$Real
    output$suitability[se:ei] <- raster::extract(spA_suitability$suitab.raster, y= df[,1:2])  #extract points from suitability file
    output$sst[se:ei] <-  raster::extract(sst, y= df[,1:2])  #extract points from suitability file
    output$chla_surface[se:ei] <-  raster::extract(chla, y= df[,1:2])
    output$zoo[se:ei] <-  raster::extract(zoo, y= df[,1:2])
    output$z[se:ei] <-  raster::extract(z, y= df[,1:2])
    
  }
  
  a_mean <- mean(pop$AnchovySSB) / 120
  a_sd <- (sd(pop$AnchovySSB) / 120) / 3 
  output$abundance <- ifelse(output$pres==1,rnorm(nrow(output),a_mean, a_sd)*output$suitability,0)

  return(output)
}
