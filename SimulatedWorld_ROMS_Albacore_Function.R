#This code is the operating model for a HMS archetype, represented by north pacific albacore tuna
#It uses average spring conditions from downscales ROMS projections using Hadley 
#IMPORTANT: Download average spring ROMS data from here: https://www.dropbox.com/sh/pj1vjz4h1n27hnb/AAD9sySDudMKmq3l4yLPvjQra?dl=0

#We include a trophic interaction between predator (albacore) and prey, but note the estmation model will use chl-a as a proxy for prey information:
#prey: distribution and suitability driven by SST and zooplankton integrated across top 200m
#predator: distribution and abundance driven by SST, MLD, and Species A

SimulateWorld_ROMS_Albacore <- function(dir){
  #dir is the local directory that points to where ROMS data is stored 

  #required libraries
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  
  #----Create output file----
  nsamples <- 500  #number of samples to use in the estimation model
  gridcells <- 21912
  output <- as.data.frame(matrix(NA, nrow=(gridcells*121), ncol=10))  #21912 non-NA grid cells in ROMS
  colnames(output) <- c("lon","lat","year","pres","suitability","sst","zoo_200","chla_surface", "mld", "sample")

  #----Load in rasters----
  #These are the average spring conditions from the downscaled gfdl earth system model
  files_sst <- list.files(paste0(dir,'sst_monthly'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_chl_surface <- list.files(paste0(dir,'chl_surface'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_mld <- list.files(paste0(dir,'ild_0.5C'), full.names = TRUE, pattern=".grd")
  files_zoo <- list.files(paste0(dir,'zoo_200m'), full.names = TRUE, pattern=".grd")
  years <- seq(1980,2100,1)
  
  #----Loop through each year----
  for (y in 1:121){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    
    #Load in environmental rasters for a specific year
    sst <- raster(files_sst[y])
    mld <- raster(files_mld[y])
    zoo <- raster(files_zoo[y])
    chla_surface <- raster(files_chl_surface[y])
    chla_surface <- log(chla_surface)
    
    #Optional: plot environmental layers
    # par(mfrow=c(1,4))
    # plot(sst, main= 'SST')
    # plot(chla_surface, main = 'Chl-a Surface')
    # plot(mld, main = 'MLD')
    # plot(zoo, main = 'Zooplankton')


    #----SPECIES A (prey): assign response curves----
    #Species A: likes high zooplankton and medium temps
    
    #Stack rasters
    spA_stack <- stack(sst, zoo)
    names(spA_stack) <- c('sst', 'zoo')
    
    #Assign preferences

    spA_parameters <- formatFunctions(sst = c(fun="dnorm",mean=14,sd=7),
                                      zoo = c(fun="logisticFun",alpha=-10,beta=45))
    spA_suitability <- generateSpFromFun(spA_stack,parameters=spA_parameters,
                                         rescale = FALSE,rescale.each.response = FALSE,
                                         species.type="multiplicative") #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability) #plot response curves

    #manually rescale
    ref_max_sst <- dnorm(spA_parameters$sst$args[1], mean=spA_parameters$sst$args[1], sd=spA_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spA_parameters$zoo$args[2])/spA_parameters$zoo$args[1]))
    ref_max <- ref_max_sst * ref_max_zoo #simple multiplication of layers.
    spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # spA_suitability$suitab.raster
    # plot(spA_suitability$suitab.raster) #plot habitat suitability

    #----SPECIES B (albacore): assign response curves----
    #Species B: likes to eat Species A, and warmer temperatues & shallow MLD
    
    #Stack rasters
    spB_stack <- stack(sst, mld, spA_suitability$suitab.raster)
    names(spB_stack) <- c('sst',"mld", "spA")
    
    #Assign preferences
    spB_parameters <- formatFunctions(sst = c(fun="dnorm",mean=17,sd=4),
                                      mld = c(fun="dnorm",mean=50,sd=30),
                                      # spA = c(fun="dnorm",mean=0.5,sd=2))
                                      spA = c(fun="logisticFun",alpha=-0.15,beta=0.4))
    spB_suitability <- generateSpFromFun(spB_stack,parameters=spB_parameters,
                                         rescale = FALSE,rescale.each.response = FALSE,
                                         species.type = "multiplicative")
    # plot(spB_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spB_suitability) #plot response curves

    #manually rescale
    ref_max_sst <- dnorm(spB_parameters$sst$args[1], mean=spB_parameters$sst$args[1], sd=spB_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_mld <- dnorm(spB_parameters$mld$args[1], mean=spB_parameters$mld$args[1], sd=spB_parameters$mld$args[2])
    ref_max_spA <- 1 / (1 + exp(((spA_suitability$suitab.raster@data@max) - spB_parameters$spA$args[2])/spB_parameters$spA$args[1]))
    ref_max <- ref_max_sst * ref_max_mld * ref_max_spA
    spB_suitability$suitab.raster <- (1/ref_max)*spB_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    # print(spB_suitability$suitab.raster)
    # plot(spB_suitability$suitab.raster) #plot habitat suitability

        #----Convert suitability to Presence-Absence----
    #Use a specific function to convert suitability (0-1) to presence or absence (1 or 0)
    
    suitability_PA <- virtualspecies::convertToPA(spB_suitability, PA.method = "probability", beta = 0.4,
                                                  alpha = -0.07, species.prevalence = NULL, plot = FALSE)
    # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    # plot(suitability_PA$pa.raster)
    
    #-----Sample Presences and Absences-----
    presence.points <- sampleOccurrences(suitability_PA,n = nsamples,type = "presence-absence",
                                         detection.probability = 1,error.probability=0, plot = FALSE,
                                         sample.prevalence = NULL)
    #convert to dataframe
    pres_df <- cbind(as.data.frame(presence.points$sample.points$x),as.data.frame(presence.points$sample.points$y))
    colnames(pres_df) <- c("x","y")
    pres_df$sample <- 1
    
    #expand dataframe to include all possible locations
    df_full <- as.data.frame(rasterToPoints(sst)[,1:2])
    df_full_2 <- left_join(df_full, pres_df, by=c('x','y'))
    df_full_2$sample <- ifelse(is.na(df_full_2$sample),0,df_full_2$sample)
    df_full_2$lon <- df_full_2$x
    df_full_2$lat <- df_full_2$y
    
    #----Extract data for each year----
    print("Extracting suitability")
    ei <- gridcells*y #end location in output grid to index to
    se <- ei - (gridcells-1) #start location in output grid to index to
    output$lat[se:ei] <- rasterToPoints(sst)[,2] #using sst as an example raster 
    output$lon[se:ei] <- rasterToPoints(sst)[,1]
    output$year[se:ei] <- rep(years[y],gridcells)
    output$pres[se:ei] <- rasterToPoints(suitability_PA$pa.raster)[,3] 
    output$suitability[se:ei] <- rasterToPoints(spB_suitability$suitab.raster)[,3]  #extract points from suitability file
    output$sst[se:ei] <-  rasterToPoints(sst)[,3]  #extract points from suitability file
    output$zoo_200[se:ei] <-  rasterToPoints(zoo)[,3]
    output$chla_surface[se:ei] <-  rasterToPoints(chla_surface)[,3]
    output$mld[se:ei] <-  rasterToPoints(mld)[,3]
    temp <-   left_join(output[se:ei,1:9], df_full_2[,c(3,4,5)], by=c('lon','lat'))
    output$sample[se:ei] <- temp$sample
  }
  
  #Average monthly biomass available to CCS is: 1.18x10^5 ± (0.13x10^5 se) mt (from Desiree Tommasi)
  mean_spatial <- log((118000/gridcells)*5) #Kind of had to make these up to get biomass pop that matches 1.18, and with error low enough that the EM does ok. 
  se_spatial <- log((13000)/10000) 
  output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),mean_spatial, se_spatial)*output$suitability,0)
  # hist(output$abundance)
  # sum(output$abundance, na.rm=T)
  # 
  print('Saving csv to working directory')
  # write.csv(output, 'Albacore_OM_Simulation.csv',row.names = FALSE)
  return(output)
}

