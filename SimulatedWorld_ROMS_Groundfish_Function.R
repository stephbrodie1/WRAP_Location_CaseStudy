#This code is the operating model for a groundfish archetype, likely represented by sablefish
#It uses average spring conditions from downscales ROMS projections using GFDL. 
#IMPORTANT: Download average spring ROMS data from here: https://www.dropbox.com/sh/pj1vjz4h1n27hnb/AAD9sySDudMKmq3l4yLPvjQra?dl=0

#We include a recruitment suitability index of a preceding year (t-1) that influences distribution of sablefish in the subsequent year (t):
#recruitment: suitability driven by temperature, oxygen, and an SSH index
#adult: distribution and abundance driven by bottom temp, bottom oxygen, habitat covariate, and recruitment at t-1.  

SimulateWorld_ROMS_Groundfish <- function(dir){
  #dir is the local directory that points to where ROMS data is stored 
  
  #required libraries
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  
  #----Create output file----
  #Assuming 400 'samples' are taken each year, from 1980-2100
  output <- as.data.frame(matrix(NA, nrow=48400,ncol=7))
  colnames(output) <- c("lon","lat","year","pres","suitability","sst", "chla")
  
  #----Load in rasters----
  #Note that these need to be updated with correct covariates. These are all contained in the dropbox link above. 
  gcm_dr <- 'gfdl'
  files_sst <- list.files(paste0(dir,'gfdl/sst_monthly'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_chl <- list.files(paste0(dir,'gfdl/chl_surface'), full.names = TRUE, pattern=".grd") #should be 1452 files
  years <- seq(1980,2100,1)
  
  #----loop through each year----
  for (y in 1:121){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    
    #Environmental at time t
    sst_t <- raster(files_sst[y])
    chla_t <- raster(files_chl[y])
    chla_t <- log(chla_t)
    # plot(sst_t)
    # plot(chla_t)
    
    #Environmental at time t-1 (previous year). 
    #If year 1, then just build recruitment at time t. 
    if (y>1){
    sst_t1 <- raster(files_sst[y-1])
    chla_t1 <- raster(files_chl[y-1])
    chla_t1 <- log(chla_t1)
    } else {
      sst_t1 <- sst_t
      chla_t1 <- chla_t
    }
    # plot(sst_t1)
    # plot(chla_t1)
    
    #insert SSH index here too. An example of how to create a raster that contains a single value. 
    ssh_index <- sst_t1
    ssh_index[] <- 1
    #Insert habitat covariate raster here, and make sure extent is the same so it can be stacked
    
    #----Recruitment Index: assign response curves----
    
    #Stack rasters
    spA_stack <- stack(sst_t1, chla_t1)
    names(spA_stack) <- c('sst_t1', 'chla_t1')
    #Assign preferences
    spA_parameters <- formatFunctions(sst_t1 = c(fun="dnorm",mean=12,sd=3),
                                      chla_t1 = c(fun="dnorm",mean=1.6,sd=9))
    spA_suitability <- generateSpFromFun(spA_stack,parameters=spA_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop. 
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spA_parameters$sst_t1$args[1], mean=spA_parameters$sst_t1$args[1], sd=spA_parameters$sst_t1$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_chl <- dnorm(spA_parameters$chla_t1$args[1], mean=spA_parameters$chla_t1$args[1], sd=spA_parameters$chla_t1$args[2]) 
    ref_max <- ref_max_sst * ref_max_chl #simple multiplication of layers. 
    spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    
    #----Adult groundfish: assign response curves----
    
    #Stack rasters
    spB_stack <- stack(sst_t, spA_suitability$suitab.raster)
    names(spB_stack) <- c('sst_t', 'spA')
    
    #Assign preferences
    spB_parameters <- formatFunctions(sst_t = c(fun="dnorm",mean=15,sd=4),
                                      spA = c(fun="dnorm",mean=0.8,sd=1))
    spB_suitability <- generateSpFromFun(spB_stack,parameters=spB_parameters, rescale = FALSE,rescale.each.response = FALSE)
    # plot(spB_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spB_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spB_parameters$sst_t$args[1], mean=spB_parameters$sst_t$args[1], sd=spB_parameters$sst_t$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_spA <- dnorm(spB_parameters$spA$args[1], mean=spB_parameters$spA$args[1], sd=spB_parameters$spA$args[2])
    ref_max <- ref_max_sst * ref_max_spA
    spB_suitability$suitab.raster <- (1/ref_max)*spB_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    # plot(spB_suitability$suitab.raster) #plot habitat suitability
    
    #----Convert suitability to Presence-Absence----
    suitability_PA <- virtualspecies::convertToPA(spB_suitability, PA.method = "probability", beta = "random",
                                                  alpha = -0.5, species.prevalence = 0.5, plot = FALSE)
    # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    # plot(suitability_PA$pa.raster)
    
    #-----Sample Presences and Absences-----
    presence.points <- sampleOccurrences(suitability_PA,n = 400,type = "presence-absence", 
                                         detection.probability = 1,error.probability=0, plot = FALSE) #default but cool options to play with                                    )
    df <- cbind(as.data.frame(round(presence.points$sample.points$x,1)),as.data.frame(round(presence.points$sample.points$y,1)))
    colnames(df) <- c("x","y")
    
    #----Extract data for each year----
    print("Extracting suitability")
    ei <- 400*y #end location in output grid to index to
    se <- ei - 399 #start location in output grid to index to
    output$lat[se:ei] <- df$y
    output$lon[se:ei] <- df$x
    output$year[se:ei] <- rep(years[y],400)
    output$pres[se:ei] <-  presence.points$sample.points$Real
    output$suitability[se:ei] <- raster::extract(spB_suitability$suitab.raster, y= df)  #extract points from suitability file
    output$sst[se:ei] <-  raster::extract(sst_t, y= df)  #extract points from suitability file
    output$chla[se:ei] <-  raster::extract(chla_t, y= df)
  }
  
  #----Create abundance as a function of the environment----
  # SB: values in Ecography paper. Initially based on EBS flounder.
  #parameters need to be updated with species specific info
  output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),2,0.1)*output$suitability,0)
  
  return(output)
}

#----example-----
test <- SimulateWorld_ROMS_Groundfish(dir = "~/Dropbox/WRAP Location^3/Rasters_2d_Spring/" )
