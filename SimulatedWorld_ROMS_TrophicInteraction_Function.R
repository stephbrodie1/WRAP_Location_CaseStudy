#This code replicates SimulatedWorld_ROMS_function, but adds in a trophic interaction, where:
#Species A: distribution and abundance drivn by SST and Chl-a
#Species B: distribution and abundance drivn by SST and Species A
#EM for Species b: only have chl-a and temp as covariates. 

#Note only using GFDL for now

SimulateWorld_ROMS_TrophicInteraction <- function(PA_shape, abund_enviro, dir){
  # # 'PA_shape' specifies how enviro suitability determines species presence-absence...
  # #... takes values of "logistic" (SB original), "logistic_prev" (JS, reduces knife-edge), "linear" (JS, reduces knife edge, encourages more absences)
  # # 'abund_enviro' specifies abundance if present, can be "lnorm_low" (SB original)...
  # #... "lnorm_high" (EW), or "poisson" (JS, increases abundance range)
  # # 'dir' points to your directory where ROMS raster are kept
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  
  #----IMPORTANT----
  #Download ROMS data from here: https://www.dropbox.com/sh/aaezimxwq3glwdy/AABHmZbmfjVJM7R4jcHCi4c9a?dl=0
  
  # #----Create output file----
  #####Needs to be modified as variables are added. Starting with sst
  #Assuming 400 'samples' are taken each year, from 1980-2100
  output <- as.data.frame(matrix(NA, nrow=48400,ncol=7))
  colnames(output) <- c("lon","lat","year","pres","suitability","sst", "chla")

    #----Load in rasters----
  gcm_dr <- 'gfdl'

  files_sst <- list.files(paste0(dir,'gfdl/sst_monthly'), full.names = TRUE, pattern=".grd") #should be 1452 files
  files_chl <- list.files(paste0(dir,'gfdl/chl_surface'), full.names = TRUE, pattern=".grd") #should be 1452 files
  months <- rep(1:12,121) 
  years <- rep(1980:2100,each=12)
  august_indexes <- which(months==8) #picking last month in summer
  
  #loop through each year
  for (y in august_indexes){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    
    sst <- raster(files_sst[y])
    chla <- raster(files_chl[y])
    chla <- log(chla)
    # plot(sst)
    # plot(chla)
    # 
    #----SPECIES A: assign response curve----
    #Species A: likes high chla and medium temps
    
    #Stack rasters
    spA_stack <- stack(sst, chla)
    names(spA_stack) <- c('sst', 'chla')
    #Assign preferences
    spA_parameters <- formatFunctions(sst = c(fun="dnorm",mean=12,sd=3),
                                                  chla = c(fun="dnorm",mean=1.6,sd=9))
    spA_suitability <- generateSpFromFun(spA_stack,parameters=spA_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop. 
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spA_parameters$sst$args[1], mean=spA_parameters$sst$args[1], sd=spA_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_chl <- dnorm(spA_parameters$chla$args[1], mean=spA_parameters$chla$args[1], sd=spA_parameters$chla$args[2]) 
    ref_max <- ref_max_sst * ref_max_chl #simple multiplication of layers. 
    spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
    # plot(spA_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spA_suitability) #plot response curves


    #----SPECIES B: assign response curve----
    #Species B: likes to eat Species A, and warmer temperatues
    
    #Stack rasters
    spB_stack <- stack(sst, spA_suitability$suitab.raster)
    names(spB_stack) <- c('sst', 'spA')
    
    #Assign preferences
    spB_parameters <- formatFunctions(sst = c(fun="dnorm",mean=15,sd=4),
                                      spA = c(fun="dnorm",mean=0.8,sd=1))
    spB_suitability <- generateSpFromFun(spB_stack,parameters=spB_parameters, rescale = FALSE,rescale.each.response = FALSE)
    # plot(spB_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spB_suitability) #plot response curves
    
    #manually rescale
    ref_max_sst <- dnorm(spB_parameters$sst$args[1], mean=spB_parameters$sst$args[1], sd=spB_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    ref_max_spA <- dnorm(spB_parameters$spA$args[1], mean=spB_parameters$spA$args[1], sd=spB_parameters$spA$args[2])
    # ref_max <- ref_max_sst * 0.5
    spB_suitability$suitab.raster <- (1/ref_max)*spB_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    # plot(spB_suitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(spB_suitability) #plot response curves
    
    #----Convert suitability to Presence-Absence----
    if (PA_shape == "logistic") {
      #SB: specifies alpha and beta of logistic - creates almost knife-edge absence -> presence
      suitability_PA <- virtualspecies::convertToPA(spB_suitability, PA.method = "probability", beta = 0.5,
                                                    alpha = -0.05, species.prevalence = NULL, plot = FALSE, )
      # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    }
    if (PA_shape == "logistic_prev") {
      #JS: relaxes logistic a little bit, by specifing reduced prevalence and fitting beta (test diff prevalence values, but 0.5 seems realistic)
      suitability_PA <- virtualspecies::convertToPA(envirosuitability, PA.method = "probability", beta = "random",
                                                    alpha = -0.3, species.prevalence = 0.5, plot = FALSE)
      # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    }
    if (PA_shape == "linear") {
      #JS: relaxes knife-edge absence -> presence further; also specifies prevalence and fits 'b'
      suitability_PA <- virtualspecies::convertToPA(envirosuitability, PA.method = "probability",
                                                    prob.method = "linear", a = NULL,
                                                    b = NULL, species.prevalence = 0.8, plot = FALSE)
      # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
    }
    # plot(suitability_PA$pa.raster)
    
    #-----Sample Presences and Absences-----
    presence.points <- sampleOccurrences(suitability_PA,n = 400,type = "presence-absence", 
                                         detection.probability = 1,error.probability=0, plot = FALSE) #default but cool options to play with                                    )
    df <- cbind(as.data.frame(round(presence.points$sample.points$x,1)),as.data.frame(round(presence.points$sample.points$y,1)))
    colnames(df) <- c("x","y")
    
    #----Extract data for each year----
    print("Extracting suitability")
    ei <- 400*which(august_indexes==y) #end location in output grid to index to
    se <- ei - 399 #start location in output grid to index to
    output$lat[se:ei] <- df$y
    output$lon[se:ei] <- df$x
    output$year[se:ei] <- rep(years[y],400)
    output$pres[se:ei] <-  presence.points$sample.points$Real
    output$suitability[se:ei] <- raster::extract(spB_suitability$suitab.raster, y= df)  #extract points from suitability file
    output$sst[se:ei] <-  raster::extract(sst, y= df)  #extract points from suitability file
    output$chla[se:ei] <-  raster::extract(chla, y= df)
    
  }
  
  #----Create abundance as a function of the environment----
  if (abund_enviro == "lnorm_low") {
    # SB: values in Ecography paper. I think initially they were based on flounder in EBS but not sure if I edited them
    output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),2,0.1)*output$suitability,0)
  }
  if (abund_enviro == "lnorm_high") {
    # EW: I'm cranking up the rlnorm parameters to make it more comparable to wc trawl survey estimates -- these new ones based on arrowtooth
    # SB: rlnorm parameters (6,1) were too large for estimation model. GAMs had explained deviance <10%. Other suggestions?
    output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),6,1)*output$suitability,0)
  }
  if (abund_enviro == "poisson") {
    # JS: sample from a Poisson distbn, with suitability proportional to mean (slower, bc it re-samples distbn for each observation)
    maxN <- 50  #max mean abundance at highest suitability
    output$abundance <- ifelse(output$pres==1,rpois(nrow(output),lambda=output$suitability*maxN),0)
  } 
  
  return(output)
}


#Example-
# test <- SimulateWorld_ROMS_multispecies(PA_shape="logisitic", abund_enviro="lnorm_low", dir <- "~/Dropbox/WRAP Location^3/Rasters_2d_monthly/")
