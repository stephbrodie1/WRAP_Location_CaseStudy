#Caution: this function uses the downscaled GCMs to simulate species distrubtion
#It has slight differences to the other SimulatedWorld_function that 'randomly' generates environmental data
#Remember: 1980-2010 are not observed data (by design)

SimulateWorld_ROMS <- function(PA_shape, abund_enviro, dir){
  # 'PA_shape' specifies how enviro suitability determines species presence-absence...
  #... takes values of "logistic" (SB original), "logistic_prev" (JS, reduces knife-edge), "linear" (JS, reduces knife edge, encourages more absences)
  # 'abund_enviro' specifies abundance if present, can be "lnorm_low" (SB original)...
  #... "lnorm_high" (EW), or "poisson" (JS, increases abundance range)
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  
  #----Create output file----
  #Needs to be modified as variables are added. Starting with sst
  #Assuming 400 'samples' are taken each year, from 1980-2100
  output <- as.data.frame(matrix(NA, nrow=48400,ncol=6))
  colnames(output) <- c("lon","lat","year","presabs","suitability","sst")
  
  #----Load in rasters----
  gcm_dr <- 'gfdl'

  files_sst <- list.files(paste0(dir,'gfdl/sst_monthly'), full.names = TRUE, pattern=".grd") #should be 1452 files
  months <- rep(1:12,121) 
  years <- rep(1980:2100,each=12)
  august_indexes <- which(months==8) #picking last month in summer

  #loop through each year
  for (y in august_indexes){
    print(paste0("Creating environmental simulation for Year ",years[y]))
    sst <- raster(files_sst[y])
    #plot(sst)
 
    #----Use Virtual Species to assign response curve----
    envir_stack <- stack(sst) #must be in raster stack format
    names(envir_stack) <- c('sst')
    
    #----assign response curve with mean at 4 degrees----
    parameters <- virtualspecies::formatFunctions(sst = c(fun="dnorm",mean=15,sd=4))
    
    #----convert temperature raster to species suitability----
    envirosuitability <- virtualspecies::generateSpFromFun(envir_stack,parameters=parameters, rescale = FALSE,rescale.each.response = FALSE)
    #rescale
    ref_max <- dnorm(parameters$sst$args[1], mean=parameters$sst$args[1], sd=parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    envirosuitability$suitab.raster <- (1/ref_max)*envirosuitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    
    #Plot suitability and response curve
    # plot(envirosuitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(envirosuitability) #plot response curves
    
    #----Convert suitability to Presence-Absence----
    if (PA_shape == "logistic") {
      #SB: specifies alpha and beta of logistic - creates almost knife-edge absence -> presence
      suitability_PA <- virtualspecies::convertToPA(envirosuitability, PA.method = "probability", beta = 0.5,
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
      output$presabs[se:ei] <-  presence.points$sample.points$Real
      output$suitability[se:ei] <- raster::extract(envirosuitability$suitab.raster, y= df)  #extract points from suitability file
      output$sst[se:ei] <-  raster::extract(sst, y= df)  #extract points from suitability file
  }
  
  #----Create abundance as a function of the environment----
  if (abund_enviro == "lnorm_low") {
    # SB: values in Ecography paper. I think initially they were based on flounder in EBS but not sure if I edited them
    output$abundance <- ifelse(output$presabs==1,rlnorm(48400,2,0.1)*output$suitability,0)
  }
  if (abund_enviro == "lnorm_high") {
    # EW: I'm cranking up the rlnorm parameters to make it more comparable to wc trawl survey estimates -- these new ones based on arrowtooth
    # SB: rlnorm parameters (6,1) were too large for estimation model. GAMs had explained deviance <10%. Other suggestions?
    output$abundance <- ifelse(output$presabs==1,rlnorm(48400,6,1)*output$suitability,0)
  }
  if (abund_enviro == "poisson") {
    # JS: sample from a Poisson distbn, with suitability proportional to mean (slower, bc it re-samples distbn for each observation)
    maxN <- 50  #max mean abundance at highest suitability
    output$abundance <- ifelse(output$presabs==1,rpois(48400,lambda=output$suitability*maxN),0)
  } 
  

  return(output)
}


#-----Example-----
# test <- SimulatedWorld_ROMS(PA_shape="logistic", abund_enviro="lnorm_low")
# head(test)
# tail(test)
# plot(aggregate(suitability~year, data=test, FUN=mean),type='b')
# plot(aggregate(presabs~year, data=test, FUN=sum),type='b')
# plot(aggregate(abundance~year, data=test, FUN=sum),type='b')
# plot(aggregate(sst~year, data=test, FUN=mean),type='b')

