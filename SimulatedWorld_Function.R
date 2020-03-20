#Function to simulate species distribution and abundance with respect to environmental covariates

SimulateWorld <- function(temp_diff, temp_spatial, PA_shape, abund_enviro){
  # 'temp_diff' specifies min and max temps at year 1 and year 100 (e.g. temp_diff=c(1,3,5,7) means year 1 varies from 1-3C and year 100 from 5-7C)
  # 'temp_spatial' specifies whether we have "simple" linear temp distbn (SB) or added "matern" variation (EW)
  # 'PA_shape' specifies how enviro suitability determines species presence-absence...
  #... takes values of "logistic" (SB original), "logistic_prev" (JS, reduces knife-edge), "linear" (JS, reduces knife edge, encourages more absences, currently throws errors)
  # 'abund_enviro' specifies abundance if present, can be "lnorm_low" (SB original)...
  #... "lnorm_high" (EW), or "poisson" (JS, increases abundance range)
  
  library(virtualspecies)
  library(scales)
  library(raster)
  library(glmmfields)
  #----Generate grid of locations through time----
  x_tot <- seq(1, 20, 1)
  y_tot <- seq(1, 20, 1)
  expand <- expand.grid(x = x_tot, y = y_tot)
  grid <- as.data.frame(cbind(x = rep(expand$x,100),y = rep(expand$y,100),year = rep(1:100,each=400)))
  years <- seq(1,100,1) #this includes both observed (20 years) and forecast years (80 years)
  
  temp_max_slope <- (temp_diff[4] - temp_diff[2])/100  # linear slope over 100 years
  temp_min_slope <- (temp_diff[3] - temp_diff[1])/100
  temp_max_int <-  temp_diff[2] - temp_max_slope
  temp_min_int <- temp_diff[1] - temp_min_slope
  
  #----Loop through each year----
  for (y in years){
    print(paste0("Creating environmental simulation for Year ",y))
    
    #----Generate Temperature Covariate----
    temp_plain <- raster(ncol=20,nrow=20)
    ex <- extent(0.5,20.5,0.5,20.5)
    extent(temp_plain) <- ex
    xy <- coordinates(temp_plain)
    min <- temp_min_slope*y + temp_min_int 
    max <- temp_max_slope*y + temp_max_int
    
    #----Decide on spatial distribution of temperature----
    if (temp_spatial=="simple"){
      temp_plain[] <- seq(min,max,length.out=400) # assign values to RasterLayer
      temp <- raster::calc (temp_plain, fun = function(x) jitter(x,amount=1))
    }
    
    if (temp_spatial=="matern"){
      # EW: simulate from matern spatial field. We could extend this by
      # (1) letting matern parameters vary over time, (2) letting latitude 
      # gradient vary over time, and (3) making the field simulation be non-independent
      # across years. Right now this generates an independent field by year
      sim_field = glmmfields::sim_glmmfields(n_knots = 40,
                                             n_draws=1, covariance="matern",
                                             g = data.frame(lon=xy[,1],lat=xy[,2]), 
                                             n_data_points=nrow(xy),
                                             B = c(0,-0.05), #making second parameter negative to invert data, and slightly higher (from 0.1) to have better latitudinal siganl 
                                             X = cbind(1,xy[,2]))
      sim_field$dat$new_y = (sim_field$dat$y + abs(min(sim_field$dat$y)))
      # make these compatible with previous range
      temp_plain[] = min + (max-min)*sim_field$dat$new_y / max(sim_field$dat$new_y)
      temp <- temp_plain
    }
    # plot(temp)
    
    #----Use Virtual Species to assign response curve----
    envir_stack <- stack(temp) #must be in raster stack format
    names(envir_stack) <- c('temp')
    
    #----assign response curve with mean at 4 degrees----
    parameters <- virtualspecies::formatFunctions(temp = c(fun="dnorm",mean=4,sd=1))
    
    #----convert temperature raster to species suitability----
    envirosuitability <- virtualspecies::generateSpFromFun(envir_stack,parameters=parameters, rescale = FALSE,rescale.each.response = FALSE)
    #rescale
    ref_max <- dnorm(parameters$temp$args[1], mean=parameters$temp$args[1], sd=parameters$temp$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
    envirosuitability$suitab.raster <- (1/ref_max)*envirosuitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum temp is encountered
    
    #Plot suitability and response curve
    # plot(envirosuitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(envirosuitability) #plot response curves

    #----Convert suitability to Presence-Absence----
    if (PA_shape == "logistic") {
      #SB: specifies alpha and beta of logistic - creates almost knife-edge absence -> presence
      suitability_PA <- virtualspecies::convertToPA(envirosuitability, PA.method = "probability", beta = 0.5,
                                                    alpha = -0.05, species.prevalence = NULL, plot = FALSE)
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
    
    #----Extract suitability for each location----
    print("Extracting suitability")
    for (i in 1:400){
      start_index <- min(which(grid$year==y))
      ii = (i + start_index) -1
      s <- raster::extract(envirosuitability$suitab.raster,grid[ii,1:2]) 
      pa <- raster::extract(suitability_PA$pa.raster,grid[ii,1:2])
      t <- raster::extract(temp,grid[ii,1:2]) 
      grid$suitability[ii] <-  s
      grid$pres[ii] <-  pa
      grid$temp[ii] <-  t
    }
  }
  # summary(grid)
  
  #----Create abundance as a function of the environment----
  if (abund_enviro == "lnorm_low") {
    # SB: values in Ecography paper. I think initially they were based on flounder in EBS but not sure if I edited them
    grid$abundance <- ifelse(grid$pres==1,rlnorm(nrow(grid),2,0.1)*grid$suitability,0)
  }
  if (abund_enviro == "lnorm_high") {
    # EW: I'm cranking up the rlnorm parameters to make it more comparable to wc trawl survey estimates -- these new ones based on arrowtooth
    # SB: rlnorm parameters (6,1) were too large for estimation model. GAMs had explained deviance <10%. Other suggestions?
    grid$abundance <- ifelse(grid$pres==1,rlnorm(nrow(grid),6,1)*grid$suitability,0)
  }
  if (abund_enviro == "poisson") {
    # JS: sample from a Poisson distbn, with suitability proportional to mean (slower, bc it re-samples distbn for each observation)
    maxN <- 20  #max mean abundance at highest suitability
    grid$abundance <- ifelse(grid$pres==1,rpois(nrow(grid),lambda=grid$suitability*maxN),0)
  } 
 
  grid$year <- ifelse(grid$year<=20,grid$year + 2000, grid$year + 2000) #give years meaning: observed years 2000-2020, forecast years 2021-2100.
  
  return(grid)
}

#----EXAMPLE----
# data <- SimulateWorld()
# #plot time-series of total catch in observed years
# plot(aggregate(abundance~year,data[data$year<=2010,],FUN="sum"),type="l",ylab="Abundance")
# #plot time-series of total catch in forecast years
# plot(aggregate(abundance~year,data[data$year>2010,],FUN="sum"),type="l",ylab="Abundance")





