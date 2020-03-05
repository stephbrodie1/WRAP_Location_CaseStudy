#Function to simulate species distribution and abundance with respect to environmental covariates

SimulateWorld <- function(){
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
  
  #----Loop through each year----
  for (y in years){
    print(paste0("Creating environmental simulation for Year ",y))
    
    #----Generate Temperature Covariate----
    temp_plain <- raster(ncol=20,nrow=20)
    ex <- extent(0.5,20.5,0.5,20.5)
    extent(temp_plain) <- ex
    xy <- coordinates(temp_plain)
    #temp_plain[] <- xy[,2] # not used
    # EW: if this is a proxy for bottom temp is 4 deg inc over 100 years too high? 
    # EW: also these values created a surface with a 4 deg change within a year -- 
    # -- basically equaivalent to 400-600m depth in the survey
    min <- 0.0404*y + 1.9596 #temp increases by 4 degrees over 100 years
    max <- 0.0404*y + 5.9596
    # EW: next lines commented out, replaced below
    # temp_plain[] <- seq(min,max,length.out=400) # assign values to RasterLayer
    # temp <- raster::calc (temp_plain, fun = function(x) jitter(x,amount=1))
    
    # EW: simulate from matern spatial field. We could extend this by
    # (1) letting matern parameters vary over time, (2) letting latitude 
    # gradient vary over time, and (3) making the field simulation be non-independent
    # across years. Right now this generates an independent field by year
    sim_field = glmmfields::sim_glmmfields(n_knots = 40,
      n_draws=1, covariance="matern",
      g = data.frame(lon=xy[,1],lat=xy[,2]), 
      n_data_points=nrow(xy),
      B = c(0,-0.05), #making second parameter negative to invert data, and slightly higher to have better latitudinal siganl 
      X = cbind(1,xy[,2]))
    sim_field$dat$new_y = (sim_field$dat$y + abs(min(sim_field$dat$y)))
    # make these compatible with previous range
    temp_plain[] = min + (max-min)*sim_field$dat$new_y / max(sim_field$dat$new_y)
    temp <- temp_plain
    # plot(temp)
    
    #----Use Virtual Species to assign response curve----
    envir_stack <- stack(temp) #must be in raster stack format for below 'virtualspecies' function
    names(envir_stack) <- c('temp')
    
    #----assign response curve with mean at 3 degrees----
    parameters <- virtualspecies::formatFunctions(temp = c(fun="dnorm",mean=4,sd=1))
    
    #----convert temperature raster to species suitability----
    envirosuitability <- virtualspecies::generateSpFromFun(envir_stack,parameters=parameters, rescale = FALSE)
    
    #Plot suitability and response curve
    # plot(envirosuitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(envirosuitability) #plot response curves
    
    #----Convert suitability to Presence-Absence----
    suitability_PA <- virtualspecies::convertToPA(envirosuitability$suitab.raster, PA.method = "probability", beta = 0.5,
      alpha = -0.05, species.prevalence = NULL, plot = FALSE)
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
  summary(grid)
  
  #----Create abundance as a function of the environment----
  # EW: I'm cranking up the rlnorm parameters to make it more comparable to wc trawl survey estimates -- these new ones based on arrowtooth
  # SB: rlnorm parameters (6,1) were too large for estimation model. GAMs had explained deviance <10%. 
  # grid$abundance <- ifelse(grid$pres==1,rlnorm(400,6,1)*grid$suitability,0)
  grid$abundance <- ifelse(grid$pres==1,rlnorm(400,2,0.1)*grid$suitability,0)
  grid$year <- ifelse(grid$year<=20,grid$year + 2000, grid$year + 2000) #give years meaning: observed years 2000-2020, forecast years 2021-2100.
  
  return(grid)
}

#----EXAMPLE----
# data <- SimulateWorld()
# #plot time-series of total catch in observed years
# plot(aggregate(abundance~year,data[data$year<=2010,],FUN="sum"),type="l",ylab="Abundance")
# #plot time-series of total catch in forecast years
# plot(aggregate(abundance~year,data[data$year>2010,],FUN="sum"),type="l",ylab="Abundance")





