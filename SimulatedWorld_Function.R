#Function to simulate species distribution and abundance with respect to environmental covariates

SimulateWorld <- function(){
  library(virtualspecies)
  library(scales)
  library(raster)
  
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
    temp_plain[] <- xy[,2]
    min <- 0.0404*y + 1.9596 #temp increases by 4 degrees over 100 years
    max <- 0.0404*y + 5.9596
    vals <- seq(min,max,0.01)
    vals <- vals[1:400]
    temp_plain[] <- vals
    temp <- raster::calc (temp_plain, fun = function(x) jitter(x,amount=1))
    # plot(temp)
    
    #----Use Virtual Species to assign response curve----
    envir_stack <- stack(temp) #must be in raster stack format for below 'virtualspecies' function
    names(envir_stack) <- c('temp')
    
    #----assign response curve with mean at 3 degrees----
    parameters <- virtualspecies::formatFunctions(temp = c(fun="dnorm",mean=3,sd=2))
    
    #----convert temperature raster to species suitability----
    envirosuitability <- virtualspecies::generateSpFromFun(envir_stack,parameters=parameters, rescale = FALSE)
    
    #Plot suitability and response curve
    # plot(envirosuitability$suitab.raster) #plot habitat suitability
    # virtualspecies::plotResponse(envirosuitability) #plot response curves
    
    #----Convert suitability to Presence-Absence----
    suitability_PA <- virtualspecies::convertToPA(envirosuitability$suitab.raster, PA.method = "probability", beta = 0.5,
                                  alpha = -0.05, species.prevalence = NULL, plot = FALSE)
    # plot(suitability_PA$pa.raster")
    
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
  grid$abundance <- ifelse(grid$pres==1,rlnorm(1,2,0.1)*grid$suitability,0)
  grid$year <- ifelse(grid$year<=20,grid$year + 2000, grid$year + 2000) #give years meaning: observed years 2000-2020, forecast years 2021-2100.
  
  return(grid)
}

#----EXAMPLE----
# data <- SimulateWorld()
# #plot time-series of total catch in observed years
# plot(aggregate(abundance~year,data[data$year<=2010,],FUN="sum"),type="l",ylab="Abundance")
# #plot time-series of total catch in forecast years
# plot(aggregate(abundance~year,data[data$year>2010,],FUN="sum"),type="l",ylab="Abundance")





