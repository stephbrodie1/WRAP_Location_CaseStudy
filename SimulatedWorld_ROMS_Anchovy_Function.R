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
  # library(glmmfields)
  library(dplyr)
  
  #----Create output file----
  #This will be the information passed to the estimation model
  gridcells <- 4012 #number of grid cells in ROMS (NAs removed)
  output <- as.data.frame(matrix(NA, nrow=(gridcells*121), ncol=10))  #21912 non-NA grid cells in ROMS
  colnames(output) <- c("lon","lat","year","pres","suitability","sst","chla_surface", "z","zoo_50", "sampled")
  
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
  pop$AnchovySSB_norm <- BBmisc::normalize(pop$AnchovySSB,method = "range",range=c(0.8,1)) #use this range to vary prevalence according to popn. dynamics
  # plot(pop$Year,pop$AnchovySSB, type='l')
  # abline(h=1172713)
  mean(pop$AnchovySSB)
  quantile(pop$AnchovySSB)
  lowYears <- which(pop$AnchovySSB<1172713)
  
  #----Load in rasters----
  #These are the average spring conditions from the downscaled Hadley earth system model
  files_sst <- list.files(paste0(dir,'sst_monthly'), full.names = TRUE, pattern=".grd") 
  files_chl <- list.files(paste0(dir,'chl_surface'), full.names = TRUE, pattern=".grd") 
  # files_zoo200 <- list.files(paste0(dir,'zoo_200m'), full.names = TRUE, pattern=".grd")
  files_zoo50 <- list.files(paste0(dir,'zoo_50m'), full.names = TRUE, pattern=".grd") 
  years <- seq(1980,2100,1)
  
  #----Loop through each year----
  for (y in 1:121){
      print(paste0("Creating environmental simulation for Year ",years[y]))
      
      #Load in environmental rasters for a specific year
      sst <- raster(files_sst[y])
      chla <- raster(files_chl[y])
      chla <- log(chla)
      zoo <- raster(files_zoo50[y])
      z <- raster('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/gfdl/bottom_layer_depth.grd')
      z <- z * -1
      lat <- z
      xy <- coordinates(lat)
      lat[] <- xy[, 2]
      lat[] <- BBmisc::normalize(lat@data@values,method = "range",range=c(1,0)) 
      
      #----create a mask of a smaller domain-----
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
      
      sst <- mask(sst,rt)
      chla <- mask(chla,rt)
      zoo <- mask(zoo,rt)
      z <- mask(z,rt)
      lat <- mask(lat,rt)
      # plot(z)
      # contour(z,add=T,levels=c(-1000))
      # # z[] <- BBmisc::normalize(z@data@values, method = "range", range = c(0, 1))
      
      # if(y %in% lowYears){
      #   #Low biomass years: force only southern habitats to be suitable
      #   
      #   #Stack rasters
      #   spA_stack <- stack(sst, z, lat)
      #   names(spA_stack) <- c('sst', 'z','lat')
      #   
      #   #Assign preferences
      #   spA_parameters <- formatFunctions(sst = c(fun="dnorm",mean=15,sd=7),
      #                                     # zoo = c(fun="logisticFun",alpha=-5,beta=20),
      #                                     z = c(fun="logisticFun",alpha=-500,beta=-2000),
      #                                     lat = c(fun="linearFun",a=1,b=0))
      #   spA_suitability <- generateSpFromFun(spA_stack,parameters=spA_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
      #   # plot(spA_suitability$suitab.raster) #plot habitat suitability
      #   # virtualspecies::plotResponse(spA_suitability) #plot response curves
      #   
      #   #manually rescale
      #   ref_max_sst <- dnorm(spA_parameters$sst$args[1], mean=spA_parameters$sst$args[1], sd=spA_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
      #   ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spA_parameters$zoo$args[2])/spA_parameters$zoo$args[1]))
      #   ref_max_z <- 1 / (1 + exp(((z@data@max) - spA_parameters$z$args[2])/spA_parameters$z$args[1]))
      #   ref_max_lat <- 1
      #   ref_max <- ref_max_sst  * ref_max_z  * ref_max_lat
      #   spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
      #   print(spA_suitability$suitab.raster)
      #   plot(spA_suitability$suitab.raster) #plot habitat suitability
      #   
      # } else{
        
      #----assign response curves----
      #Stack rasters
      spA_stack <- stack(sst, zoo, z)
      names(spA_stack) <- c('sst', 'zoo', 'z')
      
      #Assign preferences
      spA_parameters <- formatFunctions(sst = c(fun="dnorm",mean=16,sd=6),
                                        zoo = c(fun="logisticFun",alpha=-5,beta=20),
                                        z = c(fun="logisticFun",alpha=-500,beta=-2000))
      spA_suitability <- generateSpFromFun(spA_stack,parameters=spA_parameters, rescale = FALSE,rescale.each.response = FALSE) #Important: make sure rescaling is false. Doesn't work well in the 'for' loop.
      # plot(spA_suitability$suitab.raster) #plot habitat suitability
      # virtualspecies::plotResponse(spA_suitability) #plot response curves
      
      #manually rescale
      ref_max_sst <- dnorm(spA_parameters$sst$args[1], mean=spA_parameters$sst$args[1], sd=spA_parameters$sst$args[2]) #JS/BM: potential maximum suitability based on optimum temperature
      ref_max_zoo <- 1 / (1 + exp(((zoo@data@max) - spA_parameters$zoo$args[2])/spA_parameters$zoo$args[1]))
      ref_max_z <- 1 / (1 + exp(((z@data@max) - spA_parameters$z$args[2])/spA_parameters$z$args[1]))
      ref_max <- ref_max_sst * ref_max_zoo * ref_max_z
      spA_suitability$suitab.raster <- (1/ref_max)*spA_suitability$suitab.raster #JS/BM: rescaling suitability, so the max suitbaility is only when optimum is encountered
      # print(spA_suitability$suitab.raster)
      # plot(spA_suitability$suitab.raster) #plot habitat suitability
      
    # } #end else loop here
      
      #----Convert suitability to Presence-Absence----
      #Use a specific function to convert suitability (0-1) to presence or absence (1 or 0)
      # prev <- pop$AnchovySSB[pop$year==y]
      # prev <- (pop$AnchovySSB - min(pop$AnchovySSB)) / (max(pop$AnchovySSB) - min(pop$AnchovySSB))
      set.seed(y)
      suitability_PA <- virtualspecies::convertToPA(spA_suitability, PA.method = "probability", beta = 0.4,
                                                    alpha = -0.07, species.prevalence = NULL, plot = FALSE)
      # plotSuitabilityToProba(suitability_PA) #Let's you plot the shape of conversion function
      # plot(suitability_PA$pa.raster)
      
      #-----Sample Presences and generate Biomass-----
      # set.seed(y)
      # find_prevalence <- sampleOccurrences(suitability_PA,n = 500,type = "presence-absence",
      #                                      detection.probability = 1,error.probability=0, plot = FALSE)
      # suit_prev <- as.numeric(find_prevalence$sample.prevalence[1])
      # pop_prev <- pop$AnchovySSB_norm[pop$year==y]
      # sp_prev <- pop_prev * suit_prev
      
      set.seed(y)
      presence.points <- sampleOccurrences(suitability_PA,n = 500,type = "presence-absence",
                                           detection.probability = 1,error.probability=0, plot = FALSE)
      # sampling.area = test)
      # sample.prevalence = round(sp_prev,1))
      # sample.prevalence = round(prev[y],1))
      # sample.prevalence = round(sp_prev,1))
      # sample.prevalence = round(prev[y],1))
      # sum(presence.points$sample.points$Real)
      df <- cbind(as.data.frame(presence.points$sample.points$x),as.data.frame(presence.points$sample.points$y))
      colnames(df) <- c("lon","lat")
      df$sampled <- 1
    
    
    #----Extract data for each year----
    print("Extracting suitability")
    ei <- gridcells*y #end location in output grid to index to
    se <- ei - (gridcells-1) #start location in output grid to index to
    output$lat[se:ei] <- rasterToPoints(sst)[,2]
    output$lon[se:ei] <- rasterToPoints(sst)[,1]
    output$year[se:ei] <- rep(years[y],gridcells)
    output$pres[se:ei] <- rasterToPoints(suitability_PA$pa.raster)[,3] 
    output$suitability[se:ei] <- rasterToPoints(suitability_PA$suitab.raster)[,3] 
    output$sst[se:ei] <-  rasterToPoints(sst)[,3]
    output$chla_surface[se:ei] <-  rasterToPoints(chla)[,3]
    output$zoo_50[se:ei] <- rasterToPoints(zoo)[,3]
    output$z[se:ei] <-  rasterToPoints(z)[,3]
    #temp file
    temp_sampled <- left_join(output[se:ei,1:9], df, by=c('lon','lat'))
    #write to output
    output$sampled[se:ei] <- temp_sampled$sampled
    output$pop_norm[se:ei] <-  pop$AnchovySSB_norm[pop$year==y]
    output$low_biomass_year[se:ei] <- ifelse(y %in% lowYears, rep(1,4012), rep(0,4012))
  }
  
  a_mean <- log((mean(pop$AnchovySSB) / gridcells )*2) #divide by number of grid cells, but multiply by 1.7 to offset multiplying by suitablity & pop below
  a_sd <- log((sd(pop$AnchovySSB) / gridcells)/220 ) 
  # hist(rlnorm(gridcells,a_mean,a_sd))
  # sum(rlnorm(gridcells,a_mean,a_sd))
  # sum(rlnorm(nsamples,a_mean,a_sd))
  
  #add in something where abundance in low suitability years is weighted towards the south
  mean_suit_annual <- aggregate(suitability~year,data=output[output$year<=2010 & output$year>=1985,],FUN=mean)[2]
  mean_suit_annual_allyears <- aggregate(suitability~year,data=output,FUN=mean)
  suit_quantil <- round(quantile(mean_suit_annual$suitability, 0.25),2)
  low_suit <- which(mean_suit_annual_allyears$suitability<=suit_quantil)
  low_years <- mean_suit_annual_allyears$year[low_suit]
  output$low_years_binary <- ifelse(output$year %in% low_years, 1, 0)
  
  #NOW FOR LOW YEARS, weight abundance towards the south
  output$lat_norm <- BBmisc::normalize(output$lat,method = "range",range=c(1,0.8))
  output$suitability_inclowyears <- ifelse(output$low_years_binary==1,output$suitability*output$lat_norm,output$suitability)
  output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),a_mean, a_sd)*output$suitability_inclowyears*output$pop_norm,0)
  
  # output$abundance <- ifelse(output$pres==1,rlnorm(nrow(output),a_mean, a_sd)*output$suitability*output$pop_norm,0)

  
  # sum(output$abundance, na.rm=T)
  #Convert NAs in 'sampled' columns to zeros
  output$sampled <- ifelse(is.na(output$sampled),0, output$sampled)
  
  return(output)
}
