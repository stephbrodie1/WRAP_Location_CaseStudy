##---- Quantifying novel space in the environmental niche   ----##
##----  Version for Steph                                    ----##
##----  J.A.Smith - UCSC/NOAA - June 2021 (R v4.0.4)         ----##
##---- since modified by steph----

## Notes on hypervolume package:
## - Best to center and scale all variables; this requires a global mean and SD...
## ... I calculate these from my historical reference set (for each ESM) and use them for all future data.
## - For trajectories of future novelty, I use an 'inclusion test', which requires...
## ... a data.frame to build the historical hypervolume (in an SDM context, your observations), and...
## ... rasters for your future period (in an SDM context, your prediction data).
## - Hypervolumes are computationally intense, and more than 6 climate variables and ~300K observations...
## ... is asking for trouble re: fitting. For me, this meant subsetting 5% of my 30y historical...
## ... data at the native ROMS resolution (using monthly means).
## - There are a couple of algorithms to delineate the hypervolume, but I prefer support vector machine...
## ... mainly because it fits around the around the data very closely, and indicates...
## ... whether future data are 'inside or outside' historical (rather than a probability).

## *** MAKE SURE YOU DO ONE ESM/MODEL AT A TIME - the means and SDs are specific to each historical data set...
## ... and these will change if the historical data set changes


#NOTE: RECENT UPDATES USE %NEARBY FROM EXDET package

#----library-----
library(hypervolume)
library(raster)
library(glue)

##----List enviro variables which define model niche (i.e. SDM variables used in fitting)----
env_vars_hms <- c("temp", "mld")
env_vars_cps <- c("temp", "chla","z")
env_vars_gfs <- c("btemp", "O2","z")

#-----Create historical hypervolume (1985-2010)-----
#Get historical data. #Use data SDMs were trained on
#Function to scale historical data
# esm can be one of c("had","gfdl","ipsl")
# species can be one of c("Albacore EMs","Anchovy EMs","Groundfish EMs")

# make_hist_data <- function(esm,species,env_vars){
#   #load data
#   hist_data <- readRDS(paste0("~/PROJECTS/WRAP Location/",species,"/",esm,"/dat_hist_results_full.rds"))
#   #rescale
#   means <- colMeans(as.data.frame(hist_data)[env_vars])  #'global' values used to standardise all your data
#   SDs <- apply(hist_data, 2, sd)[env_vars]
#   hist_data_s <- scale(hist_data[,env_vars], center=means, scale=SDs)  #center and scale
#   hist_data_s <- as.data.frame(hist_data_s)
# }


#------loop through all years----
t1 <- Sys.time()
extrap_output <- list()
output <- as.data.frame(matrix(NA, nrow=90,ncol=5))
colnames(output) <- c("year","excluded","nearby_median","nearby_mean","extrap_percent")
counter=1
for (s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  for (e in c("had","gfdl","ipsl")){
    print(glue("Running historical hypervolume for {e} {s}"))
    #Get historical data
    if(s=="Albacore EMs"){env_vars = env_vars_hms}
    if(s=="Anchovy EMs"){env_vars = env_vars_cps}
    if(s=="Groundfish EMs"){env_vars = env_vars_gfs}
    hist_data <- readRDS(glue("~/Dropbox/PROJECTS/WRAP Location/{s}/{e}/dat_hist_results_full.rds"))
    means <- colMeans(as.data.frame(hist_data)[env_vars])  #'global' values used to standardise all your data
    SDs <- apply(hist_data, 2, sd)[env_vars]
    hist_data_s <- scale(hist_data[,env_vars], center=means, scale=SDs)  #center and scale
    hist_data_s <- as.data.frame(hist_data_s)
    hvh = hypervolume(data=hist_data_s, method='svm')
    # plot(hvh, show.3d=F)
    counter_year <- 1
    for (y in 2011:2100){
      print(glue("Running projection hypervolume for {y}"))
      
      #load in data
      fcast <- readRDS(glue('~/Dropbox/PROJECTS/WRAP Location/{s}/{e}/dat_fcast_results_full.rds'))
      #get correct year
      yi <- fcast[fcast$year==y,]

      if(s=="Albacore EMs"){
        #Create rasters
        temp <- rasterFromXYZ(yi[,c(1,2,6)])
        mld <- rasterFromXYZ(yi[,c(1,2,9)])
        rx <- stack(temp, mld)
        rx$temp <- scale(rx$temp, center=means["temp"], scale=SDs["temp"])  #center and scale using global values
        rx$mld <- scale(rx$mld, center=means["mld"], scale=SDs["mld"])
        #create hypervolume novelty
        hyp_proj <- hypervolume_project(hvh, rasters=rx,type="inclusion", fast.or.accurate="accurate")
        excl <- length(hyp_proj[hyp_proj==0])  #number of cells 'excluded'; in example ~3000 cells
        excl_prop <- length(hyp_proj[hyp_proj==0])/length(hyp_proj[hyp_proj>=0])  #proportion of cells excluded; in example ~14%
        # plot(hyp_proj)  #this is the inclusion raster; green = included (analog conditions), white = excluded (novel conditions)

        temp <- rasterFromXYZ(yi[,c(1,2,6)])
        mld <- rasterFromXYZ(yi[,c(1,2,9)])
        rx <- stack(temp, mld)
        t <- as.data.frame(rasterToPoints(rx))
        aftt_crs <- sp::CRS("+proj=longlat +datum=WGS84")
        nearby <- compute_nearby(samples = hist_data[,env_vars],
                       prediction.grid = t,
                       coordinate.system = aftt_crs,
                       covariate.names = env_vars,
                       nearby = 1)
        extrap <- compute_extrapolation(samples = hist_data[,env_vars],
                              covariate.names = env_vars,
                              prediction.grid = t,
                              coordinate.system = aftt_crs)
        extrap_percent <- 100 - extrap$summary$extrapolation$analogue.p
        nearby_median  <- median(nearby$raster@data@values, na.rm=T) #11 for y2100; 22 for 2011; 20 for 2030
        nearby_mean  <-   mean(nearby$raster@data@values, na.rm=T) #5989 for y2100; 11117 for 2011; 10039.79 for 2030
       
      }
      if(s=="Anchovy EMs"){
        #Create rasters
        temp <- rasterFromXYZ(yi[,c(1,2,6)])
        chla <- rasterFromXYZ(yi[,c(1,2,7)])
        z <- rasterFromXYZ(yi[,c(1,2,8)])
        rx <- stack(temp, chla,z)
        rx$temp <- scale(rx$temp, center=means["temp"], scale=SDs["temp"])  #center and scale using global values
        rx$chla <- scale(rx$chla, center=means["chla"], scale=SDs["chla"])
        rx$z <- scale(rx$z, center=means["z"], scale=SDs["z"])
        #create hypervolume novelty
        hyp_proj <- hypervolume_project(hvh, rasters=rx,type="inclusion", fast.or.accurate="accurate")
        excl <- length(hyp_proj[hyp_proj==0])  #number of cells 'excluded'; in example ~3000 cells
        excl_prop <- length(hyp_proj[hyp_proj==0])/length(hyp_proj[hyp_proj>=0])  #proportion of cells excluded; in example ~14%
        # plot(hyp_proj, asp=1)  #this is the inclusion raster; green = included (analog conditions), white = excluded (novel conditions)
        
        temp <- rasterFromXYZ(yi[,c(1,2,6)])
        chla <- rasterFromXYZ(yi[,c(1,2,7)])
        z <- rasterFromXYZ(yi[,c(1,2,8)])
        rx <- stack(temp, chla,z)
        t <- as.data.frame(rasterToPoints(rx))
        aftt_crs <- sp::CRS("+proj=longlat +datum=WGS84")
        nearby <- compute_nearby(samples = hist_data[,env_vars],
                                 prediction.grid = t,
                                 coordinate.system = aftt_crs,
                                 covariate.names = env_vars,
                                 nearby = 1)
        extrap <- compute_extrapolation(samples = hist_data[,env_vars],
                                        covariate.names = env_vars,
                                        prediction.grid = t,
                                        coordinate.system = aftt_crs)
        extrap_percent <- 100 - extrap$summary$extrapolation$analogue.p
        nearby_median   <- median(nearby$raster@data@values, na.rm=T) #11 for y2100; 22 for 2011; 20 for 2030
        nearby_mean  <-   mean(nearby$raster@data@values, na.rm=T) #5989 for y2100; 11117 for 2011; 10039.79 for 2030
        
      }
      if(s=="Groundfish EMs"){
        #Create rasters
        btemp <- rasterFromXYZ(yi[,c(1,2,6)])
        O2 <- rasterFromXYZ(yi[,c(1,2,7)])
        z <- rasterFromXYZ(yi[,c(1,2,8)])
        rx <- stack(btemp, O2,z)
        rx$btemp <- scale(rx$btemp, center=means["btemp"], scale=SDs["btemp"])  #center and scale using global values
        rx$O2 <- scale(rx$O2, center=means["O2"], scale=SDs["O2"])
        rx$z <- scale(rx$z, center=means["z"], scale=SDs["z"])
        #create hypervolume novelty
        hyp_proj <- hypervolume_project(hvh, rasters=rx,type="inclusion", fast.or.accurate="accurate")
        excl <- length(hyp_proj[hyp_proj==0])  #number of cells 'excluded'; in example ~3000 cells
        excl_prop <- length(hyp_proj[hyp_proj==0])/length(hyp_proj[hyp_proj>=0])  #proportion of cells excluded; in example ~14%
        # plot(hyp_proj, asp=1)  #this is the inclusion raster; green = included (analog conditions), white = excluded (novel conditions)
        
        btemp <- rasterFromXYZ(yi[,c(1,2,6)])
        O2 <- rasterFromXYZ(yi[,c(1,2,7)])
        z <- rasterFromXYZ(yi[,c(1,2,8)])
        rx <- stack(btemp, O2,z)
        t <- as.data.frame(rasterToPoints(rx))
        aftt_crs <- sp::CRS("+proj=longlat +datum=WGS84")
        nearby <- compute_nearby(samples = hist_data[,env_vars],
                                 prediction.grid = t,
                                 coordinate.system = aftt_crs,
                                 covariate.names = env_vars,
                                 nearby = 1)
        extrap <- compute_extrapolation(samples = hist_data[,env_vars],
                                        covariate.names = env_vars,
                                        prediction.grid = t,
                                        coordinate.system = aftt_crs)
        extrap_percent <- 100 - extrap$summary$extrapolation$analogue.p
        nearby_median   <- median(nearby$raster@data@values, na.rm=T) #11 for y2100; 22 for 2011; 20 for 2030
        nearby_mean  <-   mean(nearby$raster@data@values, na.rm=T) #5989 for y2100; 11117 for 2011; 10039.79 for 2030
        
      }
      output[counter_year,1] <- y
      output[counter_year,2] <- excl_prop
      output[counter_year,3] <- nearby_median
      output[counter_year,4] <- nearby_mean
      output[counter_year,5] <- extrap_percent
      counter_year <- counter_year+1
    }
    extrap_output[[counter]] <- output
    counter=counter+1
  }
}
t2 <- Sys.time()
t2-t1

#Save
saveRDS(extrap_output,'~/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet_PlusNearby.rds')

#Quick plots
# d1 <- extrap_output[[1]]
# ggplot(data=d1,aes(x=year,y=nearby_sum))+
#   geom_line()
# ggplot(data=d1,aes(x=year,y=excluded))+
#   geom_line()



#------MAKE PLOTS OF 2100-----

#enviro vars
env_vars_hms <- c("temp", "mld")
env_vars_cps <- c("temp", "chla","z")
env_vars_gfs <- c("btemp", "O2","z")

#raster mask for anchovy and groundfish
scb_coords=matrix(c(-126,48,
                    -115,48,
                    -115,30,
                    -119,30,
                    -126,40),ncol=2,byrow = T)
p=Polygon(scb_coords)
ps=Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
proj4string(sps) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
data = data.frame(f=99.9)
test = SpatialPolygonsDataFrame(sps,data)
r <- raster(ncol=185,nrow=180)
extent(r) <- extent(z)
rt <- rasterize(test,r)


for (s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  for (e in c("had","gfdl","ipsl")){
    if(s=="Albacore EMs"){env_vars = env_vars_hms}
    if(s=="Anchovy EMs"){env_vars = env_vars_cps}
    if(s=="Groundfish EMs"){env_vars = env_vars_gfs}
    
    hist_data <- readRDS(glue("~/PROJECTS/WRAP Location/{s}/{e}/dat_hist_results_full.rds"))
    means <- colMeans(as.data.frame(hist_data)[env_vars])  #'global' values used to standardise all your data
    SDs <- apply(hist_data, 2, sd)[env_vars]
    hist_data_s <- scale(hist_data[,env_vars], center=means, scale=SDs)  #center and scale
    hist_data_s <- as.data.frame(hist_data_s)
    hvh = hypervolume(data=hist_data_s, method='svm')

    for (y in 2100){
      print(glue("Running projection hypervolume for {y}"))
      
      #load in data
      fcast <- readRDS(glue('~/PROJECTS/WRAP Location/{s}/{e}/dat_fcast_results_full.rds'))
      #get correct year
      yi <- fcast[fcast$year==y,]
      
      temp <- raster(glue('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/{e}/sst_monthly/sst_monthly_{e}_SpringMean_2100.grd'))
      mld <- raster(glue('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/{e}/ild_0.5C/ild_0.5C_{e}_SpringMean_2100.grd'))
      chla <- raster(glue('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/{e}/chl_surface/chl_surface_{e}_SpringMean_2100.grd'))
      btemp<- raster(glue('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/{e}/temp_bottom/temp_bottom_{e}_SpringMean_2100.grd'))
      O2<- raster(glue('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/{e}/oxygen_bottom/oxygen_bottom_{e}_SpringMean_2100.grd'))
      z<- raster(glue('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/{e}/bottom_layer_depth.grd'))
      
      if(s=="Albacore EMs"){
        #Create rasters
        rx <- stack(temp, mld)
        names(rx) <- c("temp","mld")
        rx$temp <- scale(rx$temp, center=means["temp"], scale=SDs["temp"])  #center and scale using global values
        rx$mld <- scale(rx$mld, center=means["mld"], scale=SDs["mld"])
        #create hypervolume novelty
        hyp_proj <- hypervolume_project(hvh, rasters=rx,type="inclusion", fast.or.accurate="accurate")
        # excl <- length(hyp_proj[hyp_proj==0])  #number of cells 'excluded'; in example ~3000 cells
        # excl_prop <- length(hyp_proj[hyp_proj==0])/length(hyp_proj[hyp_proj>=0])  #proportion of cells excluded; in example ~14%
        # plot(hyp_proj, asp=1)  #this is the inclusion raster; green = included (analog conditions), white = excluded (novel conditions)
        p1 <- as.data.frame(rasterToPoints(hyp_proj))
        g1 <- ggplot(data=p1, aes(x=x, y=y))+
          geom_tile(aes(fill=temp))+
          theme_classic() +  labs(y="", x="", title = glue('Novel Habitat in 2100: {s} {e}'))+
          theme(legend.position="none",legend.title = element_blank(),
                plot.title = element_text(size = 10))+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+#makes a box
          scale_fill_gradientn(colours = c("#e74c3c","#3498db")) + 
          annotation_map(map_data("world"), colour = "black", fill="grey50")+
          coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
          scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 
        tiff(glue('~/PROJECTS/WRAP Location/Manuscript/Figures/Hypervolume/Extrapolation_{e}_{s}_2100.tiff'), res=300, units="in",width=8,height=10)
        plot(g1)
        dev.off()
      }
      if(s=="Anchovy EMs"){
        #Create rasters
        rx <- stack(temp, chla,z)
        names(rx) <- c("temp","chla","z")
        rx <- mask(rx,rt)
        rx$temp <- scale(rx$temp, center=means["temp"], scale=SDs["temp"])  #center and scale using global values
        rx$chla <- scale(rx$chla, center=means["chla"], scale=SDs["chla"])
        rx$z <- scale(rx$z, center=means["z"], scale=SDs["z"])
        #create hypervolume novelty
        hyp_proj <- hypervolume_project(hvh, rasters=rx,type="inclusion", fast.or.accurate="accurate")
        # excl <- length(hyp_proj[hyp_proj==0])  #number of cells 'excluded'; in example ~3000 cells
        # excl_prop <- length(hyp_proj[hyp_proj==0])/length(hyp_proj[hyp_proj>=0])  #proportion of cells excluded; in example ~14%
        # plot(hyp_proj, asp=1)  #this is the inclusion raster; green = included (analog conditions), white = excluded (novel conditions)
        p1 <- as.data.frame(rasterToPoints(hyp_proj))
        g1 <- ggplot(data=p1, aes(x=x, y=y))+
          geom_tile(aes(fill=temp))+
          theme_classic() +  labs(y="", x="", title = glue('Novel Habitat in 2100: {s} {e}'))+
          theme(legend.position="none",legend.title = element_blank(),
                plot.title = element_text(size = 10))+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+#makes a box
          scale_fill_gradientn(colours = c("#e74c3c","#3498db")) + 
          annotation_map(map_data("world"), colour = "black", fill="grey50")+
          coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
          scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 
        
        tiff(glue('~/PROJECTS/WRAP Location/Manuscript/Figures/Hypervolume/Extrapolation_{e}_{s}_2100.tiff'), res=300, units="in",width=8,height=10)
        plot(g1)
        dev.off()
      }
      if(s=="Groundfish EMs"){
        #Create rasters
        rx <- stack(btemp, O2,z)
        names(rx) <- c("btemp","O2","z")
        rx <- mask(rx,rt)
        rx$btemp <- scale(rx$btemp, center=means["btemp"], scale=SDs["btemp"])  #center and scale using global values
        rx$O2 <- scale(rx$O2, center=means["O2"], scale=SDs["O2"])
        rx$z <- scale(rx$z, center=means["z"], scale=SDs["z"])
        #create hypervolume novelty
        hyp_proj <- hypervolume_project(hvh, rasters=rx,type="inclusion", fast.or.accurate="accurate")
        # excl <- length(hyp_proj[hyp_proj==0])  #number of cells 'excluded'; in example ~3000 cells
        # excl_prop <- length(hyp_proj[hyp_proj==0])/length(hyp_proj[hyp_proj>=0])  #proportion of cells excluded; in example ~14%
        # plot(hyp_proj, asp=1)  #this is the inclusion raster; green = included (analog conditions), white = excluded (novel conditions)
        p1 <- as.data.frame(rasterToPoints(hyp_proj))
        g1 <- ggplot(data=p1, aes(x=x, y=y))+
          geom_tile(aes(fill=btemp))+
          theme_classic() +  labs(y="", x="", title = glue('Novel Habitat in 2100: {s} {e}'))+
          theme(legend.position="none",legend.title = element_blank(),
                plot.title = element_text(size = 10))+
          theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+#makes a box
          scale_fill_gradientn(colours = c("#e74c3c","#3498db")) + 
          annotation_map(map_data("world"), colour = "black", fill="grey50")+
          coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
          scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 
        tiff(glue('~/PROJECTS/WRAP Location/Manuscript/Figures/Hypervolume/Extrapolation_{e}_{s}_2100.tiff'), res=300, units="in",width=8,height=10)
        plot(g1)
        dev.off()
      }
    }
  }
}
