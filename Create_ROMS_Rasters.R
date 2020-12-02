#Convert ROMS GCM netcdf files into rasters
#Files provided by Mike Jacox. Use beyond the WRAP workshop requires permission by Mike Jacox (michael.jacox@noaa.gov)
#Note: 1980-2010 are not observed values. 1980-present should not be compared to observations since the interannual variability will not match up (by design). 

#8 covariates: sst, bottom temp, ild, chl-a surface, chl-a 50m, bottom oxygen, depth of 2.0 oxytherm, bottom depth layer (reference point)
#3 global climate models: Hadley, IPSL, GFDL

#Dimensions:
#1452 years: 1980 - 2100 (121 years * 12 months)
#1452 months: 1 - 12 (121 years * 12 months)
#181 lats: 30 - 48 degrees north @ 0.1 degree resolution
#186 lons: 115.5 - 134 degrees west @ 0.1 degree resolution

#-----Set Wd----
wd <- '~/Dropbox/WRAP Location^3/'
setwd(wd)

#----load library----
library(ncdf4)
library(raster)
library(tidync)

#-----Load in Dropbox files----
#note: download files locally. Netcdfs are contained in a folder called 2d_fields
files_long <- list.files('2d_fields/monthly/', full.names = TRUE)
files_short <- list.files('2d_fields/monthly', full.names = FALSE)

for (f in c(2:22)){ #skipping bottom layer depth as it is a unique (and static) 
  print(files_long[f])
  nc <- nc_open(files_long[f])
  
  # print(nc) #run if you want more details on file
  lat <- ncvar_get(nc, 'lat'); lat <- lat[1,]
  lon <- ncvar_get(nc, 'lon'); lon <- lon[,1]
  year <- ncvar_get(nc, 'year') 
  month <- ncvar_get(nc, 'month')
  name <-  names(nc$var)[5]
  tmp.array <- ncvar_get(nc,name)
  
  #Loop through every time step in ncdf. 
  for (i in 1:1452){
    r <- raster(t(tmp.array[,,i]),
                xmn=min(lon), xmx=max(lon),
                ymn=min(lat), ymx=max(lat), 
                crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    r <- flip(r,2)
    #Forcing resolution our ROMS template (a quick fix to resolve grid edge vs. mid point)
    template = raster('~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd') 
    r <- raster::resample(r, template, method="bilinear")  
    # plot(r)
    
    #create nested folders to save files
    #First create GCM folder
    gcm_options <- c("gfdl", "had", "ipsl")
    gcm_name <- unlist(strsplit(files_short[f],"_"))
    index <- which(gcm_name %in% gcm_options)
    gcm_folder <- gcm_name[index]
    dir.create("Rasters_2d_monthly",showWarnings = FALSE)
    dir.create(paste0('Rasters_2d_monthly/',gcm_folder), showWarnings = FALSE)
    
    var <- unlist(strsplit(files_short[f],"_"))[1:2]
    folder <- paste0('Rasters_2d_monthly/',gcm_folder,"/",paste(var[1],var[2],sep="_"))
    dir.create(folder, showWarnings = FALSE)
    
    month_num <- format(as.Date(paste("1666",month[i],"06",sep="-"),"%Y-%m-%d"), "%m")
    #create name of raster file
    fname <- paste0(var[1],"_",var[2],"_",gcm_folder,"_",year[i],"_month",month_num,".grd")
    
    #save raster
    writeRaster(r, paste0(folder,"/", fname), overwrite=TRUE)
  }
  nc_close(nc)
}


#----Create Spring Average Conditions----
#For each year and variable, load in monthly data from April-June and average. Output to a new folder. 
#Make sure to manually update GCM folder directory (because I'm lazy and didn't code it in)
variables <- c("chl_surface","ild_0.5C","oxygen_bottom","sst_monthly","temp_bottom","zoo_200m","zoo_50m")
for (variable in variables){
  print(variable)
  
  for (i in 1:121){
    files <- list.files(paste0('~/Dropbox/WRAP Location^3/Rasters_2d_monthly/ipsl/',variable), pattern=".grd" , full.names = T)
    month_idx <- rep(1:12,times=121)
    spring_months <- which(month_idx==4,5,6)
    
    start_indx <- spring_months[i]
    
    apr <- raster(files[start_indx])
    may <- raster(files[start_indx+1])
    jun <- raster(files[start_indx+2])
    spring_r <- mean(apr,may,jun)
    years <- seq(1980,2100,1)
    writeRaster(spring_r,paste0('~/Dropbox/WRAP Location^3/Rasters_2d_Spring/ipsl/',variable,'/',variable,'_ipsl_SpringMean_',years[i],'.grd'), overwrite=TRUE)
  }
}

# #-------Quick Data exploration-----
# #have a quick looksie at variable means etc. 
# files_long <- list.files('~/Dropbox/WRAP Location^3/2d_fields/monthly/', full.names = TRUE, )
# files_short <- list.files('~/Dropbox/WRAP Location^3/2d_fields/monthly/', full.names = FALSE)
# f=6
# nc <- nc_open(files_long[f])
# lat <- ncvar_get(nc, 'lat'); lat <- lat[1,]
# lon <- ncvar_get(nc, 'lon'); lon <- lon[,1]
# year <- ncvar_get(nc, 'year')
# month <- ncvar_get(nc, 'month')
# name <-  names(nc$var)[5]
# tmp.array <- ncvar_get(nc,name)
# 
# test <- as.data.frame(apply(tmp.array, 3, function (x) mean(x,na.rm=T)))
# colnames(test) <- c( "vals")
# test$year <- year
# test$month <- month
# 
# plot(test$vals~test$year)
# plot(aggregate(vals ~ year, data = test, FUN="mean"),type='b', pch=19, ylab="Mean Surface Chl-a") #looking at roughly 2 degrees for SST in gfdl
# points(aggregate(vals ~ year, data = test, FUN="min"),type='b', col="grey")
# points(aggregate(vals ~ year, data = test, FUN="max"),type='b', col="grey")
# # 





