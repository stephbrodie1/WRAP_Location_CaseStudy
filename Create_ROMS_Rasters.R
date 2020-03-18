
#Convert ROMS netcdf files into rasters
#Rasters need to be used for the virtualspecies simulation

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

#-----Load in Dropbox files----
#note: download files locally. Netcdfs are contained in a folder called 2d_fields
files_long <- list.files('2d_fields/monthly', full.names = TRUE)
files_short <- list.files('2d_fields/monthly', full.names = FALSE)

for (f in 5:10){ #skipping bottom layer depth as it is unique --> Steph check if it can go in the loop
  print(files_long[f])
  nc <- nc_open(files_long[f])
  # print(nc) #run if you want more details on file
  lat <- ncvar_get(nc, 'lat'); lat <- lat[1,]
  lon <- ncvar_get(nc, 'lon'); lon <- lon[,1]
  year <- ncvar_get(nc, 'year') 
  month <- ncvar_get(nc, 'month')
  name <-  names(nc$var)[5]
  tmp.array <- ncvar_get(nc,name)
  
  #Loop through every time step in ncdf
  for (i in 1:1452){
    r <- raster(t(tmp.array[,,i]),
                xmn=min(lon), xmx=max(lon),
                ymn=min(lat), ymx=max(lat), 
                crs=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
    r <- flip(r,2)
    #plot(r)
    
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
    
    #create name of raster file
    fname <- paste0(var[1],"_",var[2],"_",gcm_folder,"_",year[i],"_month",month[i],".grd")
    
    #save raster
    writeRaster(r, paste0(folder,"/", fname))
  }
  nc_close(nc)
}




