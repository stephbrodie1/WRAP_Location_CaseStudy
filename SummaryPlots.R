
#Summary Plots
#WRAP Loc^3 call

#---library----
library(tidyverse)
library(ggplot2)
library(reshape2)
library(oce)
library(gganimate)
library(gifski)
library(cmocean)
library(ncdf4)
library(raster)
library(gridExtra)

#----FIGURE ONE----
#Conceptual diagram of OM and EMs
#Will be done in powerpoint.

#----FIGURE TWO: environmental timeseries & maps----
#Show time-series of environmental covariates in each domain
# 3x4 plot
#col 1 & 3: spatial maps of historical variables (use gfdl)
#col 2 & 4: time-series of covariates for each ESM +/- SE

#Spatial maps
#For each variable in gfdl, average over years 1985-2010 and over months 4,5,6
spatial_mean_list <- list()
variables <- c("ild_0.5C","oxygen_bottom","sst","temp_bottom","zoo_200m","zoo_50m")
counter = 1
for (variable in variables){
  print(variable)
  nc <- nc_open(paste0('~/Dropbox/WRAP Location^3/2d_fields/monthly/',variable,'_monthly_roms_gfdl_1980_2100.nc'))
  lat <- ncvar_get(nc, 'lat'); lat <- lat[1,]
  lon <- ncvar_get(nc, 'lon'); lon <- lon[,1]
  year <- ncvar_get(nc, 'year') 
  month <- ncvar_get(nc, 'month')
  month <- month[61:372] #only historical select years of interest 1985-2010
  month_index <- which(month==4 | month==5 | month==6)
  var <- ncvar_get(nc,nc$var[[5]]$name)
  var <- var[,,month_index]
  var_mean <- apply(var, c(1,2),FUN = function(x) {mean(x,na.rm=T)})
  colnames(var_mean) <- lat
  rownames(var_mean) <- lon
  var_mean <- melt(var_mean)
  var_mean <- na.omit(var_mean)
  spatial_mean_list[[counter]] <- var_mean
  counter = counter +1
}
summary(spatial_mean_list)

#ggplot maps
map1 <- ggplot(data=spatial_mean_list[[1]], aes(x=Var1, y=Var2))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = "Mixed Layer Depth") +   theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("dense")(256)) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))
map3 <- ggplot(data=spatial_mean_list[[3]], aes(x=Var1, y=Var2))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = "Sea Surface Temperature") +   theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("thermal")(256), name="Sea Surface Temperature") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 
map5 <- ggplot(data=spatial_mean_list[[3]], aes(x=Var1, y=Var2))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = "Integrated Zooplankton (200 m)") +   theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("algae")(256), name="Zooplankton Integrated to 200m") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 


#Trim to inshore domain
#Build raster template
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
extent(r) <- extent(raster('~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/z_.1.grd'))
rt <- rasterize(test,r)

#Mask before plotting
var_r <- rasterFromXYZ(spatial_mean_list[[4]])
template = raster('~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd')
var_r <- raster::resample(var_r, template, method="bilinear")  
var_r <- mask(var_r,rt)
var_df <- as.data.frame(rasterToPoints(var_r))
map4 <- ggplot(data=var_df, aes(x=x, y=y))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = "Bottom Temperature") +   theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("ice")(256), name="Bottom Temperature") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 

var_r <- rasterFromXYZ(spatial_mean_list[[2]])
template = raster('~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd')
var_r <- raster::resample(var_r, template, method="bilinear")  
var_r <- mask(var_r,rt)
var_df <- as.data.frame(rasterToPoints(var_r))
map2 <- ggplot(data=var_df, aes(x=x, y=y))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = "Bottom Oxygen") +  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("rain")(256), name="Bottom Oxygen") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 


var_r <- rasterFromXYZ(spatial_mean_list[[6]])
template = raster('~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd')
var_r <- raster::resample(var_r, template, method="bilinear")  
var_r <- mask(var_r,rt)
var_df <- as.data.frame(rasterToPoints(var_r))
map6 <- ggplot(data=var_df, aes(x=x, y=y))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = "Integrated Zooplankton (50 m)") +    theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("algae")(256), name="Zooplankton Integrated to 50m") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 

#Make Time-series plots: 6 variables, 3 ESMs
timeseries_mean_list <- list()
variables <- c("ild_0.5C","oxygen_bottom","sst","temp_bottom","zoo_200m","zoo_50m")
counter = 1
for (variable in variables){
  print(variable)
  esm_count <- 1
  esm_temp_list <- list()
  for (esm in c("had","gfdl","ipsl")){
    print(esm)
    if(variable=="ild_0.5C" | variable=="sst" | variable=="zoo_200m"){
    nc <- nc_open(paste0('~/Dropbox/WRAP Location^3/2d_fields/monthly/',variable,'_monthly_roms_',esm,'_1980_2100.nc'))
    lat <- ncvar_get(nc, 'lat'); lat <- lat[1,]
    lon <- ncvar_get(nc, 'lon'); lon <- lon[,1]
    month <- ncvar_get(nc, 'month')
    month <- month[61:1452] #select all years after 1985
    month_index <- which(month==4 | month==5 | month==6)
    var <- ncvar_get(nc,nc$var[[5]]$name)
    var <- var[,,month_index]
    var_mean <- as.data.frame(apply(var, c(3),FUN = function(x){mean(x,na.rm=T)}))
    colnames(var_mean) <- "value"
    var_mean$sd <- apply(var, c(3),FUN = function(x){sd(x,na.rm=T)})
    var_mean$month <- rep(c(4,5,6),116)
    var_mean$year <- rep(1985:2100,each=3)
    var_mean_agg <- aggregate(value ~ year, data=var_mean, FUN = function(x){mean(x,na.rm=T)})
    var_sd_agg <- aggregate(sd ~ year, data=var_mean, FUN = function(x){sd(x,na.rm=T)})
    var_agg <- left_join(var_mean_agg,var_sd_agg, by="year")
    esm_temp_list[[esm_count]] <- var_agg
    esm_count <- esm_count +1
    } else {
      #Trim other variables to smaller domain
      stack <- stack(list.files(paste0('~/Dropbox/WRAP Location^3/Rasters_2d_monthly/',esm,'/',variable,'/'),full.names = TRUE, pattern = '.grd'))
      stack2 <- mask(stack,rt)
      stats <- as.data.frame(cellStats(stack2,'mean'))
      colnames(stats) <- "stats"
      stats$year <- rep(1980:2100, each=12)
      stats$month <- rep(1:12, each=121)
      stats <- stats[stats$year>=1985,]
      stats_agg <- aggregate(stats ~ year, data=stats, FUN = function(x){mean(x,na.rm=T)})
      stats_agg_sd <- aggregate(stats ~ year, data=stats, FUN = function(x){sd(x,na.rm=T)})
      var_agg <- left_join(stats_agg,stats_agg_sd, by="year")
      colnames(var_agg) <- c("year","value","sd")
      esm_temp_list[[esm_count]] <- var_agg
      esm_count <- esm_count +1
    }
  }
  df1 <- esm_temp_list[[1]] 
  df2 <- esm_temp_list[[2]]
  df3 <- esm_temp_list[[3]]
  df <- rbind(df1,df2,df3)
  df$esm <- rep(c("HAD","GFDL","IPSL"), each=116)
  df$ymin <- df$value - df$sd/sqrt(3)
  df$ymax <- df$value + df$sd/sqrt(3)
  timeseries_mean_list[[counter]] <- df
  counter = counter +1
}
summary(timeseries_mean_list)

ts1 <- ggplot(data = timeseries_mean_list[[1]],aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Mixed Layer Depth", x="") +   theme(legend.position="top")
ts2 <- ggplot(data = timeseries_mean_list[[2]],aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Bottom Oxygen", x="") +   theme(legend.position="top")
ts3 <- ggplot(data = timeseries_mean_list[[3]],aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Sea Surface Temperature", x="") +   theme(legend.position="top")
ts4 <- ggplot(data = timeseries_mean_list[[4]],aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Bottom Temperature", x="") +   theme(legend.position="top")
ts5 <- ggplot(data = timeseries_mean_list[[5]],aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Integrated Zooplankton (200 m)", x="") +   theme(legend.position="top")
ts6 <- ggplot(data = timeseries_mean_list[[6]],aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Integrated Zooplankton (50 m)", x="") +   theme(legend.position="top")

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig2.tiff',res=300,units="in",height=14,width=16)
grid.arrange(map1, ts1, map2, ts2,
             map3, ts3, map4, ts4,
             map5, ts5, map6, ts6,
              nrow=3,ncol=4)
dev.off()

#----FIGURE THREE: biomass timeseries-----
#Ensemble mean and ESM comparison
# 1x3 plot for each archetype, with each plot showing ensemble mean + error for 3 ESMs. 
# Start with one species

#Load in predictions made in EM files
species_biomass <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #create ensemble mean (remove Gam_S & Gam_EST)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                               "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                               "brt_E" , "brt_ES", "brt_EST" ,  
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"), sum)
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  all_esm_agg$min <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                              "glm_E" ,"glm_ESt", "glm_ESr",   
                                              "brt_E" , "brt_ES", "brt_EST" ,  
                                              "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){min(x,na.rm=T)})
  all_esm_agg$max <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                              "glm_E" ,"glm_ESt", "glm_ESr",   
                                              "brt_E" , "brt_ES", "brt_EST" ,  
                                              "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){max(x,na.rm=T)})
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm", "min", "max"), variable.name = "EM",value.name = "abundance")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="abundance" | all_esm_agg_longer$EM=="ens_mean",]
  species_biomass[[counter]]<- all_esm_agg_longer_filtered
  counter=counter+1
}

hms_biomass <- ggplot(data=species_biomass[[1]],aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated biomass", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps_biomass <- ggplot(data=species_biomass[[2]],aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated biomass", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs_biomass <- ggplot(data=species_biomass[[3]],aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated biomass", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig3.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_biomass, cps_biomass, gfs_biomass, nrow=3)
dev.off()

#----FIGURE 3.5: biomass timeseries for whack EMs-----

#Load in predictions made in EM files
species_biomass_whack <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #look at only remove Gam_S & Gam_EST & Glm_Sr
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_S","gam_EST", "glm_Sr",
                                          "esm")]
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("abundance","gam_S","gam_EST", "glm_Sr"),c(sum),na.rm=T)
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm"), variable.name = "EM",value.name = "abundance")
  species_biomass_whack[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

hms_biomass_whack <- ggplot(data=species_biomass_whack[[1]],aes(x=year,y=abundance))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085","#8e44ad"),labels = c("Simulated biomass", "GAM S", "GAM EST","GLM S"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps_biomass_whack <- ggplot(data=species_biomass_whack[[2]],aes(x=year,y=abundance))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085","#8e44ad"),labels = c("Simulated biomass", "GAM S", "GAM EST","GLM S"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs_biomass_whack <- ggplot(data=species_biomass_whack[[3]],aes(x=year,y=abundance))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085","#8e44ad"),labels = c("Simulated biomass", "GAM S", "GAM EST","GLM S"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig3.5.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_biomass_whack, cps_biomass_whack, gfs_biomass_whack, nrow=3)
dev.off()

#----FIGURE 4: RMSE timeseries----
#Compare average RMSE for each year to average in historical period, showing ensemble mean

#Load in predictions made in EM files
species_rmse_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #create ensemble mean (remove Gam_S & Gam_EST)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  
  #Calculate RMSE
  RMSE = function(p, o){(sqrt(mean((p - o)^2)))}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(RMSE(., o=abundance)))
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  all_esm_agg$min <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){min(x,na.rm=T)})
  all_esm_agg$max <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){max(x,na.rm=T)})
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm", "min", "max"), variable.name = "EM",value.name = "RMSE")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  species_rmse_diff[[counter]]<- all_esm_agg_longer_filtered
  counter=counter+1
}

hms_rmse_diff <- ggplot(data=species_rmse_diff[[1]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps_rmse_diff <- ggplot(data=species_rmse_diff[[2]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs_rmse_diff <- ggplot(data=species_rmse_diff[[3]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig4.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_rmse_diff, cps_rmse_diff, gfs_rmse_diff, nrow=3)
dev.off()

#----FIGURE 4.5: Correlation timeseries-----
#Compare point x point correlation for each year and EM

#Load in predictions made in EM files
species_rmse_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- as.factor(rep(c("HAD","GFDL","IPSL"),each=58000))
  #create ensemble mean (remove Gam_S & Gam_EST)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  
  #Calculate RMSE
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(COR(., y=abundance)))
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  all_esm_agg$min <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){min(x,na.rm=T)})
  all_esm_agg$max <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){max(x,na.rm=T)})
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm", "min", "max"), variable.name = "EM",value.name = "RMSE")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  species_rmse_diff[[counter]]<- all_esm_agg_longer_filtered
  counter=counter+1
}

hms_rmse_diff <- ggplot(data=species_rmse_diff[[1]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps_rmse_diff <- ggplot(data=species_rmse_diff[[2]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs_rmse_diff <- ggplot(data=species_rmse_diff[[3]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig4.5.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_rmse_diff, cps_rmse_diff, gfs_rmse_diff, nrow=3)
dev.off()



#----FIGURE 5: Latitudinal COG timeseries----
#Compare COG for each year to average in historical period, showing ensemble mean

#Load in predictions made in EM files
species_cog_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #create ensemble mean (remove Gam_S & Gam_EST)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  
  #Calculate COG
  COG = function(x, w){weighted.mean(x=x, w=w)}
  
  #aggregate over space
  all_esm_agg <- all_esm_allyears %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(weighted.mean(., x=lat)),na.rm=T)
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  all_esm_agg$min <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){min(x,na.rm=T)})
  all_esm_agg$max <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){max(x,na.rm=T)})
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm", "min", "max"), variable.name = "EM",value.name = "COG")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="abundance" | all_esm_agg_longer$EM=="ens_mean",]
  species_cog_diff[[counter]]<- all_esm_agg_longer_filtered
  counter=counter+1
}

hms_cog <- ggplot(data=species_cog_diff[[1]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Latitudinal COG")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps_cog <- ggplot(data=species_cog_diff[[2]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Latitudinal COG")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs_cog <- ggplot(data=species_cog_diff[[3]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: Latitudinal RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig5.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_cog, cps_cog, gfs_cog, nrow=3)
dev.off()


#----FIGURE 6: Longitudinal COG timeseries----
#Compare COG for each year to average in historical period, showing ensemble mean

#Load in predictions made in EM files
species_cog_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
    had_allyears$mlp_ES <- NA
    gfdl_allyears$mlp_E <- NA
    gfdl_allyears$mlp_ES <- NA
    ipsl_allyears$mlp_ES <- NA
    had_allyears <- had_allyears[,c(1:37,38,40,39)]
    gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
    ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #create ensemble mean (remove Gam_S & Gam_EST)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  
  #Calculate COG
  COG = function(x, w){weighted.mean(x=x, w=w)}
  
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                              "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                              "brt_E" , "brt_ES", "brt_EST" ,  
                                                                              "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(weighted.mean(., x=lon)),na.rm=T)
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  all_esm_agg$min <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){min(x,na.rm=T)})
  all_esm_agg$max <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){max(x,na.rm=T)})
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm", "min", "max"), variable.name = "EM",value.name = "COG")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="abundance" | all_esm_agg_longer$EM=="ens_mean",]
  species_cog_diff[[counter]]<- all_esm_agg_longer_filtered
  counter=counter+1
}

hms_cog <- ggplot(data=species_cog_diff[[1]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: longitudinal COG")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps_cog <- ggplot(data=species_cog_diff[[2]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: longitudinal COG")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs_cog <- ggplot(data=species_cog_diff[[3]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: longitudinal COG")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig6.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_cog, cps_cog, gfs_cog, nrow=3)
dev.off()

#----FIGURE 6.5: Distance to coast timeseries----

#Load in predictions made in EM files
species_cog_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #create ensemble mean (remove Gam_S & Gam_EST)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  
  #Calculate distance from coast COG
  #Create mask for distance raster #NOTE MASK ONLY NEEDED FOR CPS AND Groundfish
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
  extent(r) <- extent(raster('~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/z_.1.grd'))
  rt <- rasterize(test,r)
  
  dist <- readRDS('~/PROJECTS/WRAP Location/r_distcoast.rds')
  # dist <- raster::mask(dist,rt)
  dist_df <- as.data.frame(rasterToPoints(dist))
  colnames(dist_df) <- c("lon","lat","distX")
  dist_df$lat <- round(dist_df$lat,2) #CRITICAL STEP
  dist_df$lon <- round(dist_df$lon,2) #CRITICAL STEP

  
  all_esm_allyears_dist <- merge(all_esm_allyears, dist_df, by=c("lon","lat"), all.x=T, all.y=F)
 
  # COG = function(x, w){weighted.mean(x=x, w=w)}
  
  #aggregate over space
  all_esm_agg <-all_esm_allyears_dist %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                              "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                              "brt_E" , "brt_ES", "brt_EST" ,  
                                                                              "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(weighted.mean(., x=distX)),na.rm=T)
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  all_esm_agg$min <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){min(x,na.rm=T)})
  all_esm_agg$max <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){max(x,na.rm=T)})
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm", "min", "max"), variable.name = "EM",value.name = "COG")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="abundance" | all_esm_agg_longer$EM=="ens_mean",]
  species_cog_diff[[counter]]<- all_esm_agg_longer_filtered
  counter=counter+1
}

hms_cog <- ggplot(data=species_cog_diff[[1]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: distance to coast COG")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps_cog <- ggplot(data=species_cog_diff[[2]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: distance to coast COG")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs_cog <- ggplot(data=species_cog_diff[[3]],aes(x=year,y=COG, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated COG", "Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: distance to coast COG")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig6.5.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_cog, cps_cog, gfs_cog, nrow=3)
dev.off()






#----FIGURE 7: Forecast years Correlation (also RMSE)----
#Fit and performance metrics for each EM and OM

#Load in predictions made in EM files
species_predperf <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  # had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  # gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  # ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=4500)
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor","gam_EST", "gam_S",
                                          "glm_E" ,"glm_ESt", "glm_Sr", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  #Calculate metrics 
  COR = function(x, y){cor(x,y,method="spearman",use = "na.or.complete")}
  RMSE = function(p, o){(sqrt(mean((p - o)^2)))}
  #COR: aggregate over space and time
  all_esm_agg_cor <-all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor","gam_EST", "gam_S",
                                                                              "glm_E" ,"glm_ESt", "glm_Sr", "glm_ESr",   
                                                                              "brt_E" , "brt_ES", "brt_EST" ,  
                                                                              "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(COR(., y=abundance)))
  all_esm_agg_cor_longer <- melt(all_esm_agg_cor,id=c("esm"), variable.name = "EM",value.name = "COR")
  #RMSE: aggregate over space and time
  all_esm_agg_rmse <-all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor","gam_EST", "gam_S",
                                                                             "glm_E" ,"glm_ESt", "glm_Sr", "glm_ESr",   
                                                                             "brt_E" , "brt_ES", "brt_EST" ,  
                                                                             "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(RMSE(., o=abundance)))
  all_esm_agg_rmse_longer <- melt(all_esm_agg_rmse,id=c("esm"), variable.name = "EM",value.name = "RMSE")
  #Make one data frame
  all_esm_agg_longer <- left_join(all_esm_agg_cor_longer,all_esm_agg_rmse_longer, by =c("esm","EM"))
  species_predperf[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

hms_cor <- ggplot(data=species_predperf[[1]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="HMS cor: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)+
  ylim(c(0,1))

cps_cor <- ggplot(data=species_predperf[[2]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="CPS cor: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)+
  ylim(c(0,1))

gfs_cor <- ggplot(data=species_predperf[[3]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="Groundfish cor: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)+
  ylim(c(0,1))

hms_rmse <- ggplot(data=species_predperf[[1]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="HMS RMSE: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

cps_rmse <- ggplot(data=species_predperf[[2]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="CPS RMSE: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

gfs_rmse <- ggplot(data=species_predperf[[3]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="Groundfish RMSE: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig7.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_cor, cps_cor, gfs_cor, nrow=3)
dev.off()


#-----FIGURE 8: Fitted years Correlation (also RMSE)-------
#Load in predictions made in EM files
species_predperf <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  # had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  # gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  # ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist)
  gfdl_allyears <- rbind(gfdl_dat_hist)
  ipsl_allyears <- rbind(ipsl_dat_hist)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=13000)
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor","gam_EST", "gam_S",
                                          "glm_E" ,"glm_ESt", "glm_Sr", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  #Calculate metrics 
  COR = function(x, y){cor(x,y,method="spearman",use = "na.or.complete")}
  RMSE = function(p, o){(sqrt(mean((p - o)^2)))}
  #COR: aggregate over space and time
  all_esm_agg_cor <-all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor","gam_EST", "gam_S",
                                                                             "glm_E" ,"glm_ESt", "glm_Sr", "glm_ESr",   
                                                                             "brt_E" , "brt_ES", "brt_EST" ,  
                                                                             "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(COR(., y=abundance)))
  all_esm_agg_cor_longer <- melt(all_esm_agg_cor,id=c("esm"), variable.name = "EM",value.name = "COR")
  #RMSE: aggregate over space and time
  all_esm_agg_rmse <-all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor","gam_EST", "gam_S",
                                                                              "glm_E" ,"glm_ESt", "glm_Sr", "glm_ESr",   
                                                                              "brt_E" , "brt_ES", "brt_EST" ,  
                                                                              "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(RMSE(., o=abundance)))
  all_esm_agg_rmse_longer <- melt(all_esm_agg_rmse,id=c("esm"), variable.name = "EM",value.name = "RMSE")
  #Make one data frame
  all_esm_agg_longer <- left_join(all_esm_agg_cor_longer,all_esm_agg_rmse_longer, by =c("esm","EM"))
  species_predperf[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

hms_cor <- ggplot(data=species_predperf[[1]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="HMS cor: fitted years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)+
  ylim(c(0,1))

cps_cor <- ggplot(data=species_predperf[[2]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="CPS cor: fitted years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)+
  ylim(c(0,1))

gfs_cor <- ggplot(data=species_predperf[[3]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="Groundfish cor: fitted years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)+
  ylim(c(0,1))

hms_rmse <- ggplot(data=species_predperf[[1]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="HMS RMSE: fitted years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

cps_rmse <- ggplot(data=species_predperf[[2]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="CPS RMSE: fitted years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

gfs_rmse <- ggplot(data=species_predperf[[3]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="Groundfish RMSE: fitted years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig8.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_cor, cps_cor, gfs_cor, nrow=3)
dev.off()

#---FIGURE 9: Forecast RMSE of good models ----
#Load in predictions made in EM files
species_predperf <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  # had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  # gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  # ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_fcast)
  if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
    had_allyears$mlp_ES <- NA
    gfdl_allyears$mlp_E <- NA
    gfdl_allyears$mlp_ES <- NA
    ipsl_allyears$mlp_ES <- NA
    had_allyears <- had_allyears[,c(1:37,38,40,39)]
    gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
    ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=4500)
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  #Calculate metrics 
  COR = function(x, y){cor(x,y,method="spearman",use = "na.or.complete")}
  RMSE = function(p, o){(sqrt(mean((p - o)^2)))}
  #COR: aggregate over space and time
  all_esm_agg_cor <-all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                             "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                             "brt_E" , "brt_ES", "brt_EST" ,  
                                                                             "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(COR(., y=abundance)))
  all_esm_agg_cor_longer <- melt(all_esm_agg_cor,id=c("esm"), variable.name = "EM",value.name = "COR")
  #RMSE: aggregate over space and time
  all_esm_agg_rmse <-all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                              "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                              "brt_E" , "brt_ES", "brt_EST" ,  
                                                                              "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(RMSE(., o=abundance)))
  all_esm_agg_rmse_longer <- melt(all_esm_agg_rmse,id=c("esm"), variable.name = "EM",value.name = "RMSE")
  #Make one data frame
  all_esm_agg_longer <- left_join(all_esm_agg_cor_longer,all_esm_agg_rmse_longer, by =c("esm","EM"))
  species_predperf[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

hms_cor <- ggplot(data=species_predperf[[1]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  facet_wrap(~esm)

cps_cor <- ggplot(data=species_predperf[[2]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  facet_wrap(~esm)

gfs_cor <- ggplot(data=species_predperf[[3]],aes(x=EM,y=COR))+
  geom_bar(stat = "identity")+
  theme_classic()+
  facet_wrap(~esm)

hms_rmse <- ggplot(data=species_predperf[[1]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="HMS RMSE: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

cps_rmse <- ggplot(data=species_predperf[[2]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="CPS RMSE: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

gfs_rmse <- ggplot(data=species_predperf[[3]],aes(x=EM,y=RMSE))+
  geom_bar(stat = "identity")+
  theme_classic()+
  labs(x="",y="Groundfish RMSE: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig9.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_rmse, cps_rmse, gfs_rmse, nrow=3)
dev.off()

#-----FIGURE 10: environmental density plots-----

#Load in predictions made in EM files
species_enviro_density <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
    had_allyears$mlp_ES <- NA
    gfdl_allyears$mlp_E <- NA
    gfdl_allyears$mlp_ES <- NA
    ipsl_allyears$mlp_ES <- NA
    had_allyears <- had_allyears[,c(1:37,38,40,39)]
    gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
    ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #keep necessary columns
  if (s =="Albacore EMs"){
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","temp","chla","mld","esm")]
  all_esm_agg_longer <- melt(all_esm_allyears,id=c("lon","lat","year", "esm"), variable.name = "enviro_name",value.name = "enviro_value")
  }
  if (s =="Anchovy EMs"){
    all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","temp","chla","esm")]
    all_esm_agg_longer <- melt(all_esm_allyears,id=c("lon","lat","year", "esm"), variable.name = "enviro_name",value.name = "enviro_value")
  }
  if (s =="Groundfish EMs"){
    all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","btemp","O2","esm")]
    all_esm_agg_longer <- melt(all_esm_allyears,id=c("lon","lat","year", "esm"), variable.name = "enviro_name",value.name = "enviro_value")
  }
  
  all_esm_agg_longer$period <- ifelse(all_esm_agg_longer$year<=2010,"Historical","Forecast")
  species_enviro_density[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

hms_enviro_density <- ggplot(data=species_enviro_density[[1]],aes(x=enviro_value,fill=period))+
  geom_density(alpha=0.4)+
  theme_classic()+
  labs(x="",y="")+
  theme(legend.title = element_blank())+
  facet_wrap(c("esm","enviro_name"),scales="free")

cps_enviro_density <- ggplot(data=species_enviro_density[[2]],aes(x=enviro_value,fill=period))+
  geom_density(alpha=0.4)+
  theme_classic()+
  labs(x="",y="")+
  theme(legend.title = element_blank())+
  facet_wrap(c("enviro_name", "esm"),scales="free")

gfs_enviro_density <- ggplot(data=species_enviro_density[[3]],aes(x=enviro_value,fill=period))+
  geom_density(alpha=0.4)+
  theme_classic()+
  labs(x="",y="")+
  theme(legend.title = element_blank())+
  facet_wrap(c("enviro_name","esm"),scales="free")


tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig10_HMS.tiff',res=300,units="in",height=10,width=13)
grid.arrange(hms_enviro_density)
dev.off()

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig10_CPS.tiff',res=300,units="in",height=10,width=13)
grid.arrange(cps_enviro_density)
dev.off()

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig10_GFS.tiff',res=300,units="in",height=10,width=13)
grid.arrange(gfs_enviro_density)
dev.off()


#----FIGURE 11: Partitioning Variance-----
#Load in predictions made in EM files
variance_partition <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  # had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  # gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  # ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=4500)
  #create ensemble mean (remove Gam_S & Gam_EST)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),c(sum),na.rm=T)
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  all_esm_agg$min <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){min(x,na.rm=T)})
  all_esm_agg$max <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){max(x,na.rm=T)})
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm", "min", "max"), variable.name = "EM",value.name = "abundance")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  variance_partition[[counter]]<- all_esm_agg_longer_filtered
  counter=counter+1
}


hms_variance_partition <- ggplot(data=variance_partition[[1]],aes(x=year,y=abundance, by=esm, ymin=min,ymax=max))+
  geom_line(aes(col=esm))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085"))+
  geom_ribbon(alpha=0.2, aes(fill=esm))+
  scale_fill_manual(values=c("#2980b9","#c0392b","#16a085"))+
  theme_classic()+
  labs(x="",y="HMS Biomass")+
  theme(legend.title = element_blank())

cps_variance_partition <- ggplot(data=variance_partition[[2]],aes(x=year,y=abundance, by=esm, ymin=min,ymax=max))+
  geom_line(aes(col=esm))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085"))+
  geom_ribbon(alpha=0.2, aes(fill=esm))+
  scale_fill_manual(values=c("#2980b9","#c0392b","#16a085"))+
  theme_classic()+
  labs(x="",y="CPS Biomass")+
  theme(legend.title = element_blank())

gfs_variance_partition <- ggplot(data=variance_partition[[3]],aes(x=year,y=abundance, by=esm, ymin=min,ymax=max))+
  geom_line(aes(col=esm))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085"))+
  geom_ribbon(alpha=0.2, aes(fill=esm))+
  scale_fill_manual(values=c("#2980b9","#c0392b","#16a085"))+
  theme_classic()+
  labs(x="",y="Groundfish Biomass")+
  theme(legend.title = element_blank())


tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig11.tiff',res=300,units="in",height=10,width=8)
grid.arrange(hms_variance_partition, cps_variance_partition, gfs_variance_partition, nrow=3)
dev.off()

#---FIGURE 12: Partition variance more closely-----
variance_partition_summary <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  # if (s =="Groundfish EMs"){ #do this because MLPs weren't working correctly
  #   had_allyears$mlp_ES <- NA
  #   gfdl_allyears$mlp_E <- NA
  #   gfdl_allyears$mlp_ES <- NA
  #   ipsl_allyears$mlp_ES <- NA
  #   had_allyears <- had_allyears[,c(1:37,38,40,39)]
  #   gfdl_allyears <- gfdl_allyears[,c(1:37,39,40,38)]
  #   ipsl_allyears <- ipsl_allyears[,c(1:37,38,40,39)]
  # }
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #create ensemble mean (remove Gam_S & Gam_EST)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),c(sum),na.rm=T)
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  
  all_esm_agg$sd_esm <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST")],1,FUN=function(x){sd(x,na.rm=T)})
  variance_partition_summary[[counter]]<- all_esm_agg
  counter=counter+1
}

#HMS
all_esm_agg_esm <- aggregate(sd_esm~year,data=variance_partition_summary[[1]],FUN=mean)
all_esm_agg_em <- aggregate(abundance~year,data=variance_partition_summary[[1]],FUN="sd")
combined <- left_join(all_esm_agg_em,all_esm_agg_esm,by="year")
combined <- melt(combined,id=c("year"), variable.name = "error_source",value.name = "value")
hline <- exp(mean(log(variance_partition_summary[[1]]$sd_esm[variance_partition_summary[[1]]$year<=2010])) + sd(log(variance_partition_summary[[1]]$sd_esm[variance_partition_summary[[1]]$year<=2010]))*2)
hms_variancepartition <- ggplot(data=combined,aes(x=year, y=value, by=error_source))+
  geom_line(aes(col=error_source))+
  theme_classic()+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Error among ESM", "Error among EM"))+
  labs(x="",y="Standard deviation of HMS Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  geom_hline(yintercept = hline)
combined$year[which(combined$value[combined$error_source=="sd_esm"]>hline)]


#CPS
all_esm_agg_esm <- aggregate(sd_esm~year,data=variance_partition_summary[[2]],FUN=mean)
all_esm_agg_em <- aggregate(abundance~year,data=variance_partition_summary[[2]],FUN="sd")
combined <- left_join(all_esm_agg_em,all_esm_agg_esm,by="year")
combined <- melt(combined,id=c("year"), variable.name = "error_source",value.name = "value")
hline <- exp(mean(log(variance_partition_summary[[2]]$sd_esm[variance_partition_summary[[2]]$year<=2010])) + sd(log(variance_partition_summary[[2]]$sd_esm[variance_partition_summary[[2]]$year<=2010]))*2)
cps_variancepartition <- ggplot(data=combined,aes(x=year, y=value, by=error_source))+
  geom_line(aes(col=error_source))+
  theme_classic()+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Error among ESM", "Error among EM"))+
  labs(x="",y="Standard deviation of CPS Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  geom_hline(yintercept = hline)
combined$year[which(combined$value[combined$error_source=="sd_esm"]>hline)]

#GFS
all_esm_agg_esm <- aggregate(sd_esm~year,data=variance_partition_summary[[3]],FUN=mean)
all_esm_agg_em <- aggregate(abundance~year,data=variance_partition_summary[[3]],FUN="sd")
combined <- left_join(all_esm_agg_em,all_esm_agg_esm,by="year")
combined <- melt(combined,id=c("year"), variable.name = "error_source",value.name = "value")
hline <- exp(mean(log(variance_partition_summary[[3]]$sd_esm[variance_partition_summary[[3]]$year<=2010])) + sd(log(variance_partition_summary[[3]]$sd_esm[variance_partition_summary[[3]]$year<=2010]))*2)
gfs_variancepartition <- ggplot(data=combined,aes(x=year, y=value, by=error_source))+
  geom_line(aes(col=error_source))+
  theme_classic()+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Error among ESM", "Error among EM"))+
  labs(x="",y="Standard deviation of Groundfish Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  geom_hline(yintercept = hline)
combined$year[which(combined$value[combined$error_source=="sd_esm"]>hline)]

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig12.tiff',res=300,units="in",height=10,width=8)
grid.arrange(hms_variancepartition, cps_variancepartition, gfs_variancepartition, nrow=3)
dev.off()





#-----ANIMATED SPECIES MAPS-----
setwd('~/PROJECTS/WRAP Location/Manuscript/Figures/')

#HMS
for (esm in c("had","gfdl","ipsl")){
  print(esm)
  dat <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/",esm,"/dat_albacore_all.rds"))
  title <- ifelse(esm=="had","HAD",ifelse(esm=="gfdl","GFDL","IPSL"))
  hms_animate_map <- ggplot(data = dat, aes(x=lon,y=lat))+
    geom_tile(aes(fill=abundance))+
    theme_classic() +  labs(y="", x="") + 
    theme(legend.position="right",legend.title = element_blank())+
    theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
    scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
    annotation_map(map_data("world"), colour = "black", fill="grey50")+
    coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
    scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
    transition_time(year)+
    ease_aes("linear") +
    labs(title=paste0(title,": HMS Biomass {frame_time}")) #takes a few mins
  hms_animated <- gganimate::animate(hms_animate_map,nframes = 116, fps=3, renderer = gifski_renderer())#renders in 
  anim_save(paste0('AnimatedMap_Albacore_',esm,'.gif'), hms_animated)
}

#CPS
for (esm in c("had","gfdl","ipsl")){
  print(esm)
  dat <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/",esm,"/dat_anchovy_all.rds"))
  title <- ifelse(esm=="had","HAD",ifelse(esm=="gfdl","GFDL","IPSL"))
  hms_animate_map <- ggplot(data = dat, aes(x=lon,y=lat))+
    geom_tile(aes(fill=abundance))+
    theme_classic() +  labs(y="", x="") + 
    theme(legend.position="right",legend.title = element_blank())+
    theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
    scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
    annotation_map(map_data("world"), colour = "black", fill="grey50")+
    coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
    scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
    transition_time(year)+
    ease_aes("linear") +
    labs(title=paste0(title,": CPS Biomass {frame_time}")) #takes a few mins
  hms_animated <- gganimate::animate(hms_animate_map,nframes = 116, fps=3, renderer = gifski_renderer())#renders in 
  anim_save(paste0('AnimatedMap_Anchovy_',esm,'.gif'), hms_animated)
}

#Groundfish
for (esm in c("had","gfdl","ipsl")){
  print(esm)
  dat <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/",esm,"/dat_groundfish_all.rds"))
  title <- ifelse(esm=="had","HAD",ifelse(esm=="gfdl","GFDL","IPSL"))
  hms_animate_map <- ggplot(data = dat, aes(x=lon,y=lat))+
    geom_tile(aes(fill=abundance))+
    theme_classic() +  labs(y="", x="") + 
    theme(legend.position="right",legend.title = element_blank())+
    theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
    scale_fill_gradientn(colours = cmocean("matter")(256),limits = c(0, max(dat$abundance))) +
    annotation_map(map_data("world"), colour = "black", fill="grey50")+
    coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
    scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
    transition_time(year)+
    ease_aes("linear") +
    labs(title=paste0(title,": Groundfish Biomass {frame_time}")) #takes some time 
  hms_animated <- gganimate::animate(hms_animate_map,nframes = 116, fps=3, renderer = gifski_renderer())#renders in 
  anim_save(paste0('AnimatedMap_Sablefish_',esm,'.gif'), hms_animated)
}


#-----Plot Environmental Trends in each domain-----
#Load in cps files
cps_had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/","had","/dat_hist_results_full.rds"))
cps_had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/","had","/dat_fcast_results_full.rds"))
cps_gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/","gfdl","/dat_hist_results_full.rds"))
cps_gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/","gfdl","/dat_fcast_results_full.rds"))
cps_ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/","ipsl","/dat_hist_results_full.rds"))
cps_ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/","ipsl","/dat_fcast_results_full.rds"))
#Rbind historical and forecast years
cps_had_allyears <- rbind(cps_had_dat_hist,cps_had_dat_fcast)
cps_gfdl_allyears <- rbind(cps_gfdl_dat_hist,cps_gfdl_dat_fcast)
cps_ipsl_allyears <- rbind(cps_ipsl_dat_hist,cps_ipsl_dat_fcast)

#Load in HMS files
hms_had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","had","/dat_hist_results_full.rds"))
hms_had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","had","/dat_fcast_results_full.rds"))
hms_gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","gfdl","/dat_hist_results_full.rds"))
hms_gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","gfdl","/dat_fcast_results_full.rds"))
hms_ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","ipsl","/dat_hist_results_full.rds"))
hms_ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","ipsl","/dat_fcast_results_full.rds"))
#Rbind historical and forecast years
hms_had_allyears <- rbind(hms_had_dat_hist,hms_had_dat_fcast)
hms_gfdl_allyears <- rbind(hms_gfdl_dat_hist,hms_gfdl_dat_fcast)
hms_ipsl_allyears <- rbind(hms_ipsl_dat_hist,hms_ipsl_dat_fcast)

#Load in groundfish files
gfs_had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/","had","/dat_hist_results_full.rds"))
gfs_had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/","had","/dat_fcast_results_full.rds"))
gfs_gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/","gfdl","/dat_hist_results_full.rds"))
gfs_gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/","gfdl","/dat_fcast_results_full.rds"))
gfs_ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/","ipsl","/dat_hist_results_full.rds"))
gfs_ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/","ipsl","/dat_fcast_results_full.rds"))
#Rbind historical and forecast years
gfs_had_allyears <- rbind(gfs_had_dat_hist,gfs_had_dat_fcast)
gfs_gfdl_allyears <- rbind(gfs_gfdl_dat_hist,gfs_gfdl_dat_fcast)
gfs_ipsl_allyears <- rbind(gfs_ipsl_dat_hist,gfs_ipsl_dat_fcast)

#HMS plot
par(mfrow=c(3,1))
plot(aggregate(temp~year,hms_had_allyears,FUN="mean"),type="l",  lwd=2,ylab="Mean Temperature")
lines(aggregate(temp~year,hms_gfdl_allyears,FUN="mean"),type="l",  lwd=2, col="blue")
lines(aggregate(temp~year,hms_ipsl_allyears,FUN="mean"),type="l",  lwd=2, col="red")
abline(v=2010)

plot(aggregate(zoo_200~year,hms_had_allyears,FUN="mean"),type="l",  lwd=2,ylab="Mean Zooplankton")
lines(aggregate(zoo_200~year,hms_gfdl_allyears,FUN="mean"),type="l",  lwd=2, col="blue")
lines(aggregate(zoo_200~year,hms_ipsl_allyears,FUN="mean"),type="l",  lwd=2, col="red")
abline(v=2010)

plot(aggregate(mld~year,hms_had_allyears,FUN="mean"),type="l",  lwd=2,ylab="Mean Mixed Layer Depth")
lines(aggregate(mld~year,hms_gfdl_allyears,FUN="mean"),type="l",  lwd=2, col="blue")
lines(aggregate(mld~year,hms_ipsl_allyears,FUN="mean"),type="l",  lwd=2,col="red")
abline(v=2010)

#CPS plot
par(mfrow=c(2,1))
plot(aggregate(temp~year,cps_had_allyears,FUN="mean"),type="l",  lwd=2,ylab="Mean Temperature")
lines(aggregate(temp~year,cps_gfdl_allyears,FUN="mean"),type="l",  lwd=2, col="blue")
lines(aggregate(temp~year,cps_ipsl_allyears,FUN="mean"),type="l",  lwd=2, col="red")
abline(v=2010)

plot(aggregate(zoo_50~year,cps_had_allyears,FUN="mean"),type="l",  lwd=2,ylab="Mean Zooplankton")
lines(aggregate(zoo_50~year,cps_gfdl_allyears,FUN="mean"),type="l",  lwd=2, col="blue")
lines(aggregate(zoo_50~year,cps_ipsl_allyears,FUN="mean"),type="l",  lwd=2, col="red")
abline(v=2010)

#GFS plot
par(mfrow=c(2,1))
plot(aggregate(btemp~year,gfs_had_allyears,FUN="mean"),type="l",  lwd=2,ylab="Mean Bottom Temperature")
lines(aggregate(btemp~year,gfs_gfdl_allyears,FUN="mean"),type="l",  lwd=2, col="blue")
lines(aggregate(btemp~year,gfs_ipsl_allyears,FUN="mean"),type="l",  lwd=2, col="red")
abline(v=2010)

plot(aggregate(O2~year,gfs_had_allyears,FUN="mean"),type="l",  lwd=2,ylab="Mean Bottom Oxygen", ylim=c(35,65))
lines(aggregate(O2~year,gfs_gfdl_allyears,FUN="mean"),type="l",  lwd=2, col="blue")
lines(aggregate(O2~year,gfs_ipsl_allyears,FUN="mean"),type="l",  lwd=2, col="red")
abline(v=2010)


#-----Plot true OM abundance across ESMs----
#HMS plot
par(mfrow=c(3,1))
plot(aggregate(abundance~year,hms_had_allyears,FUN="sum"),type="l",  lwd=2,ylab="Biomass (mt)", main="HMS")
lines(aggregate(abundance~year,hms_gfdl_allyears,FUN="sum"),type="l",  lwd=2,col="blue")
lines(aggregate(abundance~year,hms_ipsl_allyears,FUN="sum"),type="l",  lwd=2,col='red')

plot(aggregate(abundance~year,cps_had_allyears,FUN="sum"),type="l",  lwd=2,ylab="Biomass (mt)", main="CPS")
lines(aggregate(abundance~year,cps_gfdl_allyears,FUN="sum"),type="l",  lwd=2,col="blue")
lines(aggregate(abundance~year,cps_ipsl_allyears,FUN="sum"),type="l",  lwd=2,col='red')

plot(aggregate(abundance~year,gfs_had_allyears,FUN="sum"),type="l",  lwd=2,ylab="Biomass (mt)", main="GFS")
lines(aggregate(abundance~year,gfs_gfdl_allyears,FUN="sum"),type="l",  lwd=2,col="blue")
lines(aggregate(abundance~year,gfs_ipsl_allyears,FUN="sum"),type="l",  lwd=2,col='red')





  
  
  
  
  
  