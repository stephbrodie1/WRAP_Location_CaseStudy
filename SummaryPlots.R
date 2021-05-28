
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
library(gratia)
library(dominanceanalysis)
# remotes::install_github("densitymodelling/dsmextra")
library(dsmextra)

#----FIGURE ONE----
#Conceptual diagram of OM and EMs
#Will be done in powerpoint.

#----FIGURE 2: environmental timeseries & maps----
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

#----FIGURE 3: biomass timeseries-----
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

#----FIGURE 3.1: biomass timeseries by region-----
#Ensemble mean and ESM comparison
# 1x3 plot for each archetype, with each plot showing ensemble mean + error for 3 ESMs. 
#For each result, plot regional 30-34.5, 34.5-40, 40-48

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
  
  #Assign region
  all_esm_allyears$region <- ifelse(all_esm_allyears$lat<=34.5,'South', 
                                    ifelse(all_esm_allyears$lat>=40, "North", "Central"))
    
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm, region) %>% summarise_at(c("abundance","gam_E", "gam_ES", "gam_ECor",
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
  # all_esm_agg$relative_biomass <- all_esm_agg$ens_mean / all_esm_agg$abundance
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm","region", "min", "max"), variable.name = "EM",value.name = "abundance")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  species_biomass[[counter]]<- all_esm_agg_longer_filtered
  counter=counter+1
}

hms_biomass <- ggplot(data=species_biomass[[1]],aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=region))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps_biomass <- ggplot(data=species_biomass[[2]],aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=region))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs_biomass <- ggplot(data= species_biomass[[3]],aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=region))+
  scale_color_manual(values=c("#2980b9","#c0392b","#16a085"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish Biomass")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig3.1.tiff',res=300,units="in",height=10,width=13)
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
  RMSE = function(p, o){(sqrt(mean((p - o)^2))) / mean(o)}
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




#----FIGURE 4.6: Maps of correlations------
#Map of correlations averaged over ~10years
#Load in predictions made in EM files
species_cor_spatial <- list()
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
  #Allocate decade 
  # all_esm_allyears$decade <- ifelse(all_esm_allyears$year<=2010,1985, #OLD CODE
  #                                   ifelse(all_esm_allyears$year>=2011 & all_esm_allyears$year<2100 , all_esm_allyears$year - all_esm_allyears$year %% 10,
  #                                          ifelse( all_esm_allyears$year==2100,2090, NA)))
  all_esm_allyears$decade <- ifelse(all_esm_allyears$year<=2010,"1985-2010",
                                    ifelse(all_esm_allyears$year>=2011 & all_esm_allyears$year<=2040, "2011-2040",
                                           ifelse(all_esm_allyears$year>=2041 & all_esm_allyears$year<=2070, "2041-2070",
                                                  ifelse(all_esm_allyears$year>=2071, "2071-2100", NA))))
  all_esm_allyears$lat_round <- round(all_esm_allyears$lat,0)
  all_esm_allyears$lon_round <- round(all_esm_allyears$lon,0)
  
  #Calculate correlation
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate by ESM and grid cell
  all_esm_agg <-all_esm_allyears %>% group_by(lon_round, lat_round, decade, esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(COR(., y=abundance)))
  all_esm_agg_longer <- melt(all_esm_agg,id=c("lat_round","lon_round","decade", "esm"), variable.name = "EM",value.name = "COR")
  # all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  species_cor_spatial[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

models <- c( "gam_E","gam_ES", "gam_ECor",
             "glm_E" ,"glm_ESt", "glm_ESr",   
             "brt_E" , "brt_ES", "brt_EST" ,  
             "mlp_E" , "mlp_ES" ,"mlp_EST","ens_mean")
for (m in models){
  print(m)
  dat <- species_cor_spatial[[1]]
  g1 <- ggplot(dat[dat$EM==m,],aes(x=lon_round,y=lat_round))+
    geom_tile(aes(fill=COR))+
    scale_fill_gradient2(
      low = 'red', mid = 'white', high = 'blue',
      midpoint = 0.5, guide = 'colourbar', aesthetics = 'fill', limits=c(-1,1))+
    theme_classic()+
    annotation_map(map_data("world"), colour = "black", fill=grey(0.7))+
    coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48))+
    facet_wrap(~esm + decade, nrow=3)+
    labs(x="",y="")+
    theme(legend.title = element_blank())
  
  dat <- species_cor_spatial[[2]]
  g2 <- ggplot(dat[dat$EM==m,],aes(x=lon_round,y=lat_round))+
    geom_tile(aes(fill=COR))+
    scale_fill_gradient2(
      low = 'red', mid = 'white', high = 'blue',
      midpoint = 0.5, guide = 'colourbar', aesthetics = 'fill', limits=c(-1,1))+
    theme_classic()+
    annotation_map(map_data("world"), colour = "black", fill=grey(0.7))+
    coord_quickmap(xlim=c(-127,-115.8),ylim=c(30,48))+
    facet_wrap(~esm + decade, nrow=3)+
    scale_x_continuous(breaks = round(seq(min(dat$lon_round), max(dat$lon_round), by =4),0))+
    theme(legend.title = element_blank())
  
  dat <- species_cor_spatial[[3]]
  g3 <- ggplot(dat[dat$EM==m,],aes(x=lon_round,y=lat_round))+
    geom_tile(aes(fill=COR))+
    scale_fill_gradient2(
      low = 'red', mid = 'white', high = 'blue',
      midpoint = 0.5, guide = 'colourbar', aesthetics = 'fill', limits=c(-1,1))+
    theme_classic()+
    annotation_map(map_data("world"), colour = "black", fill=grey(0.7))+
    coord_quickmap(xlim=c(-127,-115.8),ylim=c(30,48))+
    facet_wrap(~esm + decade, nrow=3)+
    scale_x_continuous(breaks = round(seq(min(dat$lon_round), max(dat$lon_round), by =4),0))+
    theme(legend.title = element_blank())
  
  tiff(paste0('~/PROJECTS/WRAP Location/Manuscript/Figures/Spatial Correlation Plots/HMS_Spatial_Correlation_allESM_',m,'.tiff'),
       units="in",width=13, height=10, res=200)
  plot(g1)
  dev.off()
  
  tiff(paste0('~/PROJECTS/WRAP Location/Manuscript/Figures/Spatial Correlation Plots/CPS_Spatial_Correlation_allESM_',m,'.tiff'),
       units="in",width=8, height=10, res=200)
  plot(g2)
  dev.off()
  
  tiff(paste0('~/PROJECTS/WRAP Location/Manuscript/Figures/Spatial Correlation Plots/GFS_Spatial_Correlation_allESM_',m,'.tiff'),
       units="in",width=8, height=10, res=200)
  plot(g3)
  dev.off()
}



#----FIGURE 4.7: Maps of RMSE------
#Map of correlations averaged over ~10years
#Load in predictions made in EM files
species_rmse_spatial <- list()
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
  #Allocate decade
  all_esm_allyears$decade <- ifelse(all_esm_allyears$year<=2010,"1985-2010",
                                    ifelse(all_esm_allyears$year>=2011 & all_esm_allyears$year<=2040, "2011-2040",
                                           ifelse(all_esm_allyears$year>=2041 & all_esm_allyears$year<=2070, "2041-2070",
                                                  ifelse(all_esm_allyears$year>=2071, "2071-2100", NA))))
  all_esm_allyears$lat_round <- round(all_esm_allyears$lat,0)
  all_esm_allyears$lon_round <- round(all_esm_allyears$lon,0)
  
  #Calculate correlation
  RMSE = function(p, o){(sqrt(mean((p - o)^2))) / mean(o)}
  #aggregate by ESM and grid cell
  all_esm_agg <-all_esm_allyears %>% group_by(lon_round, lat_round, decade, esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                                                    "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                                                    "brt_E" , "brt_ES", "brt_EST" ,  
                                                                                                    "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(RMSE(., o=abundance)))
  all_esm_agg_longer <- melt(all_esm_agg,id=c("lat_round","lon_round","decade", "esm"), variable.name = "EM",value.name = "RMSE")
  all_esm_agg_longer$RMSE <- as.numeric(ifelse(all_esm_agg_longer$RMSE=="Inf", "NA",all_esm_agg_longer$RMSE))
  # all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  species_rmse_spatial[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

models <- c( "gam_E","gam_ES", "gam_ECor",
             "glm_E" ,"glm_ESt", "glm_ESr",   
             "brt_E" , "brt_ES", "brt_EST" ,  
             "mlp_E" , "mlp_ES" ,"mlp_EST","ens_mean")
for (m in models){
  print(m)
  dat <- species_rmse_spatial[[1]]
  g1 <- ggplot(dat[dat$EM==m,],aes(x=lon_round,y=lat_round))+
    geom_tile(aes(fill=log(RMSE)))+
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',)+
    theme_classic()+
    annotation_map(map_data("world"), colour = "black", fill=grey(0.7))+
    coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48))+
    facet_wrap(~esm + decade, nrow=3)+
    labs(x="",y="")+
    theme(legend.title = element_blank())
  
  dat <- species_rmse_spatial[[2]]
  g2 <- ggplot(dat[dat$EM==m,],aes(x=lon_round,y=lat_round))+
    geom_tile(aes(fill=log(RMSE)))+
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',)+
    theme_classic()+
    annotation_map(map_data("world"), colour = "black", fill=grey(0.7))+
    coord_quickmap(xlim=c(-127,-115.8),ylim=c(30,48))+
    facet_wrap(~esm + decade, nrow=3)+
    scale_x_continuous(breaks = round(seq(min(dat$lon_round), max(dat$lon_round), by =4),0))+
    theme(legend.title = element_blank())
  
  dat <- species_rmse_spatial[[3]]
  g3 <- ggplot(dat[dat$EM==m,],aes(x=lon_round,y=lat_round))+
    geom_tile(aes(fill=log(RMSE)))+
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red',)+
    theme_classic()+
    annotation_map(map_data("world"), colour = "black", fill=grey(0.7))+
    coord_quickmap(xlim=c(-127,-115.8),ylim=c(30,48))+
    facet_wrap(~esm + decade, nrow=3)+
    scale_x_continuous(breaks = round(seq(min(dat$lon_round), max(dat$lon_round), by =4),0))+
    theme(legend.title = element_blank())
  
  tiff(paste0('~/PROJECTS/WRAP Location/Manuscript/Figures/Spatial RMSE Plots/HMS_Spatial_RMSE_allESM_',m,'.tiff'),
       units="in",width=13, height=10, res=200)
  plot(g1)
  dev.off()
  
  tiff(paste0('~/PROJECTS/WRAP Location/Manuscript/Figures/Spatial RMSE Plots/CPS_Spatial_RMSE_allESM_',m,'.tiff'),
       units="in",width=8, height=10, res=200)
  plot(g2)
  dev.off()
  
  tiff(paste0('~/PROJECTS/WRAP Location/Manuscript/Figures/Spatial RMSE Plots/GFS_Spatial_RMSE_allESM_',m,'.tiff'),
       units="in",width=8, height=10, res=200)
  plot(g3)
  dev.off()
}
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


#----FIGURE 8: Fitted years Correlation (also RMSE)-------
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

#----FIGURE 9: Forecast RMSE of good models ----
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

#----FIGURE 10: environmental density plots-----

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

#----FIGURE 12: Partition variance more closely: BIOMASS-----
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
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ),c(sum),na.rm=T)
  #Move EMs to long format
  all_esm_agg_long <- all_esm_agg %>% pivot_longer(cols=c("gam_E", "gam_ES", "gam_ECor",
                                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                                          "mlp_E" , "mlp_ES" ,"mlp_EST"), names_to="EM", values_to="Biomass")
  
  #get SD across EMs
  EM_SD <-all_esm_agg_long %>% group_by(year, esm) %>% summarise_at(c("Biomass" ),sd)
  
  #get SD across EMs AND ESMs
  EM_ESM_SD <-all_esm_agg_long %>% group_by(year) %>% summarise_at(c("Biomass" ),sd) 
  
  #get SD across ESMs (using observed biomass)
  #First move to long format
  ESM_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("abundance"),c(sum),na.rm=T)
  ESM_SD <-ESM_agg %>% group_by(year) %>% summarise_at(c("abundance" ),sd) 
  
  #merge three error sources
  error <- left_join(EM_ESM_SD,ESM_SD,by="year")
  colnames(error) <- c("year","biomass_error_ESM_EM","biomass_error_ESM")
  error_all <- left_join(error,EM_SD, by="year")
  colnames(error_all) <- c("year","biomass_error_ESM_EM","biomass_error_ESM", "esm", "biomass_error_EM")
  
  variance_partition_summary[[counter]]<- error_all
  counter=counter+1
}



ggplot(data = variance_partition_summary[[3]], aes(x=year))+
  geom_line(aes(y=biomass_error_ESM_EM),col="red")+
  geom_line(aes(y=biomass_error_ESM),col="blue")+
  geom_line(aes(y=biomass_error_EM, group=esm),col="green")
  


#OLD CODE
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
# combined$year[which(combined$value[combined$error_source=="sd_esm"]>hline)]


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





#----Figure 13: within model uncertainty-----
#Need to plot within model uncertainty, and hopefully it's less than between model uncertainty
#Try first for a GAM for Albacore
#read in dat files and model objects
withinmodel_uncertainty <- list()
species_counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  
  #Load in data & EMS and simulate from posterior distribution
  counter=1
  esm_sims <- as.data.frame(matrix(NA,nrow=1044, ncol=6))
  colnames(esm_sims) <- c("year","esm","model","mean","min","max")
  for (esm in c("had","ipsl","gfdl")){
    
    dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/",esm,"/dat_hist_results_full.rds"))
    dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/",esm,"/dat_fcast_results_full.rds"))
    dat_allyears <- rbind(dat_hist,dat_fcast)
    load(paste0("~/PROJECTS/WRAP Location/",s,"/",esm,"/saved_models_full.RData"))
    
    #remove unnecessary models (at this stage)
    rm( glm_E_N,glm_E_P, glm_ESr_N,glm_ESr_P, glm_ESt_N, glm_ESt_P, glm_Sr_N, glm_Sr_P,
        brt_E_N, brt_E_P, brt_ES_N, brt_ES_P, brt_EST_N, brt_EST_P,
        mlp_E_N, mlp_E_P, mlp_ES_N, mlp_ES_P, mlp_EST_N, mlp_EST_P)

    #Draw fitted values from posterior distribution of GAM
    N1 <- fitted_samples(gam_E_N, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_E_P, n=100, newdata=dat_allyears, scale="response", seed=99)  #response scale doesn't work for binomial model for some reason
    sim_data_E <- cbind(N1,P1$fitted)
    colnames(sim_data_E) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_E$biomass  <- exp(sim_data_E$bio_sim) * sim_data_E$pres_sim
    sim_data_E$year <- dat_allyears$year[sim_data_E$inputdata_row]
    sim_data_E$model <- "GAM_E"
    
    N1 <- fitted_samples(gam_ES_N, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_ES_P, n=100, newdata=dat_allyears, scale="response", seed=99)  #response scale doesn't work for binomial model for some reason
    sim_data_ES <- cbind(N1,P1$fitted)
    colnames(sim_data_ES) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_ES$biomass  <- exp(sim_data_ES$bio_sim) * sim_data_ES$pres_sim
    sim_data_ES$year <- dat_allyears$year[sim_data_ES$inputdata_row]
    sim_data_ES$model <- "GAM_ES"
    
    N1 <- fitted_samples(gam_ECor_N$gam, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_ECor_P$gam, n=100, newdata=dat_allyears, scale="response", seed=99)  #response scale doesn't work for binomial model for some reason
    sim_data_ECor <- cbind(N1,P1$fitted)
    colnames(sim_data_ECor) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_ECor$biomass  <- exp(sim_data_ECor$bio_sim) * sim_data_ECor$pres_sim
    sim_data_ECor$year <- dat_allyears$year[sim_data_ECor$inputdata_row]
    sim_data_ECor$model <- "GAM_ECor"
    
    sim_data_esm <- rbind(sim_data_E,sim_data_ES,sim_data_ECor)
    sim_data_esm$esm <- esm
    
    #Get sum of simulated biomass across grid cells
    sim_data_agg <-  sim_data_esm %>% group_by(year, draw, esm, model) %>% summarise_at(c("biomass"), funs(sum)) 
    #Get mean, min, and max of sims
    sim_data_agg_funs <-  sim_data_agg %>% group_by(year, esm, model) %>% summarise_at(c("biomass"), funs(mean, min, max))
    
    sc <- counter 
    ec <- sc + 1043
    esm_sims[sc:ec,] <- sim_data_agg_funs
    counter = counter + 1044
  }
  
  withinmodel_uncertainty[[species_counter]] <- esm_sims
  species_counter <- species_counter +1
  
}

  
#Now get true abundance for plotting
obs_dat <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
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
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("had","gfdl","ipsl"),each=58000)
  all_esm_allyears_agg <- all_esm_allyears %>% group_by(year, esm) %>% summarise_at(c("abundance"),sum)  
  obs_dat[[counter]] <- all_esm_allyears_agg
  counter=counter+1
}

hms <- ggplot(data=withinmodel_uncertainty[[1]], aes(x=year,y=mean))+
  geom_line(data=obs_dat[[1]], aes(x=year,y=abundance), col="grey")+
  geom_ribbon(aes(ymin=min,ymax=max), alpha=0.5, fill="red")+
  geom_line(aes(y=mean), col="red")+
  theme_classic()+
  facet_wrap(~esm + model)+
  labs(x="",y="HMS Biomass")+
  geom_vline(xintercept = 2010, linetype="dashed")
  
cps <- ggplot(data=withinmodel_uncertainty[[2]], aes(x=year,y=mean))+
  geom_line(data=obs_dat[[2]], aes(x=year,y=abundance), col="grey")+
  geom_ribbon(aes(ymin=min,ymax=max), alpha=0.5, fill="red")+
  geom_line(aes(y=mean), col="red")+
  theme_classic()+
  facet_wrap(~esm + model)+
  labs(x="",y="CPS Biomass")+
  geom_vline(xintercept = 2010, linetype="dashed")

gfs <- ggplot(data=withinmodel_uncertainty[[3]], aes(x=year,y=mean))+
  geom_line(data=obs_dat[[3]], aes(x=year,y=abundance), col="grey")+
  geom_ribbon(aes(ymin=min,ymax=max), alpha=0.5, fill="red")+
  geom_line(aes(y=mean), col="red")+
  theme_classic()+
  facet_wrap(~esm + model)+
  labs(x="",y="Groundfish Biomass")+
  geom_vline(xintercept = 2010, linetype="dashed")

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig13_hms_n100.tiff',res=300,units="in",height=10,width=13)
plot(hms)
dev.off()

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig13_cps_n100.tiff',res=300,units="in",height=10,width=13)
plot(cps)
dev.off()

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig13_gfs_n100.tiff',res=300,units="in",height=10,width=13)
plot(gfs)
dev.off()



# #TEMP: albacore, had, gam_e only
# dtemp <- withinmodel_uncertainty[[1]]
# dtemp2 <- dtemp[dtemp$model=="GAM_E",]
# ggplot(data=dtemp[dtemp$model=="GAM_E",], aes(x=year,y=mean))+
#   geom_line(data=obs_dat[[1]], aes(x=year,y=abundance), col="grey")+
#   geom_ribbon(aes(ymin=min,ymax=max), alpha=0.5, fill="red")+
#   geom_line(aes(y=mean), col="red")+
#   theme_classic()+
#   # facet_wrap(~esm + model)+
#   labs(x="",y="HMS Biomass")+
#   geom_vline(xintercept = 2010, linetype="dashed")


#TEMPORARY
# gam_n_pred <- as.data.frame(predict.gam(gam_E_N,newdata=dat_allyears,type="response", se.fit=TRUE ))
# gam_p_pred <- as.data.frame(predict.gam(gam_E_P,newdata=dat_allyears,type="response", se.fit=TRUE ))
# 
# gam_preds_mean <- as.data.frame((exp(gam_n_pred$fit)) * (gam_p_pred$fit))
# gam_preds_high <- as.data.frame((exp(gam_n_pred$fit + (gam_n_pred$se.fit))) * (gam_p_pred$fit + (gam_p_pred$se.fit)))
# gam_preds_low <- as.data.frame((exp(gam_n_pred$fit - (gam_n_pred$se.fit))) * (gam_p_pred$fit - (gam_p_pred$se.fit)))
# gam_preds <- cbind(gam_preds_mean, gam_preds_high, gam_preds_low)
# colnames(gam_preds) <- c("mean","high","low")
# gam_preds$year <- rep(1985:2100,each=500)
# gam_preds_agg <-  gam_preds %>% group_by(year) %>% summarise_at(c("mean","high","low"), funs(sum)) 
# ggplot(data=gam_preds_agg, aes(x=year))+
#   geom_line(aes(x=year,y=mean), col="grey")+
#   geom_ribbon(aes(ymin=low,ymax=high), alpha=0.5, fill="red")
# geom_ribbon(aes(ymin=min,ymax=max), alpha=0.5, fill="red")+


#-----FIGURE 14: delta maps------
#2x3 plot: 3 species archetypes of delta distribution maps of historical and 2100 period
#as per Fig 7 in Morley et el
#warning: this takes a few mins to plot (Should have used Faet)

delta_output <- list()
counter=1
for (s in c("hms","cps","gfs")){
  #Get data
  if (s =="hms"){
    dat_had <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","had","/dat_albacore_all.rds"))
    dat_gfdl <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","gfdl","/dat_albacore_all.rds"))
    dat_ipsl <- readRDS(paste0("~/PROJECTS/WRAP Location/Albacore EMs/","ipsl","/dat_albacore_all.rds"))
  }
  if(s=="cps"){
    dat_had <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/had/dat_anchovy_all.rds"))
    dat_gfdl <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/gfdl/dat_anchovy_all.rds"))
    dat_ipsl <- readRDS(paste0("~/PROJECTS/WRAP Location/Anchovy EMs/ipsl/dat_anchovy_all.rds"))
  }
  if(s=="gfs"){
    dat_had <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/had/dat_groundfish_all.rds"))
    dat_gfdl <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/had/dat_groundfish_all.rds"))
    dat_ipsl <- readRDS(paste0("~/PROJECTS/WRAP Location/Groundfish EMs/had/dat_groundfish_all.rds"))
  }
  #combine esm
  dat_had$esm <- "had"
  dat_gfdl$esm <- "gfdl"
  dat_ipsl$esm <- "ipsl"
  dat_allesm <- rbind(dat_had, dat_gfdl, dat_ipsl)
  #isolate years of interest
  dat_allesm_hist <- dat_allesm[dat_allesm$year<=2010,]
  dat_allesm_mid <- dat_allesm[dat_allesm$year>=2025 & dat_allesm$year<=2050,]
  dat_allesm_fcast <- dat_allesm[dat_allesm$year>=2075,]
  #average across years and ensembles
  dat_allesm_hist_agg <- dat_allesm_hist %>% group_by(lat, lon) %>% summarise_at(c("abundance"),c(mean),na.rm=T)
  dat_allesm_mid_agg <- dat_allesm_mid %>% group_by(lat, lon) %>% summarise_at(c("abundance"),c(mean),na.rm=T)
  dat_allesm_fcast_agg <- dat_allesm_fcast %>% group_by(lat, lon) %>% summarise_at(c("abundance"),c(mean),na.rm=T)
  #now merge and get delta: distant future
  dat_allesm_agg_delta <- left_join(dat_allesm_hist_agg,dat_allesm_fcast_agg ,by=c("lat","lon"))
  dat_allesm_agg_delta$delta <- dat_allesm_agg_delta$abundance.y - dat_allesm_agg_delta$abundance.x  
    #now merge and get delta: mid future
  dat_allesm_agg_delta_mid <- left_join(dat_allesm_hist_agg,dat_allesm_mid_agg ,by=c("lat","lon"))
  dat_allesm_agg_delta_mid$delta <- dat_allesm_agg_delta_mid$abundance.y - dat_allesm_agg_delta_mid$abundance.x  
  
  #saveoutput
  delta_output[[counter]] <- dat_allesm_agg_delta
  delta_output[[counter+1]] <- dat_allesm_hist_agg
  delta_output[[counter+2]] <- dat_allesm_agg_delta_mid
  counter = counter+3
}

#now plot
hms_fcast <- ggplot(data = delta_output[[1]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=delta))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.78,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#e74c3c", "white", "#3498db"))(256),
                       limits = c(min(delta_output[[1]]$delta), max(delta_output[[1]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="HMS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2075-2100", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass", size=3)

hms_mid <- ggplot(data = delta_output[[3]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=delta))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.78,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#e74c3c", "white", "#3498db"))(256),
                       limits = c(min(delta_output[[3]]$delta), max(delta_output[[3]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="HMS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2035-2060", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass", size=3)

hms_hist <- ggplot(data = delta_output[[2]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=abundance))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.78,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c"))(256),
                       limits = c(min(delta_output[[2]]$abundance), max(delta_output[[2]]$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="HMS Archetype", size=5)+
  geom_text(x=-120, y=44, label="1985-2010", size=3)+
  geom_text(x=-119.5, y=42.5, label="Biomass", size=3)

#CPS
cps_fcast <- ggplot(data = delta_output[[4]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=delta))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#e74c3c", "white", "#3498db"))(256),
                       limits = c(min(delta_output[[4]]$delta), max(delta_output[[4]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="CPS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2075-2100", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass", size=3)

cps_mid <- ggplot(data = delta_output[[6]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=delta))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#e74c3c", "white", "#3498db"))(256),
                       limits = c(min(delta_output[[6]]$delta), max(delta_output[[6]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="CPS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2035-2060", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass", size=3)

cps_hist <- ggplot(data = delta_output[[5]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=abundance))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c"))(256),
                       limits = c(min(delta_output[[5]]$abundance), max(delta_output[[5]]$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="CPS Archetype", size=5)+
  geom_text(x=-120, y=44, label="1985-2010", size=3)+
  geom_text(x=-119.5, y=42.5, label="Biomass", size=3)

#REPEAT FOR GROUNDFISH
gfs_fcast <- ggplot(data = delta_output[[7]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=delta))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#e74c3c", "white", "#3498db"))(256),
                       limits = c(min(delta_output[[7]]$delta), max(delta_output[[7]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="GFS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2075-2100", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass", size=3)

gfs_mid <- ggplot(data = delta_output[[9]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=delta))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#e74c3c", "white", "#3498db"))(256),
                       limits = c(min(delta_output[[9]]$delta), max(delta_output[[9]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="GFS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2035-2060", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass", size=3)

gfs_hist <- ggplot(data = delta_output[[8]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=abundance))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c"))(256),
                       limits = c(min(delta_output[[8]]$abundance), max(delta_output[[8]]$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="GFS Archetype", size=5)+
  geom_text(x=-120, y=44, label="1985-2010", size=3)+
  geom_text(x=-119.5, y=42.5, label="Biomass", size=3)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig14.tiff',res=300,units="in",height=10,width=12)
grid.arrange(hms_hist, cps_hist, gfs_hist,
             # hms_mid, cps_mid, gfs_mid,
             hms_fcast, cps_fcast, gfs_fcast, nrow=2)
dev.off()

#-----FIGURE 15: Partitioning uncertainty with dominance analysis and COG-----
#Morley et al. used a GLM to quantify uncertainty, where:
# y ~ RCP + ESM + SDM + (RCP * ESM * SDM), where y is the COG. 
#Partitioned sum of squares from resultant model...what does that mean? 
# Then used a dominance analysis to quantify relative importance of each variable

#Let's try it:
#Get COG for all species and EMs
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
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  
  #aggregate over space
  all_esm_agg <- all_esm_allyears %>% group_by(year,esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                               "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                               "brt_E" , "brt_ES", "brt_EST" ,  
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=lat)),na.rm=T)
  all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(3:14),names_to="EM",values_to="COG")
  species_cog_diff[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

#Now run dominance analysis for every year & archetype
da_allyears_allspecies <- list()
upper_counter <- 1
for (s in 1:3){
  dat_cog <- species_cog_diff[[s]]
  dat_cog <- dat_cog %>%  separate(EM, c("type","params"))
  dat_cog$esm <- as.factor(dat_cog$esm)
  dat_cog$type <- as.factor(dat_cog$type)
  dat_cog$params <- as.factor(dat_cog$params)

  #Anova to look at where significant differences lie:
  # m1 <- aov(COG ~ esm + type + params, data = dat_cog) #must be aov to run TukeyHSD
  # summary(m1)
  # TukeyHSD(m1)
  
  #Run dominance analysis for every year
  da_allyears <- as.data.frame(matrix(NA, nrow=348, ncol=4))
  colnames(da_allyears) <- c("r2", "source", "year", "percent")
  counter=1
  for (y in 1985:2100){
    print(y)
    m2 <- lm(COG ~ esm + type + params, data = dat_cog[dat_cog$year==y,]) #must be class lm for dominance analysis
    da <- dominanceAnalysis(m2)
    da_df <- as.data.frame(da$contribution.average)
    da_df$source <- row.names(da_df)
    da_df$year <- y
    da_df$percent <- 1 - (sum(da_df$r2) - da_df$r2)/ sum(da_df$r2)
    rownames(da_df) <- NULL
    sc <- counter
    ec <- counter + 2
    da_allyears[sc:ec, ] <- da_df
    counter=counter+3
  }
  da_allyears_allspecies[[upper_counter]] <- da_allyears
  upper_counter <- upper_counter+1
}
summary(da_allyears_allspecies)

#Find first year when esm < 50%
i1 <- da_allyears_allspecies[[1]]$year[da_allyears_allspecies[[1]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]
i2 <- da_allyears_allspecies[[2]]$year[da_allyears_allspecies[[2]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]
i3 <- da_allyears_allspecies[[3]]$year[da_allyears_allspecies[[3]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]

#Now plot
hms <- ggplot(data=da_allyears_allspecies[[1]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="HMS Archetype", size=5)
  
cps <- ggplot(data=da_allyears_allspecies[[2]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i2[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="CPS Archetype", size=5)

gfs <- ggplot(data=da_allyears_allspecies[[3]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i3[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="GFS Archetype", size=5)


tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig15.tiff', units="in", res=300, height=12, width=10)
grid.arrange(hms, cps, gfs, ncol=1)
dev.off()


#-----FIGURE 16: dominance analysis by region---------
# regions: South 30-34.5; north 40-48; central 34.5-40
#For biomass not COG though. Make it relative to historical period? 

#Let's try it: ADD REGION TO THIS STEPH
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
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  all_esm_allyears$region <- ifelse(all_esm_allyears$lat<=34.5, "South",
                                    ifelse(all_esm_allyears$lat>=40, "North", "Central"))
  #aggregate over space
  all_esm_agg <- all_esm_allyears %>% group_by(year,esm,region) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                               "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                               "brt_E" , "brt_ES", "brt_EST" ,  
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(sum(.)),na.rm=T)
  all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(4:15),names_to="EM",values_to="TotalBiomass")
  species_biomass[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

#Now run dominance analysis for every year & archetype
da_allyears_allspecies <- list()
upper_counter <- 1
for (s in 1:3){
  dat_bio <- species_biomass[[s]]
  dat_bio <- dat_bio %>%  separate(EM, c("type","params"))
  dat_bio$esm <- as.factor(dat_bio$esm)
  dat_bio$type <- as.factor(dat_bio$type)
  dat_bio$params <- as.factor(dat_bio$params)
  
  #Anova to look at where significant differences lie:
  # m1 <- aov(COG ~ esm + type + params, data = dat_cog) #must be aov to run TukeyHSD
  # summary(m1)
  # TukeyHSD(m1)
  
  for (r in c("South", "Central","North")){
    #Run dominance analysis for every year
    da_allyears <- as.data.frame(matrix(NA, nrow=348, ncol=6))
    colnames(da_allyears) <- c("r2", "source", "year", "percent", "region", "species")
    counter=1
    for (y in 1985:2100){
      print(y)
      m2 <- lm(TotalBiomass ~ esm + type + params, data = dat_bio[dat_bio$year==y & dat_bio$region==r,]) #must be class lm for dominance analysis
      da <- dominanceAnalysis(m2)
      da_df <- as.data.frame(da$contribution.average)
      da_df$source <- row.names(da_df)
      da_df$year <- y
      da_df$percent <- 1 - (sum(da_df$r2) - da_df$r2)/ sum(da_df$r2)
      da_df$region <- r
      da_df$species <- s
      rownames(da_df) <- NULL
      sc <- counter
      ec <- counter + 2
      da_allyears[sc:ec, ] <- da_df
      counter=counter+3
    }
    da_allyears_allspecies[[upper_counter]] <- da_allyears
    upper_counter <- upper_counter+1
  }
}
summary(da_allyears_allspecies)

#Find first year when esm < 50%

#Now plot

dat <- do.call(rbind, da_allyears_allspecies)
dat$species <- ifelse(dat$species==1,"HMS Archetype",
                      ifelse(dat$species==2,"CPS Archetype", "GFS Archetype"))
dat$species <- factor(dat$species, levels = c("HMS Archetype", "CPS Archetype", "GFS Archetype"))
dat$region <- factor(dat$region, levels = c("North","Central","South"))

g1 <- ggplot(data=dat, aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Biomass Uncertainty") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  facet_wrap(~region+species)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig16.tiff', units="in", res=300, height=9, width=12)
plot(g1)
dev.off()

#-----Figure 17: r2 correlation with temp only models------
#Load in predictions made in EM files
species_rmse_diff <- list()
counter=1
for(s in c("Albacore EMs TempOnly","Anchovy EMs TempOnly", "Groundfish EMs TempOnly")){ #add GFS when ready
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
  all_esm_agg <- all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(COR(., y=abundance)))
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

hms <- ggplot(data=species_rmse_diff[[1]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

cps <- ggplot(data=species_rmse_diff[[2]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

gfs <- ggplot(data=species_rmse_diff[[3]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="GFS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig17.tiff', units="in", res=300, height=12, width=10)
grid.arrange(hms, cps, gfs, ncol=1)
dev.off()


#-----Figure 17.5: r2 correlation with 2040 trained models------
#Load in predictions made in EM files
species_rmse_diff <- list()
counter=1
for(s in c("Albacore EMs 2040train","Anchovy EMs 2040train", "Groundfish EMs 2040train")){ #add GFS when ready
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

hms <- ggplot(data=species_rmse_diff[[1]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2040, linetype="dashed")+
  facet_wrap(~esm)

cps <- ggplot(data=species_rmse_diff[[2]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2040, linetype="dashed")+
  facet_wrap(~esm)

gfs <- ggplot(data=species_rmse_diff[[3]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="GFS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2040, linetype="dashed")+
  facet_wrap(~esm)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig17.5.tiff', units="in", res=300, height=12, width=10)
grid.arrange(hms, cps, gfs, ncol=1)
dev.off()


#-----Figure 18: compare uncertainty with temp only models--------
#Let's try it:
#Get COG for all species and EMs
species_cog_diff <- list()
counter=1
for(s in c("Albacore EMs TempOnly","Anchovy EMs TempOnly","Groundfish EMs TempOnly")){
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
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  
  #aggregate over space
  all_esm_agg <- all_esm_allyears %>% group_by(year,esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                               "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                               "brt_E" , "brt_ES", "brt_EST" ,  
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=lat)),na.rm=T)
  all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(3:14),names_to="EM",values_to="COG")
  species_cog_diff[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

#Now run dominance analysis for every year & archetype
da_allyears_allspecies <- list()
upper_counter <- 1
for (s in 1:3){
  dat_cog <- species_cog_diff[[s]]
  dat_cog <- dat_cog %>%  separate(EM, c("type","params"))
  dat_cog$esm <- as.factor(dat_cog$esm)
  dat_cog$type <- as.factor(dat_cog$type)
  dat_cog$params <- as.factor(dat_cog$params)
  
  #Anova to look at where significant differences lie:
  # m1 <- aov(COG ~ esm + type + params, data = dat_cog) #must be aov to run TukeyHSD
  # summary(m1)
  # TukeyHSD(m1)
  
  #Run dominance analysis for every year
  da_allyears <- as.data.frame(matrix(NA, nrow=348, ncol=4))
  colnames(da_allyears) <- c("r2", "source", "year", "percent")
  counter=1
  for (y in 1985:2100){
    print(y)
    m2 <- lm(COG ~ esm + type + params, data = dat_cog[dat_cog$year==y,]) #must be class lm for dominance analysis
    da <- dominanceAnalysis(m2)
    da_df <- as.data.frame(da$contribution.average)
    da_df$source <- row.names(da_df)
    da_df$year <- y
    da_df$percent <- 1 - (sum(da_df$r2) - da_df$r2)/ sum(da_df$r2)
    rownames(da_df) <- NULL
    sc <- counter
    ec <- counter + 2
    da_allyears[sc:ec, ] <- da_df
    counter=counter+3
  }
  da_allyears_allspecies[[upper_counter]] <- da_allyears
  upper_counter <- upper_counter+1
}
summary(da_allyears_allspecies)

#Find first year when esm < 50%
# i1 <- da_allyears_allspecies[[1]]$year[da_allyears_allspecies[[1]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]
# i2 <- da_allyears_allspecies[[2]]$year[da_allyears_allspecies[[2]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]
# i3 <- da_allyears_allspecies[[3]]$year[da_allyears_allspecies[[3]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]

#Now plot
hms <- ggplot(data=da_allyears_allspecies[[1]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="HMS Archetype", size=5)

cps <- ggplot(data=da_allyears_allspecies[[2]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="CPS Archetype", size=5)

gfs <- ggplot(data=da_allyears_allspecies[[3]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="GFS Archetype", size=5)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig18.tiff', units="in", res=300, height=12, width=10)
grid.arrange(hms, cps, gfs, ncol=1)
dev.off()

#-----Figure 19: compare uncertainty with 2040 trained models---------
#Let's try it:
species_cog_diff <- list()
counter=1
for(s in c("Albacore EMs 2040train","Anchovy EMs 2040train","Groundfish EMs 2040train")){
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
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=58000)
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  
  #aggregate over space
  all_esm_agg <- all_esm_allyears %>% group_by(year,esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                               "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                               "brt_E" , "brt_ES", "brt_EST" ,  
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=lat)),na.rm=T)
  all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(3:14),names_to="EM",values_to="COG")
  species_cog_diff[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

#Now run dominance analysis for every year & archetype
da_allyears_allspecies <- list()
upper_counter <- 1
for (s in 1:3){
  dat_cog <- species_cog_diff[[s]]
  dat_cog <- dat_cog %>%  separate(EM, c("type","params"))
  dat_cog$esm <- as.factor(dat_cog$esm)
  dat_cog$type <- as.factor(dat_cog$type)
  dat_cog$params <- as.factor(dat_cog$params)
  
  #Anova to look at where significant differences lie:
  # m1 <- aov(COG ~ esm + type + params, data = dat_cog) #must be aov to run TukeyHSD
  # summary(m1)
  # TukeyHSD(m1)
  
  #Run dominance analysis for every year
  da_allyears <- as.data.frame(matrix(NA, nrow=348, ncol=4))
  colnames(da_allyears) <- c("r2", "source", "year", "percent")
  counter=1
  for (y in 1985:2100){
    print(y)
    m2 <- lm(COG ~ esm + type + params, data = dat_cog[dat_cog$year==y,]) #must be class lm for dominance analysis
    da <- dominanceAnalysis(m2)
    da_df <- as.data.frame(da$contribution.average)
    da_df$source <- row.names(da_df)
    da_df$year <- y
    da_df$percent <- 1 - (sum(da_df$r2) - da_df$r2)/ sum(da_df$r2)
    rownames(da_df) <- NULL
    sc <- counter
    ec <- counter + 2
    da_allyears[sc:ec, ] <- da_df
    counter=counter+3
  }
  da_allyears_allspecies[[upper_counter]] <- da_allyears
  upper_counter <- upper_counter+1
}
summary(da_allyears_allspecies)

#Find first year when esm < 50%
# i1 <- da_allyears_allspecies[[1]]$year[da_allyears_allspecies[[1]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]
# i2 <- da_allyears_allspecies[[2]]$year[da_allyears_allspecies[[2]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]
# i3 <- da_allyears_allspecies[[3]]$year[da_allyears_allspecies[[3]]$percent<0.5 & da_allyears_allspecies[[1]]$source=="esm"]

#Now plot
hms <- ggplot(data=da_allyears_allspecies[[1]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2040, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="HMS Archetype", size=5)

cps <- ggplot(data=da_allyears_allspecies[[2]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2040, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="CPS Archetype", size=5)

gfs <- ggplot(data=da_allyears_allspecies[[3]], aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Centroid Uncertainty") + 
  geom_vline(xintercept = 2040, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  geom_text(x=2080, y=0.95, label="GFS Archetype", size=5)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig19.tiff', units="in", res=300, height=12, width=10)
grid.arrange(hms, cps, gfs, ncol=1)
dev.off()

#------FIGURE 20: correlation between R2 and ESM relative contribution------

#PART ONE:
#Get R2 for each model, species, ESM, and experiment (n=342)
cors_all <- as.data.frame(matrix(NA, nrow=9, ncol=3))
colnames(cors_all) <- c("species","experiment","forecast_R2")
counter=1
for(s in c("Albacore EMs", "Albacore EMs TempOnly", "Albacore EMs 2040train",
           "Anchovy EMs", "Anchovy EMs TempOnly","Anchovy EMs 2040train",
           "Groundfish EMs", "Groundfish EMs TempOnly", "Groundfish EMs 2040train")){
  print(s)
  #Only get forecast years 
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_dat_fcast,gfdl_dat_fcast,ipsl_dat_fcast)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=nrow(had_dat_fcast))
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","abundance")]
  #Get cor between observed and predicted for each model
  all_esm_agg <- all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                                                          "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(cor(., x=abundance)))
  #Get mean across EM:
  all_esm_agg$mean_r2 <- apply(all_esm_agg[,2:13],1,mean)
  
  # #Pivot longer
  # all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(2:13),names_to="EM",values_to="R2")
  
  #Make species & experiment names
  sp <- unlist(strsplit(s," "))[1]
  exp <- unlist(strsplit(s," "))[3]
  exp <- ifelse(is.na(exp),"all",exp)
  
  #Write out data
  cors_all[counter,1] <- sp
  cors_all[counter,2] <- exp
  cors_all[counter,3] <- mean(all_esm_agg$mean_r2)
  counter=counter+1
}
head(cors_all)
tail(cors_all)

#PART TWO: 
#Get %contribution of ESM via dominance analysis (first get COG)

dom_all <- as.data.frame(matrix(NA, nrow=9, ncol=3))
colnames(dom_all) <- c("species","experiment","esm_contribution")
counter=1
for(s in c("Albacore EMs", "Albacore EMs TempOnly", "Albacore EMs 2040train",
           "Anchovy EMs", "Anchovy EMs TempOnly","Anchovy EMs 2040train",
           "Groundfish EMs", "Groundfish EMs TempOnly", "Groundfish EMs 2040train")){
  print(s)
  #Forecast years only
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_dat_fcast,gfdl_dat_fcast,ipsl_dat_fcast)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=nrow(ipsl_dat_fcast))
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  
  #aggregate over space
  all_esm_agg <- all_esm_allyears %>% group_by(year,esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                               "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                               "brt_E" , "brt_ES", "brt_EST" ,  
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=lat)),na.rm=T)
  all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(3:14),names_to="EM",values_to="COG")
  
  #Make species & experiment names
  sp <- unlist(strsplit(s," "))[1]
  exp <- unlist(strsplit(s," "))[3]
  exp <- ifelse(is.na(exp),"all",exp)
  
  #Run dominance analysis across all years
  all_esm_agg_longer <- all_esm_agg_longer %>%  separate(EM, c("type","params"))
  m2 <- lm(COG ~ esm + type + params, data = all_esm_agg_longer) #must be class lm for dominance analysis
  da <- dominanceAnalysis(m2)
  da_df <- as.data.frame(da$contribution.average)
  da_df$source <- row.names(da_df)
  da_df$percent <- 1 - (sum(da_df$r2) - da_df$r2)/ sum(da_df$r2)
  
  dom_all[counter,1] <- sp
  dom_all[counter,2] <- exp
  dom_all[counter,3] <- da_df$percent[da_df$source=="esm"]
  counter = counter + 1
}


# PART THREE:
#plot
plot(cors_all$forecast_R2,dom_all$esm_contribution)


#------FIGURE 21: barplots of r2 for each ESM, Species, and experiment--------
#Get R2 for each model, species, ESM, and experiment (n=342)
cors_all <- as.data.frame(matrix(NA, nrow=9, ncol=4))
colnames(cors_all) <- c("species","experiment","forecast_mean_R2", "forecast_sd_R2")
counter=1
for(s in c("Albacore EMs", "Albacore EMs TempOnly", "Albacore EMs 2040train",
           "Anchovy EMs", "Anchovy EMs TempOnly","Anchovy EMs 2040train",
           "Groundfish EMs", "Groundfish EMs TempOnly", "Groundfish EMs 2040train")){
  print(s)
  #Only get forecast years 
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_dat_fcast,gfdl_dat_fcast,ipsl_dat_fcast)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=nrow(had_dat_fcast))
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","abundance")]
  #Get cor between observed and predicted for each model
  all_esm_agg <- all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                                                          "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(cor(., x=abundance, method="spearman")))
  #Get mean across EM:
  mean_r2 <- mean(as.matrix(all_esm_agg[,2:13]))
  sd_r2 <- sd(as.matrix(all_esm_agg[,2:13]))
  
  # #Pivot longer
  # all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(2:13),names_to="EM",values_to="R2")
  
  #Make species & experiment names
  sp <- unlist(strsplit(s," "))[1]
  exp <- unlist(strsplit(s," "))[3]
  exp <- ifelse(is.na(exp),"all",exp)
  
  #Write out data
  cors_all[counter,1] <- sp
  cors_all[counter,2] <- exp
  cors_all[counter,3] <- mean_r2
  cors_all[counter,4] <- sd_r2
  counter=counter+1
}
cors_all
cors_all$species <- ifelse(cors_all$species=="Albacore","HMS",
                           ifelse(cors_all$species=="Anchovy", "CPS", "GFS"))
cors_all$species <- factor(cors_all$species, levels = c("HMS", "CPS", "GFS"))
g1 <- ggplot(data=cors_all, aes(x=species, y=forecast_mean_R2, fill=experiment))+
  geom_bar( position="dodge",stat = "identity")+
  geom_errorbar(aes(x=species, ymin=forecast_mean_R2-forecast_sd_R2,ymax=forecast_mean_R2+forecast_sd_R2),
                width=.4, position=position_dodge(.9))+
  theme_classic()+
  scale_fill_viridis(discrete = T, labels = c("Train 2040", "Full", "Temp only"),
                     name="Experiment", option="E")+
  scale_y_continuous(expand = c(0,0), limits=c(0,1))+
  labs(x="Species Archetype",y="mean correlation coefficient")

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig21.tiff', units="in", width = 6, height = 4, res=300)
plot(g1)
dev.off()


#------FIGURE 21.5: barplots of r2 for each ESM, Species, and experiment--------
#Get R2 for each model, species, ESM, and experiment (n=342)
cors_all <- as.data.frame(matrix(NA, nrow=9, ncol=4))
colnames(cors_all) <- c("species","experiment","historical_mean_R2", "historical_sd_R2")
counter=1
for(s in c("Albacore EMs", "Albacore EMs TempOnly", "Albacore EMs 2040train",
           "Anchovy EMs", "Anchovy EMs TempOnly","Anchovy EMs 2040train",
           "Groundfish EMs", "Groundfish EMs TempOnly", "Groundfish EMs 2040train")){
  print(s)
  #Only get forecast years 
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_dat_hist,gfdl_dat_hist,ipsl_dat_hist)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=nrow(had_dat_hist))
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","abundance")]
  #Get cor between observed and predicted for each model
  all_esm_agg <- all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                                                          "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(cor(., x=abundance, method="spearman")))
  #Get mean across EM:
  mean_r2 <- mean(as.matrix(all_esm_agg[,2:13]))
  sd_r2 <- sd(as.matrix(all_esm_agg[,2:13]))
  
  # #Pivot longer
  # all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(2:13),names_to="EM",values_to="R2")
  
  #Make species & experiment names
  sp <- unlist(strsplit(s," "))[1]
  exp <- unlist(strsplit(s," "))[3]
  exp <- ifelse(is.na(exp),"all",exp)
  
  #Write out data
  cors_all[counter,1] <- sp
  cors_all[counter,2] <- exp
  cors_all[counter,3] <- mean_r2
  cors_all[counter,4] <- sd_r2
  counter=counter+1
}
cors_all
cors_all$species <- ifelse(cors_all$species=="Albacore","HMS",
                           ifelse(cors_all$species=="Anchovy", "CPS", "GFS"))
cors_all$species <- factor(cors_all$species, levels = c("HMS", "CPS", "GFS"))
g1 <- ggplot(data=cors_all, aes(x=species, y=historical_mean_R2, fill=experiment))+
  geom_bar( position="dodge",stat = "identity")+
  geom_errorbar(aes(x=species, ymin=historical_mean_R2-historical_sd_R2,ymax=historical_mean_R2+historical_sd_R2),
                width=.4, position=position_dodge(.9))+
  theme_classic()+
  scale_fill_viridis(discrete = T, labels = c("Train 2040", "Full", "Temp only"),
                     name="Experiment", option="E")+
  scale_y_continuous(expand = c(0,0), limits=c(0,1))+
  labs(x="Species Archetype",y="mean correlation coefficient")

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig21.5.tiff', units="in", width = 6, height = 4, res=300)
plot(g1)
dev.off()


#------FIGURE 22: explore r2 and environmental anomalies-----
#Get R2 for each species and each year
#Load in predictions made in EM files
species_cors <- list()
counter_large=1
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
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- as.factor(rep(c("HAD","GFDL","IPSL"),each=58000))
  
  #create ensemble mean (remove Gam_S, Gam_EST, and Glm_s)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)})

  if (s == "Albacore EMs"){
    covar_names  <- c("temp","chla","mld") #use chla not zoo because models are trained & tested with this variable (ie. you are extrpoalting on chla not zoo)
  } 
  if (s == "Anchovy EMs"){
    covar_names  <- c("temp","chla")
  } 
  if (s == "Groundfish EMs"){
    covar_names  <- c("btemp","O2")
  } 
 
  all_esm_allyears$bin <- ifelse(all_esm_allyears$year<2011, "1985-2010",
                                 ifelse(all_esm_allyears$year>=2011 & all_esm_allyears$year <=2040, "2011-2040",
                                        ifelse(all_esm_allyears$year>=2041 & all_esm_allyears$year<=2070, "2041-2070","2071-2100")))
  counter=1
  export_data <- as.data.frame(matrix(NA,nrow=270, ncol=5))
  colnames(export_data) <- c("year","esm","mean_ExDat_under0","n_uni", "n_com")
  for (e in c("HAD","IPSL","GFDL")){
    print(e)
    train_df <- as.data.frame(all_esm_allyears[all_esm_allyears$esm==e & all_esm_allyears$year<=2010,])
    colnames(train_df)[1:2] <- c("x","y")
    for (y in 2011:2100){#c("2011-2040","2041-2070","2071-2100")){
      pred_df <- all_esm_allyears[all_esm_allyears$esm==e & all_esm_allyears$year==y,]
      colnames(pred_df)[1:2] <- c("x","y")
      aftt_crs <- sp::CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      
      test <-   compute_extrapolation(samples = train_df,
                                      covariate.names = covar_names,
                                      prediction.grid = pred_df, 
                                      coordinate.system=aftt_crs)
      export_data[counter,1] <- y
      export_data[counter,2] <- e
      export_data[counter,3] <- mean(abs(test$data$all$ExDet[test$data$all$ExDet < 0]))
      export_data[counter,4] <- length(test$data$all$ExDet[test$data$all$ExDet<0]) 
      export_data[counter,5] <- length(test$data$all$ExDet[test$data$all$ExDet>1]) 
      
      counter=counter+1
    }
  }
  
  # Calculate correlation between observed and predicted
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST","ens_mean"),funs(COR(., y=abundance)))
  #Bring together model R2 and ExDat means
  all_esm_agg_enviro <- left_join(all_esm_agg[all_esm_agg$year>=2011,],export_data,by=c("year","esm"))
  
  all_esm_agg_enviro_longer <- all_esm_agg_enviro %>% pivot_longer(c(3:15),names_to="EM", values_to="R2")
  
  species_cors[[counter_large]]<- all_esm_agg_enviro_longer
  counter_large=counter_large+1
}

##########


######
#PLOT THIS DATA
hms <- species_cors[[1]]
hms$n_uni_p <- hms$n_uni/500

cps <- species_cors[[2]]
cps$n_uni_p <- cps$n_uni/500

gfs <- species_cors[[3]]
gfs$n_uni_p <- gfs$n_uni/500

plot(hms$R2,hms$n_uni)
plot(hms$R2[hms$esm=="HAD"],hms$n_uni[hms$esm=="HAD"])
plot(hms$R2[hms$esm=="GFDL"],hms$n_uni[hms$esm=="GFDL"])
plot(hms$R2[hms$esm=="IPSL"],hms$n_uni[hms$esm=="IPSL"])
cor(hms$R2,hms$n_uni, use='na.or.complete') #-0.44 when <0 inc. z00; -0.19 when mean inc. zoo

plot(hms$n_uni[hms$esm=="HAD"]~hms$year[hms$esm=="HAD"])

hms <- species_cors[[2]]
plot(hms$R2,hms$n_uni)
plot(hms$R2[hms$esm=="HAD"],hms$n_com[hms$esm=="HAD"])
plot(hms$R2[hms$esm=="GFDL"],hms$n_com[hms$esm=="GFDL"])
plot(hms$R2[hms$esm=="IPSL"],hms$n_com[hms$esm=="IPSL"])
cor(hms$R2,hms$n_uni, use='na.or.complete') #-0.07 when <0 inc. z00; -0.23 when mean inc. zoo

hms <- species_cors[[3]]
plot(hms$R2,hms$n_uni)
plot(hms$R2[hms$esm=="HAD"],hms$n_uni[hms$esm=="HAD"])
plot(hms$R2[hms$esm=="GFDL"],hms$n_uni[hms$esm=="GFDL"])
plot(hms$R2[hms$esm=="IPSL"],hms$n_uni[hms$esm=="IPSL"])
cor(hms$R2,hms$n_uni, use='na.or.complete')



# plot(all_esm_agg_enviro$brt_E[all_esm_agg_enviro$esm=="HAD"],all_esm_agg_enviro$mean_ExDat[all_esm_agg_enviro$esm=="HAD"])
# plot(all_esm_agg_enviro$brt_E[all_esm_agg_enviro$esm=="IPSL"],all_esm_agg_enviro$mean_ExDat[all_esm_agg_enviro$esm=="IPSL"])
# plot(all_esm_agg_enviro$brt_E[all_esm_agg_enviro$esm=="GFDL"],all_esm_agg_enviro$mean_ExDat[all_esm_agg_enviro$esm=="GFDL"])
# 
# cor(all_esm_agg_enviro$brt_E[all_esm_agg_enviro$esm=="HAD"],all_esm_agg_enviro$mean_ExDat[all_esm_agg_enviro$esm=="HAD"], use="na.or.complete")
# cor(all_esm_agg_enviro$brt_E[all_esm_agg_enviro$esm=="IPSL"],all_esm_agg_enviro$mean_ExDat[all_esm_agg_enviro$esm=="IPSL"], use="na.or.complete")
# cor(all_esm_agg_enviro$brt_E[all_esm_agg_enviro$esm=="GFDL"],all_esm_agg_enviro$mean_ExDat[all_esm_agg_enviro$esm=="GFDL"], use="na.or.complete")



#------Figure 23: supp figure of population dyanmics CPS & GFS----

# GFS
years <- seq(1980,2100,1)
pop <-read.csv('~/PROJECTS/WRAP Location/Groundfish_2019_AssessmentRecruits.csv', header = TRUE)
pop <- pop[,c(1,6,10,11)]
pop$year <- pop$Year-1979
lower <- subset(pop,pop$Category=="bad"|pop$Category=="average")
lower.mean <- mean(lower$Age.0.logged)
lower.sd <- sd(lower$Age.0.logged)
higher <- subset(pop,pop$Category=="good"|pop$Category=="average")
higher.mean <- mean(higher$Age.0.logged)
higher.sd <- sd(higher$Age.0.logged)
set.seed(121) #this isn't what I used for models but it looks good to plotting purposes. Didn't keep track of seed, oops. 
Age.0p <- pop$Age.0
zero <- rlnorm(21, lower.mean, lower.sd) #1980-2000
half <-rlnorm(20, higher.mean, higher.sd) #2001-2020
first <- rlnorm(20, lower.mean, lower.sd) #2021-2040
second <-rlnorm(20, higher.mean, higher.sd) #2041-2060
third <-rlnorm(20, lower.mean, lower.sd) #2061-2080
fourth <-rlnorm(20, higher.mean, higher.sd) #2061-2100
Age.0p <- c(zero,half,first,second,third,fourth)
Age.0p <- as.data.frame(cbind(years,Age.0p))
Age.0p$year <- Age.0p$years-1979
# plot(Age.0p$years,Age.0p$Age.0p, type='l',xlab="years",ylab="Age-0 abundance")
# lines(Age.0p$years[1:39],Age.0p$Age.0p[1:39],col='blue')
Age.0p$Age.0p_norm <- BBmisc::normalize(Age.0p$Age.0p,method = "range",range=c(0.7,1)) #use this range to vary prevalence according to popn. dynamics below. KA: using this range creates prevalence values that are >1 down below which throws an error. I've ifelse'd this to simply use "1" when >1, but this doesn't allow for boom years to be accounted for properly???
# plot(Age.0p$years,Age.0p$Age.0p_norm, type='l',xlab="years",ylab="Age-0 normalized abundance")
# lines(Age.0p$years[1:39],Age.0p$Age.0p_norm[1:39],col='blue')




#CPS
pop <- read.csv('~/PROJECTS/WRAP Location/SUMMARYApreycleanversion.csv', header = TRUE)
pop <- pop[pop$Year>=1880,]
pop$year <- pop$Year - 1879
error <- sd(pop$AnchovySSB) #pretty large
pop <- pop[pop$Sim==1,]
pop$AnchovySSB_norm <- BBmisc::normalize(pop$AnchovySSB,method = "range",range=c(0.8,1)) #use this range to vary prevalence according to popn. dynamics
# plot(pop$Year,pop$AnchovySSB, type='l')
# abline(h=1172713)
# mean(pop$AnchovySSB)
# quantile(pop$AnchovySSB_norm)
# lowYears <- which(pop$AnchovySSB<1172713)


gfs <- ggplot(data=Age.0p[Age.0p$years>1984,], aes(x=years,y=Age.0p_norm))+
  geom_line()+
  ylab("Normalized Biomass")+
  xlab("")+
  ggtitle("GFS")

cps <- ggplot(data=pop[pop$year>5,], aes(x=year+1979,y=AnchovySSB_norm))+
  geom_line()+
  ylab("Normalized Biomass")+
  xlab("")+
  ggtitle("CPS")+
  geom_hline(yintercept =  0.8513197, linetype="dashed")

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig23.tiff', units="in",height=8,width=6, res=300)
grid.arrange(cps,gfs,nrow=2)
dev.off()




#-----Fig 24: whats up with the high frequency noise------
#Some something to do with teh 500 random samples more than extrapolation 

species_cors <- list()
counter_large=1
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
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- as.factor(rep(c("HAD","GFDL","IPSL"),each=58000))
  
  #create ensemble mean (remove Gam_S, Gam_EST, and Glm_s)
  all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",   
                                                         "brt_E" , "brt_ES", "brt_EST" ,  
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)})
    # Calculate correlation between observed and predicted
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST","ens_mean"),funs(COR(., y=abundance)))
  
  all_esm_agg_enviro_longer <- all_esm_agg_enviro %>% pivot_longer(c(3:15),names_to="EM", values_to="R2")
  
  species_cors[[counter_large]]<- all_esm_agg_enviro_longer
  counter_large=counter_large+1
}


pres <- aggregate(abundance~year, data = all_esm_allyears[all_esm_allyears$esm=="GFDL",], mean)
cor <- all_esm_agg[all_esm_agg$esm=="GFDL",]
plot(pres$abundance~cor$ens_mean)
cor(pres$abundance,cor$ens_mean)

plot(pres$abundance~pres$year, type='b')
plot(cor$gam_E~cor$year, col="blue", type='b')





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





  
  
  
  
  
  