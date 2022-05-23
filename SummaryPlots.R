
#Summary Plots
#WRAP Loc^3 call

#---library----

library(ggplot2)
library(reshape2)
# library(sf)
# library(oce)
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
library(cowplot)
library(ggpubr)
library(mgcv)
library(geosphere)
library(ecp)
library(tidyverse)
library(glue)
library(zoo)
library(scales)

#----FIGURE ONE----
#Conceptual diagram of OM and EMs
#Will be done in powerpoint.

#----FIGURE 2: environmental timeseries & maps----
#Show time-series of environmental covariates in each domain
# 3x4 plot
#col 1 & 3: spatial maps of historical variables (use gfdl)
#col 2 & 4: time-series of covariates for each ESM +/- SE

#Spatial maps
#For each variable in ensemble mean, average over years 1985-2010 and over months 4,5,6
spatial_mean_list <- list()
variables <- c("ild_0.5C","oxygen_bottom","sst","temp_bottom","zoo_200m","zoo_50m")
counter = 1
for (variable in variables){
  print(variable)
  nc <- nc_open(paste0('~/Dropbox/WRAP Location^3/2d_fields/monthly/',variable,'_monthly_roms_ensmean_1980_2100.nc'))
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
  theme_classic() +  labs(y="", x="", title = "Mixed Layer Depth (m)") +   theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = rev(cmocean("dense")(256))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))
map3 <- ggplot(data=spatial_mean_list[[3]], aes(x=Var1, y=Var2))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = "Sea Surface Temperature (°C)") +   theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("thermal")(256), name="Sea Surface Temperature") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 
map5 <- ggplot(data=spatial_mean_list[[3]], aes(x=Var1, y=Var2))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = bquote('Integrated Zoo 200 m '(mmol~N~m^-2)))+
  theme(legend.position="right",legend.title = element_blank(),
        plot.title = element_text(size = 10))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+#makes a box
  scale_fill_gradientn(colours = rev(cmocean("algae")(256)), name="Zooplankton Integrated to 200m") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0)) 

# title = expression(paste("Integrated Zooplankton \n(200 m; mmol N"~m^-2~")")))
# labs(x = bquote('x axis'~(Å^2)), y = "y axis") +

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
  theme_classic() +  labs(y="", x="", title = "Bottom Temperature (°C)") +   theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = cmocean("ice")(256), name="Bottom Temperature") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0), breaks = pretty(var_df$x, n = 3)) +  scale_y_continuous(expand = c(0, 0))

var_r <- rasterFromXYZ(spatial_mean_list[[2]])
template = raster('~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd')
var_r <- raster::resample(var_r, template, method="bilinear")  
var_r <- mask(var_r,rt)
var_df <- as.data.frame(rasterToPoints(var_r))
map2 <- ggplot(data=var_df, aes(x=x, y=y))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = expression(paste("Bottom Oxygen (mmol"~m^3~")")))  +  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = rev(cmocean("rain")(256)), name="Bottom Oxygen") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0), breaks = pretty(var_df$x, n = 3)) +  scale_y_continuous(expand = c(0, 0)) 


var_r <- rasterFromXYZ(spatial_mean_list[[6]])
template = raster('~/Dropbox/Eco-ROMS/ROMS & Bathym Data/Bathymetry ETOPO1/template.grd')
var_r <- raster::resample(var_r, template, method="bilinear")  
var_r <- mask(var_r,rt)
var_df <- as.data.frame(rasterToPoints(var_r))
map6 <- ggplot(data=var_df, aes(x=x, y=y))+
  geom_tile(aes(fill=value))+
  theme_classic() +  labs(y="", x="", title = bquote('Integrated Zoo 50 m '(mmol~N~m^-2)))+
  theme(legend.position="right",legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1),#makes a box
         plot.title = element_text(size = 10))+
  scale_fill_gradientn(colours = rev(cmocean("algae")(256)), name="Zooplankton Integrated to 50m") +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0), breaks = pretty(var_df$x, n = 3)) +  scale_y_continuous(expand = c(0, 0)) 

#Make Time-series plots: 6 variables, 3 ESMs
#Caution this takes a while --> load file I've already made
timeseries_mean_list <- readRDS('~/Dropbox/PROJECTS/WRAP Location/esm_timeseries_dataforplot.rds')
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
# saveRDS(timeseries_mean_list,'~/PROJECTS/WRAP Location/esm_timeseries_dataforplot.rds')


temp <- timeseries_mean_list[[1]] %>% group_by(esm) %>% summarise_at(vars("value","year"), funs(rollmean(.,k=11, fill=NA)))
ts1 <- ggplot(data = temp,aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Mixed Layer Depth (m)", x="") +   theme(legend.position="top")+
  scale_x_continuous(limits=c(1989,2096), expand = c(0, 0))
temp <- timeseries_mean_list[[2]] %>% group_by(esm) %>% summarise_at(vars("value","year"), funs(rollmean(.,k=11, fill=NA)))
ts2 <- ggplot(data = temp,aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y=expression(paste("Bottom Oxygen")), x="") +   theme(legend.position="top")+
  scale_x_continuous(limits=c(1989,2096), expand = c(0, 0))
temp <- timeseries_mean_list[[3]] %>% group_by(esm) %>% summarise_at(vars("value","year"), funs(rollmean(.,k=11, fill=NA)))
ts3 <- ggplot(data = temp,aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Sea Surface Temperature", x="") +   theme(legend.position="top")+
  scale_x_continuous(limits=c(1989,2096), expand = c(0, 0))
temp <- timeseries_mean_list[[4]] %>% group_by(esm) %>% summarise_at(vars("value","year"), funs(rollmean(.,k=11, fill=NA)))
ts4 <- ggplot(data = temp,aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  labs(y="Bottom Temperature (°C)", x="") +   theme(legend.position="top")+
  scale_x_continuous(limits=c(1989,2096), expand = c(0, 0))
temp <- timeseries_mean_list[[5]] %>% group_by(esm) %>% summarise_at(vars("value","year"), funs(rollmean(.,k=11, fill=NA)))
ts5 <- ggplot(data = temp,aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        # axis.title=element_text(size=9)
        )+
  labs(y=expression(paste("Integrated Zoo 200 m")), x="") +   theme(legend.position="top")+
  scale_x_continuous(limits=c(1989,2096), expand = c(0, 0))+
  scale_y_continuous(labels = label_number(accuracy = 1))
temp <- timeseries_mean_list[[6]] %>% group_by(esm) %>% summarise_at(vars("value","year"), funs(rollmean(.,k=11, fill=NA)))
ts6 <- ggplot(data = temp,aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        # axis.title=element_text(size=9)
        )+
  labs(y=expression(paste("Integrated Zoo 50 m")), x="") +   theme(legend.position="top")+
  scale_x_continuous(limits=c(1989,2096), expand = c(0, 0))+
  scale_y_continuous(labels = label_number(accuracy = 1))

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig2_smoothed.tiff',res=300,units="in",height=11,width=13)
h = 0.25
w = 0.25
ggdraw() +
  draw_plot(map1, 0, 0.75, h, w) +
  draw_plot(ts1, .25, 0.75, h, w) +
  draw_plot(map2, 0.5, 0.75, h, w) +
  draw_plot(ts2, 0.75, 0.75, h, w) +
  
  draw_plot(map3, 0, 0.5, h, w) +
  draw_plot(ts3, .25, 0.5, h, w) +
  draw_plot(map4, 0.5, 0.5, h, w) +
  draw_plot(ts4, 0.75, 0.5, h, w) +
  
  draw_plot(map5, 0, 0.25, h, w) +
  draw_plot(ts5, .25, 0.25, h, w) +
  draw_plot(map6, 0.5, 0.25, h, w) +
  draw_plot(ts6, 0.75, 0.25, h, w) 
dev.off()

#----Figure2.1: breakpoint analysis on each variable----

timeseries_mean_list <- readRDS('~/PROJECTS/WRAP Location/esm_timeseries_dataforplot.rds')
variables <- c("ild_0.5C","oxygen_bottom","sst","temp_bottom","zoo_200m","zoo_50m")

for (v in c(1:6)){
  print(glue("running variable {v}"))
  dat_temp <- timeseries_mean_list[[v]]
  dat_temp <- dat_temp[dat_temp$year>=2011,]
  dat_temp_long <- dat_temp %>% dplyr::select(year, value, esm) %>% pivot_wider(., names_from = "esm", values_from = "value")
  dat_mat <- as.matrix(dat_temp_long)
  fit_cpm1 = e.cp3o(dat_mat, minsize=5, K=5)
  print(fit_cpm1$estimates)
}

# dat_temp$year[c(42,84)]
# dat_temp$year[c(44,89)]
# dat_temp$year[c(43,86)]
# dat_temp$year[c(42,84)]
# dat_temp$year[c(43,87)]
# dat_temp$year[c(42,84)]

dat_temp$year[c(32,68)]
dat_temp$year[c(30,62)]
dat_temp$year[c(33,66)]
dat_temp$year[c(34,67)]
dat_temp$year[c(35,71)]
dat_temp$year[c(34,69)]

#So what changes between those periods
par(mfrow=c(2,3))
for (v in c(1:6)){
dat_temp <- timeseries_mean_list[[v]]
dat_temp <- dat_temp[dat_temp$year>=2011,]
dat_temp$breaks <- ifelse(dat_temp$year>=2011 & dat_temp$year<=2043,1,
                          ifelse(dat_temp$year>2043 & dat_temp$year<2077, 2, 3))
dat_temp_agg <- dat_temp %>% select(year, value, esm, breaks) %>% group_by(breaks, esm) %>% summarise(cv =mean(value))#/mean(value))
dat_temp_agg <- dat_temp_agg %>% group_by(breaks) %>% summarise(cv =sd(cv)/mean(cv))
barplot(dat_temp_agg$cv*100, main = variables[v], ylim=c(0,15))
}


#So what changes between those periods
par(mfrow=c(2,3))
for (v in c(1:6)){
  dat_temp <- timeseries_mean_list[[v]]
  dat_temp <- dat_temp[dat_temp$year>=2011,]
  dat_temp_agg <- dat_temp %>% select(year, value, esm) %>% group_by(year) %>% summarise(cv =sd(value)/mean(value))#/mean(value))
  plot(dat_temp_agg$year,dat_temp_agg$cv*100,'b')
}

#Calculate the rate of change then get mean across time period
par(mfrow=c(2,3))
for (v in c(1:6)){
  dat_temp <- timeseries_mean_list[[v]]
  dat_temp <- dat_temp %>% group_by(esm) %>%  arrange(esm, year) %>% mutate(roc = 100*(value - lag(value))/lag(value)) %>% ungroup()
  dat_temp$breaks <- ifelse(dat_temp$year<=2025,1,
                            ifelse(dat_temp$year>2025 & dat_temp$year<2068, 2, 3))
  dat_temp_agg <- dat_temp %>% select(year, roc, esm, breaks) %>% group_by(breaks) %>% summarise(sd=mean(roc, na.rm=T))
  # barplot(dat_temp_agg$sd, main = variables[v])
  plot(dat_temp$year, dat_temp$roc, col=as.factor(dat_temp$esm), 'b')
}

#----FIGURE 2.5: environmental timeseries across region----
#col 2 & 4: time-series of covariates for each ESM +/- SE
#Make Time-series plots: 6 variables, 3 ESMs
#Caution this takes a while --> load file I've already made
timeseries_mean_list <- readRDS('~/PROJECTS/WRAP Location/esm_timeseries_dataforplot_byregion.rds')
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
      grid <- expand.grid(lon,lat)
      month <- ncvar_get(nc, 'month')
      month <- month[61:1452] #select all years after 1985
      month_index <- which(month==4 | month==5 | month==6)
      var <- ncvar_get(nc,nc$var[[5]]$name)
      
      var_n <- var[,102:181,month_index]
      var_c <- var[,47:101,month_index]
      var_s <- var[,1:46,month_index]
      
      var_mean_n <- as.data.frame(apply(var_n, c(3),FUN = function(x){mean(x,na.rm=T)}))
      colnames(var_mean_n) <- "value"
      var_mean_c <- as.data.frame(apply(var_c, c(3),FUN = function(x){mean(x,na.rm=T)}))
      colnames(var_mean_c) <- "value"
      var_mean_s <- as.data.frame(apply(var_s, c(3),FUN = function(x){mean(x,na.rm=T)}))
      colnames(var_mean_s) <- "value"
      
      var_mean_n$region <- "north"
      var_mean_c$region <- "central"
      var_mean_s$region <- "south"
      var_mean <- rbind(var_mean_n,var_mean_c,var_mean_s)

      # var_mean$sd <- apply(var, c(3),FUN = function(x){sd(x,na.rm=T)})
      var_mean$month <- rep(rep(c(4,5,6),116),3)
      var_mean$year <- rep(rep(1985:2100,each=3),3)
      var_mean_agg <- aggregate(value ~ year + region, data=var_mean, FUN = function(x){mean(x,na.rm=T)})
      # var_sd_agg <- aggregate(sd ~ year, data=var_mean, FUN = function(x){sd(x,na.rm=T)})
      # var_agg <- left_join(var_mean_agg,var_sd_agg, by="year")
      
      esm_temp_list[[esm_count]] <- var_mean_agg
      esm_count <- esm_count +1
    } else {
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
      
      #Trim other variables to smaller domain
      stack <- stack(list.files(paste0('~/Dropbox/WRAP Location^3/Rasters_2d_monthly/',esm,'/',variable,'/'),full.names = TRUE, pattern = '.grd'))
      stack2 <- raster::mask(stack,rt)
      
      stack2_n <- crop(stack2, extent(-134,-115.5,40,48))
      stack2_c <- crop(stack2, extent(-134,-115.5,34.5,40)) 
      stack2_s <- crop(stack2, extent(-134,-115.5,30,34.5))
      
      stats_n <- as.data.frame(cellStats(stack2_n,'mean'))
      stats_c <- as.data.frame(cellStats(stack2_c,'mean'))
      stats_s <- as.data.frame(cellStats(stack2_s,'mean'))
      
      colnames(stats_n) <- "stats"
      colnames(stats_c) <- "stats"
      colnames(stats_s) <- "stats"

      stats_n$region <- "north"
      stats_c$region <- "central"
      stats_s$region <- "south"
            
      stats <- rbind(stats_n,stats_c,stats_s)
      
      stats$year <- rep(rep(1980:2100, each=12),3)
      stats$month <- rep(rep(1:12, each=121),3)
      stats <- stats[stats$year>=1985,]
      stats_agg <- aggregate(stats ~ year + region, data=stats, FUN = function(x){mean(x,na.rm=T)})
      # stats_agg_sd <- aggregate(stats ~ year, data=stats, FUN = function(x){sd(x,na.rm=T)})
      # var_agg <- left_join(stats_agg,stats_agg_sd, by="year")
      # colnames(var_agg) <- c("year","value","sd")
      esm_temp_list[[esm_count]] <- stats_agg
      esm_count <- esm_count +1
    }
  }
  df1 <- esm_temp_list[[1]] 
  df2 <- esm_temp_list[[2]]
  df3 <- esm_temp_list[[3]]
  df <- rbind(df1,df2,df3)
  df$esm <- rep(c("HAD","GFDL","IPSL"), each=348)
  # df$ymin <- df$value - df$sd/sqrt(3)
  # df$ymax <- df$value + df$sd/sqrt(3)
  timeseries_mean_list[[counter]] <- df
  counter = counter +1
}
summary(timeseries_mean_list)
saveRDS(timeseries_mean_list,'~/PROJECTS/WRAP Location/esm_timeseries_dataforplot_byregion.rds')

timeseries_mean_list[[1]]$variable <- "mld"
timeseries_mean_list[[2]]$variable <- "bottom_oxygen"
colnames(timeseries_mean_list[[2]])[3] <- "value"
timeseries_mean_list[[3]]$variable <- "sst"
timeseries_mean_list[[4]]$variable <- "bottom_temp"
colnames(timeseries_mean_list[[4]])[3] <- "value"
timeseries_mean_list[[5]]$variable <- "zoo_200"
timeseries_mean_list[[6]]$variable <- "zoo_50"
colnames(timeseries_mean_list[[6]])[3] <- "value"

dat <- do.call(rbind,timeseries_mean_list)

ggplot(dat,aes(x=year,y=value,fill=esm))+
  geom_line(aes(col=esm),size=1)+
  facet_grid(region~ variable, scales="free")

# get sd of each variable in each region, then plot
test <- dat %>% group_by(variable, region, esm) %>% summarise(mean(value))
test2 <- test %>% group_by(variable, region) %>% summarise(sd(`mean(value)`))

test2$region <- as.factor(test2$region) 
levels(test2$region)
test2$region <- factor(test2$region, levels = c("north","central","south"))

test2$variable <- as.factor(test2$variable) 
# library(plyr)
test2$variable <- plyr::revalue(test2$variable, c("bottom_oxygen" = "Bottom Oxygen",
                          "bottom_temp"  = "Bottom Temperature",
                          "mld"  = "MLD",
                          "sst"   = "SST",
                          "zoo_200"   = "Zooplankton 200m",
                          "zoo_50"  = "Zooplankton 50m"))

# g <- ggplot(test2, aes(x=variable,y=`sd(\`mean(value)\`)`, fill=region))+
#   geom_bar(position="dodge", stat="identity")+
#   scale_fill_manual(values = c("#d1495b","#00798c","#edae49"))+
#   labs(x="Environmental Covariate", y="Variability among ESMs")+
#   theme_bw()+theme(legend.position = "none", legend.title = element_blank())

g <- ggplot(test2, aes(x=region,y=`sd(\`mean(value)\`)`, fill=variable))+
  geom_bar(position="dodge", stat="identity")+
  # scale_fill_manual(values = c("#d1495b","#00798c","#edae49"))+
  labs(x="Region", y="Variability among ESMs")+
  theme_bw()+theme(legend.position = "bottom", legend.title = element_blank())


a <- ggplot()+
  annotation_map(map_data("world"), colour = "black", fill="light grey")+
  coord_quickmap(xlim=c(-126,-115.5),ylim=c(30,48)) +  #Sets aspect ratio
  geom_vline(xintercept=-115.5, col="grey50", size=0.01)+ #add to force xaxis
  geom_hline(yintercept = 40, linetype="dashed")+
  geom_hline(yintercept = 34.5, linetype="dashed")+
  xlab("") + ylab("")+
  theme_classic()+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.8))+  #makes a box
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  # ggtitle("Regions")+
  annotate(geom="text", x=-122, y=45, label="north", hjust=0, fontface=2, cex=3, col="#d1495b")+
  annotate(geom="text", x=-122, y=38, label="central", hjust=0, fontface=2, cex=3, col="#00798c")+
  annotate(geom="text", x=-122, y=32, label="south", hjust=0, fontface=2, cex=3, col="#edae49")+
  theme(panel.background = element_rect(fill = "transparent"), # bg of the panel
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.text.x=element_blank())

g1 <- ggdraw() +
  draw_plot(g) +
  draw_plot(a,x=-0.36, y=0.48,height = 0.48)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig2.5.tiff', units="in",res=300,width=8,height=6)
plot_grid(g1)
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
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
  # all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  # all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  # all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  all_esm_agg <- all_esm_agg %>% group_by(esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                        "glm_E" ,"glm_ESt", "glm_ESr",
                                                                        "brt_E" , "brt_ES", "brt_EST" ,
                                                                        "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean","year"),funs(rollmean(.,k=11,fill=NA)))
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

dat <- species_biomass[[1]]
# dat <- dat %>% group_by(esm,EM) %>% summarise_at(vars("min","max","abundance","year"),funs(rollmean(.,k=11,fill=NA)))
hms_biomass <- ggplot(data=dat,aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated biomass", "Ensemble mean"))+
  geom_ribbon(alpha=0.2)+
  theme_classic()+
  labs(x="",y="HMS Biomass (kg)")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)+
  scale_x_continuous(limits = c(1990,2095), expand=c(0,0))

dat <- species_biomass[[2]]
# dat <- dat %>% group_by(esm,EM) %>% summarise_at(vars("min","max","abundance","year"),funs(rollmean(.,k=11,fill=NA)))
cps_biomass <- ggplot(data=dat,aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated biomass", "Ensemble mean"))+
  geom_ribbon(alpha=0.2)+
  theme_classic()+
  labs(x="",y="CPS Biomass (kg)")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)+
  scale_x_continuous(limits = c(1990,2095), expand=c(0,0))

dat <- species_biomass[[3]]
# dat <- dat %>% group_by(esm,EM) %>% summarise_at(vars("min","max","abundance","year"),funs(rollmean(.,k=11,fill=NA)))
gfs_biomass <- ggplot(data=dat,aes(x=year,y=abundance, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#2980b9","#c0392b"),labels = c("Simulated biomass", "Ensemble mean"))+
  geom_ribbon(alpha=0.2)+
  theme_classic()+
  labs(x="",y="GFS Biomass (kg)")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)+
  scale_x_continuous(limits = c(1990,2095), expand=c(0,0))

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig3_smoothed.tiff',res=300,units="in",height=6,width=9)
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

#----Figure 4.1: Cor timeseries of each model------
#WORK IN PROGRESS###
#Load in predictions made in EM files
species_rmse_diff <- list()
counter=1
for(s in c("Albacore EMs","Albacore EMs 2040train","Albacore EMs TempOnly",
           "Anchovy EMs","Anchovy EMs 2040train","Anchovy EMs TempOnly",
           "Groundfish EMs","Groundfish EMs 2040train","Groundfish EMs TempOnly")){
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
  # all_esm_allyears$time_bin <- ifelse(all_esm_allyears$year<=2010,"fit",
  #                                     ifelse(all_esm_allyears$year>=2056,"2056-2100","2011-2055"))
  # all_esm_allyears$time_bin <- ifelse(all_esm_allyears$year<=2010,"fit","future")
  
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(COR(., y=abundance)))
  all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
  
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm"), variable.name = "EM",value.name = "RMSE")
  # all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  
  #Rank order
  all_esm_agg_longer <- all_esm_agg_longer %>%  group_by(year,esm) %>% mutate(rank = rank(RMSE, ties.method = "first"))
  
  species_rmse_diff[[counter]]<- all_esm_agg_longer
  counter=counter+1
}



test <- rbind(species_rmse_diff[[1]],species_rmse_diff[[2]],species_rmse_diff[[3]],
               species_rmse_diff[[4]],species_rmse_diff[[5]],species_rmse_diff[[6]],
               species_rmse_diff[[7]],species_rmse_diff[[8]],species_rmse_diff[[9]])

full <- rbind(species_rmse_diff[[1]],
              species_rmse_diff[[4]],
              species_rmse_diff[[7]])
t2040 <- rbind(species_rmse_diff[[2]],
              species_rmse_diff[[5]],
              species_rmse_diff[[8]])
temp <- rbind(species_rmse_diff[[3]],
              species_rmse_diff[[6]],
              species_rmse_diff[[9]])
hms <- rbind(species_rmse_diff[[1]],
              species_rmse_diff[[2]],
              species_rmse_diff[[3]])
cps <- rbind(species_rmse_diff[[4]],
               species_rmse_diff[[5]],
               species_rmse_diff[[6]])
gfs <- rbind(species_rmse_diff[[7]],
              species_rmse_diff[[8]],
              species_rmse_diff[[9]])


plot(test$EM[test$year<2011],test$rank[test$year<2011])
plot(test$EM[test$year>=2011],test$rank[test$year>=2011])
plot(test$EM[test$year>=2011],test$RMSE[test$year>=2011])

par(mfrow=c(3,1))
plot(full$EM[full$year<2011],full$rank[full$year<2011])
plot(full$EM[full$year>=2011],full$rank[full$year>=2011])
plot(full$EM[full$year>=2011],full$RMSE[full$year>=2011])

plot(t2040$EM[t2040$year<2011],t2040$rank[t2040$year<2011])
plot(t2040$EM[t2040$year>=2011],t2040$rank[t2040$year>=2011])
plot(t2040$EM[t2040$year>=2011],t2040$RMSE[t2040$year>=2011])

plot(temp$EM[temp$year<2011],temp$rank[temp$year<2011])
plot(temp$EM[temp$year>=2011],temp$rank[temp$year>=2011])
plot(temp$EM[temp$year>=2011],temp$RMSE[temp$year>=2011])

par(mfrow=c(3,1))
plot(hms$EM[hms$year>=2011],hms$rank[hms$year>=2011])
plot(hms$EM[hms$year>=2011],hms$RMSE[hms$year>=2011])

plot(cps$EM[cps$year>=2011],cps$rank[cps$year>=2011])
plot(cps$EM[cps$year>=2011],cps$RMSE[cps$year>=2011])

plot(gfs$EM[gfs$year>=2011],gfs$rank[gfs$year>=2011])
plot(gfs$EM[gfs$year>=2011],gfs$RMSE[gfs$year>=2011])


#Group by E, ES, EST and make summary plots
#Group by mlp, gam, brt, glm and make summary plots

mean(species_rmse_diff[[1]]$RMSE)
mean(species_rmse_diff[[2]]$RMSE)
mean(species_rmse_diff[[3]]$RMSE)
mean(species_rmse_diff[[4]]$RMSE)
mean(species_rmse_diff[[5]]$RMSE)
mean(species_rmse_diff[[6]]$RMSE)
mean(species_rmse_diff[[7]]$RMSE)
mean(species_rmse_diff[[8]]$RMSE)
mean(species_rmse_diff[[9]]$RMSE)


ggplot(data=species_rmse_diff[[1]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)
ggplot(data=species_rmse_diff[[2]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)
ggplot(data=species_rmse_diff[[3]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)


ggplot(data=species_rmse_diff[[4]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)
ggplot(data=species_rmse_diff[[5]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)
ggplot(data=species_rmse_diff[[6]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)

ggplot(data=species_rmse_diff[[7]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)
ggplot(data=species_rmse_diff[[8]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)
ggplot(data=species_rmse_diff[[9]],aes(x=time_bin,y=RMSE,group=EM))+
  geom_line(aes(col=EM))+
  # scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  # geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: Average RMSE")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  facet_wrap(~esm)





#----FIGURE 4.4: Correlation timeseries-----
#Compare point x point correlation for each year and EM

#Load in predictions made in EM files
species_rmse_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","pres","suitability","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm","ens_mean")]
  
  # all_esm_allyears <- all_esm_allyears[all_esm_allyears$pres==0,]
  #Calculate RMSE
  COR = function(x,y){cor(x,y,method="spearman")}
  RMSE = function(x, y){(sqrt(mean((x - y)^2)))}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(RMSE(., y=abundance)))
  all_esm_agg <- all_esm_agg %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                     "glm_E" ,"glm_ESt", "glm_ESr",
                                                                     "brt_E" , "brt_ES", "brt_EST" ,
                                                                     "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean","year"),funs(rollmean(.,k=11,fill=NA)))
  all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(2:13), names_to =  "EM",values_to = "RMSE")
  
  # all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(3:14), names_to =  "EM",values_to = "RMSE")
  
  species_rmse_diff[[counter]]<- all_esm_agg_longer

  counter=counter+1
}


extrap <- readRDS('~/Dropbox/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet_PlusNearby.rds')
#hms
hms1 <- extrap[[1]]; hms2 <- extrap[[2]]; hms3 <- extrap[[3]]
hms1$esm = "HAD"; hms2$esm = "GFDL"; hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms <- hms %>% group_by(esm) %>% summarise_at(vars("extrap_percent","nearby_mean","excluded","year"), funs(rollmean(.,k=11,fill=NA))) 
#cps
cps1 <- extrap[[4]];cps2 <- extrap[[5]];cps3 <- extrap[[6]]
cps1$esm = "HAD";cps2$esm = "GFDL";cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps <- cps %>% group_by(esm) %>% summarise_at(vars("extrap_percent","nearby_mean","excluded","year"), funs(rollmean(.,k=11,fill=NA)))
#gfs
gfs1 <- extrap[[7]];gfs2 <- extrap[[8]];gfs3 <- extrap[[9]]
gfs1$esm = "HAD";gfs2$esm = "GFDL";gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs <- gfs %>% group_by(esm) %>% summarise_at(vars("extrap_percent","nearby_mean","excluded","year"), funs(rollmean(.,k=11,fill=NA)))

#Albacore
dat <- species_rmse_diff[[1]]
dat2 <- left_join(dat,hms,by=c("year","esm"))
dat2 <- dat2 %>%  separate(EM, c("type","params"))
# dat2 <- dat2[dat2$params=="EST",]
# ggplot(dat2, aes(x=excluded, y=RMSE, group=EM))+
#   geom_point(aes(col=EM))+
#   geom_smooth(aes(col=EM))+
#   facet_wrap(~esm, scale="free")
# ggplot(dat2, aes(x=extrap_percent, y=RMSE))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   facet_wrap(~esm+EM, scales = "free")
ggplot(dat2, aes(x=extrap_percent, y=RMSE, group=type))+
  geom_point()+
  geom_smooth(aes(col=type))+
  facet_wrap(~esm)

#Anchovy
dat <- species_rmse_diff[[2]]
dat2 <- left_join(dat,cps,by=c("year","esm"))
dat2 <- dat2 %>%  separate(EM, c("type","params"))
# ggplot(dat2, aes(x=excluded, y=RMSE, group=EM))+
#   geom_point(aes(col=EM))+
#   geom_smooth(aes(col=EM),method = "lm")+
#   facet_wrap(~esm)
# ggplot(dat2, aes(x=extrap_percent, y=RMSE))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   facet_wrap(~esm+EM, scales = "free")
ggplot(dat2, aes(x=extrap_percent, y=RMSE, group=type))+
  geom_point()+
  geom_smooth(aes(col=type))+
  facet_wrap(~esm)

dat <- species_rmse_diff[[3]]
dat2 <- left_join(dat,gfs,by=c("year","esm"))
dat2 <- dat2 %>%  separate(EM, c("type","params"))
# ggplot(dat2, aes(x=excluded, y=RMSE, group=EM))+
#   geom_point(aes(col=EM))+
#   geom_smooth(aes(col=EM),method = "lm")+
#   facet_wrap(~esm)
# ggplot(dat2, aes(x=extrap_percent, y=RMSE))+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   facet_wrap(~esm+EM, scales = "free")
ggplot(dat2, aes(x=extrap_percent, y=RMSE, group=type))+
  geom_point()+
  geom_smooth(aes(col=type))+
  facet_wrap(~esm)

timeseries_mean_list <- readRDS('~/Dropbox/PROJECTS/WRAP Location/esm_timeseries_dataforplot.rds')
zoo50 <- timeseries_mean_list[[6]]
zoo200 <- timeseries_mean_list[[3]]
zoo200 <- zoo200 %>% group_by(esm) %>%   summarise_at(vars("value","year"),funs(rollmean(.,k=11, fill=NA)))

new <- left_join(dat2,ts1,by=c("year","esm"))
head(new)
cor(new$RMSE[new$esm=="GFDL"],new$value[new$esm=="GFDL"], use='na.or.complete')
cor(new$RMSE[new$esm=="IPSL"],new$value[new$esm=="IPSL"], use='na.or.complete')
cor(new$RMSE[new$esm=="HAD"],new$value[new$esm=="HAD"], use='na.or.complete')
ggplot(new,aes(x=value,y=RMSE))+
  geom_point(aes(col=EM))+
  facet_wrap(~esm)


dat <- readRDS('~/Dropbox/PROJECTS/WRAP Location/Albacore EMs/ipsl/dat_fcast_results_full.rds')
# test <- dat %>%  group_by(year) %>% summarise_at(vars("zoo_200"),funs(cor(.,y=chla)))
test <- test %>%  summarise_at(vars("zoo_200","year"),funs(rollmean(.,k=11)))

plot(dat$year,dat$zoo_200, type='b', main="ipsl")
plot(dat2$year,dat2$RMSE)

dat <- readRDS('~/Dropbox/PROJECTS/WRAP Location/Anchovy EMs/ipsl/dat_fcast_results_full.rds')
test <- dat %>%  group_by(year) %>% summarise_at(vars("zoo_50"),funs(cor(.,y=chla)))
test <- test %>%  summarise_at(vars("zoo_50","year"),funs(rollmean(.,k=11)))
plot(test$year,test$zoo_50, type='b', main="ipsl")

# # Zoo~Chl correlation through time
timeseries_mean_list <- list()
variables <- c("chl_surface","zoo_200m")
counter = 1
for (variable in variables){
  print(variable)
  esm_count <- 1
  esm_temp_list <- list()
  for (esm in c("had","gfdl","ipsl")){
    print(esm)
    if(variable=="chl_surface" | variable=="zoo_200m"){
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

cor(timeseries_mean_list[[1]]$value, timeseries_mean_list[[2]]$value)
cor(timeseries_mean_list[[1]]$value, timeseries_mean_list[[3]]$value)

ts1 <- timeseries_mean_list[[1]] %>%  group_by(esm) %>% summarise_at(vars("value","year"),funs(rollmean(.,k=11)))
plot(ts1$year,ts1$value, col=as.factor(ts1$esm), type='b')
plot(ts1$year,ts1$value, col=as.factor(ts1$esm), type='b')


plot(timeseries_mean_list[[1]]$value[timeseries_mean_list[[1]]$esm=="HAD"]~ timeseries_mean_list[[1]]$year[timeseries_mean_list[[1]]$esm=="HAD"], type='b')
plot(timeseries_mean_list[[2]]$value[timeseries_mean_list[[2]]$esm=="HAD"]~ timeseries_mean_list[[2]]$year[timeseries_mean_list[[2]]$esm=="HAD"], type='b', col='red')

df <- leftjoin()

d1 <- timeseries_mean_list[[1]][timeseries_mean_list[[1]]$esm=="HAD",]
d2 <- timeseries_mean_list[[2]][timeseries_mean_list[[2]]$esm=="HAD",]
dj <- left_join(d1[,1:2],d2[,1:2], by='year')
ggplot(df, aes(value.x, value.y, col=year))+
         geom_line()
       

#get hist of fited data on top of response curve?       
virtualspecies::plotResponse(spA_suitability, parameters="zoo") #plot response curves
dnorm(16,mean = 16,sd = 6)
spA_parameters <- formatFunctions(sst = c(fun="dnorm",mean=16,sd=6),
                                  zoo = c(fun="logisticFun",alpha=-5,beta=20),
                                  z = c(fun="logisticFun",alpha=-500,beta=-2000))

#plot histogram of fitted data
had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Anchovy EMs/had/dat_hist_results_full.rds"))
had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Anchovy EMs/had/dat_fcast_results_full.rds"))
hist(had_dat_hist$zoo_50[had_dat_hist$pres==1])

s="Anchovy EMs"
had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
#Rbind historical and forecast years
had_allyears <- rbind(had_dat_hist,had_dat_fcast)
gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
all_esm_allyears$esm <- as.factor(rep(c("HAD","GFDL","IPSL"),each=58000))
all_esm_allyears$type <- ifelse(all_esm_allyears$year>=2011,"fcast","hcast")

ggplot(all_esm_allyears, aes(zoo_50, group=type, fill=type))+
  geom_density(alpha=0.5) +
  facet_wrap(~esm)

load('~/Dropbox/PROJECTS/WRAP Location/Anchovy EMs/had/saved_models_full.RData')  
plot(gam_E_N, pages=1, scale=0)
plot(gam_E_P, pages=1, scale=0)

#then check whether spA_suitability chanegs with ESM


#----FIGURE 4.5: Correlation timeseries-----
#Compare point x point correlation for each year and EM

#Load in predictions made in EM files
species_rmse_diff <- list()
species_thresh <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
  all_esm_agg <- all_esm_agg %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                     "glm_E" ,"glm_ESt", "glm_ESr",
                                                                     "brt_E" , "brt_ES", "brt_EST" ,
                                                                     "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean","year"),funs(rollmean(.,k=11,fill=NA)))
  # all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  # all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  # all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
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
  
  #get threshold 
  thresh <- melt(all_esm_agg,id=c("year", "esm"), variable.name = "EM",value.name = "Cor")
  th <- thresh %>% filter(year<=2010) %>% filter(EM!="ens_mean") %>% filter(EM!="min") %>% filter(EM!="max") %>% 
    group_by(esm) %>% summarise_at("Cor", funs(mean, sd))
  th$two_sd <- th$mean - th$sd*2
  th$three_sd <- th$mean - th$sd*3
  species_thresh[[counter]] <- th
  
  counter=counter+1
}

dat <- species_rmse_diff[[1]]
hms_rmse_diff <- ggplot(data=dat,aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  theme(legend.position = "none")+
  facet_wrap(~esm)+
  scale_x_continuous(limits = c(1990,2095), expand=c(0,0))

dat <- species_rmse_diff[[2]]
cps_rmse_diff <- ggplot(data=dat,aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  theme(legend.position = "none")+
  facet_wrap(~esm)+
  scale_x_continuous(limits = c(1990,2095), expand=c(0,0))

dat <- species_rmse_diff[[3]]
gfs_rmse_diff <- ggplot(data=dat,aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="GFS: Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  theme(legend.position = "none")+
  facet_wrap(~esm)+
  scale_x_continuous(limits = c(1990,2095), expand=c(0,0))

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig4.5.tiff',res=300,units="in",height=8,width=12)
grid.arrange(hms_rmse_diff, cps_rmse_diff, gfs_rmse_diff, nrow=3)
dev.off()



#----FIGURE 4.5.0: Correlation timeseries WITH EXTRAPOLATION-----
#Compare point x point correlation for each year and EM

#Load in predictions made in EM files
species_rmse_diff <- list()
species_thresh <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- as.factor(rep(c("HAD","GFDL","IPSL"),each=58000))
  #create ensemble mean (remove Gam_S & Gam_EST)
  # all_esm_allyears$ens_mean <- apply(all_esm_allyears[,c("gam_E","gam_ES", "gam_ECor",
  #                                                        "glm_E" ,"glm_ESt", "glm_ESr",   
  #                                                        "brt_E" , "brt_ES", "brt_EST" ,  
  #                                                        "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) #removing gam_S from ensemble
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  
  #Calculate RMSE
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST"),funs(COR(., y=abundance)))
  all_esm_agg <- all_esm_agg %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                     "glm_E" ,"glm_ESt", "glm_ESr",
                                                                     "brt_E" , "brt_ES", "brt_EST" ,
                                                                     "mlp_E" , "mlp_ES" ,"mlp_EST","year"),funs(rollmean(.,k=11,fill=NA)))
  all_esm_agg$ens_mean <- apply(all_esm_agg[,c("gam_E","gam_ES", "gam_ECor",
                                                         "glm_E" ,"glm_ESt", "glm_ESr",
                                                         "brt_E" , "brt_ES", "brt_EST" ,
                                                         "mlp_E" , "mlp_ES" ,"mlp_EST" )],MARGIN = 1, FUN=function(x) {mean(x,na.rm=T)}) 

  
  # all_esm_agg$mlp_E <- ifelse(all_esm_agg$mlp_E==0,NA,all_esm_agg$mlp_E)
  # all_esm_agg$mlp_ES <- ifelse(all_esm_agg$mlp_ES==0,NA,all_esm_agg$mlp_ES)
  # all_esm_agg$mlp_EST <- ifelse(all_esm_agg$mlp_EST==0,NA,all_esm_agg$mlp_EST)
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
  
  #get threshold 
  thresh <- melt(all_esm_agg,id=c("year", "esm"), variable.name = "EM",value.name = "Cor")
  th <- thresh %>% filter(year<=2010) %>% filter(EM!="ens_mean") %>% filter(EM!="min") %>% filter(EM!="max") %>% 
    group_by(esm) %>% summarise_at("Cor", funs(mean, sd))
  th$two_sd <- th$mean - th$sd*2
  th$three_sd <- th$mean - th$sd*3
  species_thresh[[counter]] <- th
  
  counter=counter+1
}

#What are the quantiles of extrapolation
extrap <- readRDS('~/Dropbox/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet_PlusNearby.rds')
#hms
hms1 <- extrap[[1]]; hms2 <- extrap[[2]]; hms3 <- extrap[[3]]
hms1$esm = "HAD"; hms2$esm = "GFDL"; hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms <- hms %>% group_by(esm) %>% summarise_at(vars("extrap_percent","nearby_mean","excluded","year"), funs(rollmean(.,k=11,fill=NA))) 
#cps
cps1 <- extrap[[4]];cps2 <- extrap[[5]];cps3 <- extrap[[6]]
cps1$esm = "HAD";cps2$esm = "GFDL";cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps <- cps %>% group_by(esm) %>% summarise_at(vars("extrap_percent","nearby_mean","excluded","year"), funs(rollmean(.,k=11,fill=NA)))
#gfs
gfs1 <- extrap[[7]];gfs2 <- extrap[[8]];gfs3 <- extrap[[9]]
gfs1$esm = "HAD";gfs2$esm = "GFDL";gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs <- gfs %>% group_by(esm) %>% summarise_at(vars("extrap_percent","nearby_mean","excluded","year"), funs(rollmean(.,k=11,fill=NA)))


myPlot <- function(data,title){
  data <- data %>% mutate_at(vars("extrap_percent"),funs("ex_scaled" = scales::rescale(.,to=c(min(min),max(max)),from=range(.)))) 
  p <- ggplot(data=data,aes(x=year, ymin=min,ymax=max))+
    geom_line(aes(y=RMSE,col=EM),size=1)+
    geom_line(aes(y=ex_scaled,col="red"),size=1)+
    geom_ribbon(alpha=0.5)+
    theme_classic()+
    labs(x="",y="")+
    ggtitle(title)+
    theme(legend.title = element_blank())+
    theme(legend.position = "none",
          # panel.border = element_rect(color = "black",
          #                                   fill = NA,
          #                                   size = 1.5)
          )+
    annotate(geom = 'segment', y = Inf, yend = Inf, color = 'black', x = -Inf, xend = Inf, size = 1) +
    scale_x_continuous(limits = c(2016,2095), expand=c(0,0))+
    scale_y_continuous( sec.axis = sec_axis(name="",~scales::rescale(., to = range(data$extrap_percent), from = range(c(min(data$min),max(data$max))))))
  p
}

#HMS
dat <- species_rmse_diff[[1]]
dat2 <- left_join(dat,hms,by=c("esm","year"))
dat2 <- na.omit(dat2)
dat2_gfdl <- dat2[dat2$esm=="GFDL",]
dat2_ipsl <- dat2[dat2$esm=="IPSL",]
dat2_had <- dat2[dat2$esm=="HAD",]
pl_h_g <- myPlot(dat2_gfdl,title=glue("HMS: GFDL"))
pl_h_i <- myPlot(dat2_ipsl,title=glue("HMS: IPSL"))
pl_h_h <- myPlot(dat2_had,title=glue("HMS: HAD"))

dat <- species_rmse_diff[[2]]
dat2 <- left_join(dat,cps,by=c("esm","year"))
dat2 <- na.omit(dat2)
dat2_gfdl <- dat2[dat2$esm=="GFDL",]
dat2_ipsl <- dat2[dat2$esm=="IPSL",]
dat2_had <- dat2[dat2$esm=="HAD",]
pl_c_g <- myPlot(dat2_gfdl,title=glue("CPS: GFDL"))
pl_c_i <- myPlot(dat2_ipsl,title=glue("CPS: IPSL"))
pl_c_h <- myPlot(dat2_had,title=glue("CPS: HAD"))

dat <- species_rmse_diff[[3]]
dat2 <- left_join(dat,gfs,by=c("esm","year"))
dat2 <- na.omit(dat2)
dat2_gfdl <- dat2[dat2$esm=="GFDL",]
dat2_ipsl <- dat2[dat2$esm=="IPSL",]
dat2_had <- dat2[dat2$esm=="HAD",]
pl_g_g <- myPlot(dat2_gfdl,title=glue("GFS: GFDL"))
pl_g_i <- myPlot(dat2_ipsl,title=glue("GFS: IPSL"))
pl_g_h <- myPlot(dat2_had,title=glue("GFS: HAD"))

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig4.5.0.tiff',res=300,units="in",height=7,width=11)
grid.arrange(pl_h_g,pl_h_h, pl_h_i,
             pl_c_g,pl_c_h, pl_c_i,
             pl_g_g,pl_g_h, pl_g_i, nrow=3)
dev.off()


# #At which year does quantile occur
# hms$species <- "hms"
# cps$species <- "cps"
# gfs$species <- "gfs"
# sp <- rbind(hms, cps, gfs)
# hms_q <- sp %>% group_by(esm,species) %>% summarise_at(vars("extrap_percent"), funs("q_percent" = quantile(., probs=c(0.25,0.5,0.75),na.rm=T)))
# hms_q$q <- rep(c(25,50,75),9)
# hms_full <- left_join(sp,hms_q,by=c("esm","species"))
# y_qs <- hms_full %>% group_by(esm,q,species) %>% summarise_at(vars("q_percent"),funs(year[which.min(abs(extrap_percent-q_percent))])) 
# 
# y_qs_hms <- y_qs[y_qs$species=="hms",]
# y_qs_cps <- y_qs[y_qs$species=="cps",]
# y_qs_gfs <- y_qs[y_qs$species=="gfs",]
# 
# dat2 <- left_join(dat,y_qs_hms,by="esm")
# # geom_vline(aes(xintercept = q_percent))+
# 
# ggplot(data=dat,aes(x=year,y=RMSE, ymin=min,ymax=max))+
#   geom_line(aes(col=EM))+
#   # geom_vline(aes(xintercept = q_percent))+
#   scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
#   geom_ribbon(alpha=0.5)+
#   theme_classic()+
#   labs(x="",y="HMS: Correlation Coefficient")+
#   theme(legend.title = element_blank())+
#   geom_vline(xintercept = 2010, linetype="dashed")+
#   theme(legend.position = "none")+
#   facet_wrap(~esm)+
#   scale_x_continuous(limits = c(1990,2095), expand=c(0,0))







#----FIGURE 4.5.1: Correlation timeseries BY REGION-----
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
  
  #Assign region
  all_esm_allyears$region <- ifelse(all_esm_allyears$lat<=34.5,'South', 
                                    ifelse(all_esm_allyears$lat>=40, "North", "Central"))
  #function
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm, region) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
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
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm",'region', "min", "max"), variable.name = "EM",value.name = "RMSE")
  all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  species_rmse_diff[[counter]]<- all_esm_agg_longer_filtered
  
  # #get threshold 
  # thresh <- melt(all_esm_agg,id=c("year", "esm"), variable.name = "EM",value.name = "Cor")
  # th <- thresh %>% filter(year<=2010) %>% filter(EM!="ens_mean") %>% filter(EM!="min") %>% filter(EM!="max") %>% 
  #   group_by(esm) %>% summarise_at("Cor", funs(mean, sd))
  # th$two_sd <- th$mean - th$sd*2
  # th$three_sd <- th$mean - th$sd*3
  # species_thresh[[counter]] <- th
  
  counter=counter+1
}

ggplot(data=species_rmse_diff[[1]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  # geom_hline(data=species_thresh[[1]], aes(yintercept = three_sd))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="HMS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  theme(legend.position = "none")+
  facet_wrap(~region + esm )

ggplot(data=species_rmse_diff[[2]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  # geom_hline(data=species_thresh[[2]], aes(yintercept = three_sd))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="CPS: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  theme(legend.position = "none")+
  facet_wrap(~region + esm )

ggplot(data=species_rmse_diff[[3]],aes(x=year,y=RMSE, ymin=min,ymax=max))+
  geom_line(aes(col=EM))+
  # geom_hline(data=species_thresh[[3]], aes(yintercept = three_sd))+
  scale_color_manual(values=c("#c0392b"),labels = c("Ensemble mean"))+
  geom_ribbon(alpha=0.5)+
  theme_classic()+
  labs(x="",y="Groundfish: Average Correlation Coefficient")+
  theme(legend.title = element_blank())+
  geom_vline(xintercept = 2010, linetype="dashed")+
  theme(legend.position = "none")+
  facet_wrap(~region + esm )

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



#----FIGURE 4.6.1: subset  of Correlation maps for final pub------
#Plot all three species, for Had only, and for gam_e, es, and est, for 2 temporal periods

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
                                          "gam_E","gam_ES","esm")]
  #Allocate decade 
  all_esm_allyears <- all_esm_allyears[all_esm_allyears$year<=2010 | all_esm_allyears$year>=2075,]
  all_esm_allyears$decade <- ifelse(all_esm_allyears$year<=2010,"1985-2010",
                                                  ifelse(all_esm_allyears$year>=2075, "2075-2100", NA))
  all_esm_allyears$lat_round <- round(all_esm_allyears$lat,0)
  all_esm_allyears$lon_round <- round(all_esm_allyears$lon,0)
  
  #Calculate correlation
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate by ESM and grid cell
  all_esm_agg <-all_esm_allyears %>% group_by(lon, lat, decade, esm) %>% summarise_at(c("gam_E", "gam_ES"),funs(COR(., y=abundance)))
  all_esm_agg_longer <- melt(all_esm_agg,id=c("lat","lon","decade", "esm"), variable.name = "EM",value.name = "COR")
  # all_esm_agg_longer_filtered <- all_esm_agg_longer[all_esm_agg_longer$EM=="ens_mean",]
  all_esm_agg_longer$species <- s
  species_cor_spatial[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

dat <- do.call(rbind,species_cor_spatial)
dat <- na.omit(dat)
dat$species2 <- ifelse(dat$species=="Albacore EMs", "HMS", 
                       ifelse(dat$species=="Anchovy EMs", "CPS", "GFS"))
dat$species2 <- as.factor(dat$species2)
dat$species2 <- factor(dat$species2, levels = c("HMS","CPS","GFS"))
g1 <- ggplot(dat[dat$esm=="HAD",],aes(x=lon,y=lat))+
    geom_tile(aes(fill=COR))+
    scale_fill_gradient2(
      low = 'red', mid = 'light grey', high = 'blue',
      midpoint = 0, guide = 'colourbar', aesthetics = 'fill', limits=c(-1,1))+
    theme_classic()+
    annotation_map(map_data("world"), colour = "black", fill=grey(0.7))+
    coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48))+
    labs(x="",y="")+
    theme(legend.title = element_blank())+
  facet_grid(species2 ~ EM + decade)

tiff(paste0('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig4.6.1.tiff'),
     units="in",width=8, height=8, res=200)
plot(g1)
dev.off()



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
  # had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  # gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  # ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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

#A better figure:
species_predperf[[1]]$species <- "HMS"
species_predperf[[2]]$species <- "CPS"
species_predperf[[3]]$species <- "GFS"
dat <- do.call(rbind,species_predperf)
dat$species <- as.factor(dat$species)
dat$species <- factor(dat$species, levels=c("HMS","CPS","GFS"))

g1 <- ggplot(data=dat,aes(x=EM,y=COR))+
  # geom_bar(stat = "identity")+
  geom_boxplot()+
  geom_point(aes(col=esm, shape=species), size=2, alpha=0.9) +
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  labs(x="",y="Correlation Coefficient: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  # facet_wrap(~species)+
  theme(legend.position = "bottom")

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig7.tiff',res=300,units="in",height=5,width=8)
plot(g1)
dev.off()


#----FIGURE 8: Fitted years Correlation (also RMSE)-------
#Load in predictions made in EM files
species_predperf <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  # had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  # gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  # ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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


#Load in predictions made in EM files
species_predperf <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  # had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  # gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  # ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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

#A better figure:
species_predperf[[1]]$species <- "HMS"
species_predperf[[2]]$species <- "CPS"
species_predperf[[3]]$species <- "GFS"
dat <- do.call(rbind,species_predperf)
dat$species <- as.factor(dat$species)
dat$species <- factor(dat$species, levels=c("HMS","CPS","GFS"))

g1 <- ggplot(data=dat,aes(x=EM,y=COR))+
  # geom_bar(stat = "identity")+
  geom_boxplot()+
  geom_point(aes(col=esm, shape=species), size=2, alpha=0.9) +
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  labs(x="",y="Correlation Coefficient: forecast years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  # facet_wrap(~species)+
  theme(legend.position = "bottom")

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig7.tiff',res=300,units="in",height=5,width=8)
plot(g1)
dev.off()

#A better figure:
species_predperf[[1]]$species <- "HMS"
species_predperf[[2]]$species <- "CPS"
species_predperf[[3]]$species <- "GFS"
dat <- do.call(rbind,species_predperf)
dat$species <- as.factor(dat$species)
dat$species <- factor(dat$species, levels=c("HMS","CPS","GFS"))

g2 <- ggplot(data=dat,aes(x=EM,y=COR))+
  # geom_bar(stat = "identity")+
  geom_boxplot()+
  geom_point(aes(col=esm, shape=species), size=2, alpha=0.9) +
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_classic()+
  labs(x="",y="Correlation Coefficient: fitted years")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  # facet_wrap(~species)+
  theme(legend.position = "none")

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig8.tiff',res=300,units="in",height=8,width=8)
grid.arrange(g2, g1, nrow=2) #g1 comes from Fig 7 code above
dev.off()

#------Figure 8.5-----

#Load in predictions made in EM files
species_predperf_hist <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  # had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  # gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  # ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist)
  gfdl_allyears <- rbind(gfdl_dat_hist)
  ipsl_allyears <- rbind(ipsl_dat_hist)
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
  species_predperf_hist[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

species_predperf_hist[[1]]$species <- "HMS"
species_predperf_hist[[2]]$species <- "CPS"
species_predperf_hist[[3]]$species <- "GFS"
dat_hist <- do.call(rbind,species_predperf_hist)
dat_hist$species <- as.factor(dat_hist$species)
dat_hist$species <- factor(dat_hist$species, levels=c("HMS","CPS","GFS"))



#Load in predictions made in EM files
species_predperf_fcast <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  # had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  # gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  # ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
  species_predperf_fcast[[counter]]<- all_esm_agg_longer
  counter=counter+1
}


#A better figure:
species_predperf_fcast[[1]]$species <- "HMS"
species_predperf_fcast[[2]]$species <- "CPS"
species_predperf_fcast[[3]]$species <- "GFS"
dat_fcast <- do.call(rbind,species_predperf_fcast)
dat_fcast$species <- as.factor(dat_fcast$species)
dat_fcast$species <- factor(dat_fcast$species, levels=c("HMS","CPS","GFS"))


dat_hist$type <- "Fitted Years"
dat_fcast$type <- "Forecast Years"
dat <- rbind(dat_hist, dat_fcast)
dat$ESM <- dat$esm #lazy
dat$Species <- dat$species

dat$EM <- as.factor(dat$EM)
levels(dat$EM) <-  list("brt_E" ="brt_E" ,   "brt_ES" = "brt_ES",  "brt_EST" = "brt_EST",
                           "gam_E"= "gam_E","gam_ECor"="gam_ECor", "gam_ES"= "gam_ES","gam_EST"="gam_EST","gam_S"="gam_S",
                           "glmm_ES"= "glm_E" ,  "glmm_EST" = "glm_ESr",  "glmm_ESTt"="glm_ESt",  "glmm_ST"= "glm_Sr",
                           "mlp_E"="mlp_E",   "mlp_ES" ="mlp_ES" ,  "mlp_EST"="mlp_EST" )
# levels(dat$EM) <-  list("brt_E" ="brt_E" ,   "gam_E"= "gam_E", "mlp_E"="mlp_E",
#                         "brt_ES" = "brt_ES", "gam_ES"= "gam_ES","glmm_ES"= "glm_E","mlp_ES" ="mlp_ES" ,
#                         "brt_EST" = "brt_EST","gam_EST"="gam_EST","glmm_EST" = "glm_ESr","glmm_ESTt"="glm_ESt", "mlp_EST"="mlp_EST" ,
#                         "gam_ECor"="gam_ECor", "gam_S"="gam_S",
#                           "glmm_ST"= "glm_Sr")

g1 <- ggplot(data=dat,aes(x=EM,y=COR))+
  # geom_bar(stat = "identity")+
  geom_boxplot()+
  geom_point(aes(col=ESM, shape=Species), size=2, alpha=0.9) +
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  theme_bw()+
  labs(x="",y="Correlation Coefficient")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  # facet_wrap(~species)+
  theme(legend.position = "bottom")+
  facet_wrap(~type, nrow=2, scales="free")+
  ylim(0.45,0.95)

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig8.5.tiff',res=300,units="in",height=8,width=8)
plot(g1)
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
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
    
    
    N1 <- fitted_samples(gam_EST_N, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_EST_P, n=100, newdata=dat_allyears, scale="response", seed=99)  #response scale doesn't work for binomial model for some reason
    sim_data_EST <- cbind(N1,P1$fitted)
    colnames(sim_data_EST) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_EST$biomass  <- exp(sim_data_EST$bio_sim) * sim_data_EST$pres_sim
    sim_data_EST$year <- dat_allyears$year[sim_data_EST$inputdata_row]
    sim_data_EST$model <- "GAM_ES"
    
    N1 <- fitted_samples(gam_ECor_N$gam, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_ECor_P$gam, n=100, newdata=dat_allyears, scale="response", seed=99)  #response scale doesn't work for binomial model for some reason
    sim_data_ECor <- cbind(N1,P1$fitted)
    colnames(sim_data_ECor) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_ECor$biomass  <- exp(sim_data_ECor$bio_sim) * sim_data_ECor$pres_sim
    sim_data_ECor$year <- dat_allyears$year[sim_data_ECor$inputdata_row]
    sim_data_ECor$model <- "GAM_ECor"
    
    sim_data_esm <- rbind(sim_data_E,sim_data_ES,sim_data_EST,sim_data_ECor)
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


#----Figure 13.5: summary plot of within-model uncertainty-----
#Simplify plot to show three species; only HAD, and four gam models
withinmodel_uncertainty <- list()
species_counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  
  #Load in data & EMS and simulate from posterior distribution
  counter=1
  esm_sims <- as.data.frame(matrix(NA,nrow=1044, ncol=6))
  colnames(esm_sims) <- c("year","esm","model","mean","min","max")
  for (esm in c("had","ipsl","gfdl")){
    
    dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/",esm,"/dat_hist_results_full.rds"))
    dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/",esm,"/dat_fcast_results_full.rds"))
    dat_allyears <- rbind(dat_hist,dat_fcast)
    load(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/",esm,"/saved_models_full.RData"))
    
    #remove unnecessary models (at this stage)
    rm( glm_E_N,glm_E_P, glm_ESr_N,glm_ESr_P, glm_ESt_N, glm_ESt_P, glm_Sr_N, glm_Sr_P,
        brt_E_N, brt_E_P, brt_ES_N, brt_ES_P, brt_EST_N, brt_EST_P,
        mlp_E_N, mlp_E_P, mlp_ES_N, mlp_ES_P, mlp_EST_N, mlp_EST_P)
    
    #Draw fitted values from posterior distribution of GAM
    N1 <- fitted_samples(gam_E_N, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_E_P, n=100, newdata=dat_allyears, scale="response", seed=99)  
    sim_data_E <- cbind(N1,P1$fitted)
    colnames(sim_data_E) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_E$biomass  <- exp(sim_data_E$bio_sim) * sim_data_E$pres_sim
    sim_data_E$year <- dat_allyears$year[sim_data_E$inputdata_row]
    sim_data_E$model <- "GAM_E"
    
    N1 <- fitted_samples(gam_ES_N, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_ES_P, n=100, newdata=dat_allyears, scale="response", seed=99)  
    sim_data_ES <- cbind(N1,P1$fitted)
    colnames(sim_data_ES) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_ES$biomass  <- exp(sim_data_ES$bio_sim) * sim_data_ES$pres_sim
    sim_data_ES$year <- dat_allyears$year[sim_data_ES$inputdata_row]
    sim_data_ES$model <- "GAM_ES"
    
    N1 <- fitted_samples(gam_EST_N, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_EST_P, n=100, newdata=dat_allyears, scale="response", seed=99)  
    sim_data_EST <- cbind(N1,P1$fitted)
    colnames(sim_data_EST) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_EST$biomass  <- exp(sim_data_EST$bio_sim) * sim_data_EST$pres_sim
    sim_data_EST$year <- dat_allyears$year[sim_data_EST$inputdata_row]
    sim_data_EST$model <- "GAM_EST"
    
    N1 <- fitted_samples(gam_ECor_N$gam, n=100, newdata=dat_allyears, scale="response", seed=99)
    P1 <- fitted_samples(gam_ECor_P$gam, n=100, newdata=dat_allyears, scale="response", seed=99)  
    sim_data_ECor <- cbind(N1,P1$fitted)
    colnames(sim_data_ECor) <- c("inputdata_row","draw","bio_sim","pres_sim")
    sim_data_ECor$biomass  <- exp(sim_data_ECor$bio_sim) * sim_data_ECor$pres_sim
    sim_data_ECor$year <- dat_allyears$year[sim_data_ECor$inputdata_row]
    sim_data_ECor$model <- "GAM_ECor"
    
    sim_data_esm <- rbind(sim_data_E,sim_data_ES,sim_data_EST,sim_data_ECor)
    sim_data_esm$esm <- esm
    
    #Get sum of simulated biomass across grid cells
    sim_data_agg <-  sim_data_esm %>% group_by(year, draw, esm, model) %>% summarise_at(c("biomass"), funs(sum)) 
    sim_data_agg <- sim_data_agg %>% group_by(esm, draw, model) %>% summarise_at(vars("biomass","year"), funs(rollmean(.,k=11, fill=NA)))
    #Get mean, min, and max of sims
    sim_data_agg_funs <-  sim_data_agg %>% group_by(year, esm, model) %>% summarise_at(c("biomass"), funs(mean, min, max))
    
    sc <- counter 
    ec <- sc + 427
    esm_sims[sc:ec,] <- sim_data_agg_funs
    counter = counter + 428
  }
  
  withinmodel_uncertainty[[species_counter]] <- esm_sims
  species_counter <- species_counter +1
  
}


#Now get true abundance for plotting
obs_dat <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- rep(c("had","gfdl","ipsl"),each=58000)
  all_esm_allyears_agg <- all_esm_allyears %>% group_by(year, esm) %>% summarise_at(c("abundance"),sum)  
  all_esm_allyears_agg <- all_esm_allyears_agg %>% group_by(esm) %>% summarise_at(vars("abundance","year"), funs(rollmean(.,k=11, fill=NA)))
  all_esm_allyears_agg$sp <- s
  obs_dat[[counter]] <- all_esm_allyears_agg
  counter=counter+1
}

withinmodel_uncertainty[[1]]$species <- "HMS"
withinmodel_uncertainty[[2]]$species <- "CPS"
withinmodel_uncertainty[[3]]$species <- "GFS"
dat <- do.call(rbind,withinmodel_uncertainty)

obs_dat[[1]]$species <- "HMS"
obs_dat[[2]]$species <- "CPS"
obs_dat[[3]]$species <- "GFS"
obs_dat <- do.call(rbind,obs_dat)

dat$species <- as.factor(dat$species)
dat$species <- factor(dat$species, levels=c("HMS","CPS","GFS"))

g1 <- ggplot(data=dat[dat$esm=="had",], aes(x=year,y=mean))+
  geom_line(data=obs_dat, aes(x=year,y=abundance), col="grey")+
  geom_ribbon(aes(ymin=min,ymax=max), alpha=0.5, fill="red")+
  geom_line(aes(y=mean), col="red")+
  theme_classic()+
  facet_wrap(~ species + model, scales = "free", nrow=3)+
  labs(x="",y="Biomass")+
  geom_vline(xintercept = 2010, linetype="dashed")+
  scale_x_continuous(limits = c(1990,2095), expand=c(0,0))


tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig13.5_smoothed.tiff',res=300,units="in",height=8,width=11)
plot(g1)
dev.off()


#-----FIGURE 14: delta maps------
#2x3 plot: 3 species archetypes of delta distribution maps of historical and 2100 period
#as per Fig 7 in Morley et el
#warning: this takes a few mins to plot (Should have used Faet)

delta_output <- list()
counter=1
for (s in c("hms","cps","gfs")){
  #Get data
  if (s =="hms"){
    dat_had <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Albacore EMs/","had","/dat_albacore_all.rds"))
    dat_gfdl <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Albacore EMs/","gfdl","/dat_albacore_all.rds"))
    dat_ipsl <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Albacore EMs/","ipsl","/dat_albacore_all.rds"))
  }
  if(s=="cps"){
    dat_had <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Anchovy EMs/had/dat_anchovy_all.rds"))
    dat_gfdl <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Anchovy EMs/gfdl/dat_anchovy_all.rds"))
    dat_ipsl <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Anchovy EMs/ipsl/dat_anchovy_all.rds"))
  }
  if(s=="gfs"){
    dat_had <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Groundfish EMs/had/dat_groundfish_all.rds"))
    dat_gfdl <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Groundfish EMs/had/dat_groundfish_all.rds"))
    dat_ipsl <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/Groundfish EMs/had/dat_groundfish_all.rds"))
  }
  #combine esm
  dat_had$esm <- "had"
  dat_gfdl$esm <- "gfdl"
  dat_ipsl$esm <- "ipsl"
  dat_allesm <- rbind(dat_had, dat_gfdl, dat_ipsl)
  #isolate years of interest
  dat_allesm_hist <- dat_allesm[dat_allesm$year<=2010,]
  dat_allesm_fcast <- dat_allesm[dat_allesm$year>=2075,]
  #average across years and ensembles
  dat_allesm_hist_agg <- dat_allesm_hist %>% group_by(lat, lon) %>% summarise_at(c("abundance"),c(mean),na.rm=T)
  dat_allesm_fcast_agg <- dat_allesm_fcast %>% group_by(lat, lon) %>% summarise_at(c("abundance"),c(mean),na.rm=T)
  #now merge and get delta: distant future
  dat_allesm_agg_delta <- left_join(dat_allesm_hist_agg,dat_allesm_fcast_agg ,by=c("lat","lon"))
  dat_allesm_agg_delta$delta <- dat_allesm_agg_delta$abundance.y - dat_allesm_agg_delta$abundance.x 

   #saveoutput
  delta_output[[counter]] <- dat_allesm_agg_delta
  delta_output[[counter+1]] <- dat_allesm_hist_agg
  counter = counter+2
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
  scale_fill_gradient2(low = "#e74c3c", mid="#ecf0f1",high = "#3498db",
    # colours = colorRampPalette(colors = c("#e74c3c", "white", "#3498db"))(256),
                       limits = c(min(delta_output[[1]]$delta), max(delta_output[[1]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-134,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="HMS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2075-2100", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass (kg)", size=3)

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
  geom_text(x=-119.5, y=42.5, label="Biomass (kg)", size=3)

#CPS
cps_fcast <- ggplot(data = delta_output[[3]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=delta))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradient2(low = "#e74c3c", mid="#ecf0f1",high = "#3498db",
                       limits = c(min(delta_output[[3]]$delta), max(delta_output[[3]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="CPS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2075-2100", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass (kg)", size=3)

cps_hist <- ggplot(data = delta_output[[4]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=abundance))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c"))(256),
                       limits = c(min(delta_output[[4]]$abundance), max(delta_output[[4]]$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="CPS Archetype", size=5)+
  geom_text(x=-120, y=44, label="1985-2010", size=3)+
  geom_text(x=-119.5, y=42.5, label="Biomass (kg)", size=3)

#REPEAT FOR GROUNDFISH
gfs_fcast <- ggplot(data = delta_output[[5]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=delta))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradient2(low = "#e74c3c", mid="#ecf0f1",high = "#3498db",
                       limits = c(min(delta_output[[5]]$delta), max(delta_output[[5]]$delta))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="GFS Archetype", size=5)+
  geom_text(x=-120, y=44, label="2075-2100", size=3)+
  geom_text(x=-119.5, y=42.5, label="Delta Biomass (kg)", size=3)

gfs_hist <- ggplot(data = delta_output[[6]], aes(x=lon,y=lat))+
  geom_tile(aes(fill=abundance))+
  theme_classic() +  labs(y="", x="") + 
  theme(legend.position=c(0.6,0.6),
        legend.direction = "horizontal",
        legend.background = element_rect(fill="grey50"),
        legend.title = element_blank())+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1)) + #makes a box
  scale_fill_gradientn(colours = colorRampPalette(colors = c("#9b59b6","#3498db", "#1abc9c", "#f1c40f", "#e74c3c"))(256),
                       limits = c(min(delta_output[[6]]$abundance), max(delta_output[[6]]$abundance))) +
  annotation_map(map_data("world"), colour = "black", fill="grey50")+
  coord_quickmap(xlim=c(-126,-115.8),ylim=c(30,48)) +  #Sets aspect ratio
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  geom_text(x=-120, y=45, label="GFS Archetype", size=5)+
  geom_text(x=-120, y=44, label="1985-2010", size=3)+
  geom_text(x=-119.5, y=42.5, label="Biomass (kg)", size=3)



tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig14.tiff',res=300,units="in",height=10,width=12)
# grid.arrange(hms_hist, cps_hist, gfs_hist,
#              # hms_mid, cps_mid, gfs_mid,
#              hms_fcast, cps_fcast, gfs_fcast, nrow=2)

h = 0.33
w = 0.5
ggdraw() +
  draw_plot(hms_hist, 0, 0.5, h, w) +
  draw_plot(cps_hist, 0.33, 0.5, h, w) +
  draw_plot(gfs_hist, 0.6, 0.5, h, w) +
  
  draw_plot(hms_fcast, 0, 0, h, w) +
  draw_plot(cps_fcast, 0.33, 0, h, w) +
  draw_plot(gfs_fcast, 0.6, 0, h, w)
  
dev.off()


#rose plots
# ggplot(d, aes(x = Angle, y = Frequency)) +
#   coord_polar(theta = "x", start = -pi/45) +
#   geom_bar(stat = "identity") +
#   scale_x_continuous(breaks = seq(0, 360, 60))

#------Figure 14.5------
#Time series of environmental novelty
#Novelty File
extrap <- readRDS('~/Dropbox/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet_PlusNearby.rds')

#ggplot for paper
#prepare data LAT cogs:
#hms
hms1 <- extrap[[1]]; hms2 <- extrap[[2]]; hms3 <- extrap[[3]]
hms1$esm = "HAD"; hms2$esm = "GFDL"; hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
#cps
cps1 <- extrap[[4]];cps2 <- extrap[[5]];cps3 <- extrap[[6]]
cps1$esm = "HAD";cps2$esm = "GFDL";cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
#gfs
gfs1 <- extrap[[7]];gfs2 <- extrap[[8]];gfs3 <- extrap[[9]]
gfs1$esm = "HAD";gfs2$esm = "GFDL";gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)

hms <- hms %>% group_by(esm) %>% summarise_at(vars("extrap_percent", "year"), funs(rollmean(., k=11, fill=NA)))
cps <- cps %>% group_by(esm) %>% summarise_at(vars("extrap_percent", "year"), funs(rollmean(., k=11, fill=NA)))
gfs <- gfs %>% group_by(esm) %>% summarise_at(vars("extrap_percent", "year"), funs(rollmean(., k=11, fill=NA)))
hms$species <- "HMS"
cps$species <- "CPS"
gfs$species <- "GFS"
data <- rbind(hms, cps, gfs)
data$species <- factor(data$species, levels = c("HMS","CPS","GFS"))
g1 <- ggplot(data=data, aes(x=year, y=extrap_percent, group=esm))+
  geom_line(aes(col=esm),size=1)+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  facet_wrap(~species)+
  ylab("Percent of environmental novelty")+xlab("")+
  theme_bw()+
  theme(legend.position = "bottom", legend.title=element_blank())


tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig14.5.tiff', res=300, units="in",height=4,width=10)
plot(g1)
dev.off()


#-----FIGURE 15: Partitioning uncertainty with dominance analysis and BIOMASS-----
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
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(sum(.)),na.rm=T)
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
    dat <- dat_cog[dat_cog$year==y,]
    m2 <- lm(COG ~ esm + type + params, data = dat) #must be class lm for dominance analysis
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

da_allyears_allspecies_fullmodels <- da_allyears_allspecies

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


#Any correlations between scenario uncertainty and variability OR magnitude? 
test <- da_allyears_allspecies
test <- do.call(rbind,test)
test$species <- rep(c("HMS","CPS","GFS"),each=348)
test <- test[test$source=="esm" & test$year>=2011,]
plot(test$year,test$percent, col=as.factor(test$species)) #90 points

extrap <- readRDS('~/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet_PlusNearby.rds')
#hms
hms1 <- extrap[[1]]; hms2 <- extrap[[2]]; hms3 <- extrap[[3]]
hms1$esm = "HAD"; hms2$esm = "GFDL"; hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms$species <- "HMS"
#cps
cps1 <- extrap[[4]];cps2 <- extrap[[5]];cps3 <- extrap[[6]]
cps1$esm = "HAD";cps2$esm = "GFDL";cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps$species <- "CPS"
#gfs
gfs1 <- extrap[[7]];gfs2 <- extrap[[8]];gfs3 <- extrap[[9]]
gfs1$esm = "HAD";gfs2$esm = "GFDL";gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs$species <- "GFS"
full_dat<- rbind(hms,cps,gfs)

full_dat <- full_dat %>% select(year, extrap_percent, esm, species) %>% 
  group_by(year,species) %>% summarise(cv = sd(extrap_percent,na.rm=T)/mean(extrap_percent,na.rm=T))
full_dat <- full_dat[full_dat$year>=2011,]



par(mfrow=c(1,3))
for (s in c("HMS","CPS","GFS")){
  extrap_cv <- full_dat$cv[full_dat$species==s]
  esm_uncertainty <- test$percent[test$source=="esm" & test$species==s]
  df <- as.data.frame(cbind(extrap_cv,esm_uncertainty))
  plot(df$extrap_cv, df$esm_uncertainty)
  cor.test(df$extrap_cv, df$esm_uncertainty)
  m <- lm(esm_uncertainty ~ extrap_cv, data=df)
  abline(m)
}

#when CVS are low, so should be esm



#Ok so uncertainty is a reflection of environmental variability among ESMs, and extent of extrapolation
#So lets plot CV of extrapolation through time and get breakpoints. 
#Then get mean CV of extrapolation at each breakpoint and we have our threshold? 

#plot CV through time
plot(full_dat$year,full_dat$sd, col=as.factor(full_dat$species))
plot(full_dat$year,full_dat$mean, col=as.factor(full_dat$species))
plot(full_dat$year,full_dat$cv, col=as.factor(full_dat$species))

#now get breakpoints for each species. #DOESN'T WORK. SKIP anduse oterh break points
# dat_temp <- full_dat
# dat_temp_long <- dat_temp %>% pivot_wider(names_from = species, values_from = cv)
# dat_mat <- as.matrix(dat_temp_long)
# fit_cpm1 = e.cp3o(dat_mat, minsize=5, K=5)
# print(fit_cpm1$estimates) #1


dat_temp <- full_dat
dat_temp$breaks <- ifelse(dat_temp$year<=2044,1,
                          ifelse(dat_temp$year>2077,3,2))
dat_temp_agg <- dat_temp %>% group_by(species, breaks) %>% summarise(mean(cv, na.rm=T))
par(mfrow=c(1,3))
barplot(dat_temp_agg$`mean(cv, na.rm = T)`[dat_temp_agg$species=="HMS"])
barplot(dat_temp_agg$`mean(cv, na.rm = T)`[dat_temp_agg$species=="CPS"])
barplot(dat_temp_agg$`mean(cv, na.rm = T)`[dat_temp_agg$species=="GFS"])




#-----Figure 15.5: changepoint analysis on uncertainty time-series-----
#Fist get uncertainty partitioned via dominance anlaysis, then calculate time-series change points

species_biom_diff <- list()
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
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(sum(.)),na.rm=T)
  all_esm_agg_longer <- all_esm_agg %>% pivot_longer(cols=c(3:14),names_to="EM",values_to="COG")
  species_biom_diff[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

#Now run dominance analysis for every year & archetype
da_allyears_allspecies <- list()
upper_counter <- 1
for (s in 1:3){
  dat_cog <- species_biom_diff[[s]]
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
    dat <- dat_cog[dat_cog$year==y,]
    m2 <- lm(COG ~ esm + type + params, data = dat) #must be class lm for dominance analysis
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

dat_da <- do.call(rbind,da_allyears_allspecies)
dat_da$species <- rep(c("HMS","CPS","GFS"),each=348)

# library(cpm) #UNIVARIATE APPROACH
# for (s in c("HMS","CPS","GFS")){
#   dat_temp <- dat_da[dat_da$species==s & dat_da$source=="params",]
#   fit_cpm = processStream(dat_temp$percent, cpmType = "Student")
#   print(fit_cpm$changePoints)
# }

#RUN ON YEARS >2011
for (s in c("HMS","CPS","GFS")){
  print(glue("running species {s}"))
  dat_temp <- dat_da[dat_da$species==s & dat_da$year>=2011,]
  dat_temp_long <- dat_temp %>% dplyr::select(year, source, percent) %>% pivot_wider(., names_from = "source", values_from = "percent")
  dat_mat <- as.matrix(dat_temp_long)
  # fit_cpm = e.cp3o(dat_mat)
  # print(fit_cpm$estimates)
  fit_cpm1 = e.cp3o(dat_mat, minsize=5, K=5)
  # fit_cpm2 = e.cp3o_delta(dat_mat, delta=20, K=3)
  # fit_cpm3 = e.divisive(dat_mat, min.size=20)
  # fit_cpm6 = e.agglo(dat_mat, alpha=1)
  print(fit_cpm1$estimates)
  # print(fit_cpm2$estimates)
  # print(fit_cpm3$estimates)
  # print(fit_cpm6$estimates)
}

# dat_mat[c(42,84)] #hms
# dat_mat[c(43,86)] #cps
# dat_mat[c(44,87)] #gfs

dat_mat[c(33,66)] #hms
dat_mat[c(34,67)] #cps
dat_mat[c(34,67)] #gfs

for (s in c("HMS","CPS","GFS")){
  print(glue("running species {s}"))
  print(mean(dat_da$percent[dat_da$source=="esm" & dat_da$species==s & dat_da$year<2027], na.rm=T))
  print(mean(dat_da$percent[dat_da$source=="esm" & dat_da$species==s & dat_da$year>2027 & dat_da<2070], na.rm=T))
  print(mean(dat_da$percent[dat_da$source=="esm" & dat_da$species==s & dat_da$year>2070], na.rm=T))
}

#Correlate year 59 breakpoint to environmental extrapolation? 
#what about running this breakpoints analysis on enviro timeseriesdata? Would it show the same thing? 
# 
extrap <- readRDS('~/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet_PlusNearby.rds')
#hms
hms1 <- extrap[[1]]; hms2 <- extrap[[2]]; hms3 <- extrap[[3]]
hms1$esm = "HAD"; hms2$esm = "GFDL"; hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms$species <- "HMS"
#cps
cps1 <- extrap[[4]];cps2 <- extrap[[5]];cps3 <- extrap[[6]]
cps1$esm = "HAD";cps2$esm = "GFDL";cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps$species <- "CPS"
#gfs
gfs1 <- extrap[[7]];gfs2 <- extrap[[8]];gfs3 <- extrap[[9]]
gfs1$esm = "HAD";gfs2$esm = "GFDL";gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs$species <- "GFS"
full_dat<- rbind(hms,cps,gfs)

mean(full_dat$extrap_percent[full_dat$year==2100])
mean(full_dat$extrap_percent[full_dat$year==2100 & full_dat$species=="GFS"])

for (s in c("HMS","CPS","GFS")){
  # for (e in c("HAD","GFDL","IPSL")){
    print(glue("running species {s}"))
    # print(glue("running esm {e}"))
    dat_temp <- full_dat[full_dat$species==s  ,]
    dat_temp_long <- dat_temp %>% dplyr::select(year, extrap_percent, esm) %>% pivot_wider(., names_from = "esm", values_from = "extrap_percent")
    dat_mat <- as.matrix(dat_temp_long)
    fit_cpm1 = e.cp3o(dat_mat, minsize=5, K=5)
    print(fit_cpm1$estimates)
  # }
}

dat_temp_long$year[c(28,56,82)]
dat_temp_long$year[c(13,38,64,83)]
dat_temp_long$year[c(32,64)]



# ggplot(data=full_dat, aes(x=year, y=extrap_percent))+
#   geom_point(aes(col=esm))+
#   scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
#   geom_vline(xintercept = 2026, linetype="dotted")+
#   geom_vline(xintercept = 2070 ,linetype="dotted")+
#   theme_bw()+
#   theme(legend.position = "bottom", legend.title = element_blank())+
#   labs(x="Year",y="Habitat extrapolation from historical conditions (%)")
#   facet_wrap(~species , scales="free")


plot(hms$year,hms$extrap_percent, col=as.factor(hms$esm), pch=19)
# abline(v=2015); abline(v=2044); abline(v=2073)
abline(v=2026, col="red"); abline(v=2070, col="red")

plot(cps$year,cps$extrap_percent, col=as.factor(cps$esm), pch=19)
# abline(v=2015); abline(v=2044); abline(v=2073)
abline(v=2026, col="red"); abline(v=2070, col="red")

plot(gfs$year,gfs$extrap_percent, col=as.factor(gfs$esm), pch=19)
# abline(v=2015); abline(v=2044); abline(v=2073)
abline(v=2026, col="red"); abline(v=2070, col="red")

#standard deviation / mean
sd(hms$extrap_percent[hms$year<2026])/mean(hms$extrap_percent[hms$year<2026])
sd(hms$extrap_percent[hms$year>2026 & hms$year<2070])/mean(hms$extrap_percent[hms$year>2026& hms$year<2070])
sd(hms$extrap_percent[hms$year<2070])/mean(hms$extrap_percent[hms$year<2070])

sd(cps$extrap_percent[cps$year<2026])/mean(cps$extrap_percent[cps$year<2026])
sd(cps$extrap_percent[cps$year>2026 & cps$year<2070])/mean(cps$extrap_percent[cps$year>2026& cps$year<2070])
sd(cps$extrap_percent[cps$year<2070])/mean(cps$extrap_percent[cps$year<2070])

sd(gfs$extrap_percent[gfs$year<2026])/mean(gfs$extrap_percent[gfs$year<2026])
sd(gfs$extrap_percent[gfs$year>2026 & gfs$year<2070])/mean(gfs$extrap_percent[gfs$year>2026& gfs$year<2070])
sd(gfs$extrap_percent[gfs$year<2070])/mean(gfs$extrap_percent[gfs$year<2070])

#-----FIGURE 16: dominance analysis by region---------
# regions: South 30-34.5; north 40-48; central 34.5-40
#For biomass not COG though. Make it relative to historical period? 

#Let's try it: 
species_biomass <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
      dat <- dat_bio[dat_bio$year==y & dat_bio$region==r,]
      m2 <- lm(TotalBiomass ~ esm + type + params, data = dat ) #must be class lm for dominance analysis
      da <- dominanceAnalysis(m2)
      # print(da)
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
  labs(y="Relative Influence on Biomass Uncertainty", x="Year") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  facet_grid(region~species)+
  theme(legend.title = element_blank(),
        legend.position = "bottom")
a1 <- ggplot()+
  annotation_map(map_data("world"), colour = "black", fill="light grey")+
  coord_quickmap(xlim=c(-126,-115.5),ylim=c(30,48)) +  #Sets aspect ratio
  geom_vline(xintercept=-115.5, col="grey50", size=0.01)+ #add to force xaxis
  geom_hline(yintercept = 40, linetype="dashed")+
  geom_hline(yintercept = 34.5, linetype="dashed")+
  xlab("") + ylab("")+
  theme_classic()+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  ggtitle("Regions")+
  annotate(geom="text", x=-120.5, y=45, label="North", hjust=0, fontface=2, cex=5)+
  annotate(geom="text", x=-120.5, y=38, label="Central", hjust=0, fontface=2, cex=5)+
  annotate(geom="text", x=-120.5, y=32, label="South", hjust=0, fontface=2, cex=5)
#   
# tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig16_partA.tiff', units="in", res=300, height=8, width=10)
# plot(g1)
# dev.off()
# 
# tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig16_partB.tiff', units="in", res=300, height=6, width=3)
# plot(a1)
# dev.off()

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig16.tiff', units="in", res=300, height=6, width=13)
plot_grid(g1, a1, nrow = 1, rel_widths = c(0.74,0.26), axis="t", align="v")
dev.off()

#SMOOTH AND RUN DOMINANCE ANALYSIS ON SMOOTHED DATA
#Now run dominance analysis for every year & archetype
da_allyears_allspecies <- list()
upper_counter <- 1
for (s in 1:3){
  dat_bio <- species_biomass[[s]]
  dat_bio <- dat_bio %>%  separate(EM, c("type","params"))
  dat_bio$esm <- as.factor(dat_bio$esm)
  dat_bio$type <- as.factor(dat_bio$type)
  dat_bio$params <- as.factor(dat_bio$params)
  #SMOOTHER
  dat_bio <- dat_bio %>% group_by(esm, region, type, params) %>% summarise_at(vars("TotalBiomass","year"), funs(rollmean(.,k=11, fill=NA)))
  
  #Anova to look at where significant differences lie:
  # m1 <- aov(COG ~ esm + type + params, data = dat_cog) #must be aov to run TukeyHSD
  # summary(m1)
  # TukeyHSD(m1)
  
  for (r in c("South", "Central","North")){
    #Run dominance analysis for every year
    da_allyears <- as.data.frame(matrix(NA, nrow=348, ncol=6))
    colnames(da_allyears) <- c("r2", "source", "year", "percent", "region", "species")
    counter=1
    for (y in 2000:2085){
      print(y)
      dat <- dat_bio[dat_bio$year==y & dat_bio$region==r,]
      m2 <- lm(TotalBiomass ~ esm + type + params, data = dat ) #must be class lm for dominance analysis
      da <- dominanceAnalysis(m2)
      # print(da)
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

dat <- do.call(rbind, da_allyears_allspecies)
dat$species <- ifelse(dat$species==1,"HMS Archetype",
                      ifelse(dat$species==2,"CPS Archetype", "GFS Archetype"))
dat$species <- factor(dat$species, levels = c("HMS Archetype", "CPS Archetype", "GFS Archetype"))
dat$region <- factor(dat$region, levels = c("North","Central","South"))
dat <- na.omit(dat)
g1 <- ggplot(data=dat, aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_classic()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Biomass Uncertainty", x="Year") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(xintercept = i1[1], linetype="solid", col="red")+
  facet_grid(region~species)+
  theme(legend.title = element_blank(),
        legend.position = "bottom")

a1 <- ggplot()+
  annotation_map(map_data("world"), colour = "black", fill="light grey")+
  coord_quickmap(xlim=c(-126,-115.5),ylim=c(30,48)) +  #Sets aspect ratio
  geom_vline(xintercept=-115.5, col="grey50", size=0.01)+ #add to force xaxis
  geom_hline(yintercept = 40, linetype="dashed")+
  geom_hline(yintercept = 34.5, linetype="dashed")+
  xlab("") + ylab("")+
  theme_classic()+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+  #makes a box
  scale_x_continuous(expand = c(0, 0)) +  scale_y_continuous(expand = c(0, 0))+
  ggtitle("Regions")+
  annotate(geom="text", x=-120.5, y=45, label="North", hjust=0, fontface=2, cex=5)+
  annotate(geom="text", x=-120.5, y=38, label="Central", hjust=0, fontface=2, cex=5)+
  annotate(geom="text", x=-120.5, y=32, label="South", hjust=0, fontface=2, cex=5)

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig16_smoothed.tiff', units="in", res=300, height=6, width=13)
plot_grid(g1, a1, nrow = 1, rel_widths = c(0.74,0.26), axis="t", align="v")
dev.off()


# 
# for (s in c("HMS Archetype","CPS Archetype","GFS Archetype")){
#   for (r in c("North","Central","South")){
#     print(glue("running species {s}"))
#   dat_temp <- dat[dat$species==s &  dat$region==r,]
#   dat_temp_long <- dat_temp %>% select(year, source, percent) %>% pivot_wider(., names_from = "source", values_from = "percent")
#   dat_mat <- as.matrix(dat_temp_long)
#   fit_cpm1 = e.cp3o(dat_mat, K=5, minsize = 5)
#   # fit_cpm2 = e.cp3o_delta(dat_mat, K=3)
#   # fit_cpm3 = e.divisive(dat_mat, min.size=20)
#   # fit_cpm4 = ks.cp3o(dat_mat, K=3)
#   # fit_cpm5 = ks.cp3o_delta(dat_mat, K=3)
#   # fit_cpm6 = e.agglo(dat_mat)
#   print(fit_cpm1$estimates)
#   # print(fit_cpm2$estimates)
#   # print(fit_cpm3$estimates)
#   # print(fit_cpm4$estimates)
#   # print(fit_cpm5$estimates)
#   # print(fit_cpm6$estimates)
#     # print(fit_cpm)
#   }
# }




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
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(sum(.)),na.rm=T)
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
    dat <-  dat_cog[dat_cog$year==y,]
    m2 <- lm(COG ~ esm + type + params, data =dat) #must be class lm for dominance analysis
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

da_allyears_allspecies_temponly <- da_allyears_allspecies


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

dat_da <- do.call(rbind,da_allyears_allspecies_temponly)
dat_da$species <- rep(c("HMS","CPS","GFS"),each=348)

for (s in c("HMS","CPS","GFS")){
  print(glue("running species {s}"))
  dat_temp <- dat_da[dat_da$species==s & dat_da$year>=2011 ,]
  dat_temp_long <- dat_temp %>% dplyr::select(year, source, percent) %>% pivot_wider(., names_from = "source", values_from = "percent")
  dat_mat <- as.matrix(dat_temp_long)
  # fit_cpm = e.cp3o(dat_mat)
  # print(fit_cpm$estimates)
  fit_cpm1 = e.cp3o(dat_mat, minsize=5, K=5)
  # fit_cpm3 = e.divisive(dat_mat, min.size=20)
  # fit_cpm6 = e.agglo(dat_mat, alpha=1)
  print(fit_cpm1$estimates)
  # print(fit_cpm3$estimates)
  # print(fit_cpm6$estimates)
}

dat_mat[c(34,67)] #2028 & 2071
dat_mat[c(33,66)] #2028 & 2071
dat_mat[c(34,67)] #2028 & 2071


#-----Figure 18.5:-----
#Mainly copies code from Figs 15 and 18

#Full models
#Get COG for all species and EMs
species_cog_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(sum(.)),na.rm=T)
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
  #SMOOTHER
  dat_cog <- dat_cog %>% group_by(esm, type, params) %>% summarise_at(vars("COG","year"), funs(rollmean(.,k=11, fill=NA))) #61
  
  #Run dominance analysis for every year
  da_allyears <- as.data.frame(matrix(NA, nrow=348, ncol=4))
  colnames(da_allyears) <- c("r2", "source", "year", "percent")
  counter=1
  for (y in 1990:2095){
    print(y)
    dat <- dat_cog[dat_cog$year==y,]
    m2 <- lm(COG ~ esm + type + params, data = dat) #must be class lm for dominance analysis
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

da_allyears_allspecies_fullmodels <- da_allyears_allspecies


#########Now run temp-only model code###

#Get COG for all species and EMs
species_cog_diff <- list()
counter=1
for(s in c("Albacore EMs TempOnly","Anchovy EMs TempOnly","Groundfish EMs TempOnly")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(sum(.)),na.rm=T)
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
  #SMOOTHER
  dat_cog <- dat_cog %>% group_by(esm, type, params) %>% summarise_at(vars("COG","year"), funs(rollmean(.,k=11, fill=NA))) #61
  
  
  #Run dominance analysis for every year
  da_allyears <- as.data.frame(matrix(NA, nrow=348, ncol=4))
  colnames(da_allyears) <- c("r2", "source", "year", "percent")
  counter=1
  for (y in 1990:2095){
    print(y)
    dat <-  dat_cog[dat_cog$year==y,]
    m2 <- lm(COG ~ esm + type + params, data =dat) #must be class lm for dominance analysis
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

da_allyears_allspecies_temponly <- da_allyears_allspecies

#Data from fig 15 (now above)
da_allyears_allspecies_fullmodels
#Data from fig 18 (now above)
da_allyears_allspecies_temponly

d_full <- do.call(rbind,da_allyears_allspecies_fullmodels)
d_temp <- do.call(rbind,da_allyears_allspecies_temponly)
d_all <- rbind(d_full,d_temp)
d_all$species <- rep(c("HMS Archetype","CPS Archetype","GFS Archetype",
                       "HMS Archetype Temp. Only","CPS Archetype Temp. Only","GFS Archetype Temp. Only"),each=348)
d_all$species <- as.factor(d_all$species)
levels(d_all$species)
d_all$species <- factor(d_all$species, levels = c("HMS Archetype","CPS Archetype","GFS Archetype",
                                                  "HMS Archetype Temp. Only","CPS Archetype Temp. Only","GFS Archetype Temp. Only"))

g1 <- ggplot(data=d_all, aes(x=year, y=percent, fill=source))+
  geom_area()+
  scale_fill_manual(labels = c("Earth System Models", "SDM Types", "SDM Parameters"),
                    values = c("#c7ecee", "#95afc0", "#30336b"))+
  theme_bw()+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y="Relative Influence on Biomass Uncertainty", x="") + 
  geom_vline(xintercept = 2010, linetype="dashed")+
  # geom_vline(aes(xintercept = changepoint1) , linetype="dotted")+
  # geom_vline(aes(xintercept = changepoint2), linetype="dotted")+
  # geom_hline(aes(yintercept = 0.33), linetype="dotted")+
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        plot.margin = margin(10, 15, 10, 10))+
  facet_wrap(~species, scales="free")

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig18.5_smoothed.tiff', units="in",res=300,width=9,height = 5)
plot(g1)
dev.off()


library(cpm)
for (s in c("HMS Archetype","CPS Archetype","GFS Archetype",
            "HMS Archetype Temp. Only","CPS Archetype Temp. Only","GFS Archetype Temp. Only")){
  dat_temp <- d_all[d_all$species==s & d_all$source=="params",]
  fit_cpm = processStream(dat_temp$percent, cpmType = "Student") 
  print(fit_cpm$changePoints)
}

library(ecp)
for (s in c("HMS Archetype","CPS Archetype","GFS Archetype",
            "HMS Archetype Temp. Only","CPS Archetype Temp. Only","GFS Archetype Temp. Only")){
  dat_temp <- d_all[d_all$species==s ,]
  dat_temp_long <- dat_temp %>% select(year, source, percent) %>% pivot_wider(., names_from = "source", values_from = "percent")
  dat_mat <- as.matrix(dat_temp_long)
  fit_cpm = e.cp3o(dat_mat)
  print(fit_cpm$estimates)
  print(fit_cpm)
}



#How much change in first 30 and 80 years?
hms <- da_allyears_allspecies[[1]]
cps <- da_allyears_allspecies[[2]]
gfs <- da_allyears_allspecies[[3]]
d <- rbind(hms,cps,gfs)
mean(hms$percent[hms$source=="esm" & hms$year>=2025 & hms$year<=2050 ], na.rm=T)
mean(cps$percent[cps$source=="esm" & cps$year>=2025 & cps$year<=2050 ], na.rm=T)
mean(hms$percent[gfs$source=="esm" & gfs$year>=2025 & gfs$year>=2050 ], na.rm=T)

mean(d$percent[d$source=="esm" & d$year>=2025 & d$year<=2050 ], na.rm=T)
mean(d$percent[d$source=="esm" & d$year>=2075 & d$year<=2100 ], na.rm=T)

mean(d$percent[d$source=="esm" & d$year>=2025 & d$year<=2050 ], na.rm=T)
mean(d$percent[d$source=="esm" & d$year>=2075 & d$year<=2100 ], na.rm=T)



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
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(sum(.)),na.rm=T)
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
    dat <- dat_cog[dat_cog$year==y,]
    m2 <- lm(COG ~ esm + type + params, data = dat) #must be class lm for dominance analysis
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


#------FIGURE 22 OLD: explore r2 and environmental anomalies: OUTDATEd-----
#Get R2 for each species and each year
#Load in predictions made in EM files
#Calculate anomalies from DSMExtra
#NOTE: NOW OUTDATED FROM JAMES HYPERVOLUME ANALYSIS
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



#------FIGURE 22: hypervolume and r2 correlations-----
#extrapolation determined from hypervolume code. Computed in another script and output loaded here

#Get correlations
species_cors <- list()
counter_large=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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

  all_esm_agg_enviro_longer <- all_esm_agg %>% pivot_longer(c(3:15),names_to="EM", values_to="R2")
  all_esm_agg_enviro_longer$species <- s
  species_cors[[counter_large]]<- all_esm_agg_enviro_longer
  counter_large=counter_large+1
}


extrap <- readRDS('~/Dropbox/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet_PlusNearby.rds')

#Prepare files for merge
#hms
hms1 <- extrap[[1]]
hms2 <- extrap[[2]]
hms3 <- extrap[[3]]
hms1$esm = "HAD"
hms2$esm = "GFDL"
hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms_cor <- rbind(species_cors[[1]])
hms_comp <- left_join(hms_cor,hms, by =c("year","esm"))

#cps
cps1 <- extrap[[4]]
cps2 <- extrap[[5]]
cps3 <- extrap[[6]]
cps1$esm = "HAD"
cps2$esm = "GFDL"
cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps_cor <- rbind(species_cors[[2]])
cps_comp <- left_join(cps_cor,cps, by =c("year","esm"))

#gfs
gfs1 <- extrap[[7]]
gfs2 <- extrap[[8]]
gfs3 <- extrap[[9]]
gfs1$esm = "HAD"
gfs2$esm = "GFDL"
gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs_cor <- rbind(species_cors[[3]])
gfs_comp <- left_join(gfs_cor,gfs, by =c("year","esm"))

#Plots
ggplot(data=hms_comp, aes(x=nearby_mean,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=hms_comp, aes(x=nearby_median,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=hms_comp, aes(x=extrap_percent,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=hms_comp, aes(x=excluded,y=R2, col=esm))+
  geom_point(aes(col=esm))

ggplot(data=cps_comp, aes(x=nearby_mean,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=cps_comp, aes(x=nearby_median,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=cps_comp, aes(x=extrap_percent,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=cps_comp, aes(x=excluded,y=R2, col=esm))+
  geom_point(aes(col=esm))

ggplot(data=gfs_comp, aes(x=nearby_mean,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=gfs_comp, aes(x=nearby_median,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=gfs_comp, aes(x=extrap_percent,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=gfs_comp, aes(x=excluded,y=R2, col=esm))+
  geom_point(aes(col=esm))



#How does each EM perform with novelty?
output <- as.data.frame(matrix(NA, nrow=117, ncol=5))
colnames(output) <- c("species","esm","EM","cor.n", "sig")
counter=1
for (s in c("hms","cps",'gfs')){
  for (e in c("HAD","GFDL","IPSL")){
    for (m in c("gam_E", "gam_ES", "gam_ECor",
                "glm_E" ,"glm_ESt", "glm_ESr",   
                "brt_E" , "brt_ES", "brt_EST" ,  
                "mlp_E" , "mlp_ES" ,"mlp_EST","ens_mean")){
      
      if (s =="hms"){dat = hms_comp}
      if (s =="cps"){dat = cps_comp}
      if (s =="gfs"){dat = gfs_comp}
      
      ch <- dat[dat$esm==e & dat$EM==m,] #get relevant data
      
      cor.n <- cor(ch$R2,ch$extrap_percent, use="na.or.complete")
      sig <- cor.test(ch$R2,ch$extrap_percent, use="na.or.complete")
      
      # near <- cor(ch$R2,ch$percent_exlcluded)
      # first_year <- ch$year[ch$nearby_sum<=near]
      
      output[counter,1] <- s
      output[counter,2] <- e
      output[counter,3] <- m
      output[counter,4] <- cor.n
      output[counter,5] <- sig$p.value
      counter=counter+1
    }
  }
}
boxplot(cor.n~EM ,data=output)
abline(h=-0.21)

output$EM <- as.factor(output$EM)
levels(output$EM) <-  list("brt_E" ="brt_E" ,   "brt_ES" = "brt_ES",  "brt_EST" = "brt_EST",
                                      "ens_mean" = "ens_mean", "gam_E"= "gam_E","gam_ECor"="gam_ECor", "gam_ES"= "gam_ES",
                                      "glmm_E"= "glm_E" ,  "glmm_ESr" = "glm_ESr",  "glmm_ESt"="glm_ESt", 
                                      "mlp_E"="mlp_E",   "mlp_ES" ="mlp_ES" ,  "mlp_EST"="mlp_EST" )

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig22.tiff', units="in",res=300,width=12,height=8)
ggplot(data=output, aes(x=EM, y=cor.n)) +
  geom_boxplot(fill="grey") +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_point(aes(col=esm, shape=species), size=3, alpha=0.9) +
  theme_bw() +
  xlab("Estimation models")+
  ylab("Correlation between model fit and degree of extrapolation")+
  geom_hline(yintercept=-0.21,linetype='dashed' )+
  theme(legend.position="bottom")+
  # scale_colour_manual(values = c("#d1495b","#00798c","#edae49"))+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  labs(col="Earth System Model", shape="Species")
dev.off()


#is there a global threshold we can calculate?
#Global wont' work becasue species $ nearby changes so much
# and then teste against experiments?
#ON HOLD FOR THE MOMENT --> NEED TO DETERMIEN A METRIC THAT SAYS WHEN A MODEL FAILS

hms_comp$species <- "hms"
cps_comp$species <- "cps"
gfs_comp$species <- "gfs"
global <- rbind(hms_comp, cps_comp, gfs_comp)

m1 <- gam(R2 ~ s(nearby_mean) + species + esm + EM + year, data = global)
summary(m1)
plot(m1,rug=T)
gam.check(m1)
  
#Changepoint Analysis for %nearby
output <- as.data.frame(matrix(NA, nrow=100, ncol=5))
colnames(output) <- c("species","esm","EM","percent_nearby")
counter=1
for (e in c("HAD","GFDL","IPSL")){
  for (m in c("gam_E", "gam_ES", "gam_ECor",
              "glm_E" ,"glm_ESt", "glm_ESr",   
              "brt_E" , "brt_ES", "brt_EST" ,  
              "mlp_E" , "mlp_ES" ,"mlp_EST","ens_mean")){
  
  ch <- cps_comp[cps_comp$esm==e & cps_comp$EM==m,] #get relevant data
  ch <- na.omit(ch) #remove NAs, which arise for historical period
  # ch <- ch[order(ch$R2),] #order by R2
  ch <- ch[order(ch$year),] #order by year
  
  fit_cpm <- processStream(ch$excluded,cpmType = "Kolmogorov-Smirnov")
  id <- fit_cpm$changePoints[1] #get index first changepoint
  near <- ch$year[id]
  
  # near <- cor(ch$R2,ch$percent_exlcluded)
  # first_year <- ch$year[ch$nearby_sum<=near]

  output[counter,1] <- "hms"
  output[counter,2] <- e
  output[counter,3] <- m
  output[counter,4] <- near
  counter=counter+1
  }
}

boxplot(output$percent_nearby~output$esm)




ch <- hms_comp[hms_comp$esm=="HAD",]
ch <- ch[!duplicated(ch$year),]
ch <- na.omit(ch)
ch <- ch[order(ch$R2),]
fit_cpm <- processStream(ch$nearby_sum,cpmType = "Kolmogorov-Smirnov")
fit_cpm$changePoints
# ch[c(29,53,72),]
ch[c(12,31),]
plot(ch$R2,ch$nearby_sum)
abline(v=0.729)
abline(v=0.772)

plot(ch$year,ch$nearby_mean, type='l')
plot(ch$R2,ch$nearby_mean)
abline(h=0.729)
abline(h=0.767)

plot(ch$year,ch$excluded,type='l')
abline(v=2032)
abline(v=2047)
abline(v=2059)
abline(v=2085)

fit_cpm <- processStream(ch$R2,cpmType = "Student")
fit_cpm$changePoints

model = list(y~1,1~1,1~1)
fit_mcp = mcp(model,data=ch, par_x="year")

#------FIGURE 22.5: hypervolume and COG anomaly correlations-----
#extrapolation determined from hypervolume code. Computed in another script and output loaded here
#Not a super good story here. 

#Get LATITUDE COGS
species_cogs <- list()
counter_large=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- as.factor(rep(c("HAD","GFDL","IPSL"),each=58000))
  
  #aggregate over space
  all_esm_agg_hist <- all_esm_allyears %>% filter(year<=2010) %>% group_by(esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                                                      "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                                                      "brt_E" , "brt_ES", "brt_EST" ,  
                                                                                                      "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=lat)),na.rm=T)
  all_esm_agg_hist$year <- 2010
  all_esm_agg_hist_longer <- all_esm_agg_hist %>% pivot_longer(c(3:14),names_to="EM", values_to="COG_hist")
  
  all_esm_agg_fcast <- all_esm_allyears %>% filter(year>=2011) %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=lat)),na.rm=T)
  all_esm_agg_fcast_longer <- all_esm_agg_fcast %>% pivot_longer(c(4:15),names_to="EM", values_to="COG_fcast")
  all_esm_agg_fcast_longer_all <- left_join(all_esm_agg_fcast_longer,all_esm_agg_hist_longer, by=c("esm","EM"))
  
  all_esm_agg_fcast_longer_all$COG_anom <- all_esm_agg_fcast_longer_all$COG_fcast - all_esm_agg_fcast_longer_all$COG_hist 
  all_esm_agg_fcast_longer_all_agg <- all_esm_agg_fcast_longer_all %>% group_by(year.x,esm) %>% summarise_at("COG_anom",mean)
  
  all_esm_agg_fcast_longer_all_agg$species <- s
  species_cogs[[counter_large]]<- all_esm_agg_fcast_longer_all_agg
  counter_large=counter_large+1
}


#Get LONGITUDE COGS
species_cogs_lon <- list()
counter_large=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_allyears,gfdl_allyears,ipsl_allyears)
  all_esm_allyears$esm <- as.factor(rep(c("HAD","GFDL","IPSL"),each=58000))
  
  #aggregate over space
  all_esm_agg_hist <- all_esm_allyears %>% filter(year<=2010) %>% group_by(esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                                                      "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                                                      "brt_E" , "brt_ES", "brt_EST" ,  
                                                                                                      "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=lon)),na.rm=T)
  all_esm_agg_hist$year <- 2010
  all_esm_agg_hist_longer <- all_esm_agg_hist %>% pivot_longer(c(3:14),names_to="EM", values_to="COG_hist")
  
  all_esm_agg_fcast <- all_esm_allyears %>% filter(year>=2011) %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                                                            "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                                                            "brt_E" , "brt_ES", "brt_EST" ,  
                                                                                                            "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=lon)),na.rm=T)
  all_esm_agg_fcast_longer <- all_esm_agg_fcast %>% pivot_longer(c(4:15),names_to="EM", values_to="COG_fcast")
  
  all_esm_agg_fcast_longer_all <- left_join(all_esm_agg_fcast_longer,all_esm_agg_hist_longer, by=c("esm","EM"))
  
  all_esm_agg_fcast_longer_all$COG_anom <- all_esm_agg_fcast_longer_all$COG_fcast - all_esm_agg_fcast_longer_all$COG_hist 
  all_esm_agg_fcast_longer_all_agg <- all_esm_agg_fcast_longer_all %>% group_by(year.x,esm) %>% summarise_at("COG_anom",mean)
  
  all_esm_agg_fcast_longer_all_agg$species <- s
  species_cogs_lon[[counter_large]]<- all_esm_agg_fcast_longer_all_agg
  counter_large=counter_large+1
}

#GET distance to coast COGS
species_cog_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
  
  dist <- readRDS('~/Dropbox/PROJECTS/WRAP Location/r_distcoast.rds')
  # dist <- raster::mask(dist,rt)
  dist_df <- as.data.frame(rasterToPoints(dist))
  colnames(dist_df) <- c("lon","lat","distX")
  dist_df$lat <- round(dist_df$lat,2) #CRITICAL STEP
  dist_df$lon <- round(dist_df$lon,2) #CRITICAL STEP
  
  all_esm_allyears_dist <- merge(all_esm_allyears, dist_df, by=c("lon","lat"), all.x=T, all.y=F)
  
  #aggregate over space
  all_esm_agg_hist <-all_esm_allyears_dist %>% filter(year<=2010) %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                                   "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                                   "brt_E" , "brt_ES", "brt_EST" ,  
                                                                                   "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=distX)),na.rm=T)
  all_esm_agg_hist$year <- 2010
  all_esm_agg_hist_longer <- all_esm_agg_hist %>% pivot_longer(c(3:14),names_to="EM", values_to="COG_hist")
  
  all_esm_agg_fcast <-all_esm_allyears_dist %>% filter(year>=2011) %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                                                               "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                                                               "brt_E" , "brt_ES", "brt_EST" ,  
                                                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(weighted.mean(., x=distX)),na.rm=T)
  all_esm_agg_fcast_longer <- all_esm_agg_fcast %>% pivot_longer(c(4:15),names_to="EM", values_to="COG_fcast")
  
  all_esm_agg_fcast_longer_all <- left_join(all_esm_agg_fcast_longer,all_esm_agg_hist_longer, by=c("esm","EM"))
  
  all_esm_agg_fcast_longer_all$COG_anom <- all_esm_agg_fcast_longer_all$COG_fcast - all_esm_agg_fcast_longer_all$COG_hist 
  all_esm_agg_fcast_longer_all_agg <- all_esm_agg_fcast_longer_all %>% group_by(year.x,esm) %>% summarise_at("COG_anom",mean, na.rm=T)
  
  all_esm_agg_fcast_longer_all_agg$species <- s
  species_cog_diff[[counter]]<- all_esm_agg_fcast_longer_all_agg
  counter=counter+1
}


#Novelty File
extrap <- readRDS('~/Dropbox/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet_PlusNearby.rds')

#ggplot for paper
#prepare data LAT cogs:
#hms
hms1 <- extrap[[1]]; hms2 <- extrap[[2]]; hms3 <- extrap[[3]]
hms1$esm = "HAD"; hms2$esm = "GFDL"; hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms_cor <- rbind(species_cogs[[1]])
colnames(hms_cor)[1] <- 'year'
hms_comp <- left_join(hms_cor,hms, by =c("year","esm"))
#cps
cps1 <- extrap[[4]];cps2 <- extrap[[5]];cps3 <- extrap[[6]]
cps1$esm = "HAD";cps2$esm = "GFDL";cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps_cor <- rbind(species_cogs[[2]])
colnames(cps_cor)[1] <- 'year'
cps_comp <- left_join(cps_cor,cps, by =c("year","esm"))
#gfs
gfs1 <- extrap[[7]];gfs2 <- extrap[[8]];gfs3 <- extrap[[9]]
gfs1$esm = "HAD";gfs2$esm = "GFDL";gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs_cor <- rbind(species_cogs[[3]])
colnames(gfs_cor)[1] <- 'year'
gfs_comp <- left_join(gfs_cor,gfs, by =c("year","esm"))

all_sp_latcog <- rbind(hms_comp,cps_comp,gfs_comp)
colnames(all_sp_latcog)[3] <- c("lat_cog_anom")

#prepare data distance to Coast
#hms
hms1 <- extrap[[1]]; hms2 <- extrap[[2]]; hms3 <- extrap[[3]]
hms1$esm = "HAD"; hms2$esm = "GFDL"; hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms_cor <- rbind(species_cog_diff[[1]])
colnames(hms_cor)[1] <- 'year'
hms_comp <- left_join(hms_cor,hms, by =c("year","esm"))
#cps
cps1 <- extrap[[4]];cps2 <- extrap[[5]];cps3 <- extrap[[6]]
cps1$esm = "HAD";cps2$esm = "GFDL";cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps_cor <- rbind(species_cog_diff[[2]])
colnames(cps_cor)[1] <- 'year'
cps_comp <- left_join(cps_cor,cps, by =c("year","esm"))
#gfs
gfs1 <- extrap[[7]];gfs2 <- extrap[[8]];gfs3 <- extrap[[9]]
gfs1$esm = "HAD";gfs2$esm = "GFDL";gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs_cor <- rbind(species_cog_diff[[3]])
colnames(gfs_cor)[1] <- 'year'
gfs_comp <- left_join(gfs_cor,gfs, by =c("year","esm"))

all_sp_distcog <- rbind(hms_comp,cps_comp,gfs_comp)
colnames(all_sp_distcog)[3] <- c("dist_cog_anom")
all_sp_distcog$dist_cog_anom <- all_sp_distcog$dist_cog_anom*-1

all_sp_all_cog <- left_join(all_sp_latcog,all_sp_distcog,by=c("year","esm","species"))
all_sp_all_cog_long <- all_sp_all_cog %>% pivot_longer(c(3,9),values_to = "value",names_to = "cog_type")


all_sp_all_cog_long$cog_type <- as.factor(all_sp_all_cog_long$cog_type)
levels(all_sp_all_cog_long$cog_type) <- c("Distance-to-Coast Centroid (km)", "Latitude Centroid (°)")
all_sp_all_cog_long$cog_type <- factor(all_sp_all_cog_long$cog_type, levels = c("Latitude Centroid (°)", "Distance-to-Coast Centroid (km)"))
all_sp_all_cog_long$species <- as.factor(all_sp_all_cog_long$species)
levels(all_sp_all_cog_long$species) <- c("HMS Archetype", "CPS Archetype", "GFS Archetype")




tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig22.5.tiff',height=7,width=10, units="in", res=300)
ggplot(data=all_sp_all_cog_long,aes(x=value,y=extrap_percent.x,group=esm))+
  geom_point(aes(col=esm))+
  # geom_smooth(aes(col=esm))+
  theme_bw()+
  facet_wrap(~cog_type +species,scales="free")+
  # facet_grid(cog_type ~species,scales="free", switch="y")+
  ylab('Habitat extrapolation from historical conditions (%)') + xlab('Habitat centroid shift from historical mean')+
  theme(legend.position = "bottom", legend.title = element_blank())+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))
dev.off()

#Investagatory plots
#Prepare files for merge
#hms
hms1 <- extrap[[1]]
hms2 <- extrap[[2]]
hms3 <- extrap[[3]]
hms1$esm = "HAD"
hms2$esm = "GFDL"
hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms_cor <- rbind(species_cogs[[1]])
colnames(hms_cor)[1] <- 'year'
hms_comp <- left_join(hms_cor,hms, by =c("year","esm"))
plot(hms_comp$COG_anom,hms_comp$extrap_percent,col=as.factor(hms_comp$esm), pch=19)
cor(hms_comp$COG_anom,hms_comp$extrap_percent)

#cps
cps1 <- extrap[[4]]
cps2 <- extrap[[5]]
cps3 <- extrap[[6]]
cps1$esm = "HAD"
cps2$esm = "GFDL"
cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps_cor <- rbind(species_cogs[[2]])
colnames(cps_cor)[1] <- 'year'
cps_comp <- left_join(cps_cor,cps, by =c("year","esm"))
plot(cps_comp$COG_anom,cps_comp$extrap_percent, col=as.factor(cps_comp$esm), pch=19)
cor(cps_comp$COG_anom,cps_comp$extrap_percent)

#gfs
gfs1 <- extrap[[7]]
gfs2 <- extrap[[8]]
gfs3 <- extrap[[9]]
gfs1$esm = "HAD"
gfs2$esm = "GFDL"
gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs_cor <- rbind(species_cogs[[3]])
colnames(gfs_cor)[1] <- 'year'
gfs_comp <- left_join(gfs_cor,gfs, by =c("year","esm"))
plot(gfs_comp$COG_anom,gfs_comp$extrap_percent,col=as.factor(gfs_comp$esm), pch=19)
cor(gfs_comp$COG_anom,gfs_comp$extrap_percent)

#LONGITUDE
#hms
hms1 <- extrap[[1]]
hms2 <- extrap[[2]]
hms3 <- extrap[[3]]
hms1$esm = "HAD"
hms2$esm = "GFDL"
hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms_cor <- rbind(species_cogs_lon[[1]])
colnames(hms_cor)[1] <- 'year'
hms_comp <- left_join(hms_cor,hms, by =c("year","esm"))
plot(hms_comp$COG_anom,hms_comp$extrap_percent)
cor(hms_comp$COG_anom,hms_comp$extrap_percent)

#cps
cps1 <- extrap[[4]]
cps2 <- extrap[[5]]
cps3 <- extrap[[6]]
cps1$esm = "HAD"
cps2$esm = "GFDL"
cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps_cor <- rbind(species_cogs_lon[[2]])
colnames(cps_cor)[1] <- 'year'
cps_comp <- left_join(cps_cor,cps, by =c("year","esm"))
plot(cps_comp$COG_anom,cps_comp$extrap_percent, col=as.factor(hms_comp$esm), pch=19)
cor(cps_comp$COG_anom,cps_comp$extrap_percent)

#gfs
gfs1 <- extrap[[7]]
gfs2 <- extrap[[8]]
gfs3 <- extrap[[9]]
gfs1$esm = "HAD"
gfs2$esm = "GFDL"
gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs_cor <- rbind(species_cogs_lon[[3]])
colnames(gfs_cor)[1] <- 'year'
gfs_comp <- left_join(gfs_cor,gfs, by =c("year","esm"))
plot(gfs_comp$COG_anom,gfs_comp$extrap_percent)
cor(gfs_comp$COG_anom,gfs_comp$extrap_percent)


#Distance to coast:
#Prepare files for merge
#hms
hms1 <- extrap[[1]]
hms2 <- extrap[[2]]
hms3 <- extrap[[3]]
hms1$esm = "HAD"
hms2$esm = "GFDL"
hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)
hms_cor <- rbind(species_cog_diff[[1]])
colnames(hms_cor)[1] <- 'year'
hms_comp <- left_join(hms_cor,hms, by =c("year","esm"))
plot(hms_comp$COG_anom,hms_comp$extrap_percent,col=as.factor(hms_comp$esm), pch=19)
cor(hms_comp$COG_anom,hms_comp$extrap_percent)

#cps
cps1 <- extrap[[4]]
cps2 <- extrap[[5]]
cps3 <- extrap[[6]]
cps1$esm = "HAD"
cps2$esm = "GFDL"
cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)
cps_cor <- rbind(species_cog_diff[[2]])
colnames(cps_cor)[1] <- 'year'
cps_comp <- left_join(cps_cor,cps, by =c("year","esm"))
plot(cps_comp$COG_anom,cps_comp$extrap_percent, col=as.factor(cps_comp$esm), pch=19)
cor(cps_comp$COG_anom,cps_comp$extrap_percent)

#gfs
gfs1 <- extrap[[7]]
gfs2 <- extrap[[8]]
gfs3 <- extrap[[9]]
gfs1$esm = "HAD"
gfs2$esm = "GFDL"
gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)
gfs_cor <- rbind(species_cog_diff[[3]])
colnames(gfs_cor)[1] <- 'year'
gfs_comp <- left_join(gfs_cor,gfs, by =c("year","esm"))
plot(gfs_comp$COG_anom,gfs_comp$extrap_percent,col=as.factor(gfs_comp$esm), pch=19)
cor(gfs_comp$COG_anom,gfs_comp$extrap_percent)




#Plots
ggplot(data=hms_comp, aes(x=nearby_mean,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=cps_comp, aes(x=nearby_mean,y=R2, col=esm))+
  geom_point(aes(col=esm))
ggplot(data=gfs_comp, aes(x=nearby_mean,y=R2, col=esm))+
  geom_point(aes(col=esm))


#How does each EM perform with novelty?
output <- as.data.frame(matrix(NA, nrow=39, ncol=5))
colnames(output) <- c("species","esm","EM","cor.n", "sig")
counter=1
for (s in c("hms","cps",'gfs')){
  for (e in c("HAD","GFDL","IPSL")){
    for (m in c("gam_E", "gam_ES", "gam_ECor",
                "glm_E" ,"glm_ESt", "glm_ESr",   
                "brt_E" , "brt_ES", "brt_EST" ,  
                "mlp_E" , "mlp_ES" ,"mlp_EST","ens_mean")){
      
      if (s =="hms"){dat = hms_comp}
      if (s =="cps"){dat = cps_comp}
      if (s =="gfs"){dat = gfs_comp}
      
      ch <- dat[dat$esm==e & dat$EM==m,] #get relevant data
      
      cor.n <- cor(ch$R2,ch$nearby_mean, use="na.or.complete")
      sig <- cor.test(ch$R2,ch$nearby_mean, use="na.or.complete")
      
      # near <- cor(ch$R2,ch$percent_exlcluded)
      # first_year <- ch$year[ch$nearby_sum<=near]
      
      output[counter,1] <- s
      output[counter,2] <- e
      output[counter,3] <- m
      output[counter,4] <- cor.n
      output[counter,5] <- sig$p.value
      counter=counter+1
    }
  }
}

boxplot(cor.n~EM ,data=output)
abline(h=0.21)

tiff('~/PROJECTS/WRAP Location/Manuscript/Figures/Fig22.tiff', units="in",res=300,width=12,height=8)
ggplot(data=output, aes(x=EM, y=cor.n)) +
  geom_boxplot(fill="grey") +
  # scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_point(aes(col=esm, shape=species), size=3, alpha=0.9) +
  theme_bw() +
  xlab("Estimation models")+
  ylab("Correlation between model fit and degree of extrapolation")+
  geom_hline(yintercept=0.21)+
  theme(legend.position="bottom")+
  scale_colour_manual(values = c("#d1495b","#00798c","#edae49"))+
  labs(col="Earth System Model", shape="Species")
dev.off()


#is there a global threshold we can calculate?
#Global wont' work becasue species $ nearby changes so much
# and then teste against experiments?
#ON HOLD FOR THE MOMENT --> NEED TO DETERMIEN A METRIC THAT SAYS WHEN A MODEL FAILS


hms_comp$species <- "hms"
cps_comp$species <- "cps"
gfs_comp$species <- "gfs"
global <- rbind(hms_comp, cps_comp, gfs_comp)

m1 <- gam(R2 ~ s(nearby_mean) + species + esm + EM + year, data = global)
summary(m1)
plot(m1,rug=T)
gam.check(m1)

#Changepoint Analysis for %nearby
output <- as.data.frame(matrix(NA, nrow=100, ncol=5))
colnames(output) <- c("species","esm","EM","percent_nearby")
counter=1
for (e in c("HAD","GFDL","IPSL")){
  for (m in c("gam_E", "gam_ES", "gam_ECor",
              "glm_E" ,"glm_ESt", "glm_ESr",   
              "brt_E" , "brt_ES", "brt_EST" ,  
              "mlp_E" , "mlp_ES" ,"mlp_EST","ens_mean")){
    
    ch <- cps_comp[cps_comp$esm==e & cps_comp$EM==m,] #get relevant data
    ch <- na.omit(ch) #remove NAs, which arise for historical period
    # ch <- ch[order(ch$R2),] #order by R2
    ch <- ch[order(ch$year),] #order by year
    
    fit_cpm <- processStream(ch$excluded,cpmType = "Kolmogorov-Smirnov")
    id <- fit_cpm$changePoints[1] #get index first changepoint
    near <- ch$year[id]
    
    # near <- cor(ch$R2,ch$percent_exlcluded)
    # first_year <- ch$year[ch$nearby_sum<=near]
    
    output[counter,1] <- "hms"
    output[counter,2] <- e
    output[counter,3] <- m
    output[counter,4] <- near
    counter=counter+1
  }
}

boxplot(output$percent_nearby~output$esm)




ch <- hms_comp[hms_comp$esm=="HAD",]
ch <- ch[!duplicated(ch$year),]
ch <- na.omit(ch)
ch <- ch[order(ch$R2),]
fit_cpm <- processStream(ch$nearby_sum,cpmType = "Kolmogorov-Smirnov")
fit_cpm$changePoints
# ch[c(29,53,72),]
ch[c(12,31),]
plot(ch$R2,ch$nearby_sum)
abline(v=0.729)
abline(v=0.772)

plot(ch$year,ch$nearby_mean, type='l')
plot(ch$R2,ch$nearby_mean)
abline(h=0.729)
abline(h=0.767)

plot(ch$year,ch$excluded,type='l')
abline(v=2032)
abline(v=2047)
abline(v=2059)
abline(v=2085)

fit_cpm <- processStream(ch$R2,cpmType = "Student")
fit_cpm$changePoints

model = list(y~1,1~1,1~1)
fit_mcp = mcp(model,data=ch, par_x="year")

#------FIGURE 22.6: When does a model fail?----
#Ok so I've been strugglng with how to compare model performance (R2) when I don't really have a threshold/definition
#of what conditions a model FAILS. E.g. Is it predictions 2 standard deviations outside of observed? for each grid cell? 
#For annual biomass? for #COG? 
#Only look at AUC from binomial component? (this might give me some more context to talk about SDMs more broadly)


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



#----Fig 25: which EMs predict biomass best?----
#Get annual estimates of FORECAST RMSE for each EM, species, and ESM.....then make boxplot

#Load in predictions made in EM files
species_predperf <- list()
counter=1
for(s in c("Albacore EMs TempOnly","Anchovy EMs TempOnly","Groundfish EMs TempOnly")){
  print(s)
    had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
    gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
    ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_dat_fcast,gfdl_dat_fcast,ipsl_dat_fcast)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=4500)
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  #Calculate metrics 
  RMSE = function(p, o){(sqrt(mean((p - o)^2)))}
  #RMSE: aggregate over space and time
  all_esm_agg_rmse <-all_esm_allyears %>% group_by(esm) %>% summarise_at(vars("gam_E", "gam_ES", "gam_ECor",
                                                                              "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                              "brt_E" , "brt_ES", "brt_EST" ,  
                                                                              "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(RMSE(., o=abundance)))
  all_esm_agg_rmse_longer <- melt(all_esm_agg_rmse,id=c("esm"), variable.name = "EM",value.name = "RMSE")
  #Make one data frame
  species_predperf[[counter]]<- all_esm_agg_rmse_longer
  counter=counter+1
}

species_predperf[[1]]$species <- "hms"
species_predperf[[2]]$species <- "cps"
species_predperf[[3]]$species <- "gfs"
all_sp <- do.call(rbind, species_predperf)

boxplot(RMSE~EM,data=all_sp)
par(mfrow=c(3,1))
boxplot(RMSE~EM,data=species_predperf[[1]])
boxplot(RMSE~EM,data=species_predperf[[2]])
boxplot(RMSE~EM,data=species_predperf[[3]])

#-----Fig 26: Timeseries of RMSE for temp only-----
#Compare average RMSE for each year to average in historical period, showing ensemble mean
#Load in predictions made in EM files
species_rmse_diff <- list()
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
  # RMSE = function(p, o){(sqrt(mean((p - o)^2))) / mean(o)}
  COR = function(p, o){cor(p,o,use="na.or.complete")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(COR(., o=abundance)))
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

#----Fig 27: compass plots------
#Bearing and magnitude of centre of gravity for each ESM & EM & Species
#Calculate x and y COG, then get distance with pointDistance (raster), and bearing with bearingRhumb (geosphere)
#timescale: native point shoudl be mean COG from 1985-2010, compared to mean from 2075-2100

#Load in predictions made in EM files
compass <- list()
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
  all_esm_agg_latCOG <- all_esm_allyears %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                               "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                               "brt_E" , "brt_ES", "brt_EST" ,  
                                                                               "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(weighted.mean(., x=lat)),na.rm=T)
  all_esm_agg_lonCOG <- all_esm_allyears %>% group_by(year,esm) %>% summarise_at(vars("abundance","gam_E", "gam_ES", "gam_ECor",
                                                                                      "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                                      "brt_E" , "brt_ES", "brt_EST" ,  
                                                                                      "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(weighted.mean(., x=lon)),na.rm=T)
  #pivot longer
  all_esm_agg_latCOG_longer <- melt(all_esm_agg_latCOG,id=c("year", "esm"), variable.name = "EM",value.name = "COG_lat")
  all_esm_agg_lonCOG_longer <- melt(all_esm_agg_lonCOG,id=c("year", "esm"), variable.name = "EM",value.name = "COG_lon")
  
  #get historical means
  all_esm_agg_latCOG_longer_histmean <- all_esm_agg_latCOG_longer[all_esm_agg_latCOG_longer$year<=2010,] %>% group_by(esm,EM) %>% summarise_at(c("COG_lat"), mean)
  all_esm_agg_lonCOG_longer_histmean <- all_esm_agg_lonCOG_longer[all_esm_agg_lonCOG_longer$year<=2010,] %>% group_by(esm,EM) %>% summarise_at(c("COG_lon"), mean)
  hist_cog <- left_join(all_esm_agg_lonCOG_longer_histmean,all_esm_agg_latCOG_longer_histmean, by =c("esm", "EM"))
  colnames(hist_cog) <- c("esm","EM","COG_lon_hist","COG_lat_hist")
  
  #get future means  
  all_esm_agg_latCOG_longer_fcastmean <- all_esm_agg_latCOG_longer[all_esm_agg_latCOG_longer$year>=2075,] %>% group_by(esm,EM) %>% summarise_at(c("COG_lat"), mean)
  all_esm_agg_lonCOG_longer_fcastmean <- all_esm_agg_lonCOG_longer[all_esm_agg_lonCOG_longer$year>=2075,] %>% group_by(esm,EM) %>% summarise_at(c("COG_lon"), mean)
  fcast_cog <- left_join(all_esm_agg_lonCOG_longer_fcastmean,all_esm_agg_latCOG_longer_fcastmean, by =c("esm", "EM"))
  colnames(fcast_cog) <- c("esm","EM","COG_lon_fcast","COG_lat_fcast")
  
  #Combine into one df, and get bearing and distance
  compass_df <- left_join(hist_cog,fcast_cog,by =c("esm", "EM"))
  dist <- compass_df %>% group_by(esm,EM) %>% summarise(dist = pointDistance(c(COG_lon_hist, COG_lat_hist),c(COG_lon_fcast, COG_lat_fcast),lonlat=T, allpairs = T))
  bearing <- compass_df %>% group_by(esm,EM) %>% summarise(bearing = bearingRhumb(c(COG_lon_hist, COG_lat_hist),c(COG_lon_fcast, COG_lat_fcast)))
  
  #Combine again into a final df
  fin_df <- left_join(dist,bearing,by=c("esm","EM"))
  fin_df$dist <- fin_df$dist/1000
  # fin2 <- left_join(fin_df,rmse_dist,by=c("esm","EM"))
  
  compass[[counter]]<- fin_df
  counter=counter+1
  
  
  #IF I wanted to include RMSE then use this:
  #get historical means
  # all_esm_agg_latCOG_longer_histmean <- all_esm_agg_latCOG_longer[all_esm_agg_latCOG_longer$year<=2010,] %>% group_by(esm,EM) %>% summarise_at(c("COG_lat", "abundance"), mean)
  # colnames(all_esm_agg_latCOG_longer_histmean)[4] <- "abundance_lat"
  # all_esm_agg_lonCOG_longer_histmean <- all_esm_agg_lonCOG_longer[all_esm_agg_lonCOG_longer$year<=2010,] %>% group_by(esm,EM) %>% summarise_at(c("COG_lon", "abundance"), mean)
  # colnames(all_esm_agg_lonCOG_longer_histmean)[4] <- "abundance_lon"
  # hist_cog <- left_join(all_esm_agg_lonCOG_longer_histmean,all_esm_agg_latCOG_longer_histmean, by =c("esm", "EM"))
  # colnames(hist_cog) <- c("esm","EM","COG_lon_hist","abundance_lon_hist","COG_lat_hist", "abundance_lat_hist")
  # 
  # #get future means  
  # all_esm_agg_latCOG_longer_fcastmean <- all_esm_agg_latCOG_longer[all_esm_agg_latCOG_longer$year>=2075,] %>% group_by(esm,EM) %>% summarise_at(c("COG_lat","abundance"), mean)
  # colnames(all_esm_agg_latCOG_longer_fcastmean)[4] <- "abundance_lat"
  # all_esm_agg_lonCOG_longer_fcastmean <- all_esm_agg_lonCOG_longer[all_esm_agg_lonCOG_longer$year>=2075,] %>% group_by(esm,EM) %>% summarise_at(c("COG_lon","abundance"), mean)
  # colnames(all_esm_agg_lonCOG_longer_fcastmean)[4] <- "abundance_lon"
  # fcast_cog <- left_join(all_esm_agg_lonCOG_longer_fcastmean,all_esm_agg_latCOG_longer_fcastmean, by =c("esm", "EM"))
  # colnames(fcast_cog) <- c("esm","EM","COG_lon_fcast","abundance_lon_fcast","COG_lat_fcast","abundance_lat_fcast")
  # 
  # #Combine into one df, and get bearing and distance
  # compass_df <- left_join(hist_cog,fcast_cog,by =c("esm", "EM"))
  # dist <- compass_df %>% group_by(esm,EM) %>% summarise(dist = pointDistance(c(COG_lon_hist, COG_lat_hist),c(COG_lon_fcast, COG_lat_fcast),lonlat=T, allpairs = T),
  #                                                       true_dist = pointDistance(c(abundance_lon_hist, abundance_lat_hist),c(abundance_lon_fcast, abundance_lat_fcast),lonlat=T, allpairs = T))
  # bearing <- compass_df %>% group_by(esm,EM) %>% summarise(bearing = bearingRhumb(c(COG_lon_hist, COG_lat_hist),c(COG_lon_fcast, COG_lat_fcast)),
  #                                                          true_bearing = bearingRhumb(c(abundance_lon_hist, abundance_lat_hist),c(abundance_lon_fcast, abundance_lat_fcast)))
  # 
  # rmse_dist <- dist %>% group_by(esm,EM) %>%  summarise(rmse_dist = RMSE(p= dist, o = true_dist))
  # rmse_bearing <- bearing %>% group_by(esm,EM) %>%  summarise(rmse_bearing = RMSE(p= bearing, o = true_bearing))
  
  
  
}


hms <- rbind(compass[[1]],compass[[4]])
hms$type <- rep(c("Full","TempOnly"),each=42)

cps <- rbind(compass[[2]],compass[[5]])
cps$type <- rep(c("Full","TempOnly"),each=42)

gfs <- rbind(compass[[3]],compass[[6]])
gfs$type <- rep(c("Full","TempOnly"),each=42)

#PLOTS
ggplot()+
  geom_bar(data=hms[hms$EM=="abundance" & hms$type=="Full",],aes(x=bearing,y=dist, color=esm, fill=esm), width=1,stat="identity")+
  coord_polar()+
  xlim(0,360)+
  theme_bw()+
  theme(legend.position = "bottom")+
  scale_fill_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))

ggplot()+
  geom_bar(data=cps[cps$EM=="abundance" & cps$type=="Full",],aes(x=bearing,y=dist, color=esm, fill=esm), width=1,stat="identity")+
  coord_polar()+
  xlim(0,360)
ggplot()+
  geom_bar(data=gfs[gfs$EM=="abundance" & gfs$type=="Full",],aes(x=bearing,y=dist, color=esm, fill=esm, alpha=0.5), width=10,stat="identity")+
  coord_polar()+
  xlim(0,360)

#PLOTS
ggplot()+
  geom_bar(data=cps[cps$esm=="HAD",],aes(x=bearing,y=dist, color=type, fill=type),stat="identity")+
  geom_bar(data=cps[cps$esm=="HAD" & cps$EM=="abundance",][1,],aes(x=bearing,y=dist, color="black", width=5),stat="identity")+
  coord_polar()+
  xlim(0,360)

#PLOTS
ggplot()+
  geom_bar(data=gfs[gfs$esm=="HAD",],aes(x=bearing,y=dist, color=type, fill=type),stat="identity")+
  geom_bar(data=gfs[gfs$esm=="HAD" & gfs$EM=="abundance",][1,],aes(x=bearing,y=dist, color="black", width=5),stat="identity")+
  coord_polar()+
  xlim(0,360)



ggplot(fin_df,aes(x=bearing,y=dist,color=EM),guide=F)+
  geom_rect(xmin=0,xmax=Inf,ymin=0,ymax=max(Inf),color=NA,alpha=.1,fill="#DCDCDC")+
  geom_segment(aes(xend = bearing, yend = 0.1),size=1)+
  geom_hline(yintercept = max(fin_df$dist),color="black")+
  geom_hline(yintercept = ybreaks[1],color="grey")+
  geom_hline(yintercept = ybreaks[2],color="grey")+
  geom_vline(xintercept = 90,color="grey")+
  geom_vline(xintercept = 0,color="grey")+
  geom_vline(xintercept = 180,color="grey")+
  geom_vline(xintercept = 270,color="grey")+
  theme_bw()+
  
  scale_color_manual("Marine heatwave event",values=c("2014"="#EFC000FF","2015"="#0073C2FF","2019"="#8a131f","2020"="darkgreen"),
                     labels=c("Blob14"="2014","Blob15"="2015","MHW19"="2019","MHW20"="2020"))+
  scale_shape_manual(values=shapes)+
  # theme_linedraw()+
  # theme_classic()+
  geom_point(size=3) +
  # scale_shape_manual(values=shapes)+
  xlab(NULL)+
  ylab(NULL)+
  geom_text(data = data.frame(x = 0, y = ybreaks, label = ybreaks),
            aes(x = x, y = y, label = label),
            inherit.aes = F,
            size = 3) +
  scale_x_continuous(limits = c(0, 360),breaks =  c(0,90,180,270),labels=c("North","East","South","West")) +
  scale_y_continuous(breaks =  ybreaks,labels=ybreaks,expand = c(0,0)) +
  coord_polar(start=0,clip="off") +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_line(size = 0),
        axis.text.x=element_text(size=10),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        # panel.grid.minor = element_line(
        #   size = 0.25, linetype = 'solid', colour = "black"),
        # panel.grid.major = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle(glue("{sp[i]} species"))#+ 



#----Figure 28: plot hypervolume novelty trends-----
#Novelty File
extrap <- readRDS('~/PROJECTS/WRAP Location/Hypervolume_ExtrapolationOutput_PlusExDet.rds')
#Prepare files for merge
#hms
hms1 <- extrap[[1]]
hms2 <- extrap[[2]]
hms3 <- extrap[[3]]
hms1$esm = "HAD"
hms2$esm = "GFDL"
hms3$esm = "IPSL"
hms <- rbind(hms1,hms2,hms3)

#cps
cps1 <- extrap[[4]]
cps2 <- extrap[[5]]
cps3 <- extrap[[6]]
cps1$esm = "HAD"
cps2$esm = "GFDL"
cps3$esm = "IPSL"
cps <- rbind(cps1,cps2,cps3)

#gfs
gfs1 <- extrap[[7]]
gfs2 <- extrap[[8]]
gfs3 <- extrap[[9]]
gfs1$esm = "HAD"
gfs2$esm = "GFDL"
gfs3$esm = "IPSL"
gfs <- rbind(gfs1,gfs2,gfs3)

hms$species <- "hms"
cps$species <- "cps"
gfs$species <- "gfs"

dat <- rbind(hms,cps,gfs)

ggplot(data=dat,aes(x=year,y=excluded,group=esm))+
  geom_line(aes(col=esm))+
  facet_wrap(~species, scales="free_y")+
  theme_classic()+
  theme(legend.position = "bottom")

#----Figure 29: comparison with temp-only----
#boxplot of projected R2 for full models and temp only

#Load in predictions made in EM files
cors_all <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs", 
           "Albacore EMs TempOnly","Anchovy EMs TempOnly","Groundfish EMs TempOnly")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/Dropbox/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
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
  
  
  all_esm_allyears <- all_esm_allyears[all_esm_allyears$year>=2011,]
  
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ,"ens_mean"),funs(cor(., abundance)))
  all_esm_agg_longer <- all_esm_agg %>%  pivot_longer(c(2:14),names_to="EM", values_to="R2",) 
  all_esm_agg_longer$species <- s
  cors_all[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

cors_all_df <- do.call(rbind,cors_all)
idx<- grep("TempOnly",cors_all_df$species)
cors_all_df$split <- "Full SDMs"
cors_all_df$split[idx] <- "Temperature-Only SDMs"

cors_all_df$species <- ifelse(cors_all_df$species=="Albacore EMs", "HMS",
                              ifelse(cors_all_df$species=="Anchovy EMs", "CPS",
                                     ifelse(cors_all_df$species=="Groundfish EMs", "GFS",
                                            ifelse(cors_all_df$species=="Albacore EMs TempOnly", "HMS",
                                                   ifelse(cors_all_df$species=="Anchovy EMs TempOnly", "CPS",
                                                          ifelse(cors_all_df$species=="Groundfish EMs TempOnly", "GFS", NA))))))

cors_all_df$EM <- as.factor(cors_all_df$EM)
levels(cors_all_df$EM) <-  list("brt_E" ="brt_E" ,   "brt_ES" = "brt_ES",  "brt_EST" = "brt_EST",
                            "gam_E"= "gam_E","gam_ECor"="gam_ECor", "gam_ES"= "gam_ES",
                           "glmm_ES"= "glm_E" ,  "glmm_EST" = "glm_ESr",  "glmm_ESTt"="glm_ESt", 
                           "mlp_E"="mlp_E",   "mlp_ES" ="mlp_ES" ,  "mlp_EST"="mlp_EST" )
cors_all_df <- na.omit(cors_all_df)

g1 <- ggplot(data=cors_all_df, aes(x=EM, y=R2))+
  geom_boxplot()+
  geom_point(aes(col=esm, shape=species), size=3, alpha=0.9) +
  theme_bw() +
  scale_color_manual(values=c("#2980b9", "#2c3e50","#c0392b"))+
  labs(col="Earth System Model", shape="Species")+
  xlab("Estimation models")+
  ylab("Correlation between observed and predicted biomass")+
  theme(legend.position = "bottom")+
  facet_wrap(~split, dir="v", scales="free_y")+
  ylim(0,1)

tiff('~/Dropbox/PROJECTS/WRAP Location/Manuscript/Figures/Fig29.tiff', res=300, units="in", width=10, height=8)
plot(g1)
dev.off()


#----Extrapolation Test----

x <- rnorm(100,mean = 1, sd = 10)
y <- rnorm(100,mean = 3, sd = 100)
df <- as.data.frame(cbind(x,y))

g1 <- gam(y~x)
b1 <- gbm.step(data = df, gbm.x="x", gbm.y="y", family= "gaussian")
m1 <- neuralnet(y~x, data = df)

df$g1 <- predict(g1, df)
df$b1 <- predict(b1, df)
df$m1 <- predict(m1, df)

plot(df$x,df$y)
points(df$x,df$g1)
points(df$x,df$b1)
points(df$x,df$m1)

ex <- df
ex$x <- ex$x+100
ex$g1 <- predict(g1, ex)
ex$b1 <- predict(b1, ex)
ex$m1 <- predict(m1, ex)
plot(ex$x,ex$y)
points(ex$x,ex$g1)
points(ex$x,ex$b1)
points(ex$x,ex$m1)


l1 <- sdmTMB(formula = y~x, family = gaussian(link = "identity"),data = df)
time_varying = NULL,
spde = spde_pos,
time = "year",
family = gaussian(link = "identity"),
data = dat_pres,
anisotropy = TRUE,
ar1_fields = FALSE,
weights = weights[which(dat_pres$pres == 1)],
spatial_only = TRUE,
spatial_trend = FALSE,
quadratic_roots = TRUE,
control = sdmTMBcontrol(step.min = 0.01, step.max = 1)


install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
pcod_spde <- make_mesh(pcod_2011, c("X", "Y"), cutoff = 25) # a coarse mesh for example speed
plot(pcod_spde)
pcod_gaus <- subset(pcod_2011, density > 0 & year >= 2013)
pcod_spde_gaus <- make_mesh(pcod_gaus, c("X", "Y"), cutoff = 30)
m_pos <- sdmTMB(log(density) ~ 0 + as.factor(year) + depth_scaled + depth_scaled2,
                data = pcod_gaus, time = "year", spde = pcod_spde_gaus)
print(m_pos)

preds <- predict(m_pos, pcod_gaus)  
plot(pcod_gaus$density,pcod_gaus$depth_scaled)
points(preds$density,pcod_gaus$depth_scaled)


#----Tables: Mean RMSE for forecast period....which EM is best under conditins?----
#Compare average RMSE for each year to average in historical period, showing ensemble mean

rmse_mean <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs",
           "Albacore EMs TempOnly","Anchovy EMs TempOnly","Groundfish EMs TempOnly",
           "Albacore EMs 2040train","Anchovy EMs 2040train","Groundfish EMs 2040train")){
  print(s)
  # had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  # gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  # ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))

  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_dat_fcast,gfdl_dat_fcast,ipsl_dat_fcast)
  all_esm_allyears$esm <- rep(c("HAD","GFDL","IPSL"),each=nrow(ipsl_dat_fcast))
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
  COR = function(p, o){cor(p,o,use="na.or.complete")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST"),funs(COR(., o=abundance)))
  all_esm_agg_longer <- melt(all_esm_agg,id=c("esm"), variable.name = "EM",value.name = "RMSE")
  all_esm_agg_longer$species <- s
  rmse_mean[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

dat <- do.call(rbind,rmse_mean)

ggplot(dat,aes(y=RMSE, x= EM))+
  geom_bar(stat="identity")+
  facet_wrap(~species+esm)

mins <- dat %>% group_by(esm,species) %>% slice(which.max(RMSE))

mean(dat$RMSE[dat$species=="Groundfish EMs"])
mean(dat$RMSE[dat$species=="Groundfish EMs TempOnly"])

hms <- rmse_mean[[8]]
hms[order(hms$RMSE),]

#Questions:
#When enviro change is large versus small? What is best EM?
#When species population dynamics are hard to predict? what is best EM?
#When you underparameterise a model? what is best EM?
#When have have strong spatial preferences. What is best EM?


#----TABLE: mean R2 for fitted period------

species_rmse_diff <- list()
counter=1
for(s in c("Albacore EMs","Anchovy EMs","Groundfish EMs",
           "Albacore EMs TempOnly","Anchovy EMs TempOnly","Groundfish EMs TempOnly")){
  print(s)
  had_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_hist_results_full.rds"))
  # had_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","had","/dat_fcast_results_full.rds"))
  gfdl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_hist_results_full.rds"))
  # gfdl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","gfdl","/dat_fcast_results_full.rds"))
  ipsl_dat_hist <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_hist_results_full.rds"))
  # ipsl_dat_fcast <- readRDS(paste0("~/PROJECTS/WRAP Location/",s,"/","ipsl","/dat_fcast_results_full.rds"))
  #Rbind historical and forecast years
  # had_allyears <- rbind(had_dat_hist,had_dat_fcast)
  # gfdl_allyears <- rbind(gfdl_dat_hist,gfdl_dat_fcast)
  # ipsl_allyears <- rbind(ipsl_dat_hist,ipsl_dat_fcast)
  #Combine all ESMs into one file
  all_esm_allyears <- rbind(had_dat_hist,gfdl_dat_hist,ipsl_dat_hist)
  all_esm_allyears$esm <- as.factor(rep(c("HAD","GFDL","IPSL"),each=nrow(ipsl_dat_hist)))
  
  #keep necessary columns
  all_esm_allyears <- all_esm_allyears[,c("lon","lat","year","abundance",
                                          "gam_E","gam_ES", "gam_ECor",
                                          "glm_E" ,"glm_ESt", "glm_ESr",   
                                          "brt_E" , "brt_ES", "brt_EST" ,  
                                          "mlp_E" , "mlp_ES" ,"mlp_EST","esm")]
  
  #Calculate RMSE
  COR = function(x,y){cor(x,y,method="spearman")}
  #aggregate over space
  all_esm_agg <-all_esm_allyears %>% group_by(year,esm) %>% summarise_at(c("gam_E", "gam_ES", "gam_ECor",
                                                                           "glm_E" ,"glm_ESt", "glm_ESr",   
                                                                           "brt_E" , "brt_ES", "brt_EST" ,  
                                                                           "mlp_E" , "mlp_ES" ,"mlp_EST" ),funs(COR(., y=abundance)))
  
  all_esm_agg_longer <- melt(all_esm_agg,id=c("year", "esm"), variable.name = "EM",value.name = "COR")
  species_rmse_diff[[counter]]<- all_esm_agg_longer
  counter=counter+1
}

mean(species_rmse_diff[[1]]$COR)^2
mean(species_rmse_diff[[2]]$COR)^2
mean(species_rmse_diff[[3]]$COR)^2

mean(species_rmse_diff[[4]]$COR)^2
mean(species_rmse_diff[[5]]$COR)^2
mean(species_rmse_diff[[6]]$COR)^2



#----TABLE: mean AUC for fitted period------

#Load in each pres-abs model object
#Then predict on 100% fitted data
#Get AUC

eval <- function(x){
  d <- cbind(dat_hist$pres, x)
  pres <- as.numeric(d[d[,1]==1,2])
  abs <- as.numeric(d[d[,1]==0,2])
  a <- evaluate(p=pres, a=abs)
  return(a@auc)
}


auc_list <- list()
counter=1
for (s in c("Albacore EMs TempOnly","Anchovy EMs TempOnly", "Groundfish EMs TempOnly")){
  for (e in c("had","gfdl","ipsl")){
    #Load data and models
    dat_hist <- readRDS(glue("~/PROJECTS/WRAP Location/{s}/{e}/dat_hist_results_full.rds"))
    load(glue("~/PROJECTS/WRAP Location/{s}/{e}/saved_models_full.RData"))  #full models
    
    # for (m in c(gam_E_P, gam_ES_P, gam_ECor_P)){
      # brt_E_P, brt_ES_P, brt_EST_P, 
      # mlp_E_P, mlp_ES_P, mlp_EST_P, 
      # glm_E_P, glm_ESr_P, glm_ESt_P)){
      print(glue('running species {s}'))
      print(e)
      
      #GAMs
      preds <- predict(gam_E_P,dat_hist, type="response")
      gam_E <- eval(x=preds)
      preds <- predict(gam_ES_P,dat_hist, type="response")
      gam_ES <- eval(x=preds)
      preds <- predict(gam_ECor_P$gam,dat_hist, type="response")
      gam_ECor <- eval(x=preds)
      
      #BRTs
      preds <- predict(brt_E_P,dat_hist,n.trees=brt_E_P$gbm.call$best.trees, type="response")
      brt_E <- eval(x=preds)
      preds <- predict(brt_ES_P,dat_hist,n.trees=brt_ES_P$gbm.call$best.trees, type="response")
      brt_ES <- eval(x=preds)
      preds <- predict(brt_EST_P,dat_hist,n.trees=brt_EST_P$gbm.call$best.trees, type="response")
      brt_EST <- eval(x=preds)
      
      #MLPs
      preds <- predict(mlp_E_P,dat_hist, type="response")
      mlp_E <- eval(x=preds)
      preds <- predict(mlp_ES_P,dat_hist, type="response")
      mlp_ES <- eval(x=preds)
      preds <- predict(mlp_EST_P,dat_hist, type="response")
      mlp_EST <- eval(x=preds)
      
      #GLMs
      # preds <- predict(glm_E_N,dat_hist, type="response")
      # glm_E <- eval(x=preds)
      # preds <- predict(glm_ES_P,dat_hist, type="response")
      # glm_ES <- eval(x=preds)
      # preds <- predict(glm_EST_P,dat_hist, type="response")
      # glm_EST <- eval(x=preds)
      
      df <- data.frame(gam_E, gam_ES, gam_ECor,
                       brt_E, brt_ES, brt_EST,
                       mlp_E, mlp_ES, mlp_EST)
      df <- df %>% pivot_longer(1:9,names_to="EM", values_to = "AUC")
      df$species <- s
      df$esm <- e

      auc_list[[counter]] <-df

      counter=counter+1
      
    }
  }
}

auc_vals <- do.call(rbind,auc_list)

auc_vals_agg <- auc_vals %>% group_by(species) %>% summarise(mean(AUC))


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





  
  
  
  
  
  