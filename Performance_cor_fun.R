### Location^3 Estimation Model Code
## Function to get correlation coefficient between observed and predicted values from SDM EMs
## Code by James Smith & Steph Brodie


sdm_cor <- function(dat_hist, dat_fcast) {
  #'dat_hist' are the 'historical' observations from the OM, and predicted abundance values from fitted models as additional columns
  #'dat_fcast' are the 'future' observations from the OM,  and predicted abundance values from fitted models as additional columns

  # Perf. metrics are:
  #       - Annual correlations between observed and predicted abundances, historical and future
  #       - Plots of correlation time-series, historical and future
  
  # There are maximum 5 GAMS, 4 BRTs, 4 MLPs, and 4 GLMs (17 total)

  
  COR = function(p, o){
    cor(p,o,method="spearman", use = "na.or.complete")
  }
  
  all_mods <- c("gam_E","gam_S","gam_ES","gam_EST","gam_ECor",
                "brt_E","brt_S","brt_ES","brt_EST",
                "mlp_E","mlp_S","mlp_ES","mlp_EST",
                "glm_E","glm_Sr","glm_ESt","glm_ESr")  #all possible models
  
  mods <- names(dat_hist)[names(dat_hist) %in% all_mods]   #which of the possible models have been run
  
  cor_abund_all <- data.frame(model=mods, cor_hist=0, cor_fcast=0)
  cor_abund_annual <- data.frame(model=mods, cor_hist=0, cor_fcast=0)
  
  cor_abund_all_y <- as.data.frame(matrix(0, nrow=length(unique(dat_fcast$year)), ncol=length(mods)+1))
  names(cor_abund_all_y) <- c("year", mods)
  cor_abund_all_y$year <- unique(dat_fcast$year)
  
  cor_abund_allh_y <- as.data.frame(matrix(0, nrow=length(unique(dat_hist$year)), ncol=length(mods)+1))
  names(cor_abund_allh_y) <- c("year", mods)
  cor_abund_allh_y$year <- unique(dat_hist$year)

  # Calculate COR, and plot time series
  
  #Hist
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {
    cor_abund_all$cor_hist[cor_abund_all$model==m]  <- round(COR(dat_hist$abundance, dat_hist[,m]),3)
    cor_abund_all$cor_fcast[cor_abund_all$model==m]  <- round(COR(dat_fcast$abundance, dat_fcast[,m]),3)
    
    # hist
    ann_obs <- aggregate(abundance~year, dat=dat_hist, FUN="sum")
    ann_mod <- aggregate(formula(paste(m,"~ year")), dat=dat_hist, FUN="sum")
    cor_abund_annual$cor_hist[cor_abund_all$model==m] <- COR(ann_obs$abundance, ann_mod[,m])
    
    plot(ann_obs, type="l", lwd=2, ylim=ylim, xlab="",
         main=paste0("Hist, ", m))
    lines(ann_mod, lwd=2, col="red")
    
    for (yy in cor_abund_allh_y$year) {  #cor of all 400 observations per year
      dat_hist_yy <- dat_hist[dat_hist$year == yy,]
      cor_abund_allh_y[cor_abund_allh_y$year == yy,m] <- round(cor(dat_hist_yy$abundance, dat_hist_yy[,m]),3)
    }
  }
  
  #Fcast
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {
    ann_obs <- aggregate(abundance~year, dat=dat_fcast, FUN="sum")
    ann_mod <- aggregate(formula(paste(m,"~ year")), dat=dat_fcast, FUN="sum")
    cor_abund_annual$cor_fcast[cor_abund_all$model==m] <- COR(ann_obs$abundance, ann_mod[,m])
    
    plot(ann_obs, type="l", lwd=2, ylim=ylim, xlab="",
         main=paste0("Fcast, ", m))
    lines(ann_mod, lwd=2, col="red")
    
    for (yy in cor_abund_all_y$year) {  #cor of all 400 observations per year
      dat_fcast_yy <- dat_fcast[dat_fcast$year == yy,]
      cor_abund_all_y[cor_abund_all_y$year == yy,m] <- round(COR(dat_fcast_yy$abundance, dat_fcast_yy[,m]),3)
    }
  }
  
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {  #plot degregdation of COR over time (if any)
    plot(cor_abund_all_y$year, cor_abund_all_y[,m], type="l", lwd=2, col="red",
         ylab="cor", xlab="", main=paste(m), ylim=c(min(cor_abund_all_y[,m]),max(cor_abund_all_y[,m])))
    abline(h=mean(cor_abund_allh_y[,m]), lty=2, col="blue")  #mean annual cor of historical period
  }

  
  return(list(cor_abund_all = cor_abund_all,
              cor_abund_annual = cor_abund_annual,
              cor_abund_all_year = cor_abund_all_y))
         
}


  
