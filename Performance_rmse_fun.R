### Location^3 Estimation Model Code
## Function to take observed and predicted values from SDM EMs and measure their performance

## Code by James Smith - April 2020


sdm_rmse <- function(dat_hist, dat_fcast) {
  #'dat_hist' are the 'historical' observations from the OM, and predicted abundance values from fitted models as additional columns
  #'dat_fcast' are the 'future' observations from the OM,  and predicted abundance values from fitted models as additional columns

  # Perf. metrics are:
  #       - RMSE and Cor between entire subset of observed and predicted abundances, historical and future
  #       - Plots and RMSE of total annual abundance, historical and future
  
  # There are maximum 5 GAMS, 4 BRTs, 4 MLPs, and 4 GLMs (17 total)

  
  RMSE = function(p, o){
    sqrt(mean((p - o)^2))
  }
  
  all_mods <- c("gam_E","gam_S","gam_ES","gam_EST","gam_ECor",
                "brt_E","brt_S","brt_ES","brt_EST",
                "mlp_E","mlp_S","mlp_ES","mlp_EST",
                "glm_E","glm_Sr","glm_ESt","glm_ESr")  #all possible models
  
  mods <- names(dat_hist)[names(dat_hist) %in% all_mods]   #which of the possible models have been run
  
  rmse_abund_all <- data.frame(model=mods, rmse_hist=0, rmse_fcast=0)
  rmse_abund_annual <- data.frame(model=mods, rmse_hist=0, rmse_fcast=0)
  
  rmse_abund_all_y <- as.data.frame(matrix(0, nrow=length(unique(dat_fcast$year)), ncol=length(mods)+1))
  names(rmse_abund_all_y) <- c("year", mods)
  rmse_abund_all_y$year <- unique(dat_fcast$year)
  
  rmse_abund_allh_y <- as.data.frame(matrix(0, nrow=length(unique(dat_hist$year)), ncol=length(mods)+1))
  names(rmse_abund_allh_y) <- c("year", mods)
  rmse_abund_allh_y$year <- unique(dat_hist$year)

  # Calculate RMSE, and plot total abundance trajectories
  
  #Hist
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {
    rmse_abund_all$rmse_hist[rmse_abund_all$model==m]  <- round(RMSE(dat_hist$abundance, dat_hist[,m]),3)
    rmse_abund_all$rmse_fcast[rmse_abund_all$model==m]  <- round(RMSE(dat_fcast$abundance, dat_fcast[,m]),3)
    
    # hist
    ann_obs <- aggregate(abundance~year, dat=dat_hist, FUN="sum")
    ann_mod <- aggregate(formula(paste(m,"~ year")), dat=dat_hist, FUN="sum")
    rmse_abund_annual$rmse_hist[rmse_abund_all$model==m] <- RMSE(ann_obs$abundance, ann_mod[,m])
    
    plot(ann_obs, type="l", lwd=2, ylim=ylim, xlab="",
         main=paste0("Hist, ", m, ", rmse=",round(RMSE(dat_hist$abundance, dat_hist[,m]),2),
                     ", cor=", round(cor(dat_hist$abundance, dat_hist[,m]),2)))
    lines(ann_mod, lwd=2, col="red")
    
    for (yy in rmse_abund_allh_y$year) {  #rmse of all 400 observations per year
      dat_hist_yy <- dat_hist[dat_hist$year == yy,]
      rmse_abund_allh_y[rmse_abund_allh_y$year == yy,m] <- round(RMSE(dat_hist_yy$abundance, dat_hist_yy[,m]),3)
    }
  }
  
  #Fcast
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {
    ann_obs <- aggregate(abundance~year, dat=dat_fcast, FUN="sum")
    ann_mod <- aggregate(formula(paste(m,"~ year")), dat=dat_fcast, FUN="sum")
    rmse_abund_annual$rmse_fcast[rmse_abund_all$model==m] <- RMSE(ann_obs$abundance, ann_mod[,m])
    
    plot(ann_obs, type="l", lwd=2, ylim=ylim, xlab="",
         main=paste0("Fcast, ", m, ", rmse=",round(RMSE(dat_fcast$abundance, dat_fcast[,m]),2),
                     ", cor=", round(cor(dat_fcast$abundance, dat_fcast[,m]),2)))
    lines(ann_mod, lwd=2, col="red")
    
    for (yy in rmse_abund_all_y$year) {  #rmse of all 400 observations per year
      dat_fcast_yy <- dat_fcast[dat_fcast$year == yy,]
      rmse_abund_all_y[rmse_abund_all_y$year == yy,m] <- round(RMSE(dat_fcast_yy$abundance, dat_fcast_yy[,m]),3)
    }
  }
  
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  
  for (m in mods) {  #plot degregdation of RMSE over time (if any)
    plot(rmse_abund_all_y$year, rmse_abund_all_y[,m], type="l", lwd=2, col="red",
         ylab="rmse", xlab="", main=paste(m), ylim=c(min(rmse_abund_all_y[,m]),max(rmse_abund_all_y[,m])))
    abline(h=mean(rmse_abund_allh_y[,m]), lty=2, col="blue")  #mean annual rmse of historical period
  }

  
  return(list(rmse_abund_all = rmse_abund_all,
              rmse_abund_annual = rmse_abund_annual,
              rmse_abund_all_year = rmse_abund_all_y))
         
}


  
