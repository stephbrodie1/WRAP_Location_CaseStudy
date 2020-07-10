### Location^3 Estimation Model Code
## Function to measure latitudinal Centre Of Gravity for each SDM EM

## Code by Steph Brodie and James Smith - April 2020


sdm_cog <- function(dat_hist, dat_fcast) {
  #'dat_hist' are the 'historical' observations from the OM, and predicted abundance values from fitted models as additional columns
  #'dat_fcast' are the 'future' observations from the OM,  and predicted abundance values from fitted models as additional columns
  #'mods' are the types of models to report: c("gam", "glm", "brt", "mlp")
  
  # Perf. metrics are:
  #       - Plots of weighted annual mean latitude (COG)
  #       - RMSE of annual COG values
  
  # There are maximum 5 GAMS, 4 BRTs, 4 MLPs, and 4 GLMs (17 total)
  
  
  RMSE = function(p, o){
    sqrt(mean((p - o)^2))
  }
  
  all_mods <- c("gam_E","gam_S","gam_ES","gam_EST","gam_ECor",
                "brt_E","brt_S","brt_ES","brt_EST",
                "mlp_E","mlp_S","mlp_ES","mlp_EST",
                "glm_E","glm_Sr","glm_ESt","glm_ESr")

  mods <- names(dat_hist)[names(dat_hist) %in% all_mods]   #which of the possible models have been run
  
  rmse_cog <- data.frame(model=mods, rmse_hist=0, rmse_fcast=0)
  
  years_hist <- unique(dat_hist$year)
  years_fcast <- unique(dat_fcast$year)
  
  cog_lat_hist <- as.data.frame(matrix(0,nrow=length(years_hist),ncol=2+length(mods)))
  cog_lat_fcast <- as.data.frame(matrix(0,nrow=length(years_fcast),ncol=2+length(mods)))
  names(cog_lat_hist) <- c("year", "truth", mods)
  names(cog_lat_fcast) <- c("year", "truth", mods)
  cog_lat_hist$year <- years_hist
  cog_lat_fcast$year <- years_fcast
  
  cog_lat_fcast_y <- as.data.frame(matrix(0, nrow=length(years_fcast), ncol=length(mods)+1))
  names(cog_lat_fcast_y) <- c("year", mods)
  cog_lat_fcast_y$year <- unique(dat_fcast$year)

  
  # Calculate COG - hist period
  for (y in cog_lat_hist$year) {
    cog_lat_truth <- weighted.mean(dat_hist[dat_hist$year==y, "lat"],w=dat_hist[dat_hist$year==y, "abundance"])
    cog_lat_hist[cog_lat_hist$year==y, "truth"] <- cog_lat_truth
    for (m in mods) {
      cog_lat_mod <- weighted.mean(dat_hist[dat_hist$year==y, "lat"],w=dat_hist[dat_hist$year==y, m])
      cog_lat_hist[cog_lat_hist$year==y, m] <- cog_lat_mod
    }
  }
  
  # Calculate COG - fcast period
  for (y in cog_lat_fcast$year) {
    cog_lat_truth <- weighted.mean(dat_fcast[dat_fcast$year==y, "lat"],w=dat_fcast[dat_fcast$year==y, "abundance"])
    cog_lat_fcast[cog_lat_fcast$year==y, "truth"] <- cog_lat_truth
    for (m in mods) {
      cog_lat_mod <- weighted.mean(dat_fcast[dat_fcast$year==y, "lat"],w=dat_fcast[dat_fcast$year==y, m])
      cog_lat_fcast[cog_lat_fcast$year==y, m] <- cog_lat_mod
      
      #calculate difference between observed and predicted COG
      cog_lat_fcast_y[cog_lat_fcast_y$year == y, m] <- cog_lat_truth - cog_lat_mod
    }
  }
  
  # Plot COG - hist period
  par(mfrow=c(ceiling(length(mods)/3),3))
  for (m in mods) {
    rmse_mod <- RMSE(cog_lat_hist$truth, cog_lat_hist[,m])
    rmse_cog$rmse_hist[rmse_cog$model==m] <- rmse_mod
    
    plot(cog_lat_hist$year,cog_lat_hist$truth, type='l', ylab="Latitude", xlim=c(min(years_hist),max(years_hist)),
         main=paste0("COG-lat, Hist, ", m, ", rmse=", round(rmse_mod,3)), lwd=2)
    lines(cog_lat_hist$year,cog_lat_hist[,m], col="blue", lwd=2)
  }
  
  # Plot COG - fcast period
  par(mfrow=c(ceiling(length(mods)/3),3))
  for (m in mods) {
    rmse_mod <- RMSE(cog_lat_fcast$truth, cog_lat_fcast[,m])
    rmse_cog$rmse_fcast[rmse_cog$model==m] <- rmse_mod
    
    plot(cog_lat_fcast$year,cog_lat_fcast$truth, type='l', ylab="Latitude", xlim=c(min(years_fcast),max(years_fcast)),
         main=paste0("COG-lat, Fcast, ", m, ", rmse=", round(rmse_mod,3)), lwd=2)
    lines(cog_lat_fcast$year,cog_lat_fcast[,m], col="blue", lwd=2)
  }
  
  #plot obs-pred difference in COG during fcast period
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  for (m in mods) {
    plot(cog_lat_fcast_y$year, cog_lat_fcast_y[,m], type="l", lwd=2, col="red",
         ylab="obs-pred cog", xlab="", main=paste(m),
         ylim=c(min(cog_lat_fcast_y[,m]),max(cog_lat_fcast_y[,m])))
    abline(h=0, lty=2, col="blue")
  }

  return(list(rmse_cog=rmse_cog,
              cog_lat_hist=cog_lat_hist,
              cog_lat_fcast=cog_lat_fcast,
              cog_lat_fcast_y=cog_lat_fcast_y))
}

  
