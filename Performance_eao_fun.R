### Location^3 Estimation Model Code
## Function to measure Effective Area Occupied

## Code by Jim Thorson and James Smith - April 2020


sdm_eao <- function(dat_hist, dat_fcast) {
  # assumes that all grid cells are equal size, such that dat$abundance is also proportional to density
  # average_density has units Biomass per area
  # total_abundance has units Biomass
  # therefore: effective_area_occupied has units area
  # ***** this should probably be done with observed/predicted data in every grid cell, not just subset of cells
  
  all_mods <- c("gam_E","gam_S","gam_ES","gam_EST","gam_ECor",
                "brt_E","brt_S","brt_ES","brt_EST",
                "mlp_E","mlp_S","mlp_ES","mlp_EST",
                "glm_E","glm_Sr","glm_ESt","glm_ESr")
  
  mods <- names(dat_hist)[names(dat_hist) %in% all_mods]   #which of the possible models have been run
  
  eao_hist_y <- as.data.frame(matrix(0, nrow=length(unique(dat_hist$year)), ncol=length(mods)+2))
  names(eao_hist_y) <- c("year", "observed", mods)
  eao_hist_y$year <- unique(dat_hist$year)
  
  eao_fcast_y <- as.data.frame(matrix(0, nrow=length(unique(dat_fcast$year)), ncol=length(mods)+2))
  names(eao_fcast_y) <- c("year", "observed", mods)
  eao_fcast_y$year <- unique(dat_fcast$year)
  
  # eao_fcast_y <- as.data.frame(matrix(0, nrow=length(unique(dat_fcast$year)), ncol=length(mods)+1))
  # names(eao_fcast_y) <- c("year", mods)
  # eao_fcast_y$year <- unique(dat_fcast$year)
  
  #observed - hist
  for (yy in eao_hist_y$year) {
    dat_hist_yy <- dat_hist[dat_hist$year == yy,]
    average_density = weighted.mean( x=dat_hist_yy$abundance, w=dat_hist_yy$abundance )
    total_abundance = sum( dat_hist_yy$abundance )
    effective_area_occupied = total_abundance / average_density
    
    eao_hist_y[eao_hist_y$year == yy, "observed"] <- effective_area_occupied
  }
  
  #observed - fcast
  for (yy in eao_fcast_y$year) {
    dat_fcast_yy <- dat_fcast[dat_fcast$year == yy,]
    average_density = weighted.mean( x=dat_fcast_yy$abundance, w=dat_fcast_yy$abundance )
    total_abundance = sum( dat_fcast_yy$abundance )
    effective_area_occupied = total_abundance / average_density
    
    eao_fcast_y[eao_fcast_y$year == yy, "observed"] <- effective_area_occupied
  }
  
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  #modelled - hist
  for (m in mods) {
    for (yy in eao_hist_y$year) {
      
      dat_hist_yy <- dat_hist[dat_hist$year == yy,]
      average_density = weighted.mean( x=dat_hist_yy[,m], w=dat_hist_yy[,m] )
      total_abundance = sum( dat_hist_yy[,m] )
      effective_area_occupied = total_abundance / average_density
      
      eao_hist_y[eao_hist_y$year == yy, m] <- effective_area_occupied
    }
    plot(eao_hist_y$year, eao_hist_y$observed, type="l", lwd=2, xlab="", ylab="EAO",
         main=paste("Hist, ",m, ", cor=", round(cor(eao_hist_y$observed, eao_hist_y[,m]),3)),
         ylim=c(0, 400), las=1)
    lines(eao_hist_y$year, eao_hist_y[,m], col="red", lwd=2)
  }
  
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  #modelled - fcast
  for (m in mods) {
    for (yy in eao_fcast_y$year) {
      
      dat_fcast_yy <- dat_fcast[dat_fcast$year == yy,]
      average_density = weighted.mean( x=dat_fcast_yy[,m], w=dat_fcast_yy[,m] )
      total_abundance = sum( dat_fcast_yy[,m] )
      effective_area_occupied = total_abundance / average_density
      
      eao_fcast_y[eao_fcast_y$year == yy, m] <- effective_area_occupied
    }
    plot(eao_fcast_y$year, eao_fcast_y$observed, type="l", lwd=2, xlab="", ylab="EAO",
         main=paste("Fcast, ",m, ", cor=", round(cor(eao_fcast_y$observed, eao_fcast_y[,m]),3)),
         ylim=c(0, 400), las=1)
    lines(eao_fcast_y$year, eao_fcast_y[,m], col="red", lwd=2)
  }
  
  par(mfrow=c(ceiling(length(mods)/3),3), mar=c(2.5,4,2.5,1))
  #plot obs-pred difference in EAO over time
  for (m in mods) {
    diffx <- eao_fcast_y$observed - eao_fcast_y[,m]
    plot(eao_fcast_y$year, diffx, type="l", lwd=2, col="red", xlab="", ylab="obs-pred EAO",
         main=paste(m), ylim=c(min(diffx), max(diffx)))
    abline(h=0, lty=2, col="blue")
    
  }
  
  return(list(eao_hist_y = eao_hist_y,
              eao_fcast_y = eao_fcast_y))
  
}
