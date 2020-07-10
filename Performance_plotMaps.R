### Location^3 Estimation Model Code
## Create maps of prediction species distribution - observed and fitted SDMs

## Code by James Smith - April 2020

performance_plotMaps_function <- function(om = om){
  
  library(raster)
  
  # prediction function for sdmTMB glmms
  if(om == "hms"){
    predict_glmm <- function(model) {
      dummy = data.frame(year = unique(dat_all$year),
                         lat=dat_all$lat[1],
                         lon=unique(dat_all$lon[1]),
                         temp_s = dat_all$temp_s[1],
                         chla_s = dat_all$chla_s[1],
                         mld_s = dat_all$mld_s[1])
      pred = predict(model,
                     newdata=rbind(dat_x[,c("year","lat","lon","temp_s","chla_s",
                                            "mld_s")],
                                   dummy), xy_cols = c("lon", "lat"))
      
      # drop dummy data
      pred = pred[-c(seq(nrow(pred)-nrow(dummy)+1,nrow(pred))),]
      pred$abundance = dat_x$abundance # add true abundance
      
      return(pred)
    }
  }
  
  if(om == "cps"){
    predict_glmm <- function(model) {
      dummy = data.frame(year = unique(dat_all$year),
                         lat=dat_all$lat[1],
                         lon=unique(dat_all$lon[1]),
                         temp_s = dat_all$temp_s[1],
                         chla_s = dat_all$chla_s[1],
                         z_s = dat_all$z_s[1])
      # z_s = dat_all$z_s[1])
      pred = predict(model,
                     newdata=rbind(dat_x[,c("year","lat","lon","temp_s","chla_s",
                                            "z_s")],
                                   dummy), xy_cols = c("lon", "lat"))
      
      # drop dummy data
      pred = pred[-c(seq(nrow(pred)-nrow(dummy)+1,nrow(pred))),]
      pred$abundance = dat_x$abundance # add true abundance
      
      return(pred)
    }
  }
  
  if(om == "grf"){
    predict_glmm <- function(model) {
      dummy = data.frame(year = unique(dat_all$year),
                         lat=dat_all$lat[1],
                         lon=unique(dat_all$lon[1]),
                         btemp_s = dat_all$btemp_s[1],
                         O2_s = dat_all$O2_s[1],
                         z_s = dat_all$z_s[1])
      # z_s = dat_all$z_s[1])
      pred = predict(model,
                     newdata=rbind(dat_x[,c("year","lat","lon","btemp_s","O2_s",
                                            "z_s")],
                                   dummy), xy_cols = c("lon", "lat"))
      
      # drop dummy data
      pred = pred[-c(seq(nrow(pred)-nrow(dummy)+1,nrow(pred))),]
      pred$abundance = dat_x$abundance # add true abundance
      
      return(pred)
    }
  }
  
  
  #GAMs
  par(mfrow=c(3,5), mar=c(2,3,3,4))
  for (yy in c(2000, 2040, 2080)) {
    
    year_x <- yy
    dat_x <- dat_all[dat_all$year==year_x,]
    
    r_obs <- rasterFromXYZ(dat_x[,c("lon","lat","abundance")])
    plot(r_obs, asp=1, main=paste0("Observed, year ", year_x))
    
    presx <- predict(gam_E_P, dat_x, type="response")
    abundx <- exp(predict(gam_E_N, dat_x, type="response"))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, gam_E, year ", year_x))
    
    presx <- predict(gam_ES_P, dat_x, type="response")
    abundx <- exp(predict(gam_ES_N, dat_x, type="response"))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, gam_ES, year ", year_x))
    
    presx <- predict(gam_EST_P, dat_x, type="response")
    abundx <- exp(predict(gam_EST_N, dat_x, type="response"))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, gam_EST, year ", year_x))
    
    if (length(ls(pattern = "gam_ECor")) > 0) {
      presx <- predict(gam_ECor_P$gam, dat_x, type="response")
      abundx <- exp(predict(gam_ECor_N$gam, dat_x, type="response"))
      dat_x$pred_x <- presx * abundx
      r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
      plot(r_pred, asp=1, main=paste0("Predicted, gam_ECor, year ", year_x))
    } else {
      plot(r_obs <- 0, main="NA")
    }
    
  }
  
  #GLMs
  par(mfrow=c(3,5), mar=c(2,3,3,4))
  for (yy in c(2000, 2040, 2080)) {
    
    year_x <- yy
    dat_x <- dat_all[dat_all$year==year_x,]
    
    r_obs <- rasterFromXYZ(dat_x[,c("lon","lat","abundance")])
    plot(r_obs, asp=1, main=paste0("Observed, year ", year_x))
    
    presx <- predict_glmm(glm_E_P)
    presx <- exp(presx$est)/(1 + exp(presx$est))
    abundx <- predict_glmm(glm_E_N)
    dat_x$pred_x <- presx * exp(abundx$est)
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, glm_E, year ", year_x))
    
    presx <- predict_glmm(glm_ESt_P)
    presx <- exp(presx$est)/(1 + exp(presx$est))
    abundx <- predict_glmm(glm_ESt_N)
    dat_x$pred_x <- presx * exp(abundx$est)
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, glm_ESt, year ", year_x))
    
    if (length(ls(pattern = "glm_ESr")) > 0) {
      presx <- predict_glmm(glm_ESr_P)
      presx <- exp(presx$est)/(1 + exp(presx$est))
      abundx <- predict_glmm(glm_ESr_N)
      dat_x$pred_x <- presx * exp(abundx$est)
      r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
      plot(r_pred, asp=1, main=paste0("Predicted, glm_ESr, year ", year_x))
    } else {
      plot(r_obs <- 0, main="NA")
    }
    
    if (length(ls(pattern = "glm_Sr")) > 0) {
      presx <- predict_glmm(glm_Sr_P)
      presx <- exp(presx$est)/(1 + exp(presx$est))
      abundx <- predict_glmm(glm_Sr_N)
      dat_x$pred_x <- presx * exp(abundx$est)
      r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
      plot(r_pred, asp=1, main=paste0("Predicted, glm_Sr, year ", year_x))
    } else {
      plot(r_obs <- 0, main="NA")
    }
  }
  
  #BRTs
  par(mfrow=c(3,4), mar=c(2,3,3,4))
  for (yy in c(2000, 2040, 2080)) {
    
    year_x <- yy
    dat_x <- dat_all[dat_all$year==year_x,]
    
    r_obs <- rasterFromXYZ(dat_x[,c("lon","lat","abundance")])
    plot(r_obs, asp=1, main=paste0("Observed, year ", year_x))
    
    presx <- predict(brt_E_P, dat_x, type="response", n.trees=brt_E_P$gbm.call$best.trees)
    abundx <- exp(predict(brt_E_N, dat_x, type="response", n.trees=brt_E_N$gbm.call$best.trees))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, brt_E, year ", year_x))
    
    presx <- predict(brt_ES_P, dat_x, type="response", n.trees=brt_ES_P$gbm.call$best.trees)
    abundx <- exp(predict(brt_ES_N, dat_x, type="response", n.trees=brt_ES_N$gbm.call$best.trees))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, brt_ES, year ", year_x))
    
    presx <- predict(brt_EST_P, dat_x, type="response", n.trees=brt_EST_P$gbm.call$best.trees)
    abundx <- exp(predict(brt_EST_N, dat_x, type="response", n.trees=brt_EST_N$gbm.call$best.trees))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, brt_EST, year ", year_x))
    
  }
  
  #MLPs
  par(mfrow=c(3,4), mar=c(2,3,3,4))
  for (yy in c(2000, 2040, 2080)) {
    
    year_x <- yy
    dat_x <- dat_all[dat_all$year==year_x,]
    
    r_obs <- rasterFromXYZ(dat_x[,c("lon","lat","abundance")])
    plot(r_obs, asp=1, main=paste0("Observed, year ", year_x))
    
    presx <- predict(mlp_E_P, dat_x, type="response")
    abundx <- exp(predict(mlp_E_N, dat_x, type="response"))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, mlp_E, year ", year_x))
    
    presx <- predict(mlp_ES_P, dat_x, type="response")
    abundx <- exp(predict(mlp_ES_N, dat_x, type="response"))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, mlp_ES, year ", year_x))
    
    presx <- predict(mlp_EST_P, dat_x, type="response")
    abundx <- exp(predict(mlp_EST_N, dat_x, type="response"))
    dat_x$pred_x <- presx * abundx
    r_pred <- rasterFromXYZ(dat_x[,c("lon","lat","pred_x")])
    plot(r_pred, asp=1, main=paste0("Predicted, mlp_EST, year ", year_x))
    
  }
}
