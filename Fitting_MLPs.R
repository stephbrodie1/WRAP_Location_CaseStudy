### Location^3 Estimation Model Code
## Function to run MLPs on data from OM
## Function returns only the fitted and predicted values
## Uses 'env_formula' as a base combination of enviro covariates to build from

## Code by James Smith, Apr 2020

#data was normalised first - IT IS ESSENTIAL THAT ALL DATA, FOR FITTING AND PREDICTING, IS NORMALISED AT SAME TIME (i.e. has the same range, or mean and SD)
#JS: I changed type of normalisation to 'range', bc lat, lon and year don't really make sense as scaled variables

#'type' is only 'delta' at this stage
#'covs' is a vector of model covariates, where:
#     - E is enviro only
#     - S is space only
#     - ES is enviro + space
#     - EST is enviro + space + time


delta1_formula <- formula(paste("pres ~", env_formula_MLP))
delta2_formula <- formula(paste("log_abundance ~", env_formula_MLP))

if (type == "delta") {
  
  if ("E" %in% covs) {
    print("Fitting MLP-E")
    mlp_E_P <- neuralnet(delta1_formula, data = dat_hist, 
                         hidden = c(3), linear.output = F, algorithm = "rprop+", threshold = 0.2)
    mlp_E_N <- neuralnet(delta2_formula, data = subset(dat_hist, !is.infinite(dat_hist$log_abundance)), 
                         hidden = c(3), linear.output = T, algorithm = "rprop+", threshold = 0.6)
    
    presx <- as.numeric(predict(mlp_E_P,dat_hist))
    abundx <- as.numeric(predict(mlp_E_N,dat_hist))
    dat_hist$mlp_E <- presx * exp(abundx)
    # par(mfrow=c(1,2))
    # plot(dat_hist$temp, presx, main="MLP-E-Pres", xlab="Temp")
    # plot(dat_hist$temp, abundx, main="MLP-E-Abund", xlab="Temp")
    
    presx <- as.numeric(predict(mlp_E_P,dat_fcast))
    abundx <- as.numeric(predict(mlp_E_N,dat_fcast))
    dat_fcast$mlp_E <- presx * exp(abundx)
  }
  
  if ("S" %in% covs) {
    print("Fitting MLP-S")
    mlp_S_P <- neuralnet(pres ~ lat_n + lon_n, data = dat_hist, 
                         hidden = c(3), linear.output = F, algorithm = "rprop+", threshold = 0.2)
    mlp_S_N <- neuralnet(log_abundance ~ lat_n + lon_n, data = subset(dat_hist, !is.infinite(dat_hist$log_abundance)), 
                         hidden = c(3), linear.output = T, algorithm = "rprop+", threshold = 0.6)
    
    presx <- as.numeric(predict(mlp_S_P,dat_hist))
    abundx <- as.numeric(predict(mlp_S_N,dat_hist))
    dat_hist$mlp_S <- presx * exp(abundx)
    # par(mfrow=c(1,2))
    # plot(dat_hist$temp, presx, main="MLP-S-Pres", xlab="Temp")
    # plot(dat_hist$temp, abundx, main="MLP-S-Abund", xlab="Temp")
    
    presx <- as.numeric(predict(mlp_S_P,dat_fcast))
    abundx <- as.numeric(predict(mlp_S_N,dat_fcast))
    dat_fcast$mlp_S <- presx * exp(abundx)
  }
  
  if ("ES" %in% covs) {
    print("Fitting MLP-ES")
    mlp_ES_P <- neuralnet(update(delta1_formula, ~. + lat_n + lon_n), data = dat_hist, 
                         hidden = c(3), linear.output = F, algorithm = "rprop+", threshold = 0.2)
    mlp_ES_N <- neuralnet(update(delta2_formula, ~. + lat_n + lon_n), data = subset(dat_hist, !is.infinite(dat_hist$log_abundance)), 
                         hidden = c(3), linear.output = T, algorithm = "rprop+", threshold = 0.6) #changing to 0.3 for anchovy temp only models
    
    presx <- as.numeric(predict(mlp_ES_P,dat_hist))
    abundx <- as.numeric(predict(mlp_ES_N,dat_hist))
    dat_hist$mlp_ES <- presx * exp(abundx)
    # par(mfrow=c(1,2))
    # plot(dat_hist$temp, presx, main="MLP-ES-Pres", xlab="Temp")
    # plot(dat_hist$temp, abundx, main="MLP-ES-Abund", xlab="Temp")
    
    presx <- as.numeric(predict(mlp_ES_P,dat_fcast))
    abundx <- as.numeric(predict(mlp_ES_N,dat_fcast))
    dat_fcast$mlp_ES <- presx * exp(abundx)
  }
  
  if ("EST" %in% covs) {
    print("Fitting MLP-EST")
    mlp_EST_P <- neuralnet(update(delta1_formula, ~. + lat_n + lon_n + year_n), data = dat_hist, 
                          hidden = c(3), linear.output = F, algorithm = "rprop+", threshold = 0.2)
    mlp_EST_N <- neuralnet(update(delta2_formula, ~. + lat_n + lon_n + year_n), data = subset(dat_hist, !is.infinite(dat_hist$log_abundance)), 
                          hidden = c(3), linear.output = T, algorithm = "rprop+", threshold = 0.5)  #changing to 0.3 for anchovy temp only models
    
    presx <- as.numeric(predict(mlp_EST_P, dat_hist))
    abundx <- as.numeric(predict(mlp_EST_N, dat_hist))
    dat_hist$mlp_EST <- presx * exp(abundx)
    # par(mfrow=c(1,2))
    # plot(dat_hist$temp, presx, main="MLP-EST-Pres", xlab="Temp")
    # plot(dat_hist$temp, abundx, main="MLP-EST-Abund", xlab="Temp")
    
    presx <- as.numeric(predict(mlp_EST_P,dat_fcast))
    abundx <- as.numeric(predict(mlp_EST_N,dat_fcast))
    dat_fcast$mlp_EST <- presx * exp(abundx)
  }
  
  ## we can adapt the code below if we get fitting issues
  # test.x <- try(predict(mlp_EST_P, dat_hist))
  # while (inherits(test.x, 'try-error') | diff(range(test.x)) < 0.1) {
  #   print("mlp_EST_P failed")
  #   mlp_EST_P <- neuralnet(update(delta1_formula, ~. + lat_n + lon_n + year_n), data = dat_hist, 
  #                          hidden = c(3), linear.output = F, algorithm = "rprop+", threshold = 0.2)
  #   test.x <- try(predict(mlp_EST_P, dat_hist))
  # }
  
  
  ## Plots - similar to partial effects
  # newdat_temp <- data.frame(temp_raw=dat$temp, temp_n=dat$temp_n,
  #                           chla_n=median(dat$chla_n),
  #                           mld_n=median(dat$mld_n),
  #                           lat_n=median(dat$lat_n),
  #                           lon_n=median(dat$lon_n),
  #                           year_n=median(dat$year_n))
  # newdat_temp <- newdat_temp[order(newdat_temp$temp_n),]

  #pres-abs; temp
  # presT1 <- predict(mlp_E_P, newdata=newdat_temp)
  # #presT2 <- predict(mlp_S_P, newdata=newdat_temp)
  # presT3 <- predict(mlp_ES_P, newdata=newdat_temp)
  # presT4 <- predict(mlp_EST_P, newdata=newdat_temp)
  # par(mfrow=c(2,2))
  # plot(newdat_temp$temp_raw, presT1, main="mlp-E, Pres-Abs", xlab="Temp", type="l")
  # plot(newdat_temp$temp_raw, presT3, main="mlp-ES, Pres-Abs", xlab="Temp", type="l")
  # plot(newdat_temp$temp_raw, presT4, main="mlp-EST, Pres-Abs", xlab="Temp", type="l")
  # 
  # #abundance; temp
  # abundT1 <- predict(mlp_E_N, newdata=newdat_temp)
  # #abundT2 <- predict(mlp_S_N, newdata=newdat_temp)
  # abundT3 <- predict(mlp_ES_N, newdata=newdat_temp)
  # abundT4 <- predict(mlp_EST_N, newdata=newdat_temp)
  # par(mfrow=c(2,2))
  # plot(newdat_temp$temp_raw, abundT1, main="mlp-E, Abund", xlab="Temp", type="l")
  # plot(newdat_temp$temp_raw, abundT3, main="mlp-ES, Abund", xlab="Temp", type="l")
  # plot(newdat_temp$temp_raw, abundT4, main="mlp-EST, Abund", xlab="Temp", type="l")
  # 
  # 
  # newdat_chla <- data.frame(chla_raw=dat$chla, chla_n=dat$chla_n,
  #                           temp_n=median(dat$temp_n),
  #                           mld_n=median(dat$mld_n),
  #                           lat_n=median(dat$lat_n),
  #                           lon_n=median(dat$lon_n),
  #                           year_n=median(dat$year_n))
  # newdat_chla <- newdat_chla[order(newdat_chla$chla_n),]
  # 
  # #pres-abs; chl
  # presC1 <- predict(mlp_E_P, newdata=newdat_chla)
  # #presC2 <- predict(mlp_S_P, newdata=newdat_chla)
  # presC3 <- predict(mlp_ES_P, newdata=newdat_chla)
  # presC4 <- predict(mlp_EST_P, newdata=newdat_chla)
  # par(mfrow=c(2,2))
  # plot(newdat_chla$chla_raw, presC1, main="mlp-E, Pres-Abs", xlab="Chla", type="l")
  # plot(newdat_chla$chla_raw, presC3, main="mlp-ES, Pres-Abs", xlab="Chla", type="l")
  # plot(newdat_chla$chla_raw, presC4, main="mlp-EST, Pres-Abs", xlab="Chla", type="l")
  # 
  # #abund; chl
  # abundC1 <- predict(mlp_E_N, newdata=newdat_chla)
  # #abundC2 <- predict(mlp_S_N, newdata=newdat_chla)
  # abundC3 <- predict(mlp_ES_N, newdata=newdat_chla)
  # abundC4 <- predict(mlp_EST_N, newdata=newdat_chla)
  # par(mfrow=c(2,2))
  # plot(newdat_chla$chla_raw, abundC1, main="mlp-E, Abund", xlab="Chla", type="l")
  # plot(newdat_chla$chla_raw, abundC3, main="mlp-ES, Abund", xlab="Chla", type="l")
  # plot(newdat_chla$chla_raw, abundC4, main="mlp-EST, Abund", xlab="Chla", type="l")
  # 
  # newdat_mld <- data.frame(mld_raw=dat$mld, mld_n=dat$mld_n,
  #                           temp_n=median(dat$temp_n),
  #                           chla_n=median(dat$chla_n),
  #                           lat_n=median(dat$lat_n),
  #                           lon_n=median(dat$lon_n),
  #                           year_n=median(dat$year_n))
  # newdat_mld <- newdat_mld[order(newdat_mld$mld_n),]
  # 
  # #pres-abs; mld
  # presM1 <- predict(mlp_E_P, newdata=newdat_mld)
  # #presM2 <- predict(mlp_S_P, newdata=newdat_mld)
  # presM3 <- predict(mlp_ES_P, newdata=newdat_mld)
  # presM4 <- predict(mlp_EST_P, newdata=newdat_mld)
  # par(mfrow=c(2,2))
  # plot(newdat_mld$mld_raw, presM1, main="mlp-E, Pres-Abs", xlab="mld", type="l")
  # plot(newdat_mld$mld_raw, presM3, main="mlp-ES, Pres-Abs", xlab="mld", type="l")
  # plot(newdat_mld$mld_raw, presM4, main="mlp-EST, Pres-Abs", xlab="mld", type="l")
  # 
  # #abund; mld
  # abundM1 <- predict(mlp_E_N, newdata=newdat_mld)
  # #abundM2 <- predict(mlp_S_N, newdata=newdat_mld)
  # abundM3 <- predict(mlp_ES_N, newdata=newdat_mld)
  # abundM4 <- predict(mlp_EST_N, newdata=newdat_mld)
  # par(mfrow=c(2,2))
  # plot(newdat_mld$mld_raw, abundM1, main="mlp-E, Abund", xlab="mld", type="l")
  # plot(newdat_mld$mld_raw, abundM3, main="mlp-ES, Abund", xlab="mld", type="l")
  # plot(newdat_mld$mld_raw, abundM4, main="mlp-EST, Abund", xlab="mld", type="l")

}

#NeuralNetTools::plotnet(mlp_E_P)  #an evaluation plot
#*** i also tried verrrry hard to get 'lekprofile' fom this pkg working, for plotting partial effects, but always got 'Error in 1:nrow(const) : argument of length 0', and i don't know why
#mod <- nnet(pres ~ X1 + X2 + X3, data = neuraldat, size = 5)  #an alternative MLP pkg




