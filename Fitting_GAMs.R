### Location^3 Estimation Model Code
## Function to run GAMs on data from OM
## Function returns only the fitted and predicted values
## Uses 'env_formula' as a base combination of enviro covariates to build from

## Code by James Smith, Apr 2020


  #'dat_hist' are the observations from the OM used to fit the models
  #'dat_fcast' are the observations from the OM used to forecast the models
  #'type' is either 'tweedie' or 'delta' ***tweedie can take longer to fit?
  #'covs' is a vector of model covariates, where:
  #         - E is enviro only
  #         - S is space only
  #         - ES is enviro + space
  #         - EST is enviro + space-time
  #         - ECor is enviro + corGaus autocorrelation
  #         - c("E", "S", "ES", "EST", "ECor") does them all
  #'env_formula' is the formula (as character) for enviro covariates
  

  tw_formula <- formula(paste("abundance ~", env_formula))
  delta1_formula <- formula(paste("pres ~", env_formula))
  delta2_formula <- formula(paste("log_abundance ~", env_formula))
  
    if ("E" %in% covs) {  #Enviro only
    print("Fitting GAM-E")
    if (type == "tweedie") {
      gam_E <- gam(tw_formula, data=dat_hist, family=tw(link="log"), method="REML")
      #plot(gam_E); AIC(gam_E); summary(gam_E); gam.check(gam_E)
      saved_theta <- gam_E$family$getTheta(TRUE)  #p=1.031
      
      dat_hist$gam_E <- predict(gam_E, dat_hist, type="response")  #save historical values
      dat_fcast$gam_E <- predict(gam_E, dat_fcast, type="response")  #save forecast values
    }
    if (type == "delta") {
      gam_E_P <- gam(delta1_formula, data=dat_hist, family=binomial)
      #plot(gam_E_P, pages=1)
      gam_E_N <- gam(delta2_formula, data=dat_hist[dat_hist$abundance>0,], family=gaussian)
      #plot(gam_E_N, pages=1)
      
      presx <- predict(gam_E_P, dat_hist, type="response")
      abundx <- exp(predict(gam_E_N, dat_hist, type="response"))
      dat_hist$gam_E <- presx * abundx
      
      presx <- predict(gam_E_P, dat_fcast, type="response")
      abundx <- exp(predict(gam_E_N, dat_fcast, type="response"))
      dat_fcast$gam_E <- presx * abundx
      
      # presx <- predict(gam_E_P, dat, type="response")
      # abundx <- exp(predict(gam_E_N, dat, type="response"))
      # dat$gam_E <- presx * abundx
    }
  }
  
  if ("S" %in% covs) {  #Space surface only
    print("Fitting GAM-S")
    if (type == "tweedie") {
      gam_S <- gam(abundance ~ s(lat,lon), data=dat_hist, family=tw(link="log"), method="REML")
      #plot(gam_ST, pages=1)
      
      dat_hist$gam_S <- predict(gam_S, dat_hist, type="response")
      dat_fcast$gam_S <- predict(gam_S, dat_fcast, type="response")
    }
    if (type == "delta") {
      gam_S_P <- gam(pres ~ s(lat,lon), data=dat_hist, family=binomial)
      #plot(gam_ST_P, pages=1)
      gam_S_N <- gam(log_abundance ~ s(lat,lon), data=dat_hist[dat_hist$abundance>0,], family=gaussian)
      #plot(gam_ST_N, pages=1)
      
      presx <- predict(gam_S_P, dat_hist, type="response")
      abundx <- exp(predict(gam_S_N, dat_hist, type="response"))
      dat_hist$gam_S <- presx * abundx
      
      presx <- predict(gam_S_P, dat_fcast, type="response")
      abundx <- exp(predict(gam_S_N, dat_fcast, type="response"))
      dat_fcast$gam_S <- presx * abundx
    }
  }
  
  if ("ES" %in% covs) {  #Enviro and Space smoother
    print("Fitting GAM-ES")
    if (type == "tweedie") {
      gam_ES <- gam(update(tw_formula, ~. + s(lat,lon)), data=dat_hist, family=tw(link="log"), method="REML")
      #plot(gam_ES, pages=1); AIC(gam_ES)
      
      dat_hist$gam_ES <- predict(gam_ES, dat_hist, type="response")
      dat_fcast$gam_ES <- predict(gam_ES, dat_fcast, type="response")
    }
    if (type == "delta") {
      gam_ES_P <- gam(update(delta1_formula, ~. + s(lat,lon)), data=dat_hist, family=binomial)
      #plot(gam_ES_P)
      gam_ES_N <- gam(update(delta2_formula, ~. + s(lat,lon)), data=dat_hist[dat_hist$abundance>0,], family=gaussian)
      #plot(gam_ES_N)
      
      presx <- predict(gam_ES_P, dat_hist, type="response")
      abundx <- exp(predict(gam_ES_N, dat_hist, type="response"))
      dat_hist$gam_ES <- presx * abundx
      
      presx <- predict(gam_ES_P, dat_fcast, type="response")
      abundx <- exp(predict(gam_ES_N, dat_fcast, type="response"))
      dat_fcast$gam_ES <- presx * abundx
    }
  }
  
  if ("EST" %in% covs) {  #Enviro and Space-Time tensor
    print("Fitting GAM-EST")
    if (type == "tweedie") {
      gam_EST <- gam(update(tw_formula, ~. + te(lat,lon,year)), data=dat_hist, family=tw(link="log"), method="REML")
      #plot(gam_EST, pages=1); AIC(gam_EST)
      
      dat_hist$gam_EST <- predict(gam_EST, dat_hist, type="response")
      dat_fcast$gam_EST <- predict(gam_EST, dat_fcast, type="response")
    }
    if (type == "delta") {
      gam_EST_P <- gam(update(delta1_formula, ~. + te(lat,lon,year)), data=dat_hist, family=binomial)
      #plot(gam_EST_P, pages=1)
      gam_EST_N <- gam(update(delta2_formula, ~. + te(lat,lon,year)), data=dat_hist[dat_hist$abundance>0,], family=gaussian)
      #plot(gam_EST_N, pages=1)
      
      presx <- predict(gam_EST_P, dat_hist, type="response")
      abundx <- exp(predict(gam_EST_N, dat_hist, type="response"))
      dat_hist$gam_EST <- presx * abundx
      
      presx <- predict(gam_EST_P, dat_fcast, type="response")
      abundx <- exp(predict(gam_EST_N, dat_fcast, type="response"))
      dat_fcast$gam_EST <- presx * abundx
    }
  }
  
  if ("ECor" %in% covs) {  #Enviro and Space-Time residual correlation
    print("Fitting GAM-ECor")
    if (type == "tweedie") {
      gam_ECor <- gamm(tw_formula, correlation=corGaus(form=~lat+lon|fYear), data=dat_hist,
                                  family=Tweedie(p=saved_theta, link="log"))
      #plot(gam_ECor$gam, pages=1)
      
      dat_hist$gam_ECor <- predict(gam_ECor, dat_hist, type="response")
      dat_fcast$gam_ECor <- predict(gam_ECor, dat_fcast, type="response")
    }
    if (type == "delta") {
      gam_ECor_P <- gamm(delta1_formula, correlation=corGaus(form=~lat+lon|fYear),
                         data=dat_hist, family=binomial)
      #plot(gam_ECor_P$gam, pages=1)
      gam_ECor_N <- gamm(delta2_formula, correlation=corGaus(form=~lat+lon|fYear),
                        data=dat_hist[dat_hist$abundance>0,], family=gaussian)
      #plot(gam_ECor_N$gam, pages=1)
      
      presx <- predict(gam_ECor_P$gam, dat_hist, type="response")
      abundx <- exp(predict(gam_ECor_N$gam, dat_hist, type="response"))
      dat_hist$gam_ECor <- presx * abundx
      
      presx <- predict(gam_ECor_P$gam, dat_fcast, type="response")
      abundx <- exp(predict(gam_ECor_N$gam, dat_fcast, type="response"))
      dat_fcast$gam_ECor <- presx * abundx
    }
  }
  
  if (type == "tweedie") {
    par(mfrow=c(3,2))
    plot(gam_E, select=1, main="gam_E, Tweedie", scale=0)
    plot(gam_S, select=1, main="gam_S, Tweedie", scheme=2, rug=F, scale=0)
    plot(gam_ES, select=1, main="gam_ES, Tweedie", scale=0)
    plot(gam_EST, select=1, main="gam_EST, Tweedie", scale=0)
    try(plot(gam_ECor$gam, select=1, main="gam_ECor, Tweedie", scale=0))
    
    par(mfrow=c(2,2))
    plot(gam_E, select=2, main="gam_E, Tweedie", scale=0)
    #plot(gam_S, select=2, main="gam_S, Tweedie", scheme=2, rug=F)
    plot(gam_ES, select=2, main="gam_ES, Tweedie", scale=0)
    plot(gam_EST, select=2, main="gam_EST, Tweedie", scale=0)
    try(plot(gam_ECor$gam, select=2, main="gam_ECor, Tweedie", scale=0))
    
    par(mfrow=c(2,2))
    plot(gam_E, select=3, main="gam_E, Tweedie", scale=0)
    #plot(gam_S, select=3, main="gam_S, Tweedie", scheme=2, rug=F)
    plot(gam_ES, select=3, main="gam_ES, Tweedie", scale=0)
    plot(gam_EST, select=3, main="gam_EST, Tweedie", scale=0)
    try(plot(gam_ECor$gam, select=3, main="gam_ECor, Tweedie", scale=0))
  }
  if (type == "delta") {
    par(mfrow=c(3,2))
    plot(gam_E_P, select=1, main="gam_E, Pres-Abs", scale=0)
    plot(gam_S_P, select=1, main="gam_S, Pres-Abs", scheme=2, rug=F, scale=0)
    plot(gam_ES_P, select=1, main="gam_ES, Pres-Abs", scale=0)
    plot(gam_EST_P, select=1, main="gam_EST, Pres-Abs", scale=0)
    try(plot(gam_ECor_P$gam, select=1, main="gam_ECor, Pres-Abs", scale=0))
    par(mfrow=c(3,2))
    plot(gam_E_N, select=1, main="gam_E, Abund", scale=0)
    plot(gam_S_N, select=1, main="gam_S, Abund", scheme=2, rug=F, scale=0)
    plot(gam_ES_N, select=1, main="gam_ES, Abund", scale=0)
    plot(gam_EST_N, select=1, main="gam_EST, Abund", scale=0)
    try(plot(gam_ECor_N$gam, select=1, main="gam_ECor, Abund", scale=0))
    
    par(mfrow=c(2,2))
    plot(gam_E_P, select=2, main="gam_E, Pres-Abs", scale=0)
    #plot(gam_S_P, select=2, main="gam_S, Pres-Abs")
    plot(gam_ES_P, select=2, main="gam_ES, Pres-Abs", scale=0)
    plot(gam_EST_P, select=2, main="gam_EST, Pres-Abs", scale=0)
    try(plot(gam_ECor_P$gam, select=2, main="gam_ECor, Pres-Abs", scale=0))
    par(mfrow=c(2,2))
    plot(gam_E_N, select=2, main="gam_E, Abund", scale=0)
    #plot(gam_S_N, select=2, main="gam_S, Abund")
    plot(gam_ES_N, select=2, main="gam_ES, Abund", scale=0)
    plot(gam_EST_N, select=2, main="gam_EST, Abund", scale=0)
    try(plot(gam_ECor_N$gam, select=2, main="gam_ECor, Abund", scale=0))
    
    par(mfrow=c(2,2))
    plot(gam_E_P, select=3, main="gam_E, Pres-Abs", scale=0)
    #plot(gam_S_P, select=3, main="gam_S, Pres-Abs")
    plot(gam_ES_P, select=3, main="gam_ES, Pres-Abs", scale=0)
    plot(gam_EST_P, select=3, main="gam_EST, Pres-Abs", scale=0)
    try(plot(gam_ECor_P$gam, select=3, main="gam_ECor, Pres-Abs", scale=0))
    par(mfrow=c(2,2))
    plot(gam_E_N, select=3, main="gam_E, Abund", scale=0)
    #plot(gam_S_N, select=3, main="gam_S, Abund", scale=0)
    plot(gam_ES_N, select=3, main="gam_ES, Abund", scale=0)
    plot(gam_EST_N, select=3, main="gam_EST, Abund", scale=0)
    try(plot(gam_ECor_N$gam, select=3, main="gam_ECor, Abund", scale=0))
  }
  par(mfrow=c(1,1))



