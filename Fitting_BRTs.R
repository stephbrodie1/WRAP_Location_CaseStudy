### Location^3 Estimation Model Code
## Function to run BRTs on data from OM
## Function returns only the fitted and predicted values
## Uses 'env_formula' as a base combination of enviro covariates to build from

## Code by James Smith, Apr 2020

#'dat_hist' are the observations from the OM used to fit the models
#'dat_fcast' are the observations from the OM used to forecast the models
#'type' is either 'tweedie' or 'delta'
#'covs' is a vector of model covariates, where:
#     - E is enviro only
#     - S is space only
#     - ES is enviro + space
#     - EST is enviro + space + time

gbm.x.env <- env_formula_BRT

if ("E" %in% covs) {
  print("Fitting BRT-E")
  if (type == "delta") {
    brt_E_P <- gbm.step(data=dat_hist,
                        gbm.x = gbm.x.env,
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.6, #making Lr0.001 for temp only anchovy
                        plot.main=FALSE, verbose = FALSE)
    
    brt_E_N <- gbm.step(data=dat_hist[dat_hist$abundance>0,], 
                       gbm.x = gbm.x.env,
                       gbm.y = 'log_abundance',
                       family = "gaussian",
                       tree.complexity = 3, learning.rate = 0.001, bag.fraction = 0.6,
                       plot.main=FALSE, verbose = FALSE)
    
    presx <- predict(brt_E_P, dat_hist, n.trees=brt_E_P$gbm.call$best.trees, type="response")
    abundx <- exp(predict(brt_E_N, dat_hist, n.trees=brt_E_N$gbm.call$best.trees, type="response"))
    dat_hist$brt_E <- presx * abundx
    
    presx <- predict(brt_E_P, dat_fcast, n.trees=brt_E_P$gbm.call$best.trees, type="response")
    abundx <- exp(predict(brt_E_N, dat_fcast, n.trees=brt_E_N$gbm.call$best.trees, type="response"))
    dat_fcast$brt_E <- presx * abundx
  }
}

if ("S" %in% covs) {
  print("Fitting BRT-S")
  if (type == "delta") {
    brt_S_P <- gbm.step(data=dat_hist,
                        gbm.x = c("lat", "lon"),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
    
    brt_S_N <- gbm.step(data=dat_hist[dat_hist$abundance>0,], 
                        gbm.x = c("lat", "lon"),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,
                        plot.main=FALSE, verbose = FALSE)
    
    presx <- predict(brt_S_P, dat_hist, n.trees=brt_S_P$gbm.call$best.trees, type="response")
    abundx <- exp(predict(brt_S_N, dat_hist, n.trees=brt_S_N$gbm.call$best.trees, type="response"))
    dat_hist$brt_S <- presx * abundx
    
    presx <- predict(brt_S_P, dat_fcast, n.trees=brt_S_P$gbm.call$best.trees, type="response")
    abundx <- exp(predict(brt_S_N, dat_fcast, n.trees=brt_S_N$gbm.call$best.trees, type="response"))
    dat_fcast$brt_S <- presx * abundx
  }
}

if ("ES" %in% covs) {
  print("Fitting BRT-ES")
  if (type == "delta") {
    brt_ES_P <- gbm.step(data=dat_hist,
                        gbm.x = c(gbm.x.env, "lat", "lon"),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,#made 0.1 for anchovy  temp only EMs
                        plot.main=FALSE, verbose = FALSE)
    
    brt_ES_N <- gbm.step(data=dat_hist[dat_hist$abundance>0,], 
                        gbm.x = c(gbm.x.env, "lat", "lon"),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,#made 0.1 for anchovy  temp only EMs
                        plot.main=FALSE, verbose = FALSE)
    
    presx <- predict(brt_ES_P, dat_hist, n.trees=brt_ES_P$gbm.call$best.trees, type="response")
    abundx <- exp(predict(brt_ES_N, dat_hist, n.trees=brt_ES_N$gbm.call$best.trees, type="response"))
    dat_hist$brt_ES <- presx * abundx
    
    presx <- predict(brt_ES_P, dat_fcast, n.trees=brt_ES_P$gbm.call$best.trees, type="response")
    abundx <- exp(predict(brt_ES_N, dat_fcast, n.trees=brt_ES_N$gbm.call$best.trees, type="response"))
    dat_fcast$brt_ES <- presx * abundx
  }
}

if ("EST" %in% covs) {
  print("Fitting BRT-EST")
  if (type == "delta") {
    brt_EST_P <- gbm.step(data=dat_hist,
                        gbm.x = c(gbm.x.env, "lat", "lon", "year"),
                        gbm.y = 'pres',
                        family = "bernoulli",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,#made 0.05 for anchovy  temp only EMs
                        plot.main=FALSE, verbose = FALSE)
    
    brt_EST_N <- gbm.step(data=dat_hist[dat_hist$abundance>0,], 
                        gbm.x = c(gbm.x.env, "lat", "lon", "year"),
                        gbm.y = 'log_abundance',
                        family = "gaussian",
                        tree.complexity = 3, learning.rate = 0.01, bag.fraction = 0.6,#made 0.05 for anchovy  temp only EMs
                        plot.main=FALSE, verbose = FALSE)
    
    presx <- predict(brt_EST_P, dat_hist, n.trees=brt_EST_P$gbm.call$best.trees, type="response")
    abundx <- exp(predict(brt_EST_N, dat_hist, n.trees=brt_EST_N$gbm.call$best.trees, type="response"))
    dat_hist$brt_EST <- presx * abundx
    
    presx <- predict(brt_EST_P, dat_fcast, n.trees=brt_EST_P$gbm.call$best.trees, type="response")
    abundx <- exp(predict(brt_EST_N, dat_fcast, n.trees=brt_EST_N$gbm.call$best.trees, type="response"))
    dat_fcast$brt_EST <- presx * abundx
  }
}

if (type == "delta") {
  gbm.plot(brt_E_P, write.title=F, main="brt_E_P", plot.layout = c(2,2))
  #gbm.plot(brt_S_P, write.title=F, main="brt_S_P", plot.layout = c(1,2))
  gbm.plot(brt_ES_P, write.title=F, main="brt_ES_P", plot.layout = c(2,3))
  gbm.plot(brt_EST_P, write.title=F, main="brt_EST_P", plot.layout = c(2,3))
  
  gbm.plot(brt_E_N, write.title=F, main="brt_E_N", plot.layout = c(2,2))
  #gbm.plot(brt_S_N, write.title=F, main="brt_S_N", plot.layout = c(1,2))
  gbm.plot(brt_ES_N, write.title=F, main="brt_ES_N", plot.layout = c(2,3))
  gbm.plot(brt_EST_N, write.title=F, main="brt_EST_N", plot.layout = c(2,3))
}
