### Location^3 Estimation Model Code
## Function to run GLMs on data from OM
## FUnction returns only the fitted and predicted values

## Code by Lewis Barnett and James Smith


  #'dat' are all the observations from the OM
  #'year_fcast' is the year used to split historical (fitting) and forecast (testing) data
  #'type' is either 'tweedie' or 'delta'
  #'covs' is a vector indicating model covariates, where:
  #         - E is enviro only
  #         - Sr is spatiotemporal random field only
  #         - ESt is enviro + spatial trend random field
  #         - ESr is enviro + spatiotemporal random field
  #         - c("E", "Sr", "ESt", "ESr") does them all
  
  ## Comments by LB:
  # note that there is a much cleaner way of iterating through models 
  # but for now simplifying for clarity and consistency with GAM fit formatting above
  # see ?sdmTMB for options



  # # project to UTM with spTransform [probably only necessary when measuring COG? otherwise this changes the original data file halfway through the analysis]
  # dat_ll <- dat
  # coordinates(dat_ll) <- c("lon", "lat")
  # proj4string(dat_ll) <- CRS("+proj=longlat +datum=WGS84")
  # dat_utm <- spTransform(dat_ll, 
  #                        CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # dat = as.data.frame(dat_utm) # convert back from sp object to data frame
  

  # need to remove future years from test data set via weights
  weights <- rep(0,nrow(dat))
  weights[which(dat$year <= year_fcast)] <- 1  #keep historical period
  
  # subset of data with only presences
  dat_pres <- dat[dat$pres == 1,]
  
  # make SPDE (stochastic partial differential equation that these INLA-based methods rely on)
  spde <- try(make_mesh(data = dat, xy_cols = c("lon","lat"), n_knots = 200), silent=TRUE) # increase knots to at least ~250 for WC data
  spde_pos <- try(make_mesh(data = dat_pres, xy_cols = c("lon","lat"), n_knots = 200), silent=TRUE)
  # plot(spde) # can vary number of knots to modify the mesh until you get what you want
  
  tw_formula <- formula(paste("abundance ~ 0 +", env_formula))
  delta1_formula <- formula(paste("pres ~ 0 +", env_formula))
  delta2_formula <- formula(paste("log_abundance ~ 0 +", env_formula_deltaN))  #JS: I find logging abundance leads to fewer fitting issues

  
  if ("E" %in% covs) {  #Enviro only; no spatiotemporal random fields (but simple space?)
    print("Fitting GLM-E")
    if (type == "tweedie") {
      glm_E <- try(sdmTMB(  #'glmm3' in original script
        formula = tw_formula,
        time_varying = NULL,
        spde = spde,
        time = "year",
        family = tweedie(link = "log"),
        data = dat,
        anisotropy = TRUE,
        ar1_fields = FALSE,
        weights = weights,
        spatial_only = TRUE,
        spatial_trend = FALSE,
        quadratic_roots = TRUE,
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      ))
      
      P_glm_E <- predict(glm_E)
      dat_hist$glm_E <- exp(P_glm_E$est[P_glm_E$year <= year_fcast])
      dat_fcast$glm_E <- exp(P_glm_E$est[P_glm_E$year > year_fcast])
    }
    if (type == "delta") {
      glm_E_P <- try(sdmTMB(  #presence-absence part
        formula = delta1_formula,
        time_varying = NULL,
        spde = spde,
        time = "year",
        family = binomial(link = "logit"),
        data = dat,
        anisotropy = TRUE,
        ar1_fields = FALSE,
        weights = weights,
        spatial_only = TRUE,
        spatial_trend = FALSE,
        quadratic_roots = TRUE,
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      ))
      glm_E_N <- try(sdmTMB(  #abundance when non-zero part
        formula = delta2_formula,
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
      ))
      
      P_glm_E_P <- predict(glm_E_P)  #binomial part, link scale
      P_glm_E_P$est_prob <- exp(P_glm_E_P$est)/(1 + exp(P_glm_E_P$est))  #as probabilities
      P_glm_E_N <- predict(glm_E_N, newdata = dat)  #abundance part
      P_glm_E_P$est_delta <- P_glm_E_P$est_prob * exp(P_glm_E_N$est)  #'exp' here if using log_abundance as response
      
      #TEMPORARY
      #plot above
      # plot(aggregate(est_delta ~ year,data=P_glm_E_P, FUN="sum"), type='b')
      #plot get_index()
      # P_glm_E_N <- predict(glm_E_N, newdata = dat,return_tmb_object = TRUE )  #abundance part
      # test <- get_index(P_glm_E_N, bias_correct = TRUE)
      # ggplot(test, aes(year, est)) + geom_line() +
      #   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
      #   xlab('Year') + ylab('Biomass estimate (metric tonnes)')
      
      dat_hist$glm_E <- P_glm_E_P$est_delta[P_glm_E_P$year <= year_fcast]
      dat_fcast$glm_E <- P_glm_E_P$est_delta[P_glm_E_P$year > year_fcast]
    }
  }
  
  if ("Sr" %in% covs) {  #only spatiotemporal random fields as AR1
    print("Fitting GLM-Sr")
    if (type == "tweedie") {
      glm_Sr <- try(sdmTMB(  #'glmm1' in original script
        formula = abundance ~ 1,
        time_varying = NULL,
        spde = spde,
        time = "year",
        family = tweedie(link = "log"),
        data = dat,
        anisotropy = TRUE,
        ar1_fields = TRUE,  #false for spatiotemporal random fields as IID, true for ar1
        weights = weights,
        spatial_only = FALSE,
        spatial_trend = FALSE,  #true for spatiotemporal random fields as IID, false for ar1
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      ))
      
      P_glm_Sr <- predict(glm_Sr)
      dat_hist$glm_Sr <- exp(P_glm_Sr$est[P_glm_Sr$year <= year_fcast])
      dat_fcast$glm_Sr <- exp(P_glm_Sr$est[P_glm_Sr$year > year_fcast])
    }
    if (type == "delta") {
      start <- Sys.time()
      print("glm_Sr_P")
      glm_Sr_P <- try(sdmTMB(  #presence-absence part
        formula = pres ~ 1,
        time_varying = NULL,
        spde = spde,
        time = "year",
        family = binomial(link = "logit"),
        data = dat,
        anisotropy = TRUE,
        ar1_fields = TRUE,
        weights = weights,
        spatial_only = FALSE,
        spatial_trend = FALSE,
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      ))
      end <- Sys.time(); round(end-start, 2)
      
      start <- Sys.time()
      print("glm_Sr_N")
      glm_Sr_N <- try(sdmTMB(  #non-zero abundance part
        formula = log_abundance ~ 1,
        time_varying = NULL,
        spde = spde_pos,
        time = "year",
        family = gaussian(link = "identity"),
        data = dat_pres,
        anisotropy = TRUE,
        ar1_fields = TRUE,
        weights = weights[which(dat_pres$pres == 1)],
        spatial_only = FALSE,
        spatial_trend = FALSE,
        control = sdmTMBcontrol(step.min = 0.01, step.max = 1)
      ))
      end <- Sys.time(); round(end-start, 2)

      P_glm_Sr_P <- predict(glm_Sr_P)  #binomial part, link scale
      P_glm_Sr_P$est_prob <- exp(P_glm_Sr_P$est)/(1 + exp(P_glm_Sr_P$est))  #as probabilities
      P_glm_Sr_N <- predict(glm_Sr_N, newdata = dat)  #abundance part
      P_glm_Sr_P$est_delta <- P_glm_Sr_P$est_prob * exp(P_glm_Sr_N$est)  #'exp' here if using log_abundance as response
      dat_hist$glm_Sr <- P_glm_Sr_P$est_delta[P_glm_Sr_P$year <= year_fcast]
      dat_fcast$glm_Sr <- P_glm_Sr_P$est_delta[P_glm_Sr_P$year > year_fcast]
    }
  }
  
  if ("ESt" %in% covs) {  # enviro and spatial trend random field, but no spatiotemporal random fields
    print("Fitting GLM-ESt")
    if (type == "tweedie") {
      glm_ESt <- try(sdmTMB(  #'glm4' in original script
        formula = tw_formula,
        time_varying = NULL,
        spde = spde,
        time = "year",
        family = tweedie(link = "log"),
        data = dat,
        anisotropy = TRUE,
        ar1_fields = FALSE,
        weights = weights,
        spatial_only = TRUE,
        spatial_trend = TRUE,
        quadratic_roots = TRUE
      ))
      
      P_glm_ESt <- predict(glm_ESt)
      dat_hist$glm_ESt <- exp(P_glm_ESt$est[P_glm_ESt$year <= year_fcast])
      dat_fcast$glm_ESt <- exp(P_glm_ESt$est[P_glm_ESt$year > year_fcast])
    }
    if (type == "delta") {
      glm_ESt_P <- try(sdmTMB(  #presence-absence part
        formula = delta1_formula,
        time_varying = NULL,
        spde = spde,
        time = "year",
        family = binomial(link = "logit"),
        data = dat,
        anisotropy = TRUE,
        ar1_fields = FALSE,
        weights = weights,
        spatial_only = TRUE,
        spatial_trend = TRUE,
        quadratic_roots = TRUE
      ))
      glm_ESt_N <- try(sdmTMB(  #non-zero abundance part
        formula = delta2_formula,
        time_varying = NULL,
        spde = spde_pos,
        time = "year",
        family = gaussian(link = "identity"),
        data = dat_pres,
        anisotropy = TRUE,
        ar1_fields = FALSE,
        weights = weights[which(dat_pres$pres == 1)],
        spatial_only = TRUE,
        spatial_trend = TRUE,
        quadratic_roots = TRUE
      ))

      P_glm_ESt_P <- predict(glm_ESt_P)  #binomial part, link scale
      P_glm_ESt_P$est_prob <- exp(P_glm_ESt_P$est)/(1 + exp(P_glm_ESt_P$est))  #as probabilities
      P_glm_ESt_N <- predict(glm_ESt_N, newdata = dat)  #abundance part
      P_glm_ESt_P$est_delta <- P_glm_ESt_P$est_prob * exp(P_glm_ESt_N$est)  #'exp' here if using log_abundance as response
      dat_hist$glm_ESt <- P_glm_ESt_P$est_delta[P_glm_ESt_P$year <= year_fcast]
      dat_fcast$glm_ESt <- P_glm_ESt_P$est_delta[P_glm_ESt_P$year > year_fcast]
    }
  }
  
  if ("ESr" %in% covs) {  # enviro and spatiotemporal random fields as AR1
    print("Fitting GLM-ESr")
    if (type == "tweedie") {
      glm_ESr <- try(sdmTMB(  #'glmm5' in original script
        formula = tw_formula,
        time_varying = NULL,
        spde = spde,
        time = "year",
        family = tweedie(link = "log"),
        data = dat,
        anisotropy = TRUE,
        ar1_fields = TRUE,
        weights = weights,
        spatial_only = FALSE,
        spatial_trend = FALSE,
        quadratic_roots = TRUE
      ))

      P_glm_ESr <- predict(glm_ESr)
      dat_hist$glm_ESr <- exp(P_glm_ESr$est[P_glm_ESr$year <= year_fcast])
      dat_fcast$glm_ESr <- exp(P_glm_ESr$est[P_glm_ESr$year > year_fcast])
    }
    if (type == "delta") {
      print("Fitting GLM-E_P")
      glm_ESr_P <- try(sdmTMB(  #presence-absence part
        formula = delta1_formula,
        time_varying = NULL,
        spde = spde,
        time = "year",
        family = binomial(link = "logit"),
        data = dat,
        anisotropy = TRUE,
        ar1_fields = TRUE,
        weights = weights,
        spatial_only = FALSE,
        spatial_trend = FALSE,
        quadratic_roots = TRUE
      ))
      print("Fitting GLM-E_N")
      glm_ESr_N <- try(sdmTMB(  #non-zero abundance part
        formula = delta2_formula,
        time_varying = NULL,
        spde = spde_pos,
        time = "year",
        family = gaussian(link = "identity"),
        data = dat_pres,
        anisotropy = TRUE,
        ar1_fields = TRUE,
        weights = weights[which(dat_pres$pres == 1)],
        spatial_only = FALSE,
        spatial_trend = FALSE,
        quadratic_roots = TRUE
      ))

      P_glm_ESr_P <- predict(glm_ESr_P)  #binomial part, link scale
      P_glm_ESr_P$est_prob <- exp(P_glm_ESr_P$est)/(1 + exp(P_glm_ESr_P$est))  #as probabilities
      P_glm_ESr_N <- predict(glm_ESr_N, newdata = dat)  #abundance part
      P_glm_ESr_P$est_delta <- P_glm_ESr_P$est_prob * exp(P_glm_ESr_N$est)  #'exp' here if using log_abundance as response
      dat_hist$glm_ESr <- P_glm_ESr_P$est_delta[P_glm_ESr_P$year <= year_fcast]
      dat_fcast$glm_ESr <- P_glm_ESr_P$est_delta[P_glm_ESr_P$year > year_fcast]
    }
  }
  
  
  # # plot 'partial effect' for base model ***for some reason this isn't working - still says mismatch of years
  # predict_glmm <- function(model, new_dat) {
  #   dummy = data.frame(year = unique(dat$year),
  #                      lat=dat$lat[1],
  #                      lon=unique(dat$lon[1]),
  #                      temp_s = dat$temp_s[1],
  #                      chla_s = dat$chla_s[1],
  #                      mld_s = dat$chla_s[1])
  #   pred = predict(model,
  #                  newdata=rbind(new_dat[,c("year","lat","lon","temp_s","chla_s", "mld_s")],
  #                                dummy), xy_cols = c("lon", "lat"))
  #   
  #   # drop dummy data
  #   pred = pred[-c(seq(nrow(pred)-nrow(dummy)+1,nrow(pred))),]
  #   pred$abundance = dat_x$abundance # add true abundance
  #   
  #   return(pred)
  # }
  # 
  # newdat_temp <- data.frame(temp_raw=dat$temp, temp_s=dat$temp_s,
  #                           chla_s=median(dat$chla_s),
  #                           mld_s=median(dat$mld_s),
  #                           lat=median(dat$lat),
  #                           lon=median(dat$lon),
  #                           year=median(dat$year))
  # newdat_temp <- newdat_temp[order(newdat_temp$temp_s),]
  # 
  # P1 <- predict_glmm(glm_E_P, newdat_temp)
  # P1 <- exp(P1$est)/(1 + exp(P1$est))
  # P2 <- predict_glmm(glm_E_N, newdat_temp)
  # P2 <- exp(P2$est)
  # 
  # plot(newdat_temp$temp_raw, P1, type="l", main="glm_E_P - Temp")
  # plot(newdat_temp$temp_raw, P2, type="l", main="glm_E_N - Temp")
  


  
