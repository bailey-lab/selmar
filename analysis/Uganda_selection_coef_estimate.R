selection_coef_estimate <- function(df, PfK13_stanmod, c469y_stanmod, a675v_stanmod, c469f_stanmod, p441l_stanmod, numofobs){
  # Posterior of the model
  df$predict_stan <- NA
  df$predict_stan[is.finite(df$lrsmed) & df$Locus == "C469Y" & df$nobs > numofobs] <-  rstanarm::posterior_predict(c469y_stanmod, type = "response") %>% colMeans
  df$predict_stan[is.finite(df$lrsmed) & df$Locus == "A675V" & df$nobs > numofobs] <-  rstanarm::posterior_predict(a675v_stanmod, type = "response") %>% colMeans
  df$predict_stan[is.finite(df$lrsmed) & df$Locus == "C469F" & df$nobs > numofobs] <-  rstanarm::posterior_predict(c469f_stanmod, type = "response") %>% colMeans
  df$predict_stan[is.finite(df$lrsmed) & df$Locus == "P441L" & df$nobs > numofobs] <-  rstanarm::posterior_predict(p441l_stanmod, type = "response") %>% colMeans
  df$predict_stan[is.finite(df$lrsmed) & df$Locus == "K13" & df$nobs > numofobs] <-  rstanarm::posterior_predict(PfK13_stanmod, type = "response") %>% colMeans

  #Selection coefficient for individual mutations
  df$conf.low <- NA
  df$conf.high <- NA
  c469y_conf <- rstanarm::posterior_predict(c469y_stanmod, type = "response") %>% apply(2, quantile, prob = (c(0.025, 0.975)))
  c469f_conf <- rstanarm::posterior_predict(c469f_stanmod, type = "response") %>% apply(2, quantile, prob = (c(0.025, 0.975)))
  a675v_conf <- rstanarm::posterior_predict(a675v_stanmod, type = "response") %>% apply(2, quantile, prob = (c(0.025, 0.975)))
  p441l_conf <- rstanarm::posterior_predict(p441l_stanmod, type = "response") %>% apply(2, quantile, prob = (c(0.025, 0.975)))
  PfK13_conf <- rstanarm::posterior_predict(PfK13_stanmod, type = "response") %>% apply(2, quantile, prob = (c(0.025, 0.975)))

  # here just add predictions  for the c469Y rows
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "C469Y" & df$nobs > numofobs] <- c469y_conf[1,]
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "A675V" & df$nobs > numofobs] <- a675v_conf[1,]
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "C469F" & df$nobs > numofobs] <- c469f_conf[1,]
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "P441L" & df$nobs > numofobs] <- p441l_conf[1,]
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "K13" & df$nobs > numofobs] <- PfK13_conf[1,]

  df$conf.high[is.finite(df$lrsmed) & df$Locus == "C469Y" & df$nobs > numofobs] <- c469y_conf[2,]
  df$conf.high[is.finite(df$lrsmed) & df$Locus == "A675V" & df$nobs > numofobs] <- a675v_conf[2,]
  df$conf.high[is.finite(df$lrsmed) & df$Locus == "C469F" & df$nobs > numofobs] <- c469f_conf[2,]
  df$conf.high[is.finite(df$lrsmed) & df$Locus == "P441L" & df$nobs > numofobs] <- p441l_conf[2,]
  df$conf.high[is.finite(df$lrsmed) & df$Locus == "K13" & df$nobs > numofobs] <- PfK13_conf[2,]

  #Prevalence scale
  df$predict_stan_prev <- NA
  df$conf.high_prev <- NA
  df$conf.low_prev <- NA

  df$predict_stan_prev <- exp(df$predict_stan)*100/(1+exp(df$predict_stan))
  df$conf.high_prev <- exp(df$conf.high)*100/(1+exp(df$conf.high))
  df$conf.low_prev <- exp(df$conf.low)*100/(1+exp(df$conf.low))
  return(df)
  }
