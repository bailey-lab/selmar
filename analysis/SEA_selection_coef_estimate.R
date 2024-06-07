#Estimating selection coefficients from GLM outputs for each mutation
#Input: GLM models for each mutation, and dataframe
#Output: Dataframe with selection coefficients and critical intervals
selection_coef_estimate <- function(df, c580y_stanmod, r539t_stanmod, k13_stanmod, mutations){

  #  get model predictions for plotting
  df$predict_stan <- NA
  df$predict_stan[is.finite(df$lrsmed) & df$Locus == "C580Y" & df$nobs > numofobs] <-  rstanarm::posterior_predict(c580y_stanmod, type = "response") %>% colMeans
  df$predict_stan[is.finite(df$lrsmed) & df$Locus == "R539T" & df$nobs > numofobs] <-  rstanarm::posterior_predict(r539t_stanmod, type = "response") %>% colMeans
  df$predict_stan[is.finite(df$lrsmed) & df$Locus == "K13" & df$nobs > numofobs] <-  rstanarm::posterior_predict(K13_stanmod, type = "response") %>% colMeans

  #confidence interval
  # get model predictions for plotting
  c580y_conf <- rstanarm::posterior_predict(c580y_stanmod, type = "response") %>% apply(2, quantile, prob = (c(0.025, 0.975)))
  r539t_conf <- rstanarm::posterior_predict(r539t_stanmod, type = "response") %>% apply(2, quantile, prob = (c(0.025, 0.975)))
  K13_conf <- rstanarm::posterior_predict(K13_stanmod, type = "response") %>% apply(2, quantile, prob = (c(0.025, 0.975)))

  #obtain lower critical interval
  df$conf.low <- NA
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "C580Y" & df$nobs > numofobs] <- c580y_conf[1,]
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "R539T" & df$nobs > numofobs] <- r539t_conf[1,]
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "K13" & df$nobs > numofobs] <- K13_conf[1,]
  #obtain higher critical interval
  df$conf.high <- NA
  df$conf.high[is.finite(df$lrsmed) & df$Locus == "C580Y" & df$nobs > numofobs] <- c580y_conf[2,]
  df$conf.high[is.finite(df$lrsmed) & df$Locus == "R539T" & df$nobs > numofobs] <- r539t_conf[2,]
  df$conf.high[is.finite(df$lrsmed) & df$Locus == "K13" & df$nobs > numofobs] <- K13_conf[2,]

  # Plot the predictions in each site with only Bayesian
  #convert to prevalence relative to WT scale
  df$predict_stan_prev <- exp(df$predict_stan)/(1+exp(df$predict_stan))
  df$conf.high_prev <- exp(df$conf.high)/(1+exp(df$conf.high))
  df$conf.low_prev <- exp(df$conf.low)/(1+exp(df$conf.low))

  return(df)
}
