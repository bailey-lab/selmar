GLM_model_mutant <- function(df, locus){
  stanmod <- rstanarm::stan_glmer(
    lrsmed ~ adj_year + (1 + adj_year | District),
    data = df %>%
      filter(Locus == locus, nobs > numofobs),
    weights = wt_med,
    adapt_delta = 0.9975
  )

  vals <- summary(stanmod$stanfit)$summary[,1]
  vals <- vals[grep("^b", names(vals))]
  ranef_stanmod <- data.frame("Intercepts" = vals[seq(1, length(vals), 2)],
                                    "Slopes" = vals[seq(2, length(vals), 2)])
  ranef_stanmod$mutation <- rep(locus, dim(ranef_stanmod)[1])

  output <- list(glmer = stanmod, ranef = ranef_stanmod)
  return(output)
}

selection_ceof_estimate <- function(df, c580y_stanmod, r539t_stanmod, k13_stanmod){
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
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "C580Y" & df$nobs > numofobs] <- c580y_conf[1,]
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "R539T" & df$nobs > numofobs] <- r539t_conf[1,]
  df$conf.low[is.finite(df$lrsmed) & df$Locus == "K13" & df$nobs > numofobs] <- K13_conf[1,]
  #obtain higher critical interval
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

preprocessing_selection_coef <- function(c580y_stanmod, r539t_stanmod, K13_stanmod){
  # Effect sizes
  stan_eff <- as.data.frame(t(data.frame(c580y_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
    mutate(method = "Bayesian (Intercept + Slope)") %>%
    setNames(c("conf.low", "estimate", "conf.high", "method"))

  effect_size <- stan_eff %>% mutate(mutation = "C580Y")

  stan_eff <- as.data.frame(t(data.frame(r539t_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
    mutate(method = "Bayesian (Intercept + Slope)") %>%
    setNames(c("conf.low", "estimate", "conf.high", "method"))

  effect_size <- rbind(effect_size, stan_eff %>% mutate(mutation = "R539T"))

  stan_eff <- as.data.frame(t(data.frame(K13_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
    mutate(method = "Bayesian (Intercept + Slope)") %>%
    setNames(c("conf.low", "estimate", "conf.high", "method"))

  effect_size <- rbind(effect_size, stan_eff %>% mutate(mutation = "K13"))

  #confidence interval per site and mutation
  stan_eff_c580y <- as.data.frame(c580y_stanmod$stan_summary[grep("year District", rownames(c580y_stanmod$stan_summary)),]) %>%
    select(c("2.5%","97.5%")) %>%
    setNames(c("conf.low", "conf.high")) %>%
    slice(1:(n()-1))
  stan_eff_r539t <- as.data.frame(r539t_stanmod$stan_summary[grep("year District", rownames(r539t_stanmod$stan_summary)),]) %>%
    select(c("2.5%","97.5%")) %>%
    setNames(c("conf.low", "conf.high")) %>%
    slice(1:(n()-1))
  stan_eff_K13 <- as.data.frame(K13_stanmod$stan_summary[grep("year District", rownames(K13_stanmod$stan_summary)),]) %>%
    select(c("2.5%","97.5%")) %>%
    setNames(c("conf.low", "conf.high")) %>%
    slice(1:(n()-1))

  coef_c580y <- ranef(c580y_stanmod)$District %>%
    mutate(District = rownames(ranef(c580y_stanmod)$District),
           mutation = "C580Y",
           estimate = ranef(c580y_stanmod)$District$adj_year,
           conf.low = stan_eff_c580y$conf.low,
           conf.high = stan_eff_c580y$conf.high) %>%
    select(conf.low, estimate, conf.high, mutation, District) %>%
    remove_rownames()

  coef_r539t <- ranef(r539t_stanmod)$District %>%
    mutate(District = rownames(ranef(r539t_stanmod)$District),
           mutation = "R539T",
           estimate = ranef(r539t_stanmod)$District$adj_year,
           conf.low = stan_eff_r539t$conf.low,
           conf.high = stan_eff_r539t$conf.high) %>%
    select(conf.low, estimate, conf.high, mutation, District) %>%
    remove_rownames()

  coef_K13 <- ranef(K13_stanmod)$District %>%
    mutate(District = rownames(ranef(K13_stanmod)$District),
           mutation = "K13",
           estimate = ranef(K13_stanmod)$District$adj_year,
           conf.low = stan_eff_K13$conf.low,
           conf.high = stan_eff_K13$conf.high) %>%
    select(conf.low, estimate, conf.high, mutation, District) %>%
    remove_rownames()

  selec_coef <- rbind(coef_c580y, coef_r539t, coef_K13)

  # plot the effect sizes for the mutation
  effect_size_bay <- effect_size %>%
    filter(method == "Bayesian (Intercept + Slope)") %>%
    mutate(District = "Southeast Asia") %>%
    select(conf.low, estimate, conf.high, mutation, District)

  selec_coef$estimate[selec_coef$mutation == "C580Y"] <- selec_coef$estimate[selec_coef$mutation == "C580Y"] + effect_size_bay$estimate[effect_size_bay$mutation == "C580Y"]
  selec_coef$estimate[selec_coef$mutation == "R539T"] <- selec_coef$estimate[selec_coef$mutation == "R539T"] + effect_size_bay$estimate[effect_size_bay$mutation == "R539T"]
  selec_coef$estimate[selec_coef$mutation == "K13"] <- selec_coef$estimate[selec_coef$mutation == "K13"] + effect_size_bay$estimate[effect_size_bay$mutation == "K13"]

  selec_coef$conf.low[selec_coef$mutation == "C580Y"] <- selec_coef$conf.low[selec_coef$mutation == "C580Y"] + effect_size_bay$conf.low[effect_size_bay$mutation == "C580Y"]
  selec_coef$conf.low[selec_coef$mutation == "R539T"] <- selec_coef$conf.low[selec_coef$mutation == "R539T"] + effect_size_bay$conf.low[effect_size_bay$mutation == "R539T"]
  selec_coef$conf.low[selec_coef$mutation == "K13"] <- selec_coef$conf.low[selec_coef$mutation == "K13"] + effect_size_bay$conf.low[effect_size_bay$mutation == "K13"]

  selec_coef$conf.high[selec_coef$mutation == "C580Y"] <- selec_coef$conf.high[selec_coef$mutation == "C580Y"] + effect_size_bay$conf.high[effect_size_bay$mutation == "C580Y"]
  selec_coef$conf.high[selec_coef$mutation == "R539T"] <- selec_coef$conf.high[selec_coef$mutation == "R539T"] + effect_size_bay$conf.high[effect_size_bay$mutation == "R539T"]
  selec_coef$conf.high[selec_coef$mutation == "K13"] <- selec_coef$conf.high[selec_coef$mutation == "K13"] + effect_size_bay$conf.high[effect_size_bay$mutation == "K13"]

  effect_size_bay$s <- 18
  effect_size_bay$sizes <- 1.1
  selec_coef$s <- 19
  selec_coef$sizes <- 0.5

  stan_eff_bay_overall <- effect_size_bay %>%
    mutate(method = "Baysian") %>%
    relocate(District, conf.low, estimate, conf.high, method, mutation, s, sizes)
  rownames(stan_eff_bay_overall) <- NULL

  # selec_coef_bay <- selec_coef %>% mutate(method = "Baysian") %>% relocate(District, conf.low, estimate, conf.high, method, mutation, s, sizes)

  effect_size_bay <- rbind(selec_coef, effect_size_bay)

  effect_size_bay <- effect_size_bay %>%
    mutate(mutation = replace(mutation, mutation == "C580Y", "580Y"),
           mutation = replace(mutation, mutation == "R539T", "539T"))
  return(effect_size_bay)
}

analysis_selection_coef <- function(effect_size_bay){
  selec_ceof_stats <- data.frame(matrix(ncol = 5, nrow = 3))
  min_580Y <- min(effect_size_bay %>% filter(mutation == "580Y" & District != "Southeast Asia") %>% pull(estimate))
  max_580Y <- max(effect_size_bay %>% filter(mutation == "580Y" & District != "Southeast Asia") %>% pull(estimate))
  selec_ceof_stats[1, ] <- cbind("580Y",
                                 round(min_580Y, 3),
                                 effect_size_bay[which(effect_size_bay$estimate == min_580Y),][5],
                                 round(max_580Y, 3),
                                 effect_size_bay[which(effect_size_bay$estimate == max_580Y),][5])

  min_539T <- min(effect_size_bay %>% filter(mutation == "539T" & District != "Southeast Asia") %>% pull(estimate))
  max_539T <- max(effect_size_bay %>% filter(mutation == "539T" & District != "Southeast Asia") %>% pull(estimate))
  selec_ceof_stats[2, ] <- cbind("539T",
                                 round(min_539T, 3),
                                 effect_size_bay[which(effect_size_bay$estimate == min_539T),][5],
                                 round(max_539T, 3),
                                 effect_size_bay[which(effect_size_bay$estimate == max_539T),][5])

  min_K13 <- min(effect_size_bay %>% filter(mutation == "K13" & District != "Southeast Asia") %>% pull(estimate))
  max_K13 <- max(effect_size_bay %>% filter(mutation == "K13" & District != "Southeast Asia") %>% pull(estimate))
  selec_ceof_stats[3, ] <- cbind("K13",
                                 round(min_K13, 3),
                                 effect_size_bay[which(effect_size_bay$estimate == min_K13),][5],
                                 round(max_K13, 3),
                                 effect_size_bay[which(effect_size_bay$estimate == max_K13),][5])
  return(selec_ceof_stats)
}
