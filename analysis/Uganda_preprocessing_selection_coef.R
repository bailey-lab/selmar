preprocessing_selection_coef <- function(PfK13_stanmod, c469y_stanmod, a675v_stanmod, c469f_stanmod, p441l_stanmod){
  stan_eff <- as.data.frame(t(data.frame(c469y_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
    mutate(method = "Bayesian (Intercept + Slope)") %>%
    setNames(c("conf.low", "estimate", "conf.high", "method"))

  effect_size <- stan_eff %>% mutate(mutation = "C469Y")

  stan_eff <- as.data.frame(t(data.frame(a675v_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
    mutate(method = "Bayesian (Intercept + Slope)") %>%
    setNames(c("conf.low", "estimate", "conf.high", "method"))

  pl <- stan_eff %>% mutate(mutation = "A675V")
  effect_size <- rbind(effect_size, pl)

  stan_eff <- as.data.frame(t(data.frame(c469f_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
    mutate(method = "Bayesian (Intercept + Slope)") %>%
    setNames(c("conf.low", "estimate", "conf.high", "method"))

  pl <- stan_eff %>% mutate(mutation = "C469F")
  effect_size <- rbind(effect_size, pl)

  stan_eff <- as.data.frame(t(data.frame(p441l_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
    mutate(method = "Bayesian (Intercept + Slope)") %>%
    setNames(c("conf.low", "estimate", "conf.high", "method"))

  pl <- stan_eff %>% mutate(mutation = "P441L")
  effect_size <- rbind(effect_size, pl)

  stan_eff <- as.data.frame(t(data.frame(PfK13_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
    mutate(method = "Bayesian (Intercept + Slope)") %>%
    setNames(c("conf.low", "estimate", "conf.high", "method"))

  pl <- stan_eff %>% mutate(mutation = "K13")
  effect_size <- rbind(effect_size, pl)


  #confidence interval per site and mutation
  stan_eff_a675v <- as.data.frame(a675v_stanmod$stan_summary[grep("year District", rownames(a675v_stanmod$stan_summary)),]) %>% select(c("2.5%","97.5%")) %>% setNames(c("conf.low", "conf.high")) %>% slice(1:(n()-1))
  stan_eff_c469y <- as.data.frame(c469y_stanmod$stan_summary[grep("year District", rownames(c469y_stanmod$stan_summary)),]) %>% select(c("2.5%","97.5%")) %>% setNames(c("conf.low", "conf.high")) %>% slice(1:(n()-1))
  stan_eff_c469f <- as.data.frame(c469f_stanmod$stan_summary[grep("year District", rownames(c469f_stanmod$stan_summary)),]) %>% select(c("2.5%","97.5%")) %>% setNames(c("conf.low", "conf.high")) %>% slice(1:(n()-1))
  stan_eff_p441l <- as.data.frame(p441l_stanmod$stan_summary[grep("year District", rownames(p441l_stanmod$stan_summary)),]) %>% select(c("2.5%","97.5%")) %>% setNames(c("conf.low", "conf.high")) %>% slice(1:(n()-1))
  stan_eff_PfK13 <- as.data.frame(PfK13_stanmod$stan_summary[grep("year District", rownames(PfK13_stanmod$stan_summary)),]) %>% select(c("2.5%","97.5%")) %>% setNames(c("conf.low", "conf.high")) %>% slice(1:(n()-1))

  coef_a675v <- ranef(a675v_stanmod)$District %>%
    mutate(District = rownames(ranef(a675v_stanmod)$District),
           mutation = "A675V",
           estimate = ranef(a675v_stanmod)$District$adj_year,
           conf.low = stan_eff_a675v$conf.low,
           conf.high = stan_eff_a675v$conf.high) %>%
    select(conf.low, estimate, conf.high, mutation, District)
  coef_c469y <- ranef(c469y_stanmod)$District %>%
    mutate(District = rownames(ranef(c469y_stanmod)$District),
           mutation = "C469Y",
           estimate = ranef(c469y_stanmod)$District$adj_year,
           conf.low = stan_eff_c469y$conf.low,
           conf.high = stan_eff_c469y$conf.high) %>%
    select(conf.low, estimate, conf.high, mutation, District)
  coef_c469f <- ranef(c469f_stanmod)$District %>%
    mutate(District = rownames(ranef(c469f_stanmod)$District),
           mutation = "C469F",
           estimate = ranef(c469f_stanmod)$District$adj_year,
           conf.low = stan_eff_c469f$conf.low,
           conf.high = stan_eff_c469f$conf.high) %>%
    select(conf.low, estimate, conf.high, mutation, District)
  coef_p441l <- ranef(p441l_stanmod)$District %>%
    mutate(District = rownames(ranef(p441l_stanmod)$District),
           mutation = "P441L",
           estimate = ranef(p441l_stanmod)$District$adj_year,
           conf.low = stan_eff_p441l$conf.low,
           conf.high = stan_eff_p441l$conf.high) %>%
    select(conf.low, estimate, conf.high, mutation, District)
  coef_PfK13 <- ranef(PfK13_stanmod)$District %>%
    mutate(District = rownames(ranef(PfK13_stanmod)$District),
           mutation = "K13",
           estimate = ranef(PfK13_stanmod)$District$adj_year,
           conf.low = stan_eff_PfK13$conf.low,
           conf.high = stan_eff_PfK13$conf.high) %>%
    select(conf.low, estimate, conf.high, mutation, District)

  selec_coef <- rbind(coef_a675v, coef_c469y, coef_c469f, coef_p441l, coef_PfK13)
  rownames(selec_coef) <- NULL
  conf_val <- as.data.frame(selec_coef$estimate[1]) %>%
    apply(2, quantile, prob = (c(0.025, 0.975)))

  # plot the effect sizes for the mutation
  effect_size_bay <- effect_size %>%
    filter(method == "Bayesian (Intercept + Slope)") %>%
    mutate(District = "Uganda") %>%
    select(conf.low, estimate, conf.high, mutation, District)

  selec_coef$estimate[selec_coef$mutation == "A675V"] <- selec_coef$estimate[selec_coef$mutation == "A675V"] + effect_size_bay$estimate[effect_size_bay$mutation == "A675V"]
  selec_coef$estimate[selec_coef$mutation == "C469Y"] <- selec_coef$estimate[selec_coef$mutation == "C469Y"] + effect_size_bay$estimate[effect_size_bay$mutation == "C469Y"]
  selec_coef$estimate[selec_coef$mutation == "C469F"] <- selec_coef$estimate[selec_coef$mutation == "C469F"] + effect_size_bay$estimate[effect_size_bay$mutation == "C469F"]
  selec_coef$estimate[selec_coef$mutation == "P441L"] <- selec_coef$estimate[selec_coef$mutation == "P441L"] + effect_size_bay$estimate[effect_size_bay$mutation == "P441L"]
  selec_coef$estimate[selec_coef$mutation == "K13"] <- selec_coef$estimate[selec_coef$mutation == "K13"] + effect_size_bay$estimate[effect_size_bay$mutation == "K13"]

  selec_coef$conf.low[selec_coef$mutation == "A675V"] <- selec_coef$conf.low[selec_coef$mutation == "A675V"] + effect_size_bay$conf.low[effect_size_bay$mutation == "A675V"]
  selec_coef$conf.low[selec_coef$mutation == "C469Y"] <- selec_coef$conf.low[selec_coef$mutation == "C469Y"] + effect_size_bay$conf.low[effect_size_bay$mutation == "C469Y"]
  selec_coef$conf.low[selec_coef$mutation == "C469F"] <- selec_coef$conf.low[selec_coef$mutation == "C469F"] + effect_size_bay$conf.low[effect_size_bay$mutation == "C469F"]
  selec_coef$conf.low[selec_coef$mutation == "P441L"] <- selec_coef$conf.low[selec_coef$mutation == "P441L"] + effect_size_bay$conf.low[effect_size_bay$mutation == "P441L"]
  selec_coef$conf.low[selec_coef$mutation == "K13"] <- selec_coef$conf.low[selec_coef$mutation == "K13"] + effect_size_bay$conf.low[effect_size_bay$mutation == "K13"]

  selec_coef$conf.high[selec_coef$mutation == "A675V"] <- selec_coef$conf.high[selec_coef$mutation == "A675V"] + effect_size_bay$conf.high[effect_size_bay$mutation == "A675V"]
  selec_coef$conf.high[selec_coef$mutation == "C469Y"] <- selec_coef$conf.high[selec_coef$mutation == "C469Y"] + effect_size_bay$conf.high[effect_size_bay$mutation == "C469Y"]
  selec_coef$conf.high[selec_coef$mutation == "C469F"] <- selec_coef$conf.high[selec_coef$mutation == "C469F"] + effect_size_bay$conf.high[effect_size_bay$mutation == "C469F"]
  selec_coef$conf.high[selec_coef$mutation == "P441L"] <- selec_coef$conf.high[selec_coef$mutation == "P441L"] + effect_size_bay$conf.high[effect_size_bay$mutation == "P441L"]
  selec_coef$conf.high[selec_coef$mutation == "K13"] <- selec_coef$conf.high[selec_coef$mutation == "K13"] + effect_size_bay$conf.high[effect_size_bay$mutation == "K13"]

  effect_size_bay$s <- 18
  effect_size_bay$sizes <- 1
  selec_coef$s <- 19
  selec_coef$sizes <- 0.5

  stan_eff_bay_overall <- effect_size_bay %>%
    mutate(method = "Baysian") %>%
    relocate(District, conf.low, estimate, conf.high, method, mutation, s, sizes)
  rownames(stan_eff_bay_overall) <- NULL

  effect_size_bay <- rbind(effect_size_bay, selec_coef)
  effect_size_bay <- effect_size_bay %>%
    group_by(mutation) %>%
    arrange(District)

  effect_size_bay <- effect_size_bay %>%
    mutate(mutation = replace(mutation, mutation == "A675V", "675V")) %>%
    mutate(mutation = replace(mutation, mutation == "C469Y", "469Y")) %>%
    mutate(mutation = replace(mutation, mutation == "C469F", "469F")) %>%
    mutate(mutation = replace(mutation, mutation == "P441L", "441L")) %>%
    mutate(mutation = replace(mutation, mutation == "K13", "K13"))

  effect_size_bay <- effect_size_bay %>%
    relocate(District, mutation) %>%
    mutate(conf.low = round(conf.low, 3),
           conf.high = round(conf.high, 3),
           estimate = round(estimate, 3))
  return(effect_size_bay)
}
