#Preprocessing for Figure 3 plotting the selection coefficients
#Input: GLM models for each mutation
#Output: Selection coefficients for each district and mutations

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
