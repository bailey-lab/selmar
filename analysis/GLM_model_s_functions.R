GLM_model_mutant <- function(df, locus, numofobs){
  stanmod <- rstanarm::stan_glmer(
    lrsmed ~ adj_year + (1 + adj_year | District),
    data = df %>%
      filter(Locus == locus,
             nobs > numofobs,
             is.finite(lrsmed)),
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
