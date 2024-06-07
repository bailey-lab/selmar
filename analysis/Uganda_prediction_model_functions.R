prediction_model <- function(df, mut, numofobs=2, pred_year){
  stan_pred <- rstanarm::stan_glmer(
    lrsmed ~ adj_year + (1 + adj_year | District),
    data = df %>%
      filter(Locus == mut & nobs > numofobs),
    weights = wt_med,
    adapt_delta = 0.9975
  )

  new_data <- df %>%
    filter(nobs > numofobs,
           is.finite(df$lrsmed),
           Locus == mut) %>% # this is for c469Y model
    group_by(District, Locus) %>% # we group by the District and Locus
    complete(year = seq(max(year), pred_year, 1)) %>% # this creates all the year values up to pred_year
    select(adj_year, year, District, Locus, min_year, max_year) %>%
    mutate(min_year = replace(min_year, is.na(min_year), 0),
           max_year = replace(max_year, is.na(max_year), 0))

  new_data <- new_data %>%
    group_by(District, Locus) %>%
    mutate(min_year = max(min_year),
           adj_year = year - min_year) %>%
    arrange(Locus, District, year) %>%
    ungroup

  new_data$predict <- rstanarm::posterior_predict(stan_pred, type = "response", newdata = new_data) %>% colMeans
  conf_val <- rstanarm::posterior_predict(stan_pred, type = "response", newdata = new_data) %>%
    apply(2, quantile, prob = (c(0.025, 0.975)))
  new_data$conf.low <- conf_val[1,]
  new_data$conf.high <- conf_val[2,]

  df_predict <- new_data
  return(df_predict)
}

plot_prediction <- function(df_predict, df, mut){
  df_true <- df_predict %>%
    filter(Locus == mut,
           year <= max_year) %>%
    mutate(label = "true")
  df_pred <- df_predict %>%
    filter(Locus == mut,
           year >= max_year) %>%
    mutate(label = "predicted")
  ex_pt <- df %>%
    filter(Locus == mut,
           nobs > numofobs,
           District %in% unique(df_pred$District))

  if (mut == "K13"){
    xlabel = "K13 allele frequency"
  }
  else{
    xlabel = paste0(str_sub(mut,2), " allele frequency")
  }

  pred_plot <- ggplot() +
    geom_line(data = df_true, aes(year, predict_freq), color = "#CC79A7", show.legend = FALSE) +
    geom_line(data = df_pred, aes(year, predict_freq), color = "#CC79A7", show.legend = FALSE) +
    geom_ribbon(data = df_pred, aes(ymin = conf.low_freq, ymax = conf.high_freq, x = year), color = "#CC79A7", fill = "#CC79A7", linetype = "dashed", alpha = 0.2, show.legend = FALSE) +
    geom_point(data = ex_pt, aes(x = year, y = freq * 100, size = n), color = "#CC79A7") +
    scale_size_binned(name = "Sample Size",range = c(0.2,1.5)) +
    theme_bw() +
    facet_wrap(~District) +
    ylab(xlabel) +
    xlab("Year") +
    scale_color_manual(
      values = c("#CC79A7"),
      labels = c("Forecasting")
    ) +
    scale_fill_manual(
      values = c("#CC79A7"),
      labels = c("Forecasting")
    ) +
    ggpubr::theme_pubclean(base_size = 12) +
    theme(
      axis.line = element_line(),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 10),
      strip.text.x = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "right",
      panel.spacing = unit(0.2, "lines"))
  return(pred_plot)
}
