prediction_model <- function(df, mut, pred_year, numofobs=2, freq_cutoff=0.5){
  data_mut <- df %>%
    filter(
      Locus == mut,
      nobs > numofobs,
      year >= min_year & year <= min_year + pred_year) %>%
    group_by(District, Locus) %>%
    filter(freq[year == min_year] < freq_cutoff)

  mut_pred <- rstanarm::stan_glmer(
    lrsmed ~ adj_year + (1 + adj_year | District),
    data = data_mut,
    weights = wt_med,
    adapt_delta = 0.9975
  )

  # forecasting
  forecast_data <- df %>%
    filter(nobs > numofobs,
           is.finite(lrsmed),
           Locus == mut) %>%
    group_by(District, Locus) %>%
    complete(year = seq(min(min_year) + pred_year, 2023, 1)) %>% # this creates all the year values up to 2023
    ungroup() %>%
    select(adj_year, year, District, country, min_year, Locus) %>%
    mutate(min_year = replace(min_year, is.na(min_year), 0))

  forecast_data <- forecast_data %>%
    group_by(District, Locus) %>%
    mutate(min_year = max(min_year),
           adj_year = year - min_year,
           country = unique(na.omit(country))) %>%
    arrange(Locus, District, year) %>%
    ungroup %>%
    group_by(Locus, District, year) %>%
    distinct(year, .keep_all = TRUE) %>%
    ungroup

  sampling <- rstanarm::posterior_predict(mut_pred, type = "response", newdata = forecast_data)
  predictions <- sampling %>% colMeans
  conf_val <- sampling %>% apply(2, quantile, prob = (c(0.025, 0.975)))

  forecast_data <- forecast_data %>%
    mutate(predict = predictions,
           conf.low = conf_val[1,],
           conf.high = conf_val[2,])

  return(forecast_data)
}

plot_prediction <- function(df_predict, df, mut) {
  #obtain dataframe that was used for prediction to then filter the Districts appropriately
  data_mut <- df %>%
    filter(
      Locus == mut,
      nobs > numofobs,
      year >= min_year & year <= min_year + pred_year) %>%
    group_by(District, Locus) %>%
    filter(freq[year == min_year] < freq_cutoff)

  df_pred <- df_predict %>%
    filter(year >= min_year + pred_y,
           Locus == mut,
           District %in% unique(data_mut$District)) %>%
    mutate(District_country = paste0(District, ", ", country, sep=""))
  ex <- df %>%
    filter(is.finite(lrsmed), nobs > numofobs, District %in% unique(data_mut$District), Locus == mut) %>%
    mutate(District_country = paste0(District, ", ", country, sep=""))
  ex_pt <- df %>%
    filter(nobs > numofobs, District %in% unique(data_mut$District), Locus == mut) %>% #, year >= min_year & year <= min_year+4
    mutate(point_col = case_when(year >= min_year & year <= min_year + 2 ~ "#E69F00",
                                 year > min_year + 2 & year <= min_year + 3 ~ "#56B4E9",
                                 year > min_year + 3 & year <= min_year + 4 ~ "#0072B2",
                                 year > min_year + 4 ~ "black",
                                 year < min_year ~ "black"),
           District_country = paste0(District, ", ", country, sep=""))

  ## plot for the mutation mut
  if (mut == "K13"){
    xlabel = "K13 allele frequency"
  }
  else{
    xlabel = paste0(str_sub(mut,2), " allele frequency")
  }

  freq_cutoff_mut_plot <-
    ggplot() +
    geom_line(data=ex, aes(year, predict_stan_freq, color = "true"), linewidth = 1) +
    geom_line(data=df_pred, aes(year, predict_freq, group = label,color = label), linetype = "dashed") +
    geom_point(data = ex_pt, aes(x=year, y=freq, size = n), color = ex_pt$point_col) +
    scale_size_binned(name = "Sample Size",range = c(0.2,2)) +
    labs(shape="Mutations", size = "Sample Size") +
    theme_bw() +
    facet_wrap(~District_country, ncol = 2) +
    ylab(xlabel) +
    xlab("Year") +
    ggpubr::theme_pubclean(base_size = 12) +
    scale_color_manual(
      values = c("#E69F00", "#56B4E9", "#0072B2", "grey"),
      labels = c("2 years", "3 years", "4 years", "true"),
      guide = "none"
    ) +
    theme(
      axis.line = element_line(),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      strip.text.x = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "right",
      panel.spacing = unit(0.2, "lines")
    )

}

correlation <- function(df, df_predict, numofobs=2){
  inital_true <- df %>%
    filter(nobs > numofobs)

  df_corr <- NULL
  for(l in unique(df_predict$label)) {
    #filter forecasting based on year (label)
    pred <- df_predict %>%
      filter(label == l) %>%
      mutate(val = paste0(year, District, Locus, min_year)) %>%
      arrange(val) %>%
      rename(p_pred = predict)
    #filter data points in true that are also in forecasted
    true <- inital_true
    true <- true %>%
      mutate(val = paste0(year, District, Locus, min_year)) %>%
      arrange(val) %>%
      filter(val %in% pred$val)
    #make sure districts match between true and forecasted
    pred <- pred %>% filter(val %in% true$val)
    df <- cbind(true, pred$p_pred, pred$label, pred$pred_y)
    colnames(df)[ncol(df)-2] <- "p_pred"
    colnames(df)[ncol(df)-1] <- "label"
    colnames(df)[ncol(df)] <- "pred_y"
    df_corr <- rbind(df_corr, df) %>%
      filter(year > min_year + pred_y)
  }

  df_corr <- df_corr %>%
    mutate(p_pred_freq = exp(p_pred)/(1 + exp(p_pred)),
           t_pred_freq = freq) %>%
    group_by(Locus, label) %>%
    mutate(mae = round(mean(abs(t_pred_freq - p_pred_freq)), 2), #mean absolute error  -> mean(abs(true-prediction)
           bias = round(mean(t_pred_freq - p_pred_freq), 2), #mean bias -> mean(true - prediction)
           corr = round(cor(t_pred_freq, p_pred_freq), 2)) %>%
    select(country, District, year, Locus, p_pred_freq, t_pred_freq, mae, bias, corr, nobs, adj_year, val , label) %>%
    arrange(Locus, label)
  return(df_corr)
}

plot_correlation <- function(df_corr, mut){
  ## Correlation plot for 50% cutoff in first year
  if (mut == "K13"){
    xlabel = "K13 allele frequency"
  }
  else{
    xlabel = paste0(str_sub(mut,2), " allele frequency")
  }

  correlation_plot <- df_corr %>%
    filter(Locus == mut) %>%
    mutate(District_country = paste0(District, ", ", country, sep="")) %>%
    ggplot(aes(x = t_pred_freq, y = p_pred_freq, color = label, group = label)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
    facet_wrap(~ District_country, ncol = 2) +
    theme_bw() +
    xlab("True") +
    ylab(xlabel) +
    scale_y_continuous(limits=c(0,1), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    scale_x_continuous(limits=c(0,1), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    scale_color_manual(
      values = c("#E69F00", "#56B4E9", "#0072B2"),
      labels = c("2 years", "3 years", "4 years")
    ) +
    theme_bw() +
    ggpubr::theme_pubclean(base_size = 12) +
    theme(
      axis.line = element_line(),
      axis.text = element_text(size = 8),
      axis.title = element_text(size = 10),
      strip.text.x = element_text(size = 10),
      legend.position = "none",
      panel.spacing = unit(0.2, "lines")
    ) +  labs(color = NULL)
  return(correlation_plot)
}

plot_forecasting_confidence <- function(df_predict, pfk7_data, mut){
  mut_df_pred <- df_predict %>%
    filter(year >= min_year + pred_y,
           pred_y == 4,
           Locus == mut) %>%
    mutate(label = "pred",
           point_col = "#0072B2")
  mut_df_true <- df_predict %>%
    filter(year <= min_year + pred_y,
           pred_y == 4,
           Locus == mut) %>%
    mutate(label = "true",
           point_col = "#0072B2")
  mut_ex_pt <- pfk7_data %>%
    filter(nobs > numofobs,
           District %in% unique(mut_df_pred$District),
           Locus == mut) %>%
    mutate(point_col = case_when(year >= min_year & year <= min_year + 4 ~ "#0072B2",
                                 year > min_year + 4 ~ "black",
                                 year < min_year ~ "black"))
  if (mut == "K13"){
    xlabel = "K13 allele frequency"
  }
  else{
    xlabel = paste0(str_sub(mut,2), " allele frequency")
  }

  fig6_plot <-
    ggplot() +
    geom_point(data = mut_ex_pt, aes(x=year, y=freq, group = Locus, size = n), color = mut_ex_pt$point_col) +
    scale_size_binned(name = "Sample \n Size", range = c(0.2,3)) +
    geom_line(data = mut_df_true, aes(year, predict_freq, color = point_col, group = Locus), show_guide=FALSE) +
    geom_line(data = mut_df_pred, aes(year, predict_freq, color = point_col, group = Locus), show_guide=FALSE) +
    geom_ribbon(data = mut_df_pred, aes(ymin = conf.low_freq, ymax = conf.high_freq, x = year, color = point_col, group = Locus, fill = point_col), linetype = "dashed", alpha = 0.2, show_guide=FALSE) +
    theme_bw() +
    facet_wrap(~District) +
    ylab(xlabel) +
    xlab("Year") +
    scale_color_manual(
      values = c("#0072B2", "black"),
      name = "Mutation") +
    scale_fill_manual(
      values = c("#0072B2", "black"),
      name = "Mutation") +
    ggpubr::theme_pubclean(base_size = 12) +
    theme(
      axis.line = element_line(),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 10),
      strip.text.x = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.position = "right",
      panel.spacing = unit(0.2, "lines")
    )
}

points_within_95CrI <- function(df_predict, pfk7_data){
  mut_df_pred <- df_predict %>%
    filter(year >= min_year + pred_y,
           pred_y == 4,
           Locus == "C580Y") %>%
    mutate(label = "pred")

  K13_df_pred <- df_predict %>%
    filter(year >= min_year + pred_y,
           pred_y == 4,
           Locus == "K13") %>%
    mutate(label = "pred")

  mut_ex_pt <- pfk7_data %>%
    filter(nobs > numofobs,
           District %in% unique(mut_df_pred$District),
           Locus == "C580Y")

  K13_ex_pt <- pfk7_data %>%
    filter(nobs > numofobs,
           District %in% unique(mut_df_pred$District),
           Locus == "K13")

  df_pred <- rbind(mut_df_pred, K13_df_pred)
  ex_pt <- rbind(mut_ex_pt, K13_ex_pt)

  pred_data <- df_pred %>% mutate(label = paste0(country, District, year, Locus))
  true_data <- ex_pt %>% mutate(label = paste0(country, District, year, Locus))

  true_freq <- true_data[true_data$label %in% pred_data$label,] %>%
    mutate(freq = freq) %>%
    select(freq)
  pred_data_filter <- pred_data[pred_data$label %in% true_data$label,]
  pred_data_filter <- cbind(pred_data_filter, true_freq$freq)
  colnames(pred_data_filter)[15] <- "freq"

  pred_data_filter <- pred_data_filter %>%
    group_by(Locus) %>%
    mutate(val = case_when(freq >= conf.low_freq & freq <= conf.high_freq ~ "yes",
                           freq < conf.low_freq | freq > conf.high_freq ~ "no"),
           perc = sum(val == "yes")/(sum(val == "yes") + sum(val == "no")),
           num = sum(val == "yes"),
           den = sum(val == "yes") + sum(val == "no")) %>%
    distinct(perc, .keep_all = T) %>%
    select(Locus, perc, num, den) %>%
    ungroup
  return(pred_data_filter)
}
