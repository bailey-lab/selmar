library(dplyr)
library(ggplot2)

##### Description of File #####
# This file takes the preprocessed PfK7 data and combines the data for each country for the forecasting analysis of SE Asia
####


#load preprocessed data for SE Asia (PfK7 dataset)
pfk7_data <- read.table("analysis/data/data-derived/pfk7_data.txt")

pfk7_country_data <-
  pfk7_data %>% select(country, District, year, n, c580y, r539t, all_kelch) %>% group_by(country, year) %>%
  mutate(c580y = sum(c580y),
    r539t = sum(r539t),
    all_kelch = sum(all_kelch),
    n = sum(n)) %>%
  distinct(year, .keep_all = TRUE) %>%
  select(!District) %>%
  mutate(prev_580y = c580y / n,
    prev_539t = r539t / n,
    prev_all_kelch = all_kelch / n) %>%
  group_by(country) %>%
  mutate(min_year = min(year)) %>%
  mutate(nobs_c580y = sum(c580y > 0),
         nobs_r539t = sum(r539t > 0),
         nobs_kelch = sum(all_kelch > 0))

#log ratio function
se_ln_ratio_noZeros <- function(x, N) {
  inds <- which(x == 0)
  x[inds] <- 0.5
  inds <- which(x == N)
  x[inds] <- x[inds] - 0.5
  se <- sqrt(1 / x + 1 / (N - x))
  return(se)
  }

numofobs <- 2

pfk7_country_data <- pfk7_country_data %>%
  mutate(lrsmed_580y = log(prev_580y / (1 - prev_580y))) %>%
  mutate(lrsmed_all_kelch = log(prev_all_kelch / (1 - prev_all_kelch))) %>%
  mutate(lrsmed_539t = log(prev_539t / (1 - prev_539t)))
pfk7_country_data <- pfk7_country_data %>%
  group_by(country) %>%
  mutate(min_year = min(year)) %>%
  ungroup %>%
  mutate(adj_year = year - min_year)

pfk7_country_data$lrsmed_539t[pfk7_country_data$lrsmed_539t == Inf] <- NA
pfk7_country_data$lrsmed_580y[pfk7_country_data$lrsmed_580y == Inf] <- NA
pfk7_country_data$lrsmed_all_kelch[pfk7_country_data$lrsmed_all_kelch == Inf] <- NA

pfk7_country_data$lrsmed_539t[pfk7_country_data$lrsmed_539t == -Inf] <- NA
pfk7_country_data$lrsmed_580y[pfk7_country_data$lrsmed_580y == -Inf] <- NA
pfk7_country_data$lrsmed_all_kelch[pfk7_country_data$lrsmed_all_kelch == -Inf] <- NA

pfk7_country_data$wt_med_539t <- NA
pfk7_country_data$wt_med_580y <- NA
pfk7_country_data$wt_med_all_kelch <- NA

x <- round(pfk7_country_data$prev_580y[pfk7_country_data$nobs_c580y > numofobs & is.finite(pfk7_country_data$lrsmed_580y)] * pfk7_country_data$n[pfk7_country_data$nobs_c580y > numofobs & is.finite(pfk7_country_data$lrsmed_580y)])
n <- pfk7_country_data$n[pfk7_country_data$nobs_c580y > numofobs & is.finite(pfk7_country_data$lrsmed_580y)]
pfk7_country_data$wt_med_580y[pfk7_country_data$nobs_c580y > numofobs & is.finite(pfk7_country_data$lrsmed_580y)] <- 1 / (se_ln_ratio_noZeros(x, n) ^ 2)

x <- round(pfk7_country_data$prev_539t[pfk7_country_data$nobs_r539t > numofobs & is.finite(pfk7_country_data$lrsmed_539t)] * pfk7_country_data$n[pfk7_country_data$nobs_r539t > numofobs & is.finite(pfk7_country_data$lrsmed_539t)])
n <- pfk7_country_data$n[pfk7_country_data$nobs_r539t > numofobs & is.finite(pfk7_country_data$lrsmed_539t)]
pfk7_country_data$wt_med_539t[pfk7_country_data$nobs_r539t > numofobs & is.finite(pfk7_country_data$lrsmed_539t)] <- 1 / (se_ln_ratio_noZeros(x, n) ^ 2)

x <- round(pfk7_country_data$prev_all_kelch[pfk7_country_data$nobs_kelch > numofobs & is.finite(pfk7_country_data$lrsmed_all_kelch)] * pfk7_country_data$n[pfk7_country_data$nobs_kelch > numofobs & is.finite(pfk7_country_data$lrsmed_all_kelch)])
n <- pfk7_country_data$n[pfk7_country_data$nobs_kelch > numofobs & is.finite(pfk7_country_data$lrsmed_all_kelch)]
pfk7_country_data$wt_med_all_kelch[pfk7_country_data$nobs_kelch > numofobs & is.finite(pfk7_country_data$lrsmed_all_kelch)] <-1 / (se_ln_ratio_noZeros(x, n) ^ 2)


write.table(pfk7_country_data,
            "analysis/data/data-derived/pfk7_country_data.txt")
