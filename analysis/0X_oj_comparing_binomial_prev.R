#set the limit of number of observations
numofobs <- 2
#load data
pfk7_data <- read.table("analysis/data/data-derived/pfk7_data.txt") %>% filter(Locus != "other_PfK13")

pred_year <- i <-  4

#------ Model prediction for C580Y ------
data_c580y <- pfk7_data %>%
  filter(
    Locus == "C580Y",
    nobs > numofobs,
    # This was not filtering to all occurrences  as you were doing ==
    District %in% c("Pailin", "Preah Vihear", "Ratanakiri", "Tak"))

selected_district <- data_c580y %>%
  group_by(Locus, country, District) %>%
  filter(adj_year == min(adj_year),
         prev <= prev_cutoff) %>%
  ungroup %>%
  select(District)

selected_district <- array(t(selected_district)) #check if the prevalence in the first year observed is less than cutoff

data_c580y <- data_c580y %>%
  filter(District %in% selected_district, #keep only the district that satisfy filtering
         year >= min_year & year <= min_year + i)

c580y_pred <- rstanarm::stan_glmer(
  cbind(x,n-x) ~ adj_year + (1 + adj_year | District),
  data = data_c580y,
  family = binomial(link = "logit"),
  adapt_delta = 0.9975
)

# forecasting
new_data <- pfk7_data %>%
  filter(nobs > numofobs,
         is.finite(pfk7_data$lrsmed),
         District %in% selected_district,
         Locus == "C580Y") %>%
  group_by(District, Locus) %>%
  complete(year = seq(min(min_year) + i, 2023, 1)) %>% # this creates all the year values up to 2023
  ungroup() %>%
  select(adj_year, year, District, country, min_year, Locus, n) %>%
  mutate(min_year = replace(min_year, is.na(min_year), 0))

new_data <- new_data %>%
  group_by(District, Locus) %>%
  mutate(min_year = max(min_year),
         adj_year = year - min_year,
         country = unique(na.omit(country))) %>%
  arrange(Locus, District, year) %>%
  ungroup %>%
  group_by(Locus, District, year) %>%
  distinct(year, .keep_all = TRUE) %>%
  ungroup

# fill this up so we get values for it
new_data$n <- 100
new_data$x <- 0

rstanarm::posterior_predict(c580y_pred, newdata = new_data %>% na.omit %>% mutate(x = 0)) -> predsx

predictions <-
  rstanarm::posterior_predict(c580y_pred, type = "response", newdata = new_data) %>% colMeans
conf_val <- rstanarm::posterior_predict(c580y_pred, type = "response", newdata = new_data) %>%
  apply(2, quantile, prob = (c(0.025, 0.975)))

new_data <- new_data %>%
  mutate(predict = predictions/n,
         conf.low = conf_val[1,]/n,
         conf.high = conf_val[2,]/n)

df_predictbin <- new_data
ggplot(df_predictbin, aes(year, predict, ymin = conf.low, ymax = conf.high)) +
  geom_point() + geom_ribbon(alpha = 0.2) + facet_wrap(~District)


# and compare against other approach
c580y_pred <- rstanarm::stan_glmer(
  lrsmed ~ adj_year + (1 + adj_year | District),
  data = data_c580y,
  weights = wt_med,
  adapt_delta = 0.9975
)

# forecasting
new_data <- pfk7_data %>%
  filter(nobs > numofobs,
         is.finite(pfk7_data$lrsmed),
         District %in% selected_district,
         Locus == "C580Y") %>%
  group_by(District, Locus) %>%
  complete(year = seq(min(min_year) + i, 2023, 1)) %>% # this creates all the year values up to 2023
  ungroup() %>%
  select(adj_year, year, District, country, min_year, Locus, n) %>%
  mutate(min_year = replace(min_year, is.na(min_year), 0))

new_data <- new_data %>%
  group_by(District, Locus) %>%
  mutate(min_year = max(min_year),
         adj_year = year - min_year,
         country = unique(na.omit(country))) %>%
  arrange(Locus, District, year) %>%
  ungroup %>%
  group_by(Locus, District, year) %>%
  distinct(year, .keep_all = TRUE) %>%
  ungroup

predictions <-
  rstanarm::posterior_predict(c580y_pred, type = "response", newdata = new_data) %>% colMeans
conf_val <- rstanarm::posterior_predict(c580y_pred, type = "response", newdata = new_data) %>%
  apply(2, quantile, prob = (c(0.025, 0.975)))

new_data <- new_data %>%
  mutate(predict = predictions,
         conf.low = conf_val[1,],
         conf.high = conf_val[2,])

df_predict <- new_data

df_predict %>% mutate(predict = exp(predict) / (1 + exp(predict))) %>%
mutate(conf.low = exp(conf.low) / (1 + exp(conf.low))) %>%
         mutate(conf.high = exp(conf.high) / (1 + exp(conf.high))) %>%

ggplot(aes(year, predict, ymin = conf.low, ymax = conf.high)) +
  geom_point() + geom_ribbon(alpha = 0.2) + facet_wrap(~District)

## compare fully


rbind(
  df_predictbin %>% mutate(method = "binomial") %>%
    select(year, predict, District, method),
  df_predict %>% mutate(method = "prev") %>% mutate(predict = exp(predict) / (1 + exp(predict))) %>%
    select(year, predict, District, method)
) %>%
  ggplot(aes(year, predict, color = method)) + geom_point() + facet_wrap( ~ District)
