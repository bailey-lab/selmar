# ------------------------------------------------
# DATA IMPORT
# ------------------------------------------------
library(tidyverse)
devtools::load_all()

# install couple of packages to help with scraping info from pdfs
#remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))
out2 <- tabulizer::extract_tables("analysis/NEJM_submission.pdf", pages = 15, guess = TRUE, output = "data.frame")

# get the data from the scraped table
df <- out2[[1]]
df[30:31,] <- df[30:31,c(ncol(df),seq_len(ncol(df)-1))]
df[df==""] <- NA
names(df) <- gsub("X","",names(df))
df$p <- NULL

df <- df %>% fill(District, .direction = "down")
for(i in 1:6){
  df[,i+2] <- gsub("\\(.*\\)", "", df[,i+2])
}

df <- df %>% pivot_longer(cols = `2016`:`2021`) %>%
  mutate(n = as.numeric(gsub("(.*)/(.*)","\\2", value))) %>%
  mutate(x = as.numeric(gsub("(.*)/(.*)","\\1", value))) %>%
  mutate(med = x/n) %>%
  na.omit %>%
  mutate(ub = Hmisc::binconf(x,n)[,3]) %>%
  rowwise() %>%
  mutate(lrsmed = log(med/(1-med))) %>%
  mutate(lrsub = log(ub/(1-ub))) %>%
  rename(year = name) %>%
  mutate(year = as.numeric(year)) %>%
  group_by(District, Locus) %>%
  mutate(nobs = sum(x>0)) %>%
  ungroup

# quick check the table seems correct
df

# lucys log ratio function
se_ln_ratio_noZeros<-function(x,N) {
  inds<-which(x==0)
  x[inds]<-0.5
  inds<-which(x==N)
  x[inds]<-x[inds]-0.5
  se<- sqrt(1/x + 1/(N-x))
  return(se)
}

df$wt_med <- 1/(se_ln_ratio_noZeros(round(df$med*df$n), df$n)^2)
df$wt_ub <- 1/(se_ln_ratio_noZeros(round(df$ub*df$n), df$n)^2)
# remove the data that is Inf here due to 0s - this will obscure it from mdel fitting
df$lrsmed[df$lrsmed == -Inf] <- NA

# ------------------------------------------------
# Fig1 Data Overview Plot
# ------------------------------------------------

binomial_smooth <- function(...) {
  geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

# quick plot of the trends
gg1 <- df %>%
  ggplot(aes(year, lrsmed, color = District)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Locus) +
  ylab("log(Kelch Prevalence / Wild-Type Prevalence)") +
  xlab("Year") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.position = "right",
        axis.text.x = element_text(size = 8))

gg2 <- df %>%
  ggplot(aes(year, med, color = District)) +
  geom_point() +
  binomial_smooth(se = FALSE) +
  facet_wrap(~Locus) +
  ylab("Kelch Mutation Prevalence") +
  xlab("Year") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.position = "right",
        axis.text.x = element_text(size = 8))

# Fig 1 Style Plot
cowplot::plot_grid(gg2 + theme(legend.position = "none"),
                   gg1 + theme(legend.position = "none"),
                   cowplot::get_legend(gg1 + theme(legend.key = element_rect(fill="white"))), ncol = 3,
                   rel_widths = c(1,1,0.2))


# ------------------------------------------------
# Modelling
# ------------------------------------------------

# Few Example Models to showcase different approaches just based on c469

# lme 4 with slope and intercept
c469y_mod <-
  lme4::lmer(
    lrsmed ~ year + (1 + year | District),
    data = df %>% filter(Locus == "C469Y" & nobs > 1) %>% mutate(year = year - 2016),
    weights = n
  )

# lme 4 without slope
c469y_mod_rint <-
  lme4::lmer(
    lrsmed ~ year + (1 |  District),
    data = df %>% filter(Locus == "C469Y" & nobs > 1) %>% mutate(year = year - 2016),
    weights = n
  )

# Bayesian stan model with slopes and intercepts
# Have a read through rstanarm vignettes and guides on how it works in comparison to lme4
c469y_stanmod <- rstanarm::stan_glmer(
  lrsmed ~ year + (1 + year | District),
  data = df %>% filter(Locus == "C469Y" & nobs > 1) %>% mutate(year = year - 2016),
  weights = wt_med, adapt_delta = 0.99
)

a675v_stanmod <- rstanarm::stan_glmer(
  lrsmed ~ year + (1 + year | District),
  data = df %>% filter(Locus == "A675V" & nobs > 1) %>% mutate(year = year - 2016),
  weights = wt_med, adapt_delta = 0.995
)

c469f_stanmod <- rstanarm::stan_glmer(
  lrsmed ~ year + (1 + year | District),
  data = df %>% filter(Locus == "C469F" & nobs > 1) %>% mutate(year = year - 2016),
  weights = wt_med, adapt_delta = 0.9975
)

# selection coefficients overall:
summary(c469y_mod)

# Random effects for lme model, which are perfectly correlated
lme4::ranef(c469y_mod)
lme4::ranef(c469y_mod)[[1]] %>%
  ggplot(aes(`(Intercept)`, year)) +
  geom_point() +
  theme_bw() +
  geom_line() +
  xlab("Random Intercept") +
  ylab("Random Slope")

# Random effects for stanmod
vals <- summary(c469y_stanmod$stanfit)$summary[,1]
vals <- vals[grep("^b", names(vals))]
data.frame("Intercepts" = vals[seq(1, length(vals), 2)],
           "Slopes" = vals[seq(2, length(vals), 2)]) %>%
  ggplot(aes(Intercepts, Slopes)) +
  geom_point() +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  xlab("Random Intercept") +
  ylab("Random Slope")


#  get model predictions for plotting
df$predict_rint <- NA
df$predict_stan <- NA

# here just add predictions  for the c469Y rows
df$predict_rint[is.finite(df$lrsmed) & df$Locus == "C469Y" & df$nobs > 1] <-  predict(c469y_mod_rint, type = "response")
df$predict_stan[is.finite(df$lrsmed) & df$Locus == "C469Y" & df$nobs > 1] <-  rstanarm::posterior_predict(c469y_stanmod, type = "response") %>% colMeans

# Plot the predictions in each site against the linear regression (geom_smooth with lm method)
df %>%
  filter(Locus == "C469Y" & is.finite(lrsmed) & nobs > 1) %>%
  ggplot(aes(year, lrsmed)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  geom_line(aes(year, predict_rint), color = "blue") +
  geom_line(aes(year, predict_stan), color = "green") +
  facet_wrap(~District) +
  ylab("log(C469Y Prevalence / Wild-Type Prevalence)") +
  xlab("Year") +
  ggpubr::theme_pubclean(base_size = 12) +
  theme(axis.line = element_line(), legend.position = "right",
        axis.text.x = element_text(size = 8),
        title = element_text(size = 10)) +
  ggtitle("Blue: Predictions from model with random intercept only \nRed: Linear regression within each site \nGreen: Bayesian Model")

# Plot the effect sizes
rint_eff <- broom.mixed::tidy(c469y_mod_rint, conf.int = 0.95) %>%
  filter(term == "year") %>%
  dplyr::select(conf.low, estimate, conf.high) %>%
  mutate(method = "Random Intercept Only")
stan_eff <- as.data.frame(t(data.frame(c469y_stanmod$stan_summary[2, c("2.5%", "50%","97.5%")]))) %>%
  mutate(method = "Bayesian (Intercept + Slope)") %>%
  setNames(c("conf.low", "estimate", "conf.high", "method"))

# plot the effect sizes for the mutation
rbind(rint_eff, stan_eff) %>% mutate(mutation = "C469Y") %>%
  ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = mutation, color = method, group = method)) +
  geom_pointrange(position = position_dodge(width = 0.1)) +
  ylab("Mutation") +
  xlab("Selection Coefficient") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  scale_color_discrete(name = "Method")



# forecasting
new_data <- df %>%
  filter(Locus == "C469Y" & nobs > 1) %>% # this is for c469Y model
  group_by(District, Locus) %>% # we group by the District and Locus
  complete(year = seq(min(year), 2023, 1)) %>% # this creates all the year values up to 2023
  mutate(year = year-2016) %>% # mutate our years to the scale of the model
  select(year, District, Locus) # select just the columns needed for the model

rstanarm::posterior_predict(c469y_stanmod, type = "response", newdata = new_data) %>% colMeans
