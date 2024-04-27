library(dplyr)
library(tidyr)
library(ggplot2)

#data files in data-raw obtained from MalariaGen PfK7 database
x <- read.csv("analysis/data/data-raw/Pf7-samples.csv")
g <- read.delim("analysis/data/data-raw/genotypes.txt")
names(g)[1] <- "sample_id"

head(x)
table(x$ARTresistant)
table(g$kelch13_349.726_ns_changes)

cam <-
  x %>% dplyr::filter(country %in% c("Cambodia", "Myanmar", "Laos", "Vietnam", "Thailand")) %>%
  dplyr::left_join(g, by = "sample_id")

cam <- cam %>%
  filter(qc_pass,
         !(kelch13_349.726_ns_changes %in% c("!", "*", "-", "!*"))) %>%
  arrange(country, site, year) %>%
  group_by(country, site, year) %>%
  summarise(C580Y = sum(grepl("580Y|580y", kelch13_349.726_ns_changes)),
            R539T = sum(grepl("539t|539T", kelch13_349.726_ns_changes)),
            K13 = sum(grepl("446|458|469|476|493|539t|539T|543|553|561|574|580Y|580y|622|675", kelch13_349.726_ns_changes)),
            other_K13 = sum(grepl("446|458|469|476|493|543|553|561|574|622|675", kelch13_349.726_ns_changes)),
            n = n()) %>%
  ungroup %>%
  pivot_longer(cols = C580Y:other_K13) %>%
  rename(Locus = name, x = value) %>%
  select(country, site, year, n, x, Locus)

pfk7_data <- cam %>%
  group_by(country, site, Locus) %>%
  filter(is.finite(x)) %>%
  mutate(min_year = min(year[x>0])) %>%
  mutate(nobs= sum(x > 0)) %>%
  ungroup %>%
  mutate(adj_year = year - min_year,
         prev = x / n) %>%
  filter(is.finite(min_year)) %>%
  arrange(Locus, country, site, year)

pfk7_data <- pfk7_data %>%
  group_by(country, site, Locus) %>%
  mutate(min_year = min(year[x>0])) %>%
  mutate(nobs= sum(x > 0)) %>%
  ungroup %>%
  mutate(adj_year = year - min_year,
         prev = x / n) %>%
  filter(is.finite(min_year))

# log ratio function
se_ln_ratio_noZeros <- function(x, N) {
  inds <- which(x == 0)
  x[inds] <- 0.5
  inds <- which(x == N)
  x[inds] <- x[inds] - 0.5
  se <- sqrt(1 / x + 1 / (N - x))
  return(se)
}

numofobs <- 2
pfk7_data$wt_med <- NA

pfk7_data <-
  pfk7_data %>% rename(District = site) %>%
  mutate(lrsmed = log(prev / (1 - prev)),
         wt_med = replace(wt_med, nobs > numofobs & is.finite(lrsmed),
                          1/(se_ln_ratio_noZeros(round(prev[nobs > numofobs & is.finite(lrsmed)]*n[nobs > numofobs & is.finite(lrsmed)]), n[nobs > numofobs & is.finite(lrsmed)])^2)
         ))

pfk7_data$lrsmed[pfk7_data$lrsmed == Inf] <- NA
pfk7_data$lrsmed[pfk7_data$lrsmed == -Inf] <- NA

write.table(pfk7_data,"analysis/data/data-derived/pfk7_data.txt")
