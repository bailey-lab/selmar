library(dplyr)
library(ggplot2)

#data files in data-raw obtained from MalariaGen PfK7 database
x <- read.csv("analysis/data/data-raw/Pf7-samples.csv")
g <- read.delim("analysis/data/data-raw/genotypes.txt")
names(g)[1] <- "sample_id"

# create our data from GMS
cam <-
  x %>% dplyr::filter(country %in% c("Cambodia", "Myanmar", "Laos", "Vietnam", "Thailand")) %>%
  dplyr::left_join(g, by = "sample_id")

# New filter to calucalte mutation prevalences
pfk7_data <- cam %>%
  # remove those that fail qc
  filter(qc_pass) %>%
  # remove those that kelch could not be reliably ascertained
  filter(!(kelch13_349.726_ns_changes %in% c("-", "!", "!*"))) %>%
  select(sample_id, year, country, site, kelch13_349.726_ns_changes) %>%
  arrange(country, site, year) %>%
  group_by(country, site, year) %>%
  summarise(C580Y = sum(grepl("580Y|580y", kelch13_349.726_ns_changes)),
            R539T = sum(grepl("539t|539T", kelch13_349.726_ns_changes)),
            PfK13 = sum(grepl("580|446|458|476|493|539|543|553|561", kelch13_349.726_ns_changes)),
            N = n()) %>%
  pivot_longer(cols = C580Y:PfK13) %>%
  rename(Locus = name, x = value) %>%
  select(country, site, year, N, x, Locus)

# group and work out our min years
pfk7_data <- pfk7_data %>%
  group_by(country, site, Locus) %>%
  filter(is.finite(x)) %>%
  mutate(min_year = min(year[x>0])) %>%
  mutate(nobs= sum(x > 0)) %>%
  ungroup %>%
  mutate(adj_year = year - min_year,
         prev = x / N) %>%
  filter(is.finite(min_year)) %>%
  arrange(Locus, country, site, year)

# log ratio function
se_ln_ratio_noZeros <- function(x, N) {
  inds <- which(x == 0)
  x[inds] <- 0.5
  inds <- which(x == N)
  x[inds] <- x[inds] - 0.5
  se <- sqrt(1 / x + 1 / (N - x))
  return(se)
}

# filter to those with more than 2 ime points
numofobs <- 2
pfk7_data <-
  pfk7_data %>% rename(District = site) %>% rename(n = N) %>%
  mutate(lrsmed = log(prev / (1 - prev)))

# calculate our weights
pfk7_data$wt_med <- NA

# Thought I would replace the following below with a
# dplyr style that might be easier to follow
# pfk7_data$wt_med[pfk7_data$nobs > numofobs & is.finite(pfk7_data$lrsmed)] <- 1/(se_ln_ratio_noZeros(round(pfk7_data$prev[pfk7_data$nobs > numofobs & is.finite(pfk7_data$lrsmed)]*pfk7_data$n[pfk7_data$nobs > numofobs & is.finite(pfk7_data$lrsmed)]), pfk7_data$n[pfk7_data$nobs > numofobs & is.finite(pfk7_data$lrsmed)])^2)

pfk7_data <- pfk7_data %>%
  mutate(wt_med = replace(
    wt_med,
    nobs > numofobs & is.finite(lrsmed),
    1/(se_ln_ratio_noZeros(x[nobs > numofobs & is.finite(lrsmed)],
                          n[nobs > numofobs & is.finite(lrsmed)]) ^ 2)
    )
  )

write.table(pfk7_data,"analysis/data/data-derived/NEW_pfk7_data.txt")


