library(dplyr)
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

table(cam$kelch13_349.726_ns_changes)
table(cam$kelch13_349.726_ns_changes[which(cam$site == "Pailin")])

cam$c580y <- NA
cam$c580y[which(cam$ARTresistant == "resistant")] <- 0
cam$c580y[grep("580Y|580y", cam$kelch13_349.726_ns_changes)] <- 1
df_c580y <- cam %>%
  select(sample_id, year, country, site, c580y) %>%
  arrange(country, site, year) %>%
  group_by(country, site, year) %>%
  summarise(N = n(),
            c580y = if(all(is.na(c580y))){NA} else{sum(c580y, na.rm = T)}) %>%
  rename(x = c580y) %>%
  mutate(Locus = "C580Y") %>%
  ungroup

cam$r539t <- NA
cam$r539t[which(cam$ARTresistant == "resistant")] <- 0
cam$r539t[grep("539t|539T", cam$kelch13_349.726_ns_changes)] <- 1
df_r539t <- cam %>%
  select(sample_id, year, country, site, r539t) %>%
  arrange(country, site, year) %>%
  group_by(country, site, year) %>%
  summarise(N = n(),
            r539t = if(all(is.na(r539t))){NA} else{sum(r539t, na.rm = T)}) %>%
  rename(x = r539t) %>%
  mutate(Locus = "R539T") %>%
  ungroup

## all validated markers.
cam$all_kelch <- NA
cam$all_kelch[which(cam$ARTresistant == "resistant")] <- 0
cam$all_kelch[grep("580Y|580y|446|458|476|493|539|543|553|561",
                   cam$kelch13_349.726_ns_changes)] <- 1
df_all_kelch <- cam %>%
  select(sample_id, year, country, site, all_kelch) %>%
  arrange(country, site, year) %>%
  group_by(country, site, year) %>%
  summarise(N = n(),
            all_kelch = if(all(is.na(all_kelch))){NA} else{sum(all_kelch, na.rm = T)}) %>%
  rename(x = all_kelch) %>%
  mutate(Locus = "PfK13") %>%
  ungroup

## all other validated markers besides 580Y and 539T
cam$other_kelch <- NA
cam$other_kelch[which(cam$ARTresistant == "resistant")] <- 0
cam$other_kelch[grep("446|458|476|493|543|553|561",
                   cam$kelch13_349.726_ns_changes)] <- 1
df_other_kelch <- cam %>%
  select(sample_id, year, country, site, other_kelch) %>%
  arrange(country, site, year) %>%
  group_by(country, site, year) %>%
  summarise(N = n(),
            other_kelch = if(all(is.na(other_kelch))){NA} else{sum(other_kelch, na.rm = T)}) %>%
  rename(x = other_kelch) %>%
  mutate(Locus = "other_PfK13") %>%
  ungroup

pfk7_data <- rbind(df_c580y, df_r539t, df_all_kelch, df_other_kelch)
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

numofobs <- 2

pfk7_data <-
  pfk7_data %>% rename(District = site) %>% rename(n = N) %>%
  mutate(lrsmed = log(prev / (1 - prev)))

pfk7_data$lrsmed[pfk7_data$lrsmed == Inf] <- NA
pfk7_data$lrsmed[pfk7_data$lrsmed == -Inf] <- NA

pfk7_data$wt_med <- NA
pfk7_data$wt_med[pfk7_data$nobs > numofobs & is.finite(pfk7_data$lrsmed)] <- 1/(se_ln_ratio_noZeros(round(pfk7_data$prev[pfk7_data$nobs > numofobs & is.finite(pfk7_data$lrsmed)]*pfk7_data$n[pfk7_data$nobs > numofobs & is.finite(pfk7_data$lrsmed)]), pfk7_data$n[pfk7_data$nobs > numofobs & is.finite(pfk7_data$lrsmed)])^2)

write.table(pfk7_data,"analysis/data/data-derived/NEW_pfk7_data.txt")
