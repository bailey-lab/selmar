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



cam$c580y <- 0
cam$c580y[grep("580Y|580y", cam$kelch13_349.726_ns_changes)] <- 1
cam$c580y[which(cam$kelch13_349.726_ns_changes == "-")] <- NA
cam$c580y[which(cam$ARTresistant == "null" |
                  cam$ARTresistant == "undetermined")] <- NA

cam$r539t <- 0
cam$r539t[grep("539t|539T", cam$kelch13_349.726_ns_changes)] <- 1
cam$r539t[which(cam$kelch13_349.726_ns_changes == "-")] <- NA
cam$r539t[which(cam$ARTresistant == "null" |
                  cam$ARTresistant == "undetermined")] <- NA

## all validated markers.
cam$all_kelch <- 0
cam$all_kelch[which(cam$kelch13_349.726_ns_changes == "-")] <- NA
cam$all_kelch[grep("580Y|580y|446|458|476|493|539|543|553|561",
                   cam$kelch13_349.726_ns_changes)] <- 1
cam$all_kelch[which(cam$ARTresistant == "null" |
                      cam$ARTresistant == "undetermined")] <- NA

pailin <- cam %>% dplyr::filter(site == "Pailin")
bago <- cam %>% dplyr::filter(site == "Bago")

pfk7_data <- cam %>%
  dplyr::group_by(country, site, year) %>%
  dplyr::summarise(
    N = n(),
    c580y = sum(c580y, na.rm = T),
    r539t = sum(r539t, na.rm = T),
    all_kelch = sum(all_kelch, na.rm = T)
  ) %>%
  ungroup %>%
  group_by(site) %>%
  mutate(min_year = min(year)) %>%
  mutate(nobs_c580y = sum(c580y > 0)) %>%
  mutate(nobs_r539t = sum(r539t > 0)) %>%
  mutate(nobs_kelch = sum(all_kelch > 0)) %>%
  ungroup %>%
  mutate(adj_year = year - min_year) %>%
  dplyr::mutate(
    prev_580y = c580y / N,
    prev_539t = r539t / N,
    prev_all_kelch = all_kelch / N,
    country_site = paste0(country, site, sep = "-")
  )

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
  mutate(lrsmed_580y = log(prev_580y / (1 - prev_580y))) %>%
  mutate(lrsmed_all_kelch = log(prev_all_kelch / (1 - prev_all_kelch))) %>%
  mutate(lrsmed_539t = log(prev_539t / (1 - prev_539t)))
pfk7_data <- pfk7_data %>%
  group_by(District) %>%
  mutate(min_year = min(year)) %>%
  ungroup %>%
  mutate(adj_year = year - min_year)

pfk7_data$lrsmed_539t[pfk7_data$lrsmed_539t == Inf] <- NA
pfk7_data$lrsmed_580y[pfk7_data$lrsmed_580y == Inf] <- NA
pfk7_data$lrsmed_all_kelch[pfk7_data$lrsmed_all_kelch == Inf] <- NA

pfk7_data$lrsmed_539t[pfk7_data$lrsmed_539t == -Inf] <- NA
pfk7_data$lrsmed_580y[pfk7_data$lrsmed_580y == -Inf] <- NA
pfk7_data$lrsmed_all_kelch[pfk7_data$lrsmed_all_kelch == -Inf] <- NA

pfk7_data$wt_med_539t <- NA
pfk7_data$wt_med_580y <- NA
pfk7_data$wt_med_all_kelch <- NA

pfk7_data$wt_med_580y[pfk7_data$nobs_c580y > numofobs & is.finite(pfk7_data$lrsmed_580y)] <- 1/(se_ln_ratio_noZeros(round(pfk7_data$prev_580y[pfk7_data$nobs_c580y > numofobs & is.finite(pfk7_data$lrsmed_580y)]*pfk7_data$n[pfk7_data$nobs_c580y > numofobs & is.finite(pfk7_data$lrsmed_580y)]), pfk7_data$n[pfk7_data$nobs_c580y > numofobs & is.finite(pfk7_data$lrsmed_580y)])^2)

pfk7_data$wt_med_539t[pfk7_data$nobs_r539t > numofobs & is.finite(pfk7_data$lrsmed_539t)] <- 1/(se_ln_ratio_noZeros(round(pfk7_data$prev_539t[pfk7_data$nobs_r539t > numofobs & is.finite(pfk7_data$lrsmed_539t)]*pfk7_data$n[pfk7_data$nobs_r539t > numofobs & is.finite(pfk7_data$lrsmed_539t)]), pfk7_data$n[pfk7_data$nobs_r539t > numofobs & is.finite(pfk7_data$lrsmed_539t)])^2)

pfk7_data$wt_med_all_kelch[pfk7_data$nobs_kelch > numofobs & is.finite(pfk7_data$lrsmed_all_kelch)] <- 1/(se_ln_ratio_noZeros(round(pfk7_data$prev_all_kelch[pfk7_data$nobs_kelch > numofobs & is.finite(pfk7_data$lrsmed_all_kelch)]*pfk7_data$n[pfk7_data$nobs_kelch > numofobs & is.finite(pfk7_data$lrsmed_all_kelch)]), pfk7_data$n[pfk7_data$nobs_kelch > numofobs & is.finite(pfk7_data$lrsmed_all_kelch)])^2)

write.table(pfk7_data,"analysis/data/data-derived/pfk7_data.txt")
