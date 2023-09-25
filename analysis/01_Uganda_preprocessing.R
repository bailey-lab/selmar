library(lme4)
library(nlme)
library(dplyr)
library(tabulizer)
library(tidyr)

#load data from raw-data
df <- read.csv("analysis/data/data-raw/conrad_raw.csv")

# get the data from the scraped table
df[df==""] <- NA
names(df) <- gsub("X","",names(df))
names(df) <- gsub("\\.\\..*","",names(df))

df <- df %>% fill(District, .direction = "down")
df$Trend.p.value <- NULL

# tidy up the \ style
for(i in 1:7){
  df[,i+2] <- gsub("\\(.*\\)", "", df[,i+2])
}

# explicitly turn - into 0/0, which is what it is.
df[df=="-"] <- "0/0"

# and turn into helpful table
df <- df %>% pivot_longer(cols = `2016`:`2022`) %>%
  mutate(n = as.numeric(gsub("(.*)/(.*)","\\2", value))) %>%
  mutate(x = as.numeric(gsub("(.*)/(.*)","\\1", value))) %>%
  mutate(med = x / n) %>%
  rowwise() %>%
  #na.omit %>%
  mutate(lrsmed = log(med/(1-med))) %>%
  mutate(x_org = x) %>%
  rename(year = name) %>%
  mutate(year = as.numeric(year)) %>%
  group_by(District, Locus) %>%
  mutate(min_year = min(year[x>0], na.rm = TRUE),
         nobs= sum(x > 0, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(adj_year = year - min_year)

# log ratio function
se_ln_ratio_noZeros<-function(x,N) {
  inds<-which(x==0)
  x[inds]<-0.5
  inds<-which(x==N)
  x[inds]<-x[inds]-0.5
  se<- sqrt(1/x + 1/(N-x))
  return(se)
}

numofobs <- 2
# remove the data that is Inf here due to 0s - this will obscure it from mdel fitting
df$lrsmed[df$lrsmed == -Inf] <- NA

df$wt_med <- NA
df$wt_med[df$nobs > numofobs & is.finite(df$lrsmed)] <- 1/(se_ln_ratio_noZeros(round(df$med[df$nobs > numofobs & is.finite(df$lrsmed)]*df$n[df$nobs > numofobs & is.finite(df$lrsmed)]), df$n[df$nobs > numofobs & is.finite(df$lrsmed)])^2)

write.table(df, "analysis/data/data-derived/Ugandan_data.txt")

## Table for supplement
df_supplement <- df %>% select(-c(value, lrsmed, x_org, wt_med)) %>%
  mutate(med = round(x / n,2))
df_supplement$adj_year[is.infinite(df_supplement$adj_year)] <- NA
df_supplement$min_year[is.infinite(df_supplement$min_year)] <- NA
write.csv(df_supplement, "analysis/data/data-derived/Ugandan_data_supplement.csv")


# Getting all mutations combined for a District
df_comb <- df %>% group_by(District, year) %>%
  mutate(n = max(n),
         x = round(med*n),
         comb = sum(x),
         comb = ifelse(is.na(comb), 0, comb)) %>%
  distinct(District, year, .keep_all = TRUE) %>%
  ungroup %>%
  select(District, year, n, comb)

# Mixed genotypes - TODO: Ask Melissa about 2022
df_comb <- as.data.frame(df_comb)

df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2019"] <-  df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2019"] - 1

df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2020"] <-  df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2020"] - 3
df_comb$comb[df_comb$District == "Kaabong" & df_comb$year == "2020"] <-  df_comb$comb[df_comb$District == "Kaabong" & df_comb$year == "2020"] - 1
df_comb$comb[df_comb$District == "Kole" & df_comb$year == "2020"] <-  df_comb$comb[df_comb$District == "Kole" & df_comb$year == "2020"] - 3
df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2020"] <-  df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2020"] - 5

df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2021"] <-  df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2021"] - 1
df_comb$comb[df_comb$District == "Katakwi" & df_comb$year == "2021"] <-  df_comb$comb[df_comb$District == "Katakwi" & df_comb$year == "2021"] - 3
df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2021"] <-  df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2021"] - 1

df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2022"] <-  df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2022"] - 2
df_comb$comb[df_comb$District == "Kapchorwa" & df_comb$year == "2022"] <-  df_comb$comb[df_comb$District == "Kapchorwa" & df_comb$year == "2022"] - 1
df_comb$comb[df_comb$District == "Kaabong" & df_comb$year == "2022"] <-  df_comb$comb[df_comb$District == "Kaabong" & df_comb$year == "2022"] - 1
df_comb$comb[df_comb$District == "Katakwi" & df_comb$year == "2022"] <-  df_comb$comb[df_comb$District == "Katakwi" & df_comb$year == "2022"] - 2
df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2022"] <-  df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2022"] - 1


# NOTE: calculating the n for the combination is not obvious.
# previously we had min(n) here. I don't think this is right as
# we are summing all the mutations observed, which often these are found
# in the loci that have the highest n. So let's use the maximum and then
# adjust expected mutations based on prevalence
df_comb2 <- df_comb %>%
  group_by(District) %>%
  mutate(nobs = sum(comb>0, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(med = comb/n,
         lrsmed = log(med/(1-med)))

df_comb2$lrsmed[df_comb2$lrsmed == -Inf] <- NA

df_comb <- df_comb2 %>%
  group_by(District) %>%
  mutate(min_year = min(year[comb>0])) %>%
  ungroup %>%
  mutate(adj_year = year - min_year)

df_comb$wt_med <- NA
df_comb$wt_med[df_comb$nobs > numofobs & is.finite(df_comb$lrsmed)] <- 1/(se_ln_ratio_noZeros(round(df_comb$med[df_comb$nobs > numofobs & is.finite(df_comb$lrsmed)]*df_comb$n[df_comb$nobs > numofobs & is.finite(df_comb$lrsmed)]), df_comb$n[df_comb$nobs > numofobs & is.finite(df_comb$lrsmed)])^2)

write.table(df_comb, "analysis/data/data-derived/Ugandan_comb_data.txt")
