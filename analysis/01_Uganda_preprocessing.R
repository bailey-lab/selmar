library(lme4)
library(nlme)
library(dplyr)
library(tabulizer)
library(tidyr)

#load data from raw-data
out2 <- tabulizer::extract_tables("analysis/data/data-raw/Conrad2023_NEJM_submission.pdf", pages = 1, guess = TRUE, output = "data.frame")

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
  mutate(med = x / n) %>%
  rowwise() %>%
  na.omit %>%
  mutate(lrsmed = log(med/(1-med))) %>%
  mutate(x_org = x) %>%
  rename(year = name) %>%
  mutate(year = as.numeric(year)) %>%
  group_by(District, Locus) %>%
  mutate(min_year = min(year[x>0]),
         nobs= sum(x > 0)) %>%
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

#COMBINED FOR OVERALL UGANDAN SELECTION COEF
df_comb <- df %>% group_by(District, year) %>% mutate(comb = sum(x)) %>% ungroup

#Mixed genotypes
df_comb <- as.data.frame(df_comb)
df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2019"] <-  df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2019"] - 1

df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2020"] <-  df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2020"] - 3
df_comb$comb[df_comb$District == "Kaabong" & df_comb$year == "2020"] <-  df_comb$comb[df_comb$District == "Kaabong" & df_comb$year == "2020"] - 1
df_comb$comb[df_comb$District == "Kole" & df_comb$year == "2020"] <-  df_comb$comb[df_comb$District == "Kole" & df_comb$year == "2020"] - 1
df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2020"] <-  df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2020"] - 5

df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2021"] <-  df_comb$comb[df_comb$District == "Agago" & df_comb$year == "2021"] - 1
df_comb$comb[df_comb$District == "Katakwi" & df_comb$year == "2021"] <-  df_comb$comb[df_comb$District == "Katakwi" & df_comb$year == "2021"] - 3
df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2021"] <-  df_comb$comb[df_comb$District == "Lamwo" & df_comb$year == "2021"] - 1
df_comb$comb[df_comb$District == "Kaabong" & df_comb$year == "2021"] <-  df_comb$comb[df_comb$District == "Kaabong" & df_comb$year == "2021"] - 1

df_comb2 <- df_comb %>% group_by(District, year) %>% mutate(n = min(n)) %>% distinct(comb, n) %>% group_by(District) %>%
  mutate(nobs = sum(comb>0)) %>% ungroup %>%mutate(med = comb/n, lrsmed = log(med/(1-med)))

df_comb2$lrsmed[df_comb2$lrsmed == -Inf] <- NA

df_comb <- df_comb2 %>%
  group_by(District) %>%
  mutate(min_year = min(year[comb>0])) %>%
  ungroup %>%
  mutate(adj_year = year - min_year)

df_comb$wt_med <- NA
df_comb$wt_med[df_comb$nobs > numofobs & is.finite(df_comb$lrsmed)] <- 1/(se_ln_ratio_noZeros(round(df_comb$med[df_comb$nobs > numofobs & is.finite(df_comb$lrsmed)]*df_comb$n[df_comb$nobs > numofobs & is.finite(df_comb$lrsmed)]), df_comb$n[df_comb$nobs > numofobs & is.finite(df_comb$lrsmed)])^2)

write.table(df_comb, "analysis/data/data-derived/Ugandan_comb_data.txt")
write.table(df, "analysis/data/data-derived/Ugandan_data.txt")
