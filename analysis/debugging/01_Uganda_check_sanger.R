#checking the number of Sanger samples
library(dplyr)
library(tidyr)
PRX2016_1 <- read.csv("analysis/data/data-raw/PRX-00.csv")
PRX2016_2 <- read.csv("analysis/data/data-raw/PRX-01.csv")
PRX2017 <- read.csv("analysis/data/data-raw/PRX-02.csv")
PRX2018 <- read.csv("analysis/data/data-raw/PRX-03.csv")
PRX2019 <- read.csv("analysis/data/data-raw/PRX-04.csv")
PRX2020 <- read.csv("analysis/data/data-raw/PRX-05.csv")
PRX2021 <- read.csv("analysis/data/data-raw/PRX-06.csv")
PRX2022 <- read.csv("analysis/data/data-raw/PRX-07.csv")

#2017
df2017 <- PRX2017 %>%
  filter(!is.finite(genotype),
         is.finite(sanger_geno),
         is.finite(merge_coverage))
df2017_mut <- df2017 %>%
  filter(merge_alt_umi_count>0)
dim(df2017)
dim(df2017_mut)

#2018
df2018 <- PRX2018 %>%
  filter(!is.finite(genotype),
         is.finite(sanger_geno),
         is.finite(merge_coverage))
df2018_mut <- df2018 %>%
  filter(merge_alt_umi_count>0)
dim(df2018)
dim(df2018_mut)

#2019
df2019 <- PRX2019 %>%
  filter(!is.finite(genotype),
         is.finite(sanger_geno),
         is.finite(merge_coverage))
df2019_mut <- df2019 %>%
  filter(merge_alt_umi_count>0)
dim(df2019)
dim(df2019_mut)

#2020
df2020 <- PRX2020 %>%
  filter(!is.finite(genotype),
         is.finite(sanger_geno),
         is.finite(merge_coverage))
df2020_mut <- df2020 %>%
  filter(merge_alt_umi_count>0)
dim(df2020)
dim(df2020_mut)

#2021
df2021 <- PRX2021 %>%
  filter(!is.finite(genotype),
         is.finite(sanger_geno),
         is.finite(merge_coverage))
df2021_mut <- df2021 %>%
  filter(merge_alt_umi_count>0)
dim(df2021)
dim(df2021_mut)

#2022
df2022 <- PRX2022 %>%
  filter(!is.finite(genotype),
         is.finite(sanger_geno),
         is.finite(merge_coverage))
df2022_mut <- df2022 %>%
  filter(merge_alt_umi_count>0)
dim(df2022)
dim(df2022_mut)
