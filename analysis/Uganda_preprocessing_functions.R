obtain_samples_low_coverage <- function(df_coverage_finite){
  #check if there are samples for which one MIP was not identified, but all other MIPs for mutations were found. We will include those samples for the overall K13 mutation calculation
  df_corrected <- data.frame() #initalize dataframe
  for (idx in unique(df_coverage_finite$id)){ #obtain unique id's
    df <- df_coverage_finite %>% filter(id == idx) #subset dataframe
    df$remove <- NA #remove will hold the mutations which we will not include for that sample in the individual mutation calculation as coverage is to low
    df$genotype_updated <- df$genotype #genotype_updated will hold the genotype of the updated
    if (sum(is.finite(df$genotype)) == 4 & sum(df$merge_coverage<=5) == 1){ #check if there is a sample where all mutations except one were MIPed
      df$remove[df$merge_coverage<=5] = df$Locus[df$merge_coverage<=5] #add mutation to remove column
      if (df$merge_ref_umi_count[df$merge_coverage<=5] <= 2 & df$merge_alt_umi_count[df$merge_coverage<=5] > 0){ #check if mutant
        df$genotype_updated[df$merge_coverage<=5] = 2
      }
      else if (df$merge_ref_umi_count[df$merge_coverage<=5] > 2 & df$merge_alt_umi_count[df$merge_coverage<=5] > 2){ #check if mixed
        df$genotype_updated[df$merge_coverage<=5] = 1
      }
      else{ #check if WT
        df$genotype_updated[df$merge_coverage<=5] = 0
      }
    }
    df_corrected  <- rbind(df_corrected, df) #append the subset
  }
  return(df_corrected)
}

calculate_ratio <- function(df_corrected, mut){
  #obtain K13 ratio of umi_alt_count to coverage using the MIPs which were not recorded due to too low coverage
  df_ratio <- df_corrected %>%
    filter(is.finite(genotype_updated)) %>%
    mutate(ratio = ifelse(merge_alt_umi_count/merge_coverage == "NaN", NA, merge_alt_umi_count/merge_coverage)) %>%
    select(site, id, ratio, Locus) %>%
    pivot_wider(names_from = Locus, values_from = ratio)

  #sum each row to obtain K13 WSAF
  df_ratio$K13 <- rowSums(df_ratio[, mut],na.rm = TRUE)
  df_ratio$WT <- 1-df_ratio$K13
  return(df_ratio)
}

add_samples_remove <- function(df_ratio, df_corrected){
  # add column to ratio dataframe which notes if we should add the sample to the individual WSAF calculation
  df_ratio$remove <- NA
  #remove all rows where MIPs are correct
  df_rem <- df_corrected %>%
    filter(!is.na(remove))
  for (r in 1:length(df_rem$id)){
    idx = df_rem$id[r]
    df_ratio$remove[df_ratio$id == idx] <- df_rem$remove[r]
  }
  return(df_ratio)
}

calc_freq_individual_mutations <- function(df_ratio, mutations){
  #calculate allele freq for each K13 mutation removing the samples that have too low coverage
  #by checking if the mutations name is in the remove column
  df_freq_site <- data.frame(matrix(ncol = length(mutations), nrow = 0))
  sites <- c()
  years <- c()
  for (y in unique(df_ratio$year)){
    for (s in unique(df_ratio$site)){ #iterate through each site
      df_site <- df_ratio %>% filter(site == s, year == y) #subsample dataset for site s
      freq_site_mutation <- c() #initialize the vector that will hold the freq for that site s
      for (m in mutations){ #iterate through all mutations
        filtered_df <- df_site %>% filter(remove != m | is.na(remove)) #remove all samples that have too low coverage
        col_idx = which(colnames(filtered_df) == m) #obtain column index for the mutation
        freq = sum(filtered_df[,col_idx], na.rm=TRUE)/(sum(filtered_df[,col_idx] + filtered_df$WT, na.rm=TRUE)) #calcualte freq = mut/(mut+WT)
        freq_site_mutation <- append(freq_site_mutation, freq) #append freq for m and s to vector
      }
      df_freq_site <- rbind(df_freq_site, freq_site_mutation)
    }
    sites <- append(sites, unique(df_ratio$site))
    years <- append(years, rep(y, length(unique(df_ratio$site))))
  }
  colnames(df_freq_site) <- mutations
  df_freq_site <- df_freq_site %>%
    mutate(site = sites,
           year = years) %>%
    select(site, year, A675V, R561H, C469Y, C469F, P441L)
  return(df_freq_site)
}

calc_freq_K13 <- function(df_ratio){
  #calculate allele freq for K13 mutation
  freq_K13 <- df_ratio %>%
    group_by(site, year) %>%
    summarise(n=n(),
              K13 = sum(PfK13, na.rm=TRUE)/(sum(PfK13 + WT, na.rm=TRUE)))
  freq_K13 <- as.data.frame(freq_K13)
  return(freq_K13)
}

calc_freq_individual_mutations_mean_year <- function(df_ratio, mutations){
  #calculate allele freq for each K13 mutation removing the samples that have too low coverage
  #by checking if the mutations name is in the remove column
  df_freq_year <- data.frame(matrix(ncol = length(mutations), nrow = 0))
  years <- c()
  for (y in unique(df_ratio$year)){
    df_year <- df_ratio %>% filter(year == y) #subsample dataset for site s
    freq_mutation <- c() #initialize the vector that will hold the freq for that site s
    for (m in mutations){ #iterate through all mutations
      filtered_df <- df_year %>% filter(remove != m | is.na(remove)) #remove all samples that have too low coverage
      col_idx = which(colnames(filtered_df) == m) #obtain column index for the mutation
      freq = mean(filtered_df[,col_idx], na.rm=TRUE)
      freq_mutation <- append(freq_mutation, freq) #append freq for m and s to vector
    }
    df_freq_year <- rbind(df_freq_year, freq_mutation)
  }
  years <- unique(df_ratio$year)
  colnames(df_freq_year) <- mutations
  df_freq_year <- df_freq_year %>%
    mutate(year = years) %>%
    select(year, A675V, R561H, C469Y, C469F, P441L)
  return(df_freq_year)
}

calc_freq_K13_mean_year <- function(df_ratio){
  #calculate allele freq for K13 mutation
  freq_K13 <- df_ratio %>%
    group_by(year) %>%
    summarise(n=n(),
              K13 = mean(PfK13, na.rm=TRUE))
  freq_K13 <- as.data.frame(freq_K13)
  return(freq_K13)
}

calc_freq_individual_mutations_mean <- function(df_ratio, mutations){
  #calculate allele freq for each K13 mutation removing the samples that have too low coverage
  #by checking if the mutations name is in the remove column
  df_freq_site <- data.frame(matrix(ncol = length(mutations), nrow = 0))
  sites <- c()
  years <- c()
  for (y in unique(df_ratio$year)){
    for (s in unique(df_ratio$site)){ #iterate through each site
      df_site <- df_ratio %>% filter(site == s, year == y) #subsample dataset for site s
      freq_site_mutation <- c() #initialize the vector that will hold the freq for that site s
      for (m in mutations){ #iterate through all mutations
        filtered_df <- df_site %>% filter(remove != m | is.na(remove)) #remove all samples that have too low coverage
        col_idx = which(colnames(filtered_df) == m) #obtain column index for the mutation
        freq = mean(filtered_df[,col_idx], na.rm=TRUE)
        freq_site_mutation <- append(freq_site_mutation, freq) #append freq for m and s to vector
      }
      df_freq_site <- rbind(df_freq_site, freq_site_mutation)
    }
    sites <- append(sites, unique(df_ratio$site))
    years <- append(years, rep(y, length(unique(df_ratio$site))))
  }
  colnames(df_freq_site) <- mutations
  df_freq_site <- df_freq_site %>%
    mutate(site = sites,
           year = years) %>%
    select(site, year, A675V, R561H, C469Y, C469F, P441L)
  return(df_freq_site)
}

calc_freq_K13_mean <- function(df_ratio){
  #calculate allele freq for K13 mutation
  freq_K13 <- df_ratio %>%
    group_by(site,year) %>%
    summarise(n=n(),
              K13 = mean(PfK13, na.rm=TRUE))
  freq_K13 <- as.data.frame(freq_K13)
  return(freq_K13)
}

calc_avg_ratio_individual_mutations <- function(df_ratio, df_corrected){
  #calculate allele freq for each K13 mutation removing the samples that have too low coverage
  #by checking if the mutations name is in the remove column
  df_mean_site <- data.frame(matrix(ncol = length(unique(df_corrected$Locus)), nrow = 0))
  for (s in unique(df_ratio$site)){
    df_site <- df_ratio %>% filter(site == s) #subsample dataset for site s
    df_site <- as.data.frame(df_site)
    mean_site_mutation <- c() #initialize the vector that will hold the freq for that site s
    for (m in unique(df_corrected$Locus)){
      filtered_df <- df_site %>%
        filter(remove != m | is.na(remove))  #remove all samples that have too low coverage
      col_idx = which(colnames(df_site) == m) #obtain column index for the mutation
      filtered_df <- filtered_df[filtered_df[,col_idx]>0 & !is.na(filtered_df[,col_idx]),]
      mean_ratio = mean(filtered_df[,col_idx], na.rm = TRUE) #calculate mean of individual muations with ratio>0
      mean_site_mutation <- append(mean_site_mutation, mean_ratio) #append mean_ratio for m and s to vector
    }
    df_mean_site <- rbind(df_mean_site, mean_site_mutation)
  }
  colnames(df_mean_site) <- unique(df_corrected$Locus)
  df_mean_site$site <- unique(df_ratio$site)
  df_mean_site <- df_mean_site[,append("site", unique(df_corrected$Locus))]
  return(df_mean_site)
}

remove_MIP_samples <- function(df_sanger){
  #remove samples which we already accounted for in MIPs using remove column
  rem_mut <- unique(df_sanger$remove)[!is.na(unique(df_sanger$remove))] #obtain mutations for which in some samples we already accounted for in the MIP calculation
  for (m in rem_mut){ #loop through mutations
    m_idx <- colnames(df_sanger)[colnames(df_sanger) == m] #obtain column index of mutation
    rem_idx <- which(df_sanger$remove == m) #obtain row indexes of mutations which need to be converted to NA
    df_sanger <- df_sanger[-rem_idx, ]  #make the sample and mutation NA, since we already accounted for it in MIP calculation
  }
  return(df_sanger)
}

sanger_calc_freq <- function(df_sanger, df_MIP_freq_avg, sanger_mutations, MIP_mutations){
  #obtain WSAF from MIP for sanger mutations
  df_MIP_freq_avg <- as.data.frame(df_MIP_freq_avg)
  df_sanger <- as.data.frame(df_sanger)
  for (m in sanger_mutations){
    if (m %in% MIP_mutations){
      m_colidx <- which(colnames(df_sanger)==m)

      #update WSAF for sanger mutant
      idx_mutant <- which(df_sanger[ ,m_colidx] == 2)
      #update WSAF for sanger mixed
      idx_mixed <- which(df_sanger[ ,m_colidx] == 1)

      # update mutant samples WSAF
      df_sanger[idx_mutant, m_colidx] <- 1
      for (r in idx_mixed){
        id <- df_sanger[r,2]
        s <- df_sanger[r,1]
        MIP_cidx <- which(colnames(df_MIP_freq_avg) == m)
        MIP_ridx <- which(df_MIP_freq_avg[,1] == s)
        if(df_MIP_freq_avg[MIP_ridx, MIP_cidx]=="NaN"){ #if we do not have WSAF from MIP
          df_sanger[r,m_colidx] <- 1 #assign WSAF of 1
        }
        else{ #if we have WSAF from MIP
          df_sanger[r, m_colidx] <- df_MIP_freq_avg[MIP_ridx, MIP_cidx] #assign MIP WSAF to sanger
        }
      }
    }
  }
  return(df_sanger)
}

updated_kelch13_ratio_greater_than_1 <- function(df){
  idx_K13_greater_1 <- which(df$K13 > 1)
  #check if the ratio for K13 > 1
  if (length(idx_K13_greater_1) != 0){
    for (idx in idx_K13_greater_1){
      df_K13_greater1 <- df[idx,]
      K13_freq <- df_K13_greater1$K13
      df_pivot <- df_K13_greater1 %>%
        select(id, A675V, R561H, C469Y, C469F, P441L, K13) %>%
        pivot_longer(cols=A675V:P441L,
                     names_to = "Locus",
                     values_to = "freq") %>%
        mutate(freq = freq/K13_freq) %>%
        pivot_wider(names_from = "Locus", values_from = "freq")
      df_pivot$K13 <- rowSums(df_pivot[,c("A675V", "R561H", "C469Y", "C469F", "P441L")], na.rm = TRUE)
      df_pivot$WT <- 1-df_pivot$K13
      df_pivot$site <- df_K13_greater1$site
      df_pivot$remove <- df_K13_greater1$remove
      df_pivot$type <- df_K13_greater1$type
      df_pivot <- df_pivot %>% select(site, id, type, A675V, R561H, C469Y, C469F, P441L, K13, WT, remove)
      df[idx,] <- df_pivot
    }
  }
  return(df)
}

add_MIP_genotype <- function(df_corrected, df_ratio_y, colnames_genotype){
  #add MIP genotype to dataframe
  df_MIP_genotype <- df_corrected %>%
    mutate(Locus = paste0(Locus, "_genotype")) %>%
    select(id, site, genotype_updated, Locus) %>%
    pivot_wider(names_from = Locus, values_from = genotype_updated)
  if (length(colnames(df_MIP_genotype)) != length(colnames_genotype)){
    missing_colnames <- colnames_genotype[!colnames_genotype %in% colnames(df_MIP_genotype)]
    for (m_missing in missing_colnames){
      df_MIP_genotype[, m_missing] <- NA
    }
    df_MIP_genotype <- df_MIP_genotype %>%
      select(site, id, A675V_genotype, R561H_genotype, C469Y_genotype, C469F_genotype, P441L_genotype)
  }
  id_NA_coverage <- df_MIP_genotype$id[!df_MIP_genotype$id %in% df_ratio_y$id]
  df_MIP_genotype <- df_MIP_genotype %>% filter(!(id %in% id_NA_coverage)) %>% select(!c(site, id))
  df_ratio_y_geno <- cbind(df_ratio_y, df_MIP_genotype)
  return(df_ratio_y_geno)
}

add_sanger_MIP_genotype <- function(df_corrected, df_ratio_y_geno, colnames_genotype){
  #add MIP and Sanger genotype to dataframe
  df_MIP_genotype <- df_corrected %>%
    mutate(Locus = paste0(Locus, "_genotype")) %>%
    select(id, site, genotype_updated, Locus) %>%
    pivot_wider(names_from = Locus, values_from = genotype_updated)
  if (length(colnames(df_MIP_genotype)) != length(colnames_genotype)){
    missing_colnames <- colnames_genotype[!colnames_genotype %in% colnames(df_MIP_genotype)]
    for (m_missing in missing_colnames){
      df_MIP_genotype[, m_missing] <- NA
    }
    df_MIP_genotype <- df_MIP_genotype %>%
      select(site, id, A675V_genotype, R561H_genotype, C469Y_genotype, C469F_genotype, P441L_genotype)
  }
  id_NA_coverage <- df_MIP_genotype$id[!df_MIP_genotype$id %in% df_MIP_ratio$id]
  df_MIP_genotype <- df_MIP_genotype %>% filter(!(id %in% id_NA_coverage)) %>% select(!site)

  #add Sanger genotype
  df_sanger_genotype <- df_sanger %>%
    rename(A675V_genotype = A675V,
           R561H_genotype = R561H,
           C469F_genotype = C469F,
           C469Y_genotype = C469Y,
           P441L_genotype = P441L) %>%
    select(!remove)
  if (length(colnames(df_sanger_genotype)) != length(colnames_genotype)){
    missing_colnames <- colnames_genotype[!colnames_genotype %in% colnames(df_MIP_genotype)]
    for (m_missing in missing_colnames){
      df_sanger_genotype[, m_missing] <- NA
    }
    df_sanger_genotype <- df_sanger_genotype %>%
      select(site, id, A675V_genotype, R561H_genotype, C469Y_genotype, C469F_genotype, P441L_genotype)
  }
  id_NA_coverage <- df_sanger_genotype$id[!df_sanger_genotype$id %in% df_sanger_ratio$id]
  df_sanger_genotype <- df_sanger_genotype %>% filter(!(id %in% id_NA_coverage)) %>% select(!c(site))

  df_genotype <- rbind(df_MIP_genotype, df_sanger_genotype)
  df_genotype <- df_genotype %>% select(!id)
  df_ratio_y_geno <- cbind(df_ratio_y_geno, df_genotype)

  return(df_ratio_y_geno)
}
