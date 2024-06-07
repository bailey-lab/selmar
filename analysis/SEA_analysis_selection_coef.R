analysis_selection_coef <- function(effect_size_bay){

  selec_ceof_stats <- data.frame()
  min_580Y <- min(effect_size_bay %>% filter(mutation == "580Y" & District != "Southeast Asia") %>% pull(estimate))
  max_580Y <- max(effect_size_bay %>% filter(mutation == "580Y" & District != "Southeast Asia") %>% pull(estimate))
  stats_580Y <- cbind("580Y",
                     round(min_580Y, 3),
                     effect_size_bay[which(effect_size_bay$estimate == min_580Y),][5],
                     round(max_580Y, 3),
                     effect_size_bay[which(effect_size_bay$estimate == max_580Y),][5])

  min_539T <- min(effect_size_bay %>% filter(mutation == "539T" & District != "Southeast Asia") %>% pull(estimate))
  max_539T <- max(effect_size_bay %>% filter(mutation == "539T" & District != "Southeast Asia") %>% pull(estimate))
  stats_539T <- cbind("539T",
                     round(min_539T, 3),
                     effect_size_bay[which(effect_size_bay$estimate == min_539T),][5],
                     round(max_539T, 3),
                     effect_size_bay[which(effect_size_bay$estimate == max_539T),][5])

  min_K13 <- min(effect_size_bay %>% filter(mutation == "K13" & District != "Southeast Asia") %>% pull(estimate))
  max_K13 <- max(effect_size_bay %>% filter(mutation == "K13" & District != "Southeast Asia") %>% pull(estimate))
  stats_K13 <- cbind("K13",
                     round(min_K13, 3),
                     effect_size_bay[which(effect_size_bay$estimate == min_K13),][5],
                     round(max_K13, 3),
                     effect_size_bay[which(effect_size_bay$estimate == max_K13),][5])

  colnames(stats_580Y) <- c("Mutation", "Min", "District", "Max", "District")
  colnames(stats_539T) <- c("Mutation", "Min", "District", "Max", "District")
  colnames(stats_K13) <- c("Mutation", "Min", "District", "Max", "District")

  selec_ceof_stats <- rbind(stats_539T, stats_580Y, stats_K13)

  return(selec_ceof_stats)
}
