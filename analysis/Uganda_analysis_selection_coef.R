analysis_selection_coef <- function(effect_size_bay){
  selec_ceof_stats <- data.frame()
  min_675V <- min(effect_size_bay %>%
                    filter(mutation == "675V" & District != "Uganda") %>%
                    pull(estimate))
  max_675V <- max(effect_size_bay %>%
                    filter(mutation == "675V" & District != "Uganda") %>%
                    pull(estimate))
  stats_675V <- cbind("675V", min_675V, effect_size_bay[which(effect_size_bay$estimate == min_675V),][1], max_675V, effect_size_bay[which(effect_size_bay$estimate == max_675V),][1])

  min_469Y <- min(effect_size_bay %>%
                    filter(mutation == "469Y" & District != "Uganda") %>%
                    pull(estimate))
  max_469Y <- max(effect_size_bay %>%
                    filter(mutation == "469Y" & District != "Uganda") %>%
                    pull(estimate))
  stats_469Y <- cbind("469Y", min_469Y, effect_size_bay[which(effect_size_bay$estimate == min_469Y),][1], max_469Y, effect_size_bay[which(effect_size_bay$estimate == max_469Y),][1])

  min_469F <- min(effect_size_bay %>%
                    filter(mutation == "469F" & District != "Uganda") %>%
                    pull(estimate))
  max_469F <- max(effect_size_bay %>%
                    filter(mutation == "469F" & District != "Uganda") %>%
                    pull(estimate))
  stats_469F <- cbind("469F", min_469F, effect_size_bay[which(effect_size_bay$estimate == min_469F),][1], max_469F, effect_size_bay[which(effect_size_bay$estimate == max_469F),][1])

  min_441L <- min(effect_size_bay %>%
                    filter(mutation == "441L" & District != "Uganda") %>%
                    pull(estimate))
  max_441L <- max(effect_size_bay %>%
                    filter(mutation == "441L" & District != "Uganda") %>%
                    pull(estimate))
  stats_441L <- cbind("441L", min_441L, effect_size_bay[which(effect_size_bay$estimate == min_441L),][1], max_441L, effect_size_bay[which(effect_size_bay$estimate == max_441L),][1])

  min_PfK13 <- min(effect_size_bay %>%
                     filter(mutation == "K13" & District != "Uganda") %>%
                     pull(estimate))
  max_PfK13 <- max(effect_size_bay %>%
                     filter(mutation == "K13" & District != "Uganda") %>%
                     pull(estimate))
  stats_PfK13 <- cbind("K13", min_PfK13, effect_size_bay[which(effect_size_bay$estimate == min_PfK13),][1], max_PfK13, effect_size_bay[which(effect_size_bay$estimate == max_PfK13),][1])

  colnames(stats_675V) <- c("Mutation", "Min Select Ceof", "District", "Max Select Ceof", "District")
  colnames(stats_469Y) <- c("Mutation", "Min Select Ceof", "District", "Max Select Ceof", "District")
  colnames(stats_469F) <- c("Mutation", "Min Select Ceof", "District", "Max Select Ceof", "District")
  colnames(stats_441L) <- c("Mutation", "Min Select Ceof", "District", "Max Select Ceof", "District")
  colnames(stats_PfK13) <- c("Mutation", "Min Select Ceof", "District", "Max Select Ceof", "District")
  selec_coef_stats <- rbind(stats_675V, stats_469Y, stats_469F, stats_441L, stats_PfK13)
}
