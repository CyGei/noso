#From Paper1
MoRatio <- function(trees, to = NULL, linelist_group) {
  require(dplyr)
  require(tidyr)

  n_pairs <- vapply(trees, function(tree) sum(!is.na(tree[[1]])),
                    FUN.VALUE = integer(1)) %>% max()

  # Get number of infectors per group
  n_infectors <- trees %>%
    lapply(function(tree) {
      if (!is.null(to)) {
        tree <- tree %>%
          filter(to_group %in% to)
      }
      ttab <- table(tree[[1]], tree[[2]])
      rowSums(ttab)
    }) %>%
    bind_rows()

  # Observed frequency and probability of infectors
  observed <-
    n_infectors %>%
    tidyr::pivot_longer(everything(), names_to = "group", values_to = "n_infectors") %>%
    group_by(group, n_infectors) %>%
    summarize(count = n()) %>%
    mutate(
      prob_y = count / sum(count),
      freq_x = n_infectors / n_pairs,
      source = "observed"
    ) %>%
    select(group, source, freq_x, prob_y)

  # Expected frequency and probability of infectors
  f_expected <- table(linelist_group) / length(linelist_group)
  expected <- purrr::map_dfr(names(f_expected), ~ {
    group_name <- .x
    group_prob <- f_expected[[.x]]
    tibble(
      n_infectors = 0:n_pairs,
      prob_y = dbinom(0:n_pairs, size = n_pairs, prob = group_prob),
      freq_x = n_infectors / n_pairs,
      group = group_name,
      source = "expected"
    )
  }) %>%
    select(group, source, freq_x, prob_y)


  # Compute mean and 95% CI
  mean_ci <- bind_rows(observed, expected) %>% # can be used to plot
    group_by(source, group) %>%
    summarize(
      mean = weighted.mean(freq_x, prob_y),
      lwr = cNORM::weighted.quantile(x = freq_x, weights = prob_y, 0.025),
      upr = cNORM::weighted.quantile(x = freq_x, weights = prob_y, 0.975)
    )

  # Calculate observed:expected ratio with confidence intervals
  obs_exp_ratio <- mean_ci %>%
    tidyr::pivot_wider(names_from = source, values_from = c(mean, lwr, upr)) %>%
    mutate(
      ratio = mean_observed / mean_expected,
      ratio_lwr = lwr_observed / upr_expected,
      ratio_upr = upr_observed / lwr_expected
    ) %>%
    select(group, ratio, ratio_lwr, ratio_upr)

  return(list(mean_ci = mean_ci, obs_exp_ratio = obs_exp_ratio))
}

trees2 <- lapply(trees, function(tree) {
  tree <- tree %>%
    select(from_group, to_group)
  tree
})
linelist_group <- linelist$group
MoRatio(trees2, to = NULL, linelist_group)






# MoRatioCI <- function(trees, to = NULL, linelist_group) {
#   require(dplyr)
#   require(tidyr)
#
#   n_pairs <- nrow(trees[[length(trees)]])
#
#   # Get number of infectors per group
#   n_infectors <- trees %>%
#     lapply(function(tree) {
#       if (!is.null(to)) {
#         tree <- tree %>%
#           filter(to_group %in% to)
#       }
#       ttab <- table(tree[[1]], tree[[2]])
#       rowSums(ttab)
#     }) %>%
#     bind_rows()
#
#   # Observed frequency and probability of infectors
#   observed <-
#     n_infectors %>%
#     tidyr::pivot_longer(everything(), names_to = "group", values_to = "n_infectors") %>%
#     group_by(group, n_infectors) %>%
#     mutate(
#       observed_frequency = n_infectors  / n_pairs,
#     ) %>%
#     ungroup() %>%
#     group_by(group) %>%
#     summarize(
#       mean = mean(observed_frequency, na.rm = TRUE),
#       lwr = quantile(observed_frequency, 0.025, na.rm = TRUE),
#       upr = quantile(observed_frequency, 0.975, na.rm = TRUE),
#       source = "observed"
#     )
#
#
#   # Expected frequency and probability of infectors
#   f_expected <- table(linelist_group) / length(linelist_group)
#   expected <- purrr::map_dfr(names(f_expected), ~ {
#     group_name <- .x
#     group_prob <- f_expected[[.x]]
#     binomtest <- prop.test(x = group_prob*n_pairs , n = n_pairs,
#                            conf.level=0.95, correct = FALSE)
#     tibble(
#       group = group_name,
#       mean = binomtest$estimate[[1]],
#       lwr = binomtest$conf.int[1],
#       upr = binomtest$conf.int[2],
#       source = "expected"
#     )
#   })
#
#   # Compute mean and 95% CI
#   mean_ci <- bind_rows(observed, expected)
#
#   # Calculate observed:expected ratio with confidence intervals
#   obs_exp_ratio <- mean_ci %>%
#     tidyr::pivot_wider(names_from = source, values_from = c(mean, lwr, upr)) %>%
#     mutate(
#       ratio = mean_observed / mean_expected,
#       ratio_lwr = lwr_observed / upr_expected,
#       ratio_upr = upr_observed / lwr_expected
#     ) %>%
#     select(group, ratio, ratio_lwr, ratio_upr)
#
#   return(list(
#     mean_ci = mean_ci,
#     obs_exp_ratio = obs_exp_ratio
#   ))
# }


get_prop <- function (i, type_from_m, infector_type, id_info) {

  obs_table <- table(factor(type_from_m[, i], levels = unique(id_info$patient)))

  obs_prop <- prop.table(obs_table)[[infector_type]]

  exp_table_random <- table(sample(factor(id_info$patient, levels = unique(id_info$patient)), sum(obs_table), replace = TRUE))

  exp_prop <- prop.table(exp_table_random)[[infector_type]]

  c(obs = obs_prop, exp = ifelse(is.na(exp_prop), 0, exp_prop))

}



p_infector <- as.data.frame(t(vapply(1:ncol(type_from_m), function (i) get_prop(i, type_from_m, infector_type, id_info), numeric(2))))



## calculating ratios

p_infector$ratio <- p_infector$obs/p_infector$exp

mean_ratio <- mean(p_infector$ratio[is.finite(p_infector$ratio)])

median_ratio <- median(p_infector$ratio[is.finite(p_infector$ratio)])

ci_ratio <- quantile(p_infector$ratio[is.finite(p_infector$ratio)], c(0.025, 0.975))

p_ratio <- 1 - sum(p_infector$ratio[is.finite(p_infector$ratio)] > 1)/length(p_infector$ratio[is.finite(p_infector$ratio)])

p_ratio <- ifelse(p_ratio > 0.5, 1 - p_ratio, p_ratio)

p_ratio_2 <- 1 - sum(p_infector$ratio > 1)/length(p_infector$ratio[is.finite(p_infector$ratio)])

p_ratio_2 <- ifelse(p_ratio_2 > 0.5, 1 - p_ratio_2, p_ratio_2)

##


