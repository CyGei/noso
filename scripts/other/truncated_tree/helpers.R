# # Chi-Square Test ---------------------------------------------------------
# get_chisq <- function(tabx, taby) {
#   tab <- inner_join(tabx, taby, by = c("from", "to")) %>%
#     filter(!(n.x == 0 & n.y == 0)) %>%
#     select(starts_with("n"))
#
#   # Calculate expected values
#   total_sum <- sum(tab)
#   row_sums <- rowSums(tab)
#   col_sums <- colSums(tab)
#   expected_values <- outer(row_sums, col_sums, "*") / total_sum
#   #the expected values represent the values that we would expect to see in each cell of the table
#   #if the null hypothesis of independence between the two variables (rows and columns) is true.
#   #(Row Total * Column Total) / Grand Total
#
#   # Determine whether to simulate p-values
#   simulate_p <- any(expected_values < 5)
#
#   # Check if there is only one non-zero column
#   only_one_non_zero_col <- sum(col_sums > 0) == 1
#
#   # Perform appropriate test
#   if (nrow(tab) < 2) {
#     p_value <- NA
#   } else if (simulate_p && !only_one_non_zero_col) {
#     # Perform Fisher's exact test with simulated p-value
#     res <- fisher.test(tab, simulate.p.value = TRUE)
#     p_value <- res$p.value
#   } else if (!simulate_p && !only_one_non_zero_col) {
#     # Perform Chi-squared test
#     res <- chisq.test(tab)
#     p_value <- res$p.value
#   } else if (only_one_non_zero_col) {
#     # Perform Fisher's exact test without simulating p-value
#     res <- fisher.test(tab)
#     p_value <- res$p.value
#   } else {
#     warning("Invalid input data. Cannot perform any test.")
#     p_value <- NA
#   }
#
#   return(p_value)
# }



get_chisq<- function(tabx, taby) {

  tab <- inner_join(tabx, taby, by = c("from", "to")) %>%
    filter(!(n.x == 0 & n.y == 0)) %>%
    select(starts_with("n"))

  # Check input
  if (!is_tibble(tab) || ncol(tab) != 2) {
    stop("Input must be a tibble with exactly 2 columns.")
  }

  # Calculate key metrics
  total_obs <- sum(tab)
  expected <- outer(rowSums(tab), colSums(tab)) / total_obs
  min_expected <- min(expected)
  prop_small_expected <- sum(expected < 5) / length(expected) #Proportion of cells with expected frequency < 5

  # Condition 1: No data
  if (nrow(tab) <= 1 || total_obs == 0) {
    return(NA_real_)
  }

  # Condition 2: Use Fisher's exact test
  use_fisher <- total_obs < 30 || min_expected < 1 || prop_small_expected > 0.2

  if (use_fisher) {
    tryCatch({
      # Try regular Fisher's test first
      fisher_test <- fisher.test(as.matrix(tab))
      return(fisher_test$p.value)
    }, error = function(e) {
      # If regular Fisher's test fails, use simulation
      fisher_test <- fisher.test(as.matrix(tab), simulate.p.value = TRUE, B = 10000)
      return(fisher_test$p.value)
    })
  }

  # Condition 3: Use Chi-square test
  # Determine if simulation is needed for Chi-square test
  use_chisq_simulation <- total_obs >= 1000 || nrow(tab) > 100

  if (use_chisq_simulation) {
    chisq_test <- chisq.test(as.matrix(tab), simulate.p.value = TRUE, B = 10000)
  } else {
    chisq_test <- chisq.test(as.matrix(tab), simulate.p.value = FALSE)
  }

  return(chisq_test$p.value)
}
