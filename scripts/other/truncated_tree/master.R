source(here::here("scripts/outbreaker", "helpers.R"))
paper <- c("JHI2021", "eLife2022")[2]

source(here::here("scripts", "helpers.R"))
load_libraries()
load_data(paper)
cutoff_dates

# Retrospective Analysis -------------------------------------------------
outA <- out[[length(out)]] %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  filter_alpha_by_kappa(1L) %>%
  get_trees(ids = linelist$case_id,
            group = linelist$group,
            date = linelist$onset)

trunc_ancestries <- lapply(cutoff_dates, function(cutoff_date) {
  trees_trunc <- bind_rows(outA, .id = "tree") %>%
    filter(
      (is.na(from_date) & to_date <= cutoff_date) |
        (
          !is.na(from_date) &
            from_date <= cutoff_date & to_date <= cutoff_date
        )
    )
  ances <- linktree:::ttable(from = trees_trunc$from,
                             to = trees_trunc$to,
                             level = linelist$case_id) %>%
    as_tibble() %>%
    filter(from != to)
  return(ances)
})



# Real-Time Analysis ------------------------------------------------------

RT_ancestries <- lapply(seq_along(cutoff_dates), function(cutoff_date) {
  linelistRT <- linelist %>%
    filter(onset <= cutoff_dates[cutoff_date])

  trees_RT <- out[[cutoff_date]] %>%
    filter(step > 500) %>%
    identify(ids = linelistRT$case_id) %>%
    filter_alpha_by_kappa(1L) %>%
    get_trees(
      ids = linelistRT$case_id,
      group = linelistRT$group,
      date = linelistRT$onset
    ) %>%
    bind_rows(.id = "tree")
  ances <- linktree:::ttable(from = trees_RT$from,
                             to = trees_RT$to,
                             level = linelist$case_id) %>%
    as_tibble() %>%
    filter(from != to)
  return(ances)
})


# Chi-Square Test ---------------------------------------------------------
get_chisq <- function(tabx, taby) {
  tab <- left_join(tabx, taby, by = c("from", "to")) %>%
    filter(!(n.x == 0 & n.y == 0)) %>%
    select(starts_with("n"))

  # Calculate expected values
  total_sum <- sum(tab)
  row_sums <- rowSums(tab)
  col_sums <- colSums(tab)
  expected_values <- outer(row_sums, col_sums, "*") / total_sum
  #the expected values represent the values that we would expect to see in each cell of the table
  #if the null hypothesis of independence between the two variables (rows and columns) is true.
  #(Row Total * Column Total) / Grand Total

   # Determine whether to simulate p-values
  simulate_p <- any(expected_values < 5)

  # Check if there is only one non-zero column
  only_one_non_zero_col <- sum(col_sums > 0) == 1

  # Perform appropriate test
  if (simulate_p && !only_one_non_zero_col) {
    # Perform Fisher's exact test with simulated p-value
    res <- fisher.test(tab, simulate.p.value = TRUE)
    p_value <- res$p.value
  } else if (!simulate_p && !only_one_non_zero_col) {
    # Perform Chi-squared test
    res <- chisq.test(tab)
    p_value <- res$p.value
  } else if (only_one_non_zero_col) {
    # Perform Fisher's exact test without simulating p-value
    res <- fisher.test(tab)
    p_value <- res$p.value
  } else {
    warning("Invalid input data. Cannot perform any test.")
    p_value <- NA
  }

  return(p_value)
}


plan(multisession, workers = availableCores() - 1)
chisq <- furrr::future_map(seq_along(cutoff_dates), function(i) {
  tibble(
    cutoff_date = cutoff_dates[i],
    p_value = get_chisq(tabx = RT_ancestries[[i]], taby = trunc_ancestries[[i]])
  )
}, .options = furrr_options(seed = TRUE)) %>%
  bind_rows()


# Plot --------------------------------------------------------------------

p_chisq <- ggplot(chisq, aes(x = cutoff_date, y = p_value)) +
  geom_point(size = 3,
             fill = "black") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "", y =  paste0("\u03C7\u00B2", " p-value")) +
  theme_noso(day_break = 4)


p_chisq_JHI2021 <- p_chisq
epicurve_JHI2021 <- epicurve()

p_chisq_eLife2022 <- p_chisq
epicurve_eLife2022 <- epicurve(day_break = 4)
grid <- cowplot::plot_grid(p_chisq_JHI2021,
                   p_chisq_eLife2022,
                   epicurve_JHI2021 + theme(legend.position = "none"),
                   epicurve_eLife2022 + theme(legend.position = "none"),
                   ncol = 2, align = "v")
cowplot::plot_grid(grid,
                   cowplot::get_legend(epicurve()),
                   nrow = 2, rel_heights = c(1, 0.075))
