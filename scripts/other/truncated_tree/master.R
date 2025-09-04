load_helpers()
load_libraries()
load_paper()

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
  geom_point(size = 3, fill = "black") +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_log10(
    breaks = c(1e-04, 1e-03, 0.05, 1e+00),
    # Set the breaks at 1e-04, 0.05, and 1e+00
    labels = c(
      expression(10 ^ -4),
      expression(10 ^ -3),
      "0.05",
      expression(10 ^ 0)
    )  # Custom labels if needed
  ) +
  labs(x = "", y =  paste0("\u03C7\u00B2", " p-value")) +
  theme_noso(day_break = 2)

epic <-
  epicurve() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.key.size = unit(0.4, "cm"),                  # Reduce the size of the legend keys
    legend.margin = margin(t = 0, r = 1, b = 0, l = 1),
    legend.spacing.y = unit(0, "cm")
  ) +
  labs(x = "")+
  guides(fill = guide_legend(title = NULL),
         color = guide_legend(title = NULL),
         linetype = guide_legend(title = NULL))

cowplot::plot_grid(
  epic,
  NULL,
  p_chisq + labs(x = "Onset"),
  ncol = 1,
  align = "v",
  rel_heights = c(1, -0.06, 0.45),
  labels = c("A", "", "B"),
  label_x = c(0, NA, 0),
  label_y = c(1, NA, 1.25)
)
ggsave(
  filename = "figs/chisq2.png",
  width = 8,
  height = 4.5,
  dpi = 500
)
