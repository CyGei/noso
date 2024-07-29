load_helpers()
load_paper()
epicurve()

future::plan("multisession", workers = availableCores() - 2)
trees <- furrr::future_map(seq_along(cutoff_dates), function(x) {
  linelist_cut <- linelist %>%
    filter(onset <= cutoff_dates[x])

  out[[x]] %>%
    filter(step > 500) %>%
    identify(ids = linelist_cut$case_id) %>%
    get_trees(
      ids = linelist_cut$case_id,
      group = linelist_cut$group,
      date = linelist_cut$onset,
      kappa = TRUE
    )
})




grouped_kappa_df <- furrr::future_map(seq_along(trees), function(i) {
  trees <- trees[[i]]
  tab <- expand_grid(
    from_group = unique(linelist$group),
    to_group = unique(linelist$group)
  )

  data <- map(trees, ~ {
    g_kappa <- .x %>%
      drop_na(kappa) %>%
      group_by(from_group, to_group, kappa) %>%
      summarise(n = n(), .groups = "drop") %>%
      # group_by(from_group, to_group) %>%
      # summarise(prop_kappa_over1 = sum(n[kappa > 1]) / sum(n),
      #           .groups = "drop")
      group_by(from_group) %>%
      mutate(total_from = sum(n)) %>%
      group_by(from_group, to_group) %>%
      summarise(
        prop_kappa_over1 = sum(n[kappa > 1]) / first(total_from),
        .groups = "drop"
      )
    left_join(tab, g_kappa, by = c("from_group", "to_group"))
  }) %>% bind_rows(.id = "tree")

  data %>%
    drop_na(prop_kappa_over1) %>%
    group_by(from_group, to_group) %>%
    summarise(
      mean = mean(prop_kappa_over1),
      lwr = quantile(prop_kappa_over1, 0.025),
      upr = quantile(prop_kappa_over1, 0.975),
      .groups = "drop"
    )
})

#
# bind_rows(grouped_kappa_df, .id = "cutoff_date") %>%
#   mutate(from_label = as.factor("from"), to_label = as.factor("to")) %>%
#   mutate(cutoff_date = cutoff_dates[as.integer(cutoff_date)]) %>%
#   ggplot(aes(x = cutoff_date, y = mean, color = from_group)) +
#   # facet_grid(rows = vars(from_group),
#   #            cols = vars(to_group),
#   #            switch = "y")+
#   ggh4x::facet_nested(
#     rows = vars(from_label, from_group),
#     cols = vars(to_label, to_group),
#     switch = "y"
#   ) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lwr, ymax = upr)) +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   theme_noso(day_break = 4) +
#   labs(x = "", y = expression(paste(
#     "Proportion of indirect transmissions: ", kappa > 1
#   ))) +
#   theme(legend.position = "none")
#





overall_kappa_df <- furrr::future_map(seq_along(trees), function(i) {
  trees <- trees[[i]]

  data <- map(trees, ~ {
    .x %>%
      drop_na(kappa) %>%
      group_by(from_group, kappa) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(from_group) %>%
      summarise(prop_kappa_over1 = sum(n[kappa > 1]) / sum(n))
  }) %>% bind_rows(.id = "tree")

  data %>%
    drop_na(prop_kappa_over1) %>%
    group_by(from_group) %>%
    summarise(
      mean = mean(prop_kappa_over1),
      lwr = quantile(prop_kappa_over1, 0.025),
      upr = quantile(prop_kappa_over1, 0.975),
      .groups = "drop"
    ) %>%
    mutate(to_group = "overall")
})


kappa_df <-
  bind_rows(
    bind_rows(overall_kappa_df, .id = "cutoff_date"),
    bind_rows(grouped_kappa_df, .id = "cutoff_date")
  ) %>%
  mutate(
    cutoff_date = cutoff_dates[as.integer(cutoff_date)],
    from_label = as.factor("from"), to_label = as.factor("to"),
    to_group = factor(to_group, levels = c("overall", "hcw", "patient"))
  )

kappa_df %>%
  ggplot(aes(x = cutoff_date, y = mean, color = from_group)) +
  ggh4x::facet_nested(
    rows = vars(from_label, from_group),
    cols = vars(to_label, to_group),
    switch = "y"
  ) +
  geom_point() +
  geom_errorbar(aes(ymin = lwr, ymax = upr)) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    position = "right"
  ) +
  theme_noso(day_break = 4) +
  labs(x = "", y = expression(paste(
    "Proportion of indirect transmissions: ", kappa > 1
  ))) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = NA)
  )
