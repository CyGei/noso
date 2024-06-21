# The purpose of this script is to estimate the reproduction number
# using transmission chain data.
# The effective reproductive number describes how many new infections an individual
# infected  causes on average in a population which is subject to a certain degree of immunity and intervention measures.
# This is not applicable in real-time.
load_helpers()
load_paper()
epicurve()

#Final outbreaker2 run
trees <- out[[length(out)]] %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  get_trees(
    ids = linelist$case_id,
    group = linelist$group,
    date = linelist$onset,
    kappa = TRUE
  )


# Re: estimating R using 5/7 days time windows ------------------------------
window <- ifelse(paper == "JHI2021", 5, 7)
start <- min(linelist$onset)
end <- max(linelist$onset)
cutoff_breaks <- c(seq.Date(start, end, by = window), end) %>% unique()

plan(multisession, workers = length(cutoff_breaks))
Re_df_group <- future_map_dfr(seq_len(length(cutoff_breaks) - 1),
                              ~ process_window(.x, trees, cutoff_breaks, TRUE))
Re_df_global <- future_map_dfr(seq_len(length(cutoff_breaks) - 1),
                               ~ process_window(.x, trees, cutoff_breaks, FALSE))


Re_df <- bind_rows(Re_df_group, Re_df_global) %>%
  mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)

Re_summary <- Re_df %>%
  group_by(from_group, window_median, window_start, window_end) %>%
  summarise(
    mean = mean(Re),
    lwr = quantile(Re, 0.025),
    upr = quantile(Re, 0.975),
    .groups = "drop"
  )

p_Re <- Re_df %>%
  ggplot(aes(
    x = window_median,
    y = Re,
    fill = from_group,
    group = interaction(from_group, window_median)
  )) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbarh(
    data = Re_summary %>% filter(from_group == "Global"),
    aes(
      y = mean,
      xmin = window_start - 0.5,
      xmax = window_end + 0.5
    ),
    height = 0.1,
    color = "black"
  ) +
  geom_violin(
    col = NA,
    adjust = 1.8,
    alpha = 0.7,
    scale = "width",
    position = position_dodge(width = 3)
  ) +
  geom_pointrange(
    data = Re_summary,
    aes(
      x = window_median,
      y = mean,
      ymin = lwr,
      ymax = upr,
    ),
    position = position_dodge(width = 3),
    pch = 21,
    stroke = 0.5,
    size = 0.75
  ) +
  scale_y_continuous(breaks = seq(0, 10, 0.5)) +
  theme_noso(day_break = window, date = TRUE) +
  coord_cartesian(ylim = c(0, 3.5)) +
  scale_fill_manual(
    values = c(
      Global = "#767676",
      hcw = "#FF9300",
      patient = "#d800d1"
    ),
    labels = c(
      "Global" = "Global",
      "hcw" = "HCW",
      "patient" = "Patient"
    )
  ) +
  labs(x = "", y = expression(R["e"]), fill = "")

p_main <-
  cowplot::plot_grid(
    epicurve(day_break = window) +
      theme(legend.position = "none") +
      labs(x = ""),
    NULL,
    p_Re + labs(x = "Onset") + theme(legend.position = "none"),
    ncol = 1,
    rel_heights = c(1, -0.075, 1),
    align = "v",
    labels = c("A", "", "B")
  )

#add legends
p_legends <-
  cowplot::plot_grid(
  peak_legend(),
  NULL,
  cowplot::get_plot_component(p_Re, 'guide-box-bottom', return_all = TRUE),
  nrow = 1,
  rel_widths = c(1, -0.75, 1)
)

cowplot::plot_grid(
  p_main,
  p_legends,
  nrow = 2,
  rel_heights = c(1, 0.1)
)
