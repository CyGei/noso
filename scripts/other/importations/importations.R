load_helpers()
load_paper()
epicurve()
out <- out[[length(out)]]
# Total Introductions
alpha_cols <- grep("alpha_", colnames(out))
fixed_imports <- sum(is.na(out[sample(1:nrow(out), 1), alpha_cols]))
fixed_imports

trees <- out %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  get_trees(
    ids = linelist$case_id,
    group = linelist$group,
    date = linelist$onset,
    kappa = TRUE,
    t_inf = TRUE
  )

# Distribution of  imports per group
total_imports <-
  lapply(trees, function(tree) {
  tree %>%
    filter(is.na(from)) %>%
    group_by(to_group) %>%
    summarise(n_imports = n(), .groups = "drop")
}) %>% bind_rows(.id = "tree") %>%
  group_by(to_group) %>%
  summarise(
    mean = mean(n_imports),
    lwr = quantile(n_imports, 0.025),
    upr = quantile(n_imports, 0.975),
    .groups = "drop"
  )
# as frequencies
total_imports %>%
  mutate(across(where(is.numeric), \(x) x /fixed_imports))

# Estimate the average number of imports per time window ------------
window <- ifelse(paper == "JHI2021", 5, 7)
start <- min(linelist$onset)
end <- max(linelist$onset)
cutoff_breaks <- c(seq.Date(start, end, by = window), end) %>% unique()

plan(multisession, workers = length(cutoff_breaks))
# Compute imports by group
imports_df_group <- future_map_dfr(seq_len(length(cutoff_breaks) - 1),
                                   ~ process_window(.x, trees, cutoff_breaks, TRUE, levels = unique(linelist$group)))

# Compute overall imports
imports_df_overall <- future_map_dfr(
  seq_len(length(cutoff_breaks) - 1),
  ~ process_window(.x, trees, cutoff_breaks, FALSE, levels = unique(linelist$group))
) %>%
  mutate(to_group = "global")

imports_df <- bind_rows(imports_df_group, imports_df_overall)
imports_summary <- imports_df %>%
  group_by(to_group, window_median, window_start, window_end) %>%
  summarise(
    mean = mean(n_imports),
    lwr = quantile(n_imports, 0.025),
    upr = quantile(n_imports, 0.975),
    .groups = "drop"
  )


p_imports <-
  imports_df %>%
  ggplot(aes(
    x = window_median,
    y = n_imports,
    fill = to_group,
    group = interaction(to_group, window_median)
  )) +
  geom_errorbarh(
    data = imports_summary %>%
      filter(to_group == "global"),
    aes(
      xmin = window_start - 0.5,
      xmax = window_end + 0.5,
      y = mean,
      height = 0
    ),
    #size = 0.3,
    height = 0.1
  ) +
  geom_violin(
    position = position_dodge(width = 3),
    kernel = "rectangular",
    adjust = 0.6,
    alpha = 0.7,
    col = NA,
    scale = "width"
  ) +
  geom_pointrange(
    data = imports_summary,
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
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  theme_noso(day_break = window, date = TRUE) +
  scale_fill_manual(
    values = c(
      global = "#767676",
      hcw = "#FF9300",
      patient = "#d800d1"
    ),
    labels = c(
      "global" = "Global",
      "hcw" = "HCW",
      "patient" = "Patient"
    )
  ) +
  coord_cartesian(ylim = c(0, 6)) +
  labs(x = "", y = "Importations", fill = "")

p_main <- cowplot::plot_grid(
  epicurve(day_break = window) + labs(x = "") + theme(legend.position = "none"),
  NULL,
  p_imports + theme(legend.position = "none"),
  ncol = 1,
  rel_heights = c(1.5, -0.12, 2),
  align = "v",
  labels = c("A", "", "B")
)
p_legends <-
  cowplot::plot_grid(
    peak_legend(),
    NULL,
    cowplot::get_plot_component(p_imports, 'guide-box-bottom', return_all = TRUE),
    nrow = 1,
    rel_widths = c(1, -0.75, 1)
  )

cowplot::plot_grid(p_main,
                   p_legends,
                   nrow = 2,
                   rel_heights = c(1, 0.1))




# #eLife2022: no introductions in the 1st 2 weeks!
# ids_to_check <- linelist %>%
#   arrange(onset) %>%
#   filter(onset < "2020-03-19") %>%
#   pull(case_id)
#
# bind_rows(trees) %>%
#   filter(to %in% ids_to_check) %>%
#   mutate(SI = to_date - from_date) %>%
#   ggplot(aes(x = SI)) +
#   geom_histogram(binwidth = 1)


# With posterior support --------------------------------------------------

# How many times is each case_id identified as an introduction?
bind_rows(trees) %>%
  filter(is.na(from)) %>%
  group_by(to) %>%
  summarise(n = n(),
            support = n / length(trees),
            .groups = "drop") %>%
  filter(support > 0.1) %>%
  arrange(to)
