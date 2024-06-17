source(here::here("scripts", "helpers.R"))
load_libraries()
paper <- c("JHI2021", "eLife2022")[2] #chose which paper to run
load_data(paper)
output_path <- here::here("data", paper, "output")
epicurve()
out <- out[[length(out)]]
# Total Introductions
alpha_cols <- grep("alpha_", colnames(out))
sum(is.na(out[sample(2:999, 1), alpha_cols]))
trees <- out %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  get_trees(ids = linelist$case_id,
            group = linelist$group,
            date = linelist$onset,
            kappa = TRUE,
            t_inf = TRUE)

window <- ifelse(paper == "JHI2021", 5, 7)
start <- min(linelist$onset)
end <- max(linelist$onset)
cutoff_breaks <- c(seq.Date(start, end, by = window), end) %>% unique()


# Estimate the average number of introductions per time window ------------
intro_df <- furrr::future_map(seq_along(cutoff_breaks), function(i) {
  window_start <- cutoff_breaks[i]
  window_end <- cutoff_breaks[i + 1] - 1

  intros <- lapply(trees, function(tree) {
    intro_window <- tree %>%
      filter(is.na(from)) %>%
      arrange(to_date) %>%
      filter(to_date >= window_start & to_date <= window_end) %>%
      group_by(to_group) %>%
      summarise(n_intros = sum(is.na(from)), .groups = "drop") %>%
      mutate(window_start = window_start, window_end = window_end)

    # Create a data frame with all possible groups
    all_groups <- data.frame(
      to_group = unique(linelist$group),
      window_start = window_start,
      window_end = window_end
    )

    # Merge the all_groups data frame with the intro_window data frame
    intro_window <- merge(all_groups, intro_window, all.x = TRUE)

    # Replace NA values with 0
    intro_window$n_intros[is.na(intro_window$n_intros)] <- 0

    return(intro_window)
  }) %>%
    bind_rows()
  return(intros)
}) %>%
  bind_rows(.id = "window_id") %>%
  mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)

intro_summary <- intro_df %>%
  group_by(to_group, window_id) %>%
  summarise(
    mean = mean(n_intros),
    lwr = quantile(n_intros, 0.025),
    upr = quantile(n_intros, 0.975),
    window_start = first(window_start),
    window_end = first(window_end),
    window_median = first(window_median),
    .groups = "drop"
  ) %>%
  group_by(window_id) %>%
  mutate(global_intro = sum(mean)) %>%
  ungroup()

p_intros <- intro_df %>%
  ggplot(aes(
    x = window_median,
    y = n_intros,
    fill = to_group,
    group = interaction(to_group, window_id)
  )) +
  geom_errorbarh(
    data = intro_summary,
    aes(xmin = window_start, xmax = window_end, y = global_intro),
    height = 0.05,
    linetype = "solid",
    col = "#636363"
  ) +
  geom_violin(
    adjust = 1.8,
    col = NA,
    position = position_dodge(width = 4),
    alpha = 0.7
  ) +
  geom_pointrange(
    data = intro_summary,
    aes(
      x = window_median,
      y = mean,
      ymin = lwr,
      ymax = upr,
    ),
    position = position_dodge(width = 4),
    pch = 21,
    stroke = 0.5,
    size = 0.75
  ) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  theme_noso(day_break = 7, date = TRUE) +
  coord_cartesian(ylim = c(0, 14)) +
  labs(x = "", y = "Importations")


cowplot::plot_grid(
  p_intros + theme(legend.position = "none"),
  epicurve(day_break = 7),
  ncol = 1,
  rel_heights = c(2, 1.5),
  align = "v",
  labels = "AUTO"
)
intro_summary %>% pull(mean) %>% sum()


#eLife2022: no introductions in the 1st 2 weeks!
ids_to_check <- linelist %>%
  arrange(onset) %>%
  filter(onset < "2020-03-19") %>%
  pull(case_id)

bind_rows(trees) %>%
  filter(to %in% ids_to_check) %>%
  mutate(SI = to_date - from_date) %>%
  ggplot(aes(x = SI)) +
  geom_histogram(binwidth = 1)


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
