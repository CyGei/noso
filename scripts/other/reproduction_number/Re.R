# The purpose of this script is to estimate the reproduction number
# using transmission chain data.
# The effective reproductive number describes how many new infections an individual
# infected  causes on average in a population which is subject to a certain degree of immunity and intervention measures.
# This is not applicable in real-time.

source(here::here("scripts", "helpers.R"))
load_libraries()
paper <- c("JHI2021", "eLife2022")[1] #chose which paper to run
load_data(paper)
output_path <- here::here("data", paper, "output")
epicurve()

#complete tree
trees <- out[[length(out)]] %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  get_trees(ids = linelist$case_id,
            group = linelist$group,
            date = linelist$onset,
            kappa = TRUE)
#C218   C203  H2028
get_Re <- function(linelist, tree, by_group = FALSE) {
  tree <- tree[!is.na(tree$from), ]

  individualR <- count(tree, from)

  df <- merge(
    x = linelist[, c("case_id", "group")],
    y = individualR,
    by.x = "case_id",
    by.y = "from",
    all.x = TRUE
  )
  df$n <- ifelse(is.na(df$n), 0, df$n)

  if (by_group) {
    Re <- df %>%
      rename(from_group = group) %>%
      group_by(from_group) %>%
      summarise(Re = mean(n),
                cases = n(),
                .groups = "drop")
  } else {
    Re <- data.frame(Re = mean(df$n), cases = nrow(df))
  }

  return(Re)
}


# Re: estimating R using 5/7 days time windows ------------------------------
window <- ifelse(paper == "JHI2021", 5, 7)
start <- min(linelist$onset) + window
end <- max(linelist$onset)
cutoff_breaks <- c(seq.Date(start, end, by = window), end) %>% unique()

pacman::p_load(furrr)
plan(multisession, workers = length(cutoff_breaks))
R_df <- furrr::future_map(cutoff_breaks, function(cutoff_break) {
  window_start <- cutoff_break - window
  window_end <- cutoff_break

  R <- lapply(trees, function(tree) {
    cut_tree <-
      tree %>%
      drop_na(from) %>%
      filter(from_date >= window_start &
               from_date <= window_end)

    cut_linelist <- linelist %>%
      filter(onset >= window_start &
               onset <= window_end)

    get_Re(cut_linelist, cut_tree, by_group = TRUE)
  }) %>%
    bind_rows() %>%
    dplyr::mutate(window_end = window_end, window_start = window_start)
}) %>%
  bind_rows(.id = "window_id") %>%
  mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)

Rsummary <- R_df %>%
  group_by(from_group, window_id) %>%
  summarise(
    mean = mean(Re),
    lwr = quantile(Re, 0.025),
    upr = quantile(Re, 0.975),
    window_start = first(window_start),
    window_end = first(window_end),
    window_median = first(window_median),
    cases = first(cases),
    .groups = "drop"
  ) %>%
  group_by(window_id) %>%
  mutate(global_Re = mean(rep(mean, cases))) %>%
  ungroup()

p_Re <- R_df %>%
  ggplot(aes(
    x = window_median,
    y = Re,
    fill = from_group,
    group = interaction(from_group, window_id)
  )) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbarh(
    data = Rsummary,
    aes(xmin = window_start, xmax = window_end, y = global_Re),
    height = 0.1,
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
    data = Rsummary,
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
  theme_noso(day_break = 7, date = TRUE) +
  labs(x = "", y = expression(R["e"]))

cowplot::plot_grid(
  p_Re + theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, 4)),
  epicurve(day_break = 7),
  ncol = 1,
  rel_heights = c(2, 1.5),
  align = "v",
  labels = "AUTO"
)

