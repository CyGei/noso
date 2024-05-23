# The purpose of this script is to estimate the reproduction number
# using transmission chain data.
# The effective reproductive number describes how many new infections an individual
# infected  causes on average in a population which is subject to a certain degree of immunity and intervention measures.

source(here::here("scripts", "helpers.R"))
load_libraries()
load_data()
outbreaker_chains <- outbreaker_chains %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id)
outbreaker_chains <- filter_alpha_by_kappa(outbreaker_chains, 1)

trees <- get_trees(
  out = outbreaker_chains,
  ids = linelist$case_id,
  group = linelist$group,
  date = linelist$onset_inferred
)

cutoff_dates <- unique(c(seq(min(
  as.Date(linelist$onset_inferred)
), max(
  as.Date(linelist$onset_inferred)
), by = 7), max(linelist$onset_inferred)))





# get_R functions -------------------------------------------------------------
# get_Re <- function(linelist, tree, by_group = FALSE) {
#   tree <- tree[!is.na(tree$from),]
#   if (by_group) {
#     cases <- table(linelist$group)
#     transmissions <- tapply(tree$from, tree$from_group, length)
#
#     #FIXED: transmissions should have as many elements as cases
#     transmissions <- transmissions[names(cases)]
#     transmissions[is.na(transmissions)] <- 0
#     names(transmissions) <- names(cases)
#
#     R <- data.frame(
#       from_group = names(transmissions),
#       transmissions = as.vector(transmissions),
#       cases = as.vector(cases),
#       R = as.vector(transmissions / cases)
#     )
#   } else {
#     transmissions <- nrow(tree)
#     cases <- nrow(linelist)
#     R <- data.frame(R = transmissions / cases)
#   }
#
#   return(R)
# }


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



# Re: estimating R using 7 days time windows ------------------------------
pacman::p_load(furrr)
plan(multisession, workers = length(cutoff_dates))

R_df <- furrr::future_map(cutoff_dates[-1], function(cutoff_date) {
  window <- 7
  cutoff_date <- as.Date(cutoff_date)
  R <- lapply(trees, function(tree) {
    cut_tree <-
      tree %>%
      drop_na(from) %>%
      mutate(from_date = as.Date(from_date)) %>%
      filter(from_date >= (cutoff_date - window) &
               from_date <= cutoff_date)

    cut_linelist <- linelist %>%
      filter(onset_inferred >= (cutoff_date - window) &
               onset_inferred <= cutoff_date)

    get_Re(cut_linelist, cut_tree, by_group = TRUE)
  }) %>%
    bind_rows() %>%
    dplyr::mutate(window_end = cutoff_date, window_start = cutoff_date - window)
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
  theme_noso() +
  labs(x = "", y = expression(R["e"]))

cowplot::plot_grid(
  p_Re + theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, 4)),
  epicurve(),
  ncol = 1,
  rel_heights = c(2, 1.5),
  align = "v",
  labels = "AUTO"
)



# Estimating Re from start of outbreak up to a cutoff ------------------------------

R_df <- furrr::future_map(cutoff_dates[-1], function(cutoff_date) {
  cutoff_date <- as.Date(cutoff_date)
  R <- lapply(trees, function(tree) {
    cut_tree <-
      tree %>%
      drop_na(from) %>%
      mutate(from_date = as.Date(from_date)) %>%
      filter(from_date <= cutoff_date)

    cut_linelist <- linelist %>%
      filter(onset_inferred <= cutoff_date)

    get_Re(cut_linelist, cut_tree, by_group = TRUE)
  }) %>%
    bind_rows() %>%
    dplyr::mutate(window_end = cutoff_date,
                  window_start = cutoff_date - as.numeric(cutoff_date -
                                                            min(as.Date(
                                                              linelist$onset_inferred
                                                            ))))
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
    x = window_end,
    y = Re,
    fill = from_group,
    group = interaction(from_group, window_id)
  )) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbarh(
    data = Rsummary,
    aes(
      xmin = window_start,
      xmax = window_end,
      y = global_Re,
      group = window_id
    ),
    height = 0.1,
    linewidth = 0.2,
    linetype = "solid",
    col = "#636363"
  ) +
  geom_violin(
    adjust = 1.8,
    col = NA,
    position = position_dodge(width = 1),
    alpha = 0.7
  ) +
  geom_pointrange(
    data = Rsummary,
    aes(
      x = window_end,
      y = mean,
      ymin = lwr,
      ymax = upr,
    ),
    position = position_dodge(width = 1),
    pch = 21,
    stroke = 0.5,
    size = 0.75
  ) +
  theme_noso() +
  labs(x = "", y = expression(R["e"]))

cowplot::plot_grid(
  p_Re + theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, 4)),
  epicurve(),
  ncol = 1,
  rel_heights = c(2, 1.5),
  align = "v",
  labels = "AUTO"
)

#R0
Rsummary %>%
  mutate(across(where(is.numeric), \(x) round(x, 2))) %>%
  filter(window_id == 2)

#remaining susceptible at cutoff
linelist %>%
  filter(onset_inferred <= cutoff_dates[3]) %>%
  {100 - (nrow(.) / 148 * 100)}

