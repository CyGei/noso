source(here::here("scripts", "helpers.R"))
load_libraries()
paper <- c("JHI2021", "eLife2022")[2]

input_path <- here::here("data", paper, "input/")
output_path <- here::here("data", paper, "output/")
load_data(input_path)

out <- out %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id)
out <- filter_alpha_by_kappa(out, 1L)

trees <- get_trees(
  out = out,
  ids = linelist$case_id,
  group = linelist$group,
  date = linelist$onset
)
epicurve()
metrics <- list("delta" = draw_delta,
                "rho" = draw_rho)

pacman::p_load(furrr)
future::plan(list(
  future::tweak(future::multisession, workers = length(metrics)),
  future::tweak(
    future::multisession,
    workers = future::availableCores() - (length(metrics) + 2)
  )
))
fHCW <- 0.33
set.seed(123)
system.time({
  array <- furrr::future_map(names(metrics), function(metric_name) {
    draw_array(
      from_col = "from_group",
      to_col = "to_group",
      levels = c("hcw", "patient"),
      trees = trees,
      draw_function = metrics[[metric_name]],
      args = list(
        f = c("hcw" = fHCW, "patient" = 1 - fHCW),
        from_id = "from",
        to_id = "to",
        diag = TRUE
      ),
      cutoff_dates = cutoff_dates,
      n_samples = ifelse(metric_name == "rho", 1, 1000)
    )
  }, .options = furrr_options(
    seed = TRUE,
    globals = c(
      "fHCW",
      "metrics",
      "trees",
      "cutoff_dates",
      "cut_tree_by_date",
      grep("draw_*", names(.GlobalEnv), value = TRUE),
      grep("*_formula", names(.GlobalEnv), value = TRUE)
    )
  ))
})
saveRDS(array, here::here(output_path, "array.rds"))
array <- readRDS(here::here(output_path, "array.rds"))

#rescale rho values with gamma2delta()
rho_idx <- which(names(metrics) == "rho")
array[[rho_idx]] <- apply(array[[rho_idx]], 1:4, gamma2delta)

CrI <- lapply(array, function(x)
  draw_CrI(x, c(2, 3)))

df <- lapply(CrI, function(x) {
  x %>%
    reshape2::melt() %>%
    pivot_wider(names_from = Var1, values_from = value)
}) %>%
  bind_rows(.id = "metric") %>%
  mutate(
    metric = factor(
      metric,
      levels = 1:length(metrics),
      labels = names(metrics)
    ),
    target = case_when(
      metric == "gamma" ~ 1,
      metric == "delta" ~ 0,
      metric == "rho" ~ 0,
      metric == "R0" ~ 1
    )
  )
df %>%
  mutate(cutoff = as.Date(cutoff)) %>%
  arrange(cutoff) %>%
  mutate(across(where(is.numeric), \(x) round(x, 2))) %>%
  print(n = 20)

p_metrics <- ggplot(df, aes(
  x = as.Date(cutoff),
  y = mean,
  ymin = lwr,
  ymax = upr,
  col = level
)) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1,
             labeller = label_parsed) +
  geom_pointrange(position =
                    position_dodge(width = 0.65)) +
  geom_hline(aes(yintercept = target), linetype = "dashed") +
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "delta" ~ scale_y_continuous(limit = c(-1, 1), breaks = seq(-1, 1, 0.5)),
      metric == "gamma" ~ scale_y_continuous(limit = c(0, 5), breaks = seq(0, 20, 1)),
      metric == "rho" ~ scale_y_continuous(limit = c(-1, 1), breaks = seq(-1, 1, 0.5)),
      metric == "R0" ~ scale_y_continuous(limit = c(0, 5), breaks = seq(0, 5, 1))
    )
  ) +
  labs(x = "", y = "metric value") +
  theme_noso() +
  theme(legend.position = "none")


# ncases <- table(linelist$group)
# p_sus <- linelist %>%
#   arrange(onset) %>%
#   group_by(group, onset) %>%
#   summarise(cases = n(), .groups = "drop") %>%
#   group_by(group) %>%
#   mutate(
#     cum_cases = cumsum(cases),
#     prop_susceptible = 1 - cum_cases / as.integer(ncases[group]),
#     cum_cases_scaled = cum_cases / max(ncases)
#   ) %>%
#   ggplot(aes(x = onset, y = prop_susceptible, col = group)) +
#   facet_wrap(~ group, ncol = 1)+
#   geom_col(aes(y = cum_cases_scaled), col = NA, alpha = 0.5) +
#   geom_line() +
#   geom_point() +
#   geom_hline(yintercept = 0.5, linetype = "dashed") +
#   scale_y_continuous(sec.axis = sec_axis(~ . * max(ncases),
#                                          name = "Cumulative cases")) +
#   labs(x = "", y = "Proportion susceptible") +
#   theme_noso() +
#   theme(legend.position = "none")


cowplot::plot_grid(
  p_metrics,
  epicurve(),
  align = "v",
  ncol = 1,
  rel_heights = c(1,0.7),
  labels = "AUTO"
)

# TTABLE ------------------------------------------------------------------
peaks <- get_peak(
  linelist$onset,
  linelist$group
) %>%
  arrange(observed_peak)

#Let's get the average transmission table
info <- tibble(
  group = c(peaks$group,"outbreak"),
  peak = c(
    peaks$observed,
    as.Date("2100-01-01") # arbitrary, post-outbreak
  )
)
ttabs <- future_map(trees, function(tree) {
  cut_trees <- lapply(info$peak, function(peak) {
    cut_tree_by_date(tree, peak) %>%
      mutate(cut_at = peak, peak_group = info$group[info$peak == peak])
  })

  ttabs <- lapply(cut_trees, function(cut_tree) {
    cut_at <- unique(cut_tree$cut_at)
    peak_group <- unique(cut_tree$peak_group)

    ttab <- linktree:::ttable(
      from = cut_tree$from_group ,
      to = cut_tree$to_group,
      levels = c("hcw", "patient")
    )
    ttab <- as_tibble(ttab) %>%
      mutate(cut_at = cut_at, peak_group = peak_group)
  }) %>% bind_rows()
}) %>% bind_rows(.id = "tree")

peak_ttabs <- split(ttabs, ttabs$peak_group) %>%
  map(~ {
    .x %>% select(from, to, n) %>%
      group_by(from, to) %>%
      summarise(n = round(mean(n)), .groups = "drop")
  })
# Average Xtabs
sjtabs <- map(peak_ttabs, ~ {
  .x %>%
    uncount(n) %>%
    select(from, to) %>%
    sjPlot::sjtab(
      fun = "xtab",
      var.labels = c("infector", "infectee"),
      show.row.prc = T,
      show.col.prc = T,
      show.summary = T,
      show.exp = T,
      show.legend = T
    )
})
names(peak_ttabs)
sjtabs[[1]]
sjtabs[[2]]
sjtabs[[3]]

linktree:::gamma_formula(0.57, 0.8)


# trees %>%
#   future_map( ~ {
#     linktree:::ttable(
#       from = .x$from_group,
#       to = .x$to_group,
#       levels = c("hcw", "patient")
#     ) %>%
#       as_tibble()
#   }) %>%
#   bind_rows() %>%
#   group_by(from, to) %>%
#   summarise(mean = round(mean(n)), .groups = "drop")




# R0 ----------------------------------------------------------------------
set.seed(123)
array <- draw_array(
  from_col = "from_group",
  to_col = "to_group",
  levels = c("hcw", "patient"),
  trees = trees,
  draw_function = draw_R0,
  args = list(
    f = c("hcw" = fHCW, "patient" = 1-fHCW),
    from_id = "from",
    to_id = "to"
  ),
  cutoff_dates = cutoff_dates,
  n_samples = 1000
)

CrI <- draw_CrI(array, c(2, 3)) %>%
  reshape2::melt() %>%
  pivot_wider(names_from = Var1, values_from = value)

ggplot(CrI, aes(
  x = as.Date(cutoff),
  y = mean,
  ymin = lwr,
  ymax = upr,
  col = level
)) +
  geom_pointrange(position =
                    position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  theme_noso() +
  coord_cartesian(ylim = c(0, 4)) +
  labs(x = "cutoff date", y = expression(R[0]))
