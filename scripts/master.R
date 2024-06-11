source(here::here("scripts", "helpers.R"))
load_libraries()
paper <- c("JHI2021", "eLife2022")[2] #chose which paper to run
load_data(paper)
output_path <- here::here("data", paper, "output")

epicurve()
plan(multisession, workers = availableCores())
trees <- furrr::future_map(seq_along(cutoff_dates), function(x) {
  linelist_cut <- linelist %>%
    filter(onset <= cutoff_dates[x])

  out[[x]] %>%
    filter(step > 500) %>%
    identify(ids = linelist_cut$case_id) %>%
    filter_alpha_by_kappa(1L) %>%
    get_trees(
      ids = linelist_cut$case_id,
      group = linelist_cut$group,
      date = linelist_cut$onset
    )
})

linelist %>%
  group_by(onset, group) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(group) %>%
  mutate(cum_n = cumsum(n))

# Compute array -----------------------------------------------------------
metrics <- list("delta" = draw_delta, "rho" = draw_rho)

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
      n_samples = ifelse(metric_name == "rho", 1, 1000)
    )
  }, .options = furrr_options(
    seed = TRUE,
    globals = c(
      "fHCW",
      "metrics",
      "trees",
      grep("draw_*", names(.GlobalEnv), value = TRUE),
      grep("*_formula", names(.GlobalEnv), value = TRUE)
    )
  ))
})
saveRDS(array, here::here(output_path, "array.rds"))

# Compute credible intervals ----------------------------------------------
array <- readRDS(here::here(output_path, "array.rds"))

#rescale rho values with gamma2delta()
rho_idx <- which(names(metrics) == "rho")
array[[rho_idx]] <- apply(array[[rho_idx]], 1:4, gamma2delta)

df <-
  reshape2::melt(array) %>%
  group_by(level, cutoff, L1) %>%
  drop_na() %>%
  summarise(
    mean = mean(value),
    lwr = quantile(value, 0.025),
    upr = quantile(value, 0.975),
    .groups = "drop"
  ) %>%
  mutate(
    metric = factor(
      L1,
      levels = 1:length(metrics),
      labels = names(metrics)
    ),
    target = case_when(
      metric == "gamma" ~ 1,
      metric == "delta" ~ 0,
      metric == "rho" ~ 0,
      metric == "R0" ~ 1
    ),
    cutoff_date = cutoff_dates[cutoff],
    #if target is contained in the CrI then it is significant
    significant = ifelse(target >= lwr &
                           target <= upr, "not significant", "significant")
  ) %>%
  arrange(cutoff_date)

df %>%
  dplyr::mutate(across(where(is.numeric), \(x) round(x, 2))) %>%
  print(n = 20)

# Plot --------------------------------------------------------------------
p_metrics <- ggplot(df,
                    aes(
                      x = cutoff_date,
                      y = mean,
                      ymin = lwr,
                      ymax = upr,
                      col = level,
                      group = interaction(metric, level, cutoff_date)
                    )) +
  facet_wrap(
    ~ metric,
    ncol = 1,
    labeller = label_parsed,
    strip.position = "left"
  ) +
  geom_hline(aes(yintercept = target), linetype = "dashed") +
  geom_errorbar(position = position_dodge(width = 0.5),
                width = 0.5,
                linewidth = 0.5) +
  geom_point(aes(shape = significant),
             position = position_dodge(width = 0.5),
             size = 2.5,
             fill = "white",
             alpha = 0.9) +
  scale_shape_manual(values =
                       c("significant" = 16, "not significant" = 21)) +
  theme_noso() +
  labs(x = "", y = "") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 16,
      margin = margin(r = -0.5, l = -1, 0, 0, "pt")
    ),
    panel.spacing.y = unit(0.5, "lines")
  )

cowplot::plot_grid(
  p_metrics + theme(plot.margin = unit(c(0.1, 0.2, -1, 0), "lines")),
  epicurve(facet = TRUE),
  align = "v",
  ncol = 1,
  rel_heights = c(1, 0.7),
  labels = "AUTO"
)
