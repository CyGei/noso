source(here::here("scripts", "helpers.R"))
load_libraries()
pacman::p_load(o2groups, dplyr, here)


moVSdelta <- function(o2groups.params) {
  sim <- o2groups::simulate_groups(
    duration = o2groups.params$duration,
    group_n = o2groups.params$group_n,
    size = o2groups.params$size,
    name = o2groups.params$name,
    gamma = o2groups.params$gamma,
    intro_n = o2groups.params$intro_n,
    r0 = o2groups.params$r0,
    generation_time = o2groups.params$generation_time,
    incubation_period = o2groups.params$incubation_period
  )
  true_delta <- linktree::gamma2delta(o2groups.params$gamma)
  f <- setNames(o2groups.params$size / sum(o2groups.params$size),
                o2groups.params$name)
  # reformat data to allow cut_tree_by_date()
  sim <- sim %>%
    mutate(from_date = sim$date_onset[match(source, sim$id)]) %>%
    rename(to_date = date_onset)
  peak_dates <-
    get_peak(date = sim$to_date,
             group = sim$group)
  cases <- sim %>%
    group_by(group) %>%
    summarise(cases = n(), .groups = "drop")

  info <- left_join(peak_dates, cases, by = "group") %>%
    mutate(true_delta = true_delta) %>%
    rename(level = group)

  ttab <-
    linktree:::ttable(
      from = sim$source_group,
      to = sim$group,
      levels = o2groups.params$name
    )

  delta <- lapply(o2groups.params$name, function(level) {
    cut_tree <- cut_tree_by_date(tree = sim,
                                 cutoff_date = info$observed_peak[info$level == level])
    tryCatch({
      delta <- get_delta(
        from = cut_tree$source_group,
        to = cut_tree$group,
        f = f,
        alpha = 0.05
      ) %>%
        filter(group == level) %>%
        select(-group)
      return(delta)
    }, error = function(e) {
      return(data.frame(
        est = NA_real_,
        lwr = NA_real_,
        upr = NA_real_
      ))
    })
  }) %>% bind_rows(.id = "level") %>%
    mutate(metric = "delta",
           level = factor(level, labels = o2groups.params$name))

  mo <- draw_mo(
    from = sim$source_group,
    to = sim$group,
    level = o2groups.params$name,
    args = list(
      diag = TRUE,
      from_id = sim$source,
      to_id = sim$id
    )
  ) %>%
    as_tibble(rownames = "level") %>%
    rename(est = value) %>%
    mutate(
      est = gamma2delta(est),
      lwr = NA_real_,
      upr = NA_real_,
      metric = "mo"
    )

  out <- bind_rows(delta, mo) %>%
    left_join(info, by = "level")

  return(list(out = out,
              ttab = ttab))

}

set.seed(123)
o2groups.params <- list(
  duration = 365,
  group_n = 2,
  size = c(100, 300),
  name = c("HCW", "patient"),
  intro_n = c(1, 1),
  r0 = c(2, 2),
  generation_time = c(0, 0.1, 0.2, 0.4, 0.2, 0.1, 0),
  incubation_period = sample(1:14, 1000, replace = TRUE),
  gamma = c(14, 0.6)
)
set.seed(123)
results <- moVSdelta(o2groups.params)
out <- results$out
ttab <- results$ttab

out %>%
  ggplot(aes(
    x = level,
    y = est,
    ymin = lwr,
    ymax = upr,
    col = metric
  )) +
  geom_hline(aes(yintercept = 0),
             linetype = "dotted") +
  geom_point(
    aes(y = true_delta,
        shape = "truth"),
    size = 5,
    col = "#5c5e60"
  ) +
  geom_pointrange(size = 0.5, position = position_dodge(width = 0.25)) +
  scale_color_manual("metric",
                     values =
                       c("delta" = "red", "mo" = "blue")) +
  scale_shape_manual("", values = c("truth" = 15)) +
  theme_bw() +
  labs(x = "",
       y = "Estimate") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(-1, 1))



ttab
ttab %>%
  as_tibble() %>%
  uncount(n) %>%
  select(from, to) %>%
  sjPlot::sjtab(fun = "xtab",
                var.labels = c("infector", "infectee"),
                show.row.prc=T,
                show.col.prc=T,
                show.summary=T,
                show.exp=T,
                show.legend=T)





# Simulation --------------------------------------------------------------
pacman::p_load(furrr)
plan(multisession, workers = availableCores() - 1)
#generate list of randomly generated parameters
set.seed(123)
o2groups.params_list <- lapply(1:100, function(i) {
    duration = 365
    group_n = round(runif(1, 2, 10))
    size = round(runif(group_n, 100, 1000))
    name = letters[1:group_n]
    intro_n = rep(1, group_n)
    r0 = rnorm(group_n, 2.3, 0.5)
    generation_time = c(0, 0.1, 0.2, 0.4, 0.2, 0.1, 0)
    incubation_period = sample(1:14, 1000, replace = TRUE)
    gamma = delta2gamma(runif(group_n, -1, 1))

    o2groups.params <- list(
        duration = duration,
        group_n = group_n,
        size = size,
        name = name,
        intro_n = intro_n,
        r0 = r0,
        generation_time = generation_time,
        incubation_period = incubation_period,
        gamma = gamma
    )
})

#run simulation
set.seed(123)
results_list <- future_map(o2groups.params_list, moVSdelta,
                           .options = furrr_options(seed = TRUE))
out_list <- future_map(results_list, function(x) x$out,
                       .options = furrr_options(seed = TRUE))
ttab_list <- future_map(results_list, function(x) x$tta,
                        .options = furrr_options(seed = TRUE))

params_df <- lapply(o2groups.params_list, function(x) {
    data.frame(
        group_n = x$group_n,
        size = x$size,
        f = x$size / sum(x$size),
        r0 = x$r0,
        gamma = x$gamma,
        level = x$name
    )
}) %>%
    bind_rows(.id = "sim")

out_list %>%
  bind_rows(.id = "sim") %>%
  left_join(params_df, by = c("sim", "level")) %>%
  mutate(error = est - true_delta) %>%
  ggplot(aes(x = true_delta, y = est, col = f)) +
  geom_point()+
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  facet_wrap(~metric)+
  scale_y_continuous(limits = c(-1, 1))+
  scale_colour_gradientn(colours = c("red", "black", "red"),
                         values = c(0, 0.5, 1))+
  coord_fixed()
