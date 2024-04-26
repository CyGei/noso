# This script aims at understanding the relationship between gamma and mo.
# We build a funciton to simulate an outbreak with group specific assortativity coefficients (gamma)
# and plot the estimated gamma and mo values.

source(here::here("scripts", "helpers.R"))
load_libraries()
pacman::p_load(o2groups, furrr, dplyr, here)


moVSdelta <- function(o2groups.params) {
  # simulate outbreak
  simulations <- o2groups::simulate_groups_furrr(
    core = 10,
    sim_n = 100,
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
  #plan(sequential)

  true_delta <- linktree::gamma2delta(o2groups.params$gamma)
  f <- setNames(o2groups.params$size / sum(o2groups.params$size),
                o2groups.params$name)
  # reformat data to allow cut_tree_by_date()
  simulations <-
    lapply(simulations, function(x) {
      x %>%
        mutate(from_date = x$date_onset[match(source, x$id)]) %>%
        rename(to_date = date_onset)
    })
  sim_df <-  bind_rows(simulations, .id = "sim_n")

  # estimate the average peak date (the date with the highest number of cases) across simulations
  peak_dates <-
    sim_df %>%
    group_by(sim_n, to_date, group) %>%
    summarise(count = n(), .groups = "drop") %>%
    ungroup() %>%
    group_by(group, to_date) %>%
    summarise(average_count = mean(count)) %>%
    slice(which.max(average_count)) %>%
    select(group, peak_date = to_date)

  # average number of cases in each group
  cases <- sim_df %>%
    group_by(sim_n, group) %>%
    count() %>%
    ungroup() %>%
    group_by(group) %>%
    summarise(cases = mean(n))

  #averge ttable
  ttab <- lapply(simulations, function(x) {
    linktree:::ttable(
      from = x$source_group,
      to = x$group,
      levels = o2groups.params$name
    ) %>%
      as.data.frame()
  }) %>%
    bind_rows(.id = "sim_n") %>%
    group_by(from, to) %>%
    summarise(Freq = mean(Freq), .groups = "drop")


  info <- left_join(peak_dates, cases, by = "group") %>%
    ungroup() %>%
    mutate(true_delta = true_delta) %>%
    rename(level = group)

  # Estimate gamma and Mo
  metrics <- list("delta" = draw_delta,
                  "Mo" = draw_mo)

  future::plan(list(
    future::tweak(future::multisession,
                  workers = 2),
    future::tweak(future::multisession,
                  workers = availableCores() - 4)
  ))


  array <- future_map(names(metrics),
                      function(metric_name) {
                        draw_array(
                          from_col = "source_group",
                          to_col = "group",
                          levels = o2groups.params$name,
                          trees = simulations,
                          draw_function = metrics[[metric_name]],
                          args = list(
                            f = f,
                            from_id = "source",
                            to_id = "id",
                            diag = TRUE
                          ),
                          cutoff_dates = round(peak_dates$peak_date),
                          n_samples = ifelse(metric_name == "Mo", 1, 1000)
                        )
                      },
                      .options = furrr_options(
                        seed = TRUE,
                        globals = list(
                          simulations = simulations,
                          metrics = metrics,
                          f = f,
                          peak_dates = peak_dates,
                          o2groups.params = o2groups.params,
                          draw_CrI = draw_CrI,
                          draw_array = draw_array,
                          draw_delta = draw_delta,
                          draw_mo = draw_mo,
                          draw_gamma = draw_gamma,
                          draw_pi = draw_pi,
                          mo_formula = mo_formula,
                          cut_tree_by_date = cut_tree_by_date
                        )
                      ))
  #scale Mo
  array[[2]] <-
    apply(array[[2]], MARGIN = 1:length(dim(array[[2]])), FUN = gamma2delta)

  CrI <- lapply(array, \(x) {
    draw_CrI(x, c(2, 3)) %>%
      reshape2::melt() %>%
      pivot_wider(names_from = Var1, values_from = value)
  }) %>%
    bind_rows(.id = "metric") %>%
    mutate(metric = factor(metric, labels = c("delta", "Mo"))) %>%
    left_join(info, by = "level")

  return(list(
    CrI = CrI,
    ttab = ttab
  ))
}


set.seed(123)
o2groups.params <- list(
  duration = 365,
  group_n = 2,
  size = c(300, 300),
  name = c("HCW", "patient"),
  intro_n = c(1, 1),
  r0 = c(2, 2),
  generation_time = c(0, 0.1, 0.2, 0.4, 0.2, 0.1, 0),
  incubation_period = sample(1:14, 1000, replace = TRUE),
  gamma = c(1, 1)
)

x = moVSdelta(o2groups.params)$CrI %>%
  filter(cutoff == round(peak_date))

x %>%
  ggplot(aes(
    x = level,
    y = mean,
    ymin = lwr,
    ymax = upr,
    col = metric
  )) +
  geom_hline(aes(yintercept = 0),
             linetype = "dotted") +
  geom_point(aes(y = true_delta),
             shape = 18,
             size = 5,
             col = "black") +
  geom_pointrange(size = 0.5, position = position_dodge(width = 0.1)) +
  scale_color_manual("Metric",
                     values =
                       c("delta" = "red", "Mo" = "blue")) +
  theme_bw() +
  labs(x = "",
       y = "Estimate") +
  theme(legend.position = "bottom")+
  scale_y_continuous(limits = c(-1, 1))

x




#use moVSdelta for all combinations of gamma
gamma_grid <- expand.grid(HCW = c(1/4, 1, 4),
                          patient = c(1/4, 1, 4)) %>%
  #remove duplicated pairs regrdaless of order
  filter(HCW >= patient)

o2groups.params <- list(
  duration = 365,
  group_n = 2,
  size = c(300, 300),
  name = c("HCW", "patient"),
  intro_n = c(1, 1),
  r0 = c(2, 2),
  generation_time = c(0, 0.1, 0.2, 0.4, 0.2, 0.1, 0),
  incubation_period = sample(1:14, 1000, replace = TRUE)
)

moVSdelta_grid <- lapply(1:nrow(gamma_grid), function(i) {
  o2groups.params$gamma <- as.numeric(gamma_grid[i,])
  moVSdelta(o2groups.params)
})

moVSdelta_grid_CrI <- lapply(moVSdelta_grid, function(x) {
  x$CrI %>%
    rowwise() %>%
    mutate(across(c(mean, lwr, upr), ~ unlist(.x)[1])) %>%
    ungroup() %>%
    filter(cutoff == round(peak_date))
})

bind_rows(moVSdelta_grid_CrI, .id = "scenario") %>%
  mutate(scenario = factor(scenario,
                           labels = lapply(1:nrow(gamma_grid), function(i) {
                             paste0("\u03b3: ",
                                    "HCW=", round(gamma_grid[i, 1], 2),
                                    ", patient=", round(gamma_grid[i, 2], 2))
                           }) %>% unlist())) %>%
  ggplot(aes(
    x = level,
    y = mean,
    ymin = lwr,
    ymax = upr,
    col = metric
  )) +
  facet_wrap(~ scenario) +
  geom_text(aes(label = paste0("N=", round(cases)),
                group = level),
            col = "black",
            y = 0.9) +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  geom_point(aes(y = true_delta),
             shape = 18,
             size = 3,
             col = "black") +
  geom_pointrange(size = 0.5, position = position_dodge(width = 0.1)) +
  scale_color_manual("Metric", values = c("delta" = "red", "Mo" = "blue")) +
  theme_bw() +
  labs(x = "",
       y = "Delta est") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(-1, 1))


lapply(moVSdelta_grid, \(x) x$ttab)[[2]] %>%
  mutate(Freq = round(Freq)) %>%
  xtabs(Freq ~ from + to, data = .)

lapply(moVSdelta_grid, \(x) x$ttab)[[2]] %>%
  mutate(Freq = round(Freq)) %>%
  uncount(Freq) %>%
  select(from, to) %>%
  sjPlot::sjtab(fun = "xtab",
                var.labels = c("from", "to"),
                show.row.prc=T,
                show.col.prc=T,
                show.summary=T,
                show.exp=T,
                show.legend=T)
