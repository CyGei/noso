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
  #warning if sim has less than 50 rows
  if (nrow(sim) < 50) {
    warning("simulation has less than 50 cases, re-run with another seed")
  }

  true_delta <- linktree::gamma2delta(o2groups.params$gamma)
  f <- setNames(o2groups.params$size / sum(o2groups.params$size),
                o2groups.params$name)
  # reformat data to allow cut_tree_by_date()
  sim <- sim %>%
    mutate(from_date = sim$date_onset[match(source, sim$id)]) %>%
    rename(to_date = date_onset)

  #randomly remove half of the pairs where the infectee is a hcw
  sim <- sim %>%
    group_by(group) %>%
    mutate(keep = ifelse(group == "patient", runif(n()) < 0.5, TRUE)) %>%
    ungroup() %>%
    filter(keep) %>%
    select(-keep)

  peak_dates <-
    get_peak(date = sim$to_date, group = sim$group)
  cases <- sim %>%
    group_by(group) %>%
    summarise(cases = n(), .groups = "drop")

  info <- left_join(peak_dates, cases, by = "group") %>%
    mutate(true_delta = true_delta) %>%
    rename(level = group)

  ttab_all <-
    linktree:::ttable(
      from = sim$source_group,
      to = sim$group,
      levels = o2groups.params$name
    )

  ttab_cut <- lapply(o2groups.params$name, function(level) {
    cut_tree <- cut_tree_by_date(tree = sim, cutoff_date = info$observed_peak[info$level == level])

    ttab_cut <- linktree:::ttable(
      from = cut_tree$source_group,
      to = cut_tree$group,
      levels = o2groups.params$name
    )
    return(as.data.frame(ttab_cut))
  }) %>% bind_rows(.id = "level")


  delta <- lapply(o2groups.params$name, function(level) {
    cut_tree <- cut_tree_by_date(tree = sim, cutoff_date = info$observed_peak[info$level == level])

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

  return(list(
    out = out,
    ttab_all = ttab_all,
    ttab_cut = ttab_cut
  ))
}



# Params ------------------------------------------------------------------
incubation_period <-
  distcrete::distcrete(
    "gamma",
    shape = epitrix::gamma_mucv2shapescale(5.95, (4.31 / 5.95))$shape,
    scale = epitrix::gamma_mucv2shapescale(5.95, (4.31 / 5.95))$scale,
    w = 0.5,
    interval = 1
  )

generation_time <-
  distcrete::distcrete(
    "gamma",
    shape = epitrix::gamma_mucv2shapescale(5, 0.6)$shape,
    scale =  epitrix::gamma_mucv2shapescale(5, 0.6)$scale,
    w = 0.5,
    interval = 1
  )
generation_time$d(1:100) %>% plot(y = ., type = "l")
generation_time = c(0.1, 0.2, 0.4, 0.2, 0.1) / sum(c(0.1, 0.2, 0.4, 0.2, 0.1))

f <- setNames(c(0.23, 0.77), c("HCW", "patient"))
set.seed(123)
o2groups.params <- list(
  duration = 365,
  group_n = 2,
  size = as.vector((f * 1000)),
  name = c("HCW", "patient"),
  intro_n = c(1, 1),
  r0 = c(1.32, 1.61),
  generation_time = generation_time,
  incubation_period = incubation_period$r(1000),
  gamma = c(15, 0.5)
)


set.seed(123)
results <- moVSdelta(o2groups.params)
out <- results$out
ttab_all <- results$ttab_all

out %>%
  mutate(
    lwr = if_else(metric == "mo" & is.na(lwr), est, lwr),
    upr = if_else(metric == "mo" & is.na(upr), est, upr)
  ) %>%
  mutate(across(c(est, lwr, upr, true_delta), ~ delta2gamma(.))) %>%
  ggplot(aes(
    x = level,
    y = est,
    ymin = lwr,
    ymax = upr,
    col = metric
  )) +
  geom_hline(aes(yintercept = 1), linetype = "dotted") +
  geom_point(aes(y = true_delta, shape = "truth"),
             size = 5,
             col = "#5c5e60") +
  geom_pointrange(size = 0.5, position = position_dodge(width = 0.25)) +
  scale_color_manual("metric", values =
                       c("delta" = "red", "mo" = "blue")) +
  scale_shape_manual("", values = c("truth" = 15)) +
  theme_bw() +
  labs(x = "", y = "Estimate") +
  theme(legend.position = "bottom") +
  #scale_y_continuous(limits = c(-1, 1))
  scale_y_continuous()



sjtabs <- split(results$ttab_cut, results$ttab_cut$level) %>%
  map( ~ {
    .x %>%
      uncount(Freq) %>%
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
sjtabs[[1]] # hcw
sjtabs[[2]] # patient

results$ttab_all %>%
  as.data.frame() %>%
  uncount(Freq) %>%
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
