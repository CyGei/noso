explorer <- function(o2groups_params){

  sims <- furrr::future_map(1:100, function(i) {
    o2groups::simulate_groups(
      duration = o2groups_params$duration,
      group_n = o2groups_params$group_n,
      size = o2groups_params$size,
      name = o2groups_params$name,
      gamma = o2groups_params$gamma,
      intro_n = o2groups_params$intro_n,
      r0 = o2groups_params$r0,
      generation_time = o2groups_params$generation_time,
      incubation_period = o2groups_params$incubation_period
    )
  }, .options = furrr_options(seed = TRUE)) %>%
    bind_rows(.id = "sim")

  # No of Cases
  cases <- sims %>%
    group_by(sim, group) %>%
    summarise(case = n(), .groups = "drop") %>%
    group_by(group) %>%
    summarise(mean = mean(case), .groups = "drop")
  cat("average cases: ", cases$mean, "\n")

  # Average Peak Date
  peaks <- sims %>%
    group_by(sim, group, date_onset) %>%
    summarise(case = n(), .groups = "drop") %>%
    group_by(sim, group) %>%
    arrange(sim, desc(case)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    group_by(group) %>%
    summarise(mean = mean(date_onset), .groups = "drop")
  cat("average peak date: ", peaks$mean, "\n")


  #EpiCurve
  p <- sims %>%
    ggplot(aes(date_onset, fill = sim)) +
    geom_bar(position = "stack", show.legend = FALSE) +
    geom_vline(data = peaks,
               aes(xintercept = mean, group = group),
               linetype = "dashed") +
    facet_wrap( ~ group, ncol = 1, scales = "free_y") +
    theme_bw() +
    labs(title = "Epidemic Curve", x = "Date", y = "Cases")
  print(p)


  # ttable: total outbreak
  ttab <- split(sims, sims$sim) %>%
    future_map( ~ {
      linktree:::ttable(
        from = .x$source_group,
        to = .x$group,
        levels = c("HCW", "patient")
      ) %>%
        as_tibble()
    }) %>%
    bind_rows() %>%
    group_by(from, to) %>%
    summarise(mean = round(mean(n)), .groups = "drop")

  #ttable : HCW peak
  ttab_HCW <- split(sims, sims$sim) %>%
    future_map( ~ {
      peak <- .x %>%
        filter(group == "HCW") %>%
        group_by(date_onset) %>%
        summarise(n = n()) %>%
        arrange(desc(n)) %>%
        slice_head(n = 1) %>%
        pull(date_onset)

      tree <- .x %>% filter(date_onset <= peak)

      linktree:::ttable(
        from = tree$source_group,
        to = tree$group,
        levels = c("HCW", "patient")
      ) %>%
        as_tibble()
    }) %>%
    bind_rows() %>%
    group_by(from, to) %>%
    summarise(mean = round(mean(n)), .groups = "drop")

  #ttable : patient peak
  ttab_patient <- split(sims, sims$sim) %>%
    future_map( ~ {
      peak <- .x %>%
        filter(group == "patient") %>%
        group_by(date_onset) %>%
        summarise(n = n()) %>%
        arrange(desc(n)) %>%
        slice_head(n = 1) %>%
        pull(date_onset)

      tree <- .x %>% filter(date_onset <= peak)

      linktree:::ttable(
        from = tree$source_group,
        to = tree$group,
        levels = c("HCW", "patient")
      ) %>%
        as_tibble()
    }) %>%
    bind_rows() %>%
    group_by(from, to) %>%
    summarise(mean = round(mean(n)), .groups = "drop")


  return(list(
    cases = cases,
    peaks = peaks,
    p = p,
    ttab = ttab,
    ttab_HCW = ttab_HCW,
    ttab_patient = ttab_patient
  ))

}

get_mo <- function(ttab_df){
  ttab <- ttab_df %>%
    uncount(mean) %>%
    select(from, to) %>%
    xtabs(~ from + to, data = .)

  prop_received <- prop.table(ttab, margin = 2)
  prop_cases <- colSums(ttab) / sum(ttab)

  return(diag(prop_received) / prop_cases)
}


pacman::p_load(tidyverse, linktree, o2groups, furrr)
plan(multisession, workers = 20)


incubation_period <-
  distcrete::distcrete(
    "gamma",
    shape = epitrix::gamma_mucv2shapescale(5.95, (4.31/5.95))$shape,
    scale = epitrix::gamma_mucv2shapescale(5.95, (4.31/5.95))$scale,
    w = 0.5,
    interval = 1
  )

generation_time = c(0.1, 0.2, 0.4, 0.2, 0.1)/sum(c(0.1, 0.2, 0.4, 0.2, 0.1))


f = c(0.5, 0.5)
popsize = 1000
r0 = c(2, 2)
gamma = c(1,1)

o2groups_params <- list(
  duration = 365,
  group_n = 2,
  f = f,
  size = f * popsize,
  name = c("HCW", "patient"),
  intro_n = c(1,1),
  r0 = r0,
  generation_time = generation_time,
  incubation_period = incubation_period$r(1000),
  gamma = gamma
)
set.seed(123)
explo = explorer(o2groups_params)

explo$ttab %>%
  uncount(mean) %>%
  select(from, to) %>%
  sjPlot::sjtab(fun = "xtab",
                var.labels = c("infector", "infectee"),
                show.row.prc=T,
                show.col.prc=T,
                show.summary=T,
                show.exp=T,
                show.legend=T)

explo$ttab_HCW %>%
  uncount(mean) %>%
  select(from, to) %>%
  sjPlot::sjtab(fun = "xtab",
                var.labels = c("infector", "infectee"),
                show.row.prc=T,
                show.col.prc=T,
                show.summary=T,
                show.exp=T,
                show.legend=T)



explo$ttab_patient %>%
  uncount(mean) %>%
  select(from, to) %>%
  sjPlot::sjtab(fun = "xtab",
                var.labels = c("infector", "infectee"),
                show.row.prc=T,
                show.col.prc=T,
                show.summary=T,
                show.exp=T,
                show.legend=T)


#MO
get_mo(explo$ttab) %>% round(2)


# Gamma for HCW
explo$ttab_HCW %>%
  uncount(mean) %>%
  select(from, to) %>%
  {
    linktree::get_gamma(
      from = .$from,
      to = .$to,
      f = setNames(o2groups_params$f, o2groups_params$name),
      alpha = 0.05
    )
  } %>%
  mutate(across(where(is.numeric),\(x) round(x, 2)))
# Gamma for patient
explo$ttab_patient %>%
  uncount(mean) %>%
  select(from, to) %>%
  {
    linktree::get_gamma(
      from = .$from,
      to = .$to,
      f = setNames(o2groups_params$f, o2groups_params$name),
      alpha = 0.05
    )
  } %>%
  mutate(across(where(is.numeric),\(x) round(x, 2)))


linktree:::gamma_formula(0.67, 0.5)

# From actual results -----------------------------------------------------

# gammaf <- data.frame(
#   fHCW = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
#   hcw = c(33.77, 15.13, 8.90, 5.79, 3.91, 2.65, 1.74, 1.05, 0.49),
#   patient = c(0.21, 0.42, 0.67, 0.97, 1.37, 1.94, 2.85, 4.59, 9.63)
# )
#
# params_list <- lapply(1:nrow(gammaf), function(i){
#   list(
#     duration = 70,
#     group_n = 2,
#     f = as.vector(unlist(c(gammaf[i, "fHCW"], 1 - gammaf[i, "fHCW"]))),
#     name = c("HCW", "patient"),
#     intro_n = c(1,1),
#     r0 = c(1.29, 1.61),
#     generation_time = generation_time,
#     incubation_period = incubation_period$r(1000),
#     gamma = as.vector(unlist((gammaf[i, c("hcw", "patient")])))
#   )
# })
#
# o2groups_params <- params_list[[1]]
# o2groups_params$size <- o2groups_params$f * 1000


