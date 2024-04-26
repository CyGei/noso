# The purpose of this script is to explore the relationship between gamma and mo's ratio.
# We will simulate multiple outbreaks with different gamma values and calculate the mo's ratio for each outbreak.
# We will then plot the relationship between gamma and mo's ratio.


# Scenarios ---------------------------------------------------------------
pacman::p_load(o2groups, furrr, dplyr, here)
duration <- 365
group_n <- 2
size <- c(150, 30)
name <- c("HCW", "patient")
intro_n <- c(1, 1)
r0 <- c(2, 2)
generation_time <- c(0, 0.1, 0.2, 0.4, 0.2, 0.1, 0)
incubation_period <- sample(1:14, sum(size), replace = TRUE)
# delta <- sort(c(0, seq(-0.9, 0.9, 0.2))) %>% round(2)
# gamma <- split(x = rep(linktree::delta2gamma(delta), each = 2),
#                f = rep(LETTERS[1:length(delta)], each = 2))
# n_scenarios <- length(gamma)
n_scenarios <- 100

# Generate a list of scenarios with different gamma values
set.seed(123)
scenarios <- lapply(1:n_scenarios, function(i) {
  list(
    duration = duration,
    group_n = group_n,
    size = size,
    name = name,
    gamma = linktree::delta2gamma(runif(2, -1, 1)),
    intro_n = intro_n,
    r0 = r0,
    generation_time = generation_time,
    incubation_period = incubation_period
  )
})

scenario_df <- lapply(scenarios, function(scenario) {
  data.frame(
    level = scenario$name,
    input_gamma = scenario$gamma
  )
}) %>%
  bind_rows(.id = "scenario")

# Simulations -------------------------------------------------------------
# Run simulate_groups for each scenario
set.seed(123)
timing <- system.time(
  simulations <- lapply(scenarios, function(scenario) {
    o2groups::simulate_groups_furrr(
      sim_n = 100,
      duration = scenario$duration,
      group_n = scenario$group_n,
      size = scenario$size,
      name = scenario$name,
      gamma = scenario$gamma,
      intro_n = scenario$intro_n,
      r0 = scenario$r0,
      generation_time = scenario$generation_time,
      incubation_period = scenario$incubation_period
    )
  })
)
timing / 60
dir.create(here::here("data"), showWarnings = FALSE)
saveRDS(simulations, here::here("data", "simulations.rds"))
simulations <- readRDS(here::here("data", "simulations.rds"))

# Analysis ----------------------------------------------------------------
source(here::here("scripts", "helpers.R"))
load_libraries()
plan(multisession, workers = availableCores() - 2)
peaks <- furrr::future_map(1:length(scenarios), function(i) {
  scenario <- i
  map(1:length(simulations[[i]]), function(j) {
    sim_n <- j
    x <- simulations[[i]][[j]]
    peak <- get_peak(date = x$date_onset, group = x$group) %>%
      mutate(scenario = as.character(scenario), sim_n = sim_n) %>%
      rename(level = group)
  })
}, .options = furrr_options(seed = TRUE)) %>%
  bind_rows()

# Peaks
left_join(bind_rows(peaks),
          scenario_df,
          by = c("scenario", "level")) %>%
  ggplot(aes(x = observed_peak, fill = scenario)) +
  facet_wrap(~level, ncol = 1) +
  geom_density(alpha = 0.5, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0, 365, 7)) +
  labs(x = "Peak date", y = "Density") +
  theme_bw()
#~day 20

# total case count
simulations %>%
  lapply(bind_rows, .id = "sim_n") %>%
  bind_rows(.id = "scenario") %>%
  group_by(scenario, sim_n) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = n)) +
  geom_histogram(aes(fill = scenario),
                 binwidth = 1,
                 show.legend = FALSE) +
  labs(x = "Total case count", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 200, 10)) +
  theme_bw()


# reformat data to allow cut_tree_by_date()
simulations <- lapply(simulations, function(scenario) {
  lapply(scenario, function(x) {
    x %>%
      mutate(from_date = x$date_onset[match(source, x$id)]) %>%
      rename(to_date = date_onset)
  })
})

# Metrics -----------------------------------------------------------------
metrics <- list(
  "gamma" = draw_gamma,
  "Mo" = draw_mo
)

arrays <- lapply(seq_along(scenarios), function(i) {
  set.seed(123)
  array <- lapply(
    names(metrics),
    function(metric_name) {
      draw_array(
        from_col = "source_group",
        to_col = "group",
        levels = c("HCW", "patient"),
        trees = simulations[[i]], # simulation corresponding to the i-th scenario
        draw_function = metrics[[metric_name]],
        args = list(
          f = c("HCW" = 0.5, "patient" = 0.5),
          from_id = "source",
          to_id = "id",
          diag = TRUE
        ),
        cutoff_dates = 20,
        n_samples = ifelse(metric_name == "Mo", 1, 1000)
      )
    }
  )
})
saveRDS(arrays, here::here("data", "sim_arrays.rds"))
arrays <- readRDS(here::here("data", "sim_arrays.rds"))


CrI_list <- lapply(arrays, function(scenario) {
  lapply(scenario, function(x)
    draw_CrI(x, c(2, 3)))
})
df <- lapply(CrI_list, function(scenario) {
  lapply(scenario, function(x) {
    x %>%
      reshape2::melt() %>%
      pivot_wider(names_from = Var1, values_from = value)
  }) %>%
    bind_rows(.id = "metric") %>%
    mutate(metric = factor(
      metric,
      levels = 1:length(metrics),
      labels = names(metrics)
    ),
    reference = 1)
}) %>%
  bind_rows(.id = "scenario")

df <- left_join(x = df, y = scenario_df, by = c("scenario", "level"))

# True gamma vs estimated gamma
df %>%
  filter(level == "HCW" & metric == "gamma") %>%
  ggplot(aes(x = input_gamma,
             y = mean,
             ymin = lwr,
             ymax = upr)) +
  geom_pointrange() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "true gamma", y = "estimated gamma")+
  coord_fixed(ylim = c(0,10),
                  xlim = c(0,10)) +
  scale_x_continuous(breaks = seq(0, 100, 1)) +
  scale_y_continuous(breaks = seq(0, 100, 1))

df %>%
  ggplot(aes(
    x = input_gamma,
    y = mean,
    ymin = lwr,
    ymax = upr,
    col = metric
  )) +
  facet_wrap(~ level) +
  geom_pointrange(size = 0.1) +
  geom_hline(aes(yintercept = reference), linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_color_manual("Metric",
                     values =
                       c("gamma" = "red", "Mo" = "blue")) +
  coord_fixed(ylim = c(0,10),
              xlim = c(0,10)) +
  scale_x_continuous(breaks = seq(0, 100, 1)) +
  scale_y_continuous(breaks = seq(0, 100, 1))+
  theme_bw() +
  labs(x = "True Gamma",
       y = "Estimate")


