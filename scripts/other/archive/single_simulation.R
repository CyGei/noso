source(here::here("scripts", "helpers.R"))
load_libraries()
pacman::p_load(o2groups, furrr, dplyr, here)
duration <- 365
group_n <- 2
size <- c(150, 30)
name <- c("HCW", "patient")
intro_n <- c(1, 1)
r0 <- c(2, 2)
generation_time <- c(0, 0.1, 0.2, 0.4, 0.2, 0.1, 0)
incubation_period <- sample(1:14, sum(size), replace = TRUE)

input_gamma <- c(10, 1)
f <- size / sum(size)


set.seed(123)
simulations <- o2groups::simulate_groups_furrr(
  sim_n = 100,
  duration = duration,
  group_n = group_n,
  size = size,
  name = name,
  gamma = input_gamma,
  intro_n = intro_n,
  r0 = r0,
  generation_time = generation_time,
  incubation_period = incubation_period
)
# reformat data to allow cut_tree_by_date()
simulations <-
  lapply(simulations, function(x) {
    x %>%
      mutate(from_date = x$date_onset[match(source, x$id)]) %>%
      rename(to_date = date_onset)
  })

# estimate peak date
peak_date <- bind_rows(simulations, .id = "sim_n") %>%
  group_by(sim_n) %>%
  count(to_date) %>%
  filter(n == max(n)) %>%
  slice(1) %>%
  pull(to_date) %>%
  mean() - 1
# lapply(simulations, function(sim) {
#   get_peak(date = sim$to_date,
#            group = sim$group)
# }) %>%
#   bind_rows(.id = "sim_n") %>%
#   summarise(mean_peak = mean(observed_peak))
#peak day 20

#average number of cases
lapply(simulations, function(sim){
  sim %>%
    group_by(group) %>%
    summarise(n = n())
}) %>%
  bind_rows(.id = "sim_n") %>%
  ggplot(aes(x = n, fill = group)) +
  geom_density(alpha = 0.5, show.legend = FALSE) +
  scale_x_continuous(breaks = seq(0, 200, 20))


# Metrics -----------------------------------------------------------------
metrics <- list(
  "gamma" = draw_gamma,
  "Mo" = draw_mo
)

set.seed(123)
array <- lapply(
  names(metrics),
  function(metric_name) {
    draw_array(
      from_col = "source_group",
      to_col = "group",
      levels = c("HCW", "patient"),
      trees = simulations,
      draw_function = metrics[[metric_name]],
      args = list(
        #HCW to patient ratio is 9:30
        f = c("HCW" = f[1], "patient" = f[2]),
        from_id = "source",
        to_id = "id",
        diag = TRUE
      ),
      cutoff_dates = as.integer(peak_date),
      n_samples = ifelse(metric_name == "Mo", 1, 1000)
    )
  }
)

CrI <- lapply(array, \(x) {draw_CrI(x, 2) %>%
                    reshape2::melt() %>%
                    pivot_wider(names_from = Var1, values_from = value) %>%
                    mutate(input_gamma = input_gamma)
  }) %>%
  bind_rows(.id = "metric") %>%
  mutate(metric = factor(metric, labels = c("gamma", "Mo")))


CrI %>%
  ggplot(aes(
    x = level,
    y = mean,
    ymin = lwr,
    ymax = upr,
    col = metric
  )) +
  geom_pointrange(size = 0.5, position = position_dodge(width = 0.1)) +
  geom_hline(aes(yintercept = input_gamma, linetype = "true gamma")) +
  geom_hline(aes(yintercept = 1, linetype = "null value")) +
  scale_linetype_manual("Reference",
                        values = c("true gamma" = "solid", "null value" = "dashed"),
                        breaks = c("true gamma", "null value")) +
  scale_color_manual("Metric",
                     values =
                       c("gamma" = "red", "Mo" = "blue")) +
  theme_bw() +
  labs(x = "",
       y = "Estimate")



#Mo's expected:
#mo_formula
n_cases <- c(120, 30)
names(n_cases) <- c("HCW", "patient")
infector <- "patient"
infectee <- "patient"
expected_numerator <- n_cases[[infector]] - 1
expected_denominator <- sum(n_cases) - 1
expected_ratio <- expected_numerator / expected_denominator
expected_ratio
# we expect 0.8 of the infectors to be HCW
# we expect 0.2 of the infectors to be patients

#Mo's observed:
#mo_formula
lapply(simulations, function(sim){
  ttab <- linktree:::ttable(from = sim$source_group, to = sim$group, levels = c("HCW", "patient"))
  levels <- rownames(ttab)
  # Observed
  observed_numerator <- ttab[infector, infectee]
  observed_denominator <- sum(ttab[, infectee])
  observed_ratio <- observed_numerator / observed_denominator
  return(observed_ratio)
}) %>% as.numeric() %>% na.omit() %>%  mean() #/ expected_ratio


# Gamma
#gamma_formula
linktree::get_gamma
lapply(simulations, function(sim){
  get_gamma(from = sim$source_group, to = sim$group, f = c("HCW" = 0.23, "patient" = 0.77))
}) %>%
  bind_rows(.id = "sim_n") %>%
  group_by(group) %>%
  filter(is.finite(est)) %>%
  summarise(mean = mean(est, na.rm = TRUE))


#for each simulation count the average number of group to group transmissions
ttab <- lapply(simulations, function(sim){
  ttab <- linktree:::ttable(from = sim$source_group, to = sim$group, levels = c("HCW", "patient"))
  as.data.frame(ttab)
}) %>%
  bind_rows(.id = "sim_n") %>%
  group_by(from, to) %>%
  summarise(mean = mean(Freq, na.rm = TRUE))
ttab
n_cases
#pi
n_cases <- c(120, 30)
names(n_cases) <- c("HCW", "patient")
infector <- "HCW"
infectee <- "HCW"
pi <- ttab$mean[ttab$from == infector & ttab$to == infectee] / sum(ttab$mean[ttab$from == infector])
g = (pi * (1 - 0.5)) / (0.5 * (1 - pi))
g


# array %>%
#   reshape2::melt() %>%
#   rename(metric = L1) %>%
#   mutate(metric = ifelse(metric == 1, "gamma", "Mo")) %>%
#   group_by(metric, level) %>%
#   drop_na() %>%
#   filter(is.finite(value)) %>%
#   summarise(mean = mean(value),
#             median = median(value),
#             lwr = quantile(value, 0.025),
#             upr = quantile(value, 0.975)) %>%
#   mutate(input_gamma = ifelse(level == "HCW", input_gamma[1], input_gamma[2])) %>%
#   ggplot(aes(
#     x = input_gamma,
#     y = median,
#     ymin = lwr,
#     ymax = upr,
#     col = metric
#   )) +
#   facet_wrap(~ level, scales = "free") +
#   geom_pointrange(size = 0.1) +
#   geom_hline(aes(yintercept = 1), linetype = "dotted") +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   scale_color_manual("Metric",
#                      values =
#                        c("gamma" = "red", "Mo" = "blue")) +
#   theme_bw() +
#   labs(x = "Input Gamma",
#        y = "Estimate")
