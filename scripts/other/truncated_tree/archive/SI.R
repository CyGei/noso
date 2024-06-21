source(here::here("scripts", "helpers.R"))
load_libraries()
paper <- c("JHI2021", "eLife2022")[1] #chose which paper to run
load_data(paper)
output_path <- here::here("data", paper, "output")
epicurve()

#out <- readRDS(here("data", paper, "input", "out.rds"))
# out <- readRDS(here::here("data", paper, "output", "out_cyril4.rds"))

#final outbreaker run
trees <-
  out[[length(out)]] %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  get_trees(
    ids = linelist$case_id,
    group = linelist$group,
    date = linelist$onset,
    kappa = TRUE
  ) %>%
  map(~ .x %>% filter(kappa == 1))

# Serial interval ----------------------------------------------------------------
si_ecdf <- lapply(trees, function(x){
  si <- x %>%
    mutate(
      serial_interval = as.numeric(to_date - from_date)
    ) %>%
    pull(serial_interval)

  f <- ecdf(si)
  return(
    tibble(
      x = -50:50,
      y = f(-50:50)
    )
  )
}) %>% bind_rows(.id = "tree")

si_ecdf %>%
  group_by(x) %>%
  summarise(
    mean = mean(y, na.rm = T),
    lwr = quantile(y, 0.025, na.rm = T),
    upr = quantile(y, 0.975, na.rm = T)
  ) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "grey",
              alpha = 0.5) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Serial Interval (days)", y = "Cumulative Probability") +
  theme_noso(fill = F, col = F, date = F)




# At peak# ----------------------------------------------------------------

trunc_trees <- readRDS(here("data/other/truncated_tree", "trunc_trees.rds"))
realtime_trees <- readRDS(here("data/other/truncated_tree", "realtime_trees.rds"))

peak <- as.Date("2020-04-02")
peak_idx <- which(cutoff_dates == peak)


si_ecdf <- lapply(trunc_trees[[peak_idx]], function(x){
  si <- x %>%
    mutate(
      serial_interval = as.numeric(to_date - from_date)
    ) %>%
    pull(serial_interval)

  f <- ecdf(si)
  return(
    tibble(
      x = -50:50,
      y = f(-50:50)
    )
  )
}) %>% bind_rows(.id = "tree")

si_ecdf %>%
  group_by(x) %>%
  summarise(
    mean = mean(y, na.rm = T),
    lwr = quantile(y, 0.025, na.rm = T),
    upr = quantile(y, 0.975, na.rm = T)
  ) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "grey",
              alpha = 0.5) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Serial Interval (days)", y = "Cumulative Probability") +
  theme_noso(fill = F, col = F, date = F) +
  scale_x_continuous(limits = c(-20, 30), breaks = seq(-20, 30, 5))



si_ecdf <- lapply(realtime_trees[[peak_idx]], function(x){
  si <- x %>%
    mutate(
      serial_interval = as.numeric(to_date - from_date)
    ) %>%
    pull(serial_interval)

  f <- ecdf(si)
  return(
    tibble(
      x = -50:50,
      y = f(-50:50)
    )
  )
}) %>% bind_rows(.id = "tree")

si_ecdf %>%
  group_by(x) %>%
  summarise(
    mean = mean(y, na.rm = T),
    lwr = quantile(y, 0.025, na.rm = T),
    upr = quantile(y, 0.975, na.rm = T)
  ) %>%
  ggplot(aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "grey",
              alpha = 0.5) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Serial Interval (days)", y = "Cumulative Probability") +
  theme_noso(fill = F, col = F, date = F) +
  scale_x_continuous(limits = c(-20, 30), breaks = seq(-20, 30, 5))
