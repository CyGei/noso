load_helpers()
load_libraries()
load_paper()

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
  geom_line(linetype = "dotted") +
  #geom_segment(aes(xend = x, yend = 0)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Serial Interval (days)", y = "Cumulative Probability") +
  scale_y_continuous(
    breaks = seq(0, 1, 0.1),
    labels = scales::percent_format(accuracy = 1)) +
  scale_x_continuous(
    limits = c(-25, 25),
    breaks = seq(-50, 50, 5)) +
  theme_noso(fill = F, col = F, date = F)


cowplot::plot_grid(
  p1,
  p2 + labs(y = ""),
  ncol = 2,
  align = "h",
  labels = c("A", "B")
)
