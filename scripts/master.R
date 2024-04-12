source(here::here("scripts", "helpers.R"))
load_libraries()
load_data()
linelist$group <- ifelse(grepl("^C", linelist$case_id), "patient", "hcw")
outbreaker_chains <- outbreaker_chains %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id)
outbreaker_chains <- filter_alpha_by_kappa(outbreaker_chains, 1)

trees <- get_trees(
  out = outbreaker_chains,
  ids = linelist$case_id,
  group = linelist$group,
  date = linelist$onset_inferred
)

cutoff_dates <- unique(c(seq(min(as.Date(linelist$onset_inferred)),
                             max(as.Date(linelist$onset_inferred)),
                             by = 7),
                         max(linelist$onset_inferred)))


metrics <- list(
  "gamma" = draw_gamma,
  "delta" = draw_delta,
  "Mo" = draw_mo,
  "R0" = draw_R0
)



set.seed(123)
array <- lapply(names(metrics),
                function(metric_name) {
                  draw_array(
                    from_col = "from_group",
                    to_col = "to_group",
                    levels = c("hcw", "patient"),
                    trees = trees,
                    draw_function = metrics[[metric_name]],
                    args = list(
                      f = c("hcw" = 0.5, "patient" = 0.5),
                      from_id = "from",
                      to_id = "to",
                      diag = TRUE
                    ),
                    cutoff_dates = cutoff_dates,
                    n_samples = ifelse(metric_name == "Mo", 1, 1000)
                  )
                })
#saveRDS(array, here::here("data", "array.rds"))
CrI <- lapply(array, function(x) draw_CrI(x, c(2,3)))

df <- lapply(CrI, function(x) {
  x %>%
    reshape2::melt() %>%
    pivot_wider(names_from = Var1, values_from = value)
}) %>%
  bind_rows(.id = "metric") %>% #
  mutate(metric = factor(metric, levels = 1:4, labels = names(metrics)),
         target = case_when(
           metric == "gamma" ~ 1,
           metric == "delta" ~ 0,
           metric == "Mo" ~ 1,
           metric == "R0" ~ 1
         ))


p_metrics <- ggplot(df,
                    aes(
                      x = cutoff,
                      y = mean,
                      ymin = lwr,
                      ymax = upr,
                      col = level
                    )) +
  facet_wrap(~ metric,
             scales = "free_y",
             ncol = 1) +
  geom_pointrange(position =
                    position_dodge(width = 0.5)) +
  geom_hline(aes(yintercept = target), linetype = "dashed") +
  scale_color_manual("Group", values = c("hcw" = "orange", "patient" = "purple"))+
  scale_x_discrete(labels = format(cutoff_dates, "%d\n%b")) +
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "delta" ~ scale_y_continuous(limit = c(-1, 1),
                                             breaks = seq(-1, 1, 0.5)),
      metric == "gamma" ~ scale_y_continuous(limit = c(0, 20),
                                             breaks = seq(0, 20, 5)),
      metric == "Mo" ~ scale_y_continuous(),
      metric == "R0" ~ scale_y_continuous(limit = c(0, 5),
                                          breaks = seq(0, 5, 1))
    )
  )+
  labs(x = "",
       y = "metric value") +
  theme_bw() +
  theme(legend.position = "none")



cowplot::plot_grid(p_metrics,
                   epicurve(),
                   ncol = 1,
                   rel_heights = c(1, 0.4),
                   labels = "AUTO"
)


# R0 ----------------------------------------------------------------------
set.seed(123)
array <- draw_array(
  from_col = "from_group",
  to_col = "to_group",
  levels = c("hcw", "patient"),
  trees = trees,
  draw_function = draw_R0,
  args = list(
    f = c("hcw" = 0.5, "patient" = 0.5),
    from_id = "from",
    to_id = "to"
  ),
  cutoff_dates = "2020-03-20",
  n_samples = 1000
)

CrI <- draw_CrI(array, c(2,3)) %>%
  reshape2::melt() %>%
  pivot_wider(names_from = Var1, values_from = value)

ggplot(CrI,
       aes(
         x = cutoff,
         y = mean,
         ymin = lwr,
         ymax = upr,
         col = level
       )) +
  geom_pointrange(position =
                    position_dodge(width = 0.5))+
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual("Group", values = c("orange", "purple")) +
  scale_x_discrete(
    #breaks = seq_along(cutoff_dates),
    labels = format(cutoff_dates, "%d\n%b")
  )+
  labs(
    x = "",
    y = "Delta"
  )+
  theme_bw()+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 20))


array %>%
  reshape2::melt() %>%
  ggplot(aes(x = level, y = value)) +
  geom_violin(adjust=2, trim = TRUE)+
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "R0") +
  theme_bw()+
  scale_y_continuous(limits = c(0, 10))

