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

set.seed(123)
array <- draw_array(
  from_col = "from_group",
  to_col = "to_group",
  levels = c("hcw", "patient"),
  trees = trees,
  draw_function = draw_delta,
  draw_args = list(f = c("hcw" = 0.5, "patient" = 0.5)),
  cutoff_dates = cutoff_dates,
  n_samples = 1000
)

CrI_df <- draw_CrI(array, c(2,3)) %>%
  reshape2::melt() %>%
  rename(metric = Var1) %>%
  pivot_wider(names_from = metric, values_from = value)

p_delta <-
  ggplot(CrI_df,
         aes(
           x = cutoff,
           y = mean,
           ymin = lwr,
           ymax = upr,
           col = level
         )) +
  geom_pointrange(position =
                    position_dodge(width = 0.5))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual("Group", values = c("orange", "purple")) +
  scale_x_discrete(
    #breaks = seq_along(cutoff_dates),
    labels = format(cutoff_dates, "%d\n%b")
    )+
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5),
    limits = c(-1, 1)
  ) +
  labs(
    x = "",
    y = "Delta"
  )+
  theme_bw()+
  theme(legend.position = "none")



cowplot::plot_grid(
  plotlist = list(p_delta, epicurve()),
  ncol = 1,
  rel_heights = c(1, 1),
  labels = "AUTO"
)

# library(gghalves)
# p_delta <-
#   reshape2::melt(array) %>%
#   ggplot() +
#   gghalves::geom_half_violin(
#     aes(
#       x = cutoff,
#       y = value,
#       fill = level,
#       split = level
#     ),
#     position = position_dodge(width = .25),
#     col = NA,
#     na.rm = TRUE,
#     alpha = 0.5
#   ) +
#   stat_summary(
#     fun.data = function(x) {
#       data.frame(y = mean(x),
#                  ymin = quantile(x, 0.025),
#                  ymax = quantile(x, 0.975))
#     },
#     geom = "pointrange",
#     aes(x = cutoff,
#         y = value,
#         col = level,
#         group = level),
#     position = position_dodge(width = .25)
#   ) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_fill_manual(values = c("hcw" = "orange", "patient" = "purple")) +
#   scale_color_manual(values = c("hcw" = "orange", "patient" = "purple")) +
#   scale_x_discrete(
#     labels = format(cutoff_dates, "%d\n%b")
#   )+
#   labs(
#   x = "",
#   y = "Delta") +
#   theme_bw() +
#   theme(legend.position = "none")



# Sensitivity -------------------------------------------------------------
#
# # No uncertainty in Transmission Tree & uncertainty in delta --------------------------------------------------------------------
# max_tree <- trees[[which.max(outbreaker_chains$post)]]
#
#
# CrI_df <- draw_array(
#   from_col = "from_group",
#   to_col = "to_group",
#   f = c("hcw" = 0.5, "patient" = 0.5),
#   trees = list(max_tree),
#   cutoff_dates = cutoff_dates
#   ) %>%
#   draw_CrI(dims = c(2,4), alpha = 0.05)
#
# p_delta <-
#   ggplot(CrI_df,
#          aes(
#            x = as.factor(cutoff),
#            y = mean,
#            ymin = lwr,
#            ymax = upr,
#            col = group
#          )) +
#   geom_pointrange(position =
#                     position_dodge(width = 0.5))+
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_color_manual(values = c("orange", "purple")) +
#   scale_x_discrete(
#     breaks = seq_along(cutoff_dates),
#     labels = format(cutoff_dates, "%d\n%b")
#   )+
#   scale_y_continuous(
#     breaks = seq(-1, 1, 0.25),
#     limits = c(-1, 1)
#   ) +
#   labs(
#     x = "",
#     y = "Delta",
#     col = "Group"
#   )+
#   theme_bw()
# p_delta
#
#
# # No uncertainty in Transmission Tree & no uncertainty in delta --------------------------------------------------------------------
#
# linktree::get_delta(
#   from = max_tree$from_group,
#   to = max_tree$to_group,
#   f = c("hcw" = 0.5, "patient" = 0.5)
# ) %>% plot()



# Sensitivity in f --------------------------------------------------------------------
# f_seq <- seq(0.1, 0.9, 0.1)
# set.seed(123)
# CrI_f <- lapply(f_seq, function(f) {
#   array <- draw_array(
#     from_col = "from_group",
#     to_col = "to_group",
#     levels = c("hcw", "patient"),
#     trees = trees,
#     draw_function = draw_delta,
#     draw_args = list(f = c("hcw" = f, "patient" = 1 - f)),
#     cutoff_dates = cutoff_dates,
#     n_samples = 1000
#   )
#   CrI_df <- draw_CrI(array, c(2,3)) %>%
#     reshape2::melt() %>%
#     rename(metric = Var1) %>%
#     mutate(f = paste0("hcw:", f, " | patient:", 1 - f))
#   return(CrI_df)
# })
# saveRDS(CrI_f, here("data", "CrI_f_delta.rds"))

# bind_rows(CrI_f) %>%
#   pivot_wider(names_from = metric, values_from = value) %>%
#   drop_na() %>%
#   ggplot(aes(
#            x = cutoff,
#            y = mean,
#            ymin = lwr,
#            ymax = upr,
#            col = level,
#            group = interaction(cutoff, level, f)
#          )) +
#   facet_wrap(~f) +
#   geom_pointrange(position = position_dodge(width = 0.25))+
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_color_manual(values = c("orange", "purple")) +
#   scale_x_discrete(
#     labels =  format(cutoff_dates, "%d\n%b")
#   )+
#   scale_y_continuous(
#     breaks = seq(-1, 1, 0.5),
#     limits = c(-1, 1)
#   ) +
#   labs(
#     x = "",
#     y = "Delta",
#     col = "Group"
#   )+
#   theme_bw()

