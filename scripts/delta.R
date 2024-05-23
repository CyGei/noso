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




# Sensitivity in fHCW --------------------------------------------------------------------
f_seq <- seq(0.1, 0.9, 0.1)
pacman::p_load(furrr)
future::plan(list(
  future::tweak(future::multisession,
                workers = length(f_seq)),
  future::tweak(future::multisession,
                workers = future::availableCores() - (length(f_seq) + 2))
))

set.seed(123)
CrI_f <- furrr::future_map(f_seq, function(fHCW) {
  array <- draw_array(
    from_col = "from_group",
    to_col = "to_group",
    levels = c("hcw", "patient"),
    trees = trees,
    draw_function = draw_delta,
    args = list(f = c(
      "hcw" = fHCW, "patient" = 1 - fHCW
    )),
    cutoff_dates = cutoff_dates,
    n_samples = 1000
  )
  CrI_df <- draw_CrI(array, c(2, 3)) %>%
    reshape2::melt() %>%
    dplyr::rename(metric = Var1) %>%
    dplyr::mutate(f = paste0("hcw:", fHCW, " | patient:", 1 - fHCW))
  return(CrI_df)
},
.options = furrr_options(
  seed = TRUE,
  globals = c(
    "trees",
    "cutoff_dates",
    "cut_tree_by_date",
    grep("draw_*", names(.GlobalEnv),
         value = TRUE),
    grep("*_formula", names(.GlobalEnv),
         value = TRUE)
  )
))
saveRDS(CrI_f, here("data", "CrI_f_delta.rds"))

#CrI_f <- readRDS(here("data", "CrI_f_delta.rds"))

bind_rows(CrI_f) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  drop_na() %>%
  ggplot(aes(
           x = cutoff,
           y = mean,
           ymin = lwr,
           ymax = upr,
           col = level,
           group = interaction(cutoff, level, f)
         )) +
  facet_wrap(~f) +
  geom_pointrange(
    size = 0.25,
    position = position_dodge(width = 0.25))+
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(
    data = tibble(
      level = c("hcw", "patient"),
      cutoff = c(6, 4)
    ),
    aes(xintercept = cutoff, col = level),
    linetype = "dotted",
    linewidth = 1
  )+
  scale_color_manual(values = c("orange", "purple")) +
  scale_x_discrete(
    labels =  format(cutoff_dates, "%d\n%b")
  )+
  scale_y_continuous(
    breaks = seq(-1, 1, 0.5),
    limits = c(-1, 1)
  ) +
  labs(
    x = "",
    y = "Delta",
    col = "Group"
  )+
  theme_bw()



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


