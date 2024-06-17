source(here::here("scripts", "helpers.R"))
load_libraries()
paper <- c("JHI2021", "eLife2022")[2] #chose which paper to run
load_data(paper)
output_path <- here::here("data", paper, "output")

epicurve()
plan(multisession, workers = availableCores())
trees <- furrr::future_map(seq_along(cutoff_dates), function(x) {
  linelist_cut <- linelist %>%
    filter(onset <= cutoff_dates[x])

  out[[x]] %>%
    filter(step > 500) %>%
    identify(ids = linelist_cut$case_id) %>%
    filter_alpha_by_kappa(1L) %>%
    get_trees(
      ids = linelist_cut$case_id,
      group = linelist_cut$group,
      date = linelist_cut$onset
    )
})


# Sensitivity in fHCW --------------------------------------------------------------------
f_seq <- seq(0.1, 0.9, 0.1)
pacman::p_load(furrr)
future::plan(future::multisession, workers = length(f_seq))
set.seed(123)
df_sensitivity <- furrr::future_map(f_seq, function(fHCW) {
  array <- draw_array(
    from_col = "from_group",
    to_col = "to_group",
    levels = c("hcw", "patient"),
    trees = trees,
    draw_function = draw_delta,
    args = list(
      f = c("hcw" = fHCW, "patient" = 1 - fHCW),
      from_id = "from",
      to_id = "to"
    ),
    n_samples = 1000
  )

  df <-
    reshape2::melt(array) %>%
    dplyr::group_by(level, cutoff) %>%
    tidyr::drop_na() %>%
    dplyr::summarise(
      mean = mean(value),
      lwr = quantile(value, 0.025),
      upr = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      cutoff_date = cutoff_dates[cutoff],
      f = paste0("HCW:", fHCW, " | patient:", 1 - fHCW),
      significant = ifelse(0 >= lwr &
                             0 <= upr, "not significant", "significant")
    )
  return(df)
}, .options = furrr_options(
  seed = TRUE,
  globals = c(
    "trees",
    "cutoff_dates",
    grep("draw_*", names(.GlobalEnv), value = TRUE),
    grep("*_formula", names(.GlobalEnv), value = TRUE)
  )
))


saveRDS(df_sensitivity,
        here::here(output_path, "df_sensitivity.rds"))
df_sensitivity <- readRDS(here::here(output_path, "df_sensitivity.rds"))




# Plot --------------------------------------------------------------------
cutoff_breaks <- cutoff_dates[seq(1, length(cutoff_dates), 4)]
day_month <- format(cutoff_breaks, "%d\n%B")
day <- format(cutoff_breaks, "%d")
dup <- which(duplicated(format(cutoff_breaks, "%b")))
day_month[dup] <- day[dup]
day_month[1] <- "11 \n    March"
day_month[12] <- "06 \nMay    "


bind_rows(df_sensitivity) %>%
  mutate(
    fHCW = str_extract(f, "(?<=HCW:)[0-9.]+") %>% as.numeric(),
    fPatient = str_extract(f, "(?<=patient:)[0-9.]+") %>% as.numeric(),
    f_label = paste0(
      '<span style="color:orange;">',
      fHCW * 100,
      # add percent here
      '%',
      '</span> / ',
      '<span style="color:purple;">',
      fPatient * 100,
      '%'
    )
  ) %>%
  ggplot(aes(
    x = cutoff_date,
    y = mean,
    ymin = lwr,
    ymax = upr,
    col = level,
    group = interaction(cutoff, level, f)
  )) +
  facet_wrap(~ f_label) +
  geom_vline(
    data = get_peak(date = linelist$onset, group = linelist$group),
    aes(xintercept = observed_peak,
        col = group),
    linetype = "dotted",
    show.legend = FALSE,
    linewidth = 0.65
  ) +
  geom_errorbar(position = position_dodge(width = 0.5),
                width = 0.5,
                linewidth = 0.5) +
  geom_point(
    aes(shape = significant),
    position = position_dodge(width = 0.5),
    size = 2.5,
    fill = "white",
    alpha = 0.9
  ) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  scale_shape_manual(values =
                       c("significant" = 16, "not significant" = 21)) +
  theme_noso() +
  labs(x = "",
       y = "",
       col = "",
       shape = "") +
  theme(
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      size = 1,
      linetype = "solid"
    ),
    strip.text = ggtext::element_markdown(
      # Use element_markdown for HTML formatted text
      angle = 0,
      face = "bold",
      size = 12,
      margin = margin(2, 0, .1, 0, "pt")
    ),
    legend.margin = margin(-10, 0, 0, 0, "pt")
  )+
  scale_x_date(
    # break 1 every other 2 cutoff_date
    breaks = cutoff_breaks,
    labels = day_month,
    expand = c(0.01, 0.5)
  )


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
