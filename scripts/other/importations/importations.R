source(here::here("scripts", "helpers.R"))
source(here::here("scripts/other/importations", "helpers.R"))
load_libraries()
paper <- c("JHI2021", "eLife2022")[1] #chose which paper to run
load_data(paper)
output_path <- here::here("data", paper, "output")
epicurve()
out <- out[[length(out)]]
# Total Introductions
alpha_cols <- grep("alpha_", colnames(out))
sum(is.na(out[sample(1:nrow(out), 1), alpha_cols]))


trees <- out %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  get_trees(
    ids = linelist$case_id,
    group = linelist$group,
    date = linelist$onset,
    kappa = TRUE,
    t_inf = TRUE
  )


# Estimate the average number of imports per time window ------------
window <- ifelse(paper == "JHI2021", 5, 7)
start <- min(linelist$onset)
end <- max(linelist$onset)
cutoff_breaks <- c(seq.Date(start, end, by = window), end) %>% unique()

plan(multisession, workers = length(cutoff_breaks))
# Compute imports by group
imports_df_group <- future_map_dfr(seq_len(length(cutoff_breaks) - 1),
                                   ~ process_window(.x, trees, cutoff_breaks, TRUE, levels = unique(linelist$group)))

# Compute overall imports
imports_df_overall <- future_map_dfr(
  seq_len(length(cutoff_breaks) - 1),
  ~ process_window(.x, trees, cutoff_breaks, FALSE, levels = unique(linelist$group))
) %>%
  mutate(to_group = "global")

imports_df <- bind_rows(imports_df_group, imports_df_overall)
imports_summary <- imports_df %>%
  group_by(to_group, window_median, window_start, window_end) %>%
  summarise(
    mean = mean(n_imports),
    lwr = quantile(n_imports, 0.025),
    upr = quantile(n_imports, 0.975),
    .groups = "drop"
  )


p_imports <-
imports_df %>%
  ggplot(aes(
    x = window_median,
    y = n_imports,
    fill = to_group,
    group = interaction(to_group, window_median)
  )) +
  # geom_errorbarh(
  #   data = imports_summary %>%
  #     filter(to_group == "global"),
  #   aes(xmin = window_start,
  #       xmax = window_end,
  #       y = mean,
  #       height = 0),
  #   size = 0.25
  # ) +
  geom_violin(position = position_dodge(width = 3),
              col = NA,
              scale = "width") +
  geom_pointrange(
    data = imports_summary,
    aes(
      x = window_median,
      y = mean,
      ymin = lwr,
      ymax = upr,
    ),
    position = position_dodge(width = 3),
    pch = 21,
    stroke = 0.5,
    size = 0.75
  ) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  theme_noso(day_break = window, date = TRUE) +
  scale_fill_manual(
    values = c(
      global = "#767676",
      hcw = "#FF9300",
      patient = "#d800d1"
    ),
    labels = c(
      "global" = "Global",
      "hcw" = "HCW",
      "patient" = "Patient"
    )
  ) +
  coord_cartesian(ylim = c(0, 8)) +
  labs(x = "", y = "Importations",
       fill = "")


cowplot::plot_grid(
  p_imports + theme(legend.position = "none"),
  epicurve(day_break = window),
  ncol = 1,
  rel_heights = c(2, 1.5),
  align = "v",
  labels = "AUTO"
)
intro_summary %>% pull(mean) %>% sum()


# #eLife2022: no introductions in the 1st 2 weeks!
# ids_to_check <- linelist %>%
#   arrange(onset) %>%
#   filter(onset < "2020-03-19") %>%
#   pull(case_id)
#
# bind_rows(trees) %>%
#   filter(to %in% ids_to_check) %>%
#   mutate(SI = to_date - from_date) %>%
#   ggplot(aes(x = SI)) +
#   geom_histogram(binwidth = 1)


# With posterior support --------------------------------------------------

# How many times is each case_id identified as an introduction?
bind_rows(trees) %>%
  filter(is.na(from)) %>%
  group_by(to) %>%
  summarise(n = n(),
            support = n / length(trees),
            .groups = "drop") %>%
  filter(support > 0.1) %>%
  arrange(to)
