# The purpose of this script is to estimate the reproduction number
# using transmission chain data.
# For each type of transmission.
load_helpers()
load_paper()
epicurve()
trees <- out[[length(out)]] %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  get_trees(
    out = .,
    ids = linelist$case_id,
    group = linelist$group,
    date = linelist$onset,
    kappa = TRUE
  ) %>%
  map(~ .x %>% filter(kappa == 1))


# Analysis ----------------------------------------------------------------
window <- ifelse(paper == "JHI2021", 5, 7)
start <- min(linelist$onset)
end <- max(linelist$onset)
cutoff_breaks <- c(seq.Date(start, end, by = window), end) %>% unique()

if (paper == "eLife2022") {
  day_month <- format(cutoff_breaks, "%d\n%B")
  day <- format(cutoff_breaks, "%d")
  dup <- which(duplicated(format(cutoff_breaks, "%b")))
  day_month[dup] <- day[dup]
  day_month[1] <- "05 \n    March"
  day_month[length(day_month)] <- "06 \nMay    "

} else{
  day_month <- format(cutoff_breaks, "%d\n%B")
  day <- format(cutoff_breaks, "%d")
  dup <- which(duplicated(format(cutoff_breaks, "%b")))
  day_month[dup] <- day[dup]
}

plan(multisession, workers = length(cutoff_breaks))
Rematrix_df <- future_map_dfr(
  seq_len(length(cutoff_breaks) - 1),
  ~ process_window(
    .x,
    trees = trees,
    cutoff_breaks = cutoff_breaks,
    as_matrix = TRUE
  )
) %>%
  mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)

Reglobal_df <- future_map_dfr(
  seq_len(length(cutoff_breaks) - 1),
  ~ process_window(
    .x,
    trees = trees,
    cutoff_breaks = cutoff_breaks,
    by_group = TRUE,
    as_matrix = FALSE
  )
) %>%
  mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)

Re_df <- bind_rows(
  Reglobal_df %>% rename(from = from_group) %>%
    mutate(to = "overall") %>% select(from, to, window_median, window_start, window_end, Re),
  Rematrix_df %>% select(from, to, window_median, window_start, window_end, Re)
) %>%
  mutate(to = factor(to, levels = c("overall", "hcw", "patient")),
         from_label = as.factor("from"),
         to_label = as.factor("to"))


Re_summary <- Re_df %>%
  group_by(from, to, window_median, window_start, window_end) %>%
  drop_na(Re) %>%
  summarise(
    mean = mean(Re),
    lwr = quantile(Re, 0.025),
    upr = quantile(Re, 0.975),
    .groups = "drop"
  ) %>%
  mutate(from_label = as.factor("from"),
         to_label = as.factor("to"))

p_Rematrix <- Re_df %>%
  ggplot(aes(
    x = window_median,
    y = Re,
    group = interaction(from, to, from_label, to_label, window_median),
  )) +
  # facet_grid(rows = vars(from),
  #            cols = vars(as.factor(to)),
  #            switch = "y") +
  ggh4x::facet_nested(rows = vars(from_label, from),
                      cols = vars(to_label,as.factor(to)),
                      switch = "y") +
  geom_hline(
    data = tibble(
      from = c("hcw", "patient", "hcw", "patient" , "hcw", "patient"),
      to = c("overall", "overall", "hcw", "patient" , "patient", "hcw"),
      yintercept = c(1, 1, 0.5, 0.5, 0.5, 0.5)
    ) %>%
      mutate(to = factor(to, levels = c("overall", "hcw", "patient")),
             from_label = as.factor("from"),
             to_label = as.factor("to")),
    aes(yintercept = yintercept),
    linetype = "dashed"
  ) +
  geom_errorbarh(
    data = Re_summary,
    aes(xmin = window_start, xmax = window_end, y = mean),
    height = 0.1
  ) +
  geom_violin(
    aes(fill = from),
    col = NA,
    adjust = 1.8,
    position = position_dodge(width = 4),
    alpha = 0.5,
    linewidth = 0.1
  ) +
  geom_pointrange(data = Re_summary,
                  aes(
                    x = window_median,
                    y = mean,
                    ymin = lwr,
                    ymax = upr,
                  ),
                  size = 0.4) +
  scale_y_continuous(breaks = seq(0, 10, 0.5), position = "right") +
  coord_cartesian(ylim = c(0, 2.5), clip = "off") +
  labs(
    x = "",
    fill = "Infector",
    y = expression(R["e"]),
    color = expression(paste("Group level ", R["e"]))
  ) +
  theme_noso(date = TRUE)

p_Rematrix +
  theme(legend.position = "none")+
  theme(
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      size = 1,
      linetype = "solid"
    )
  )+
  scale_x_date(
    breaks = cutoff_breaks,
    labels = day_month,
    expand = c(0.01, 0.5),
    limits = c(min(cutoff_breaks)-0.5, max(cutoff_breaks)+0.5)
  )

Re_matrix_ratios(Re_summary = Re_summary) %>%
  ggplot(aes(x = window_median, y = ratio, group = from)) +
  geom_point()+
  geom_line(aes(color = from), size = 1)
