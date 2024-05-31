# The purpose of this script is to estimate the reproduction number
# using transmission chain data.
# For each type of transmission.
source(here::here("scripts", "helpers.R"))
load_libraries()
paper <- c("JHI2021", "eLife2022")[1]
input_path <- here::here("data", paper, "input/")
load_data(input_path)
out <- out %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id)
out <- filter_alpha_by_kappa(out, 1L)

trees <- get_trees(
  out = out,
  ids = linelist$case_id,
  group = linelist$group,
  date = linelist$onset
)

# get_R functions -------------------------------------------------------------

get_Rematrix <- function(linelist, tree) {
  tree <- tree[!is.na(tree$from), ]

  cases <- table(linelist$group)
  ttab <- linktree:::ttable(
    from = tree$from_group,
    to = tree$to_group,
    levels = c("hcw", "patient")
  )

  Re <- as_tibble(ttab) %>%
    left_join(as.data.frame(cases), by = c("from" = "Var1")) %>%
    rename(cases = Freq) %>%
    mutate(Re = n / cases)

  #sweep(ttab, 1, cases, "/")
  return(Re)
}




# Analysis ----------------------------------------------------------------
cutoff_dates_init <- cutoff_dates
window <- 5
cutoff_dates <- seq.Date(min(cutoff_dates_init),
                         max(cutoff_dates_init),
                         by = window) %>%
  c(., max(cutoff_dates_init)) %>% unique()

day_month <- format(cutoff_dates, "%d\n%B")
day <- format(cutoff_dates, "%d")
dup <- which(duplicated(format(cutoff_dates, "%b")))
day_month[dup] <- day[dup]

pacman::p_load(furrr)
plan(multisession, workers = length(cutoff_dates))

R_df <- furrr::future_map(cutoff_dates[-1], function(cutoff_date) {
  cutoff_date <- as.Date(cutoff_date)
  R <- lapply(trees, function(tree) {
    cut_tree <-
      tree %>%
      drop_na(from) %>%
      mutate(from_date = as.Date(from_date)) %>%
      filter(from_date >= (cutoff_date - window) &
               from_date <= cutoff_date)

    cut_linelist <- linelist %>%
      filter(onset >= (cutoff_date - window) &
               onset <= cutoff_date)

    get_Rematrix(cut_linelist, cut_tree)
  }) %>%
    bind_rows() %>%
    dplyr::mutate(window_end = cutoff_date, window_start = cutoff_date - window)
}) %>%
  bind_rows(.id = "window_id") %>%
  mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)



Rsummary <- R_df %>%
  group_by(from, to, window_id) %>%
  drop_na(Re) %>%
  summarise(
    mean = mean(Re),
    lwr = quantile(Re, 0.025),
    upr = quantile(Re, 0.975),
    window_start = first(window_start),
    window_end = first(window_end),
    window_median = first(window_median),
    cases = first(cases),
    .groups = "drop"
  ) %>%
  group_by(from, window_id) %>%
  mutate(from_Re = sum(mean)) %>%
  ungroup()

p_Re <- R_df %>%
  ggplot(aes(
    x = window_median,
    y = Re,
    group = interaction(from, to, window_id),
  )) +
  facet_grid(rows = vars(from),
             cols = vars(to),
             switch = "y") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_errorbarh(
    data = Rsummary,
    aes(
      xmin = window_start,
      xmax = window_end,
      y = from_Re,
      col = from
    ),
    height = 0.1,
    linetype = "solid",
  ) +
  geom_violin(
    adjust = 1.8,
    fill = "#636363",
    position = position_dodge(width = 4),
    alpha = 0.5,
    linewidth = 0.1
  ) +
  geom_pointrange(data = Rsummary,
                  aes(
                    x = window_median,
                    y = mean,
                    ymin = lwr,
                    ymax = upr,
                  ),
                  size = 0.4) +
  scale_y_continuous(breaks = seq(0, 2, 0.5), position = "right") +
  coord_cartesian(ylim = c(0, 2)) +
  labs(x = "",
       y = expression(R["e"]),
       color = expression(paste("Group level ", R["e"]))) +
  theme_noso(date = TRUE)+
  scale_x_date(
    breaks = cutoff_dates,
    labels = day_month,
    expand = c(0.01, 0.01)
  )+
  theme(panel.spacing.x = unit(4, "mm")) +
  theme(axis.text.x = element_text(
    hjust = c(
      0,
      rep(0.5, length(cutoff_dates) - 2),
      0.5
    )
  ))
p_Re


p_Re2021 <- p_Re
p_Re2020 <- p_Re

legend <- cowplot::get_legend(p_Re2021)
p_Re<- cowplot::plot_grid(p_Re2020 +
                     theme(legend.position = "none") +
                     theme(axis.text.y = element_blank(),
                           axis.title.y = element_blank(),
                           axis.ticks.y = element_blank()),
                   p_Re2021 +
                     theme(legend.position = "none")+
                     theme(strip.text.y = element_blank()),
                   ncol = 2,
                   labels = "AUTO")

p_Re_lgd <- cowplot::plot_grid(p_Re, legend, ncol = 1, rel_heights = c(1, 0.1))
p_Re_lgd
