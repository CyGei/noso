# The purpose of this script is to estimate the basic reproduction number R0
# using transmission chain data.

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

get_R0 <- function(linelist, tree, by_group = FALSE){
  tree <- tree %>% drop_na(from)
  if (by_group) {
    transmissions <- tree %>%
      group_by(from_group) %>%
      summarise(transmissions = n(), .groups = "drop")
    cases <- linelist %>%
      group_by(group) %>%
      summarise(cases = n(), .groups = "drop")
    R0 <- transmissions %>%
      left_join(cases, by = c("from_group" = "group")) %>%
      mutate(R0 = transmissions / cases)
  } else {
    transmissions <- nrow(tree)
    cases <- nrow(linelist)
    R0 <- data.frame(R0 = transmissions / cases)
  }

  return(R0)
}

pacman::p_load(furrr)
plan(multisession, workers = length(cutoff_dates))
R0_df <- furrr::future_map(cutoff_dates, function(cutoff_date) {
  R0 <- lapply(trees, function(tree) {
    # truncate tree to only include infectors reporting onset before cutoff_date
    # this will retain transmissions even past the cutoff_date
    # What kind of censoring is this? Right?
    cut_tree <- filter(tree, from_date <= cutoff_date) %>% drop_na(from)
    cut_linelist <- filter(linelist, onset_inferred <= cutoff_date)
    get_R0(cut_linelist, cut_tree, by_group = TRUE)
  }) %>%
    bind_rows() %>%
    dplyr::mutate(cutoff_date = cutoff_date)
}) %>%
  bind_rows()


p_R0 <- R0_df %>%
  group_by(from_group, cutoff_date) %>%
  summarise(mean = mean(R0),
            lwr = quantile(R0, 0.025),
            upr = quantile(R0, 0.975)) %>%
  ggplot(aes(
    x = cutoff_date,
    y = mean,
    ymin = lwr,
    ymax = upr,
    color = from_group
  )) +
  geom_pointrange(position = position_dodge(width = 1), size = 0.5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual("Group",
                     values = c("hcw" = "orange", "patient" = "purple"),
                     guide = "none")+
  scale_x_date(breaks = cutoff_dates, date_labels = "%d\n%b") +
  theme_bw()+
  labs(
    x = "",
    y = "R0"
  )

cowplot::plot_grid(p_R0,
                   epicurve(),
                   ncol = 1,
                   labels = "AUTO"
)
