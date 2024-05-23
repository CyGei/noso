# Re.R gave the estimates for the basic reproduciton number as follow:
# Global R0 = 1.36
# Patient R0 = 1.61
# HCW R0 = 1.29
source(here::here("scripts", "helpers.R"))
load_libraries()
load_data()
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
cutoff_date = cutoff_dates[2]
cut_linelist <- linelist %>%
  filter(onset_inferred <= cutoff_date)
gammaf <- data.frame(
  fHCW = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
  hcw = c(33.77, 15.13, 8.90, 5.79, 3.91, 2.65, 1.74, 1.05, 0.49),
  patient = c(0.21, 0.42, 0.67, 0.97, 1.37, 1.94, 2.85, 4.59, 9.63)
)
cases <- table(cut_linelist$group)
trueR0 <- data.frame(
  level = c("hcw", "patient"),
  R0 = c(1.29, 1.61)
)


Rgamma <- lapply(trees, function(tree) {
  cut_tree <-
    tree %>%
    drop_na(from) %>%
    mutate(from_date = as.Date(from_date)) %>%
    filter(from_date <= cutoff_date)

  ttab <- linktree:::ttable(from = cut_tree[["from_group"]],
                            to = cut_tree[["to_group"]],
                            levels = c("hcw", "patient"))

  f_effect <- lapply(1:nrow(gammaf), function(i) {
    R0_formula(
      tau = diag(ttab),
      gamma = as.numeric(gammaf[i, c("hcw", "patient")]),
      n_cases = as.numeric(cases),
      f = as.numeric(c(gammaf[i, "fHCW"], 1 - gammaf[i, "fHCW"]))
    )
  }) %>% bind_rows(.id = "f")

}) %>% bind_rows(.id = "tree_id") %>%
  mutate(fHCW = factor(f, labels = gammaf$fHCW))


pal <- scales::seq_gradient_pal(low = "#f1dddf", high = "#9d0759")

Rgamma %>%
  pivot_longer(
    cols = c("hcw", "patient"),
    names_to = "level",
    values_to = "R0"
  ) %>%
  ggplot(aes(x = level, y = R0)) +
  gghalves::geom_half_dotplot(aes(group = interaction(fHCW, level), fill = fHCW),
                              binaxis = "y",
                              method="histodot",
                              binwidth = 0.1,
                              stackdir="up",
                              dotsize = 0.01,
                              stackratio = 0.5) +
  gghalves::geom_half_violin(aes(group = interaction(fHCW, level), fill = fHCW),
                             adjust = 1,
                             col = NA,
                             side = "l") +
  geom_errorbar(
    data = trueR0,
    aes(
      x = as.factor(level),
      y = R0 ,
      ymin = R0,
      ymax = R0
    ),
    linewidth = 0.7,
    col = "#66696f"
  ) +
  stat_summary(fun.y = mean,
               fun.ymin = function(x) quantile(x, 0.025),
               fun.ymax = function(x) quantile(x, 0.975),
               aes(group = interaction(fHCW, level)),
               position = position_dodge(width = 0.75),
               geom = "pointrange") +
  stat_summary(fun.y = median,
               aes(group = interaction(fHCW, level)),
               position = position_dodge(width = 0.75),
               size = 4,
               shape = 18,
               col = "#76cc8a",
               geom = "point") +
  scale_fill_manual("% of HCW amongst suceptibles", values = pal(seq(0, 1, 0.1))) +
  theme_noso(col = F, fill = F, date = F) +
  labs(x = "", color = "True R0") +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(0.8, "cm"),
      title.position = "top",
      label.position = "bottom"
    ),
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      keywidth = unit(1.5, "cm"),
      title.position = "top",
      label.position = "bottom"
    )
  )


Rgamma %>%
  pivot_longer(
    cols = c("hcw", "patient"),
    names_to = "level",
    values_to = "R0"
  ) %>%
  group_by(fHCW, level) %>%
  summarise(
    mean = mean(R0),
    lwr = quantile(R0, 0.025),
    upr = quantile(R0, 0.975),
    .groups = "drop"
  ) %>%
  ggplot() +
  geom_point(
    data = trueR0,
    aes(x = level, y = R0),
    shape = 15,
    size = 4,
    col = "green"
  )+
  geom_pointrange(
    aes(
    x = level,
    y = mean,
    ymin = lwr,
    ymax = upr,
    group = interaction(fHCW, level),
    col = as.numeric(as.character(fHCW))
  ),
  position = position_dodge(width = 0.5)) +
  scale_colour_gradientn(name = '% of HCWs amongst suceptibles', colours = c('blue', 'red'),
                         breaks = seq(0, 1, 0.1),
                         labels = scales::percent(seq(0, 1, 0.1))) +
  #scale_color_viridis_d("% of HCW" ,option = "inferno")+
  theme_noso(col = F, fill = F, date = F)+
  guides(colour = guide_colorbar(barwidth = 20, barheight = 1,
                                 title.position = "top"))
