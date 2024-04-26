#https://www.journalofhospitalinfection.com/article/S0195-6701(21)00308-X/fulltext
source(here::here("scripts", "helpers.R"))
load_libraries()
load(here("data","paper1", "data_for_Thibaut_20210225.RData"))
outbreaker_chains <- out
linelist <- ids
rm(list = c("out", "ids"))

linelist$group <- ifelse(grepl("^C", linelist$case_id), "patient", "hcw")

outbreaker_chains <- outbreaker_chains %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id)
outbreaker_chains <- filter_alpha_by_kappa(outbreaker_chains, 1)

#onset unavailable, use average date of infection from outbreaker
linelist <- outbreaker_chains %>%
  select(starts_with("t_inf")) %>%
  pivot_longer(cols = everything(), names_to = "case_id", values_to = "t_inf") %>%
  mutate(case_id = str_remove(case_id, "t_inf_")) %>%
  group_by(case_id) %>%
  summarise(t_inf =round(mean(t_inf)), .groups = "drop") %>%
  left_join(x = linelist, y = ., by = "case_id") %>%
  mutate(
    onset_inferred = t_inf + min(t_inf),
    onset_inferred = as.Date(t_inf, origin = "2020-01-01"))

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
epicurve()

metrics <- list(
  "delta" = draw_delta,
  "Mo" = draw_mo
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
                      f = c("hcw" = 0.23, "patient" = 0.77),
                      from_id = "from",
                      to_id = "to",
                      diag = TRUE
                    ),
                    cutoff_dates = cutoff_dates,
                    n_samples = ifelse(metric_name == "Mo", 1, 1000)
                  )
                })

#Rescale Mo
array[[2]] <- apply(array[[2]], 1:length(dim(array[[2]])), gamma2delta)

CrI <- lapply(array, function(x) draw_CrI(x, c(2,3)))

df <- lapply(CrI, function(x) {
  x %>%
    reshape2::melt() %>%
    pivot_wider(names_from = Var1, values_from = value)
}) %>%
  bind_rows(.id = "metric") %>% #
  mutate(metric = factor(metric, levels = 1:2, labels = names(metrics)),
         target = 0)
df




p_metrics <- ggplot(df,
                    aes(
                      x = as.Date(cutoff),
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
  scale_y_continuous(limit = c(-1, 1),
                     breaks = seq(-1, 1, 0.5)) +
  scale_x_date(breaks = cutoff_dates, date_labels = "%d\n%b") +
  labs(x = "",
       y = "metric value") +
  theme_bw() +
  theme(legend.position = "none")



cowplot::plot_grid(p_metrics,
                   epicurve(),
                   ncol = 1,
                   rel_heights = c(1, 0.5),
                   labels = "AUTO"
)

#number of cases
table(linelist$group)

#averge ttable
ttab <- lapply(trees, function(x) {
  linktree:::ttable(
    from = x$from_group ,
    to = x$to_group,
    levels = c("hcw", "patient")
  ) %>%
    as.data.frame()
}) %>%
  bind_rows(.id = "sim_n") %>%
  group_by(from, to) %>%
  summarise(Freq = mean(Freq), .groups = "drop")

ttab %>%
  mutate(Freq = round(Freq)) %>%
  uncount(Freq) %>%
  select(from, to) %>%
  sjPlot::sjtab(fun = "xtab",
                var.labels = c("infector", "infectee"),
                show.row.prc=T,
                show.col.prc=T,
                show.summary=T,
                show.exp=T,
                show.legend=T)


