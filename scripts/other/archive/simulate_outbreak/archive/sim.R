source(here::here("scripts", "helpers.R"))
load_libraries()
pacman::p_load(o2groups, dplyr, here, distcrete)

incubation_period <-
  distcrete::distcrete(
    "gamma",
    shape = epitrix::gamma_mucv2shapescale(5.95, (4.31/5.95))$shape,
    scale = epitrix::gamma_mucv2shapescale(5.95, (4.31/5.95))$scale,
    w = 0.5,
    interval = 1
  )

generation_time = c(0.1, 0.2, 0.4, 0.2, 0.1)/sum(c(0.1, 0.2, 0.4, 0.2, 0.1))
gamma_f <- readRDS(here("data", "CrI_f_delta.rds")) %>%
  bind_rows() %>%
  drop_na(value) %>%
  mutate(value = delta2gamma(value),
         fHCW = as.numeric(sub(".*hcw:(.*?)\\|.*", "\\1", f))) %>%
  filter(
    metric == "mean",
    (level == "hcw" & cutoff == "2020-04-09") |
      (level == "patient" & cutoff == "2020-03-26")
  ) %>%
  select(fHCW, value, level) %>%
  pivot_wider(names_from = level, values_from = value) %>%
  relocate(fHCW, hcw, patient)

o2groups.paramsL <- lapply(1:nrow(gamma_f), function(i){
  list(
    duration = 70,
    group_n = 2,
    size = as.vector(unlist(c(gamma_f[i, "fHCW"], 1 - gamma_f[i, "fHCW"]))),
    name = c("HCW", "patient"),
    intro_n = c(1,1),
    r0 = c(1.29, 1.61),
    generation_time = generation_time,
    incubation_period = incubation_period$r(1000),
    gamma = as.vector(unlist((gamma_f[i, c("hcw", "patient")])))
  )
})

pop <- 700
params <- o2groups.paramsL[[9]]
params$size <- params$size * pop

set.seed(123)
sims <- lapply(1:50, function(i){
  sim <- o2groups::simulate_groups(
    duration = params$duration,
    group_n = params$group_n,
    size = params$size,
    name = params$name,
    gamma = params$gamma,
    intro_n = params$intro_n,
    r0 = params$r0,
    generation_time = params$generation_time,
    incubation_period = params$incubation_period
  )
  return(sim)
})
sims_df <-bind_rows(sims, .id = "sim")
target_cases <- c(105, 43)
target_peak <- c(34, 21)
# Average number of cases
cases <- sims_df %>%
  group_by(sim, group) %>%
  summarise(case = n(), .groups = "drop") %>%
  group_by(group) %>%
  summarise(mean = mean(case), .groups = "drop") %>%
  mutate(target = target_cases)
cases
colSums(cases[, c("mean", "target")])

# Average peak date
peaks <- sims_df %>%
  group_by(sim, group, date_onset) %>%
  summarise(case = n(), .groups = "drop") %>%
  group_by(sim, group) %>%
  arrange(sim, desc(case)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  group_by(group) %>%
  summarise(mean = mean(date_onset), .groups = "drop") %>%
  mutate(target = target_peak)
peaks

# Average ttable up to the group's peak
source(here::here("scripts", "other", "simulate_outbreak", "moVSdelta.R"))
set.seed(123)
results <- moVSdelta(params)
results$out
#plot
results$out %>%
  ggplot(aes(
    x = level,
    y = est,
    ymin = lwr,
    ymax = upr,
    col = metric
  )) +
  geom_hline(aes(yintercept = 0),
             linetype = "dotted") +
  geom_point(
    aes(y = true_delta,
        shape = "truth"),
    size = 5,
    col = "#5c5e60"
  ) +
  geom_pointrange(size = 0.5, position = position_dodge(width = 0.25)) +
  scale_color_manual("metric",
                     values =
                       c("delta" = "red", "mo" = "blue")) +
  scale_shape_manual("", values = c("truth" = 15)) +
  theme_bw() +
  labs(x = "",
       y = "Estimate") +
  theme(legend.position = "bottom") +
  scale_y_continuous(limits = c(-1, 1))

sjtabs <- split(results$ttab_cut, results$ttab_cut$level) %>%
  map(~{
    .x %>%
      uncount(Freq) %>%
      select(from, to) %>%
      sjPlot::sjtab(fun = "xtab",
                    var.labels = c("infector", "infectee"),
                    show.row.prc=T,
                    show.col.prc=T,
                    show.summary=T,
                    show.exp=T,
                    show.legend=T)
  })
sjtabs[[1]] # hcw
sjtabs[[2]] # patient

results$ttab_all %>%
  as.data.frame() %>%
  uncount(Freq) %>%
  select(from, to) %>%
  sjPlot::sjtab(fun = "xtab",
                var.labels = c("infector", "infectee"),
                show.row.prc=T,
                show.col.prc=T,
                show.summary=T,
                show.exp=T,
                show.legend=T)
