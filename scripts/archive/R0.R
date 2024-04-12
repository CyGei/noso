source(here::here("scripts", "helpers.R"))
load_libraries()
load_data()
linelist$group <-
  ifelse(grepl("^C", linelist$case_id), "patient", "hcw")
outbreaker_chains <- outbreaker_chains %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id)

trees <- get_trees(
  out = outbreaker_chains,
  id = linelist$case_id,
  group = linelist$group,
  date = linelist$onset_inferred
)

#draw R0 for each tree
R0s <- lapply(trees, function(tree) {
  draw_R0(from = tree$from_group,
          to = tree$to_group,
          f = c("hcw" = 0.8, "patient" = 0.2),
          n_cases = as.vector(table(linelist$group))) %>%
    as_tibble()
}) %>%
  bind_rows(.id = "tree")

R0_summary <- R0s %>%
  select(-tree) %>%
  summarise(across(everything(),
                   list(
                     mean = mean,
                     lwr = ~ quantile(.x, 0.025),
                     upr = ~ quantile(.x, 0.975)
                   ))) %>%
  pivot_longer(cols = everything(),
               names_to = "metric") %>%
  separate(metric, into = c("group", "metric")) %>%
  pivot_wider(names_from = metric, values_from = value)

p_R0 <- R0s %>%
  select(-tree) %>%
  pivot_longer(cols = everything(),
               names_to = "group") %>%
  ggplot(aes(x = value, fill = group, group = group)) +
  #stat_bin(binwidth = 0.1, aes(y = ..count../sum(..count..))) +
  geom_density(alpha = 0.5, col = NA) +
  geom_pointrange(
    data = R0_summary,
    aes(x = mean,
        xmin = lwr,
        xmax = upr,
        y = c(0.9, 1.1),
        color = group),
        )+
  theme_bw()+
  scale_fill_manual("Group", values = c(hcw = "orange", patient = "purple")) +
  scale_color_manual("Group", values = c(hcw = "orange", patient = "purple")) +
  labs(x = "R0",
       y = "Density")



f_seq <- seq(0.1, 0.9, 0.1)
R0_f <- lapply(f_seq, function(f) {

  #1 calculate R0 for each tree given f
  R0s <- lapply(trees, function(tree) {
    draw_R0(from = tree$from_group,
            to = tree$to_group,
            f = c("hcw" = f, "patient" = 1 - f),
            n_cases = as.vector(table(linelist$group))
            ) %>%
      as_tibble()
  }) %>%
    bind_rows(.id = "tree") %>%
    mutate(f_hcw = f,
           f_patient = 1 - f)
}) %>%
  bind_rows(.id = "f_iter") %>%
  mutate(f_label = paste0("f_hcw = ", f_hcw, ", f_patient = ", f_patient))


R0_summary <- R0_f %>%
  group_by(f_label) %>%
  summarise(across(c("hcw", "patient"),
                   list(
                     mean = mean,
                     lwr = ~ quantile(.x, 0.025),
                     upr = ~ quantile(.x, 0.975)
                   ))) %>%
  pivot_longer(cols = -f_label,
               names_to = "metric") %>%
  separate(metric, into = c("group", "metric")) %>%
  pivot_wider(names_from = metric, values_from = value)


p_R0_f <- R0_f %>%
  pivot_longer(cols = c("hcw", "patient"),
               names_to = "group") %>%
  ggplot(aes(x = value, fill = group, group = group)) +
  facet_wrap(~f_label, scales = "free") +
  geom_density(alpha = 0.5, col = NA) +
  geom_pointrange(
    data = R0_summary,
    aes(x = mean,
        xmin = lwr,
        xmax = upr,
        y = 1,
        color = group),
  )+
  theme_bw()+
  scale_fill_manual("Group", values = c(hcw = "orange", patient = "purple")) +
  scale_color_manual("Group", values = c(hcw = "orange", patient = "purple")) +
  labs(x = "R0",
       y = "Density")
p_R0_f







