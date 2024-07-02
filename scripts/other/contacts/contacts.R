load_helpers()
load_paper()
pacman::p_load(linktree)

# Realtime contacts from linelist ---------------------------------------------------------
contacts_rt <- loadcontacts(paper)


delta_rt <- map(contacts_rt, ~ {
  ttab <- linktree:::ttable(
    from = .x$from_group,
    to = .x$to_group,
    levels = c("hcw", "patient")
  )
  pi_df <- contact_pi(ttab)
  delta_df <-  tibble(hcw = contact_delta(pi_df$hcw, 0.23),
                      patient = contact_delta(pi_df$patient, 0.77)) %>%
    pivot_longer(cols = everything(),
                 names_to = "level",
                 values_to = "delta")
  return(delta_df)
}) %>%
  bind_rows(.id = "cutoff") %>%
  mutate(cutoff_date = cutoff_dates[as.integer(cutoff)])

summary_rt <- delta_rt %>%
  group_by(cutoff_date, level) %>%
  summarise(
    mean = mean(delta),
    lwr = quantile(delta, 0.025),
    upr = quantile(delta, 0.975),
    .groups = "drop"
  ) %>%
  mutate(significant = ifelse(lwr > 0 |
                                upr < 0, "significant", "not significant"))

ggplot(delta_rt,
       aes(
         x = cutoff_date,
         y = delta,
         fill = level,
         group = interaction(cutoff_date, level)
       )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(
    col = NA,
    adjust = 0.5,
    position = position_dodge(width = 0.2),
    scale = "width"
  ) +
  geom_point(
    data = summary_rt,
    aes(y = mean, shape = significant),
    position = position_dodge(width = 0.2),
    size = 1.5
  ) +
  geom_errorbar(
    data = summary_rt,
    aes(y = mean, ymin = lwr, ymax = upr),
    position = position_dodge(width = 0.2),
    width = 0.1
  ) +
  theme_noso(date = T) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_shape_manual(values = c("significant" = 16, "not significant" = 1)) +
  labs(x = "",
       y = "Assortativity",
       fill = "Group",
       shape = "")



#########################################################################
# Individual pi ---------------------------------------------------------
# Function to count contacts for a single ID
contact_counts = map(contacts_rt[-1], ~ count_contacts(.x)) %>%
  bind_rows(.id = "cutoff") %>%
  mutate(
    cutoff_date = cutoff_dates[as.integer(cutoff)],
    level = ifelse(startsWith(id, "C"), "patient", "hcw"),
    delta = contact_delta(freq, ifelse(level == "hcw", 0.23, 0.77))
  ) %>%
  filter(level == contact_level)


stat_contacts <- contact_counts %>%
  group_by(level, cutoff_date) %>%
  summarise(
    mean = mean(delta),
    lwr = quantile(delta, 0.025),
    upr = quantile(delta, 0.975),
    .groups = "drop"
  ) %>%
  mutate(significant = ifelse(lwr > 0 |
                                upr < 0, "significant", "not significant"))

ggplot(contact_counts,
       aes(
         x = cutoff_date,
         y = delta,
         fill = level,
         group = interaction(cutoff_date, level)
       )) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin(
    col = NA,
    adjust = 0.5,
    position = position_dodge(width = 0.2),
    scale = "width"
  ) +
  geom_point(
    data = stat_contacts,
    aes(y = mean, shape = significant),
    position = position_dodge(width = 0.2),
    size = 1.5
  ) +
  geom_errorbar(
    data = stat_contacts,
    aes(y = mean, ymin = lwr, ymax = upr),
    position = position_dodge(width = 0.2),
    width = 0.1
  ) +
  theme_noso(date = T) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  scale_shape_manual(values = c("significant" = 16, "not significant" = 1)) +
  labs(x = "",
       y = "Assortativity",
       fill = "Group",
       shape = "")





# STATIC ANALYSIS ---------------------------------------------------------

contacts <- readRDS(here::here("data", paper, "input", "contacts.rds")) %>%
  mutate(
    from_group = ifelse(startsWith(from_id, "C"), "patient", "hcw"),
    to_group   = ifelse(startsWith(to_id, "C"), "patient", "hcw")
  ) %>%
  filter(from_id != to_id) %>%
  mutate(pair = paste(pmin(from_id, to_id), pmax(from_id, to_id))) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)


id_tab <- tibble(id = unique(c(contacts$from_id, contacts$to_id)),
                 level = ifelse(startsWith(id, "C"), "patient", "hcw"))
table(id_tab$level) %>% prop.table()

ttab <- linktree:::ttable(
  from = contacts$from_group,
  to = contacts$to_group,
  levels = c("hcw", "patient")
)
pi_df <- contact_pi(ttab)
delta_df <-  tibble(hcw = contact_delta(pi_df$hcw, 0.23),
                    patient = contact_delta(pi_df$patient, 0.77)) %>%
  pivot_longer(cols = everything(),
               names_to = "level",
               values_to = "delta")

summary_df <- delta_df %>%
  group_by(level) %>%
  summarise(
    mean = mean(delta),
    lwr = quantile(delta, 0.025),
    upr = quantile(delta, 0.975),
    .groups = "drop"
  ) %>%
  mutate(significant = ifelse(lwr > 0 |
                                upr < 0, "significant", "not significant"))



ggplot(delta_df, aes(x = delta, fill = level)) +
  geom_density(col = NA) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = summary_df, aes(x = mean, y = 0), size = 3) +
  geom_errorbarh(data = summary_df,
                 aes(
                   x = mean,
                   xmin = lwr,
                   xmax = upr,
                   y = 0
                 ),
                 height = 0.05) +
  theme_noso(date = F) +
  labs(x = "Assortativity", y = "Density", shape = "") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25))



# Static analysis with individual pi ------------------------------------

delta_df <- count_contacts(contacts) %>%
  mutate(level = ifelse(startsWith(id, "C"), "patient", "hcw"),
       delta = contact_delta(freq, ifelse(level == "hcw", 0.23, 0.77))) %>%
  filter(level == contact_level)


summary_df <- delta_df %>%
  group_by(level) %>%
  summarise(
    mean = mean(delta),
    lwr = quantile(delta, 0.025),
    upr = quantile(delta, 0.975),
    .groups = "drop"
  ) %>%
  mutate(significant = ifelse(lwr > 0 |
                                upr < 0, "significant", "not significant"))



ggplot(delta_df, aes(x = delta, fill = level)) +
  geom_histogram(binwidth = 0.01, col = NA) +
  #geom_density(col = NA, alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = summary_df, aes(x = mean, y = 0), size = 3) +
  geom_errorbarh(data = summary_df,
                 aes(
                   x = mean,
                   xmin = lwr,
                   xmax = upr,
                   y = 0
                 ),
                 height = 0.05) +
  theme_noso(date = F) +
  labs(x = "Assortativity", y = "Count", shape = "") +
  scale_x_continuous(breaks = seq(-1, 1, 0.25))+
  coord_cartesian(xlim = c(-1, 1))
