source(here::here("scripts", "helpers.R"))
load_libraries()

# Mo by transmission type -------------------------------------------------
#Here we calculate the Mo metric for each type of transmission

mo_df <- lapply(c("eLife2022", "JHI2021"), function(paper) {
  input_path <- here::here(paste0("data/", paper, "/input/"))
  output_path <- here::here(paste0("data/", paper, "/output/"))
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


  plan(multisession, workers = availableCores() - 1)
  mo <- future_map(trees, ~ {
    draw_mo(
      from = .x$from_group,
      to = .x$to_group,
      levels = c("hcw", "patient"),
      args = list(
        from_id = .x$from,
        to_id = .x$to,
        diag = FALSE,
        scale = FALSE
      )
    ) %>%
      as.data.frame() %>%
      rownames_to_column("from") %>%
      pivot_longer(-from, names_to = "to", values_to = "mo")
  }, options = furrr_options(seed = TRUE)) %>%
    bind_rows(.id = "tree") %>%
    mutate(paper = paper)

  mo
}) %>%
  bind_rows()

pacman::p_load(ggplot2, gghalves)

mo_df %>%
  mutate(mo = linktree::gamma2delta(mo)) %>%
  ggplot(aes(
    x = "",
    y = mo,
    group = paper,
    fill = paper,
    split = paper
  )) +
  facet_grid(rows = vars(from),
             cols = vars(to),
             switch = "y") +
  geom_hline(yintercept = 0, linetype = 2) +
  gghalves::geom_half_violin(
    color = "black",
    linewidth = 0.1,
    alpha = 0.7,
    position = position_dodge(width = 0.1)
  ) +
  # gghalves::geom_half_dotplot(
  #   binaxis = "y",
  #   stackdir = "down",
  #   stackratio = 0.6,
  #   dotsize = 0.35,
  #   binwidth = 0.04
  # )+
  stat_summary(
    fun = mean,
    fun.min = function(x)
      quantile(x, 0.025),
    fun.max = function(x)
      quantile(x, 0.975),
    geom = "pointrange",
    position = position_dodge(width = 0.1)
  ) +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    position = "right"
  ) +
  scale_fill_manual(values = c("#ff0085", "#00a087")) +
  theme_noso(fill = FALSE, col = FALSE, date = FALSE) +
  labs(x = NULL) +
  theme(axis.ticks = element_blank())



# Global Mo by case type --------------------------------------------------
#Here we calculate the Mo excess at the case type level
mo_df <- lapply(c("eLife2022", "JHI2021"), function(paper) {
  input_path <- here::here(paste0("data/", paper, "/input/"))
  output_path <- here::here(paste0("data/", paper, "/output/"))
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

  f_cases <- prop.table(table(linelist$group)) %>%
    as.data.frame() %>%
    rename(level = Var1,
           f_case = Freq)

  mo <- lapply(trees, function(tree) {
    prop.table(table(tree$from_group)) %>%
      as.data.frame() %>%
      rename(level = Var1)
  }) %>%
    bind_rows() %>%
    # group_by(level) %>%
    # summarise(
    #   mean = mean(Freq),
    #   lower = quantile(Freq, 0.025),
    #   upper = quantile(Freq, 0.975),
    #   .groups = "drop"
    # ) %>%
    left_join(y = f_cases, by = "level") %>%
    mutate(paper = paper)

  return(mo)
}) %>%
  bind_rows()


mo_df %>%
  ggplot(aes(
    x = "",
    y = Freq,
    group = level,
    fill = level,
    split = level
  )) +
  facet_wrap( ~ paper) +
  geom_hline(aes(yintercept = f_case, color = level,
                 linetype = "Proportion of cases")) +
  gghalves::geom_half_violin(
    color = "black",
    linewidth = 0.1,
    alpha = 0.7,
    position = position_dodge(width = 0.1)
  ) +
  stat_summary(
    fun = mean,
    fun.min = function(x)
      quantile(x, 0.025),
    fun.max = function(x)
      quantile(x, 0.975),
    geom = "pointrange",
    position = position_dodge(width = 0.1)
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_linetype_manual(values = c("Proportion of cases" = 2)) +
  theme_noso(fill = TRUE, col = TRUE, date = FALSE) +
  labs(x = NULL,
       y = "Proportion of infectors",
       fill = "",
       color = "",
       linetype = "") +
  theme(axis.ticks = element_blank())



# Does f_hcw change over time ---------------------------------------------

#Investigate the relative prevalence of each case type over time
papers <- c("eLife2022", "JHI2021")

plots <- lapply(papers, function(paper) {
  input_path <- here::here(paste0("data/", paper, "/input/"))
  output_path <- here::here(paste0("data/", paper, "/output/"))

  load_data(input_path)

  # use cutoff_dates to count the % of cases in each group at each timepoint
  p_prevalence <- lapply(cutoff_dates, function(cutoff) {
    linelist %>%
      filter(onset <= cutoff) %>%
      group_by(group) %>%
      summarise(n = n()) %>%
      mutate(pct = n / sum(n),
             cutoff = cutoff)
  }) %>%
    bind_rows() %>%
    ggplot(aes(x = cutoff, y = pct, fill = group)) +
    geom_bar(stat = "identity") +
    labs(x = "Cutoff", y = "Proportion of cases") +
    theme_noso(fill = TRUE, col = FALSE, date = TRUE)

  calculate_ci_wilson <- function(df, cutoff) {
    z <- 1.96
    df %>%
      filter(onset <= cutoff) %>%
      group_by(group) %>%
      summarise(n = n()) %>%
      mutate(
        pct = n / sum(n),
        lower = (pct + z^2 / (2 * n) - z * sqrt((pct * (1 - pct) + z^2 / (4 * n)) / n)) / (1 + z^2 / n),
        upper = (pct + z^2 / (2 * n) + z * sqrt((pct * (1 - pct) + z^2 / (4 * n)) / n)) / (1 + z^2 / n),
        cutoff = cutoff
      )
  }

  p_prevalence <- lapply(cutoff_dates, function(cutoff) {
    calculate_ci_wilson(linelist, cutoff)
  }) %>%
    bind_rows() %>%
    ggplot(aes(x = cutoff, y = pct, fill = group)) +
    geom_line(aes(group = group, col = group)) +
    geom_point(aes(col = group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = group), alpha = 0.2) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    theme_noso(fill = TRUE, col = TRUE, date = TRUE) +
    theme(legend.position = "none")

  # if(paper == "JHI2021") {
  #   final_plot <- cowplot::plot_grid(
  #     (p_prevalence + labs(x = "", y = "")),
  #     (epicurve() + labs(x = "", y = "")),
  #     align = "v",
  #     ncol = 1,
  #     rel_heights = c(1, 1)
  #   )
  # }
  # final_plot <- cowplot::plot_grid(
  #   p_prevalence + labs(x = "", y = "Case prevalence (%)"),
  #   epicurve(),
  #   align = "v",
  #   ncol = 1,
  #   rel_heights = c(1, 1)
  # )

  return(list(
    p_prevalence,
    epicurve()
  )
  )
})

p_col1 <- cowplot::plot_grid(
  plots[[1]][[1]] + theme(legend.position = "none")+ labs(x = "", y = "Case prevalence (%)"),
  plots[[1]][[2]] + theme(legend.position = "none"),
  ncol = 1,
  labels = c("A1", "B1")
)
p_col2 <- cowplot::plot_grid(
  plots[[2]][[1]] + theme(legend.position = "none") + labs(x = "", y = ""),
  plots[[2]][[2]] + theme(legend.position = "none") + labs(y = ""),
  ncol = 1,
  labels = c("A2", "B2")
)

p_cols <- cowplot::plot_grid(
  p_col1,
  p_col2
)

p_final <- cowplot::plot_grid(
  p_cols,
  cowplot::get_legend(epicurve()),
  nrow = 2,
  rel_heights = c(1, 0.1)
)
p_final
