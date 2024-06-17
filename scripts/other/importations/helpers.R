intro_df <- furrr::future_map(seq_along(cutoff_breaks), function(i) {
  window_start <- cutoff_breaks[i]
  window_end <- cutoff_breaks[i + 1] - 1

  intros <- lapply(trees, function(tree) {
    intro_window <- tree %>%
      filter(is.na(from)) %>%
      arrange(to_date) %>%
      filter(to_date >= window_start & to_date <= window_end) %>%
      summarise(n_intros = sum(is.na(from)), .groups = "drop") %>%
      mutate(window_start = window_start, window_end = window_end)

    # Create a data frame with all possible groups
    all_groups <- data.frame(window_start = window_start, window_end = window_end)

    # Merge the all_groups data frame with the intro_window data frame
    intro_window <- merge(all_groups, intro_window, all.x = TRUE)

    # Replace NA values with 0
    intro_window$n_intros[is.na(intro_window$n_intros)] <- 0

    return(intro_window)
  }) %>%
    bind_rows()
  return(intros)
}) %>%
  bind_rows(.id = "window_id") %>%
  mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)

bind_rows(intro_df$data) %>%
  ggplot(aes(x = window_median, y = n_intros)) +
  geom_violin(aes(group = window_id), #adjust = 1.8,
              col = "black", alpha = 0.7) +
  scale_y_continuous(breaks = seq(0, 14, 2)) +
  theme_noso(day_break = 7, date = TRUE) +
  coord_cartesian(ylim = c(0, 14)) +
  labs(x = "", y = "Importations")








tree <- trees[[1]]
window_start <- cutoff_breaks[1]
window_end <- cutoff_breaks[2] - 1
levels <- unique(linelist$group)

get_imports <- function(tree,
                        window_start,
                        window_end,
                        by_group = FALSE,
                        levels = NULL) {
  data <- tree %>%
    filter(is.na(from), to_date >= window_start, to_date <= window_end)

  tab <- tibble(
    window_start = window_start,
    window_end = window_end,
    window_median = window_start + as.numeric(window_end - window_start) / 2,
    n_imports = nrow(data)
  )

  if (by_group) {
    #check that levels is not NULL
    if (is.null(levels)) {
      stop("levels must be provided when by_group = TRUE")
    }
    y <- tree %>%
      filter(is.na(from), to_date >= window_start, to_date <= window_end) %>%
      group_by(to_group) %>%
      summarise(n_imports = n(), .groups = "drop")

    tab <- tibble(
      to_group = levels,
      window_start = window_start,
      window_end = window_end,
      window_median = window_start + as.numeric(window_end - window_start) / 2
    ) %>%
      left_join(y, by = "to_group") %>%
      mutate(n_imports = ifelse(is.na(n_imports), 0, n_imports))
  }

  return(tab)
}
get_imports(tree, window_start, window_end, by_group = TRUE)










get_impots <- function(trees, cutoff_breaks) {
  # Create window pairs
  windows <- tibble(
    window_start = cutoff_breaks[-length(cutoff_breaks)],
    window_end = cutoff_breaks[-1] - 1,
    window_id = seq_along(window_start)
  )

  # Function to process a single tree
  process_tree <- function(tree, window) {
    tree %>%
      filter(is.na(from),
             to_date >= window$window_start,
             to_date <= window$window_end) %>%
      summarise(n_imports = n(), .groups = "drop") %>%
      right_join(window, by = character()) %>%
      mutate(n_imports = coalesce(n_imports, 0L))
  }

  # Process all trees for all windows
  intro_df <- windows %>%
    mutate(data = future_map(row_number(), ~ map(trees, process_tree, window = windows[.x, ]))) %>%
    unnest(data) %>%
    mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)

  return(intro_df)
}

# Usage
intro_df <- get_impots(trees, cutoff_breaks)
