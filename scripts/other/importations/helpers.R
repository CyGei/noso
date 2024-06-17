
# Function to compute imports
get_imports <- function(tree, window_start, window_end, by_group = FALSE, levels = NULL) {
  data <- tree %>%
    filter(is.na(from), to_date >= window_start, to_date <= window_end)

  if (by_group) {
    if (is.null(levels)) {
      stop("levels must be provided when by_group = TRUE")
    }
    data_grouped <- data %>%
      group_by(to_group) %>%
      summarise(n_imports = n(), .groups = "drop")

    result <- tibble(to_group = levels) %>%
      left_join(data_grouped, by = "to_group") %>%
      mutate(n_imports = replace_na(n_imports, 0),
             window_start = window_start,
             window_end = window_end)
  } else {
    result <- tibble(window_start = window_start, window_end = window_end, n_imports = nrow(data))
  }

  return(result)
}


# Function to process each window
process_window <- function(i, trees, cutoff_breaks, by_group = FALSE, levels = NULL) {
  if (i >= length(cutoff_breaks)) return(NULL)
  window_start <- cutoff_breaks[i]
  if (i == length(cutoff_breaks) - 1) {
    window_end <- cutoff_breaks[i + 1]
  } else {
  window_end <- cutoff_breaks[i + 1] - 1
  }

  imports <- trees %>%
    map_dfr(~ get_imports(.x, window_start, window_end, by_group, levels))

  imports %>%
    mutate(window_median = window_start + as.numeric(window_end - window_start) / 2)
}
