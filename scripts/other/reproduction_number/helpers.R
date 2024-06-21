get_Re <- function(linelist, tree, by_group = FALSE) {
  #tree <- tree[!is.na(tree$from), ]

  individualR <- count(tree, from)

  df <- merge(
    x = linelist[, c("case_id", "group")],
    y = individualR,
    by.x = "case_id",
    by.y = "from",
    all.x = TRUE
  )
  df$n <- ifelse(is.na(df$n), 0, df$n)

  if (by_group) {
    Re <- df %>%
      rename(from_group = group) %>%
      group_by(from_group) %>%
      summarise(Re = mean(n),
                cases = n(),
                .groups = "drop")
  } else {
    Re <- data.frame(from_group = "Global",
                     Re = mean(df$n),
                     cases = nrow(df))
  }

  return(Re)
}


get_Rematrix <- function(linelist, tree) {
 # tree <- tree[!is.na(tree$from), ]

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





# Function to process each window
process_window <- function(i, trees, cutoff_breaks, by_group = FALSE, as_matrix = FALSE) {
  if (i >= length(cutoff_breaks))
    return(NULL)
  window_start <- cutoff_breaks[i]
  if (i == length(cutoff_breaks) - 1) {
    window_end <- cutoff_breaks[i + 1]
  } else {
    window_end <- cutoff_breaks[i + 1] - 1
  }

  trees <- trees %>%
    map( ~ {
      .x %>%
        drop_na(from) %>%
        filter(from_date >= window_start &
                 from_date <= window_end)
    })

  window_linelist <- linelist %>%
    filter(onset >= window_start &
             onset <= window_end)

  if (as_matrix) {
    result <- trees %>%
      map_dfr( ~ get_Rematrix(
        tree = .x,
        linelist = window_linelist
      ))
  } else {
    result <- trees %>%
      map_dfr( ~ get_Re(
        tree = .x,
        linelist = window_linelist,
        by_group = by_group
      ))
  }

  result %>%
    mutate(window_start = window_start, window_end = window_end)
}
