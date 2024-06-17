#######################################################
# Load libraries, data and draw_* functions
#######################################################
load_libraries <- function() {
  # Load libraries
  if (!require("pacman")) {
    install.packages("pacman")
  }
  pacman::p_load(here, outbreaker2, tidyverse, purrr, furrr, reshape2, cowplot)
  pacman::p_load_gh("CyGei/linktree")
}

load_data <- function(paper) {
  linelist <- readRDS(here("data", paper, "input/linelist.rds"))
  cutoff_dates <- sort(unique(linelist$onset))

  #outbreaker
  outbreaker_path <- here("data", paper, "output/outbreaker/")
  outbreaker_files <- list.files(outbreaker_path, full.names = FALSE)
  cutoff_days <- sort(as.integer(str_extract(outbreaker_files, "\\d+")))
  cutoff_map <- setNames(cutoff_dates, cutoff_days)

  # Load and process outbreaker files
  out <- lapply(cutoff_days, function(cutoff_day) {
    outbreaker_data <- readRDS(here(
      outbreaker_path,
      paste0("outbreaker_", cutoff_day, ".rds")
    ))
    outbreaker_data$cutoff <- cutoff_day
    outbreaker_data$cutoff_date <- cutoff_map[[as.character(cutoff_day)]]
    outbreaker_data
  })

  if (paper == "eLife2022") {
    #remove 1st cutoff off of everything since there's only one case
    out <- out[-1]
    cutoff_dates <- cutoff_dates[-1]
    cutoff_days <- cutoff_days[-1]
    cutoff_map <- cutoff_map[-1]
  }

  list2env(
    list(
      linelist = linelist,
      out = out,
      cutoff_dates = cutoff_dates,
      cutoff_days = cutoff_days,
      cutoff_map = cutoff_map
    ),
    envir = .GlobalEnv
  )
}

load_draw_functions <- function() {
  # Load draw_* functions
  source(here::here("scripts", "draw_functions.R"))
}
load_draw_functions()

#######################################################
# outbreaker2 output processing
#######################################################

#' Replace integers from the outbreaker2 output with unique identifiers.
#'
#' @param out A data frame of class "outbreaker_chains".
#' @param ids A vector of IDs from the original linelist, see outbreaker2::outbreaker_data()[["ids"]].
#'
#' @return The input data frame with the integers replaced by the corresponding IDs.
#'
#' @export
#' @examples
#' \dontrun{
#' identify(out, ids)
#' }
#'
identify <- function(out, ids) {
  # Check inputs
  stopifnot(is.data.frame(out))
  stopifnot(inherits(out, "outbreaker_chains"))
  stopifnot(is.character(ids))

  alpha_cols <- grep("^alpha_", names(out), value = TRUE)
  stopifnot(length(ids) == length(alpha_cols))

  out[alpha_cols] <- lapply(out[alpha_cols], function(x) {
    ids[x]
  })
  cols_idx <- grep("(\\d+)$", names(out))
  cols <- names(out)[cols_idx]
  digits <- as.numeric(gsub("[^0-9]", "", cols))
  cols_replaced <- paste0(gsub("\\d+$", "", cols), ids[digits])
  names(out)[cols_idx] <- cols_replaced

  return(out)
}

#' Filter alpha values based on kappa threshold
#'
#' This function filters the alpha values in an outbreaker2 data frame based on a given kappa threshold.
#'
#' @param out The outbreaker2 data frame.
#' @param kappa_threshold The threshold value for kappa.
#'
#' @return The modified outbreaker2 data frame with filtered alpha values.
#'
#'
#' @export
filter_alpha_by_kappa <- function(out, kappa_threshold) {
  if (!is.data.frame(out)) {
    stop("The 'out' argument must be an outbreaker2 data frame.")
  }
  if (!is.integer(kappa_threshold) ||
      length(kappa_threshold) != 1) {
    stop("The 'kappa_threshold' argument must be a single integer value.")
  }

  alpha_cols <- grep("^alpha_", names(out), value = TRUE)
  kappa_cols <- grep("^kappa_", names(out), value = TRUE)

  if (length(alpha_cols) != length(kappa_cols)) {
    stop("The number of 'alpha_' and 'kappa_' columns do not match.")
  }

  for (i in seq_along(alpha_cols)) {
    alpha_col <- alpha_cols[i]
    kappa_col <- kappa_cols[i]
    out[[alpha_col]][out[[kappa_col]] > kappa_threshold] <- NA
  }

  return(out)
}


#' Return a list of posterior transmission trees from an outbreaker2 object.
#'
#' This function takes an outbreaker2 object and returns a list of data frames with the 'from' and 'to' columns, and additional columns if provided.
#' The additional arguments should be vectors of values, and the name of the argument will be used as the name of the additional column. The additional columns will be named 'from_' and 'to_' followed by the name of the argument.
#'
#' @param out  A data frame of class "outbreaker_chains".
#' @param ids  A vector of IDs from the original linelist.
#' @param kappa A logical value indicating whether to include the kappa values in the output.
#' @param ... Additional columns from the original linelist to include. Each argument should be an atomic vector of values, and the name of the argument will be used as the name of the additional column.
#'
#' @return A list of data frames. Each data frame has 'from' and 'to' columns,
#'         and additional columns based on the additional arguments.
#' @export
#'
#' @examples
#' \dontrun{
#' get_trees(out, ids = ids, group = group, dates = dates)
#' }
#'
get_trees <- function(out, ids, kappa = FALSE, t_inf = FALSE, ...) {
  # Check inputs
  stopifnot(is.data.frame(out))
  args <- list(...)
  stopifnot(all(sapply(args, is.atomic)))

  # Create a mapping of ids to their original indices
  id_map <- setNames(seq_along(ids), ids)

  # Retrieve all columns starting with alpha_
  alpha_cols <- grep("^alpha_", names(out), value = TRUE)
  to <- gsub("alpha_", "", alpha_cols)

  kappa_cols <- grep("^kappa_", names(out), value = TRUE)
  t_inf_cols <- grep("^t_inf_", names(out), value = TRUE)

  # Create a list of data frames, each with 'from' and 'to' columns
  tree_list <- lapply(seq_len(nrow(out)), function(i) {
    from <- unlist(out[i, alpha_cols], use.names = FALSE)
    if (is.integer(from)) {
      to <- as.integer(to)
    }

    df <- data.frame(from = from,
                     to = to,
                     stringsAsFactors = FALSE)

    if (kappa) {
      df$kappa <- unlist(out[i, kappa_cols], use.names = FALSE)
    }
    if (t_inf) {
      df$t_inf <- unlist(out[i, t_inf_cols], use.names = FALSE)
    }

    # Add additional columns
    for (arg in names(args)) {
      map <- args[[arg]]
      df[paste0("from_", arg)] <- map[id_map[df$from]]
      df[paste0("to_", arg)] <- map[id_map[df$to]]
    }

    df
  })

  return(tree_list)
}


#' Estimate the peak of incidence
#'
#' This function estimates the peak of incidence for a given date. If a group is specified,
#' the function will estimate the peak for each unique group.
#'
#' @param date A vector of dates for which to estimate the peak of incidence.
#' @param group An optional vector of group identifiers. If specified, the function will
#' estimate the peak for each unique group.
#'
#' @return A data frame with the group (if specified) and the estimated peak of incidence.
#'
#'
#' @export
get_peak <- function(date, group = NULL) {
  pacman::p_load(incidence)
  if (is.null(group)) {
    incid <- incidence(date)
    p <- estimate_peak(incid, n = 100, alpha = 0.05)
    result <- data.frame(group = NA, observed_peak = p$observed)
  } else {
    data <- tibble(date = date, group = group)
    result <- lapply(unique(group), function(g) {
      g_data <- data[data$group == g, ]
      incid <- incidence(g_data$date)
      p <- estimate_peak(incid, n = 100, alpha = 0.05)
      data.frame(group = g, observed_peak = p$observed)
    })
    result <- do.call(rbind, result)
  }
  return(result)
}


# computes pi from gamma and f
pi_formula <- function(gamma, f) {
  numerator <- gamma * f
  denominator <- 1 - f + gamma * f
  return(numerator / denominator)
}

# GGPLOT ------------------------------------------------------------------

theme_noso <- function(day_break = 2,
                       fill = TRUE,
                       col = TRUE,
                       date = TRUE,
                       legend_position = "bottom") {
  list(
    theme(
      strip.background = element_rect(
        colour = "black",
        fill = "#EEEEEE",
        linewidth = 1.2
      ),
      strip.text = element_text(
        size = 14,
        colour = "black",
        face = "bold"
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 0.5
      ),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_line(linewidth = 0.3, color = "grey"),
      panel.grid.minor = element_line(linewidth = 0.075, color = "grey"),
      axis.line = element_line(linewidth = 0.6, color = "black"),
      legend.position = legend_position,
      text = element_text(size = 12)
    ),
    if (fill) {
      scale_fill_manual(
        values = c(hcw = "#FF9300", patient = "#d800d1"),
        labels = c("HCW", "Patient")
      )
    },
    if (col) {
      scale_color_manual(
        values = c(hcw = "#FF9300", patient = "#d800d1"),
        labels = c("HCW", "Patient")
      )
    },
    if (date) {
      cutoff_breaks <- seq(min(linelist$onset), max(linelist$onset), by = day_break)
      day_month <- format(cutoff_breaks, "%d\n%B")
      day <- format(cutoff_breaks, "%d")
      dup <- which(duplicated(format(cutoff_breaks, "%b")))
      day_month[dup] <- day[dup]
      scale_x_date(
        limits = c(min(linelist$onset) - 1, max(linelist$onset) + 1),
        breaks = cutoff_breaks,
        labels = day_month,
        expand = c(0.01, 0.01)
      )
    }
  )
}



# epicurve
epicurve <- function(day_break = 2,
                     legend_position = "bottom",
                     facet = FALSE) {
  case_counts <-
    linelist %>%
    mutate(col = as.character(row_number())) %>%
    group_by(group, onset, col) %>%
    summarise(n = n())

  peaks <- get_peak(date = linelist$onset, group = linelist$group) %>%
    rowwise() %>%
    mutate(max_n = sum(case_counts$n[case_counts$onset == observed_peak]))

  gg <- case_counts %>%
    ggplot(aes(x = onset, y = n, fill = group)) +
    geom_col(
      position = "stack",
      width = 1,
      col = "#494949",
      linewidth = 0.35
    ) +
    geom_segment(
      data = peaks,
      aes(
        x = observed_peak,
        y = c(10, 10),
        xend = observed_peak,
        yend = max_n,
        linetype = "Peak"
      ),
      col = "black",
      linewidth = 1.75,
      linejoin = 'mitre',
      lineend = 'square',
      arrow = arrow(length = unit(0.25, 'cm'), type = 'closed')
    ) +
    geom_segment(
      data = peaks,
      aes(
        x = observed_peak,
        y = c(10, 10),
        xend = observed_peak,
        yend = max_n,
        linetype = "Peak",
        col = group
      ),
      show.legend = FALSE,
      linewidth = 0.5,
      linejoin = 'mitre',
      lineend = 'square',
      arrow = arrow(length = unit(0.25, 'cm'), type = 'closed')
    ) +
    scale_y_continuous(breaks = seq(0, 50, by = 2)) +
    labs(
      x = "Onset",
      y = "Cases",
      fill = "",
      col = "",
      linetype = ""
    ) +
    theme_noso(day_break) +
    theme(
      legend.position = legend_position,
      legend.background = element_rect(color = "black"),
      legend.key=element_blank(),
      legend.key.size = unit(0.5, "cm"),
      legend.margin = margin(0.1, 0.1, 0.1, 0, "cm"),
    ) +
    guides(fill = guide_legend(override.aes = list(size = 0.25)),
           linetype = guide_legend(override.aes = list(linewidth = 0.2,
                                                       size = 0.1)))

  if (facet) {
    gg <- gg +
      facet_wrap(
        ~ "Cases",
        ncol = 1,
        labeller = label_parsed,
        strip.position = "left"
      ) +
      theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(
          angle = 90,
          face = "bold",
          size = 14,
          margin = margin(r = -0.5, l = -1, 0, 0, "pt")
        ),
        panel.spacing.y = unit(0.5, "lines")
      ) +
      labs(y = "")
  }

  return(gg)
}
