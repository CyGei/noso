# Description: The uncertainty in the assorativity coefficient is derived from the uncertainty in `pi`.
# The draw_* functions (draw_gamma, draw_delta, draw_R0) are used to draw random values of the metric given the uncertainty in `pi` for a unique tree.
# The draw_array function is used to draw random values of the metric for multiple trees. This accounts for the uncertainty in the metric and the uncertainty in the tree.
# Note that draw_rho has no uncertainty, as it is a deterministic function of the observed and expected ratios.
# Functions workflow:
# *_formula: Compute the metric for a unique tree
# draw_*: Draw random values of the metric for a unique tree
# draw_array: Draw random values of the metric for multiple trees


rho_formula <- function(from, to, levels = NULL, from_id, to_id, infector, infectee) {
  linktree:::check_fromto(from, to)
  ttab <- linktree:::ttable(from, to, levels)
  levels <- rownames(ttab)

  # Observed
  observed_numerator <- ttab[infector, infectee]
  observed_denominator <- sum(ttab[, infectee])
  observed_ratio <- observed_numerator / observed_denominator

  # Expected
  n_cases <- table(factor(
    unique(na.omit(cbind(
      c(from, to), c(from_id, to_id)
    )))[, 1],
    levels = levels
  ))
  expected_numerator <- n_cases[[infector]] - 1
  expected_denominator <- sum(n_cases) - 1
  expected_ratio <- expected_numerator / expected_denominator

  return(observed_ratio / expected_ratio)
}

# Note there is no uncertainty here: n_samples is ignored (only to be consistent with the other draw_* functions)
# diag is a logical indicating whether the diagonal elements should be returned
# This is for compatibility with the other draw_* functions (comparison with gamma/delta)
draw_rho <- function(from, to, levels = NULL, n_samples = 1, args) {
  from_id <- args$from_id
  to_id <- args$to_id
  diag <- args$diag

  levels <- rownames(linktree:::ttable(from, to, levels))
  grid_levels <- expand.grid(from = levels, to = levels)

  rho <- matrix(
    apply(grid_levels, 1, function(x) {
      rho_formula(from, to, levels, from_id, to_id, x[1], x[2])
    }),
    nrow = length(levels), ncol = length(levels),
    dimnames = list(levels, levels)
  )
  if (diag) {
    rho <- diag(rho)
  }
  return(rho)
}

# This is where the uncertainty lies, all subsequent draw_* functions derive their values from draw_pi
draw_pi <- function(from, to, levels = NULL, n_samples = 1000, args = NULL) {
  linktree:::check_fromto(from, to)
  ttab <- linktree:::ttable(from, to, levels)
  levels <- rownames(ttab)
  success <- diag(ttab)
  trial <- rowSums(ttab)

  pi <- vapply(levels, function(l) {
    i <- which(levels == l)
    if (trial[i] < 1) {
      rep(NA, n_samples)
    } else {
      rbeta(n_samples, success[i] + 1, trial[i] - success[i] + 1)
    }
  }, FUN.VALUE = numeric(n_samples))
  # matrix(pi_values, nrow = n_samples, ncol = length(levels), byrow = FALSE)
  return(pi)
}

draw_gamma <- function(from, to, levels = NULL, n_samples = 1000, args) {
  f <- args$f

  levels <- names(f)
  pi <- draw_pi(from, to, levels, n_samples)
  linktree:::check_flevels(f, levels)
  f <- f[levels]
  f <- f / sum(f, na.rm = TRUE)
  gamma <- vapply(levels, function(i) {
    linktree:::gamma_formula(pi[, i], f[i])
  }, FUN.VALUE = numeric(n_samples))
  return(gamma)
}

draw_delta <- function(from, to, levels = NULL, n_samples = 1000, args) {
  gamma <- draw_gamma(from, to, levels, n_samples, args)
  delta <- apply(gamma, 2, function(x) {
    linktree::gamma2delta(x)
  })
  return(delta)
}

R0_formula <- function(tau, gamma, f, n_cases) {
  numerator <- tau * gamma * f + (1 - f)
  denominator <- gamma * f * n_cases
  return(numerator / denominator)
}

draw_R0 <- function(from, to, levels = NULL, n_samples = 1000, args) {
  f <- args$f
  from_id <- args$from_id
  to_id <- args$to_id

  ttab <- linktree:::ttable(from, to, levels)
  levels <- rownames(ttab)

  n_cases <- table(factor(
    unique(na.omit(cbind(
      c(from, to), c(from_id, to_id)
    )))[, 1],
    levels = levels
  ))

  gamma <- draw_gamma(from, to, levels, n_samples, args["f"])
  diag <- diag(ttab)
  R0 <- vapply(levels, function(l) {
    R0_formula(tau = diag[[l]], gamma = gamma[, l], f = f[[l]], n_cases = n_cases[[l]])
  }, FUN.VALUE = numeric(n_samples))

  return(R0)
}


draw_array <- function(from_col, to_col, levels, trees, draw_function, args = list(), cutoff_dates = NULL, n_samples = 1000) {
  pacman::p_load("furrr")
  # Check inputs
  if (!is.function(draw_function)) stop("draw_function must be a function")
  if (!is.list(args)) stop("args must be a list")
  if (length(levels) < 2) stop("There must be at least two group levels in the data")

  N <- n_samples
  L <- length(levels)
  C <- if (is.null(cutoff_dates)) 1 else length(cutoff_dates)
  Tr <- length(trees)

  if (!is.null(cutoff_dates)) {
    dimnames <- list("sample" = NULL, "level" = levels, "cutoff" = as.character(cutoff_dates), "tree" = NULL)
  } else {
    dimnames <- list("sample" = NULL, "level" = levels, "cutoff" = NULL, "tree" = NULL)
  }
  result <- array(NA_real_, dim = c(N, L, C, Tr), dimnames = dimnames)

 # plan(multisession, workers = future::availableCores() - 10)

  cut_trees <- furrr::future_map(trees, function(tree) {
    if (is.null(cutoff_dates)) {
      list(tree)
    } else {
      return(lapply(cutoff_dates, function(cutoff_date) cut_tree_by_date(tree, cutoff_date)))
    }
  }, .options = furrr::furrr_options(seed = FALSE))

  draw_results <-
    furrr::future_map(unlist(cut_trees, recursive = FALSE),
                      function(cut_tree) {
                        if ("from_id" %in% names(args) && "to_id" %in% names(args)) {
                          args$from_id <- cut_tree[[args$from_id]]
                          args$to_id <- cut_tree[[args$to_id]]
                        }

                        arg_list <- list(
                          from = cut_tree[[from_col]],
                          to = cut_tree[[to_col]],
                          levels = levels,
                          n_samples = N,
                          args = args
                        )

                        do.call(draw_function, arg_list)
                      },
                      .options = furrr::furrr_options(seed = TRUE)
    )
  #plan(sequential)
  result[] <- unlist(draw_results)

  return(result)
}


# dims refers to the grouping dimensions we want to compute the CrIs for
draw_CrI <- function(array, dims, alpha = 0.05) {
  compute_CrI <- function(x, alpha = 0.05) {
    x <- x[!is.na(x)] # Remove NAs
    mean_x <- mean(x)
    quantiles <- quantile(x, probs = c(alpha / 2, 1 - alpha / 2))
    return(c(mean = mean_x, lwr = quantiles[[1]], upr = quantiles[[2]]))
  }
  CrI_data <- apply(array, dims, compute_CrI, alpha = alpha)
  return(CrI_data)
}


# ARCHIVE -------------------------------------------------------------------------------------


# draw_array <- function(from_col, to_col, trees, levels, draw_function, draw_args, cutoff_dates = NULL, n_samples = 1000) {
#   # check if the draw_function is valid
#   if (!is.function(draw_function)) {
#     stop("draw_function must be a function")
#   }
#   # check if the draw_args is a list
#   if (!is.list(draw_args)) {
#     stop("draw_args must be a list")
#   }
#   # check levels
#   if (length(levels) < 2) {
#     stop("There must be at least two group levels in the data")
#   }

#   N <- n_samples
#   L <- length(levels)
#   Tr <- length(trees)
#   C <- if (is.null(cutoff_dates)) 1 else length(cutoff_dates)
#   result <- array(NA, dim = c(N, L, Tr, C), dimnames = list("sample" = NULL, "level" = levels, "tree" = NULL, "cutoff" = NULL))

#   for (t in 1:Tr) {
#     tree <- trees[[t]]
#     for (c in 1:C) {
#       if (is.null(cutoff_dates)) {
#         cut_tree <- tree
#       } else {
#         cutoff_date <- cutoff_dates[c]
#         cut_tree <- cut_tree_by_date(tree, cutoff_date)
#       }
#       draw_result <- do.call(draw_function, c(list(from = cut_tree[[from_col]], to = cut_tree[[to_col]], levels = names(f), n_samples = N), draw_args))
#       result[, , t, c] <- draw_result
#     }
#   }
#   return(result)
# }


# draw_density <- function(array, dims) {
#   compute_density <- function(data) {
#     clean_data <- na.omit(data)
#     if (length(clean_data) < 2) {
#       return(data.frame(x = rep(NA, 512), y = rep(NA, 512)))
#     } else {
#       dens <- density(clean_data, from = -1, to = 1, na.rm = TRUE)
#       return(data.frame(x = dens$x, y = dens$y))
#     }
#   }
#   replace_indices_with_values <- function(indices, values) {
#     if (!is.null(values)) {
#       return(values[indices])
#     } else {
#       return(indices)
#     }
#   }

#   density_data <- apply(array, dims, compute_density)
#   dim_meta <- dimnames(array)[dims]
#   dim_values <- lapply(dim(array)[dims], seq_len)
#   dim_grid <- expand.grid(Map(replace_indices_with_values, dim_values, dim_meta), stringsAsFactors = FALSE)

#   # Set names for dimensions
#   names(dim_grid) <- ifelse(sapply(names(dim_meta), is.null), paste0("dim_", dims), names(dim_meta))
#   # names(dim_grid) <- sapply(seq_along(dim_meta), function(i) {
#   #   if (is.null(dim_meta[[i]])) {
#   #     return(paste0("dim_", dims[i]))
#   #   } else {
#   #     return(dim_meta[[i]])
#   #   }
#   # })

#   density_df <- data.frame(
#     density_x = unlist(lapply(density_data, "[[", "x")),
#     density_y = unlist(lapply(density_data, "[[", "y")),
#     dim_grid
#   )

#   return(density_df)
# }


#
# draw_CrI <- function(array, dims, alpha = 0.05) {
#   compute_CrI <- function(x, alpha = 0.05) {
#     x <- x[!is.na(x)]  # Remove NAs
#     mean_x <- mean(x)
#     quantiles <- quantile(x, probs = c(alpha/2, 1 - alpha/2))
#     return(c(mean = mean_x, lwr = quantiles[[1]], upr = quantiles[[2]]))
#   }
#
#   replace_indices_with_values <- function(indices, values) {
#     if (!is.null(values)) {
#       return(values[indices])
#     } else {
#       return(indices)
#     }
#   }
#
#   CrI_data <- apply(array, dims, compute_CrI, alpha = alpha)
#   dim_meta <- dimnames(array)[dims]
#   dim_values <- lapply(dim(array)[dims], seq_len)
#   dim_grid <- expand.grid(Map(replace_indices_with_values, dim_values, dim_meta), stringsAsFactors = FALSE)
#
#   # Set names for dimensions
#   names(dim_grid) <- ifelse(sapply(names(dim_meta), is.null), paste0("dim_", dims), names(dim_meta))
#
#   CrI_df <- data.frame(
#     mean = unlist(lapply(CrI_data, "[[", "mean")),
#     lwr = unlist(lapply(CrI_data, "[[", "lwr")),
#     upr = unlist(lapply(CrI_data, "[[", "upr")),
#     dim_grid
#   )
#
#   return(CrI_df)
# }
