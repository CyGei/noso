# Description:
# For a unique transmission tree:
# The uncertainty in rho is derived from the uncertainty in `phi` i.e. the the proportion of cases in the infectee group relative to all cases.
# The uncertainty in the assorativity coefficient (gamma/delta) is derived from the uncertainty in `pi`.
# The draw_* functions (draw_gamma, draw_delta, draw_R0) are used to draw random values of the metric given the uncertainty in `pi` for a unique tree.
# The draw_array function is used to draw random values of the metric for multiple trees. This accounts for the uncertainty in the metric and the uncertainty in the tree.
# Functions workflow:
# *_formula: Compute the metric for a unique tree
# draw_*: Draw random values of the metric for a unique tree
# draw_array: Draw random values of the metric for multiple trees


draw_phi <- function(from,
                     to,
                     levels = NULL,
                     from_id,
                     to_id,
                     n_samples = 1000) {
  linktree:::check_fromto(from, to)
  ttab <- linktree:::ttable(from, to, levels)
  levels <- rownames(ttab)
  # Expected
  n_cases <- table(factor(unique(na.omit(cbind(
    c(from, to), c(from_id, to_id)
  )))[, 1], levels = levels))
  expected_numerator <- n_cases - 1
  expected_denominator <- sum(n_cases) - 1
  expected_ratio <- expected_numerator / expected_denominator

  phi <- vapply(seq_along(levels), function(i) {
    if (expected_denominator < 1 || expected_ratio[i] < 0) {
      rep(NA, n_samples)
    } else {
      rbinom(n = n_samples,
             size = expected_denominator,
             prob = expected_ratio[i]) / expected_denominator
    }
  }, FUN.VALUE = numeric(n_samples))
  dimnames(phi) <- list(NULL, levels)
  return(phi)
}


# We draw uncertainty in rho from the uncertainty in expected_ratio using binomial distribution
draw_rho <- function(from,
                     to,
                     levels = NULL,
                     n_samples = 1000,
                     args) {
  from_id <- args$from_id
  to_id <- args$to_id
  linktree:::check_fromto(from, to)
  ttab <- linktree:::ttable(from, to, levels)
  levels <- rownames(ttab)

  observed_numerator <- diag(ttab)
  observed_denominator <- colSums(ttab)
  observed_ratio <- observed_numerator / observed_denominator
  phi <- draw_phi(from, to, levels, from_id, to_id, n_samples)

  rho <- vapply(levels, function(i) {
    linktree::gamma2delta(observed_ratio[i] / phi[, i])
  }, FUN.VALUE = numeric(n_samples), USE.NAMES = TRUE)

  return(rho)
}

# This is where the uncertainty lies, all subsequent draw_* functions derive their values from draw_pi
draw_pi <- function(from,
                    to,
                    levels = NULL,
                    n_samples = 1000,
                    args = NULL) {
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
  return(pi)
}

draw_gamma <- function(from,
                       to,
                       levels = NULL,
                       n_samples = 1000,
                       args) {
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

draw_delta <- function(from,
                       to,
                       levels = NULL,
                       n_samples = 1000,
                       args) {
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

draw_R0 <- function(from,
                    to,
                    levels = NULL,
                    n_samples = 1000,
                    args) {
  f <- args$f
  from_id <- args$from_id
  to_id <- args$to_id

  ttab <- linktree:::ttable(from, to, levels)
  levels <- rownames(ttab)

  n_cases <- table(factor(unique(na.omit(cbind(
    c(from, to), c(from_id, to_id)
  )))[, 1], levels = levels))

  gamma <- draw_gamma(from, to, levels, n_samples, args["f"])
  diag <- diag(ttab)
  R0 <- vapply(levels, function(l) {
    R0_formula(
      tau = diag[[l]],
      gamma = gamma[, l],
      f = f[[l]],
      n_cases = n_cases[[l]]
    )
  }, FUN.VALUE = numeric(n_samples))

  return(R0)
}


draw_array <- function(from_col,
                       to_col,
                       levels,
                       trees,
                       draw_function,
                       args = list(),
                       n_samples = 1000) {
  # Check inputs
  if (!is.function(draw_function))
    stop("draw_function must be a function")
  if (!is.list(args))
    stop("args must be a list")
  if (length(levels) < 2)
    stop("There must be at least two group levels in the data")
  if (length(unique(vapply(trees, length, FUN.VALUE = integer(1)))) > 1)
    stop("All trees must have the same length")

  N <- n_samples
  L <- length(levels)
  Tr <- length(trees[[1]])
  C <- length(trees)

  dimnames <- list(
    "sample" = NULL,
    "level" = levels,
    "tree" = NULL,
    "cutoff" = NULL
  )
  result <- array(NA_real_, dim = c(N, L, Tr, C), dimnames = dimnames)
  draw_results <-
    furrr::future_map(unlist(trees, recursive = FALSE), function(tree) {
      if ("from_id" %in% names(args) && "to_id" %in% names(args)) {
        args$from_id <- tree[[args$from_id]]
        args$to_id <- tree[[args$to_id]]
      }

      arg_list <- list(
        from = tree[[from_col]],
        to = tree[[to_col]],
        levels = levels,
        n_samples = N,
        args = args
      )

      do.call(draw_function, arg_list)
    }, .options = furrr::furrr_options(seed = TRUE))
  result[] <- unlist(draw_results)

  return(result)
}


# dims refers to the grouping dimensions we want to compute the CrIs for
draw_CrI <- function(array, dims, alpha = 0.05) {
  compute_CrI <- function(x, alpha = 0.05) {
    x <- x[!is.na(x)] # Remove NAs
    mean_x <- mean(x)
    quantiles <- quantile(x, probs = c(alpha / 2, 1 - alpha / 2))
    return(c(
      mean = mean_x,
      lwr = quantiles[[1]],
      upr = quantiles[[2]]
    ))
  }
  CrI_data <- apply(array, dims, compute_CrI, alpha = alpha)
  return(CrI_data)
}



# test <- future_map(seq_along(cutoff_dates), function(cutoff) {
#   lapply(seq_along(trees[[1]]), function(step) {
#     draw_delta(
#       from = trees[[cutoff]][[step]]$from_group,
#       to = trees[[cutoff]][[step]]$to_group,
#       args = list(f = c(
#         "hcw" = 0.33, "patient" = 0.67
#       ))
#     ) %>% as.data.frame()
#   }) %>% bind_rows(.id = "step") %>% mutate(cutoff_date = cutoff_dates[cutoff])
# }, .options = furrr_options(seed = TRUE)) %>%
#   bind_rows(.id = "cutoff")


# test %>%
#   pivot_longer(cols = c("hcw", "patient"),
#                names_to = "level", values_to = "value") %>%
#   ggplot(aes(x = cutoff_date,
#              y = value,
#              color = level, group = interaction(level, cutoff_date))) +
#   geom_violin()

# draw_delta(
#   from = trees[[1]][[1]]$from_group,
#   to = trees[[1]][[1]]$to_group,
#   args = list(
#     f = c("hcw" = 0.33, "patient" = 0.67)
#   )
# )
