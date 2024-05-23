pacman::p_load(outbreaker2)
furrrbreaker <- function(data = outbreaker_data(),
                       config = create_config(),
                       priors = custom_priors(),
                       likelihoods = custom_likelihoods(),
                       moves = custom_moves(),
                       n_sim = 1,
                       plan = "multisession",
                       workers = NULL) {
  pacman::p_load(outbreaker2, furrr)

  if (is.null(workers)) {
    workers <- future::availableCores() - 1
  }

  future::plan(plan, workers = workers)

  results <-
    furrr::future_map(
      1:n_sim,
      ~ outbreaker2::outbreaker(data, config, priors, likelihoods, moves),
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )

  future::plan("sequential")

  return(results)
}
