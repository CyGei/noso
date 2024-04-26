pacman::p_load(outbreaker2)
furrrbreaker <- function(data = outbreaker_data(),
                         config = create_config(),
                         priors = custom_priors(),
                         likelihoods = custom_likelihoods(),
                         moves = custom_moves(),
                         n_sim = 1,
                         plan = "multisession",
                         workers = NULL) {
  pacman::p_load(outbreaker2, furrr, progressr)

  if (is.null(workers)) {
    workers <- future::availableCores() - 1
  }

  future::plan(plan, workers = workers)

  p <- progressr::progressor(steps = n_sim)


  results <-
    furrr::future_map(
      1:n_sim,
      ~ {
        p()
        outbreaker2::outbreaker(data, config, priors, likelihoods, moves)
      },
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )

  future::plan("sequential")

  return(results)
}


furrrplot <- function(furrrbreaker_results, y) {
  pacman::p_load(dplyr, ggplot2)
  ggplot(
    data = bind_rows(furrrbreaker_results, .id = "sim"),
    aes(x = step, y = {{y}}, col = sim)) +
    geom_line()
}


#
#
# source(here("scripts", "other", "outbreaker.R"))
# set.seed(123)
# progressr::with_progress({
#   out <- furrrbreaker(
#     data = o2_data,
#     config = config_list,
#     priors = custom_priors(pi = prior_list$pi$fn),
#     n_sim = 10
#   )
# })
# furrrplot(out, post)
# furrrplot(out, pi)
# furrrplot(out, like)
# furrrplot(out, mu)
