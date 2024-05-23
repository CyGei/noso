


# Run multiple iterations of the outbreaker2 model in parallel
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
  source(here::here("scripts", "helpers.R"))
  pacman::p_load(dplyr, ggplot2)
  ggplot(data = bind_rows(furrrbreaker_results, .id = "sim"),
         aes(x = step, y = {{y}}, col = sim)) +
    geom_line(alpha = 0.5, linewidth = 0.5) +
    theme_noso(fill = FALSE, col = FALSE) +
    theme(legend.position = "none")
}




# Load the data/parameters to be used in the outbreaker2 model
loadbreaker <- function() {
  pacman::p_load(ape, outbreaker2, distcrete, epitrix, here)

  # Parameters ---------------------------------------------------------
  distribution_list <- list(
    serial_interval = list(mu = 3.0, sd = 4.1),
    incubation_period = list(mu = 5.95, sd = 4.31)
  ) %>%
    lapply(., function(x) {
      x$cv <- x$sd / x$mu
      x$shape_scale <- epitrix::gamma_mucv2shapescale(x$mu, x$cv)
      x$dist <- distcrete::distcrete(
        "gamma",
        shape = x$shape_scale$shape,
        scale = x$shape_scale$scale,
        w = 0.5,
        interval = 1
      )
      return(x)
    })


  # Prior -------------------------------------------------------------------
  prior_list <- list(pi = list(
    lower_bound_prior = 0.55,
    fn = function(x) {
      ifelse(x$pi > 0.55, log(1 / (1 - 0.55)), log(0))
    }
  ))


  # Config ------------------------------------------------------------------
  config_list <- list(
    n_iter = 10000,
    sample_every = 50,
    max_kappa = 3,
    # missed no cases
    init_kappa = 1,
    # number of generations before the last sampled ancestor
    move_kappa = TRUE,
    init_pi = 0.6,
    # 60% reporting rate
    move_pi = TRUE,
    init_tree = "star",
    outlier_threshold = 2
  )

  # Data -------------------------------------------------------------------
  data_list <- list(
    dna = read.dna(
      here("data/eLife2022", "sequences.fasta"),
      skip = 0,
      format = "fasta"
    ),
    contacts = read.csv(here("data/eLife2022", "contacts.csv")),
    linelist = read.csv(here("data/eLife2022", "linelist.csv"))
  )

  # Reorder DNA sequences
  matching_indices <-
    match(data_list$linelist$case_id, rownames(data_list$dna))
  data_list$dna <- data_list$dna[matching_indices, ]
  # Date as integer
  data_list$linelist$onset_inferred <-
    as.Date(data_list$linelist$onset_inferred, format = "%Y-%m-%d")
  data_list$linelist$date_integer <-
    as.integer(difftime(
      data_list$linelist$onset_inferred,
      min(data_list$linelist$onset_inferred),
      units = "days"
    ))
  data_list$linelist$group <- ifelse(grepl("^C", data_list$linelist$case_id), "patient", "hcw")


  # outbreaker format ----------------------------------------------------------------
  o2_data <- outbreaker_data(
    dates = data_list$linelist$onset_inferred,
    dna = data_list$dna,
    ids = data_list$linelist$case_id,
    w_dens = distribution_list$incubation_period$dist$d(1:50),
    f_dens = distribution_list$serial_interval$dist$d(1:50)
  )

  cutoff_dates <- unique(c(
    seq(min(as.Date(data_list$linelist$onset_inferred)),
        max(as.Date(data_list$linelist$onset_inferred)),
        by = 7
    ),
    max(data_list$linelist$onset_inferred)
  ))

  # assign o2_data, config_list and prior_list to the global environment
  list2env(list(
    linelist = data_list$linelist,
    o2_data = o2_data,
    config_list = config_list,
    prior_list = prior_list,
    cutoff_dates = cutoff_dates
  ),
  envir = .GlobalEnv)

}

# Chi-square test ---------------------------------------------------------

get_chisq <- function(x, y, ids, ntrees = 190) {
  if(is.null(x)) {
    x <- data.frame(from = character(0), to = character(0))
  }
  if(is.null(y)) {
    y <- data.frame(from = character(0), to = character(0))
  }
  tabx <- linktree:::ttable(
    from = x$from,
    to = x$to,
    levels = ids
  ) %>%
    as.data.frame() %>%
    select(from, to, Freq) %>%
    rename(Freq_x = Freq)
  taby <- linktree:::ttable(
    from = y$from,
    to = y$to,
    levels = ids
  ) %>%
    as.data.frame() %>%
    select(from, to, Freq) %>%
    rename(Freq_y = Freq)
  tab <- merge(tabx, taby, by = c("from", "to")) %>%
    filter(from != to) %>%
    filter(!(Freq_x == 0 & Freq_y == 0)) %>%
    select(Freq_x, Freq_y) %>%
    mutate(across(everything(), \(x) (x / ntrees) * 100))

  tryCatch(
    {
      if (nrow(tab) <= 5) {
        test <- stats::fisher.test(tab)$p.value
      } else {
        test <- chisq.test(tab,
                           simulate.p.value = TRUE)$p.value
      }
      return(test)
    },
    error = function(e) {
      return(NA)
    }
  )
}



