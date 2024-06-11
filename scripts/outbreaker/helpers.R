# Load the data/parameters to be used in the outbreaker2 model
# specified in the respective paper
pacman::p_load(ape, outbreaker2, distcrete, epitrix, here, dplyr, purrr, furrr)

loadcutoff_dates <- function(paper) {
  input_path <- here::here("data", paper, "input/")
  linelist <- readRDS(here(input_path, "linelist.rds"))
  cutoff_dates <- sort(unique(linelist$onset))
  return(cutoff_dates)
}

loadbreaker <- function(paper) {
  cutoff_dates <- loadcutoff_dates(paper)
  lapply(cutoff_dates, \(x) loadbreaker_helper(paper, x))
}


loadbreaker_helper <- function(paper, cutoff_date) {
  # Parameters ---------------------------------------------------------
  distribution_list <- list(
    serial_interval = list(mu = 3.0, sd = 4.1),
    incubation_period = list(mu = 5.95, sd = 4.31)
  ) %>%
    purrr::imap(function(x, name) {
      if (name == "serial_interval") {
        x$cv <- x$sd / x$mu
        x$shape_scale <- epitrix::gamma_mucv2shapescale(x$mu, x$cv)
        x$dist <- distcrete::distcrete(
          "gamma",
          shape = x$shape_scale$shape,
          scale = x$shape_scale$scale,
          w = 0.5,
          interval = 1
        )
      } else if (name == "incubation_period") {
        meanlog <- log(x$mu^2 / sqrt(x$sd^2 + x$mu^2))
        sdlog <- sqrt(log(1 + (x$sd^2 / x$mu^2)))
        x$dist <- distcrete::distcrete(
          "lnorm",
          meanlog = meanlog,
          sdlog = sdlog,
          w = 0.5,
          interval = 1
        )
      }

      x
    })


  # Prior -------------------------------------------------------------------
  min_pi <- ifelse(paper == "eLife2022", 0.55, 0.75)
  prior_list <- list(pi = function(x) {
    ifelse(x$pi > min_pi, log(1 / (1 - min_pi)), log(0))
  })

  # Config ------------------------------------------------------------------
  max_kappa <- ifelse(paper == "eLife2022", 3, 2)
  config_list <- list(
    n_iter = 5e5,
    sample_every = 500,
    max_kappa = max_kappa,
    # missed no cases
    init_kappa = 1,
    # number of generations before the last sampled ancestor
    move_kappa = TRUE,
    init_pi = 0.8,
    move_pi = TRUE,
    init_tree = "star",
    outlier_threshold = 2
  )


  # Data -------------------------------------------------------------------
  input_path <- here::here("data", paper, "input/")
  data_list <- list(
    dna = read.dna(
      here(input_path, "dna.fasta"),
      skip = 0,
      format = "fasta"
    ),
    linelist = readRDS(here(input_path, "linelist.rds")) %>%
      filter(onset <= cutoff_date) %>%
      arrange(onset)
  )

  # Reorder DNA sequences
  matching_indices <-
    match(data_list$linelist$case_id, rownames(data_list$dna))
  # Remove unmatched sequences
  data_list$dna <- data_list$dna[matching_indices, ]
  # Date as integer
  # data_list$linelist$onset <-
  #   as.Date(data_list$linelist$onset, format = "%Y-%m-%d")
  data_list$linelist$date_integer <-
    as.integer(difftime(
      data_list$linelist$onset,
      min(data_list$linelist$onset),
      units = "days"
    ))


  # outbreaker format ----------------------------------------------------------------
  return(
    list(
      data = outbreaker_data(
        dates = data_list$linelist$date_integer,
        dna = data_list$dna,
        ids = data_list$linelist$case_id,
        w_dens = distribution_list$incubation_period$dist$d(1:50),
        f_dens = distribution_list$serial_interval$dist$d(1:50)
      ),
      config = create_config(config_list),
      priors = custom_priors(pi = prior_list$pi)
    )
  )
}

