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
        meanlog <- log(x$mu ^ 2 / sqrt(x$sd ^ 2 + x$mu ^ 2))
        sdlog <- sqrt(log(1 + (x$sd ^ 2 / x$mu ^ 2)))
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
  prior_list <- list(
    pi = function(x) {
      ifelse(x$pi > min_pi, log(1 / (1 - min_pi)), log(0))
    }
  )

  # Config ------------------------------------------------------------------
  max_kappa <- ifelse(paper == "eLife2022", 3, 2)
  init_pi <- ifelse(paper == "eLife2022", 0.60, 0.80)
  config_list <- list(
    n_iter = 5e5,
    sample_every = 500,
    max_kappa = max_kappa,
    init_kappa = 1,
    move_kappa = TRUE,
    init_pi = init_pi,
    move_pi = TRUE,
    init_tree = "star",
    outlier_threshold = 2
  )


  # Data -------------------------------------------------------------------
  input_path <- here::here("data", paper, "input/")

  linelist <- readRDS(here(input_path, "linelist.rds")) %>%
    filter(onset <= cutoff_date) %>%
    arrange(onset) %>%
    mutate(date_integer = as.integer(difftime(onset, min(onset), units = "days")))

  ctd <- readRDS(here(input_path, "contacts.rds")) %>%
    filter(from_id %in% linelist$case_id & to_id %in% linelist$case_id) %>%
    filter(from_id != to_id) %>%
    #unique combinations of contacts
    mutate(combination = paste(pmin(from_id, to_id), pmax(from_id, to_id))) %>%
    distinct(combination, .keep_all = TRUE) %>%
    select(-combination)

  dna <- read.dna(here(input_path, "dna.fasta"), skip = 0, format = "fasta")
  dna <- dna[match(linelist$case_id, rownames(dna)), ]

  # outbreaker format ----------------------------------------------------------------
  return(list(
    data = outbreaker_data(
      dates = linelist$date_integer,
      dna = dna,
      ids = linelist$case_id,
      ctd = if (nrow(ctd) > 1) ctd else NULL,
      w_dens = distribution_list$incubation_period$dist$d(1:50),
      f_dens = distribution_list$serial_interval$dist$d(1:50)
    ),
    config = create_config(config_list),
    priors = custom_priors(pi = prior_list$pi)
  ))
}
