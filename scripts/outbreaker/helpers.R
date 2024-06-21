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
  si_mu <- 3.0
  si_sd <- 4.1
  si_cv <- si_sd / si_mu
  si_params <- gamma_mucv2shapescale(si_mu, si_cv)
  si <- distcrete::distcrete(
    "gamma",
    shape = si_params$shape,
    scale = si_params$scale,
    w = 0.5,
    interval = 1
  )


  incub_mu <- 5.95
  incub_sd <- 4.31
  incub_cv <- incub_sd / incub_mu
  incub_params <- gamma_mucv2shapescale(incub_mu, incub_cv)
  incub <- distcrete::distcrete(
    "gamma",
    shape = incub_params$shape,
    scale = incub_params$scale,
    w = 0.5,
    interval = 1
  )

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
    filter(onset <= cutoff_date)
  # %>%
  #   arrange(onset)
  # %>%
  #   mutate(date_integer = as.integer(difftime(onset, min(onset), units = "days")))

  ctd <- readRDS(here(input_path, "contacts.rds")) %>%
    filter(from_id %in% linelist$case_id &
             to_id %in% linelist$case_id) %>%
    filter(from_id != to_id)
  # %>%
  #   #unique combinations of contacts
  #   mutate(combination = paste(pmin(from_id, to_id), pmax(from_id, to_id))) %>%
  #   distinct(combination, .keep_all = TRUE) %>%
  #   select(-combination)

  # if eLife paper remove id C123 in the contacts
  if (paper == "eLife2022") {
    ctd <- ctd %>% filter(from_id != "C123" & to_id != "C123")
  }

  if (nrow(ctd) <= 1) {
    ctd <- NULL
  } else{
    ctd <- epicontacts::make_epicontacts(
      linelist = linelist,
      contacts = ctd,
      id = "case_id",
      from = "from_id",
      to = "to_id",
      directed = TRUE
    )
  }


  dna <- read.dna(here(input_path, "dna.fasta"),
                  skip = 0,
                  format = "fasta")
  dna <- dna[match(linelist$case_id, rownames(dna)), ]

  # outbreaker format ----------------------------------------------------------------
  return(
    list(
      data = outbreaker_data(
        ids = linelist$case_id,
        dates = linelist$onset,
        dna = dna,
        ctd = ctd,
        w_dens = si$d(1:20),
        f_dens = incub$d(1:20)
      ),
      config = create_config(config_list),
      priors = custom_priors(pi = prior_list$pi),
      likelihoods = custom_likelihoods()
    )
  )
}
