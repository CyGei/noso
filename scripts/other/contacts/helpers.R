# Contacts are loaded in realtime based on the IDs present in the linelist up to time t

pacman::p_load(here, dplyr, purrr)


# Load Data ---------------------------------------------------------------
loadcutoff_dates <- function(paper) {
  input_path <- here::here("data", paper, "input/")
  linelist <- readRDS(here(input_path, "linelist.rds"))
  cutoff_dates <- sort(unique(linelist$onset))
  return(cutoff_dates)
}

loadcontacts <- function(paper) {
  cutoff_dates <- loadcutoff_dates(paper)
  lapply(cutoff_dates, \(x) loadcontacts_helper(paper, x))
}


loadcontacts_helper <- function(paper, cutoff_date) {
  input_path <- here::here("data", paper, "input/")
  linelist <- readRDS(here(input_path, "linelist.rds")) %>%
    filter(onset <= cutoff_date)

  contacts <- readRDS(here(input_path, "contacts.rds")) %>%
    filter(from_id %in% linelist$case_id &
             to_id %in% linelist$case_id) %>%
    mutate(
      from_group = ifelse(startsWith(from_id, "C"), "patient", "hcw"),
      to_group   = ifelse(startsWith(to_id, "C"), "patient", "hcw")
    ) %>%
    filter(from_id != to_id) %>%
    mutate(pair = paste(pmin(from_id, to_id), pmax(from_id, to_id))) %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(-pair)

}



# Formula Helpers ---------------------------------------------------------
contact_delta <- function(pi, f) {
  numerator <- pi - f
  denominator <- pi + (f * (1 - 2 * pi))
  return(numerator / denominator)
}

contact_pi <- function(ttab, n_samples = 1000) {
  tibble(
    level = c("hcw", "patient"),
    success = c(diag(ttab)[["hcw"]], diag(ttab)[["patient"]]),
    trial = c(colSums(ttab)[["hcw"]], rowSums(ttab)[["patient"]])
  ) %>%
    mutate(draws = map2(success, trial, ~ rbeta(n_samples, .x + 1, .y - .x + 1))) %>%
    select(level, draws) %>%
    pivot_wider(names_from = level, values_from = draws) %>%
    unnest(cols = c(hcw, patient))
}



count_contacts <- function(contacts) {
  unique_ids <- unique(c(contacts$from_id, contacts$to_id))

  # Count contacts type for a single ID
  count_for_id <- function(id) {
    contacts %>%
      filter(from_id == id | to_id == id) %>%
      mutate(contact_level = if_else(from_id == id, to_group, from_group)) %>%
      count(contact_level) %>%
      complete(contact_level = c("hcw", "patient"), fill = list(n = 0)) %>%
      mutate(id = id)

  }
  # Calculate frequencies
  contact_counts <- map_df(unique_ids, count_for_id) %>%
    group_by(id) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup()

  return(contact_counts)
}

