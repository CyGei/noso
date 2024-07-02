source(here::here("scripts/outbreaker", "helpers.R"))
paper = "eLife2022"
cutoff_dates = loadcutoff_dates(paper)
input <- loadbreaker_helper(paper, cutoff_date = cutoff_dates[length(cutoff_dates)])

x = input
system.time({
  set.seed(123)
  out <- outbreaker2::outbreaker(
    data = x$data,
    config = x$config,
    priors = x$priors,
    likelihoods = x$likelihoods
  )
})
#92 min

saveRDS(out, here::here("data", paper, "output", "out_cyril4.rds"))

out <- readRDS(here::here("data", paper, "output", "out_cyril.rds"))



# Final Check -------------------------------------------------------------
source(here::here("scripts", "helpers.R"))
load_libraries()
load_data(paper)
# Check that our final iteration of outbreaker is the same as Abbas et al.
abbas_out <- readRDS(here::here("data", paper, "input", "out.rds"))
geismar_out <- readRDS(here::here("data", paper, "output", "out_cyril4.rds"))

# Compare the two outputs
abbas_trees <-
  abbas_out %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  filter_alpha_by_kappa(1L) %>%
  get_trees(
    ids = linelist$case_id,
    group = linelist$group,
    date = linelist$onset,
    kappa = TRUE
  ) %>%
  bind_rows(.id = "tree")
abbas_ances <- linktree:::ttable(from = abbas_trees$from,
                           to = abbas_trees$to,
                           level = linelist$case_id) %>%
  as_tibble() %>%
  filter(from != to)

geismar_trees <-
  geismar_out %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id) %>%
  filter_alpha_by_kappa(1L) %>%
  get_trees(
    ids = linelist$case_id,
    group = linelist$group,
    date = linelist$onset,
    kappa = TRUE
  ) %>%
  bind_rows(.id = "tree")

geismar_ances <- linktree:::ttable(from = geismar_trees$from,
                           to = geismar_trees$to,
                           level = linelist$case_id) %>%
  as_tibble() %>%
  filter(from != to)

# Compare the two outputs using the chi-squared test

get_chisq <- function(tabx, taby) {
  tab <- left_join(tabx, taby, by = c("from", "to")) %>%
    filter(!(n.x == 0 & n.y == 0)) %>%
    select(starts_with("n"))

  # Calculate expected values
  total_sum <- sum(tab)
  row_sums <- rowSums(tab)
  col_sums <- colSums(tab)
  expected_values <- outer(row_sums, col_sums, "*") / total_sum
  #the expected values represent the values that we would expect to see in each cell of the table
  #if the null hypothesis of independence between the two variables (rows and columns) is true.
  #(Row Total * Column Total) / Grand Total

  # Determine whether to simulate p-values
  simulate_p <- any(expected_values < 5)

  # Check if there is only one non-zero column
  only_one_non_zero_col <- sum(col_sums > 0) == 1

  # Perform appropriate test
  if (simulate_p && !only_one_non_zero_col) {
    # Perform Fisher's exact test with simulated p-value
    res <- fisher.test(tab, simulate.p.value = TRUE)
    p_value <- res$p.value
  } else if (!simulate_p && !only_one_non_zero_col) {
    # Perform Chi-squared test
    res <- chisq.test(tab)
    p_value <- res$p.value
  } else if (only_one_non_zero_col) {
    # Perform Fisher's exact test without simulating p-value
    res <- fisher.test(tab)
    p_value <- res$p.value
  } else {
    warning("Invalid input data. Cannot perform any test.")
    p_value <- NA
  }

  return(p_value)
}


get_chisq(abbas_ances, geismar_ances)




# test inputs -------------------------------------------------------------
# check that inputs are the same

# Contacts
abbas_contacts <- read.csv(
  "C:\\Users\\cg1521\\OneDrive - Imperial College London\\PROJECTS\\OUTBREAKER2\\NOSO\\data\\eLife2022\\input\\contacts.csv"
) %>% select(-X) %>%
  as_tibble()
geismar_contacts <- readRDS(here(input_path, "contacts.rds"))

# check that the two dataframes are the same
all.equal(abbas_contacts, geismar_contacts)
# CONTACTS OK


# Linelist
abbas_linelist <- readRDS(
  "C:\\Users\\cg1521\\OneDrive - Imperial College London\\PROJECTS\\OUTBREAKER2\\NOSO\\data\\eLife2022\\input\\linelist.rds"
)
geismar_linelist <- readRDS(here(input_path, "linelist.rds"))
all.equal(abbas_linelist, geismar_linelist)
# LINELIST OK


# DNA
abbas_dna <- ape::read.dna(
  "C:\\Users\\cg1521\\OneDrive - Imperial College London\\PROJECTS\\OUTBREAKER2\\NOSO\\data\\eLife2022\\input\\dna.fasta",
  format = "fasta",
  skip = 0
)
geismar_dna <- ape::read.dna(here(input_path, "dna.fasta"), format = "fasta", skip = 0)
all.equal(abbas_dna, geismar_dna)
# DNA OK
all(rownames(geismar_dna) == abbas_linelist$case_id)
