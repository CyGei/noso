# Assess chi-square -------------------------------------------------------
# Given that the chains have all converged, we expect that there is no significant difference
# between the trees. To test this hypothesis, we will use a chi-square test on the ancestry frequencies
# of any two chains.
source(here::here("scripts", "helpers.R"))
load_data()
load_libraries()
out <- readRDS(here::here("data/other/truncated_tree", "out.rds"))
trees <- lapply(out, \(chain){
  chain <- chain %>%
    filter(step > 500) %>%
    identify(., ids = linelist$case_id)
  get_trees(chain, ids = linelist$case_id)
})

plan(multisession, workers = availableCores() - 1)
# Test the chi-square on the same chain using bootstrapping
boot_chain <- function(chain, ntrees = 190) {
  sample(chain, size = length(chain), replace = TRUE)
}
set.seed(123)
future_map(1:length(trees), \(i){
  get_chi.test(
    x = bind_rows(trees[[i]]),
    y = bind_rows(boot_chain(trees[[i]])),
    ids = linelist$case_id,
    ntrees = 190
  )
}, .options = furrr_options(seed = TRUE)) %>%
  unlist() %>% summary()

# Test the chi-square on different chains
grid <- expand.grid(1:length(trees), 1:length(trees)) %>%
  filter(Var1 != Var2) %>%
  filter(Var1 < Var2)
set.seed(123)
grid <- sample_n(grid, 50)
grid_results <- future_map(1:nrow(grid), \(i){
  get_chi.test(
    x = bind_rows(trees[[grid$Var1[i]]]),
    y = bind_rows(trees[[grid$Var2[i]]]),
    ids = linelist$case_id,
    ntrees = 190
  )
}, .options = furrr_options(seed = TRUE))
grid_results %>%
  unlist() %>%
  summary()
