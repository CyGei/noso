source(here::here("scripts", "helpers.R"))
source(here::here("scripts/other/truncated_tree", "helpers.R"))
loadbreaker()

cutoff_dates <- unique(c(seq(min(as.Date(linelist$onset_inferred)),
                             max(as.Date(linelist$onset_inferred)),
                             by = 7),
                         max(linelist$onset_inferred)))
burnin <- 500

############################################################################
# trunc_trees
############################################################################
# In the below we 1st run the tree on the entire dataset and then we cut the tree at different dates
set.seed(123)
out <- outbreaker2::outbreaker(
  data = o2_data,
  config = config_list,
  priors = custom_priors(
    pi = prior_list$pi$fn
  )
)
out <- out %>%
  filter(step > burnin) %>%
  identify(ids = o2_data$ids)

trees <- get_trees(out = out,
                   ids = linelist$case_id,
                   date = linelist$onset_inferred)
trunc_trees <- lapply(cutoff_dates, function(cutoff_date) {
  lapply(trees, function(tree) {
    tree_cut <- cut_tree_by_date(tree, cutoff_date)
    if (nrow(tree_cut) == 0) {
      return(NULL)
    }
    return(tree_cut)
  })
})

saveRDS(trunc_trees, here("data/other/truncated_tree", "trunc_trees.rds"))

# Reconstruction ---------------------------------------------------------
# here we reconstruct the tree using the o2_data up to the cutoff date
# set.seed(123)
realtime_trees <- lapply(cutoff_dates, function(cutoff_date) {
  filtered_indices <- linelist$onset_inferred <= cutoff_date
  if(sum(filtered_indices) < 2) {
    return(NULL)
  }

  o2_data_cut <- outbreaker2::outbreaker_data(
    dates = linelist$onset_inferred[filtered_indices],
    ids = linelist$case_id[filtered_indices],
    dna = o2_data$dna[filtered_indices, ],
    w_dens = distribution_list$incubation_period$dist$d(1:50),
    f_dens = distribution_list$serial_interval$dist$d(1:50)
  )
  out_cut <-  outbreaker2::outbreaker(
    data = o2_data_cut,
    config = config_list,
    priors =  outbreaker2::custom_priors(
      pi = prior_list$pi$fn
    )
  )
  out_cut <- out_cut %>% filter(step > burnin)
  out_cut <- identify(out_cut, ids = o2_data_cut$ids)
  trees_cut <- get_trees(out_cut,
                         ids = o2_data_cut$ids,
                         date = o2_data_cut$dates)
  return(trees_cut)
})
saveRDS(realtime_trees, here("data/other/truncated_tree", "realtime_trees.rds"))


