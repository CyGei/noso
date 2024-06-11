source(
  here::here(
    "scripts/other/truncated_tree",
    "helpers.R"
  )
)
#load the data
loadbreaker()
set.seed(123)
system.time({
  progressr::with_progress({
    out <- furrrbreaker(
      data = o2_data,
      config = config_list,
      priors = custom_priors(pi = prior_list$pi$fn),
      n_sim = 100
    )
  })
})
#saveRDS(out, here::here("data/other/truncated_tree", "out.rds"))
out <- readRDS(here::here("data/other/truncated_tree", "out.rds"))

cowplot::plot_grid(
  furrrplot(out, post) + labs(x = ""),
  furrrplot(out, like) + labs(x = ""),
  furrrplot(out, pi),
  furrrplot(out, mu),
  nrow = 2,
  align = "v",
  labels = "AUTO"
)

