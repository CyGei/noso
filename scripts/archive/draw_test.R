source(here::here("scripts", "helpers.R"))
load_libraries()
load_data()
linelist$group <- ifelse(grepl("^C", linelist$case_id), "patient", "hcw")
outbreaker_chains <- outbreaker_chains %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id)
outbreaker_chains <- filter_alpha_by_kappa(outbreaker_chains, 1)

trees <- get_trees(
  out = outbreaker_chains,
  ids = linelist$case_id,
  group = linelist$group,
  date = linelist$onset_inferred
)

cutoff_dates <- unique(c(seq(min(as.Date(linelist$onset_inferred)),
                             max(as.Date(linelist$onset_inferred)),
                             by = 7),
                         max(linelist$onset_inferred)))

tree = trees[[200]]

mo_formula(
  from = tree$from_group,
  to = tree$to_group,
  levels = NULL,
  from_id = tree$from,
  to_id = tree$to,
  infector = "hcw",
  infectee = "patient"
)

draw_mo(
  from = tree$from_group,
  to = tree$to_group,
  levels = NULL,
  n_samples = 1,
  args = list(
    from_id = tree$from,
    to_id = tree$to,
    diag = TRUE
  )
)

draw_pi(
  from = tree$from_group,
  to = tree$to_group,
  levels = c("hcw", "patient"),
  n_samples = 4
)


draw_gamma(
  from = tree$from_group,
  to = tree$to_group,
  levels = NULL,
  n_samples = 4,
  args = list(
    f = c("hcw" = 0.5, "patient" = 0.5)
  )
)

draw_delta(
  from = tree$from_group,
  to = tree$to_group,
  levels = NULL,
  n_samples = 4,
  args = list(
    f = c("hcw" = 0.5, "patient" = 0.5)
  )
)

draw_R0(
  from = tree$from_group,
  to = tree$to_group,
  levels = c("hcw", "patient"),
  n_samples = 4,
  args = list(
    f = c("hcw" = 0.5, "patient" = 0.5),
    from_id = tree$from,
    to_id = tree$to
  )
)

array <- draw_array(
  from_col = "from_group",
  to_col = "to_group",
  levels = c("hcw", "patient"),
  trees = trees[1:100],
  draw_function = draw_gamma,
  args = list(
    f = c("hcw" = 0.5, "patient" = 0.5),
    from_id = "from",
    to_id = "to",
    diag = TRUE
  ),
  cutoff_dates = cutoff_dates,
  n_samples = 100
)

CrI <- draw_CrI(
  array = array,
  dims = c(2,3)
)
CrI
