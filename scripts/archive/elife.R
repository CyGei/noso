# The purpose of this script is to reproduce Mo's
# results from e-life paper using the monte-carlo approach.

# Random Tree Method ------------------------------------------------------
# Function to sample a random infector for each infectee
sample_infector <- function(infectee, linelist) {
  possible_infectors <- linelist$case_id[linelist$case_id != infectee]
  sample(possible_infectors, size = 1)
}

# Function to build a random tree
# It takes a posterior tree and a linelist as input, and returns a new tree
# where for each infectee, a random infector is sampled from the linelist
build_random_tree <- function(tree, linelist) {

  # exclude importations
  tree <- tree[complete.cases(tree),]

  infectees <- tree$to

  # Sample randomly from the potential infectors for each infectee (excluding the infectee themselves)
  random_infectors <- vapply(X = infectees,
                             FUN = function(e) sample_infector(e, linelist),
                             FUN.VALUE = character(1),
                             USE.NAMES = FALSE)

  # Create a new data frame with randomly sampled infectors
  random_tree <- data.frame(
    from = random_infectors,
    to = infectees,
    from_group = linelist[match(random_infectors, linelist$case_id), "group"],
    to_group = tree$to_group
  )

  return(random_tree)
}

# Function to compute the proportion of infectors for each infectee group
# It takes a tree as input and returns a data frame with the proportion of infectors for each infectee group
# If from.group is specified, it computes the proportion of infectors for that group relative to the tree
compute_infector_proportions <- function(tree,
                                         from.group = NULL) {
  tab <- table(tree$from_group, tree$to_group)

  #from.group is to compute for the total excess in the tree
  if (!is.null(from.group)) {
    tab <- rowSums(tab) / sum(tab)
    tab <- tab[from.group]
    names(tab) <- "proportion"
    return(tab)
  }

  # The % of infectors for each infectee group
  ptab <- as.data.frame(prop.table(tab, margin = 2))
  names(ptab) <- c("from", "to", "proportion")
  return(ptab)
}


# Analysis ----------------------------------------------------------------
source(here::here("scripts", "helpers.R"))
load_libraries()
load_data()

#linelist$group <- ifelse(grepl("^C", linelist$case_id), "patient", "hcw")
linelist <- linelist %>%
  mutate(group = case_when(
    patient_cat  == "patient-HA" ~ "patient_noso",
    patient_cat == "patient-CA" ~ "patient_community",
    patient_cat == "hcw covid" ~ "hcw_covid",
    patient_cat == "hcw non-covid" ~ "hcw_outbreak",
  ))

outbreaker_chains <- outbreaker_chains %>%
  filter(step > 500) %>%
  identify(ids = linelist$case_id)
outbreaker_chains <- filter_alpha_by_kappa(outbreaker_chains, 1)

trees <- get_trees(
  out = outbreaker_chains,
  id = linelist$case_id,
  group = linelist$group,
  date = linelist$onset_inferred
)

#make sure trees are complete cases/pairs (i.e. exclude importations)
trees <- lapply(trees, function(tree) tree[complete.cases(tree),])


#1: build random trees (i.e. randomly sample an infector for each infectee)
random_trees <- lapply(trees, \(x) build_random_tree(x, linelist))

#2: compute the proportion of infectors for each infectee group
expected <-
  lapply(random_trees,
         \(x) compute_infector_proportions(x, from.group = NULL)) %>%
  bind_rows(.id = "tree") %>%
  mutate(source = "expected")

observed <- lapply(trees,
                   \(x) compute_infector_proportions(x, from.group = NULL)
                   ) %>%
  bind_rows(.id = "tree") %>%
  mutate(source = "observed")



# Plotting ---------------------------------------------------------------
# Plot if you don't specify `from.group` in `compute_infector_proportions`
modat <- bind_rows(expected, observed) %>%
  filter(from == to) %>%
  filter(from != "patient_community" &
           to != "patient_community") %>%
  mutate(from_to = paste(from, to, sep = " -> "))

modat %>%
  ggplot(aes(x = proportion, fill = source)) +
  facet_wrap( ~ from_to, scales = "free_y",
              ncol = 1) +
  geom_histogram(
    aes(y = ..density..),
    binwidth = 0.025,
    bins = 40,
    position = "identity",
    alpha = 0.5,
    col = "black"
  ) +
  scale_fill_manual(values = c(
    "expected" = "#4399AE",
    "observed" = "#AE4356"
  )) +
  scale_x_continuous(breaks = seq(0, 1, 0.25),
                     limits = c(0,1)) +
  theme_bw() +
  labs(x = "Proportion of infectors", fill = "Source")


#ratios MC
mo_ratio_MC <- modat %>%
  select(tree, from_to, source, proportion) %>%
  pivot_wider(names_from = source, values_from = proportion) %>%
  mutate(ratio = observed / expected) %>%
  filter(ratio != Inf) %>% # remove Inf values
  group_by(from_to) %>%
  summarise(
    mean = mean(ratio),
    lwr = quantile(ratio, 0.025),
    upr = quantile(ratio, 0.975)
  )

# Analytical Method -------------------------------------------------------
unique_groups <- unique(linelist$group)
unique_groups

ratios <- lapply(trees, \(x) {
  draw_mo(
    from = x$from_group,
    to = x$to_group,
    levels = NULL,
    n_samples = 1000,
    args = list(
      from_id = x$from,
      to_id = x$to,
      diag = FALSE
    )
  ) %>%
    as.data.frame() %>%
    rownames_to_column(var = "from")
}) %>%
  bind_rows() %>%
  pivot_longer(cols = -from, names_to = "to", values_to = "ratio")

mo_ratio_AN <- ratios %>%
  filter(from == to) %>%
  filter(from != "patient_community" &
           to != "patient_community") %>%
  mutate(from_to = paste(from, to, sep = " -> ")) %>%
  group_by(from_to) %>%
  summarise(
    mean = mean(ratio),
    lwr = quantile(ratio, 0.025),
    upr = quantile(ratio, 0.975)
  )

bind_rows(mo_ratio_MC, mo_ratio_AN,
          .id = "source") %>%
  mutate(source = ifelse(source == "1", "MonteCarlo", "Analytical")) %>%
  ggplot(aes(
    x = mean,
    y = from_to,
    xmin = lwr,
    xmax = upr,
    col = source
  )) +
  geom_pointrange(position=ggstance::position_dodgev(height=0.3)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 8, 0.5)) +
  theme_bw()
