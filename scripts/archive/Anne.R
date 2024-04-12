# Random Tree Method ------------------------------------------------------
# Function to sample an infector for each infectee
sample_infector <- function(infectee, linelist) {
  possible_infectors <- linelist$case_id[linelist$case_id != infectee]
  sample(possible_infectors, size = 1)
}

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

linelist$group <- ifelse(grepl("^C", linelist$case_id), "patient", "hcw")
#linelist$group <- linelist$patient_cat

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
random_trees <- lapply(trees, build_random_tree, linelist = linelist)

#2: compute the proportion of infectors for each infectee group (in the observed and random trees)
expected <-
  lapply(
    random_trees,
    compute_infector_proportions,
    from.group = NULL
  ) %>%
  bind_rows(.id = "tree") %>%
  mutate(source = "expected")

observed <- lapply(trees,
                  compute_infector_proportions,
                  from.group = NULL
                   ) %>%
  bind_rows(.id = "tree") %>%
  mutate(source = "observed")



# Plotting ---------------------------------------------------------------
# Plot if you don't specify `from.group` in `compute_infector_proportions`
bind_rows(expected, observed) %>%
  mutate(
    from_to = paste(from, to, sep = " -> ")
  ) %>%
  ggplot(aes(x = proportion, fill = source)) +
  facet_wrap(~from_to, scales = "free_y") +
  geom_histogram(aes(y=..density..),
                 binwidth = 0.025,
                 bins = 40,
                 position = "identity",
                 alpha = 0.5,
                 col = NA) +
  scale_fill_manual(values = c("expected" = "#4399AE", "observed" = "#AE4356")) +
  theme_bw() +
  labs(x = "Proportion of infectors", fill = "Source")



bind_rows(expected, observed) %>%
  mutate(
    from_to = paste(from, to, sep = " -> ")
  ) %>%
  select(tree, from_to, source, proportion) %>%
  pivot_wider(names_from = source, values_from = proportion) %>%
  mutate(ratio = observed / expected) %>%
  filter(!is.infinite(ratio)) %>%
  group_by(from_to) %>%
  summarise(across(ratio, .fns = list(
    mean = ~mean(., na.rm = TRUE),
    lwr = ~quantile(., 0.025, na.rm = TRUE),
    upr = ~quantile(., 0.975,  na.rm = TRUE)
  ))) %>%
  ggplot(aes(x = from_to, y = ratio_mean, ymin = ratio_lwr, ymax = ratio_upr)) +
  geom_pointrange() +
  theme_bw() +
  coord_flip()


# Plot if you specify `from.group` in `compute_infector_proportions`
bind_rows(expected, observed) %>%
  ggplot(aes(x =  proportion, fill = source)) +
  geom_histogram(aes(y=..density..),
                 binwidth = 0.01,
                 position = "identity",
                 alpha = 0.5,
                 col = "grey") +
  theme_bw() +
  scale_fill_manual(values = c("expected" = "#4399AE", "observed" = "#AE4356")) +
  labs(x = "Proportion of infectors", fill = "Source")

bind_rows(expected, observed) %>%
  pivot_wider(names_from = source, values_from = proportion) %>%
  mutate(ratio = observed / expected) %>%
  summarise(across(ratio, .fns = list(
    mean = ~mean(.),
    lwr = ~quantile(., 0.025),
    upr = ~quantile(., 0.975)
  )))




# Analytical Method -------------------------------------------------------
unique_groups <- unique(linelist$group)
ttabs <- lapply(trees, function(tree){
  linktree:::ttable(from = tree$from_group, to = tree$to_group, levels = unique_groups)
})

calculate_ratio <- function(ttab, from, to) {

  # Observed
  observed_numerator <- ttab[from, to]
  observed_denominator <- sum(ttab[, to])
  observed_ratio <- observed_numerator / observed_denominator

  # Expected
  expected_numerator <- sum(ttab[, from]) - 1
  expected_denominator <- sum(ttab) - 1
  expected_ratio <- expected_numerator / expected_denominator

  ratio <- observed_ratio / expected_ratio
  return(ratio)
}

group_grid <- expand.grid(from = unique_groups, to = unique_groups)

calculate_all_ratios <- function(ttab, group_grid) {
  results <- apply(group_grid, 1, function(row) {
    ratio <- calculate_ratio(ttab, row['from'], row['to'])
    data.frame(from = row['from'], to = row['to'], ratio = ratio,
               row.names = NULL)
  })
  do.call(rbind, results)
}
all_ratios <- lapply(ttabs, calculate_all_ratios, group_grid) %>%
  bind_rows(.id = "tree")
all_ratios %>%
  group_by(from, to) %>%
  drop_na() %>%
  filter(!is.infinite(ratio)) %>%
  summarise(across(ratio, .fns = list(
    mean = ~mean(.),
    lwr = ~quantile(., 0.025),
    upr = ~quantile(., 0.975)
  ))) %>%
  ggplot(aes(x = from, y = ratio_mean, ymin = ratio_lwr, ymax = ratio_upr)) +
  facet_wrap(~to) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(limits = c(0, 3)) +
  theme_bw() +
  coord_flip()+
  labs(x = "From group", y = "Ratio of observed to expected infector proportion")
