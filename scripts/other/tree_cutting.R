# This scripts aims at assessing whether truncating the trees at given dates yields different results compared to reconstructing the tree up to a given date.
source(here::here("scripts", "helpers.R"))
load_libraries()
load_data()
library(ape)
library(outbreaker2)

# Parameters ---------------------------------------------------------
distribution_list <- list(
  serial_interval = list(
    mu = 3.0,
    sd = 4.1
  ),
  incubation_period = list(
    mu = 5.95,
    sd = 4.31
  )
) %>%
  lapply(., function(x) {
    x$cv <- x$sd / x$mu
    x$shape_scale <- epitrix::gamma_mucv2shapescale(x$mu, x$cv)
    x$dist <- distcrete::distcrete(
      "gamma",
      shape = x$shape_scale$shape,
      scale = x$shape_scale$scale,
      w = 0.5,
      interval = 1
    )
    return(x)
  })


# Prior -------------------------------------------------------------------
prior_list <- list(pi = list(
  lower_bound_prior = 0.55,
  fn = function(x) {
    ifelse(x$pi > 0.55, log(1 / (1 - 0.55)), log(0))
  }
))


# Config ------------------------------------------------------------------
config_list <- list(
  n_iter = 10000,
  sample_every = 50,
  max_kappa = 3, # missed no cases
  init_kappa = 1, # number of generations before the last sampled ancestor
  move_kappa = TRUE,
  init_pi = 0.6, # 60% reporting rate
  move_pi = TRUE,
  init_tree = "star",
  outlier_threshold = 2
)

# Data -------------------------------------------------------------------
data_list <- list(
  dna = read.dna(here("data", "sequences.fasta"), skip = 0, format = "fasta"),
  contacts = read.csv(here("data", "contacts.csv")),
  linelist = read.csv(here("data", "linelist.csv"))
)

# Reorder DNA sequences
matching_indices <-
  match(data_list$linelist$case_id, rownames(data_list$dna))
data_list$dna <- data_list$dna[matching_indices, ]
# Date as integer
data_list$linelist$onset_inferred <-
  as.Date(data_list$linelist$onset_inferred, format = "%Y-%m-%d")
data_list$linelist$date_integer <-
  as.integer(difftime(
    data_list$linelist$onset_inferred,
    min(data_list$linelist$onset_inferred),
    units = "days"
  ))
data_list$linelist$group <- ifelse(
  grepl("^C", data_list$linelist$case_id),
  "patient",
  "hcw"
)


# outbreaker format ----------------------------------------------------------------
o2_data <- outbreaker_data(
  dates = data_list$linelist$onset_inferred,
  dna = data_list$dna,
  ids = data_list$linelist$case_id,
  w_dens = distribution_list$incubation_period$dist$d(1:50),
  f_dens = distribution_list$serial_interval$dist$d(1:50)
)


############################################################################
# TREE CUTTING
############################################################################
# In the below we 1st run the tree on the entire dataset and then we cut the tree at different dates
cutoff_dates <- unique(c(seq(min(as.Date(linelist$onset_inferred)),
                             max(as.Date(linelist$onset_inferred)),
                             by = 7),
                         max(linelist$onset_inferred)))
burnin <- 500

# Whole Run ---------------------------------------------------------
# set.seed(123)
# out <- outbreaker(
#   data = o2_data,
#   config = config_list,
#   priors = custom_priors(
#     pi = prior_list$pi$fn
#   )
# )
# out <- out %>%
#   filter(step > burnin) %>%
#   identify(ids = o2_data$ids)
#
# trees <- get_trees(out = out,
#                    ids = data_list$linelist$case_id,
#                    date = data_list$linelist$onset_inferred)
# cut_trees <- lapply(cutoff_dates, function(cutoff_date) {
#   lapply(trees, function(tree) {
#     tree_cut <- cut_tree_by_date(tree, cutoff_date)
#     if (nrow(tree_cut) == 0) {
#       return(NULL)
#     }
#     return(tree_cut)
#   })
# })
#
# saveRDS(cut_trees, here("data", "cut_trees.rds"))

# Reconstruction ---------------------------------------------------------
# here we reconstruct the tree using the o2_data up to the cutoff date
# set.seed(123)
# reconstructed_cut_trees <- lapply(cutoff_dates, function(cutoff_date) {
#   cat(which(cutoff_dates == cutoff_date), "out of", length(cutoff_dates), "\n")
#
#
#   filtered_indices <- data_list$linelist$onset_inferred <= cutoff_date
#   if(sum(filtered_indices) < 2) {
#     return(NULL)
#   }
#
#   o2_data_cut <- outbreaker2::outbreaker_data(
#     dates = data_list$linelist$onset_inferred[filtered_indices],
#     ids = data_list$linelist$case_id[filtered_indices],
#     dna = data_list$dna[filtered_indices, ],
#     w_dens = distribution_list$incubation_period$dist$d(1:50),
#     f_dens = distribution_list$serial_interval$dist$d(1:50)
#   )
#   out_cut <-  outbreaker2::outbreaker(
#     data = o2_data_cut,
#     config = config_list,
#     priors =  outbreaker2::custom_priors(
#       pi = prior_list$pi$fn
#     )
#   )
#   out_cut <- out_cut %>% filter(step > burnin)
#   out_cut <- identify(out_cut, ids = o2_data_cut$ids)
#   trees_cut <- get_trees(out_cut,
#                          ids = o2_data_cut$ids,
#                          date = o2_data_cut$dates)
#   return(trees_cut)
# })
# saveRDS(reconstructed_cut_trees, here("data", "reconstructed_cut_trees.rds"))


# Epicurve ----------------------------------------------------------------
linelist <- data_list$linelist
epicurve(cutoff_dates)

# Data --------------------------------------------------------------------
cut_trees <- readRDS(here("data", "cut_trees.rds"))
reconstructed_cut_trees <- readRDS(here("data", "reconstructed_cut_trees.rds")) # [-1]

# chi-square method -------------------------------------------------------
get_chi.test <- function(x, y) {
  if(is.null(x)) {
    x <- data.frame(from = character(0), to = character(0))
  }
  if(is.null(y)) {
    y <- data.frame(from = character(0), to = character(0))
  }
  tabx <- linktree:::ttable(
    from = x$from,
    to = x$to,
    levels = data_list$linelist$case_id
  ) %>%
    as.data.frame() %>%
    select(from, to, Freq) %>%
    rename(Freq_x = Freq)
  taby <- linktree:::ttable(
    from = y$from,
    to = y$to,
    levels = data_list$linelist$case_id
  ) %>%
    as.data.frame() %>%
    select(from, to, Freq) %>%
    rename(Freq_y = Freq)
  tab <- merge(tabx, taby, by = c("from", "to")) %>%
    filter(from != to) %>%
    filter(!(Freq_x == 0 & Freq_y == 0)) %>%
    select(Freq_x, Freq_y)

  tryCatch(
    {
      if (nrow(tab) <= 5) {
        test <- stats::fisher.test(tab)$p.value
      } else {
        test <- chisq.test(tab)$p.value
      }
      return(test)
    },
    error = function(e) {
      return(NA)
    }
  )
}
#@thibaut: ramener a une frequence d'ancestries
get_chi.test(x = bind_rows(reconstructed_cut_trees[[10]]),
             y = bind_rows(cut_trees[[10]])
)
# chidf <- furrr::future_map(2:length(cutoff_dates), function(i) {
#   data.frame(
#     cutoff_date = cutoff_dates[i],
#     p_value = get_chi.test(x = bind_rows(cut_trees[[i]]),
#                            y = bind_rows(reconstructed_cut_trees[[i]]))
#   )
# }) %>%
#   bind_rows() #!returns signficant p-values! @thibaut?

library(furrr)
plan(multisession, workers = 10)

#chisq at the sample level
chidf <- furrr::future_map(2:length(cutoff_dates), function(i) {
  lapply(1:190, function(j) {
    data.frame(
      cutoff_date = cutoff_dates[i],
      p_value = get_chi.test(x = cut_trees[[i]][[j]],
                             y = reconstructed_cut_trees[[i]][[j]]),
      sim = j
    )
  })
}) %>%
  bind_rows()


# compare any random posterior sample in each set 50 times
chidf <- furrr::future_map(2:length(cutoff_dates), function(i) {
  lapply(1:50, function(j) {
    sample <- sample(1:190, 2, replace = TRUE) #190 is length of the posterior samples
    data.frame(
      cutoff_date = cutoff_dates[i],
      p_value = get_chi.test(x = cut_trees[[i]][[sample[1]]],
                             y = reconstructed_cut_trees[[i]][[sample[2]]])
    )
  })
}) %>%
  bind_rows()


p_chi <- chidf %>%
  ggplot(aes(x = as.factor(cutoff_date), y = p_value)) +
  geom_violin() +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_x_discrete(
    limits = as.factor(cutoff_dates),
    breaks = as.factor(cutoff_dates),
    labels = format(cutoff_dates, "%d\n%b")) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "",
       y = "chi-square test p-value")

cowplot::plot_grid(p_chi, epicurve(), nrow = 2)



get_chi.test(x = bind_rows(reconstructed_cut_trees[[10]][[1]]),
             y = bind_rows(reconstructed_cut_trees[[10]][[10]])
)

# Compute Negative SIs ----------------------------------------------------
#
# six <- lapply(cut_trees[-1], \(x) bind_rows(x) %>%
#   mutate(
#     serial_interval = to_date - from_date,
#     serial_interval = as.integer(serial_interval)
#   ) %>%
#   select(serial_interval))
# siy <- lapply(reconstructed_cut_trees[-1], \(x) bind_rows(x) %>%
#   mutate(
#     serial_interval = to_date - from_date,
#     serial_interval = as.integer(serial_interval)
#   ) %>%
#   select(serial_interval))
#
# si <- lapply(1:length(six), function(i) {
#   bind_rows(
#     six[[i]] %>% mutate(set = "x"),
#     siy[[i]] %>% mutate(set = "y")
#   ) %>%
#     mutate(cut_date = cutoff_dates[i + 1])
# }) %>%
#   bind_rows()
#
# library(gghalves)
# p_si <- ggplot() +
#   aes(
#     x = cut_date,
#     y = serial_interval,
#     # fill = set,
#     col = set,
#     split = set,
#     group = interaction(set, cut_date)
#   ) +
#   gghalves::geom_half_boxplot(
#     data = si %>% filter(set == "x"),
#     side = "r",
#     outlier.shape = NA
#   ) +
#   gghalves::geom_half_boxplot(
#     data = si %>% filter(set == "y"),
#     side = "l",
#     outlier.shape = NA
#   ) +
#   # gghalves::geom_half_dotplot(
#   #   method = "histodot",
#   #   position = "identity",
#   #   stackratio = 0.5,
#   #   binwidth = 1,
#   #   dotsize = 0.01,
#   #   alpha = 0.5
#   # )+
#   # gghalves::geom_half_violin(aes(split = set),
#   #                            nudge = -1.2,
#   #                            side = c("l", "r")) +
#   geom_hline(yintercept = 0, lty = "dashed") +
#   scale_x_date(
#     breaks = cutoff_dates, date_labels = "%d\n%b",
#     limits = c(min(cutoff_dates), max(cutoff_dates) + 3)
#   ) +
#   coord_cartesian() +
#   theme_bw() +
#   theme(
#     legend.position = c(.01, .99),
#     legend.justification = c("left", "top"),
#     legend.box.just = "left",
#     legend.box.background = element_rect(colour = "black")
#   )
#
# cowplot::plot_grid(p_si, epicurve(), nrow = 2)


#
# si <- lapply(1:length(six), function(i){
#   x_si <- as.data.frame(table(six[[i]])) %>%
#     mutate(set = "x")
#   y_si <- as.data.frame(table(siy[[i]])) %>%
#     mutate(set = "y")
#   bind_rows(x_si, y_si) %>%
#     mutate(cut_date = cutoff_dates[i+1])
# }) %>%
#   bind_rows() %>%
#   group_by(cut_date, set) %>%
#   mutate(freq = Freq / sum(Freq),
#          si = as.numeric(as.character(Var1)))
#
# p_si <- si %>%
#   ggplot(aes(
#     x = si,
#     y = freq,
#     fill = set,
#     group = interaction(set, cut_date)
#   )) +
#   geom_col(position = "identity", alpha = 0.5) +
#   geom_vline(xintercept = 0, col = "black", linetype = "dashed") +
#   coord_flip()+
#   facet_wrap( ~ cut_date, scales = "free_x", nrow = 1)+
#   theme_bw()+
#   theme(axis.text.x = element_blank(),
#         axis.ticks.x = element_blank())+
#   labs(x = "",
#        y = "Serial Interval (days)")
# p_si











# Jaccard Similarity --------------------------------------------------------------
# pmin_table: Contains the minimum number of transmissions between pairs of individuals that both trees agree on.
# pmax_table: Contains the maximum number of transmissions between pairs of individuals considering all transmissions recorded in both trees.
# Jaccard similarity coefficient:
# jaccard <- function(A, B, levels = NULL) {
#   # Compute tables from dataframes
#   tabA <- linktree:::ttable(A$from, A$to, levels)
#   tabB <- linktree:::ttable(B$from, B$to, levels)
#   return(sum(pmin(tabA, tabB)) / sum(pmax(tabA, tabB)))
# }
#
# # Monte Carlo function
# mc_jaccard <- function(A, B, levels = NULL, n = 999) {
#   data <- rbind(A, B)
#   data$set <- c(rep("A", nrow(A)), rep("B", nrow(B)))
#
#   # reshuffle the set membership using sample()
#   jaccard_values <- replicate(n, {
#     data$set <- sample(data$set, replace = FALSE)
#     jaccard(
#       A = data[data$set == "A", ], B = data[data$set == "B", ],
#       levels = levels
#     )
#   })
#   return(jaccard_values)
# }
#
# library(furrr)
# plan(multisession, workers = length(listx))
# set.seed(123)
# null_jaccard <- future_map2(listx, listy, \(x, y) mc_jaccard(x, y, levels = data_list$linelist$case_id, n = 999),
#                             .options = furrr_options(seed = TRUE)
# )
# plan(sequential)
# est_jaccard <- lapply(1:length(listx), function(i) {
#   jaccard(listx[[i]], listy[[i]], levels = data_list$linelist$case_id)
# })
#
# # Proportion of null Jaccard values that are greater than or equal to the estimated Jaccard value.
# p_values <- lapply(1:length(est_jaccard), function(i) {
#   mean(null_jaccard[[i]] >= est_jaccard[[i]])
# }) %>%
#   unlist() %>%
#   data.frame(cutoff_date = cutoff_dates[-1], p_value = round(., 3))
#
#
# # make a dataframe with cutoff, type (est vs null) and jaccard
# jaccard_df <- bind_rows(
#   data.frame(cutoff_date = cutoff_dates[-1], jaccard = unlist(est_jaccard), type = "est"),
#   data.frame(cutoff_date = cutoff_dates[-1], jaccard = unlist(null_jaccard), type = "null")
# )
#
# p_jaccard <- ggplot() +
#   aes(
#     x = cutoff_date,
#     y = jaccard,
#     group = interaction(type, cutoff_date)
#   ) +
#   gghalves::geom_half_violin(
#     data = jaccard_df %>% filter(type == "null"),
#     aes(col = "null"), fill = "#cdcdcd", col = NA
#   ) +
#   gghalves::geom_half_boxplot(
#     data = jaccard_df %>% filter(type == "null"),
#     width = 0.5,
#     side = "r"
#   ) +
#   geom_point(
#     data = jaccard_df %>% filter(type == "est"),
#     aes(col = "est"), size = 3
#   ) +
#   geom_label(
#     data = p_values,
#     aes(
#       label = paste0("p = ", p_value),
#       group = cutoff_date,
#       x = cutoff_date,
#       y = 0.6
#     )
#   ) +
#   scale_x_date(
#     breaks = cutoff_dates, date_labels = "%d\n%b",
#     limits = c(min(cutoff_dates), max(cutoff_dates) + 3)
#   ) +
#   scale_color_manual(values = c("null" = "#3b3b3b", "est" = "red")) +
#   labs(
#     x = "",
#     y = "Jaccard Similarity"
#   ) +
#   theme_bw() +
#   theme(
#     legend.position = c(.01, .99),
#     legend.justification = c("left", "top"),
#     legend.box.just = "left",
#     legend.box.background = element_rect(colour = "black")
#   )
#
# # 70% of the counts of the elements in the two sets overlap
#
# cowplot::plot_grid(
#   plotlist = list(p_jaccard, p_epicurve),
#   ncol = 1,
#   rel_heights = c(0.5, 1),
#   labels = "AUTO"
# )
