source(
    here::here("scripts", "helpers.R")
)
source(
    here::here("scripts/other/truncated_tree", "helpers.R")
)
loadbreaker()

# Data --------------------------------------------------------------------
trunc_trees <- readRDS(here("data/other/truncated_tree", "trunc_trees.rds"))
realtime_trees <- readRDS(here("data/other/truncated_tree", "realtime_trees.rds"))

plan(multisession, workers = availableCores() - 1)
chidf <- furrr::future_map(2:length(cutoff_dates), function(i) {
  data.frame(
    cutoff_date = cutoff_dates[i],
    p_value = get_chisq(
      x = bind_rows(realtime_trees[[i]]),
      y = bind_rows(trunc_trees[[i]]),
      ids = linelist$case_id
    )
  )
},
.options = furrr_options(seed = TRUE)) %>%
  bind_rows()


p_chi <- chidf %>%
  add_row(cutoff_date = cutoff_dates[1], p_value = NA) %>%
  ggplot(aes(x = as.Date(cutoff_date), y = p_value)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "",
       y =  paste0("\u03C7\u00B2"," p-value"))+
  theme_noso(fill=F, col=F, date = T)

cowplot::plot_grid(p_chi, epicurve(), nrow = 2)












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
