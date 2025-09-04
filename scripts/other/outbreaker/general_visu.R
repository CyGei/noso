
# Network -------------------------------------------------------------
library(epicontacts)
library(igraph)

load_helpers()
load_paper()
epicurve()

input <- loadbreaker(paper)
ctd = input[[length(input)]][["data"]][["ctd"]]

# contact_ids <- unique(c(ctd$contacts$from, ctd$contacts$to))
# linelist_ids <- ctd$linelist$id
# common_ids <- intersect(linelist_ids, contact_ids)
# ctd_trimmed <- subset(ctd, node_attribute = list("id" %in% common_ids))
g <- epicontacts:::as.igraph.epicontacts(ctd)
g <- as.undirected(g, mode = "collapse")
g <- delete_vertices(g, which(igraph::degree(g) == 0))

V(g)$color <- ifelse(V(g)$group == "patient", "purple", "orange")
set.seed(123)  # for reproducible layout
plot(
  g,
  layout = layout_nicely(g),
  vertex.label = NA,
  # remove labels for clarity
  vertex.size = 6,
  # node size
  edge.arrow.size = 0.3,
  # smaller arrow heads
  vertex.frame.color = "black"  # no border around nodes
)
group_colors <- c("patient" = "purple", "hcw" = "orange")

legend("topright",
       legend = names(group_colors),
       fill = group_colors,
       title = "Group",
       bty = "n",
       cex = 0.8)

proportions <- ctd$contact |>
  mutate(
    i = ifelse(str_starts(from, "C"), "patient", "hcw"),
    j = ifelse(str_starts(to, "C"), "patient", "hcw")
  ) %>%
  select(i, j) %>%
  summarise(
    # Count contacts where at least one patient is involved
    total_patient_involved = sum(i == "patient" | j == "patient"),
    # Count patient-patient contacts
    within_patient = sum(i == "patient" & j == "patient"),
    # Count contacts where at least one hcw is involved
    total_hcw_involved = sum(i == "hcw" | j == "hcw"),
    # Count hcw-hcw contacts
    within_hcw = sum(i == "hcw" & j == "hcw")
  ) |>
  mutate(
    # Calculate the proportion of within-group contacts for patients
    prop_patient_within = if_else(
      total_patient_involved > 0,
      within_patient / total_patient_involved,
      0
    ),
    # Calculate the proportion of within-group contacts for hcw
    prop_hcw_within = if_else(total_hcw_involved > 0, within_hcw / total_hcw_involved, 0)
  )
proportions



# Transmission tree ------------------------------------------------
#https://fontawesome.com/v4/icons/
load_helpers()
load_paper()
epicurve()
library(tidygraph)
library(ggraph)

x = out[[length(out)]]
x <- o2ools::identify(x, linelist$case_id)
df <- o2ools::get_trees(x,
                        t_inf = TRUE,
                        group = linelist$group,
                        onset = linelist$onset) %>%
  bind_rows(.id = "iteration")

# for each infectee, what is the most common infector across all 999 iterations
most_common_infectors <- df %>%
  group_by(to) %>%
  count(from, sort = TRUE) %>%
  # as frequency
  mutate(frequency = n / 999) %>%
  slice(1) %>%
  ungroup() %>%
  select(from, to, frequency)

epi <- make_epicontacts(linelist = linelist,
                        contacts = most_common_infectors,
                        directed = TRUE)
g <- epicontacts:::as.igraph.epicontacts(epi)

g_tidy <- as_tbl_graph(g)
g_tidy <- g_tidy %>%
  activate(edges) %>%
  mutate(from_group = .N()$group[from])

layout_df <- create_layout(g_tidy, layout = "stress") # Use stress as a base for initial positioning if needed
layout_df$x <- V(g)$onset # Set x-coordinate to onset date

ggraph(g_tidy, layout = layout_df) +
  geom_edge_link(
    aes(edge_colour = from_group, width = frequency, #label = round(frequency, 2)),
        # angle_calc = 'along',
        # label_dodge = unit(0.1, 'mm'),
        arrow = arrow(length = unit(2, 'mm')),
        end_cap = circle(1.5, 'mm'),
        alpha = 0.6,
        show.legend = c(edge_colour = FALSE, width = TRUE)
    ) +
      geom_node_point(aes(colour = group), size = 3) +
      scale_colour_manual(
        "",
        values = c("hcw" = "orange", "patient" = "purple"),
        breaks = c("hcw", "patient"),
        labels = c("HCW", "Patient")
      ) +
      scale_edge_colour_manual(
        values = c("hcw" = "orange", "patient" = "purple"),
        breaks = c("hcw", "patient")
      ) +
      scale_edge_width("Posterior support",
                       range = c(0.1, 3),
                       breaks = c(0.1, 0.2, 0.5)) + # Adjust the range to control min and max line thickness
      theme_graph() +
      theme(
        axis.text.x = element_text(
          size = 8,
          angle = 45,
          hjust = 1
        ),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_line(colour = "black"),
        axis.title.x = element_text(margin = margin(t = 10), colour = "black"),
        legend.position = "bottom",
      ) +
      xlab("Onset Date")



    ############
    # ids <- unique(ctd$linelist$id)
    # group_colours <- ifelse(grepl("^C", ids), "purple", "orange")
    # names(group_colours) <- ids
    #
    # vis_epicontacts(x = make_epicontacts(linelist = ctd$linelist,
    #                 contacts = most_common_infectors,
    #                 directed = TRUE),
    #                 node_shape = "group",
    #                 shapes = c("patient" = "bed", "hcw" = "user-md"),
    #                 node_colour = "id",  # use id as colouring variable
    #                 col_pal = scales::manual_pal(group_colours),
    #                 edge_colour = "from_group")
    # linelist <- ctd$linelist
    # linelist$onset_numeric <- as.numeric(as.Date(linelist$onset))
    # set.seed(123)  # for reproducibility
    # layout_mat <- data.frame(
    #   id = linelist$id,
    #   x = linelist$onset_numeric,
    #   y = runif(nrow(linelist), min = 0, max = 1)
    # )
    # layout_coords <- as.matrix(layout_mat[, c("x", "y")])
    # rownames(layout_coords) <- layout_mat$id
    # epi <- make_epicontacts(
    #   linelist = linelist,
    #   contacts = most_common_infectors,
    #   directed = TRUE
    # )
    #
    # vis_epicontacts(
    #   x = epi,
    #   layout = layout_coords,
    #   node_colour = "id",
    #   col_pal = scales::manual_pal(group_colours),
    #   node_shape = "group",
    #   shapes = c("patient" = "bed", "hcw" = "user-md"),
    #   edge_colour = "from_group"
    # )
    #
