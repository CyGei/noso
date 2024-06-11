# This script will run outbreaker2 on the truncated linelist data
# at every time point for real-time analysis
source(here::here("scripts/outbreaker", "helpers.R"))

# Run outbreaker2 for both papers
papers <- c("JHI2021", "eLife2022")

for (paper in papers) {
  # Create output directory if it doesn't exist
  output_path <- here::here("data", paper, "output/outbreaker")
  dir.create(output_path, showWarnings = FALSE)

  input <- loadbreaker(paper)
  workers <- future::availableCores() - 1
  future::plan(future::multisession, workers = workers)

  f_null <- function(data, param) {
    return(0.0)
  }

  set.seed(123)

  furrr::future_walk(input, function(x) {
    outbreaker_output <- outbreaker2::outbreaker(
      data = x$data,
      config = x$config,
      priors = x$priors,
      likelihoods = custom_likelihoods(
        contact = f_null,
        timing_sampling = f_null
      )
    )

    # Save output for the current cutoff date
    cutoff_date <- unique(max(x$data$dates))
    saveRDS(outbreaker_output, here::here(output_path, paste0("outbreaker_", cutoff_date, ".rds")))
  }, .options = furrr::furrr_options(seed = TRUE))
}


# Trace plots for the outbreaker2 results
paper <- papers[2]
input_path <- here::here("data", paper, "output/outbreaker")
cutoff_dates <- loadcutoff_dates(paper)
results <- furrr::future_map(
  .x = list.files(input_path, full.names = TRUE),
  ~ readRDS(.x) %>%
    mutate(cutoff = as.integer(gsub(".*outbreaker_(.*).rds", "\\1", .x)),
           cutoff_date = cutoff_dates[cutoff+1])
)
pacman::p_load(ggplot2)
trace <- function(x) {
  ggplot(data = bind_rows(results),
         aes(x = step, y = !!sym(x),
             col = cutoff_date, group = cutoff_date)) +
    geom_line(alpha = 0.5, linewidth = 0.5)+
    scale_colour_viridis_c("Date Cutoff", trans = "date")+
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    theme(legend.position="bottom",
          legend.spacing.x = unit(0.25, 'cm'),
          legend.key.width = unit(1.5, "cm"),
          legend.title = element_text(hjust = 0.5),
          legend.justification = "center")
}
grid <- cowplot::plot_grid(
  trace("post") + labs(x = "") + theme(legend.position = "none"),
  trace("like") + labs(x = "") + theme(legend.position = "none"),
  trace("mu") + theme(legend.position = "none"),
  trace("pi") + theme(legend.position = "none"),
  nrow = 2,
  align = "v"
)
cowplot::plot_grid(
  grid,
  cowplot::get_legend(trace("post")),
  nrow = 2,
  rel_heights = c(1, 0.2)
)
