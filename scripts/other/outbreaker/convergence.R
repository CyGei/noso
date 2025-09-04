# convergence

# This will run outbreaker on 10 chains for each paper
# and save the results in data/<paper>/output/convergence

source(here::here("scripts/other/outbreaker", "helpers.R"))
papers <- c("JHI2021", "eLife2022")

# --- Pre-load Data ---
inputs <- list()
for (paper in papers) {
  output_path <- here::here("data", paper, "output/convergence")
  dir.create(output_path, showWarnings = FALSE)

  input_list <- loadbreaker(paper)
  inputs[[paper]] <- input_list[[length(input_list)]]
}
# --- End Data Loading ---


# --- Run Chains in Parallel ---
future::plan(future::multisession, workers = 21)
params <- expand.grid(paper = papers, chain_id = 1:10)

furrr::future_walk(1:nrow(params), function(row_idx) {
  paper <- params$paper[row_idx]
  chain_id <- params$chain_id[row_idx]
  input <- inputs[[paper]]

  out <- outbreaker2::outbreaker(
    data = input$data,
    config = input$config,
    priors = input$priors
  )

  saveRDS(out,
          file = here::here(
            "data",
            paper,
            "output/convergence",
            paste0("outbreaker_chain_", chain_id, ".rds")
          ))

}, .options = furrr::furrr_options(seed = TRUE))
# --- End Parallel Processing ---

future::plan(future::sequential)

# --- Combine Chains ---
all_chains <- lapply(papers, function(paper) {
  output_path <- here::here("data", paper, "output/convergence")
  files <- list.files(output_path, pattern = "outbreaker_chain_\\d+\\.rds", full.names = TRUE)
  chains <- lapply(files, readRDS)
  chains
})
names(all_chains) <- papers
# --- End Combine Chains ---

# --- Convergence Diagnostics ---

all_chains$eLife2022
trace <- function(x, paper, burnin = 500) {
  data <- bind_rows(all_chains[[paper]][2:10], .id = "chain")
  data <- data %>%
    filter(step > burnin) %>%
    mutate(chain = as.factor(chain))
  ggplot(data = data,
         aes(x = step, y = !!sym(x),
             col = chain, group = chain)) +
    geom_line(alpha = 0.5, linewidth = 0.5)+
    scale_colour_viridis_d("Chain ID")+
    theme(legend.position="bottom",
          legend.spacing.x = unit(0.25, 'cm'),
          legend.key.width = unit(1.5, "cm"),
          legend.title = element_text(hjust = 0.5),
          legend.justification = "center")
}

grid <- function(paper){
  cowplot::plot_grid(
    trace("post",paper) + labs(x = "") + theme(legend.position = "none"),
    trace("like", paper) + labs(x = "") + theme(legend.position = "none"),
    trace("mu", paper) + theme(legend.position = "none"),
    trace("pi", paper) + theme(legend.position = "none"),
    trace("eps", paper) + theme(legend.position = "none"),
    trace("lambda", paper) + theme(legend.position = "none"),
    nrow = 3,
    align = "v"
  )
}

grid("eLife2022")
grid("JHI2021")

library(coda)
rhat <- function(paper, params = c("mu", "pi", "eps", "lambda")) {
  chains <- all_chains[[paper]][-1]
  mcmc_list <- lapply(chains, function(chain) {
    chain %>%
      filter(step>500) %>%
      select(all_of(params)) %>%
      as.matrix() %>%
      mcmc()
  }) %>%
    mcmc.list()
  gelman.diag(mcmc_list, autoburnin = FALSE)

}
rhat("eLife2022")
rhat("JHI2021")

# --- End Convergence Diagnostics ---
