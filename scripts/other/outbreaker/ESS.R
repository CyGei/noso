source(here::here("scripts/outbreaker", "helpers.R"))

paper <- "eLife2022"

# Create output directory if it doesn't exist
output_path <- here::here("data", paper, "output/outbreaker")
dir.create(output_path, showWarnings = FALSE)
input <- loadbreaker(paper)
input <- input[[length(input)]]
set.seed(123)
#runs on 2k iterations with thinnig of 1
outbreaker_output <- outbreaker2::outbreaker(
  data = input$data,
  config = input$config,
  priors = input$priors)

outbreaker_output
plot(outbreaker_output, "post")
plot(outbreaker_output, "like")
plot(outbreaker_output, "mu")
plot(outbreaker_output, "lambda")
plot(outbreaker_output, "eps")

ess <- outbreaker_output %>%
  filter(step > 200) %>% # Burn-in
  select(post, like, pi, mu, lambda, eps) %>%
  coda::as.mcmc() %>%
  coda::effectiveSize()
#       pi        mu    lambda       eps
# 43.033798 11.976024  4.832339  4.166159
#1.  Required Total Iterations
get_total_iter <- function(desired_eff_size, post_burnin_samples, ess) {
  (post_burnin_samples * desired_eff_size) / ess
}
get_total_iter(desired_eff_size = 1000,
               post_burnin_samples = 2000 - 200,
               min(ess))
#432052.7

get_thinning_interval <- function(total_iter, desired_eff_size) {
  total_iter / desired_eff_size
}
get_thinning_interval(500000, 1000)
#500


# Conclusion
#total_iter = 500000 = 5e5
#thin = 500



#####################
# Convergence diagnostics
######################
library(coda)
# read all runs:
all_chains <- bind_rows(out, .id = "time_step") %>%
  mutate(time_step = factor(time_step, levels = seq_along(cutoff_dates)))

all_chains %>%
  ggplot(aes(x = step, y = post, colour = time_step)) +
  geom_line() +
  scale_colour_viridis_d()

params <- c("mu", "pi", "eps")
# Split into a list by chain (i.e., by time_step)
chain_list <- all_chains %>%
  select(time_step, step, all_of(params)) %>%
  group_split(time_step) %>%
  lapply(function(df) {
    mcmc(df %>% filter(step>500) %>% select(all_of(params)))
  })

# Convert to mcmc.list object
mcmc_list <- mcmc.list(chain_list)
gelman.diag(mcmc_list, autoburnin = FALSE)
