# The purpose of this script is to explore the relationship between gamma and mo's ratio.
# We will simulate multiple outbreaks with different gamma values and calculate the mo's ratio for each outbreak.
# We will then plot the relationship between gamma and mo's ratio.

#devtools::install_github("CyGei/o2groups")
library(o2groups)
duration = 365
group_n = 2
size = c(500, 500)
name = c("HCW", "patient")
intro_n = c(1, 3)
r0 = c(2, 2)
generation_time = c(0, 0.1, 0.2, 0.4, 0.2, 0.1, 0)
incubation_period = sample(1:14, sum(size), replace = TRUE)


# Generate a list of scenarios with different gamma values
scenarios <- lapply(1:100, function(i) {
  list(
    duration = duration,
    group_n = group_n,
    size = size,
    name = name,
    gamma = runif(2, 0, 10),
    intro_n = intro_n,
    r0 = r0,
    generation_time = generation_time,
    incubation_period = incubation_period
  )
})

# Run simulate_groups for each scenario
set.seed(123)
timing <- system.time(
  outbreaks <- lapply(scenarios, function(scenario) {
    o2groups::simulate_groups_furrr(
      sim_n = 100,
      duration = scenario$duration,
      group_n = scenario$group_n,
      size = scenario$size,
      name = scenario$name,
      gamma = scenario$gamma,
      intro_n = scenario$intro_n,
      r0 = scenario$r0,
      generation_time = scenario$generation_time,
      incubation_period = scenario$incubation_period
    )
  })
)
print(timing)
