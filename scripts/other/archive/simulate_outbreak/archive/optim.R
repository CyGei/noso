source(here::here("scripts", "helpers.R"))
load_libraries()
pacman::p_load(o2groups, dplyr, here, distcrete)


# Fixed Params ------------------------------------------------------------

# Parameters
incubation_period <-
  distcrete::distcrete(
    "gamma",
    shape = epitrix::gamma_mucv2shapescale(5.95, (4.31 / 5.95))$shape,
    scale = epitrix::gamma_mucv2shapescale(5.95, (4.31 / 5.95))$scale,
    w = 0.5,
    interval = 1
  )

generation_time <-
  distcrete::distcrete(
    "gamma",
    shape = epitrix::gamma_mucv2shapescale(4.5, (4.1 / 3.0))$shape,
    scale =  epitrix::gamma_mucv2shapescale(4.5, (4.1 / 3.0))$scale,
    w = 0.5,
    interval = 1
  )
generation_time$d(1:100) %>% plot()
generation_time = c(0.05, 0.2, 0.35, 0.3, 0.1)/sum(c(0.05, 0.2, 0.35, 0.3, 0.1))
plot(generation_time, type = "h",  lwd = 20)

# Loss Function  ----------------------------------------------------------
loss_function <- function(o2groups.params) {
  sims <- lapply(1:50, function(i) {
    sim <- o2groups::simulate_groups(
      duration = o2groups.params$duration,
      group_n = o2groups.params$group_n,
      size = o2groups.params$size,
      name = o2groups.params$name,
      gamma = o2groups.params$gamma,
      intro_n = o2groups.params$intro_n,
      r0 = o2groups.params$r0,
      generation_time = o2groups.params$generation_time,
      incubation_period = o2groups.params$incubation_period
    )
    return(sim)
  })

  sims_df <- dplyr::bind_rows(sims, .id = "sim")

  # Targets
  target_cases <- c(105, 43)
  target_peak <- c(34, 21)

  # Average number of cases
  cases <- sims_df %>%
    dplyr::group_by(sim, group) %>%
    dplyr::summarise(case = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(mean = mean(case), .groups = "drop") %>%
    dplyr::mutate(target = target_cases)

  # Average peak date
  peaks <- sims_df %>%
    dplyr::group_by(sim, group, date_onset) %>%
    dplyr::summarise(case = n(), .groups = "drop") %>%
    dplyr::group_by(sim, group) %>%
    dplyr::arrange(sim, desc(case)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(mean = mean(date_onset), .groups = "drop") %>%
    dplyr::mutate(target = target_peak)

  # error <- c(sum((cases$mean - cases$target)^2),
  #            sum((peaks$mean - peaks$target)^2)
  #            )

  error <- sum(abs(cases$mean - cases$target))

  return(sum(error))
}

# Objective Function ------------------------------------------------------
obj_func <- function(params) {
  fHCW <- params[1]
  f <- setNames(c(fHCW, 1 - fHCW), c("HCW", "patient"))
  pop <-round(params[2])
  size <- round(as.vector(f * pop))


  generation_time = c(0.05, 0.2, 0.35, 0.3, 0.1)/sum(c(0.05, 0.2, 0.35, 0.3, 0.1))
  incubation_period <-
    distcrete::distcrete(
      "gamma",
      shape = epitrix::gamma_mucv2shapescale(5.95, (4.31 / 5.95))$shape,
      scale = epitrix::gamma_mucv2shapescale(5.95, (4.31 / 5.95))$scale,
      w = 0.5,
      interval = 1
    )

  o2groups.params <- list(
    duration = 100,
    group_n = 2,
    size = size,
    name = c("HCW", "patient"),
    intro_n = c(12, 4),
    r0 = c(1.32, 1.61),
    generation_time = generation_time,
    incubation_period = incubation_period$r(1000),
    gamma = linktree::delta2gamma(params[3:4])
  )

  loss_function(o2groups.params)
}

# Optimization ------------------------------------------------------------
pacman::p_load(optimParallel)
cl <- makeCluster(20)
setDefaultCluster(cl=cl)
clusterEvalQ(cl, {library(dplyr); library(magrittr)})
clusterExport(cl, c("obj_func", "loss_function"))

system.time({
  result <- optimParallel(
    par = c(0.2, 500, 0.5, 0.5),
    fn = obj_func,
    method = "L-BFGS-B",
    lower = c(0.1, 145, -0.9, -0.9),
    upper = c(0.4, 1000, 0.9, 0.9)
  )
})
stopCluster(cl)

result

# Grid style --------------------------------------------------------------

library(tidyverse)
f_grid <- seq(0.1, 0.9, by = 0.1)  # Range of f values to search
pop_grid <- seq(150, 1000, by = 200)  # Range of pop values to search

# Perform grid search
results <-
  expand.grid(f = f_grid, pop = pop_grid) %>%
  rowwise() %>%
  mutate(loss = loss_function(c(f, pop))) %>%
  arrange(loss)

# View results
print(results)
loss_function()
obj_func(c(0.5, 500))
