pacman::p_load(o2groups, dplyr, here, GA)
obj_fn <- function(x) {
  size <- x[1:2]
  size <- round(size)

  o2groups.params <- list(
    duration = 100,
    group_n = 2,
    size = size,
    name = c("HCW", "patient"),
    gamma = c(12.7, 0.5),
    intro_n = c(1, 1),
    r0 = c(1.32, 1.61),
    generation_time = generation_time$d(1:100),
    incubation_period = incubation_period$r(1000)
  )

  output <- moVSdelta(o2groups.params)
  info <- output$out

  hcw_cases <- info$cases[info$level == "HCW"]
  patient_cases <- info$cases[info$level == "patient"]
  hcw_peak <- info$observed_peak[info$level == "HCW"]
  patient_peak <- info$observed_peak[info$level == "patient"]

  target_hcw_cases <- 105
  target_patient_cases <- 43
  target_hcw_peak <- 34
  target_patient_peak <- 21

  error1 <- (hcw_cases - target_hcw_cases)^2
  error2 <- (patient_cases - target_patient_cases)^2
  error3 <- (hcw_peak - target_hcw_peak)^2
  error4 <- (patient_peak - target_patient_peak)^2

  -sum(error1, error2, error3, error4)
}

# Set up the GA
ga_params <- ga(
  type = "real-valued",
  fitness = obj_fn,
  lower = c(105, 43),  # Lower bounds for gamma and size
  upper = c(1000, 1000),  # Upper bounds for gamma and size
  maxiter = 50,
  monitor = TRUE,
  run = 10
)

ga_params@solution

