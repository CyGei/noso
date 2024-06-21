# if (interactive()) {
#   message("Loading project-specific settings...")
  
#   # Source the helpers.R file
#   source(here::here("scripts", "helpers.R"))
  
#   # Run load_helpers()
#   load_helpers()
  
#   # Set up an event handler for when a new file is opened, if RStudio is available
#   if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
#     setHook("rstudio.sessionInit", function(newSession) {
#       if (!newSession) {
#         load_helpers()
#       }
#     }, action = "append")
#   }
# }


if (!requireNamespace("here", quietly = TRUE)) {
  install.packages("here")
}
source(here::here("scripts", "helpers.R"))
