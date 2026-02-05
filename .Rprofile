# .Rprofile - Auto-load treefarmr from local development
# This file is automatically sourced when R starts in this directory

# Try to load treefarmr from sibling directory
if (file.exists("../treefarmr")) {
  # Check if devtools is available
  if (requireNamespace("devtools", quietly = TRUE)) {
    tryCatch({
      devtools::load_all("../treefarmr")
      message("Loaded treefarmr from local development directory")
    }, error = function(e) {
      message("Note: Could not load treefarmr from ../treefarmr")
      message("Error: ", conditionMessage(e))
    })
  } else {
    message("Note: devtools not available, skipping treefarmr auto-load")
    message("Install devtools to auto-load treefarmr: install.packages('devtools')")
  }
} else {
  message("Note: treefarmr not found in ../treefarmr")
  message("If treefarmr is elsewhere, load it manually or update this .Rprofile")
}
