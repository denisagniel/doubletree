# .Rprofile - Auto-load optimaltrees from local development
# This file is automatically sourced when R starts in this directory

# Try to load optimaltrees from sibling directory
if (file.exists("../optimaltrees")) {
  # Check if devtools is available
  if (requireNamespace("devtools", quietly = TRUE)) {
    tryCatch({
      devtools::load_all("../optimaltrees")
      message("Loaded optimaltrees from local development directory")
    }, error = function(e) {
      message("Note: Could not load optimaltrees from ../optimaltrees")
      message("Error: ", conditionMessage(e))
    })
  } else {
    message("Note: devtools not available, skipping optimaltrees auto-load")
    message("Install devtools to auto-load optimaltrees: install.packages('devtools')")
  }
} else {
  message("Note: optimaltrees not found in ../optimaltrees")
  message("If optimaltrees is elsewhere, load it manually or update this .Rprofile")
}
