# Debug Config 4: Why is continuous, n=800, tree-DML hanging?

library(optimaltrees)

cat("\n=== DEBUGGING CONFIG 4 HANG ===\n\n")

# Source dmltree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")

# Reproduce exact config that's hanging
cat("Generating continuous DGP, n=800...\n")
flush.console()

set.seed(100000 + 4*10000 + 1)  # Same seed as Config 4, Rep 1
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = 100000 + 4*10000 + 1)

cat(sprintf("Data generated: %d rows, %d columns\n", nrow(d$X), ncol(d$X)))
cat(sprintf("Treatment: %d treated, %d control\n", sum(d$A), sum(1-d$A)))
cat(sprintf("Outcome range: [%.3f, %.3f]\n", min(d$Y), max(d$Y)))
flush.console()

# Try fitting with verbose output
cat("\nAttempting DML-ATT fit with verbose=TRUE...\n")
flush.console()

start_time <- Sys.time()

# Set a timeout alarm
result <- tryCatch({
  fit <- estimate_att(
    X = d$X,
    A = d$A,
    Y = d$Y,
    K = 5,
    outcome_type = "binary",
    regularization = log(800) / 800,
    cv_regularization = FALSE,
    use_rashomon = FALSE,
    worker_limit = 4,
    verbose = TRUE  # Enable verbose output
  )
  list(success = TRUE, theta = fit$theta, error = NA)
}, error = function(e) {
  list(success = FALSE, theta = NA, error = conditionMessage(e))
})

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\n=== RESULT ===\n"))
cat(sprintf("Time elapsed: %.1f seconds\n", elapsed))

if (result$success) {
  cat(sprintf("SUCCESS: theta = %.4f\n", result$theta))
} else {
  cat(sprintf("FAILED: %s\n", result$error))
}

cat("\n=== DONE ===\n")
