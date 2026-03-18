# Test COMPLETE estimate_att() call - the EXACT call from production run

library(optimaltrees)

cat("\n=== TESTING COMPLETE estimate_att() - Config 4, Rep 1 ===\n\n")

# Source dmltree exactly as production run does
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")

# EXACT same config and seed
seed <- 100000 + 4*10000 + 1
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = seed)

cat(sprintf("Data: n=%d, treated=%d, control=%d\n\n",
            nrow(d$X), sum(d$A), sum(1-d$A)))

# Call dml_att with EXACT production parameters
cat("Calling estimate_att() with production parameters...\n")
flush.console()

start_time <- Sys.time()

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
    verbose = FALSE  # Disable verbose to avoid C++ spam
  )
  list(success = TRUE, theta = fit$theta, sigma = fit$sigma)
}, error = function(e) {
  cat(sprintf("\nERROR: %s\n", conditionMessage(e)))
  list(success = FALSE, error = conditionMessage(e))
})

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\n=== RESULT ===\n"))
cat(sprintf("Elapsed: %.2f seconds\n", elapsed))

if (result$success) {
  cat(sprintf("SUCCESS: theta = %.4f, sigma = %.4f\n", result$theta, result$sigma))
} else {
  cat(sprintf("FAILED: %s\n", result$error))
}

cat("\n")
