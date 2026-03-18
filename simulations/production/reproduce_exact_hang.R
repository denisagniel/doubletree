# Reproduce EXACT Config 4, Rep 1 that's hanging

library(optimaltrees)

cat("\n=== REPRODUCING EXACT HANG: Config 4, Rep 1 ===\n\n")

# Source everything exactly as production run does
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")

# EXACT same config and seed as production run
config_i <- 4
rep <- 1
SEED_OFFSET <- 100000
K_FOLDS <- 5

seed <- SEED_OFFSET + config_i*10000 + rep

cat(sprintf("Config: 4 (continuous, n=800, tree-DML)\n"))
cat(sprintf("Rep: %d\n", rep))
cat(sprintf("Seed: %d\n\n", seed))

# Generate data with EXACT same seed
cat("Generating data...\n")
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = seed)

cat(sprintf("Data: n=%d, treated=%d, control=%d\n",
            nrow(d$X), sum(d$A), sum(1-d$A)))
cat(sprintf("Y range: [%.3f, %.3f]\n\n", min(d$Y), max(d$Y)))

# Call dml_att with EXACT same parameters as production run
cat("Calling estimate_att() with production parameters...\n")
cat("(will timeout after 15 seconds if hanging)\n\n")
flush.console()

start_time <- Sys.time()

# Run in separate process with timeout
result <- tryCatch({
  fit <- estimate_att(
    X = d$X,
    A = d$A,
    Y = d$Y,
    K = K_FOLDS,
    outcome_type = "binary",
    regularization = log(800) / 800,
    cv_regularization = FALSE,
    use_rashomon = FALSE,
    worker_limit = 4,
    verbose = TRUE  # Enable to see where it hangs
  )
  list(success = TRUE, theta = fit$theta, sigma = fit$sigma)
}, error = function(e) {
  list(success = FALSE, error = conditionMessage(e))
})

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\n=== RESULT ===\n"))
cat(sprintf("Time: %.1f seconds\n", elapsed))

if (result$success) {
  cat(sprintf("SUCCESS: theta = %.4f, sigma = %.4f\n", result$theta, result$sigma))
} else {
  cat(sprintf("FAILED: %s\n", result$error))
}

cat("\n")
