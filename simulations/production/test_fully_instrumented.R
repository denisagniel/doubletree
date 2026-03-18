# Test with both dml_att and fit_nuisances_fold instrumented

library(optimaltrees)

cat("\n=== FULLY INSTRUMENTED TEST ===\n\n")

# Source dmltree functions
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")

# Source instrumented versions
source("fit_nuisances_fold_instrumented.R")  # Instrumented fit_nuisances_fold

# Override the original function with instrumented version
fit_nuisances_fold <- fit_nuisances_fold_instrumented

# Now source the rest (which will use our instrumented version)
source("dml_att_instrumented.R")
source("dgps/dgps_smooth.R")

# EXACT same problematic data
seed <- 100000 + 4*10000 + 1
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = seed)

cat(sprintf("Data: n=%d, treated=%d, control=%d\n\n",
            nrow(d$X), sum(d$A), sum(1-d$A)))

cat("Running fully instrumented estimate_att()...\n\n")
flush.console()

start_time <- Sys.time()

result <- tryCatch({
  fit <- dml_att_instrumented(
    X = d$X,
    A = d$A,
    Y = d$Y,
    K = 5,
    outcome_type = "binary",
    regularization = log(800) / 800,
    cv_regularization = FALSE,
    use_rashomon = FALSE,
    worker_limit = 4,
    verbose = FALSE
  )
  list(success = TRUE, theta = fit$theta, sigma = fit$sigma)
}, error = function(e) {
  cat(sprintf("\n[ERROR] %s\n", conditionMessage(e)))
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
