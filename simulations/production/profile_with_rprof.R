# Use Rprof to profile actual estimate_att() call

library(optimaltrees)

cat("\n=== PROFILING estimate_att() WITH Rprof (n=800) ===\n\n")

# Source dmltree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")

# Generate data
set.seed(100000 + 4*10000 + 1)
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = 100000 + 4*10000 + 1)

cat(sprintf("Data: n=%d, p=%d features\n\n", nrow(d$X), ncol(d$X)))

cat("Starting profiled fit (this may take a while if there's a bottleneck)...\n")
flush.console()

# Profile with fine-grained timing
Rprof("dml_att_profile.out", interval = 0.01, memory.profiling = TRUE)

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
    verbose = FALSE
  )
  list(success = TRUE, theta = fit$theta, sigma = fit$sigma)
}, error = function(e) {
  list(success = FALSE, error = conditionMessage(e))
})

Rprof(NULL)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

cat(sprintf("\n=== RESULTS ===\n"))
cat(sprintf("Elapsed time: %.2f seconds\n\n", elapsed))

if (result$success) {
  cat(sprintf("SUCCESS: theta = %.4f, sigma = %.4f\n\n", result$theta, result$sigma))
} else {
  cat(sprintf("FAILED: %s\n\n", result$error))
}

# Analyze profile
cat("=== PROFILE ANALYSIS ===\n")
profile <- summaryRprof("dml_att_profile.out", memory = "both")

cat("\nTop 10 functions by total time:\n")
print(head(profile$by.total, 10))

cat("\nTop 10 functions by self time:\n")
print(head(profile$by.self, 10))

cat("\nMemory usage:\n")
print(profile$by.total[1:5, c("total.time", "mem.total")])

# Clean up
unlink("dml_att_profile.out")

cat("\n=== DONE ===\n")
