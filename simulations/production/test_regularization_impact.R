# Test if stronger regularization speeds up tree optimization at n=800

library(optimaltrees)

cat("\n=== TESTING REGULARIZATION IMPACT ON SPEED (n=800) ===\n\n")

source("dgps/dgps_smooth.R")

# Generate data
set.seed(100000 + 4*10000 + 1)
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = 100000 + 4*10000 + 1)

# Test different regularization values
reg_values <- c(
  log(800) / 800,        # Current (weak)
  0.01,                  # Moderate
  0.05,                  # Strong
  0.10,                  # Very strong
  0.20                   # Extremely strong
)

results <- data.frame()

for (reg in reg_values) {
  cat(sprintf("Testing regularization = %.4f... ", reg))
  flush.console()

  start_time <- Sys.time()

  fit <- optimaltrees(
    X = d$X,
    y = d$A,
    loss_function = "log_loss",
    regularization = reg,
    worker_limit = 4,
    verbose = FALSE
  )

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  cat(sprintf("%.3fs", elapsed))

  # Check tree complexity
  tree_text <- capture.output(print(fit$model))
  n_lines <- length(tree_text)

  cat(sprintf(" (tree lines: %d)\n", n_lines))
  flush.console()

  results <- rbind(results, data.frame(
    regularization = reg,
    time_sec = elapsed,
    tree_lines = n_lines
  ))
}

cat("\n=== RESULTS ===\n")
print(results, row.names = FALSE)

cat("\n=== ANALYSIS ===\n")
baseline_time <- results$time_sec[1]
cat(sprintf("Baseline (reg=%.4f): %.3fs\n", reg_values[1], baseline_time))

for (i in 2:nrow(results)) {
  speedup <- baseline_time / results$time_sec[i]
  cat(sprintf("reg=%.4f: %.3fs (%.1fx speedup)\n",
              results$regularization[i], results$time_sec[i], speedup))
}

# Estimate impact on full DML
cat("\n=== IMPACT ON FULL DML ===\n")
cat("Full DML requires 15 tree fits per estimate\n\n")

for (i in 1:nrow(results)) {
  dml_time <- 15 * results$time_sec[i]
  reps_per_hour <- 3600 / dml_time
  cat(sprintf("reg=%.4f: ~%.1fs per DML estimate, ~%.0f reps/hour\n",
              results$regularization[i], dml_time, reps_per_hour))
}

cat("\n=== RECOMMENDATION ===\n")
if (results$time_sec[nrow(results)] < 0.5 * baseline_time) {
  cat("✓ Higher regularization significantly speeds up optimization\n")
  cat(sprintf("  Recommend using reg >= %.2f for n=800\n", reg_values[which.min(results$time_sec)]))
} else {
  cat("⚠️  Regularization has limited impact on speed\n")
  cat("  Slowness is inherent to larger sample sizes\n")
}

cat("\n")
