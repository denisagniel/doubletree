# Test: Does memory accumulate when running many configs?

library(optimaltrees)
library(dplyr)

cat("\n=== MEMORY ACCUMULATION TEST ===\n\n")

# Source doubletree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")
source("methods/method_forest_dml.R")
source("methods/method_linear_dml.R")

# Test 20 consecutive runs
N_RUNS <- 20
K_FOLDS <- 5

cat("Running", N_RUNS, "consecutive forest-DML fits (n=400)\n\n")

memory_usage <- numeric(N_RUNS)

for (i in 1:N_RUNS) {
  # Generate data
  d <- generate_dgp_binary_att(n = 400, tau = 0.10, seed = 90000 + i)

  # Fit
  fit <- att_forest(
    X = d$X, A = d$A, Y = d$Y,
    K = K_FOLDS,
    num.trees = 500,
    seed = 90000 + i
  )

  # Check memory
  mem <- gc(verbose = FALSE)
  memory_usage[i] <- sum(mem[,2])

  cat(sprintf("[%d/%d] theta=%.4f, memory=%.1f MB\n", i, N_RUNS, fit$theta, memory_usage[i]))
}

cat("\n=== MEMORY ANALYSIS ===\n\n")
cat(sprintf("Initial memory: %.1f MB\n", memory_usage[1]))
cat(sprintf("Final memory: %.1f MB\n", memory_usage[N_RUNS]))
cat(sprintf("Increase: %.1f MB (%.1f%%)\n",
            memory_usage[N_RUNS] - memory_usage[1],
            100 * (memory_usage[N_RUNS] - memory_usage[1]) / memory_usage[1]))

# Check for linear growth
if (length(memory_usage) > 5) {
  lm_fit <- lm(memory_usage ~ seq_along(memory_usage))
  slope <- coef(lm_fit)[2]
  cat(sprintf("\nMemory growth rate: %.3f MB/run\n", slope))

  if (slope > 1) {
    cat("⚠️  MEMORY LEAK: Memory grows >1 MB per run\n")
    cat("    Estimated memory at 240 runs: %.1f MB\n", memory_usage[1] + slope*240)
  } else if (slope > 0.1) {
    cat("⚠️  SLOW LEAK: Memory grows slightly each run\n")
  } else {
    cat("✓ NO LEAK: Memory stable across runs\n")
  }
}

cat("\n")
