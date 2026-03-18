# Debug: Why was verification script so slow?
# Run 5 reps of a single config with detailed logging

library(optimaltrees)
library(dplyr)

cat("\n=== DEBUGGING VERIFICATION SLOWNESS ===\n\n")

# Source doubletree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")

# Source DGPs
source("dgps/dgps_smooth.R")

# Source baseline methods
source("methods/method_forest_dml.R")
source("methods/method_linear_dml.R")

K_FOLDS <- 5
N_REPS <- 5

# Test problematic config: binary, n=800, forest (slowest combo)
cat("Testing: binary DGP, n=800, forest-DML\n")
cat("Reps:", N_REPS, "\n\n")

start_total <- Sys.time()

for (rep in 1:N_REPS) {
  cat(sprintf("[%d/%d] ", rep, N_REPS))

  start_rep <- Sys.time()

  # Data generation
  t1 <- Sys.time()
  d <- generate_dgp_binary_att(n = 800, tau = 0.10, seed = 90000 + rep)
  t2 <- Sys.time()
  data_time <- as.numeric(difftime(t2, t1, units = "secs"))

  # Model fit
  t1 <- Sys.time()
  fit <- att_forest(
    X = d$X, A = d$A, Y = d$Y,
    K = K_FOLDS,
    num.trees = 500,
    seed = 90000 + rep
  )
  t2 <- Sys.time()
  fit_time <- as.numeric(difftime(t2, t1, units = "secs"))

  rep_time <- as.numeric(difftime(Sys.time(), start_rep, units = "secs"))

  cat(sprintf("total=%.2fs (data=%.2fs, fit=%.2fs) theta=%.4f\n",
              rep_time, data_time, fit_time, fit$theta))

  # Check memory usage
  mem <- gc(verbose = FALSE)
  cat(sprintf("    Memory: %.1f MB used\n", sum(mem[,2])))
}

total_time <- as.numeric(difftime(Sys.time(), start_total, units = "secs"))

cat(sprintf("\nTotal time: %.2f sec (%.2f sec/rep)\n", total_time, total_time/N_REPS))
cat("\n=== MEMORY CHECK ===\n")
mem_final <- gc(verbose = FALSE)
print(mem_final)

cat("\n=== CONCLUSION ===\n")
if (total_time/N_REPS > 5) {
  cat("⚠️  SLOW: Average time > 5 sec/rep\n")
  cat("    Investigate: forest-DML with n=800 is bottleneck\n")
} else {
  cat("✓ NORMAL: Performance as expected\n")
}

cat("\n")
