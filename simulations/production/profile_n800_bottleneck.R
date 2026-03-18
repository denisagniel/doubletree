# Profile n=800 to find bottleneck
# Tree optimization is fast (0.04s), but something else is slow

library(optimaltrees)

cat("\n=== PROFILING n=800 BOTTLENECK ===\n\n")

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

cat("Data: n=800, p=4 features\n\n")

# Profile with detailed timing
cat("Starting profiled fit...\n")
flush.console()

start_total <- Sys.time()

# Create fold indices
K <- 5
set.seed(42)
fold_ids <- sample(rep(1:K, length.out = nrow(d$X)))

# Time each fold
fold_times <- numeric(K)

for (k in 1:K) {
  cat(sprintf("Fold %d/%d... ", k, K))
  flush.console()

  start_fold <- Sys.time()

  # Get train/test split
  test_idx <- which(fold_ids == k)
  train_idx <- which(fold_ids != k)

  X_train <- d$X[train_idx, ]
  A_train <- d$A[train_idx]
  Y_train <- d$Y[train_idx]
  X_test <- d$X[test_idx, ]

  # Fit propensity model
  t1 <- Sys.time()
  prop_fit <- optimaltrees(
    X = X_train,
    y = A_train,
    loss_function = "log_loss",
    regularization = log(length(A_train)) / length(A_train),
    worker_limit = 4,
    verbose = FALSE
  )
  t2 <- Sys.time()
  time_prop <- as.numeric(difftime(t2, t1, units = "secs"))

  # Fit outcome models
  t1 <- Sys.time()
  outcome_fit_treated <- optimaltrees(
    X = X_train[A_train == 1, ],
    y = Y_train[A_train == 1],
    loss_function = "log_loss",
    regularization = log(sum(A_train)) / sum(A_train),
    worker_limit = 4,
    verbose = FALSE
  )
  t2 <- Sys.time()
  time_outcome_treated <- as.numeric(difftime(t2, t1, units = "secs"))

  t1 <- Sys.time()
  outcome_fit_control <- optimaltrees(
    X = X_train[A_train == 0, ],
    y = Y_train[A_train == 0],
    loss_function = "log_loss",
    regularization = log(sum(1 - A_train)) / sum(1 - A_train),
    worker_limit = 4,
    verbose = FALSE
  )
  t2 <- Sys.time()
  time_outcome_control <- as.numeric(difftime(t2, t1, units = "secs"))

  fold_time <- as.numeric(difftime(Sys.time(), start_fold, units = "secs"))
  fold_times[k] <- fold_time

  cat(sprintf("%.2fs (prop=%.2fs, out_trt=%.2fs, out_ctl=%.2fs)\n",
              fold_time, time_prop, time_outcome_treated, time_outcome_control))
  flush.console()
}

total_time <- as.numeric(difftime(Sys.time(), start_total, units = "secs"))

cat(sprintf("\n=== RESULTS ===\n"))
cat(sprintf("Total time: %.2f seconds\n", total_time))
cat(sprintf("Average per fold: %.2f seconds\n", mean(fold_times)))
cat(sprintf("Fold times: %s\n", paste(sprintf("%.2f", fold_times), collapse=", ")))

cat("\n=== ANALYSIS ===\n")
if (mean(fold_times) > 10) {
  cat("⚠️  SLOW: > 10 seconds per fold\n")
  cat("Expected for full fit: ~50 seconds (5 folds)\n")
} else {
  cat("✓ NORMAL: < 10 seconds per fold\n")
  cat(sprintf("Expected for full fit: ~%.1f seconds\n", mean(fold_times) * 5))
}

cat("\n")
