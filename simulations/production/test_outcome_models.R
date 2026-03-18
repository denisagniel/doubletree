# Test outcome models specifically - these might be the slow ones

library(optimaltrees)

cat("\n=== TESTING OUTCOME MODEL SPEED (n=800) ===\n\n")

source("dgps/dgps_smooth.R")

# Generate data
set.seed(100000 + 4*10000 + 1)
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = 100000 + 4*10000 + 1)

cat(sprintf("Full data: n=%d\n", nrow(d$X)))
cat(sprintf("Treated: n=%d (%.1f%%)\n", sum(d$A), 100*mean(d$A)))
cat(sprintf("Control: n=%d (%.1f%%)\n", sum(1-d$A), 100*mean(1-d$A)))
cat("\n")

# Test propensity model
cat("1. Propensity model (full data, n=800)...\n")
t1 <- Sys.time()
prop_fit <- optimaltrees(
  X = d$X,
  y = d$A,
  loss_function = "log_loss",
  regularization = log(800) / 800,
  worker_limit = 4,
  verbose = FALSE
)
t2 <- Sys.time()
time_prop <- as.numeric(difftime(t2, t1, units = "secs"))
cat(sprintf("   Time: %.3fs\n\n", time_prop))

# Test outcome model on treated
cat(sprintf("2. Outcome model | treated (n=%d)...\n", sum(d$A)))
t1 <- Sys.time()
outcome_treated <- optimaltrees(
  X = d$X[d$A == 1, ],
  y = d$Y[d$A == 1],
  loss_function = "log_loss",
  regularization = log(sum(d$A)) / sum(d$A),
  worker_limit = 4,
  verbose = FALSE
)
t2 <- Sys.time()
time_treated <- as.numeric(difftime(t2, t1, units = "secs"))
cat(sprintf("   Time: %.3fs\n\n", time_treated))

# Test outcome model on control
cat(sprintf("3. Outcome model | control (n=%d)...\n", sum(1-d$A)))
t1 <- Sys.time()
outcome_control <- optimaltrees(
  X = d$X[d$A == 0, ],
  y = d$Y[d$A == 0],
  loss_function = "log_loss",
  regularization = log(sum(1-d$A)) / sum(1-d$A),
  worker_limit = 4,
  verbose = FALSE
)
t2 <- Sys.time()
time_control <- as.numeric(difftime(t2, t1, units = "secs"))
cat(sprintf("   Time: %.3fs\n\n", time_control))

# Simulate full cross-fitting
cat("4. Simulating full cross-fitting (5 folds, 15 trees total)...\n")
t1 <- Sys.time()

K <- 5
fold_ids <- sample(rep(1:K, length.out = nrow(d$X)))
total_trees <- 0

for (k in 1:K) {
  train_idx <- which(fold_ids != k)

  # Propensity
  prop <- optimaltrees(d$X[train_idx, ], d$A[train_idx],
                      loss_function = "log_loss",
                      regularization = log(length(train_idx)) / length(train_idx),
                      worker_limit = 4, verbose = FALSE)

  # Outcome | treated
  treated_idx <- train_idx[d$A[train_idx] == 1]
  out_trt <- optimaltrees(d$X[treated_idx, ], d$Y[treated_idx],
                          loss_function = "log_loss",
                          regularization = log(length(treated_idx)) / length(treated_idx),
                          worker_limit = 4, verbose = FALSE)

  # Outcome | control
  control_idx <- train_idx[d$A[train_idx] == 0]
  out_ctl <- optimaltrees(d$X[control_idx, ], d$Y[control_idx],
                          loss_function = "log_loss",
                          regularization = log(length(control_idx)) / length(control_idx),
                          worker_limit = 4, verbose = FALSE)

  total_trees <- total_trees + 3
}

t2 <- Sys.time()
time_full_cf <- as.numeric(difftime(t2, t1, units = "secs"))

cat(sprintf("   Total time: %.3fs for %d trees\n", time_full_cf, total_trees))
cat(sprintf("   Average per tree: %.3fs\n\n", time_full_cf / total_trees))

cat("=== SUMMARY ===\n")
cat(sprintf("Single propensity model: %.3fs\n", time_prop))
cat(sprintf("Single outcome|treated: %.3fs\n", time_treated))
cat(sprintf("Single outcome|control: %.3fs\n", time_control))
cat(sprintf("Full cross-fitting (15 trees): %.3fs\n", time_full_cf))
cat(sprintf("Average per tree in cross-fitting: %.3fs\n\n", time_full_cf / total_trees))

cat("Expected time for 500 DML estimates:\n")
cat(sprintf("  500 × %.3fs = %.0f seconds (%.1f minutes)\n",
            time_full_cf, 500 * time_full_cf, 500 * time_full_cf / 60))

cat("\n")
