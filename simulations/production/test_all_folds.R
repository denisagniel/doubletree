# Test all 5 folds to find which one hangs

library(optimaltrees)

cat("\n=== TESTING ALL 5 FOLDS FOR Config 4, Rep 1 ===\n\n")

source("dgps/dgps_smooth.R")

# Exact same data as production run
seed <- 100000 + 4*10000 + 1
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = seed)

K <- 5
set.seed(42)  # doubletree uses seed 42 for folds
fold_ids <- sample(rep(1:K, length.out = nrow(d$X)))

cat(sprintf("Data: n=%d, K=%d folds\n\n", nrow(d$X), K))

for (k in 1:K) {
  cat(sprintf("=== Fold %d/%d ===\n", k, K))

  train_idx <- which(fold_ids != k)
  X_train <- d$X[train_idx, ]
  y_train <- d$A[train_idx]

  cat(sprintf("  Train n=%d, Test n=%d\n", length(train_idx), sum(fold_ids == k)))

  # Try propensity model
  cat("  Fitting propensity... ")
  flush.console()
  t1 <- Sys.time()
  prop_fit <- optimaltrees(X_train, y_train,
                          loss_function = "log_loss",
                          regularization = log(length(y_train)) / length(y_train),
                          worker_limit = 4, verbose = FALSE)
  t2 <- Sys.time()
  time_prop <- as.numeric(difftime(t2, t1, units = "secs"))
  cat(sprintf("%.3fs\n", time_prop))

  # Try outcome | control
  control_idx <- train_idx[d$A[train_idx] == 0]
  X_control <- d$X[control_idx, ]
  y_control <- d$Y[control_idx]

  cat(sprintf("  Fitting outcome|control (n=%d)... ", length(control_idx)))
  flush.console()
  t1 <- Sys.time()
  out_fit <- optimaltrees(X_control, y_control,
                         loss_function = "log_loss",
                         regularization = log(length(control_idx)) / length(control_idx),
                         worker_limit = 4, verbose = FALSE)
  t2 <- Sys.time()
  time_out <- as.numeric(difftime(t2, t1, units = "secs"))
  cat(sprintf("%.3fs\n", time_out))

  cat(sprintf("  Fold %d total: %.3fs\n\n", k, time_prop + time_out))
  flush.console()
}

cat("=== ALL FOLDS COMPLETED ===\n")
cat("No hang detected!\n\n")
