# Complete profiling of estimate_att() to find where the hang occurs

library(optimaltrees)

cat("\n=== COMPLETE DML-ATT PROFILING (n=800) ===\n\n")

# Source dmltree with instrumentation
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")

# Instrument nuisance_trees.R by adding timing
cat("Reading nuisance_trees.R to add timing instrumentation...\n")
nuisance_code <- readLines("../../R/nuisance_trees.R")

# Add timing checkpoint after each major step
instrumented_file <- "../../R/nuisance_trees_instrumented.R"
writeLines(nuisance_code, instrumented_file)
source(instrumented_file)

source("../../R/estimate_att.R")
source("dgps/dgps_smooth.R")

# Generate data
set.seed(100000 + 4*10000 + 1)
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = 100000 + 4*10000 + 1)

cat(sprintf("Data: n=%d, p=%d features\n\n", nrow(d$X), ncol(d$X)))

# Wrap dml_att with detailed timing
cat("Starting instrumented estimate_att()...\n")
flush.console()

# Time checkpoints
times <- list()
times$start <- Sys.time()

# Call with tryCatch to capture any hangs
result <- tryCatch({

  # We'll manually call parts of dml_att to time each step
  K <- 5
  n <- nrow(d$X)

  # 1. Create folds
  cat("Step 1: Creating folds... ")
  flush.console()
  t1 <- Sys.time()
  fold_ids <- sample(rep(1:K, length.out = n))
  times$fold_creation <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("%.3fs\n", times$fold_creation))
  flush.console()

  # 2. Fit nuisance functions (this should be fast based on previous test)
  cat("Step 2: Fitting nuisance functions (5 folds)... ")
  flush.console()
  t1 <- Sys.time()

  nuisances <- list()
  for (k in 1:K) {
    test_idx <- which(fold_ids == k)
    train_idx <- which(fold_ids != k)

    # Fit trees
    prop_fit <- optimaltrees(
      X = d$X[train_idx, ],
      y = d$A[train_idx],
      loss_function = "log_loss",
      regularization = log(length(train_idx)) / length(train_idx),
      worker_limit = 4,
      verbose = FALSE
    )

    outcome_treated <- optimaltrees(
      X = d$X[train_idx, ][d$A[train_idx] == 1, ],
      y = d$Y[train_idx][d$A[train_idx] == 1],
      loss_function = "log_loss",
      regularization = log(sum(d$A[train_idx])) / sum(d$A[train_idx]),
      worker_limit = 4,
      verbose = FALSE
    )

    outcome_control <- optimaltrees(
      X = d$X[train_idx, ][d$A[train_idx] == 0, ],
      y = d$Y[train_idx][d$A[train_idx] == 0],
      loss_function = "log_loss",
      regularization = log(sum(1 - d$A[train_idx])) / sum(1 - d$A[train_idx]),
      worker_limit = 4,
      verbose = FALSE
    )

    # Get predictions
    prop_pred <- get_probabilities(prop_fit, newdata = d$X[test_idx, ])[, 2]
    mu1_pred <- get_probabilities(outcome_treated, newdata = d$X[test_idx, ])[, 2]
    mu0_pred <- get_probabilities(outcome_control, newdata = d$X[test_idx, ])[, 2]

    nuisances[[k]] <- list(
      prop = prop_pred,
      mu1 = mu1_pred,
      mu0 = mu0_pred,
      idx = test_idx
    )
  }

  times$nuisance_fitting <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("%.3fs\n", times$nuisance_fitting))
  flush.console()

  # 3. Compute scores
  cat("Step 3: Computing DML scores... ")
  flush.console()
  t1 <- Sys.time()

  scores <- rep(NA, n)
  for (k in 1:K) {
    idx <- nuisances[[k]]$idx
    A_k <- d$A[idx]
    Y_k <- d$Y[idx]
    prop_k <- nuisances[[k]]$prop
    mu1_k <- nuisances[[k]]$mu1
    mu0_k <- nuisances[[k]]$mu0

    # ATT score
    treated_idx <- A_k == 1
    n_treated <- sum(treated_idx)

    score_k <- rep(NA, length(idx))
    score_k[treated_idx] <- (Y_k[treated_idx] - mu0_k[treated_idx]) / n_treated
    score_k[!treated_idx] <- ((1 - prop_k[!treated_idx]) / prop_k[!treated_idx]) *
                              (Y_k[!treated_idx] - mu0_k[!treated_idx]) / n_treated

    scores[idx] <- score_k
  }

  times$score_computation <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("%.3fs\n", times$score_computation))
  flush.console()

  # 4. Compute theta
  cat("Step 4: Computing theta... ")
  flush.console()
  t1 <- Sys.time()
  theta <- mean(scores, na.rm = TRUE)
  times$theta_computation <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("%.3fs (theta=%.4f)\n", times$theta_computation, theta))
  flush.console()

  # 5. Compute variance
  cat("Step 5: Computing variance... ")
  flush.console()
  t1 <- Sys.time()
  sigma <- sd(scores, na.rm = TRUE) / sqrt(n)
  times$variance_computation <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("%.3fs (sigma=%.4f)\n", times$variance_computation, sigma))
  flush.console()

  list(success = TRUE, theta = theta, sigma = sigma, times = times)

}, error = function(e) {
  list(success = FALSE, error = conditionMessage(e), times = times)
})

times$total <- as.numeric(difftime(Sys.time(), times$start, units = "secs"))

cat("\n=== TIMING BREAKDOWN ===\n")
if (result$success) {
  cat(sprintf("Fold creation:        %.3fs\n", result$times$fold_creation))
  cat(sprintf("Nuisance fitting:     %.3fs\n", result$times$nuisance_fitting))
  cat(sprintf("Score computation:    %.3fs\n", result$times$score_computation))
  cat(sprintf("Theta computation:    %.3fs\n", result$times$theta_computation))
  cat(sprintf("Variance computation: %.3fs\n", result$times$variance_computation))
  cat(sprintf("TOTAL:                %.3fs\n", result$times$total))
  cat(sprintf("\nResult: theta = %.4f, sigma = %.4f\n", result$theta, result$sigma))
} else {
  cat(sprintf("FAILED: %s\n", result$error))
  cat("Completed steps:\n")
  print(names(result$times))
}

# Clean up
unlink(instrumented_file)

cat("\n=== DONE ===\n")
