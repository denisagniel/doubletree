# diagnose_nuisance_mse.R
# Diagnostic: Measure nuisance estimation quality vs n and lambda
#
# Goal: Understand if trees are estimating nuisances well
# - Compute MSE of hat{e}, hat{m0} vs true nuisances
# - Test different lambda values
# - Check if MSE improves with n (as theory predicts)

# Setup ------------------------------------------------------------------------

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("devtools is required")
}
devtools::load_all()

if (!requireNamespace("treefarmr", quietly = TRUE)) {
  stop("treefarmr is required")
}

library(cli)

# Load DGP3 function
source("simulations/run_simulations_extended.R", local = TRUE)

# Diagnostic configuration -----------------------------------------------------

sample_sizes <- c(400, 800, 1600, 3200)
lambda_values <- c(0.01, 0.05, 0.1, 0.2, 0.5)  # Test range of regularization
n_reps <- 50  # Fewer reps for speed

cli_h1("Nuisance MSE Diagnostic")
cli_text("Goal: Measure nuisance estimation quality")
cli_text("Sample sizes: {paste(sample_sizes, collapse=', ')}")
cli_text("Lambda values: {paste(lambda_values, collapse=', ')}")
cli_text("Replications: {n_reps} per config")
cli_text("")

# Run diagnostic ---------------------------------------------------------------

results <- expand.grid(
  n = sample_sizes,
  lambda = lambda_values,
  rep = 1:n_reps,
  stringsAsFactors = FALSE
)
results$mse_e <- NA_real_
results$mse_m0 <- NA_real_
results$mse_m1 <- NA_real_
results$num_leaves_e <- NA_integer_
results$num_leaves_m0 <- NA_integer_
results$num_leaves_m1 <- NA_integer_

cli_progress_bar("Diagnostic runs", total = nrow(results))

for (i in 1:nrow(results)) {
  n <- results$n[i]
  lambda <- results$lambda[i]
  rep <- results$rep[i]
  seed <- 1000 * i + rep

  # Generate data with true nuisances
  data <- generate_data_dgp3(n = n, seed = seed)

  # Fit trees with this lambda (no cross-fitting, just train/test split)
  set.seed(seed)
  train_idx <- sample(1:n, size = floor(0.8 * n))
  test_idx <- setdiff(1:n, train_idx)

  X_train <- data$X[train_idx, , drop = FALSE]
  A_train <- data$A[train_idx]
  Y_train <- data$Y[train_idx]

  X_test <- data$X[test_idx, , drop = FALSE]
  true_e_test <- data$true_e[test_idx]
  true_m0_test <- data$true_m0[test_idx]
  true_m1_test <- data$true_m1[test_idx]

  # Fit propensity tree
  tryCatch({
    e_model <- treefarmr::fit_tree(X_train, A_train,
                                    loss_function = "log_loss",
                                    regularization = lambda,
                                    verbose = FALSE)
    e_pred <- predict(e_model, X_test, type = "prob")
    e_hat <- if (is.matrix(e_pred)) e_pred[, 2] else rep(0.5, length(test_idx))
    results$mse_e[i] <- mean((e_hat - true_e_test)^2)
    results$num_leaves_e[i] <- length(unique(e_model$leaf_ids))
  }, error = function(e) {
    results$mse_e[i] <<- NA_real_
  })

  # Fit m0 tree (control outcomes)
  idx0_train <- train_idx[A_train[train_idx %in% train_idx] == 0]
  if (length(idx0_train) > 10) {
    tryCatch({
      m0_model <- treefarmr::fit_tree(X_train[A_train == 0, , drop = FALSE],
                                       Y_train[A_train == 0],
                                       loss_function = "log_loss",
                                       regularization = lambda,
                                       verbose = FALSE)
      m0_pred <- predict(m0_model, X_test, type = "prob")
      m0_hat <- if (is.matrix(m0_pred)) m0_pred[, 2] else rep(0.5, length(test_idx))
      results$mse_m0[i] <- mean((m0_hat - true_m0_test)^2)
      results$num_leaves_m0[i] <- length(unique(m0_model$leaf_ids))
    }, error = function(e) {
      results$mse_m0[i] <<- NA_real_
    })
  }

  # Fit m1 tree (treated outcomes)
  idx1_train <- train_idx[A_train[train_idx %in% train_idx] == 1]
  if (length(idx1_train) > 10) {
    tryCatch({
      m1_model <- treefarmr::fit_tree(X_train[A_train == 1, , drop = FALSE],
                                       Y_train[A_train == 1],
                                       loss_function = "log_loss",
                                       regularization = lambda,
                                       verbose = FALSE)
      m1_pred <- predict(m1_model, X_test, type = "prob")
      m1_hat <- if (is.matrix(m1_pred)) m1_pred[, 2] else rep(0.5, length(test_idx))
      results$mse_m1[i] <- mean((m1_hat - true_m1_test)^2)
      results$num_leaves_m1[i] <- length(unique(m1_model$leaf_ids))
    }, error = function(e) {
      results$mse_m1[i] <<- NA_real_
    })
  }

  cli_progress_update()
}

cli_progress_done()

# Save results
saveRDS(results, "simulations/nuisance_mse_diagnostic.rds")
write.csv(results, "simulations/nuisance_mse_diagnostic.csv", row.names = FALSE)

# Summarize --------------------------------------------------------------------

cli_h1("Diagnostic Results")

# Aggregate by n and lambda
agg <- aggregate(cbind(mse_e, mse_m0, mse_m1, num_leaves_e, num_leaves_m0, num_leaves_m1) ~ n + lambda,
                 data = results, FUN = mean, na.rm = TRUE)

cli_h2("MSE by Sample Size (averaged over lambda)")
mse_by_n <- aggregate(cbind(mse_e, mse_m0, mse_m1) ~ n, data = agg, FUN = mean)
for (i in 1:nrow(mse_by_n)) {
  cli_text("n={mse_by_n$n[i]}: mse_e={round(mse_by_n$mse_e[i], 4)}, mse_m0={round(mse_by_n$mse_m0[i], 4)}, mse_m1={round(mse_by_n$mse_m1[i], 4)}")
}

cli_h2("MSE by Lambda (averaged over n)")
mse_by_lambda <- aggregate(cbind(mse_e, mse_m0, mse_m1) ~ lambda, data = agg, FUN = mean)
for (i in 1:nrow(mse_by_lambda)) {
  cli_text("lambda={mse_by_lambda$lambda[i]}: mse_e={round(mse_by_lambda$mse_e[i], 4)}, mse_m0={round(mse_by_lambda$mse_m0[i], 4)}, mse_m1={round(mse_by_lambda$mse_m1[i], 4)}")
}

cli_h2("Number of Leaves by Lambda")
leaves_by_lambda <- aggregate(cbind(num_leaves_e, num_leaves_m0, num_leaves_m1) ~ lambda,
                               data = agg, FUN = mean)
for (i in 1:nrow(leaves_by_lambda)) {
  cli_text("lambda={leaves_by_lambda$lambda[i]}: e={round(leaves_by_lambda$num_leaves_e[i], 1)}, m0={round(leaves_by_lambda$num_leaves_m0[i], 1)}, m1={round(leaves_by_lambda$num_leaves_m1[i], 1)}")
}

# Check if MSE decreases with n (theory prediction)
cli_h2("Rate Check: Does MSE decrease with n?")
mse_rate_e <- mse_by_n$mse_e[1] / mse_by_n$mse_e[nrow(mse_by_n)]
mse_rate_m0 <- mse_by_n$mse_m0[1] / mse_by_n$mse_m0[nrow(mse_by_n)]
cli_text("MSE(e) ratio (n=400 / n=3200): {round(mse_rate_e, 2)}x")
cli_text("MSE(m0) ratio (n=400 / n=3200): {round(mse_rate_m0, 2)}x")
if (mse_rate_e > 2 && mse_rate_m0 > 2) {
  cli_alert_success("✓ MSE decreases with n (consistent with theory)")
} else {
  cli_alert_warning("⚠ MSE not decreasing as expected")
}

# Best lambda
best_lambda_e <- mse_by_lambda$lambda[which.min(mse_by_lambda$mse_e)]
best_lambda_m0 <- mse_by_lambda$lambda[which.min(mse_by_lambda$mse_m0)]
cli_h2("Optimal Lambda (minimizes MSE)")
cli_alert_info("Best for e(X): lambda = {best_lambda_e}")
cli_alert_info("Best for m0(X): lambda = {best_lambda_m0}")

if (best_lambda_e < 0.1 || best_lambda_m0 < 0.1) {
  cli_alert_warning("⚠ Current default (0.1) may be too high!")
  cli_text("  Recommendation: Try lambda = {min(best_lambda_e, best_lambda_m0)} or use CV")
}

cli_h2("Saved")
cli_text("- simulations/nuisance_mse_diagnostic.rds")
cli_text("- simulations/nuisance_mse_diagnostic.csv")
