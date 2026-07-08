# Test Lambda Sensitivity for Complex DGP
# Date: 2026-05-22
#
# Hypothesis: Fixed lambda=0.1 is too large (over-regularizes)
# Test: Run complex DGP at n=2000 with different lambda values
# Compare: Coverage rates, tree complexity, bias/SE

library(optimaltrees)
library(doubletree)

source("code/dgps.R")
source("code/estimators.R")

# Configuration
set.seed(2026522)
n <- 2000
dgp <- 4  # Complex DGP
n_reps <- 100  # Enough to detect coverage differences

# Lambda values to test
lambda_values <- c(
  0.1,                    # Current (baseline)
  0.05,                   # Intermediate
  log(n) / n              # Theory-consistent ≈ 0.0038
)
lambda_names <- c("lambda_0.1", "lambda_0.05", "lambda_theory")

# Approaches to test (fast ones)
approaches <- 1:3
approach_names <- c("full_sample", "crossfit", "doubletree")

cat("Lambda Sensitivity Test\n")
cat("======================\n\n")
cat("Configuration:\n")
cat("  DGP: Complex (dgp=4)\n")
cat("  n: ", n, "\n")
cat("  Replications: ", n_reps, "\n")
cat("  Lambda values:\n")
for (i in seq_along(lambda_values)) {
  cat(sprintf("    %s = %.6f\n", lambda_names[i], lambda_values[i]))
}
cat("\n")

# Storage
results <- vector("list", length(lambda_values) * length(approaches) * n_reps)
idx <- 1

# Run simulations
start_time <- Sys.time()

for (lam_idx in seq_along(lambda_values)) {
  lam <- lambda_values[lam_idx]
  lam_name <- lambda_names[lam_idx]

  cat(sprintf("\n--- Testing %s = %.6f ---\n", lam_name, lam))

  for (rep in 1:n_reps) {
    if (rep %% 10 == 0) cat("  Rep", rep, "/", n_reps, "\n")

    # Generate data
    seed <- 2026522000 + lam_idx * 10000 + rep
    set.seed(seed)
    data <- generate_dgp_complex(n = n)

    # True ATT
    theta_true <- data$theta_true

    # Test each approach
    for (app in approaches) {
      tryCatch({
        # Fit model with this lambda
        if (app == 1) {
          fit <- estimate_att_fullsample(data$X, data$A, data$Y, regularization = lam)
        } else if (app == 2) {
          fit <- estimate_att_crossfit(data$X, data$A, data$Y, K = 5, regularization = lam)
        } else if (app == 3) {
          fit <- estimate_att_doubletree(data$X, data$A, data$Y, K = 5, regularization = lam)
        }

        # Check for errors/failures
        if (!is.null(fit$error) || is.na(fit$theta) || is.na(fit$se)) {
          # Record failure
          results[[idx]] <- list(
            lambda = lam,
            lambda_name = lam_name,
            approach = app,
            approach_name = approach_names[app],
            rep = rep,
            seed = seed,
            theta_hat = NA,
            se = NA,
            theta_true = theta_true,
            coverage = NA,
            ci_width = NA,
            n_leaves_e = NA,
            n_leaves_m0 = NA,
            error = if (!is.null(fit$error)) fit$error else "NA result"
          )
        } else {
          # Compute coverage
          ci_lower <- fit$theta - 1.96 * fit$se
          ci_upper <- fit$theta + 1.96 * fit$se
          coverage <- (ci_lower <= theta_true) & (ci_upper >= theta_true)
          ci_width <- ci_upper - ci_lower

          # Extract tree complexity if available
          n_leaves_e <- NA
          n_leaves_m0 <- NA

          if (!is.null(fit$trees)) {
            if (!is.null(fit$trees$e)) {
              n_leaves_e <- length(fit$trees$e$leaf_indices)
            }
            if (!is.null(fit$trees$m0)) {
              n_leaves_m0 <- length(fit$trees$m0$leaf_indices)
            }
          } else if (!is.null(fit$structures)) {
            if (!is.null(fit$structures$e)) {
              n_leaves_e <- fit$structures$e@n_leaves
            }
            if (!is.null(fit$structures$m0)) {
              n_leaves_m0 <- fit$structures$m0@n_leaves
            }
          }

          # Record success
          results[[idx]] <- list(
            lambda = lam,
            lambda_name = lam_name,
            approach = app,
            approach_name = approach_names[app],
            rep = rep,
            seed = seed,
            theta_hat = fit$theta,
            se = fit$se,
            theta_true = theta_true,
            bias = fit$theta - theta_true,
            coverage = as.numeric(coverage),
            ci_width = ci_width,
            n_leaves_e = n_leaves_e,
            n_leaves_m0 = n_leaves_m0,
            error = NA
          )
        }
      }, error = function(e) {
        # Record error
        results[[idx]] <- list(
          lambda = lam,
          lambda_name = lam_name,
          approach = app,
          approach_name = approach_names[app],
          rep = rep,
          seed = seed,
          theta_hat = NA,
          se = NA,
          theta_true = theta_true,
          coverage = NA,
          ci_width = NA,
          n_leaves_e = NA,
          n_leaves_m0 = NA,
          error = e$message
        )
      })

      idx <- idx + 1
    }
  }
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n\nSimulations complete!\n")
cat("Elapsed time:", round(elapsed, 1), "minutes\n\n")

# Convert to data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))

# Save results
saveRDS(results_df, "results/lambda_sensitivity_test.rds")
cat("Results saved to: results/lambda_sensitivity_test.rds\n\n")

# Summary statistics
cat("======================\n")
cat("Summary by Lambda\n")
cat("======================\n\n")

for (lam_name in lambda_names) {
  cat("--- ", lam_name, " ---\n")
  subset_data <- results_df[results_df$lambda_name == lam_name, ]

  # By approach
  for (app_name in approach_names) {
    app_data <- subset_data[subset_data$approach_name == app_name, ]

    # Success rate
    n_success <- sum(!is.na(app_data$coverage))
    success_rate <- n_success / nrow(app_data)

    if (n_success > 0) {
      # Among successful fits
      coverage_rate <- mean(app_data$coverage, na.rm = TRUE)
      mean_bias <- mean(app_data$bias, na.rm = TRUE)
      mean_se <- mean(app_data$se, na.rm = TRUE)
      mean_ci_width <- mean(app_data$ci_width, na.rm = TRUE)
      mean_leaves_e <- mean(app_data$n_leaves_e, na.rm = TRUE)
      mean_leaves_m0 <- mean(app_data$n_leaves_m0, na.rm = TRUE)

      cat(sprintf("  %s:\n", app_name))
      cat(sprintf("    Success rate: %.1f%% (%d/%d)\n",
                  100 * success_rate, n_success, nrow(app_data)))
      cat(sprintf("    Coverage: %.3f\n", coverage_rate))
      cat(sprintf("    Mean bias: %.4f\n", mean_bias))
      cat(sprintf("    Mean SE: %.4f\n", mean_se))
      cat(sprintf("    Mean CI width: %.4f\n", mean_ci_width))
      cat(sprintf("    Mean leaves (e): %.1f\n", mean_leaves_e))
      cat(sprintf("    Mean leaves (m0): %.1f\n", mean_leaves_m0))
    } else {
      cat(sprintf("  %s: 0%% success (%d errors)\n", app_name, nrow(app_data)))
    }
  }
  cat("\n")
}

cat("======================\n")
cat("Coverage Comparison\n")
cat("======================\n\n")

# Create comparison table
comparison <- aggregate(coverage ~ lambda_name + approach_name,
                       data = results_df,
                       FUN = function(x) mean(x, na.rm = TRUE))
comparison_wide <- reshape(comparison,
                          idvar = "approach_name",
                          timevar = "lambda_name",
                          direction = "wide")
print(comparison_wide)

cat("\n\nTarget coverage: 0.95 (acceptable: 0.93-0.97)\n")

# Highlight improvements
cat("\n======================\n")
cat("Coverage Improvement\n")
cat("======================\n\n")

for (app_name in approach_names) {
  baseline <- comparison[comparison$lambda_name == "lambda_0.1" &
                         comparison$approach_name == app_name, "coverage"]
  theory <- comparison[comparison$lambda_name == "lambda_theory" &
                       comparison$approach_name == app_name, "coverage"]

  if (length(baseline) > 0 && length(theory) > 0 && !is.na(baseline) && !is.na(theory)) {
    improvement <- theory - baseline
    cat(sprintf("%s:\n", app_name))
    cat(sprintf("  Baseline (λ=0.1): %.3f\n", baseline))
    cat(sprintf("  Theory (λ≈0.004): %.3f\n", theory))
    cat(sprintf("  Improvement: %+.3f\n", improvement))

    if (theory >= 0.93 && theory <= 0.97) {
      cat("  ✓ Theory lambda achieves target coverage!\n")
    } else if (improvement > 0) {
      cat("  → Theory lambda improves coverage (but not to target yet)\n")
    }
    cat("\n")
  }
}

cat("\nTest complete! See results/lambda_sensitivity_test.rds for full data.\n")
