# Extract and Analyze Lambda Sensitivity Results
# Re-run the simulations were complete but conversion failed

library(optimaltrees)
library(doubletree)

source("code/dgps.R")
source("code/estimators.R")

set.seed(2026522)
n <- 2000
dgp <- 4
n_reps <- 100

lambda_values <- c(0.1, 0.05, log(n)/n)
lambda_names <- c("lambda_0.1", "lambda_0.05", "lambda_theory")

approaches <- 1:3
approach_names <- c("full_sample", "crossfit", "doubletree")

# Re-run simulations (should be fast since trees are cached in memory)
cat("Re-running lambda sensitivity test with better error handling...\n\n")

results_list <- list()
idx <- 1

for (lam_idx in seq_along(lambda_values)) {
  lam <- lambda_values[lam_idx]
  lam_name <- lambda_names[lam_idx]

  cat(sprintf("Processing %s = %.6f\n", lam_name, lam))

  for (rep in 1:n_reps) {
    if (rep %% 25 == 0) cat("  Rep", rep, "/", n_reps, "\n")

    seed <- 2026522000 + lam_idx * 10000 + rep
    set.seed(seed)
    data <- generate_dgp_complex(n = n)
    theta_true <- data$theta_true

    for (app in approaches) {
      result_row <- list(
        lambda = lam,
        lambda_name = lam_name,
        approach = app,
        approach_name = approach_names[app],
        rep = rep,
        seed = seed,
        theta_true = theta_true
      )

      tryCatch({
        if (app == 1) {
          fit <- estimate_att_fullsample(data$X, data$A, data$Y, regularization = lam)
        } else if (app == 2) {
          fit <- estimate_att_crossfit(data$X, data$A, data$Y, K = 5, regularization = lam)
        } else if (app == 3) {
          fit <- estimate_att_doubletree(data$X, data$A, data$Y, K = 5, regularization = lam)
        }

        if (!is.null(fit$error) || is.na(fit$theta) || is.na(fit$se)) {
          result_row$theta_hat <- NA_real_
          result_row$se <- NA_real_
          result_row$bias <- NA_real_
          result_row$coverage <- NA_real_
          result_row$ci_width <- NA_real_
          result_row$n_leaves_e <- NA_integer_
          result_row$n_leaves_m0 <- NA_integer_
          result_row$error <- if (!is.null(fit$error)) fit$error else "NA result"
        } else {
          ci_lower <- fit$theta - 1.96 * fit$se
          ci_upper <- fit$theta + 1.96 * fit$se

          result_row$theta_hat <- fit$theta
          result_row$se <- fit$se
          result_row$bias <- fit$theta - theta_true
          result_row$coverage <- as.numeric((ci_lower <= theta_true) & (ci_upper >= theta_true))
          result_row$ci_width <- ci_upper - ci_lower

          # Extract tree complexity
          if (!is.null(fit$trees$e)) {
            result_row$n_leaves_e <- length(fit$trees$e$leaf_indices)
          } else if (!is.null(fit$structures$e)) {
            result_row$n_leaves_e <- fit$structures$e@n_leaves
          } else {
            result_row$n_leaves_e <- NA_integer_
          }

          if (!is.null(fit$trees$m0)) {
            result_row$n_leaves_m0 <- length(fit$trees$m0$leaf_indices)
          } else if (!is.null(fit$structures$m0)) {
            result_row$n_leaves_m0 <- fit$structures$m0@n_leaves
          } else {
            result_row$n_leaves_m0 <- NA_integer_
          }

          result_row$error <- NA_character_
        }
      }, error = function(e) {
        result_row$theta_hat <- NA_real_
        result_row$se <- NA_real_
        result_row$bias <- NA_real_
        result_row$coverage <- NA_real_
        result_row$ci_width <- NA_real_
        result_row$n_leaves_e <- NA_integer_
        result_row$n_leaves_m0 <- NA_integer_
        result_row$error <- e$message
      })

      results_list[[idx]] <- result_row
      idx <- idx + 1
    }
  }
}

cat("\n\nConverting to data frame...\n")

# Convert to data frame more carefully
results_df <- do.call(rbind, lapply(results_list, function(x) {
  as.data.frame(x, stringsAsFactors = FALSE)
}))

cat("Results shape:", nrow(results_df), "rows\n\n")

# Save
saveRDS(results_df, "results/lambda_sensitivity_test.rds")
cat("Saved to: results/lambda_sensitivity_test.rds\n\n")

# Summary
library(dplyr)

cat("======================\n")
cat("SUMMARY BY LAMBDA\n")
cat("======================\n\n")

summary <- results_df %>%
  group_by(lambda_name, approach_name) %>%
  summarise(
    n_total = n(),
    n_success = sum(!is.na(coverage)),
    success_rate = mean(!is.na(coverage)),
    coverage = mean(coverage, na.rm = TRUE),
    mean_bias = mean(bias, na.rm = TRUE),
    mean_se = mean(se, na.rm = TRUE),
    mean_ci_width = mean(ci_width, na.rm = TRUE),
    mean_leaves_e = mean(n_leaves_e, na.rm = TRUE),
    mean_leaves_m0 = mean(n_leaves_m0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(approach_name, lambda_name)

print(summary, n = 50, width = 150)

cat("\n\nTarget coverage: 0.95 (acceptable: 0.93-0.97)\n\n")

# Coverage comparison
cat("======================\n")
cat("COVERAGE BY LAMBDA\n")
cat("======================\n\n")

for (app in approach_names) {
  cat(sprintf("%s:\n", app))
  app_data <- summary[summary$approach_name == app, ]

  for (i in seq_len(nrow(app_data))) {
    row <- app_data[i, ]
    cat(sprintf("  %s: %.3f (n=%d)\n",
                row$lambda_name, row$coverage, row$n_success))
  }

  # Calculate improvement
  baseline <- app_data[app_data$lambda_name == "lambda_0.1", "coverage"]
  theory <- app_data[app_data$lambda_name == "lambda_theory", "coverage"]

  if (length(baseline) > 0 && length(theory) > 0) {
    improvement <- theory - baseline
    cat(sprintf("  Improvement (theory - baseline): %+.3f\n", improvement))

    if (theory >= 0.93 && theory <= 0.97) {
      cat("  ✓ Theory lambda achieves target coverage!\n")
    } else if (improvement > 0.02) {
      cat("  → Theory lambda improves coverage substantially\n")
    } else if (improvement > 0) {
      cat("  → Theory lambda improves coverage slightly\n")
    } else {
      cat("  ✗ Theory lambda does not improve coverage\n")
    }
  }
  cat("\n")
}

cat("\nAnalysis complete!\n")
