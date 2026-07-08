# Quick Lambda Sensitivity Test - Minimal Version
# Just enough reps to detect coverage differences

library(optimaltrees)
library(doubletree)

source("code/dgps.R")
source("code/estimators.R")

set.seed(20265220)
n <- 2000
n_reps <- 30  # Enough to estimate coverage ±0.18

lambda_values <- c(0.1, 0.05)  # Skip theory lambda for now
lambda_names <- c("lambda_0.1", "lambda_0.05")

approaches <- 1:3
approach_names <- c("full_sample", "crossfit", "doubletree")

cat("Quick Lambda Test (30 reps)\n")
cat("===========================\n\n")

results_list <- list()
idx <- 1

for (lam_idx in seq_along(lambda_values)) {
  lam <- lambda_values[lam_idx]
  lam_name <- lambda_names[lam_idx]

  cat(sprintf("%s = %.3f\n", lam_name, lam))

  for (rep in 1:n_reps) {
    seed <- 20265220 + lam_idx * 1000 + rep
    set.seed(seed)
    data <- generate_dgp_complex(n = n)
    theta_true <- data$true_att  # Fixed: was data$theta_true

    for (app in approaches) {
      tryCatch({
        if (app == 1) {
          fit <- estimate_att_fullsample(data$X, data$A, data$Y, regularization = lam)
        } else if (app == 2) {
          fit <- estimate_att_crossfit(data$X, data$A, data$Y, K = 5, regularization = lam)
        } else {
          fit <- estimate_att_doubletree(data$X, data$A, data$Y, K = 5, regularization = lam)
        }

        if (!is.na(fit$theta) && !is.na(fit$se)) {
          ci_lower <- fit$theta - 1.96 * fit$se
          ci_upper <- fit$theta + 1.96 * fit$se
          coverage <- (ci_lower <= theta_true) & (ci_upper >= theta_true)

          results_list[[idx]] <- data.frame(
            lambda = lam,
            lambda_name = lam_name,
            approach = app,
            approach_name = approach_names[app],
            rep = rep,
            theta_hat = fit$theta,
            se = fit$se,
            theta_true = theta_true,
            bias = fit$theta - theta_true,
            coverage = as.numeric(coverage),
            ci_width = ci_upper - ci_lower,
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) {
        # Skip errors
      })

      idx <- idx + 1
    }
  }
}

cat("\n\nAnalyzing...\n")

results_df <- do.call(rbind, results_list)

library(dplyr)

summary <- results_df %>%
  group_by(lambda_name, approach_name) %>%
  summarise(
    n = n(),
    coverage = mean(coverage),
    mean_bias = mean(bias),
    mean_se = mean(se),
    .groups = "drop"
  )

cat("\n======================\n")
cat("RESULTS\n")
cat("======================\n\n")

print(summary, n=50)

cat("\n\nTarget: 0.95\n")
cat("Current (λ=0.1): from full simulation\n\n")

# Compare
cat("Coverage Change (λ=0.05 vs λ=0.1):\n")
for (app in approach_names) {
  baseline <- summary$coverage[summary$lambda_name == "lambda_0.1" & summary$approach_name == app]
  new_val <- summary$coverage[summary$lambda_name == "lambda_0.05" & summary$approach_name == app]

  if (length(baseline) > 0 && length(new_val) > 0) {
    cat(sprintf("  %s: %.3f → %.3f (%+.3f)\n", app, baseline, new_val, new_val - baseline))
  }
}

saveRDS(results_df, "results/lambda_quick_test.rds")
cat("\n\nSaved to: results/lambda_quick_test.rds\n")
