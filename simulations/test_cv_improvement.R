#!/usr/bin/env Rscript
# Quick test: Compare fixed λ vs CV λ for coverage

suppressMessages(library(dplyr))
suppressMessages(library(optimaltrees))

# Source package functions
invisible(sapply(c(
  "R/estimate_att.R",
  "R/nuisance_trees.R",
  "R/score_att.R",
  "R/inference.R",
  "R/utils.R",
  "simulations/production/dgps/dgps_smooth.R"
), function(f) source(f, local = FALSE)))

set.seed(42)
n <- 800
n_sims <- 100

cat("═══════════════════════════════════════════════════════════════\n")
cat("Coverage Comparison: Fixed λ vs CV λ\n")
cat("═══════════════════════════════════════════════════════════════\n\n")
cat(sprintf("DGP: Binary features (DGP1)\n"))
cat(sprintf("Sample size: n = %d\n", n))
cat(sprintf("Simulations: %d\n", n_sims))
cat(sprintf("True ATT: τ = 0.10\n\n"))

# Storage
results_fixed <- data.frame(
  theta = numeric(n_sims),
  ci_lower = numeric(n_sims),
  ci_upper = numeric(n_sims),
  lambda_used = numeric(n_sims)
)

results_cv <- data.frame(
  theta = numeric(n_sims),
  ci_lower = numeric(n_sims),
  ci_upper = numeric(n_sims),
  lambda_used = numeric(n_sims)
)

cat("Running simulations...\n")
pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

for (i in 1:n_sims) {
  # Generate data
  d <- generate_dgp_binary_att(n = n, tau = 0.10, seed = 1000 + i)

  # Fixed λ = log(n)/n
  fit_fixed <- tryCatch({
    estimate_att(
      X = d$X, A = d$A, Y = d$Y, K = 5,
      regularization = log(n) / n,
      use_rashomon = FALSE,
      verbose = FALSE
    )
  }, error = function(e) NULL)

  if (!is.null(fit_fixed)) {
    results_fixed[i, ] <- c(
      fit_fixed$theta,
      fit_fixed$ci_95[1],
      fit_fixed$ci_95[2],
      log(n) / n
    )
  }

  # CV λ (auto-selected)
  fit_cv <- tryCatch({
    estimate_att(
      X = d$X, A = d$A, Y = d$Y, K = 5,
      cv_regularization = TRUE,
      cv_K = 5,
      use_rashomon = FALSE,
      verbose = FALSE
    )
  }, error = function(e) NULL)

  if (!is.null(fit_cv)) {
    # Extract median lambda from nuisance fits (approximate)
    # In practice, each fold may select different lambda
    results_cv[i, ] <- c(
      fit_cv$theta,
      fit_cv$ci_95[1],
      fit_cv$ci_95[2],
      NA  # Lambda varies by fold
    )
  }

  setTxtProgressBar(pb, i)
}
close(pb)

# Remove failed simulations
results_fixed <- results_fixed[results_fixed$theta != 0, ]
results_cv <- results_cv[results_cv$theta != 0, ]

cat("\n\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("RESULTS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Calculate metrics
true_att <- 0.10

metrics_fixed <- data.frame(
  method = "Fixed λ",
  n_valid = nrow(results_fixed),
  bias = mean(results_fixed$theta - true_att),
  rmse = sqrt(mean((results_fixed$theta - true_att)^2)),
  coverage = mean(results_fixed$ci_lower <= true_att & results_fixed$ci_upper >= true_att),
  mean_ci_width = mean(results_fixed$ci_upper - results_fixed$ci_lower),
  lambda = mean(results_fixed$lambda_used, na.rm = TRUE)
)

metrics_cv <- data.frame(
  method = "CV λ",
  n_valid = nrow(results_cv),
  bias = mean(results_cv$theta - true_att),
  rmse = sqrt(mean((results_cv$theta - true_att)^2)),
  coverage = mean(results_cv$ci_lower <= true_att & results_cv$ci_upper >= true_att),
  mean_ci_width = mean(results_cv$ci_upper - results_cv$ci_lower),
  lambda = NA
)

combined <- rbind(metrics_fixed, metrics_cv)

cat(sprintf("%-10s %8s %8s %8s %10s %10s %10s\n",
            "Method", "N Valid", "Bias", "RMSE", "Coverage", "CI Width", "Lambda"))
cat(strrep("-", 75), "\n")

for (i in 1:nrow(combined)) {
  cat(sprintf("%-10s %8d %8.4f %8.4f %9.1f%% %10.4f %10.6f\n",
              combined$method[i],
              combined$n_valid[i],
              combined$bias[i],
              combined$rmse[i],
              combined$coverage[i] * 100,
              combined$mean_ci_width[i],
              combined$lambda[i]))
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
cat("INTERPRETATION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

improvement_coverage <- (metrics_cv$coverage - metrics_fixed$coverage) * 100
improvement_bias <- abs(metrics_cv$bias) - abs(metrics_fixed$bias)

cat("Coverage improvement: ", sprintf("%+.1f%%", improvement_coverage), "\n")
cat("Bias change: ", sprintf("%+.4f", improvement_bias), " (negative = better)\n\n")

if (metrics_cv$coverage >= 0.93) {
  cat("✓ SUCCESS: CV achieves nominal coverage (≥93%)\n")
} else if (metrics_cv$coverage > metrics_fixed$coverage) {
  cat("✓ PARTIAL: CV improves coverage but not yet at 93-95% target\n")
  cat("  → May need more flexible models (BART, boosting)\n")
} else {
  cat("✗ NO IMPROVEMENT: CV doesn't help\n")
  cat("  → Issue may be in DML framework, not just λ choice\n")
}

cat("\n")
cat("═══════════════════════════════════════════════════════════════\n")
