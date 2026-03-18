# Comprehensive baseline test: 100 replications on DGP 1
# Validates that forest-DML and linear-DML achieve reasonable coverage

library(dplyr)

cat("Loading functions...\n")
source("../production/dgps/dgps_smooth.R")
source("../production/methods/method_forest_dml.R")
source("../production/methods/method_linear_dml.R")

# Parameters
n <- 400
tau <- 0.10
K <- 5
n_reps <- 100
seed_offset <- 10000

cat(sprintf("\nRunning %d replications...\n", n_reps))
cat(sprintf("  DGP: Binary (4 features)\n"))
cat(sprintf("  n = %d, τ = %.2f\n", n, tau))
cat(sprintf("  Methods: Forest-DML, Linear-DML\n\n"))

# Storage
results <- data.frame(
  rep = integer(),
  method = character(),
  theta = numeric(),
  sigma = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  true_att = numeric(),
  stringsAsFactors = FALSE
)

# Run simulations
pb_update <- seq(10, n_reps, by = 10)
for (i in 1:n_reps) {
  # Generate data
  d <- generate_dgp_binary_att(n, tau = tau, seed = seed_offset + i)

  # Forest-DML
  res_forest <- tryCatch({
    dml_att_forest(d$X, d$A, d$Y, K = K, seed = seed_offset + i,
                   num.trees = 500, verbose = FALSE)
  }, error = function(e) NULL)

  if (!is.null(res_forest)) {
    results <- rbind(results, data.frame(
      rep = i,
      method = "Forest",
      theta = res_forest$theta,
      sigma = res_forest$sigma,
      ci_lower = res_forest$ci[1],
      ci_upper = res_forest$ci[2],
      true_att = d$true_att
    ))
  }

  # Linear-DML
  res_linear <- tryCatch({
    dml_att_linear(d$X, d$A, d$Y, K = K, seed = seed_offset + i,
                   verbose = FALSE)
  }, error = function(e) NULL)

  if (!is.null(res_linear)) {
    results <- rbind(results, data.frame(
      rep = i,
      method = "Linear",
      theta = res_linear$theta,
      sigma = res_linear$sigma,
      ci_lower = res_linear$ci[1],
      ci_upper = res_linear$ci[2],
      true_att = d$true_att
    ))
  }

  # Progress
  if (i %in% pb_update) {
    cat(sprintf("  Completed %d/%d\n", i, n_reps))
  }
}

cat("\nAnalyzing results...\n\n")

# Compute summary statistics
summary_stats <- results %>%
  group_by(method) %>%
  summarize(
    n_valid = n(),
    bias = mean(theta - true_att),
    rmse = sqrt(mean((theta - true_att)^2)),
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att),
    mean_ci_width = mean(ci_upper - ci_lower),
    .groups = "drop"
  )

# Print results
cat(strrep("=", 70), "\n")
cat("Summary Statistics (100 replications)\n")
cat(strrep("=", 70), "\n\n")
print(summary_stats, n = Inf)

# Quality checks
cat("\n")
cat(strrep("=", 70), "\n")
cat("Quality Checks\n")
cat(strrep("=", 70), "\n\n")

# Check coverage
for (method_name in unique(results$method)) {
  cov <- summary_stats %>% filter(method == method_name) %>% pull(coverage)
  bias_val <- summary_stats %>% filter(method == method_name) %>% pull(bias)

  if (cov >= 0.90 && cov <= 0.98) {
    cat(sprintf("✓ %s coverage: %.1f%% (good)\n", method_name, 100 * cov))
  } else {
    cat(sprintf("✗ %s coverage: %.1f%% (outside 90-98%%)\n", method_name, 100 * cov))
  }

  if (abs(bias_val) < 0.03) {
    cat(sprintf("✓ %s bias: %.4f (small)\n", method_name, bias_val))
  } else {
    cat(sprintf("⚠ %s bias: %.4f (moderate)\n", method_name, bias_val))
  }
  cat("\n")
}

# Overall assessment
cat(strrep("=", 70), "\n")
forest_ok <- summary_stats %>%
  filter(method == "Forest", coverage >= 0.90, coverage <= 0.98) %>%
  nrow() > 0

linear_ok <- summary_stats %>%
  filter(method == "Linear", coverage >= 0.90, coverage <= 0.98) %>%
  nrow() > 0

if (forest_ok && linear_ok) {
  cat("✓ BASELINE METHODS VALIDATED\n")
  cat("  Both methods achieve 90-98% coverage on DGP 1\n")
  cat("  Ready for production simulations\n")
} else {
  cat("⚠ ISSUES DETECTED\n")
  if (!forest_ok) cat("  - Forest-DML coverage outside acceptable range\n")
  if (!linear_ok) cat("  - Linear-DML coverage outside acceptable range\n")
}
cat(strrep("=", 70), "\n")
