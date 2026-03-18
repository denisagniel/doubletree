# Test baseline methods (forest-DML and linear-DML)
# Quick validation before full simulation runs

# Load packages and functions
devtools::load_all("../../../optimaltrees")
devtools::load_all("../../")  # Load doubletree
source("../production/dgps/dgps_smooth.R")
source("../production/methods/method_forest_dml.R")
source("../production/methods/method_linear_dml.R")

# Generate test data
set.seed(123)
n <- 400
tau <- 0.10
d <- generate_dgp_binary_att(n, tau = tau, seed = 123)

cat("Test data generated:\n")
cat(sprintf("  n = %d\n", n))
cat(sprintf("  n_treated = %d (%.1f%%)\n", sum(d$A), 100 * mean(d$A)))
cat(sprintf("  True ATT = %.4f\n\n", d$true_att))

# Test 1: Tree-DML (reference)
cat(strrep("=", 60), "\n")
cat("Test 1: Tree-DML (reference implementation)\n")
cat(strrep("=", 60), "\n")

result_tree <- doubletree::estimate_att(
  X = d$X,
  A = d$A,
  Y = d$Y,
  K = 5,
  regularization = log(n) / n,
  cv_regularization = FALSE,
  verbose = FALSE
)

cat(sprintf("Estimate: %.4f\n", result_tree$theta))
cat(sprintf("Std err:  %.4f\n", result_tree$sigma))
cat(sprintf("95%% CI:   [%.4f, %.4f]\n", result_tree$ci[1], result_tree$ci[2]))
cat(sprintf("Covers true ATT: %s\n\n",
            ifelse(d$true_att >= result_tree$ci[1] && d$true_att <= result_tree$ci[2],
                   "YES", "NO")))

# Test 2: Forest-DML
cat(strrep("=", 60), "\n")
cat("Test 2: Forest-DML (ranger baseline)\n")
cat(strrep("=", 60), "\n")

result_forest <- att_forest(
  X = d$X,
  A = d$A,
  Y = d$Y,
  K = 5,
  seed = 123,
  num.trees = 500,
  verbose = FALSE
)

cat(sprintf("Estimate: %.4f\n", result_forest$theta))
cat(sprintf("Std err:  %.4f\n", result_forest$sigma))
cat(sprintf("95%% CI:   [%.4f, %.4f]\n", result_forest$ci[1], result_forest$ci[2]))
cat(sprintf("Covers true ATT: %s\n\n",
            ifelse(d$true_att >= result_forest$ci[1] && d$true_att <= result_forest$ci[2],
                   "YES", "NO")))

# Test 3: Linear-DML
cat(strrep("=", 60), "\n")
cat("Test 3: Linear-DML (GLM baseline)\n")
cat(strrep("=", 60), "\n")

result_linear <- att_linear(
  X = d$X,
  A = d$A,
  Y = d$Y,
  K = 5,
  seed = 123,
  interactions = FALSE,
  verbose = FALSE
)

cat(sprintf("Estimate: %.4f\n", result_linear$theta))
cat(sprintf("Std err:  %.4f\n", result_linear$sigma))
cat(sprintf("95%% CI:   [%.4f, %.4f]\n", result_linear$ci[1], result_linear$ci[2]))
cat(sprintf("Covers true ATT: %s\n\n",
            ifelse(d$true_att >= result_linear$ci[1] && d$true_att <= result_linear$ci[2],
                   "YES", "NO")))

# Summary comparison
cat(strrep("=", 60), "\n")
cat("Summary: Method Comparison\n")
cat(strrep("=", 60), "\n\n")

comparison <- data.frame(
  Method = c("Tree-DML", "Forest-DML", "Linear-DML"),
  Estimate = c(result_tree$theta, result_forest$theta, result_linear$theta),
  SE = c(result_tree$sigma, result_forest$sigma, result_linear$sigma),
  CI_width = c(
    diff(result_tree$ci),
    diff(result_forest$ci),
    diff(result_linear$ci)
  ),
  Bias = c(
    result_tree$theta - d$true_att,
    result_forest$theta - d$true_att,
    result_linear$theta - d$true_att
  )
)

print(comparison, row.names = FALSE, digits = 4)

cat("\nExpected behavior:\n")
cat("- All methods should be reasonably close to true ATT (0.10)\n")
cat("- Tree and forest should be similar (bias within 20%)\n")
cat("- Linear may differ slightly (DGP 1 has modest nonlinearity)\n")
cat("- All CIs should cover true ATT (single replication, so not guaranteed)\n\n")

cat("Test complete. If all methods ran without error, baseline methods are ready.\n")
