# Quick test of baseline methods (suppresses compilation output)
# Test that forest-DML and linear-DML run without errors

cat("Loading packages...\n")
suppressMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../../", quiet = TRUE)
})

source("../production/dgps/dgps_smooth.R")
source("../production/methods/method_forest_dml.R")
source("../production/methods/method_linear_dml.R")

cat("Packages loaded successfully.\n\n")

# Generate small test dataset
set.seed(123)
n <- 400
tau <- 0.10
d <- generate_dgp_binary_att(n, tau = tau, seed = 123)

cat("Test Data:\n")
cat(sprintf("  n = %d, n_treated = %d (%.1f%%)\n", n, sum(d$A), 100 * mean(d$A)))
cat(sprintf("  True ATT = %.4f\n\n", d$true_att))

# Test 1: Tree-DML
cat("Testing Tree-DML...\n")
result_tree <- tryCatch({
  suppressMessages({
    capture.output({
      doubletree::estimate_att(
        X = d$X, A = d$A, Y = d$Y,
        K = 5,
        regularization = log(n) / n,
        cv_regularization = FALSE,
        verbose = FALSE
      )
    })
  })
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  NULL
})

if (!is.null(result_tree)) {
  cat(sprintf("  ✓ Tree-DML: θ̂ = %.4f, SE = %.4f, CI = [%.4f, %.4f]\n",
              result_tree$theta, result_tree$sigma,
              result_tree$ci[1], result_tree$ci[2]))
  covers <- result_tree$ci[1] <= d$true_att && result_tree$ci[2] >= d$true_att
  cat(sprintf("  %s Coverage: %s\n\n",
              if(covers) "✓" else "✗",
              if(covers) "YES" else "NO"))
} else {
  cat("  ✗ Tree-DML FAILED\n\n")
}

# Test 2: Forest-DML
cat("Testing Forest-DML (ranger)...\n")
result_forest <- tryCatch({
  att_forest(
    X = d$X, A = d$A, Y = d$Y,
    K = 5,
    seed = 123,
    num.trees = 500,
    verbose = FALSE
  )
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  NULL
})

if (!is.null(result_forest)) {
  cat(sprintf("  ✓ Forest-DML: θ̂ = %.4f, SE = %.4f, CI = [%.4f, %.4f]\n",
              result_forest$theta, result_forest$sigma,
              result_forest$ci[1], result_forest$ci[2]))
  covers <- result_forest$ci[1] <= d$true_att && result_forest$ci[2] >= d$true_att
  cat(sprintf("  %s Coverage: %s\n\n",
              if(covers) "✓" else "✗",
              if(covers) "YES" else "NO"))
} else {
  cat("  ✗ Forest-DML FAILED\n\n")
}

# Test 3: Linear-DML
cat("Testing Linear-DML (GLM)...\n")
result_linear <- tryCatch({
  att_linear(
    X = d$X, A = d$A, Y = d$Y,
    K = 5,
    seed = 123,
    interactions = FALSE,
    verbose = FALSE
  )
}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  NULL
})

if (!is.null(result_linear)) {
  cat(sprintf("  ✓ Linear-DML: θ̂ = %.4f, SE = %.4f, CI = [%.4f, %.4f]\n",
              result_linear$theta, result_linear$sigma,
              result_linear$ci[1], result_linear$ci[2]))
  covers <- result_linear$ci[1] <= d$true_att && result_linear$ci[2] >= d$true_att
  cat(sprintf("  %s Coverage: %s\n\n",
              if(covers) "✓" else "✗",
              if(covers) "YES" else "NO"))
} else {
  cat("  ✗ Linear-DML FAILED\n\n")
}

# Summary
cat(strrep("=", 60), "\n")
cat("Summary\n")
cat(strrep("=", 60), "\n\n")

if (!is.null(result_tree) && !is.null(result_forest) && !is.null(result_linear)) {
  comparison <- data.frame(
    Method = c("Tree", "Forest", "Linear"),
    Estimate = c(result_tree$theta, result_forest$theta, result_linear$theta),
    SE = c(result_tree$sigma, result_forest$sigma, result_linear$sigma),
    Bias = c(result_tree$theta - d$true_att,
             result_forest$theta - d$true_att,
             result_linear$theta - d$true_att)
  )
  print(comparison, row.names = FALSE, digits = 4)

  cat("\n✓ All methods completed successfully\n")
  cat(sprintf("✓ Tree vs Forest bias difference: %.4f\n",
              abs(result_tree$theta - result_forest$theta)))
  cat(sprintf("✓ All estimates within [%.3f, %.3f] range\n",
              min(comparison$Estimate), max(comparison$Estimate)))

} else {
  cat("\n✗ One or more methods failed\n")
}
