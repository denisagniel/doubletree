#!/usr/bin/env Rscript

#' Test Phase 2 DGPs (7, 8, 9)
#'
#' Verify that:
#' - DGP7 (deep interaction): Tree should beat linear
#' - DGP8 (threshold): Tree should beat linear
#' - DGP9 (weak overlap): Stress test for both methods

suppressMessages({
  library(doubletree)
  library(dplyr)
})

# Source DGPs and methods
source("dgps/dgps_phase2.R")
source("methods/method_forest.R")
source("methods/method_linear.R")

cat("Testing Phase 2 DGPs\n")
cat(strrep("=", 70), "\n\n", sep = "")

# Test parameters
n <- 800  # Larger sample for stress test
K <- 5
seed_data <- 123
seed_method <- 456

# Function to test one DGP
test_dgp <- function(dgp_func, dgp_name) {
  cat(sprintf("Testing %s\n", dgp_name))
  cat(strrep("-", 70), "\n", sep = "")

  # Generate data
  data <- dgp_func(n, tau = 0.10, seed = seed_data)

  cat(sprintf("Data characteristics:\n"))
  cat(sprintf("  n = %d\n", n))
  cat(sprintf("  Y range: [%.2f, %.2f]\n", min(data$Y), max(data$Y)))
  cat(sprintf("  Treatment rate: %.1f%%\n", 100 * mean(data$A)))
  cat(sprintf("  Propensity range: [%.3f, %.3f]\n", min(data$true_e), max(data$true_e)))

  # Overlap diagnostic
  overlap_count <- sum(data$true_e > 0.1 & data$true_e < 0.9)
  overlap_pct <- 100 * overlap_count / n
  cat(sprintf("  Overlap (e in [0.1, 0.9]): %.1f%%\n", overlap_pct))

  cat(sprintf("  True ATT: %.4f\n\n", data$true_att))

  # Test tree method
  cat("Tree method:\n")

  # Auto-detect outcome type
  outcome_type <- if (all(data$Y %in% c(0, 1))) "binary" else "continuous"

  fit_tree <- tryCatch({
    doubletree::estimate_att(
      X = data$X, A = data$A, Y = data$Y,
      outcome_type = outcome_type,
      K = K, seed = seed_method
    )
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    NULL
  })

  if (!is.null(fit_tree)) {
    tree_rmse <- sqrt((fit_tree$theta - data$true_att)^2)
    tree_covers <- data$true_att >= fit_tree$ci[1] && data$true_att <= fit_tree$ci[2]
    cat(sprintf("  θ̂ = %.4f (RMSE: %.4f)\n", fit_tree$theta, tree_rmse))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", fit_tree$ci[1], fit_tree$ci[2]))
    cat(sprintf("  Covers truth: %s\n", ifelse(tree_covers, "✓", "✗")))
  } else {
    tree_rmse <- NA
    tree_covers <- FALSE
  }

  cat("\n")

  # Test linear method
  cat("Linear method:\n")
  fit_linear <- tryCatch({
    att_linear(data$X, data$A, data$Y, K = K, seed = seed_method)
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    NULL
  })

  if (!is.null(fit_linear)) {
    linear_rmse <- sqrt((fit_linear$theta - data$true_att)^2)
    linear_covers <- data$true_att >= fit_linear$ci[1] && data$true_att <= fit_linear$ci[2]
    cat(sprintf("  θ̂ = %.4f (RMSE: %.4f)\n", fit_linear$theta, linear_rmse))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n", fit_linear$ci[1], fit_linear$ci[2]))
    cat(sprintf("  Covers truth: %s\n", ifelse(linear_covers, "✓", "✗")))
  } else {
    linear_rmse <- NA
    linear_covers <- FALSE
  }

  cat("\n")

  # Comparison
  if (!is.null(fit_tree) && !is.null(fit_linear)) {
    rmse_ratio <- tree_rmse / linear_rmse
    cat("Comparison:\n")
    cat(sprintf("  Tree RMSE / Linear RMSE: %.2fx\n", rmse_ratio))
    if (rmse_ratio < 1.0) {
      cat("  → Tree WINS (lower RMSE) ✓\n")
    } else if (rmse_ratio < 2.0) {
      cat("  → Tree comparable (within 2x)\n")
    } else {
      cat("  → Linear wins (tree RMSE > 2x)\n")
    }
  }

  cat("\n", strrep("=", 70), "\n\n", sep = "")

  list(
    dgp = dgp_name,
    tree_rmse = tree_rmse,
    linear_rmse = linear_rmse,
    tree_covers = tree_covers,
    linear_covers = linear_covers,
    overlap_pct = overlap_pct
  )
}

# Run tests
results <- list()

results[[1]] <- test_dgp(generate_dgp7, "DGP7: Deep 3-way interaction")
results[[2]] <- test_dgp(generate_dgp8, "DGP8: Double nonlinearity (sin/cos)")
results[[3]] <- test_dgp(generate_dgp9, "DGP9: Weak overlap")

# Summary
cat(strrep("=", 70), "\n", sep = "")
cat("SUMMARY\n")
cat(strrep("=", 70), "\n", sep = "")

summary_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(
    DGP = r$dgp,
    Tree_RMSE = sprintf("%.4f", r$tree_rmse),
    Linear_RMSE = sprintf("%.4f", r$linear_rmse),
    Ratio = sprintf("%.2fx", r$tree_rmse / r$linear_rmse),
    Tree_Covers = ifelse(r$tree_covers, "✓", "✗"),
    Linear_Covers = ifelse(r$linear_covers, "✓", "✗"),
    Overlap = sprintf("%.1f%%", r$overlap_pct)
  )
}))

print(summary_df, row.names = FALSE)

cat("\n")
cat("Expected outcomes:\n")
cat("  DGP7 (deep interaction):   Tree RMSE < Linear RMSE (tree wins)\n")
cat("  DGP8 (double nonlinear):   Tree RMSE < Linear RMSE (tree wins)\n")
cat("  DGP9 (weak overlap):       Both maintain coverage (stress test)\n")
cat("\n")
cat("Note: DGP8 uses sin/cos in BOTH e(X) and m0(X) to break double robustness.\n")
cat("Linear misspecifies both models → cannot rely on DML's double robustness.\n")
cat("\n")

# Check if expectations met
dgp7_result <- results[[1]]
dgp8_result <- results[[2]]
dgp9_result <- results[[3]]

expectations_met <- TRUE

if (!is.na(dgp7_result$tree_rmse) && !is.na(dgp7_result$linear_rmse)) {
  if (dgp7_result$tree_rmse >= dgp7_result$linear_rmse) {
    cat("⚠ DGP7: Expected tree to beat linear, but didn't\n")
    expectations_met <- FALSE
  }
}

if (!is.na(dgp8_result$tree_rmse) && !is.na(dgp8_result$linear_rmse)) {
  if (dgp8_result$tree_rmse >= dgp8_result$linear_rmse) {
    cat("⚠ DGP8: Expected tree to beat linear, but didn't\n")
    expectations_met <- FALSE
  }
}

if (!dgp9_result$tree_covers || !dgp9_result$linear_covers) {
  cat("⚠ DGP9: Expected both methods to maintain coverage under weak overlap\n")
  expectations_met <- FALSE
}

if (expectations_met) {
  cat("✓ All expectations met!\n")
  cat("✓ Phase 2 DGPs ready for full simulation study\n")
}
