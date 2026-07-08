# Detailed Diagnostic: Show How Each Approach Works
# Run one replication with full tree details

library(doubletree)
library(optimaltrees)

# Source estimators
source("code/estimators.R")
source("code/dgps.R")

set.seed(12345)

# Generate data (DGP 1, n=500)
cat("\n=== GENERATING DATA ===\n")
cat("DGP: Simple confounding\n")
cat("n: 500\n")
cat("True ATT: 0.15\n\n")

dgp_result <- generate_dgp_simple(n = 500)
X <- dgp_result$X
A <- dgp_result$A
Y <- dgp_result$Y

cat("Sample sizes:\n")
cat("  Treated:", sum(A == 1), "\n")
cat("  Control:", sum(A == 0), "\n")
cat("Sample ATT (mean(Y|A=1) - mean(Y|A=0)):",
    mean(Y[A == 1]) - mean(Y[A == 0]), "\n")

cat("\n" , rep("=", 80), "\n\n", sep = "")

# ============================================================================
# APPROACH 1: FULLSAMPLE
# ============================================================================

cat("APPROACH 1: FULLSAMPLE\n")
cat(rep("-", 80), "\n", sep = "")
cat("Strategy: Fit nuisance functions on full sample, compute EIF\n")
cat("Trees: 2 (propensity e(X), outcome m0(X))\n\n")

result1 <- estimate_fullsample(X, A, Y, regularization = NULL, seed = 1)

cat("\n1. Propensity tree e(X):\n")
cat("   Lambda:", result1$trees$propensity$model@selected_regularization, "\n")
cat("   Number of leaves:", result1$trees$propensity$n_leaves, "\n")
cat("   Tree structure:\n")
print(result1$trees$propensity$splits)

cat("\n2. Outcome tree m0(X) [for control units]:\n")
cat("   Lambda:", result1$trees$outcome$model@selected_regularization, "\n")
cat("   Number of leaves:", result1$trees$outcome$n_leaves, "\n")
cat("   Tree structure:\n")
print(result1$trees$outcome$splits)

cat("\nEstimate computation:\n")
cat("  1. Predict e(X) for all units\n")
cat("  2. Predict m0(X) for all units\n")
cat("  3. Compute EIF: psi_i = A_i/e(X_i) * [Y_i - m0(X_i)] + m0(X_i) - theta\n")
cat("  4. theta_hat = mean(psi) among treated\n\n")

cat("RESULT: theta =", result1$estimate, "\n")
cat("        SE    =", result1$se, "\n")
cat("        CI    = [", result1$ci[1], ",", result1$ci[2], "]\n")

cat("\n", rep("=", 80), "\n\n", sep = "")

# ============================================================================
# APPROACH 2: CROSSFIT
# ============================================================================

cat("APPROACH 2: CROSSFIT\n")
cat(rep("-", 80), "\n", sep = "")
cat("Strategy: K-fold cross-fitting with CV per fold\n")
cat("Trees: 2K = 10 (K propensity + K outcome)\n\n")

result2 <- estimate_crossfit(X, A, Y, K = 5, regularization = NULL, seed = 2)

cat("Cross-fitting with K =", result2$K, "folds\n\n")

for (fold in 1:result2$K) {
  cat("Fold", fold, ":\n")
  cat("  Training size:", result2$fold_details[[fold]]$n_train, "\n")
  cat("  Test size:", result2$fold_details[[fold]]$n_test, "\n")

  cat("  Propensity tree:\n")
  cat("    Lambda:", result2$fold_details[[fold]]$propensity_lambda, "\n")
  cat("    Leaves:", result2$fold_details[[fold]]$propensity_leaves, "\n")

  cat("  Outcome tree:\n")
  cat("    Lambda:", result2$fold_details[[fold]]$outcome_lambda, "\n")
  cat("    Leaves:", result2$fold_details[[fold]]$outcome_leaves, "\n")
  cat("\n")
}

cat("Estimate computation:\n")
cat("  For each fold k:\n")
cat("    1. Train e_k(X) and m0_k(X) on fold k training data\n")
cat("    2. Predict on fold k test data\n")
cat("    3. Compute EIF for fold k test units\n")
cat("  4. Pool EIF scores across all folds\n")
cat("  5. theta_hat = mean(psi) among treated\n\n")

cat("RESULT: theta =", result2$estimate, "\n")
cat("        SE    =", result2$se, "\n")
cat("        CI    = [", result2$ci[1], ",", result2$ci[2], "]\n")

cat("\n", rep("=", 80), "\n\n", sep = "")

# ============================================================================
# APPROACH 3: DOUBLETREE (RASHOMON INTERSECTION)
# ============================================================================

cat("APPROACH 3: DOUBLETREE (RASHOMON INTERSECTION)\n")
cat(rep("-", 80), "\n", sep = "")
cat("Strategy: Single tree structure via Rashomon intersection + cross-fitting\n")
cat("Trees: Many collected, 1 used for inference\n\n")

result3 <- estimate_doubletree(X, A, Y,
                               rashomon_tol = 0.05,
                               cv_regularization = TRUE,
                               seed = 3)

cat("Rashomon tolerance:", result3$rashomon_tol, "\n")
cat("Use Rashomon:", result3$use_rashomon, "\n")
cat("Intersection succeeded:", result3$intersection_succeeded, "\n\n")

if (result3$intersection_succeeded) {
  cat("Final tree structure (from intersection):\n")
  cat("  Number of leaves:", length(result3$tree_structure$leaf_paths), "\n")
  cat("  Partition hash:", result3$tree_structure$partition_hash, "\n")
  cat("  Splits:\n")
  print(result3$tree_structure$splits)
} else {
  cat("Intersection failed, using fold-specific tree from fold 1\n")
  cat("  Number of leaves:", length(result3$tree_structure$leaf_paths), "\n")
  cat("  Splits:\n")
  print(result3$tree_structure$splits)
}

cat("\nEstimate computation:\n")
cat("  1. Collect Rashomon sets across K folds (within tolerance of best loss)\n")
cat("  2. Find common structure via intersection\n")
cat("  3. Refit leaf values using cross-fitting with common structure\n")
cat("  4. Compute EIF with single tree structure\n\n")

cat("RESULT: theta =", result3$estimate, "\n")
cat("        SE    =", result3$se, "\n")
cat("        CI    = [", result3$ci[1], ",", result3$ci[2], "]\n")

cat("\n", rep("=", 80), "\n\n", sep = "")

# ============================================================================
# APPROACH 4: DOUBLETREE AVERAGED (RASHOMON + AVERAGING)
# ============================================================================

cat("APPROACH 4: DOUBLETREE AVERAGED (RASHOMON + TREE AVERAGING)\n")
cat(rep("-", 80), "\n", sep = "")
cat("Strategy: Average K fold-specific trees with common structure\n")
cat("Trees: K trees with matched structure\n\n")

result4 <- estimate_doubletree_averaged(X, A, Y,
                                        rashomon_tol = 0.05,
                                        seed = 4)

cat("Method:", result4$method, "\n")
cat("Tier used:", result4$tier, "\n")
cat("Number of trees averaged:", result4$n_trees, "\n\n")

if (!is.null(result4$common_structure)) {
  cat("Common tree structure:\n")
  cat("  Number of leaves:", length(result4$common_structure$leaf_paths), "\n")
  cat("  Partition hash:", result4$common_structure$partition_hash, "\n")
  cat("  Splits:\n")
  print(result4$common_structure$splits)

  cat("\nFold-specific trees (before averaging):\n")
  for (k in 1:min(3, length(result4$fold_trees))) {
    cat("  Fold", k, ":\n")
    cat("    Leaves:", result4$fold_trees[[k]]$n_leaves, "\n")
    cat("    Matches common structure: TRUE\n")
  }
  if (length(result4$fold_trees) > 3) {
    cat("  ... and", length(result4$fold_trees) - 3, "more\n")
  }
}

cat("\nEstimate computation:\n")
cat("  1. Find common structure via Rashomon intersection (or fallback tiers)\n")
cat("  2. Collect K fold-specific trees matching this structure\n")
cat("  3. Average leaf values across K trees\n")
cat("  4. Compute EIF using averaged tree\n\n")

cat("RESULT: theta =", result4$estimate, "\n")
cat("        SE    =", result4$se, "\n")
cat("        CI    = [", result4$ci[1], ",", result4$ci[2], "]\n")

cat("\n", rep("=", 80), "\n\n", sep = "")

# ============================================================================
# APPROACH 5: M-SPLIT (MODAL STRUCTURE)
# ============================================================================

cat("APPROACH 5: M-SPLIT (MODAL STRUCTURE)\n")
cat(rep("-", 80), "\n", sep = "")
cat("Strategy: Repeated sample splitting, find modal tree structure\n")
cat("Trees: M×K trees, find most common structure\n\n")

result5 <- estimate_msplit(X, A, Y,
                           M = 10,
                           K = 5,
                           regularization = 0.1,
                           seed = 5)

cat("Number of splits (M):", result5$M, "\n")
cat("Folds per split (K):", result5$K, "\n")
cat("Total trees considered:", result5$M * result5$K, "\n\n")

cat("Modal structure found:", result5$modal_structure_found, "\n")
if (result5$modal_structure_found) {
  cat("Modal structure frequency:", result5$modal_structure_freq, "out of", result5$M * result5$K, "\n")
  cat("  (", round(100 * result5$modal_structure_freq / (result5$M * result5$K), 1), "%)\n\n", sep = "")

  cat("Modal tree structure:\n")
  cat("  Number of leaves:", length(result5$modal_structure$leaf_paths), "\n")
  cat("  Partition hash:", result5$modal_structure$partition_hash, "\n")
  cat("  Splits:\n")
  print(result5$modal_structure$splits)
}

cat("\nEstimate computation:\n")
cat("  1. For m = 1 to M:\n")
cat("       Split data into K folds\n")
cat("       For each fold, fit propensity and outcome trees\n")
cat("  2. Find most common tree structure across M×K trees\n")
cat("  3. Refit leaf values for modal structure using K-fold cross-fitting\n")
cat("  4. Compute EIF with modal tree structure\n\n")

cat("RESULT: theta =", result5$estimate, "\n")
cat("        SE    =", result5$se, "\n")
cat("        CI    = [", result5$ci[1], ",", result5$ci[2], "]\n")

cat("\n", rep("=", 80), "\n\n", sep = "")

# ============================================================================
# APPROACH 6: M-SPLIT AVERAGED
# ============================================================================

cat("APPROACH 6: M-SPLIT AVERAGED (MODAL STRUCTURE + TREE AVERAGING)\n")
cat(rep("-", 80), "\n", sep = "")
cat("Strategy: Average all trees matching modal structure\n")
cat("Trees: All M×K trees with modal structure\n\n")

result6 <- estimate_msplit_averaged(X, A, Y,
                                    M = 10,
                                    K = 5,
                                    seed = 6)

cat("Number of splits (M):", result6$M, "\n")
cat("Folds per split (K):", result6$K, "\n")
cat("Total trees considered:", result6$M * result6$K, "\n\n")

cat("Modal structure found:", result6$modal_structure_found, "\n")
if (result6$modal_structure_found) {
  cat("Modal structure frequency:", result6$modal_structure_freq, "out of", result6$M * result6$K, "\n")
  cat("Trees averaged:", result6$n_trees_averaged, "\n\n")

  cat("Modal tree structure:\n")
  cat("  Number of leaves:", length(result6$modal_structure$leaf_paths), "\n")
  cat("  Partition hash:", result6$modal_structure$partition_hash, "\n")
  cat("  Splits:\n")
  print(result6$modal_structure$splits)

  cat("\nTree averaging details:\n")
  cat("  All", result6$n_trees_averaged, "trees with modal structure averaged\n")
  cat("  Each tree's leaf values weighted equally (1/", result6$n_trees_averaged, ")\n", sep = "")
}

cat("\nEstimate computation:\n")
cat("  1. For m = 1 to M, k = 1 to K:\n")
cat("       Fit propensity and outcome trees with CV\n")
cat("  2. Find most common structure across M×K trees\n")
cat("  3. Filter trees to those matching modal structure\n")
cat("  4. Average leaf values across all matching trees\n")
cat("  5. Compute EIF with averaged tree\n\n")

cat("RESULT: theta =", result6$estimate, "\n")
cat("        SE    =", result6$se, "\n")
cat("        CI    = [", result6$ci[1], ",", result6$ci[2], "]\n")

cat("\n", rep("=", 80), "\n\n", sep = "")

# ============================================================================
# SUMMARY COMPARISON
# ============================================================================

cat("SUMMARY COMPARISON\n")
cat(rep("=", 80), "\n", sep = "")

results_df <- data.frame(
  Approach = c("1. Fullsample", "2. Crossfit", "3. Doubletree",
               "4. DT Averaged", "5. M-split", "6. MS Averaged"),
  Estimate = c(result1$estimate, result2$estimate, result3$estimate,
               result4$estimate, result5$estimate, result6$estimate),
  SE = c(result1$se, result2$se, result3$se,
         result4$se, result5$se, result6$se),
  CI_lower = c(result1$ci[1], result2$ci[1], result3$ci[1],
               result4$ci[1], result5$ci[1], result6$ci[1]),
  CI_upper = c(result1$ci[2], result2$ci[2], result3$ci[2],
               result4$ci[2], result5$ci[2], result6$ci[2]),
  Trees_fitted = c(2, 10, "Many", "K", "M×K", "M×K"),
  Trees_used = c(2, 10, 1, "K", 1, "Many")
)

print(results_df, row.names = FALSE)

cat("\nTrue ATT: 0.15\n\n")

cat("Key differences:\n")
cat("  - Approaches 1-2: Independent trees, no structure sharing\n")
cat("  - Approach 3: Single structure via intersection\n")
cat("  - Approach 4: K structures averaged (matched via intersection)\n")
cat("  - Approach 5: Single modal structure via voting\n")
cat("  - Approach 6: Many structures averaged (all matching modal)\n")

cat("\n", rep("=", 80), "\n", sep = "")
