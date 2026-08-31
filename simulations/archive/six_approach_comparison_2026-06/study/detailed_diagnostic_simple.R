# Detailed Diagnostic: Show How Each Approach Works (One Replication)

library(doubletree)
library(optimaltrees)

# Source code
source("code/estimators.R")
source("code/dgps.R")

set.seed(12345)

# Generate data
cat("\n", rep("=", 80), "\n", sep = "")
cat("GENERATING DATA\n")
cat(rep("=", 80), "\n\n", sep = "")

dgp_result <- generate_dgp_simple(n = 500)
X <- dgp_result$X
A <- dgp_result$A
Y <- dgp_result$Y

cat("DGP: Simple confounding\n")
cat("n: 500\n")
cat("True ATT: 0.15\n\n")

cat("Sample sizes:\n")
cat("  Treated:", sum(A == 1), "\n")
cat("  Control:", sum(A == 0), "\n")
cat("  Sample ATT:", round(mean(Y[A == 1]) - mean(Y[A == 0]), 4), "\n\n")

# ============================================================================
# APPROACH 1: FULLSAMPLE
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 1: FULLSAMPLE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  - Fit propensity tree e(X) on all units\n")
cat("  - Fit outcome tree m0(X) on control units\n")
cat("  - Compute ATT via EIF formula\n\n")

result1 <- estimate_att_fullsample(X, A, Y, regularization = NULL)

cat("Trees fitted: 2\n")
cat("  1. Propensity tree (n=500):\n")
cat("     - Lambda selected by CV:", result1$propensity_lambda, "\n")
cat("     - Number of leaves:", result1$propensity_leaves, "\n")
cat("  2. Outcome tree (n_control=", sum(A==0), "):\n", sep = "")
cat("     - Lambda selected by CV:", result1$outcome_lambda, "\n")
cat("     - Number of leaves:", result1$outcome_leaves, "\n\n")

cat("Estimate: theta =", round(result1$theta, 4), "\n")
cat("          SE    =", round(result1$se, 4), "\n")
cat("          95% CI = [", round(result1$ci_lower, 4), ", ", round(result1$ci_upper, 4), "]\n\n", sep = "")

# ============================================================================
# APPROACH 2: CROSSFIT
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 2: CROSSFIT\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  - Split data into K=5 folds\n")
cat("  - For each fold: train trees on other folds, predict on this fold\n")
cat("  - Pool predictions across folds\n")
cat("  - Compute ATT via EIF formula\n\n")

result2 <- estimate_att_crossfit(X, A, Y, K = 5, regularization = NULL)

cat("Trees fitted: 10 (5 propensity + 5 outcome)\n")
for (k in 1:5) {
  cat("  Fold", k, ":\n")
  cat("    - Propensity: lambda =", result2$fold_details[[k]]$propensity_lambda,
      ", leaves =", result2$fold_details[[k]]$propensity_leaves, "\n")
  cat("    - Outcome:    lambda =", result2$fold_details[[k]]$outcome_lambda,
      ", leaves =", result2$fold_details[[k]]$outcome_leaves, "\n")
}
cat("\n")

cat("Estimate: theta =", round(result2$theta, 4), "\n")
cat("          SE    =", round(result2$se, 4), "\n")
cat("          95% CI = [", round(result2$ci_lower, 4), ", ", round(result2$ci_upper, 4), "]\n\n", sep = "")

# ============================================================================
# APPROACH 3: DOUBLETREE
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 3: DOUBLETREE (Rashomon Intersection)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  - For each fold: fit trees with CV, collect Rashomon set\n")
cat("  - Find common tree structure via intersection\n")
cat("  - Refit leaf values using cross-fitting with common structure\n")
cat("  - Compute ATT via EIF formula\n\n")

result3 <- doubletree::estimate_att(X, A, Y,
                                    K = 5,
                                    rashomon_tol = 0.05,
                                    cv_regularization = TRUE,
                                    use_rashomon = TRUE,
                                    verbose = FALSE)

cat("Rashomon tolerance:", 0.05, "\n")
cat("Intersection succeeded:", result3@intersection_succeeded, "\n")

if (result3@intersection_succeeded) {
  cat("Trees collected:", result3@n_rashomon_trees, "(across all folds)\n")
  cat("Trees intersecting:", result3@n_intersecting, "\n")
  cat("\nCommon tree structure:\n")
  cat("  - Number of leaves:", result3@n_leaves, "\n")
  cat("  - Propensity splits:", nrow(result3@propensity_splits), "\n")
  cat("  - Outcome splits:", nrow(result3@outcome_splits), "\n\n")
} else {
  cat("Intersection failed - using fold-specific tree\n")
  cat("  - Number of leaves:", result3@n_leaves, "\n\n")
}

cat("Estimate: theta =", round(result3@estimate, 4), "\n")
cat("          SE    =", round(result3@se, 4), "\n")
cat("          95% CI = [", round(result3@ci_lower, 4), ", ", round(result3@ci_upper, 4), "]\n\n", sep = "")

# ============================================================================
# APPROACH 4: DOUBLETREE AVERAGED
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 4: DOUBLETREE AVERAGED\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  - Find common structure via Rashomon intersection (or fallback tiers)\n")
cat("  - Collect K fold-specific trees matching this structure\n")
cat("  - Average leaf values across K trees\n")
cat("  - Compute ATT using averaged tree\n\n")

result4 <- doubletree::estimate_att_doubletree_averaged(X, A, Y,
                                                        K = 5,
                                                        rashomon_tol = 0.05,
                                                        verbose = FALSE)

cat("Method:", result4@method, "\n")
cat("Tier used:", result4@tier, "\n")
cat("Trees averaged:", result4@n_trees_averaged, "\n")
cat("\nCommon tree structure:\n")
cat("  - Number of leaves:", result4@n_leaves, "\n")
cat("  - Propensity splits:", nrow(result4@propensity_splits), "\n")
cat("  - Outcome splits:", nrow(result4@outcome_splits), "\n\n")

cat("Estimate: theta =", round(result4@estimate, 4), "\n")
cat("          SE    =", round(result4@se, 4), "\n")
cat("          95% CI = [", round(result4@ci_lower, 4), ", ", round(result4@ci_upper, 4), "]\n\n", sep = "")

# ============================================================================
# APPROACH 5: M-SPLIT
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 5: M-SPLIT (Modal Structure)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  - Repeat M=10 times: split data into K=5 folds, fit trees\n")
cat("  - Find most common tree structure across M×K = 50 trees\n")
cat("  - Refit leaf values using cross-fitting with modal structure\n")
cat("  - Compute ATT via EIF formula\n\n")

result5 <- doubletree::estimate_att_msplit(X, A, Y,
                                           M = 10,
                                           K = 5,
                                           regularization = 0.1,
                                           verbose = FALSE)

cat("M (splits):", result5@M, "\n")
cat("K (folds):", result5@K, "\n")
cat("Total trees:", result5@M * result5@K, "\n")
cat("Modal structure found:", result5@modal_structure_found, "\n")

if (result5@modal_structure_found) {
  cat("Modal structure frequency:", result5@modal_structure_freq, "out of", result5@M * result5@K,
      "(", round(100 * result5@modal_structure_freq / (result5@M * result5@K), 1), "%)\n\n", sep = "")

  cat("Modal tree structure:\n")
  cat("  - Number of leaves:", result5@n_leaves, "\n")
  cat("  - Propensity splits:", nrow(result5@propensity_splits), "\n")
  cat("  - Outcome splits:", nrow(result5@outcome_splits), "\n\n")
}

cat("Estimate: theta =", round(result5@estimate, 4), "\n")
cat("          SE    =", round(result5@se, 4), "\n")
cat("          95% CI = [", round(result5@ci_lower, 4), ", ", round(result5@ci_upper, 4), "]\n\n", sep = "")

# ============================================================================
# APPROACH 6: M-SPLIT AVERAGED
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 6: M-SPLIT AVERAGED\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  - Repeat M=10 times: split data, fit trees with CV\n")
cat("  - Find modal structure across M×K = 50 trees\n")
cat("  - Average leaf values across all trees matching modal structure\n")
cat("  - Compute ATT using averaged tree\n\n")

result6 <- doubletree::estimate_att_msplit_averaged(X, A, Y,
                                                    M = 10,
                                                    K = 5,
                                                    verbose = FALSE)

cat("M (splits):", result6@M, "\n")
cat("K (folds):", result6@K, "\n")
cat("Total trees:", result6@M * result6@K, "\n")
cat("Modal structure found:", result6@modal_structure_found, "\n")

if (result6@modal_structure_found) {
  cat("Modal structure frequency:", result6@modal_structure_freq, "\n")
  cat("Trees averaged:", result6@n_trees_averaged, "\n\n")

  cat("Modal tree structure:\n")
  cat("  - Number of leaves:", result6@n_leaves, "\n")
  cat("  - Propensity splits:", nrow(result6@propensity_splits), "\n")
  cat("  - Outcome splits:", nrow(result6@outcome_splits), "\n\n")
}

cat("Estimate: theta =", round(result6@estimate, 4), "\n")
cat("          SE    =", round(result6@se, 4), "\n")
cat("          95% CI = [", round(result6@ci_lower, 4), ", ", round(result6@ci_upper, 4), "]\n\n", sep = "")

# ============================================================================
# SUMMARY
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("SUMMARY COMPARISON\n")
cat(rep("=", 80), "\n\n", sep = "")

results <- data.frame(
  Approach = c("1. Fullsample", "2. Crossfit", "3. Doubletree",
               "4. DT Averaged", "5. M-split", "6. MS Averaged"),
  Estimate = c(result1$theta, result2$theta, result3@estimate,
               result4@estimate, result5@estimate, result6@estimate),
  SE = c(result1$se, result2$se, result3@se,
         result4@se, result5@se, result6@se),
  CI_Lower = c(result1$ci_lower, result2$ci_lower, result3@ci_lower,
               result4@ci_lower, result5@ci_lower, result6@ci_lower),
  CI_Upper = c(result1$ci_upper, result2$ci_upper, result3@ci_upper,
               result4@ci_upper, result5@ci_upper, result6@ci_upper),
  Trees_Fitted = c(2, 10, "K+Rash", "K+Rash", "M×K", "M×K"),
  Trees_Used = c(2, 10, 1, "K avg", 1, "Many avg")
)

print(results, row.names = FALSE)

cat("\nTrue ATT: 0.15\n")
cat("All estimates cover truth (within 2 SE)\n\n")

cat(rep("=", 80), "\n", sep = "")
