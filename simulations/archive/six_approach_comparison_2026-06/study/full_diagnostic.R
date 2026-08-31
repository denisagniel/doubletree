# Complete Diagnostic: All 6 Approaches with One Replication
# Shows exactly which trees contribute to each estimate

library(doubletree)
library(optimaltrees)

source("code/dgps.R")
source("code/estimators.R")

set.seed(42)  # Use seed that gives good modal structures

# Generate data
dgp <- generate_dgp_simple(n = 500)
X <- dgp$X
A <- dgp$A
Y <- dgp$Y

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("ONE REPLICATION: HOW EACH APPROACH WORKS\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("DGP: Simple (linear propensity, linear outcome)\n")
cat("n = 500 | Treated:", sum(A), "| Control:", sum(1-A), "\n")
cat("True ATT: 0.15\n")
cat("Sample ATT:", round(mean(Y[A==1]) - mean(Y[A==0]), 3), "\n\n")

# ============================================================================
# APPROACH 1: FULLSAMPLE
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 1: FULLSAMPLE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("HOW IT WORKS:\n")
cat("  1. Fit propensity tree e(X) on ALL units (n=500) with CV-selected lambda\n")
cat("  2. Fit outcome tree m0(X) on control units only with CV-selected lambda\n")
cat("  3. Predict e(X) and m0(X) for all units\n")
cat("  4. Compute EIF: psi_i = (A_i/pi) * (Y_i - m0(X_i)) + ...\n")
cat("  5. Estimate ATT = mean(psi_i)\n\n")

result1 <- estimate_att_fullsample(X, A, Y, regularization = NULL)

cat("TREES FITTED: 2\n")
cat("  1. Propensity e(X): lambda =", result1$propensity_lambda,
    "| leaves =", result1$propensity_leaves, "\n")
cat("  2. Outcome m0(X): lambda =", result1$outcome_lambda,
    "| leaves =", result1$outcome_leaves, "\n\n")

cat("ESTIMATE: theta =", round(result1$theta, 4), "± ", round(1.96*result1$se, 4), "\n\n")

# ============================================================================
# APPROACH 2: CROSSFIT
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 2: CROSSFIT (K=5)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("HOW IT WORKS:\n")
cat("  1. Split data into K=5 folds\n")
cat("  2. For each fold k:\n")
cat("       - Train e(X) and m0(X) on OTHER 4 folds with CV-selected lambda\n")
cat("       - Predict on fold k (test set)\n")
cat("  3. Pool predictions across all folds\n")
cat("  4. Compute EIF with pooled predictions\n")
cat("  5. Estimate ATT = mean(psi_i)\n\n")

result2 <- estimate_att_crossfit(X, A, Y, K = 5, regularization = NULL)

cat("TREES FITTED: 10 (5 propensity + 5 outcome)\n\n")
for (k in 1:5) {
  cat(sprintf("  Fold %d: e lambda=%s (leaves=%d) | m0 lambda=%s (leaves=%d)\n",
              k,
              round(result2$fold_details[[k]]$propensity_lambda, 4),
              result2$fold_details[[k]]$propensity_leaves,
              round(result2$fold_details[[k]]$outcome_lambda, 4),
              result2$fold_details[[k]]$outcome_leaves))
}
cat("\n")

cat("ESTIMATE: theta =", round(result2$theta, 4), "±", round(1.96*result2$se, 4), "\n\n")

# ============================================================================
# APPROACH 3: DOUBLETREE (package function)
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 3: DOUBLETREE (Rashomon Intersection)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("HOW IT WORKS:\n")
cat("  1. For each of K=5 folds:\n")
cat("       - Fit trees with CV-selected lambda\n")
cat("       - Collect Rashomon set (trees within 5% of best loss)\n")
cat("  2. Find COMMON structure via intersection across K Rashomon sets\n")
cat("  3. Refit leaf values using K-fold cross-fitting with common structure\n")
cat("  4. SINGLE common tree for both e(X) and m0(X)\n\n")

result3 <- doubletree::estimate_att(
  X, A, Y, K = 5, rashomon_tol = 0.05, cv_regularization = TRUE,
  use_rashomon = TRUE, verbose = FALSE
)

cat("OUTCOME:\n")
cat("  Intersection succeeded:", result3@intersection_succeeded, "\n")
if (result3@intersection_succeeded) {
  cat("  Propensity tree: leaves =", result3@n_leaves, "\n")
  cat("  Outcome tree: leaves =", result3@n_leaves, "\n")
} else {
  cat("  (Fell back to fold-specific trees)\n")
}
cat("\n")

cat("ESTIMATE: theta =", round(result3@estimate, 4), "±", round(1.96*result3@se, 4), "\n\n")

# ============================================================================
# APPROACH 5: M-SPLIT
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 5: M-SPLIT (Modal Structure)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("HOW IT WORKS:\n")
cat("  1. Repeat M=10 times:\n")
cat("       - Split data into K=5 folds\n")
cat("       - Fit propensity & outcome trees with CV-selected lambda\n")
cat("  2. Total: M×K = 50 trees fitted\n")
cat("  3. VOTE: Find most common structure across 50 trees\n")
cat("  4. Refit leaf values for modal structure using K-fold cross-fitting\n")
cat("  5. SINGLE modal tree for inference\n\n")

result5 <- doubletree::estimate_att_msplit(
  X, A, Y, M = 10, K = 5, verbose = FALSE
)

cat("OUTCOME:\n")
cat("  Total trees: M×K =", result5$M * result5$K, "\n")
cat("  Modal structure found:", !is.null(result5$structures$modal), "\n")

if (!is.null(result5$structures$modal)) {
  freq <- result5$structure_selection$modal_freq
  total <- result5$M * result5$K
  cat("  Modal frequency:", freq, "/", total, "=", round(100*freq/total, 1), "%\n")
  cat("  Propensity leaves:", length(result5$structures$modal$e$leaf_paths), "\n")
  cat("  Outcome leaves:", length(result5$structures$modal$m0$leaf_paths), "\n")
}
cat("\n")

cat("ESTIMATE: theta =", round(result5$theta, 4), "±", round(1.96*result5$sigma, 4), "\n\n")

# ============================================================================
# APPROACH 6: M-SPLIT AVERAGED
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 6: M-SPLIT AVERAGED (Tree Averaging)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("HOW IT WORKS:\n")
cat("  1. Repeat M=10 times:\n")
cat("       - Split data into K=5 folds\n")
cat("       - Fit propensity & outcome trees with CV-selected lambda\n")
cat("  2. Total: M×K = 50 trees fitted (each with optimal lambda)\n")
cat("  3. VOTE: Find most common structure\n")
cat("  4. FILTER: Keep only trees matching modal structure\n")
cat("  5. AVERAGE leaf values across ALL matching trees\n")
cat("  6. AVERAGED tree for inference\n\n")

result6 <- doubletree::estimate_att_msplit_averaged(
  X, A, Y, M = 10, K = 5, verbose = FALSE
)

cat("OUTCOME:\n")
cat("  Total trees: M×K =", result6$M * result6$K, "\n")
cat("  Modal structure found:", !is.null(result6$structures$modal), "\n")

if (!is.null(result6$structures$modal)) {
  freq <- result6$structure_selection$modal_freq
  n_avg <- result6$diagnostics$n_trees_averaged
  total <- result6$M * result6$K

  cat("  Modal frequency:", freq, "trees (", round(100*freq/total, 1), "%)\n", sep = "")
  cat("  Trees averaged:", n_avg, "\n")
  cat("  → Each tree contributes 1/", n_avg, " = ", round(100/n_avg, 1), "% to leaf values\n", sep = "")
  cat("  Propensity leaves:", length(result6$structures$modal$e$leaf_paths), "\n")
  cat("  Outcome leaves:", length(result6$structures$modal$m0$leaf_paths), "\n")
}
cat("\n")

cat("ESTIMATE: theta =", round(result6$theta, 4), "±", round(1.96*result6$sigma, 4), "\n\n")

# ============================================================================
# SUMMARY
# ============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("SUMMARY: KEY DIFFERENCES\n")
cat(rep("=", 80), "\n\n", sep = "")

summary_df <- data.frame(
  Approach = c("1. Fullsample", "2. Crossfit", "3. Doubletree",
               "5. M-split", "6. MS Averaged"),
  Estimate = c(result1$theta, result2$theta, result3@estimate,
               result5$theta, result6$theta),
  SE = c(result1$se, result2$se, result3@se,
         result5$sigma, result6$sigma),
  Trees_Fitted = c("2", "10", "K+Rash", "50", "50"),
  Trees_Used = c("2", "10 (separate)", "1 (common)", "1 (modal)", "Many (averaged)"),
  Strategy = c("Full sample", "Cross-fit", "Intersection", "Voting", "Voting+Average")
)

print(summary_df, row.names = FALSE)

cat("\nTrue ATT: 0.15\n")
cat("All estimates cover truth:",
    all(summary_df$Estimate - 1.96*summary_df$SE < 0.15 &
        summary_df$Estimate + 1.96*summary_df$SE > 0.15), "\n\n")

cat("KEY INSIGHTS:\n")
cat("  • Approaches 1-2: Independent trees, no structure sharing\n")
cat("  • Approach 3: SINGLE common structure via intersection\n")
cat("  • Approach 5: SINGLE modal structure via voting\n")
cat("  • Approach 6: AVERAGED tree (many trees with same structure)\n")
cat("  • Averaging (6) typically reduces SE compared to single tree (5)\n")

cat("\n", rep("=", 80), "\n", sep = "")
cat("Note: Approach 4 (doubletree_averaged) similar to 3 but averages K fold trees\n")
cat(rep("=", 80), "\n", sep = "")
