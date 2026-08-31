# Run One Replication - Show How Each Approach Works

library(doubletree)
library(optimaltrees)

source("code/dgps.R")

set.seed(888)

# Generate data
dgp <- generate_dgp_simple(n = 500)
X <- dgp$X
A <- dgp$A
Y <- dgp$Y

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("DATA: Simple DGP, n=500\n")
cat(rep("=", 80), "\n", sep = "")
cat("Treated:", sum(A), "| Control:", sum(1-A), "\n")
cat("True ATT: 0.15 | Sample ATT:", round(mean(Y[A==1]) - mean(Y[A==0]), 3), "\n\n")

# ==================================================================
# APPROACH 4: DOUBLETREE AVERAGED
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 4: DOUBLETREE AVERAGED\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  1. Find common tree structure (Rashomon intersection or fallback)\n")
cat("  2. Collect K=5 trees matching this structure\n")
cat("  3. Average leaf values across K trees\n")
cat("  4. Compute ATT using averaged tree\n\n")

result4 <- doubletree::estimate_att_doubletree_averaged(
  X, A, Y,
  K = 5,
  rashomon_tol = 0.05,
  verbose = FALSE
)

cat("Result:\n")
cat("  Method:", result4@method, "\n")
cat("  Tier:", result4@tier, "\n")
cat("  Trees averaged:", result4@n_trees_averaged, "\n")
cat("  Leaves:", result4@n_leaves, "\n\n")

if (result4@tier <= 2) {
  cat("✓ Rashomon intersection succeeded\n")
  cat("  → Common structure found across", result4@n_trees_averaged, "folds\n")
  cat("  → Each fold contributed one tree with this structure\n")
  cat("  → Averaged leaf values across all", result4@n_trees_averaged, "trees\n\n")
} else {
  cat("✗ Rashomon intersection failed\n")
  cat("  → Fallback tier", result4@tier, "\n")
  cat("  → Averaged", result4@n_trees_averaged, "fold-specific trees\n\n")
}

cat("Propensity tree structure (", nrow(result4@propensity_splits), " splits):\n", sep = "")
print(result4@propensity_splits)
cat("\nOutcome tree structure (", nrow(result4@outcome_splits), " splits):\n", sep = "")
print(result4@outcome_splits)

theta4 <- result4@estimate
se4 <- result4@se
cat("\nEstimate:", round(theta4, 4), "±", round(1.96*se4, 4),
    "[", round(theta4 - 1.96*se4, 4), ",", round(theta4 + 1.96*se4, 4), "]\n\n")

# ==================================================================
# APPROACH 5: M-SPLIT
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 5: M-SPLIT (Modal Structure)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  1. Repeat M=10 times: split data, fit K=5 trees\n")
cat("  2. Total M×K = 50 trees (with fixed lambda = 0.1)\n")
cat("  3. Find most common structure (voting)\n")
cat("  4. Refit leaf values using cross-fitting\n")
cat("  5. Compute ATT with single modal tree\n\n")

result5 <- doubletree::estimate_att_msplit(
  X, A, Y,
  M = 10,
  K = 5,
  regularization = 0.1,
  verbose = FALSE
)

cat("Result:\n")
cat("  M × K:", result5@M, "×", result5@K, "=", result5@M * result5@K, "trees\n")
cat("  Modal structure found:", result5@modal_structure_found, "\n")

if (result5@modal_structure_found) {
  freq <- result5@modal_structure_freq
  total <- result5@M * result5@K
  pct <- round(100 * freq / total, 1)
  cat("  Modal frequency:", freq, "/", total, "(", pct, "%)\n", sep = "")
  cat("  Leaves:", result5@n_leaves, "\n\n")

  cat("✓ Modal structure selected\n")
  cat("  →", freq, "out of", total, "trees had this structure\n")
  cat("  → Used cross-fitting to refit leaf values\n")
  cat("  → Single tree used for inference\n\n")

  cat("Propensity tree structure (", nrow(result5@propensity_splits), " splits):\n", sep = "")
  print(result5@propensity_splits)
  cat("\nOutcome tree structure (", nrow(result5@outcome_splits), " splits):\n", sep = "")
  print(result5@outcome_splits)
}

theta5 <- result5@estimate
se5 <- result5@se
cat("\nEstimate:", round(theta5, 4), "±", round(1.96*se5, 4),
    "[", round(theta5 - 1.96*se5, 4), ",", round(theta5 + 1.96*se5, 4), "]\n\n")

# ==================================================================
# APPROACH 6: M-SPLIT AVERAGED
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 6: M-SPLIT AVERAGED\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Strategy:\n")
cat("  1. Repeat M=10 times: split data, fit K=5 trees with CV\n")
cat("  2. Total M×K = 50 trees (each with CV-selected lambda)\n")
cat("  3. Find most common structure (voting)\n")
cat("  4. Filter to ALL trees matching modal structure\n")
cat("  5. Average leaf values across all matching trees\n")
cat("  6. Compute ATT with averaged tree\n\n")

result6 <- doubletree::estimate_att_msplit_averaged(
  X, A, Y,
  M = 10,
  K = 5,
  verbose = FALSE
)

cat("Result:\n")
cat("  M × K:", result6@M, "×", result6@K, "=", result6@M * result6@K, "trees\n")
cat("  Modal structure found:", result6@modal_structure_found, "\n")

if (result6@modal_structure_found) {
  freq <- result6@modal_structure_freq
  avg_n <- result6@n_trees_averaged
  total <- result6@M * result6@K
  pct <- round(100 * freq / total, 1)

  cat("  Modal frequency:", freq, "(", pct, "%)\n", sep = "")
  cat("  Trees averaged:", avg_n, "\n")
  cat("  Leaves:", result6@n_leaves, "\n\n")

  cat("✓ Modal structure selected and averaged\n")
  cat("  →", freq, "out of", total, "trees had modal structure\n")
  cat("  → Averaged leaf values across", avg_n, "of these trees\n")
  cat("  → Each tree weighted equally (1/", avg_n, ")\n\n", sep = "")

  cat("Propensity tree structure (", nrow(result6@propensity_splits), " splits):\n", sep = "")
  print(result6@propensity_splits)
  cat("\nOutcome tree structure (", nrow(result6@outcome_splits), " splits):\n", sep = "")
  print(result6@outcome_splits)
}

theta6 <- result6@estimate
se6 <- result6@se
cat("\nEstimate:", round(theta6, 4), "±", round(1.96*se6, 4),
    "[", round(theta6 - 1.96*se6, 4), ",", round(theta6 + 1.96*se6, 4), "]\n\n")

# ==================================================================
# SUMMARY
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("True ATT: 0.15\n\n")

summary_df <- data.frame(
  Approach = c("4. DT Averaged", "5. M-split", "6. MS Averaged"),
  Estimate = c(theta4, theta5, theta6),
  SE = c(se4, se5, se6),
  Trees_Used = c(paste(result4@n_trees_averaged, "averaged"),
                 "1 modal",
                 paste(result6@n_trees_averaged, "averaged")),
  Covers_Truth = c(
    theta4 - 1.96*se4 < 0.15 & theta4 + 1.96*se4 > 0.15,
    theta5 - 1.96*se5 < 0.15 & theta5 + 1.96*se5 > 0.15,
    theta6 - 1.96*se6 < 0.15 & theta6 + 1.96*se6 > 0.15
  )
)

print(summary_df, row.names = FALSE)

cat("\nKey Point:\n")
cat("  • Approaches 4 & 6 AVERAGE multiple trees (more stable)\n")
cat("  • Approach 5 uses SINGLE modal tree (less stable)\n")
cat("  • Averaging reduces variance from tree structure selection\n")

cat("\n", rep("=", 80), "\n", sep = "")
