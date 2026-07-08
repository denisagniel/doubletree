# Show M-Split Approaches (5 and 6) in Detail

library(doubletree)
library(optimaltrees)

source("code/dgps.R")

set.seed(777)

# Generate data
dgp <- generate_dgp_simple(n = 500)
X <- dgp$X
A <- dgp$A
Y <- dgp$Y

cat("\n================================================================================\n")
cat("DATA\n")
cat("================================================================================\n\n")
cat("DGP: Simple (linear propensity, linear outcome)\n")
cat("n = 500 | Treated:", sum(A), "| Control:", sum(1-A), "\n")
cat("True ATT: 0.15\n")
cat("Sample ATT:", round(mean(Y[A==1]) - mean(Y[A==0]), 3), "\n\n")

# ==================================================================
# APPROACH 5: M-SPLIT (Single Modal Tree)
# ==================================================================

cat("================================================================================\n")
cat("APPROACH 5: M-SPLIT (Modal Structure, Single Tree)\n")
cat("================================================================================\n\n")

cat("HOW IT WORKS:\n")
cat("  Step 1: Repeat M=10 times ->\n")
cat("            For each repetition, randomly split data into K=5 folds\n")
cat("            Fit propensity & outcome trees on each fold\n")
cat("            Total: M×K = 50 trees fitted\n\n")
cat("  Step 2: Vote across all 50 trees -> Find MODAL (most common) structure\n\n")
cat("  Step 3: Using the modal structure:\n")
cat("            Refit leaf values via K-fold cross-fitting\n")
cat("            Compute ATT using SINGLE tree with modal structure\n\n")

result5 <- doubletree::estimate_att_msplit(
  X, A, Y,
  M = 10,
  K = 5,
  regularization = 0.1,  # Fixed lambda for approach 5
  verbose = FALSE
)

cat("OUTCOME:\n")
cat("  Total trees fitted:", result5$M * result5$K, "\n")
cat("  Modal structure found:", !is.null(result5$structures$modal), "\n")

if (!is.null(result5$structures$modal)) {
  modal <- result5$structures$modal
  freq <- result5$structure_selection$modal_freq
  total <- result5$M * result5$K
  pct <- round(100 * freq / total, 1)

  cat("  Modal structure appeared in:", freq, "/", total, "trees (", pct, "%)\n\n", sep = "")

  cat("MODAL TREE STRUCTURE:\n")
  cat("  Number of leaves:", length(modal$e$leaf_paths), "\n\n")

  cat("  Propensity splits (", nrow(modal$e$splits), "):\n", sep = "")
  print(modal$e$splits, row.names = FALSE)
  cat("\n")

  cat("  Outcome splits (", nrow(modal$m0$splits), "):\n", sep = "")
  print(modal$m0$splits, row.names = FALSE)
  cat("\n")

  cat("INFERENCE:\n")
  cat("  Uses SINGLE tree with modal structure\n")
  cat("  Leaf values fitted via cross-fitting on full sample\n")
  cat("  ATT computed from this single tree\n\n")
}

theta5 <- result5$theta
se5 <- result5$sigma
ci5_lower <- theta5 - 1.96 * se5
ci5_upper <- theta5 + 1.96 * se5

cat("RESULT:\n")
cat("  Estimate: ", round(theta5, 4), "\n", sep = "")
cat("  SE:       ", round(se5, 4), "\n", sep = "")
cat("  95% CI:   [", round(ci5_lower, 4), ", ", round(ci5_upper, 4), "]\n", sep = "")
cat("  Covers truth (0.15):", ci5_lower < 0.15 & ci5_upper > 0.15, "\n\n")

# ==================================================================
# APPROACH 6: M-SPLIT AVERAGED (Average Across Modal Trees)
# ==================================================================

cat("================================================================================\n")
cat("APPROACH 6: M-SPLIT AVERAGED (Modal Structure, Averaged Trees)\n")
cat("================================================================================\n\n")

cat("HOW IT WORKS:\n")
cat("  Step 1: Repeat M=10 times ->\n")
cat("            For each repetition, randomly split data into K=5 folds\n")
cat("            Fit trees with CV-SELECTED lambda on each fold\n")
cat("            Total: M×K = 50 trees fitted (each with optimal lambda)\n\n")
cat("  Step 2: Vote across all 50 trees -> Find MODAL structure\n\n")
cat("  Step 3: AVERAGE across ALL trees matching modal structure:\n")
cat("            Filter 50 trees to those with modal structure\n")
cat("            Average leaf values across ALL matching trees\n")
cat("            Each tree contributes equally to the average\n")
cat("            Compute ATT using AVERAGED tree\n\n")

result6 <- doubletree::estimate_att_msplit_averaged(
  X, A, Y,
  M = 10,
  K = 5,
  verbose = FALSE
)

cat("OUTCOME:\n")
cat("  Total trees fitted:", result6$M * result6$K, "\n")
cat("  Modal structure found:", !is.null(result6$structures$modal), "\n")

if (!is.null(result6$structures$modal)) {
  modal <- result6$structures$modal
  freq <- result6$structure_selection$modal_freq
  n_avg <- result6$diagnostics$n_trees_averaged
  total <- result6$M * result6$K
  freq_pct <- round(100 * freq / total, 1)

  cat("  Modal structure appeared in:", freq, "/", total, "trees (", freq_pct, "%)\n", sep = "")
  cat("  Trees averaged:", n_avg, "\n\n")

  cat("MODAL TREE STRUCTURE:\n")
  cat("  Number of leaves:", length(modal$e$leaf_paths), "\n\n")

  cat("  Propensity splits (", nrow(modal$e$splits), "):\n", sep = "")
  print(modal$e$splits, row.names = FALSE)
  cat("\n")

  cat("  Outcome splits (", nrow(modal$m0$splits), "):\n", sep = "")
  print(modal$m0$splits, row.names = FALSE)
  cat("\n")

  cat("INFERENCE:\n")
  cat("  Uses AVERAGED tree (", n_avg, " trees combined)\n", sep = "")
  cat("  Each of the ", n_avg, " trees contributes 1/", n_avg, " = ",
      round(100/n_avg, 1), "% to leaf values\n", sep = "")
  cat("  Averaging reduces variance from structure selection\n")
  cat("  ATT computed from averaged tree\n\n")
}

theta6 <- result6$theta
se6 <- result6$sigma
ci6_lower <- theta6 - 1.96 * se6
ci6_upper <- theta6 + 1.96 * se6

cat("RESULT:\n")
cat("  Estimate: ", round(theta6, 4), "\n", sep = "")
cat("  SE:       ", round(se6, 4), "\n", sep = "")
cat("  95% CI:   [", round(ci6_lower, 4), ", ", round(ci6_upper, 4), "]\n", sep = "")
cat("  Covers truth (0.15):", ci6_lower < 0.15 & ci6_upper > 0.15, "\n\n")

# ==================================================================
# COMPARISON
# ==================================================================

cat("================================================================================\n")
cat("KEY DIFFERENCES: APPROACH 5 vs 6\n")
cat("================================================================================\n\n")

cat("APPROACH 5 (M-split):\n")
cat("  • Fixed lambda = 0.1 for all trees\n")
cat("  • Voting across 50 trees -> modal structure\n")
cat("  • REFIT leaf values using ONE cross-fitted tree\n")
cat("  • Uses SINGLE tree for inference\n")
cat("  • More variable (depends on one structure fit)\n\n")

cat("APPROACH 6 (M-split Averaged):\n")
cat("  • CV-selected lambda for each tree (adaptive)\n")
cat("  • Voting across 50 trees -> modal structure\n")
cat("  • AVERAGE leaf values across ALL trees with modal structure\n")
cat("  • Uses AVERAGED tree for inference\n")
cat("  • More stable (averages out noise)\n\n")

cat("COMPARISON:\n")
comp_df <- data.frame(
  Metric = c("Estimate", "SE", "CI Width", "Covers Truth"),
  Approach_5 = c(
    round(theta5, 4),
    round(se5, 4),
    round(ci5_upper - ci5_lower, 4),
    ci5_lower < 0.15 & ci5_upper > 0.15
  ),
  Approach_6 = c(
    round(theta6, 4),
    round(se6, 4),
    round(ci6_upper - ci6_lower, 4),
    ci6_lower < 0.15 & ci6_upper > 0.15
  )
)
print(comp_df, row.names = FALSE)

cat("\n✓ Approach 6 typically has similar or smaller SE due to averaging\n")
cat("✓ Both use same modal structure but differ in how leaf values are computed\n")
cat("✓ Averaging (Approach 6) reduces sensitivity to individual tree fits\n\n")

cat("================================================================================\n")
