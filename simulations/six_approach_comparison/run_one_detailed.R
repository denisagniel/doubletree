# Run One Detailed Replication - Show Trees Contributing to Each Estimate

library(doubletree)
library(optimaltrees)

source("code/dgps.R")

set.seed(999)

# Generate data
dgp <- generate_dgp_simple(n = 500)
X <- dgp$X
A <- dgp$A
Y <- dgp$Y

cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("DATA GENERATION\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("DGP: Simple (linear propensity, linear outcome)\n")
cat("n =", nrow(X), "\n")
cat("Treated:", sum(A), "| Control:", sum(1-A), "\n")
cat("True ATT: 0.15\n")
cat("Sample ATT:", round(mean(Y[A==1]) - mean(Y[A==0]), 3), "\n\n")

# Helper to show result
show_result <- function(name, result, show_trees = TRUE) {
  cat(rep("=", 80), "\n", sep = "")
  cat(name, "\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Estimate
  theta <- result@estimate
  se <- result@se
  ci_lower <- theta - 1.96 * se
  ci_upper <- theta + 1.96 * se

  cat("Estimate: ", round(theta, 4), "\n", sep = "")
  cat("SE:       ", round(se, 4), "\n", sep = "")
  cat("95% CI:   [", round(ci_lower, 4), ", ", round(ci_upper, 4), "]\n\n", sep = "")

  if (show_trees) {
    # Tree structure
    cat("Tree Structure:\n")
    cat("  Leaves:", result@n_leaves, "\n")

    if (!is.null(result@propensity_splits) && nrow(result@propensity_splits) > 0) {
      cat("  Propensity splits:", nrow(result@propensity_splits), "\n")
      print(result@propensity_splits)
      cat("\n")
    }

    if (!is.null(result@outcome_splits) && nrow(result@outcome_splits) > 0) {
      cat("  Outcome splits:", nrow(result@outcome_splits), "\n")
      print(result@outcome_splits)
      cat("\n")
    }
  }
}

# ==================================================================
# APPROACH 3: DOUBLETREE (baseline - single tree from intersection)
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 3: DOUBLETREE (Rashomon Intersection)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Method:\n")
cat("  1. For each of K=5 folds, fit trees with CV-selected lambda\n")
cat("  2. Collect Rashomon sets (trees within 5% of best loss)\n")
cat("  3. Find common structure via intersection\n")
cat("  4. Refit leaf values using cross-fitting\n")
cat("  5. Compute ATT with single common tree\n\n")

result3 <- doubletree::estimate_att(
  X, A, Y,
  K = 5,
  rashomon_tol = 0.05,
  cv_regularization = TRUE,
  use_rashomon = TRUE,
  verbose = TRUE
)

cat("\nRashomon Details:\n")
cat("  Tolerance: 0.05\n")
cat("  Trees collected:", result3@n_rashomon_trees, "\n")
cat("  Trees intersecting:", result3@n_intersecting, "\n")
cat("  Intersection succeeded:", result3@intersection_succeeded, "\n\n")

show_result("APPROACH 3 RESULT", result3)

# ==================================================================
# APPROACH 4: DOUBLETREE AVERAGED
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 4: DOUBLETREE AVERAGED (Tree Averaging)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Method:\n")
cat("  1. Find common structure (via Rashomon intersection or 5-tier fallback)\n")
cat("  2. Collect K=5 fold-specific trees matching this structure\n")
cat("  3. Average leaf values across K trees\n")
cat("  4. Compute ATT with averaged tree\n\n")

result4 <- doubletree::estimate_att_doubletree_averaged(
  X, A, Y,
  K = 5,
  rashomon_tol = 0.05,
  verbose = TRUE
)

cat("\nAveraging Details:\n")
cat("  Method:", result4@method, "\n")
cat("  Tier used:", result4@tier, "\n")
cat("  Trees averaged:", result4@n_trees_averaged, "\n")

if (result4@tier <= 2) {
  cat("  ✓ Rashomon intersection succeeded\n")
  cat("  ✓ Averaged", result4@n_trees_averaged, "fold-specific trees with common structure\n\n")
} else if (result4@tier == 5) {
  cat("  ✗ Rashomon intersection failed\n")
  cat("  → Fallback tier", result4@tier, ": Fold-specific trees\n")
  cat("  → Averaged", result4@n_trees_averaged, "trees (one per fold)\n\n")
}

show_result("APPROACH 4 RESULT", result4)

# ==================================================================
# APPROACH 5: M-SPLIT (modal structure, single tree)
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 5: M-SPLIT (Modal Structure)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Method:\n")
cat("  1. Repeat M=10 times: split data, fit K=5 trees per split\n")
cat("  2. Total M×K = 50 trees fitted\n")
cat("  3. Find most common structure (modal)\n")
cat("  4. Refit leaf values using cross-fitting\n")
cat("  5. Compute ATT with single modal tree\n\n")

result5 <- doubletree::estimate_att_msplit(
  X, A, Y,
  M = 10,
  K = 5,
  regularization = 0.1,  # Fixed lambda (no CV in approach 5)
  verbose = TRUE
)

cat("\nModal Structure Details:\n")
cat("  M (splits):", result5@M, "\n")
cat("  K (folds per split):", result5@K, "\n")
cat("  Total trees:", result5@M * result5@K, "\n")
cat("  Modal structure found:", result5@modal_structure_found, "\n")

if (result5@modal_structure_found) {
  freq_pct <- round(100 * result5@modal_structure_freq / (result5@M * result5@K), 1)
  cat("  Modal frequency:", result5@modal_structure_freq, "/", result5@M * result5@K,
      "(", freq_pct, "%)\n\n", sep = "")
}

show_result("APPROACH 5 RESULT", result5)

# ==================================================================
# APPROACH 6: M-SPLIT AVERAGED
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("APPROACH 6: M-SPLIT AVERAGED (Tree Averaging)\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Method:\n")
cat("  1. Repeat M=10 times: split data, fit K=5 trees with CV per split\n")
cat("  2. Total M×K = 50 trees fitted (each with CV-selected lambda)\n")
cat("  3. Find most common structure (modal)\n")
cat("  4. Filter to trees matching modal structure\n")
cat("  5. Average leaf values across ALL matching trees\n")
cat("  6. Compute ATT with averaged tree\n\n")

result6 <- doubletree::estimate_att_msplit_averaged(
  X, A, Y,
  M = 10,
  K = 5,
  verbose = TRUE
)

cat("\nAveraging Details:\n")
cat("  M (splits):", result6@M, "\n")
cat("  K (folds per split):", result6@K, "\n")
cat("  Total trees:", result6@M * result6@K, "\n")
cat("  Modal structure found:", result6@modal_structure_found, "\n")

if (result6@modal_structure_found) {
  freq_pct <- round(100 * result6@modal_structure_freq / (result6@M * result6@K), 1)
  cat("  Modal frequency:", result6@modal_structure_freq, "(", freq_pct, "%)\n", sep = "")
  cat("  Trees averaged:", result6@n_trees_averaged, "\n")
  cat("  ✓ Each tree contributes equally: 1/", result6@n_trees_averaged, "\n\n", sep = "")
}

show_result("APPROACH 6 RESULT", result6)

# ==================================================================
# SUMMARY
# ==================================================================

cat(rep("=", 80), "\n", sep = "")
cat("SUMMARY: APPROACHES 3-6 COMPARISON\n")
cat(rep("=", 80), "\n\n", sep = "")

summary_df <- data.frame(
  Approach = c("3. Doubletree", "4. DT Averaged", "5. M-split", "6. MS Averaged"),
  Estimate = c(result3@estimate, result4@estimate, result5@estimate, result6@estimate),
  SE = c(result3@se, result4@se, result5@se, result6@se),
  Trees_Contrib = c("1", paste0(result4@n_trees_averaged, " avg"), "1",
                    paste0(result6@n_trees_averaged, " avg")),
  Method = c("Intersection", "Intersection+Avg", "Modal", "Modal+Avg")
)

print(summary_df, row.names = FALSE)

cat("\nTrue ATT: 0.15\n")
cat("All estimates within 2 SE of truth:", all(abs(summary_df$Estimate - 0.15) < 2 * summary_df$SE), "\n")

cat("\nKey Differences:\n")
cat("  • Approach 3: Single tree from intersection\n")
cat("  • Approach 4: Averages K fold-specific trees (same structure)\n")
cat("  • Approach 5: Single modal tree from M×K trees\n")
cat("  • Approach 6: Averages ALL trees matching modal structure\n")

cat("\n", rep("=", 80), "\n", sep = "")
