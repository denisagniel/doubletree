#!/usr/bin/env Rscript

# Master Script: Run All Diagnostic Scripts
# Created: 2026-05-27
#
# Executes all diagnostic scripts in sequence and generates summary report

# ============================================================================
# Setup
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("  Complex DGP Diagnostic Suite\n")
cat("  Systematic analysis of tree-based nuisance estimation\n")
cat("================================================================\n\n")

start_time <- Sys.time()

# Track which diagnostics succeeded
diagnostic_status <- list()

# Output directory for summary
summary_dir <- "diagnostics/results/summary"
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE)
}

# ============================================================================
# Phase 1.1: Propensity Score Diagnostics
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("Phase 1.1: Propensity Score Diagnostics\n")
cat("================================================================\n\n")

ps_status <- tryCatch({
  source("diagnostics/01_propensity_diagnostics.R", echo = FALSE)
  cat("\n✓ Propensity diagnostics complete\n")
  "SUCCESS"
}, error = function(e) {
  cat(sprintf("\n✗ ERROR in propensity diagnostics:\n%s\n", conditionMessage(e)))
  "FAILED"
})

diagnostic_status$propensity <- ps_status

# ============================================================================
# Phase 1.2: Outcome Model Diagnostics
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("Phase 1.2: Outcome Model Diagnostics\n")
cat("================================================================\n\n")

outcome_status <- tryCatch({
  source("diagnostics/02_outcome_diagnostics.R", echo = FALSE)
  cat("\n✓ Outcome diagnostics complete\n")
  "SUCCESS"
}, error = function(e) {
  cat(sprintf("\n✗ ERROR in outcome diagnostics:\n%s\n", conditionMessage(e)))
  "FAILED"
})

diagnostic_status$outcome <- outcome_status

# ============================================================================
# Phase 5: EIF Decomposition
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("Phase 5: EIF Component Decomposition\n")
cat("================================================================\n\n")

eif_status <- tryCatch({
  source("diagnostics/05_eif_decomposition.R", echo = FALSE)
  cat("\n✓ EIF decomposition complete\n")
  "SUCCESS"
}, error = function(e) {
  cat(sprintf("\n✗ ERROR in EIF decomposition:\n%s\n", conditionMessage(e)))
  "FAILED"
})

diagnostic_status$eif <- eif_status

# ============================================================================
# Generate Summary Report
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("Generating Summary Report\n")
cat("================================================================\n\n")

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

# Load all results
results_summary <- list()

if (diagnostic_status$propensity == "SUCCESS") {
  results_summary$propensity <- tryCatch({
    readRDS("diagnostics/results/propensity/propensity_diagnostics.rds")
  }, error = function(e) NULL)
}

if (diagnostic_status$outcome == "SUCCESS") {
  results_summary$outcome <- tryCatch({
    readRDS("diagnostics/results/outcome/outcome_diagnostics.rds")
  }, error = function(e) NULL)
}

if (diagnostic_status$eif == "SUCCESS") {
  results_summary$eif <- tryCatch({
    readRDS("diagnostics/results/eif_decomposition/eif_decomposition.rds")
  }, error = function(e) NULL)
}

# Create summary report
report_file <- file.path(summary_dir, "diagnostic_summary.txt")
sink(report_file)

cat("================================================================\n")
cat("  DIAGNOSTIC SUMMARY REPORT\n")
cat("================================================================\n\n")

cat(sprintf("Run date: %s\n", Sys.Date()))
cat(sprintf("Total time: %.1f minutes\n\n", as.numeric(elapsed)))

# Status
cat("Diagnostic Status:\n")
cat(sprintf("  Propensity score: %s\n", diagnostic_status$propensity))
cat(sprintf("  Outcome model: %s\n", diagnostic_status$outcome))
cat(sprintf("  EIF decomposition: %s\n", diagnostic_status$eif))

n_success <- sum(sapply(diagnostic_status, function(x) x == "SUCCESS"))
n_total <- length(diagnostic_status)

cat(sprintf("\nOverall: %d / %d diagnostics completed successfully\n\n", n_success, n_total))

# ============================================================================
# Key Findings
# ============================================================================

cat("================================================================\n")
cat("KEY FINDINGS\n")
cat("================================================================\n\n")

# Focus on complex DGP (dgp == 3)
if (!is.null(results_summary$propensity)) {
  ps_complex <- results_summary$propensity[results_summary$propensity$dgp == 3, ]

  cat("PROPENSITY SCORE (Complex DGP):\n")
  cat(sprintf("  Mean RMSE: %.4f\n", mean(ps_complex$rmse_e)))
  cat(sprintf("  Mean calibration slope: %.3f (should be ~1.0)\n",
              mean(ps_complex$calibration_slope, na.rm = TRUE)))
  cat(sprintf("  Mean tree size: %.1f leaves\n", mean(ps_complex$n_leaves)))
  cat(sprintf("  Extreme weights: %.1f%%\n", 100 * mean(ps_complex$extreme_total)))

  if (mean(ps_complex$rmse_e) > 0.10) {
    cat("  → ⚠ HIGH RMSE: Propensity trees struggling\n")
  } else if (mean(ps_complex$rmse_e) > 0.05) {
    cat("  → ⚠ MODERATE RMSE: Some PS estimation error\n")
  } else {
    cat("  → ✓ LOW RMSE: PS estimation adequate\n")
  }

  cat("\n")
}

if (!is.null(results_summary$outcome)) {
  outcome_complex <- results_summary$outcome[results_summary$outcome$dgp == 3, ]

  cat("OUTCOME MODEL (Complex DGP):\n")
  cat(sprintf("  Mean RMSE (control): %.4f\n", mean(outcome_complex$rmse_control)))
  cat(sprintf("  Mean RMSE (treated): %.4f\n", mean(outcome_complex$rmse_treated, na.rm = TRUE)))
  cat(sprintf("  Mean oracle RMSE: %.4f\n", mean(outcome_complex$oracle_rmse)))
  cat(sprintf("  Mean tree size: %.1f leaves\n", mean(outcome_complex$n_leaves)))

  oracle_rmse <- mean(outcome_complex$oracle_rmse)
  if (oracle_rmse > 0.10) {
    cat("  → ⚠ HIGH ORACLE RMSE: Trees cannot represent function\n")
    cat("     (Even with infinite data, trees struggle)\n")
  } else if (oracle_rmse > 0.05) {
    cat("  → ⚠ MODERATE ORACLE RMSE: Some expressiveness issues\n")
  } else {
    cat("  → ✓ LOW ORACLE RMSE: Trees can represent function\n")
  }

  extrapolation_ratio <- mean(outcome_complex$extrapolation_rmse_ratio, na.rm = TRUE)
  if (extrapolation_ratio > 1.5) {
    cat(sprintf("  → ⚠ EXTRAPOLATION ISSUE: RMSE %.1fx higher on treated units\n",
                extrapolation_ratio))
  }

  cat("\n")
}

if (!is.null(results_summary$eif)) {
  eif_complex <- results_summary$eif[results_summary$eif$dgp == 3, ]

  cat("EIF BIAS DECOMPOSITION (Complex DGP):\n")
  cat(sprintf("  Total bias: %.4f\n", mean(eif_complex$theta_bias)))
  cat(sprintf("  Component 1 (Outcome): %.4f\n", mean(eif_complex$comp1_bias)))
  cat(sprintf("  Component 2 (Propensity): %.4f\n", mean(eif_complex$comp2_bias)))

  comp1_abs <- mean(abs(eif_complex$comp1_bias))
  comp2_abs <- mean(abs(eif_complex$comp2_bias))

  if (comp1_abs > 2 * comp2_abs) {
    cat("\n  → PRIMARY ISSUE: Outcome model (Component 1)\n")
    cat("     Outcome trees are the main source of bias\n")
  } else if (comp2_abs > 2 * comp1_abs) {
    cat("\n  → PRIMARY ISSUE: Propensity score (Component 2)\n")
    cat("     Propensity trees are the main source of bias\n")
  } else {
    cat("\n  → BOTH components contribute to bias\n")
    cat("     Both nuisance functions need improvement\n")
  }

  cat("\n")
}

# ============================================================================
# Recommendations
# ============================================================================

cat("================================================================\n")
cat("RECOMMENDATIONS\n")
cat("================================================================\n\n")

# Automated recommendations based on findings
if (!is.null(results_summary$outcome)) {
  outcome_complex <- results_summary$outcome[results_summary$outcome$dgp == 3, ]
  oracle_rmse <- mean(outcome_complex$oracle_rmse)

  if (oracle_rmse > 0.10) {
    cat("1. HIGH PRIORITY: Trees cannot represent functions\n")
    cat("   - Consider ensemble methods (random forests)\n")
    cat("   - Or test much weaker regularization (allow deeper trees)\n")
    cat("   - Or increase sample size substantially\n\n")
  } else if (oracle_rmse > 0.05) {
    cat("1. MODERATE PRIORITY: Tree expressiveness\n")
    cat("   - Test weaker regularization (λ × 0.25 or λ × 0.5)\n")
    cat("   - Allow trees to grow deeper\n\n")
  } else {
    cat("1. Trees CAN represent functions (oracle RMSE low)\n")
    cat("   - Issue is estimation noise, not expressiveness\n")
    cat("   - Focus on regularization calibration\n\n")
  }
}

if (!is.null(results_summary$eif)) {
  eif_complex <- results_summary$eif[results_summary$eif$dgp == 3, ]
  comp1_abs <- mean(abs(eif_complex$comp1_bias))
  comp2_abs <- mean(abs(eif_complex$comp2_bias))

  if (comp1_abs > 1.5 * comp2_abs) {
    cat("2. Focus on OUTCOME MODEL first\n")
    cat("   - Test weaker regularization for outcome trees\n")
    cat("   - Investigate extrapolation to treated units\n\n")
  } else if (comp2_abs > 1.5 * comp1_abs) {
    cat("2. Focus on PROPENSITY SCORE first\n")
    cat("   - Test weaker regularization for PS trees\n")
    cat("   - Check overlap quality\n\n")
  } else {
    cat("2. Improve BOTH nuisance functions\n")
    cat("   - Test weaker regularization globally\n")
    cat("   - Or use CV to select lambda\n\n")
  }
}

cat("3. Next steps:\n")
cat("   - Review plots in diagnostics/results/*/\n")
cat("   - Run calibration experiments (06_calibration_experiments.R)\n")
cat("   - Test interventions based on findings above\n\n")

# ============================================================================
# File Locations
# ============================================================================

cat("================================================================\n")
cat("OUTPUT FILES\n")
cat("================================================================\n\n")

cat("Results:\n")
if (diagnostic_status$propensity == "SUCCESS") {
  cat("  - diagnostics/results/propensity/propensity_diagnostics.rds\n")
  cat("  - diagnostics/results/propensity/*.png\n")
}
if (diagnostic_status$outcome == "SUCCESS") {
  cat("  - diagnostics/results/outcome/outcome_diagnostics.rds\n")
  cat("  - diagnostics/results/outcome/*.png\n")
}
if (diagnostic_status$eif == "SUCCESS") {
  cat("  - diagnostics/results/eif_decomposition/eif_decomposition.rds\n")
  cat("  - diagnostics/results/eif_decomposition/*.png\n")
}

cat(sprintf("\nThis report: %s\n", report_file))

sink()

# ============================================================================
# Print Summary to Console
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("ALL DIAGNOSTICS COMPLETE\n")
cat("================================================================\n\n")

cat(sprintf("Total time: %.1f minutes\n\n", as.numeric(elapsed)))

cat("Status:\n")
for (name in names(diagnostic_status)) {
  status <- diagnostic_status[[name]]
  symbol <- ifelse(status == "SUCCESS", "✓", "✗")
  cat(sprintf("  %s %s\n", symbol, name))
}

cat(sprintf("\nSummary report saved to:\n  %s\n\n", report_file))

cat("Review:\n")
cat("  1. Summary report (above file)\n")
cat("  2. Plots in diagnostics/results/*/\n")
cat("  3. Full results in *.rds files\n\n")

cat("Next steps:\n")
cat("  - Review findings in summary report\n")
cat("  - Examine plots for visual confirmation\n")
cat("  - Proceed to calibration experiments\n\n")

# Exit with appropriate code
if (n_success == n_total) {
  cat("✓ All diagnostics successful\n\n")
  quit(status = 0)
} else {
  cat(sprintf("⚠ %d / %d diagnostics failed\n\n", n_total - n_success, n_total))
  quit(status = 1)
}
