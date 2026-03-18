#' Verification Script for Beta Continuous DGPs
#'
#' Tests that the new dgps_beta_continuous.R DGPs are correctly implemented.
#' Checks:
#' 1. X is truly continuous (not manually discretized)
#' 2. Nuisances vary smoothly
#' 3. ATT ≈ 0.10 across all regimes
#' 4. Theoretical s_n predictions are correct
#' 5. No extreme overlap violations

source("dgps/dgps_beta_continuous.R")

cat(strrep("=", 70), "\n")
cat("Verification: β Continuous DGPs\n")
cat(strrep("=", 70), "\n\n")

n_test <- 800
tau_test <- 0.10
seed_test <- 12345

# Generate all three DGPs
cat("Generating DGPs (n=800, tau=0.10)...\n\n")
d_high <- generate_dgp_beta_high(n = n_test, tau = tau_test, seed = seed_test)
d_boundary <- generate_dgp_beta_boundary(n = n_test, tau = tau_test, seed = seed_test)
d_low <- generate_dgp_beta_low(n = n_test, tau = tau_test, seed = seed_test)

dgps <- list(
  "beta_high (β=3)" = d_high,
  "beta_boundary (β=2)" = d_boundary,
  "beta_low (β=1)" = d_low
)

#------------------------------------------------------------------------------
# Check 1: X is continuous (not manually discretized)
#------------------------------------------------------------------------------

cat("Check 1: Features are continuous\n")
cat(strrep("-", 70), "\n")

for (dgp_name in names(dgps)) {
  d <- dgps[[dgp_name]]
  n_unique <- sapply(d$X, function(x) length(unique(x)))

  cat(sprintf("%s:\n", dgp_name))
  cat(sprintf("  X1: %d unique values\n", n_unique[1]))
  cat(sprintf("  X2: %d unique values\n", n_unique[2]))
  cat(sprintf("  X3: %d unique values\n", n_unique[3]))
  cat(sprintf("  X4: %d unique values\n", n_unique[4]))

  if (all(n_unique > 100)) {
    cat("  ✓ All features are continuous (>100 unique values)\n")
  } else {
    cat("  ✗ FAILED: Features appear discretized\n")
  }
  cat("\n")
}

#------------------------------------------------------------------------------
# Check 2: ATT is stable across regimes
#------------------------------------------------------------------------------

cat("\nCheck 2: ATT stability across regimes\n")
cat(strrep("-", 70), "\n")

for (dgp_name in names(dgps)) {
  d <- dgps[[dgp_name]]
  cat(sprintf("%s:\n", dgp_name))
  cat(sprintf("  Target tau: %.3f\n", tau_test))
  cat(sprintf("  True ATT:   %.3f\n", d$true_att))
  cat(sprintf("  Difference: %.4f\n", abs(d$true_att - tau_test)))

  if (abs(d$true_att - tau_test) < 0.02) {
    cat("  ✓ ATT close to target (within 0.02)\n")
  } else {
    cat("  ⚠️  ATT differs from target by >0.02\n")
  }
  cat("\n")
}

#------------------------------------------------------------------------------
# Check 3: Propensity overlap
#------------------------------------------------------------------------------

cat("\nCheck 3: Propensity overlap\n")
cat(strrep("-", 70), "\n")

for (dgp_name in names(dgps)) {
  d <- dgps[[dgp_name]]
  e_min <- min(d$true_e)
  e_max <- max(d$true_e)
  e_mean <- mean(d$true_e)

  cat(sprintf("%s:\n", dgp_name))
  cat(sprintf("  e(X) range: [%.3f, %.3f]\n", e_min, e_max))
  cat(sprintf("  e(X) mean:  %.3f\n", e_mean))

  if (e_min >= 0.05 && e_max <= 0.95) {
    cat("  ✓ Good overlap (e ∈ [0.05, 0.95])\n")
  } else {
    cat("  ⚠️  Potential overlap issues\n")
  }
  cat("\n")
}

#------------------------------------------------------------------------------
# Check 4: Theoretical s_n predictions
#------------------------------------------------------------------------------

cat("\nCheck 4: Theoretical s_n predictions (n=800)\n")
cat(strrep("-", 70), "\n")

for (dgp_name in names(dgps)) {
  d <- dgps[[dgp_name]]
  diag <- d$diagnostics

  cat(sprintf("%s:\n", dgp_name))
  cat(sprintf("  β = %d, d = %d\n", diag$beta, diag$d))
  cat(sprintf("  Formula: s_n = n^(d/(2β+d)) = n^(%d/(%d*%d+%d)) = n^%.3f\n",
              diag$d, 2, diag$beta, diag$d,
              diag$d / (2*diag$beta + diag$d)))
  cat(sprintf("  Theoretical s_n at n=800: %.1f leaves\n", diag$theoretical_sn))
  cat(sprintf("  Convergence rate: n^%.3f\n", diag$rate_exponent))
  cat(sprintf("  DML threshold: n^%.3f\n", diag$dml_threshold))
  cat(sprintf("  Condition: %s\n", diag$description))

  # Manual verification
  expected_sn <- n_test^(diag$d / (2*diag$beta + diag$d))
  if (abs(diag$theoretical_sn - expected_sn) < 0.01) {
    cat("  ✓ Theoretical s_n formula correct\n")
  } else {
    cat(sprintf("  ✗ FAILED: Expected %.1f, got %.1f\n", expected_sn, diag$theoretical_sn))
  }
  cat("\n")
}

#------------------------------------------------------------------------------
# Check 5: Nuisance smoothness (visual check)
#------------------------------------------------------------------------------

cat("\nCheck 5: Nuisance smoothness characterization\n")
cat(strrep("-", 70), "\n")

for (dgp_name in names(dgps)) {
  d <- dgps[[dgp_name]]
  diag <- d$diagnostics

  cat(sprintf("%s:\n", dgp_name))

  # For beta_high: should see cubic variation
  # For beta_boundary: should see quadratic variation
  # For beta_low: should see kinks (non-differentiable)

  # Sample 10 observations to check nuisance values
  idx <- seq(1, n_test, length.out = 10)
  cat(sprintf("  Sample e(X) values (10 obs): %.3f, %.3f, %.3f, ..., %.3f\n",
              d$true_e[idx[1]], d$true_e[idx[2]], d$true_e[idx[3]], d$true_e[idx[10]]))
  cat(sprintf("  Sample m0(X) values: %.3f, %.3f, %.3f, ..., %.3f\n",
              d$true_m0[idx[1]], d$true_m0[idx[2]], d$true_m0[idx[3]], d$true_m0[idx[10]]))

  # Check variation
  e_sd <- sd(d$true_e)
  m0_sd <- sd(d$true_m0)
  cat(sprintf("  SD(e):  %.3f\n", e_sd))
  cat(sprintf("  SD(m0): %.3f\n", m0_sd))

  if (e_sd > 0.05 && m0_sd > 0.05) {
    cat("  ✓ Sufficient variation in nuisances\n")
  } else {
    cat("  ⚠️  Nuisances may be too flat\n")
  }
  cat("\n")
}

#------------------------------------------------------------------------------
# Summary
#------------------------------------------------------------------------------

cat(strrep("=", 70), "\n")
cat("Verification Summary\n")
cat(strrep("=", 70), "\n\n")

cat("All checks complete. If all ✓, proceed with full simulation.\n")
cat("If any ✗ or ⚠️, review DGP implementation before running.\n\n")

cat("Expected theoretical s_n at n=800:\n")
cat(sprintf("  β=3 (high):     ~%.1f leaves (small trees, fast convergence)\n",
            d_high$diagnostics$theoretical_sn))
cat(sprintf("  β=2 (boundary): ~%.1f leaves (medium trees)\n",
            d_boundary$diagnostics$theoretical_sn))
cat(sprintf("  β=1 (low):      ~%.1f leaves (large trees, slow convergence)\n",
            d_low$diagnostics$theoretical_sn))

cat("\nNext steps:\n")
cat("  1. If verification passes, run: Rscript run_beta_study.R\n")
cat("  2. Monitor fitted s_n in results (should match theoretical predictions)\n")
cat("  3. If fitted s_n << theoretical s_n, regularization is too strong\n\n")
