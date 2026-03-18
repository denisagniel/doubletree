# check_dgp_att.R
# Check if the DGP's true ATT actually equals the specified tau

suppressMessages({
  devtools::load_all()
  source("simulations/dgps_realistic.R")
})

cat("Checking if true ATT = specified tau...\n\n")

# Generate large sample to estimate true ATT
n_large <- 100000
tau_specified <- 0.15
d <- generate_dgp_simple(n_large, tau = tau_specified, seed = 123)

cat(sprintf("DGP: simple (4 binary features)\n"))
cat(sprintf("n = %d (large sample for precise estimation)\n", n_large))
cat(sprintf("Specified tau: %.4f\n\n", tau_specified))

# Compute true ATT using potential outcomes
# ATT = E[Y1 - Y0 | A=1]
treated_indices <- which(d$A == 1)
n_treated <- length(treated_indices)

true_att <- mean(d$true_m1[treated_indices] - d$true_m0[treated_indices])
true_ate <- mean(d$true_m1 - d$true_m0)

cat("True quantities:\n")
cat(sprintf("  True ATT:  %.4f\n", true_att))
cat(sprintf("  True ATE:  %.4f\n", true_ate))
cat(sprintf("  Specified: %.4f\n", tau_specified))
cat(sprintf("\nDiscrepancy (ATT - specified): %.4f\n", true_att - tau_specified))

# Check propensity distribution
cat(sprintf("\nPropensity e(X):\n"))
cat(sprintf("  Min:  %.4f\n", min(d$true_e)))
cat(sprintf("  Mean: %.4f\n", mean(d$true_e)))
cat(sprintf("  Max:  %.4f\n", max(d$true_e)))

# Check if ATT ≈ specified tau
if (abs(true_att - tau_specified) < 0.01) {
  cat("\n✓ True ATT matches specified tau\n")
  cat("DGP is correctly specified for ATT estimation\n")
} else {
  cat("\n✗ True ATT DOES NOT match specified tau\n")
  cat("This explains the bias in estimation!\n")
  cat("\nPossible causes:\n")
  cat("  - DGP constructs ATE, not ATT\n")
  cat("  - Treatment assignment depends on covariates\n")
  cat("  - Need to adjust DGP or estimand\n")
}
