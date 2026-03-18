# verify_att_dgp.R
# Verify that the corrected DGPs have the right true ATT

suppressMessages({
  devtools::load_all()
  source("simulations/dgps_att_correct.R")
})

cat("Verifying corrected DGPs...\n\n")

# Test parameters
n_large <- 100000
tau_target <- 0.10

# Test each DGP
dgps <- list(
  binary = generate_dgp_binary_att,
  continuous = generate_dgp_continuous_att,
  moderate = generate_dgp_moderate_att
)

for (dgp_name in names(dgps)) {
  cat(sprintf("DGP: %s\n", dgp_name))
  cat(rep("-", 40), "\n", sep = "")

  # Generate large sample
  d <- dgps[[dgp_name]](n_large, tau = tau_target, seed = 123)

  # Compute true ATT
  treated_idx <- which(d$A == 1)
  true_att <- mean(d$true_m1[treated_idx] - d$true_m0[treated_idx])
  true_ate <- mean(d$true_m1 - d$true_m0)

  cat(sprintf("  Target tau:  %.4f\n", tau_target))
  cat(sprintf("  True ATT:    %.4f\n", true_att))
  cat(sprintf("  True ATE:    %.4f\n", true_ate))
  cat(sprintf("  Discrepancy: %.4f\n", abs(true_att - tau_target)))

  # Check propensity
  cat(sprintf("  e(X) range:  [%.3f, %.3f]\n", min(d$true_e), max(d$true_e)))

  # Success?
  if (abs(true_att - tau_target) < 0.01) {
    cat("  ✓ CORRECT\n")
  } else {
    cat("  ✗ INCORRECT\n")
  }
  cat("\n")
}

cat("If all DGPs show ✓ CORRECT, proceed to coverage testing.\n")
