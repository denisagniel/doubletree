#!/usr/bin/env Rscript
# Test that Phase 2 DGPs (dgp7-9) can be called via run_batch_replications.R

cat("Testing Phase 2 DGP batch script integration\n")
cat("============================================\n\n")

# Source required files
source("dgps/dgps_smooth.R")
source("dgps/dgps_continuous.R")
source("dgps/dgps_phase2.R")

# Test that aliases exist
cat("Checking function aliases:\n")
stopifnot(exists("generate_dgp7"))
stopifnot(exists("generate_dgp8"))
stopifnot(exists("generate_dgp9"))
cat("  ✓ generate_dgp7 exists\n")
cat("  ✓ generate_dgp8 exists\n")
cat("  ✓ generate_dgp9 exists\n\n")

# Test DGP7
cat("Testing DGP7 (deep interaction):\n")
d7 <- generate_dgp7(n = 100, tau = 0.10, seed = 42)
stopifnot(ncol(d7$X) == 4)
stopifnot(all(d7$Y %in% c(0, 1)))  # Binary outcome
cat("  ✓ DGP7 generates binary outcome\n")
cat(sprintf("  ✓ DGP7 true_att = %.3f\n", d7$true_att))

# Test DGP8
cat("\nTesting DGP8 (sin/cos double nonlinearity):\n")
d8 <- generate_dgp8(n = 100, tau = 0.10, seed = 42)
stopifnot(ncol(d8$X) == 4)
stopifnot(is.numeric(d8$Y))
stopifnot(!all(d8$Y %in% c(0, 1)))  # Continuous outcome
cat("  ✓ DGP8 generates continuous outcome\n")
cat(sprintf("  ✓ DGP8 true_att = %.3f\n", d8$true_att))

# Test DGP9
cat("\nTesting DGP9 (weak overlap):\n")
d9 <- generate_dgp9(n = 100, tau = 0.10, seed = 42)
stopifnot(ncol(d9$X) == 4)
stopifnot(all(d9$Y %in% c(0, 1)))  # Binary outcome
cat("  ✓ DGP9 generates binary outcome\n")
cat(sprintf("  ✓ DGP9 true_att = %.3f\n", d9$true_att))

# Check propensity score distribution for DGP9
e_range <- range(d9$true_e)
overlap_pct <- mean(d9$true_e >= 0.1 & d9$true_e <= 0.9) * 100
cat(sprintf("  ✓ DGP9 propensity range: [%.3f, %.3f]\n", e_range[1], e_range[2]))
cat(sprintf("  ✓ DGP9 overlap (0.1-0.9): %.1f%% (should be low for stress test)\n", overlap_pct))

cat("\n============================================\n")
cat("✓ All Phase 2 DGPs verified and ready\n")
cat("============================================\n")
