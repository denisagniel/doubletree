#!/usr/bin/env Rscript

#' Test Phase 1 DGPs (DGP4-6) locally before O2 deployment
#'
#' Validates:
#' - DGP4: Continuous features, binary outcome
#' - DGP5: Continuous features, continuous outcome
#' - DGP6: Mixed features, binary outcome

suppressMessages({
  library(optimaltrees)
  library(doubletree)
})

source("dgps/dgps_continuous.R")

cat("Testing Phase 1 DGPs (DGP4-6)\n")
cat("=============================\n\n")

n <- 400
seed <- 10001
lambda <- log(n) / n

# Test DGP4: Continuous features, binary outcome
cat("DGP4: Continuous features, binary outcome\n")
cat(strrep("-", 40), "\n", sep = "")
d4 <- generate_dgp_continuous_binary(n = n, tau = 0.10, seed = seed)
cat(sprintf("  n = %d, true ATT = %.4f\n", nrow(d4$X), d4$true_att))
cat(sprintf("  Feature types: %s\n", paste(sapply(d4$X, class), collapse = ", ")))
cat(sprintf("  Outcome: binary, mean = %.3f\n", mean(d4$Y)))
cat(sprintf("  Treatment rate: %.1f%%\n", 100*mean(d4$A)))

cat("  Testing tree method... ")
start <- Sys.time()
fit4 <- tryCatch({
  estimate_att(
    X = d4$X, A = d4$A, Y = d4$Y,
    K = 5, outcome_type = "binary",
    regularization = lambda,
    cv_regularization = FALSE,
    use_rashomon = FALSE,
    verbose = FALSE
  )
}, error = function(e) {
  cat(sprintf("\n    ERROR: %s\n", conditionMessage(e)))
  NULL
})
time4 <- difftime(Sys.time(), start, units = "secs")

if (!is.null(fit4)) {
  cat(sprintf("SUCCESS (%.2fs)\n", time4))
  cat(sprintf("    θ = %.4f, σ = %.4f, bias = %.4f\n",
              fit4$theta, fit4$sigma, fit4$theta - d4$true_att))
} else {
  cat("FAILED\n")
}

# Test DGP5: Continuous features, continuous outcome
cat("\n\nDGP5: Continuous features, continuous outcome\n")
cat(strrep("-", 40), "\n", sep = "")
d5 <- generate_dgp_continuous_continuous(n = n, tau = 0.10, seed = seed)
cat(sprintf("  n = %d, true ATT = %.4f\n", nrow(d5$X), d5$true_att))
cat(sprintf("  Feature types: %s\n", paste(sapply(d5$X, class), collapse = ", ")))
cat(sprintf("  Outcome: continuous, mean = %.3f, sd = %.3f\n",
            mean(d5$Y), sd(d5$Y)))
cat(sprintf("  Treatment rate: %.1f%%\n", 100*mean(d5$A)))

cat("  Testing tree method... ")
start <- Sys.time()
fit5 <- tryCatch({
  estimate_att(
    X = d5$X, A = d5$A, Y = d5$Y,
    K = 5, outcome_type = "continuous",
    regularization = lambda,
    cv_regularization = FALSE,
    use_rashomon = FALSE,
    verbose = FALSE
  )
}, error = function(e) {
  cat(sprintf("\n    ERROR: %s\n", conditionMessage(e)))
  NULL
})
time5 <- difftime(Sys.time(), start, units = "secs")

if (!is.null(fit5)) {
  cat(sprintf("SUCCESS (%.2fs)\n", time5))
  cat(sprintf("    θ = %.4f, σ = %.4f, bias = %.4f\n",
              fit5$theta, fit5$sigma, fit5$theta - d5$true_att))
} else {
  cat("FAILED\n")
}

# Test DGP6: Mixed features, binary outcome
cat("\n\nDGP6: Mixed features (2 binary + 2 continuous), binary outcome\n")
cat(strrep("-", 40), "\n", sep = "")
d6 <- generate_dgp_mixed(n = n, tau = 0.10, seed = seed)
cat(sprintf("  n = %d, true ATT = %.4f\n", nrow(d6$X), d6$true_att))
cat(sprintf("  Feature types: %s\n", paste(sapply(d6$X, class), collapse = ", ")))
cat(sprintf("  Outcome: binary, mean = %.3f\n", mean(d6$Y)))
cat(sprintf("  Treatment rate: %.1f%%\n", 100*mean(d6$A)))

cat("  Testing tree method... ")
start <- Sys.time()
fit6 <- tryCatch({
  estimate_att(
    X = d6$X, A = d6$A, Y = d6$Y,
    K = 5, outcome_type = "binary",
    regularization = lambda,
    cv_regularization = FALSE,
    use_rashomon = FALSE,
    verbose = FALSE
  )
}, error = function(e) {
  cat(sprintf("\n    ERROR: %s\n", conditionMessage(e)))
  NULL
})
time6 <- difftime(Sys.time(), start, units = "secs")

if (!is.null(fit6)) {
  cat(sprintf("SUCCESS (%.2fs)\n", time6))
  cat(sprintf("    θ = %.4f, σ = %.4f, bias = %.4f\n",
              fit6$theta, fit6$sigma, fit6$theta - d6$true_att))
} else {
  cat("FAILED\n")
}

# Summary
cat("\n\n")
cat(strrep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n", sep = "")

success_count <- sum(c(
  !is.null(fit4),
  !is.null(fit5),
  !is.null(fit6)
))

if (success_count == 3) {
  cat("✅ All 3 DGPs working correctly!\n")
  cat("\nReady for O2 deployment:\n")
  cat("  - DGP4: Continuous features ✓\n")
  cat("  - DGP5: Continuous outcomes ✓\n")
  cat("  - DGP6: Mixed features ✓\n")
  cat("\nNext steps:\n")
  cat("  1. Update launch_batched_simulations.sh to include dgp4-6\n")
  cat("  2. Deploy to O2 (36,000 new replications)\n")
} else {
  cat(sprintf("⚠️  Only %d/3 DGPs working\n", success_count))
  cat("Fix issues before O2 deployment\n")
}
