#!/usr/bin/env Rscript
# Minimal test to diagnose background execution issue

cat("=== Background Execution Test ===\n")
flush.console()

# Test 1: Can we load packages?
cat("Test 1: Loading packages...\n")
flush.console()
suppressMessages({
  library(dplyr)
  library(optimaltrees)
})
cat("  ✓ Packages loaded\n")
flush.console()

# Test 2: Can we source functions?
cat("\nTest 2: Loading functions...\n")
flush.console()
source("../simulation_helpers.R")
invisible(sapply(c(
  "../../R/estimate_att.R",
  "../../R/nuisance_trees.R",
  "../../R/score_att.R",
  "../../R/inference.R",
  "../../R/utils.R",
  "dgps/dgps_smooth.R"
), safe_source))
cat("  ✓ Functions loaded\n")
flush.console()

# Test 3: Can we run ONE simulation?
cat("\nTest 3: Running single simulation...\n")
flush.console()
set.seed(123)
d <- generate_dgp_binary_att(n = 400, tau = 0.1)
cat("  Data generated (n=400)\n")
flush.console()

fit <- estimate_att(
  X = d$X, A = d$A, Y = d$Y,
  K = 2,
  regularization = 0.01,
  use_rashomon = FALSE,
  verbose = FALSE
)
cat(sprintf("  ✓ Fit complete: theta = %.4f\n", fit$theta))
flush.console()

# Test 4: Can we run 5 simulations in a loop?
cat("\nTest 4: Running 5 simulations in loop...\n")
flush.console()
results <- list()
for (i in 1:5) {
  cat(sprintf("  Sim %d...", i))
  flush.console()

  set.seed(1000 + i)
  d <- generate_dgp_binary_att(n = 400, tau = 0.1)
  fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = 2,
                 regularization = 0.01, use_rashomon = FALSE, verbose = FALSE)
  results[[i]] <- data.frame(sim = i, theta = fit$theta)

  cat(sprintf(" theta = %.4f\n", fit$theta))
  flush.console()
}
cat("  ✓ All 5 sims complete\n")
flush.console()

# Test 5: Can we use lapply?
cat("\nTest 5: Running 5 simulations with lapply...\n")
flush.console()
run_one <- function(i) {
  cat(sprintf("    lapply sim %d\n", i))
  flush.console()

  set.seed(2000 + i)
  d <- generate_dgp_binary_att(n = 400, tau = 0.1)
  fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = 2,
                 regularization = 0.01, use_rashomon = FALSE, verbose = FALSE)
  data.frame(sim = i, theta = fit$theta)
}

results_lapply <- lapply(1:5, run_one)
results_df <- dplyr::bind_rows(results_lapply)
cat(sprintf("  ✓ lapply complete: %d results\n", nrow(results_df)))
flush.console()

# Test 6: Can we save results?
cat("\nTest 6: Saving results...\n")
flush.console()
dir.create("results/test_background", recursive = TRUE, showWarnings = FALSE)
saveRDS(results_df, "results/test_background/test_results.rds")
cat("  ✓ Results saved\n")
flush.console()

cat("\n=== ALL TESTS PASSED ===\n")
cat("Background execution should work!\n")
flush.console()
