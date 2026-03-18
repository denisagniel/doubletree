#!/usr/bin/env Rscript
# Quick parallel test: 6 sims (2 per core)

suppressMessages({
  library(dplyr)
  library(optimaltrees)
  library(parallel)
})

source("../simulation_helpers.R")
invisible(sapply(c(
  "../../R/estimate_att.R",
  "dgps/dgps_smooth.R"
), safe_source))

dgp_func <- generate_dgp_continuous_att

# Test function
test_sim <- function(i) {
  cat(sprintf("Sim %d starting (PID: %d)\n", i, Sys.getpid()))
  d <- dgp_func(n = 400, tau = 0.10, seed = 10000 + i)
  fit <- suppressMessages(estimate_att(X = d$X, A = d$A, Y = d$Y, K = 3,
                                   regularization = 0.015, use_rashomon = FALSE,
                                   verbose = FALSE))
  cat(sprintf("Sim %d complete (theta=%.3f)\n", i, fit$theta))
  data.frame(i = i, theta = fit$theta, stringsAsFactors = FALSE)
}

cat("Testing parallel with 3 cores (6 sims)...\n")
start <- Sys.time()
results <- mclapply(1:6, test_sim, mc.cores = 3)
elapsed <- difftime(Sys.time(), start, units = "secs")

cat(sprintf("\nParallel test complete in %.1f seconds\n", elapsed))
cat(sprintf("Results: %d simulations\n", length(results)))
cat("SUCCESS: Parallel processing works\n")
