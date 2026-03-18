#!/usr/bin/env Rscript
# Final test: 10 simulations to confirm everything works

cat("Testing 10 simulations...\n\n")
flush.console()

library(dplyr)
library(optimaltrees)

source("../simulation_helpers.R")
invisible(sapply(c(
  "../../R/estimate_att.R",
  "../../R/nuisance_trees.R",
  "../../R/score_att.R",
  "../../R/inference.R",
  "../../R/utils.R",
  "dgps/dgps_smooth.R"
), safe_source))

# 10 simulations: dgp1, tree method only
run_sim <- function(i) {
  suppressWarnings(suppressMessages({
    d <- generate_dgp_binary_att(n = 400, tau = 0.1, seed = 1000 + i)
    fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = 2,
                   regularization = 0.01, use_rashomon = FALSE, verbose = FALSE)
  }))

  data.frame(sim = i, theta = fit$theta, sigma = fit$sigma,
             ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2])
}

cat("Running 10 simulations...\n")
start <- Sys.time()
results <- lapply(1:10, run_sim)
results_df <- dplyr::bind_rows(results)
elapsed <- difftime(Sys.time(), start, units = "secs")

cat(sprintf("\n✓ Complete in %.1f seconds\n", elapsed))
cat(sprintf("Mean theta: %.4f\n", mean(results_df$theta)))
cat("\nResults:\n")
print(results_df)

# Save
dir.create("results/test_10", showWarnings = FALSE, recursive = TRUE)
saveRDS(results_df, "results/test_10/results.rds")
cat("\nSaved to results/test_10/results.rds\n")
