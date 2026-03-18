#!/usr/bin/env Rscript
# Minimal Test: Single simulation for DGP2 and DGP3
# Purpose: Quick diagnostic to verify basic functionality

suppressMessages({
  library(dplyr)
  library(optimaltrees)
})

cat("Loading sources...\n")
source("../simulation_helpers.R")
invisible(sapply(c(
  "../../R/estimate_att.R",
  "../../R/dml_att_repeated.R",
  "../../R/nuisance_trees.R",
  "../../R/score_att.R",
  "../../R/inference.R",
  "../../R/utils.R",
  "dgps/dgps_smooth.R",
  "methods/method_forest_dml.R",
  "methods/method_linear_dml.R"
), safe_source))

cat("Sources loaded successfully.\n\n")

# Test DGP2
cat("=== Testing DGP2 (Continuous) ===\n")
set.seed(10001)
d2 <- generate_dgp_continuous_att(n = 400, tau = 0.10, seed = 10001)
cat(sprintf("Data generated: n=%d, treatment prop=%.2f\n", nrow(d2$X), mean(d2$A)))

cat("Running tree method...\n")
fit2 <- estimate_att(X = d2$X, A = d2$A, Y = d2$Y, K = 5,
                regularization = log(400) / 400,
                use_rashomon = FALSE, verbose = TRUE)
cat(sprintf("Result: theta=%.4f, sigma=%.4f\n", fit2$theta, fit2$sigma))

# Test DGP3
cat("\n=== Testing DGP3 (Moderate) ===\n")
set.seed(10002)
d3 <- generate_dgp_moderate_att(n = 400, tau = 0.10, seed = 10002)
cat(sprintf("Data generated: n=%d, treatment prop=%.2f\n", nrow(d3$X), mean(d3$A)))

cat("Running tree method...\n")
fit3 <- estimate_att(X = d3$X, A = d3$A, Y = d3$Y, K = 5,
                regularization = log(400) / 400,
                use_rashomon = FALSE, verbose = TRUE)
cat(sprintf("Result: theta=%.4f, sigma=%.4f\n", fit3$theta, fit3$sigma))

cat("\n=== Test Complete ===\n")
cat("Both DGPs work successfully.\n")
