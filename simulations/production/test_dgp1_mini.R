# Quick test of DGP1 batch (2 reps only)
suppressMessages({
  library(dplyr)
  library(optimaltrees)
})

source("../simulation_helpers.R")
invisible(sapply(c(
  "../../R/estimate_att.R",
  "../../R/nuisance_trees.R",
  "../../R/score_att.R",
  "../../R/inference.R",
  "../../R/utils.R",
  "dgps/dgps_smooth.R"
), safe_source))

cat("Quick test: DGP1 with 2 reps\n\n")

dgp_func <- generate_dgp_binary_att

for (i in 1:2) {
  cat("Rep", i, "...")
  d <- dgp_func(n = 400, tau = 0.1, seed = 10000 + i)
  fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = 2,
                 regularization = 0.01, use_rashomon = FALSE, verbose = FALSE)
  cat(" theta =", round(fit$theta, 4), "\n")
}

cat("\n✓ Test successful - batch script should work\n")
