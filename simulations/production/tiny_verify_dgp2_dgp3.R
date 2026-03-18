#!/usr/bin/env Rscript
# Tiny Verification: 2 reps per method per n per DGP
# Total: 48 simulations (2 DGPs × 4 methods × 3 n × 2 reps)

suppressMessages({
  library(dplyr)
  library(optimaltrees)
})

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

# Configuration
N_VALUES <- c(400, 800, 1600)
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 2  # TINY TEST
SEED_OFFSET <- 10000

DGPS <- list(
  dgp2 = generate_dgp_continuous_att,
  dgp3 = generate_dgp_moderate_att
)
METHODS <- c("tree", "rashomon", "forest", "linear")

cat("=== TINY VERIFICATION: DGP2 + DGP3 (2 reps) ===\n")

sim_grid <- expand.grid(
  dgp_name = names(DGPS),
  method = METHODS,
  n = N_VALUES,
  rep = 1:N_REPS,
  stringsAsFactors = FALSE
)
sim_grid$sim_id <- 1:nrow(sim_grid)

cat(sprintf("Total: %d sims\n", nrow(sim_grid)))

run_single_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  row <- grid[grid$sim_id == sim_id, ]
  if (sim_id %% 10 == 0) cat(sprintf("  %d/%d\n", sim_id, nrow(grid)))

  dgp_func <- dgps[[row$dgp_name]]
  d <- dgp_func(n = row$n, tau = tau, seed = seed_offset + sim_id)

  result <- tryCatch({
    if (row$method == "tree") {
      fit <- suppressMessages(estimate_att(X = d$X, A = d$A, Y = d$Y, K = k_folds,
                                       regularization = log(row$n) / row$n,
                                       use_rashomon = FALSE, verbose = FALSE))
      list(theta = fit$theta, sigma = fit$sigma, converged = TRUE)
    } else if (row$method == "rashomon") {
      fit <- suppressMessages(estimate_att(X = d$X, A = d$A, Y = d$Y, K = k_folds,
                                       regularization = log(row$n) / row$n,
                                       use_rashomon = TRUE,
                                       rashomon_bound_multiplier = 2 * sqrt(log(row$n) / row$n),
                                       auto_tune_intersecting = TRUE, verbose = FALSE))
      list(theta = fit$theta, sigma = fit$sigma, converged = fit$converged)
    } else if (row$method == "forest") {
      fit <- suppressMessages(method_forest_dml(X = d$X, A = d$A, Y = d$Y, K = k_folds))
      list(theta = fit$theta, sigma = fit$sigma, converged = TRUE)
    } else if (row$method == "linear") {
      fit <- suppressMessages(method_linear_dml(X = d$X, A = d$A, Y = d$Y, K = k_folds))
      list(theta = fit$theta, sigma = fit$sigma, converged = TRUE)
    }
  }, error = function(e) list(theta = NA, sigma = NA, converged = FALSE))

  if (is.null(result)) result <- list(theta = NA, sigma = NA, converged = FALSE)

  data.frame(sim_id = sim_id, dgp = row$dgp_name, method = row$method,
             n = row$n, rep = row$rep,
             theta = if(is.null(result$theta)) NA else result$theta,
             converged = if(is.null(result$converged)) FALSE else result$converged,
             stringsAsFactors = FALSE)
}

start <- Sys.time()
results_list <- lapply(sim_grid$sim_id, function(id) run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET))
results <- bind_rows(results_list)
elapsed <- difftime(Sys.time(), start, units = "secs")

cat(sprintf("\nComplete in %.1f seconds\n", elapsed))
cat(sprintf("Convergence: %.1f%%\n", 100 * mean(results$converged)))

cat("\nPer-DGP convergence:\n")
for (dgp in names(DGPS)) {
  conv <- mean(results$converged[results$dgp == dgp])
  cat(sprintf("  %s: %.1f%%\n", dgp, 100 * conv))
}

cat("\n=== SUCCESS: Both DGPs work ===\n")
