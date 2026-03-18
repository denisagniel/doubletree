#!/usr/bin/env Rscript
# Production DGP2 Batch: 6,000 simulations
# Uses lessons from background execution debugging

cat("======================================================================\n")
cat("DGP2 (Continuous Outcome) - Production Simulations\n")
cat("======================================================================\n\n")
flush.console()

# Load packages
suppressMessages({
  library(dplyr)
  library(optimaltrees)
})

# Load functions
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

cat("✓ Packages and functions loaded\n\n")
flush.console()

# Configuration
N_VALUES <- c(400, 800, 1600)
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 500
SEED_OFFSET <- 10000

DGPS <- list(dgp2 = generate_dgp_continuous_att)
METHODS <- c("tree", "rashomon", "forest", "linear")

# Create grid
sim_grid <- expand.grid(
  dgp_name = names(DGPS),
  method = METHODS,
  n = N_VALUES,
  rep = 1:N_REPS,
  stringsAsFactors = FALSE
)
sim_grid$sim_id <- 1:nrow(sim_grid)

cat(sprintf("Total simulations: %d\n", nrow(sim_grid)))
cat(sprintf("Expected runtime: ~1.5 hours\n\n"))
flush.console()

# Simulation function
run_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  if (sim_id %% 50 == 0) {
    cat(sprintf("Progress: %d/%d (%.1f%%)\n", sim_id, nrow(grid), 100*sim_id/nrow(grid)))
    flush.console()
  }

  row <- grid[grid$sim_id == sim_id, ]
  dgp_func <- dgps[[row$dgp_name]]
  seed <- seed_offset + sim_id
  d <- dgp_func(n = row$n, tau = tau, seed = seed)

  result <- tryCatch({
    if (row$method == "tree") {
      fit <- suppress_all(estimate_att(X = d$X, A = d$A, Y = d$Y, K = k_folds,
                                   regularization = log(row$n)/row$n,
                                   use_rashomon = FALSE, verbose = FALSE))
      list(theta = fit$theta, sigma = fit$sigma, ci_lower = fit$ci_95[1],
           ci_upper = fit$ci_95[2], converged = TRUE, epsilon_n = NA)

    } else if (row$method == "rashomon") {
      fit <- suppress_all(estimate_att(X = d$X, A = d$A, Y = d$Y, K = k_folds,
                                   regularization = log(row$n)/row$n,
                                   use_rashomon = TRUE,
                                   rashomon_bound_multiplier = 2*sqrt(log(row$n)/row$n),
                                   auto_tune_intersecting = TRUE, verbose = FALSE))
      list(theta = fit$theta, sigma = fit$sigma, ci_lower = fit$ci_95[1],
           ci_upper = fit$ci_95[2], converged = fit$converged, epsilon_n = fit$epsilon_n)

    } else if (row$method == "forest") {
      fit <- suppress_all(method_forest_dml(X = d$X, A = d$A, Y = d$Y, K = k_folds))
      list(theta = fit$theta, sigma = fit$sigma, ci_lower = fit$ci_lower,
           ci_upper = fit$ci_upper, converged = TRUE, epsilon_n = NA)

    } else if (row$method == "linear") {
      fit <- suppress_all(method_linear_dml(X = d$X, A = d$A, Y = d$Y, K = k_folds))
      list(theta = fit$theta, sigma = fit$sigma, ci_lower = fit$ci_lower,
           ci_upper = fit$ci_upper, converged = TRUE, epsilon_n = NA)
    }
  }, error = function(e) {
    list(theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
         converged = FALSE, epsilon_n = NA)
  })

  data.frame(
    sim_id = sim_id, dgp = row$dgp_name, method = row$method,
    n = row$n, rep = row$rep, true_att = tau,
    theta = result$theta, sigma = result$sigma,
    ci_lower = result$ci_lower, ci_upper = result$ci_upper,
    converged = result$converged, epsilon_n = result$epsilon_n,
    stringsAsFactors = FALSE
  )
}

# Run simulations
cat("Starting simulations...\n")
flush.console()
start_time <- Sys.time()

results_list <- lapply(sim_grid$sim_id, function(id) {
  run_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET)
})

results <- dplyr::bind_rows(results_list)
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "hours"))

cat(sprintf("\n✓ Complete in %.2f hours\n", elapsed))
cat(sprintf("Convergence: %.1f%%\n", 100*mean(results$converged, na.rm=TRUE)))
flush.console()

# Save
output_dir <- sprintf("results/dgp2_batch_%s", Sys.Date())
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(results, file.path(output_dir, "dgp2_results.rds"))

cat(sprintf("\nSaved: %s/dgp2_results.rds\n", output_dir))
flush.console()

# Summary
summary_stats <- results %>%
  filter(converged) %>%
  group_by(method, n) %>%
  summarize(
    n_valid = n(),
    bias = mean(theta - true_att, na.rm = TRUE),
    rmse = sqrt(mean((theta - true_att)^2, na.rm = TRUE)),
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_stats, file.path(output_dir, "dgp2_summary.csv"), row.names = FALSE)

cat("\nSummary Statistics:\n")
print(summary_stats, n = Inf)
flush.console()

cat("\n======================================================================\n")
cat("DGP2 Batch Complete\n")
cat("======================================================================\n")
