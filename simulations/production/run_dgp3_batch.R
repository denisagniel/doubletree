#!/usr/bin/env Rscript
# Batch: DGP1 (Moderate) - All methods, all sample sizes
# Total: 6,000 simulations (4 methods × 3 n × 500 reps)

suppressMessages({
  library(dplyr)
  library(optimaltrees)
  library(parallel)
})

source("../simulation_helpers.R")
invisible(sapply(c(
  "../../R/estimate_att.R",
  "../../R/att_repeated.R",
  "../../R/nuisance_trees.R",
  "../../R/score_att.R",
  "../../R/inference.R",
  "../../R/utils.R",
  "dgps/dgps_smooth.R",
  "methods/method_forest.R",
  "methods/method_linear.R"
), safe_source))

# Configuration
N_VALUES <- c(400, 800, 1600)
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 500
SEED_OFFSET <- 10000

# Single DGP
DGPS <- list(dgp3 = generate_dgp_moderate_att)
METHODS <- c("tree", "rashomon", "forest", "linear")

cat(strrep("=", 70), "\n")
cat("BATCH: DGP1 (Moderate Outcome)\n")
cat(strrep("=", 70), "\n\n")

# Create output directory
output_dir <- sprintf("results/dgp3_batch_%s", Sys.Date())
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Create simulation grid
sim_grid <- expand.grid(
  dgp_name = names(DGPS),
  method = METHODS,
  n = N_VALUES,
  rep = 1:N_REPS,
  stringsAsFactors = FALSE
)
sim_grid$sim_id <- 1:nrow(sim_grid)

cat(sprintf("Total simulations: %d\n", nrow(sim_grid)))
cat(sprintf("Estimated runtime: %.1f hours\n\n", nrow(sim_grid) * 0.08 / 3600))

# Simulation function
run_single_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  row <- grid[grid$sim_id == sim_id, ]

  # Progress every 50 sims
  if (sim_id %% 50 == 0) {
    cat(sprintf("  Progress: %d/%d (%.1f%%)\n", sim_id, nrow(grid), 100*sim_id/nrow(grid)))
  }

  dgp_func <- dgps[[row$dgp_name]]
  seed <- seed_offset + sim_id
  d <- dgp_func(n = row$n, tau = tau, seed = seed)

  result <- tryCatch({
    if (row$method == "tree") {
      fit <- suppressWarnings(suppressMessages({
        estimate_att(X = d$X, A = d$A, Y = d$Y, K = k_folds,
                regularization = log(row$n) / row$n,
                use_rashomon = FALSE, verbose = FALSE)
      }))
      list(theta = fit$theta, sigma = fit$sigma,
           ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
           converged = TRUE, epsilon_n = NA)

    } else if (row$method == "rashomon") {
      fit <- suppressWarnings(suppressMessages({
        estimate_att(X = d$X, A = d$A, Y = d$Y, K = k_folds,
                regularization = log(row$n) / row$n,
                use_rashomon = TRUE,
                rashomon_bound_multiplier = 2 * sqrt(log(row$n) / row$n),
                auto_tune_intersecting = TRUE, verbose = FALSE)
      }))
      list(theta = fit$theta, sigma = fit$sigma,
           ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
           converged = fit$converged, epsilon_n = fit$epsilon_n)

    } else if (row$method == "forest") {
      fit <- suppressWarnings(suppressMessages({
        att_forest(X = d$X, A = d$A, Y = d$Y, K = k_folds)
      }))
      list(theta = fit$theta, sigma = fit$sigma,
           ci_lower = fit$ci[1], ci_upper = fit$ci[2],
           converged = TRUE, epsilon_n = NA)

    } else if (row$method == "linear") {
      fit <- suppressWarnings(suppressMessages({
        att_linear(X = d$X, A = d$A, Y = d$Y, K = k_folds)
      }))
      list(theta = fit$theta, sigma = fit$sigma,
           ci_lower = fit$ci[1], ci_upper = fit$ci[2],
           converged = TRUE, epsilon_n = NA)
    }
  }, error = function(e) {
    list(theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
         converged = FALSE, epsilon_n = NA)
  })

  data.frame(
    sim_id = sim_id,
    dgp = row$dgp_name,
    method = row$method,
    n = row$n,
    rep = row$rep,
    true_att = tau,
    theta = result$theta,
    sigma = result$sigma,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    converged = result$converged,
    epsilon_n = result$epsilon_n,
    stringsAsFactors = FALSE
  )
}

# Run simulations (parallelized with 3 cores)
cat("Starting simulations...\n")
cat("Using 3 cores for parallel processing\n")
start_time <- Sys.time()

results_list <- mclapply(
  sim_grid$sim_id,
  function(id) run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET),
  mc.cores = 3
)

results <- dplyr::bind_rows(results_list)
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "hours"))

cat(sprintf("\n✓ Simulations complete in %.2f hours\n", elapsed))
cat(sprintf("Convergence rate: %.1f%%\n", 100 * mean(results$converged, na.rm = TRUE)))

# Save results
saveRDS(results, file.path(output_dir, "dgp3_results.rds"))
cat(sprintf("Results saved to: %s\n", file.path(output_dir, "dgp3_results.rds")))

# Summary stats
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

write.csv(summary_stats, file.path(output_dir, "dgp3_summary.csv"), row.names = FALSE)
cat("\nSummary Statistics:\n")
print(summary_stats, n = Inf)

cat("\n", strrep("=", 70), "\n")
cat("DGP1 Batch Complete\n")
cat(strrep("=", 70), "\n")
