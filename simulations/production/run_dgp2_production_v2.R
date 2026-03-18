#!/usr/bin/env Rscript
# Batch: DGP1 (Binary) - All methods, all sample sizes
# Total: 6,000 simulations (4 methods × 3 n × 500 reps)

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
N_REPS <- 500
SEED_OFFSET <- 10000

# Single DGP
DGPS <- list(dgp2 = generate_dgp_continuous_att)
METHODS <- c("tree", "rashomon", "forest", "linear")

cat(strrep("=", 70), "\n")
cat("BATCH: DGP2 (Continuous Outcome)\n")
cat(strrep("=", 70), "\n\n")

# Create output directory
output_dir <- sprintf("results/dgp2_batch_%s", Sys.Date())
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
      suppressWarnings(suppressMessages({
        fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = k_folds,
                       regularization = log(row$n) / row$n,
                       use_rashomon = FALSE, verbose = FALSE)
      }))
      list(theta = fit$theta, sigma = fit$sigma,
           ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
           converged = TRUE, epsilon_n = NA)

    } else if (row$method == "rashomon") {
      suppressWarnings(suppressMessages({
        fit <- estimate_att(X = d$X, A = d$A, Y = d$Y, K = k_folds,
                       regularization = log(row$n) / row$n,
                       use_rashomon = TRUE,
                       rashomon_bound_multiplier = 2 * sqrt(log(row$n) / row$n),
                       auto_tune_intersecting = TRUE, verbose = FALSE)
      }))
      list(theta = fit$theta, sigma = fit$sigma,
           ci_lower = fit$ci_95[1], ci_upper = fit$ci_95[2],
           converged = fit$converged, epsilon_n = fit$epsilon_n)

    } else if (row$method == "forest") {
      suppressWarnings(suppressMessages({
        fit <- method_forest_dml(X = d$X, A = d$A, Y = d$Y, K = k_folds)
      }))
      list(theta = fit$theta, sigma = fit$sigma,
           ci_lower = fit$ci[1], ci_upper = fit$ci[2],
           converged = TRUE, epsilon_n = NA)

    } else if (row$method == "linear") {
      suppressWarnings(suppressMessages({
        fit <- method_linear_dml(X = d$X, A = d$A, Y = d$Y, K = k_folds)
      }))
      list(theta = fit$theta, sigma = fit$sigma,
           ci_lower = fit$ci[1], ci_upper = fit$ci[2],
           converged = TRUE, epsilon_n = NA)
    }
  }, error = function(e) {
    cat(sprintf("    Error in sim %d: %s\n", sim_id, conditionMessage(e)))
    list(theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
         converged = FALSE, epsilon_n = NA)
  })

  # Ensure result is a proper list
  if (is.null(result) || !is.list(result)) {
    result <- list(theta = NA, sigma = NA, ci_lower = NA, ci_upper = NA,
                   converged = FALSE, epsilon_n = NA)
  }

  # Create result row
  data.frame(
    sim_id = sim_id,
    dgp = row$dgp_name,
    method = row$method,
    n = row$n,
    rep = row$rep,
    true_att = tau,
    theta = if(is.null(result$theta)) NA else result$theta,
    sigma = if(is.null(result$sigma)) NA else result$sigma,
    ci_lower = if(is.null(result$ci_lower)) NA else result$ci_lower,
    ci_upper = if(is.null(result$ci_upper)) NA else result$ci_upper,
    converged = if(is.null(result$converged)) FALSE else result$converged,
    epsilon_n = if(is.null(result$epsilon_n)) NA else result$epsilon_n,
    stringsAsFactors = FALSE
  )
}

# Run simulations
cat("Starting simulations...\n")
start_time <- Sys.time()

results_list <- lapply(
  sim_grid$sim_id,
  function(id) run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET)
)

results <- dplyr::bind_rows(results_list)
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "hours"))

cat(sprintf("\n✓ Simulations complete in %.2f hours\n", elapsed))
cat(sprintf("Convergence rate: %.1f%%\n", 100 * mean(results$converged, na.rm = TRUE)))

# Save results
saveRDS(results, file.path(output_dir, "dgp1_results.rds"))
cat(sprintf("Results saved to: %s\n", file.path(output_dir, "dgp1_results.rds")))

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

write.csv(summary_stats, file.path(output_dir, "dgp1_summary.csv"), row.names = FALSE)
cat("\nSummary Statistics:\n")
print(summary_stats, n = Inf)

cat("\n", strrep("=", 70), "\n")
cat("DGP1 Batch Complete\n")
cat(strrep("=", 70), "\n")
