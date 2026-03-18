#' Primary Simulations for Manuscript Table 1
#'
#' Comprehensive comparison of 4 DML-ATT methods across 3 smooth DGPs
#' and 3 sample sizes.
#'
#' **Grid:**
#' - DGPs: 1-3 (smooth, validated)
#' - Methods: tree-DML, rashomon-DML, forest-DML, linear-DML
#' - Sample sizes: n ∈ {400, 800, 1600}
#' - Replications: 500 per config
#' - Total: 3 × 4 × 3 × 500 = 18,000 runs
#'
#' **Output:**
#' - results/primary_YYYY-MM-DD/simulation_results.rds (full replication data)
#' - results/primary_YYYY-MM-DD/summary_stats.csv (aggregated metrics)
#'
#' **Runtime:** ~15 hours single-threaded; ~4 hours with 4-core parallelization
#'
#' **Three-way fidelity:**
#' - Paper Section 4, Table 1 will cite this script
#' - Methods match paper descriptions (lines 157-162)
#' - DGPs match paper specifications
#'
#' **Tuning parameters:**
#' - Tree-DML: λ = log(n)/n (minimax-optimal rate, theory-justified)
#' - Rashomon-DML: ε_n = 2√(log(n)/n) (theory-justified Rashomon bound)
#' - Forest-DML: 500 trees, mtry = sqrt(p) (standard ranger defaults)
#' - Linear-DML: no interactions (main effects only)
#' - All methods: K = 5 folds (standard for DML)

# Load packages and functions
suppressMessages({
  library(parallel)
  library(dplyr)
  library(optimaltrees)  # Use installed package instead of devtools::load_all()
})

# Load simulation helpers (prevents log bloat)
source("../simulation_helpers.R")

# Source dmltree functions SILENTLY
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
N_REPS <- 500  # For production; use 50-100 for testing
SEED_OFFSET <- 10000

# DGP functions
DGPS <- list(
  dgp1 = generate_dgp_binary_att,
  dgp2 = generate_dgp_continuous_att,
  dgp3 = generate_dgp_moderate_att
)

# Method functions
METHODS <- c("tree", "rashomon", "forest", "linear")

# Sequential mode (mclapply has stability issues on macOS)
N_CORES <- 1
cat("Running in sequential mode\n\n")

# Create output directory
output_dir <- sprintf("results/primary_%s", Sys.Date())
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
cat(sprintf("Estimated runtime: %.1f hours (sequential mode)\n\n",
            nrow(sim_grid) * 3 / 3600))

# Simulation function for single replication
run_single_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  row <- grid[grid$sim_id == sim_id, ]

  # Progress tracking: console only, infrequent (prevents log bloat)
  progress_msg(sim_id, nrow(grid), every = 100, label = "Simulations")

  # Memory monitoring every 10 simulations
  if (sim_id %% 10 == 0) {
    monitor_memory(threshold_mb = 12000, verbose = FALSE)
  }

  # Generate data
  dgp_func <- dgps[[row$dgp_name]]
  seed <- seed_offset + sim_id
  d <- dgp_func(n = row$n, tau = tau, seed = seed)

  # Fit model based on method - ALL wrapped in suppress_all()
  result <- tryCatch({

    if (row$method == "tree") {
      # Tree-DML (fold-specific regularization)
      fit <- suppress_all({
        estimate_att(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          regularization = log(row$n) / row$n,
          cv_regularization = FALSE,
          use_rashomon = FALSE,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = TRUE,
        epsilon_n = NA  # Not applicable for tree-DML
      )

    } else if (row$method == "rashomon") {
      # Rashomon-DML (structure intersection)
      # Use theory-justified epsilon_n = 2 * sqrt(log(n)/n)
      epsilon_n <- 2 * sqrt(log(row$n) / row$n)

      fit <- suppress_all({
        estimate_att(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          regularization = log(row$n) / row$n,
          cv_regularization = FALSE,
          use_rashomon = TRUE,
          rashomon_bound_multiplier = epsilon_n,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = TRUE,
        epsilon_n = epsilon_n  # Store for diagnostics
      )

    } else if (row$method == "forest") {
      # Forest-DML (ranger)
      fit <- suppress_all({
        att_forest(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          seed = seed,
          num.trees = 500,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = fit$convergence == "converged",
        epsilon_n = NA  # Not applicable for forest-DML
      )

    } else if (row$method == "linear") {
      # Linear-DML (logistic regression)
      fit <- suppress_all({
        att_linear(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          seed = seed,
          interactions = FALSE,
          verbose = FALSE
        )
      })

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = fit$convergence == "converged",
        epsilon_n = NA  # Not applicable for linear-DML
      )

    } else {
      stop("Unknown method: ", row$method)
    }

  }, error = function(e) {
    list(
      theta = NA,
      sigma = NA,
      ci_lower = NA,
      ci_upper = NA,
      converged = FALSE,
      epsilon_n = NA
    )
  })

  # Return results with metadata
  data.frame(
    sim_id = sim_id,
    dgp = row$dgp_name,
    method = row$method,
    n = row$n,
    rep = row$rep,
    true_att = d$true_att,
    theta = result$theta,
    sigma = result$sigma,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    converged = result$converged,
    epsilon_n = result$epsilon_n,  # Rashomon bound (NA for non-Rashomon methods)
    stringsAsFactors = FALSE
  )
}

# Run simulations (sequential - mclapply has fork issues on macOS)
cat("Starting simulations...\n")
cat("Running in sequential mode due to mclapply stability issues\n\n")
start_time <- Sys.time()

results_list <- lapply(
  sim_grid$sim_id,
  function(id) run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET)
)

# Combine results
results <- do.call(rbind, results_list)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "hours"))

cat(sprintf("\nSimulations complete in %.2f hours\n", elapsed))
cat(sprintf("Convergence rate: %.1f%%\n", 100 * mean(results$converged, na.rm = TRUE)))

# Save full results (atomic write)
safe_save(results, file.path(output_dir, "simulation_results.rds"))

# Compute summary statistics
summary_stats <- results %>%
  filter(converged) %>%
  group_by(dgp, method, n) %>%
  summarize(
    n_valid = n(),
    bias = mean(theta - true_att, na.rm = TRUE),
    rmse = sqrt(mean((theta - true_att)^2, na.rm = TRUE)),
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att, na.rm = TRUE),
    mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE),
    median_ci_width = median(ci_upper - ci_lower, na.rm = TRUE),
    .groups = "drop"
  )

# Save summary statistics
write.csv(summary_stats,
          file.path(output_dir, "summary_stats.csv"),
          row.names = FALSE)
cat(sprintf("Summary stats saved to: %s\n", file.path(output_dir, "summary_stats.csv")))

# Print summary
cat("\n")
cat("=" %R% 70, "\n")
cat("Summary Statistics (by DGP, Method, Sample Size)\n")
cat("=" %R% 70, "\n\n")
print(summary_stats, n = Inf)

# Flag any concerning results
cat("\n")
cat("=" %R% 70, "\n")
cat("Quality Checks\n")
cat("=" %R% 70, "\n\n")

# Check coverage
bad_coverage <- summary_stats %>%
  filter(coverage < 0.90 | coverage > 0.98)

if (nrow(bad_coverage) > 0) {
  cat("⚠️  Configurations with poor coverage (<90% or >98%):\n")
  print(bad_coverage[, c("dgp", "method", "n", "coverage")])
} else {
  cat("✓ All methods achieve 90-98% coverage\n")
}

# Check bias
bad_bias <- summary_stats %>%
  filter(abs(bias) > 0.03)  # >30% of tau = 0.10

if (nrow(bad_bias) > 0) {
  cat("\n⚠️  Configurations with large bias (>0.03):\n")
  print(bad_bias[, c("dgp", "method", "n", "bias")])
} else {
  cat("✓ All methods have bias < 0.03 (30% of truth)\n")
}

cat("\n")
cat("=" %R% 70, "\n")
cat("Simulation complete. Results ready for manuscript Table 1.\n")
cat("=" %R% 70, "\n")

# Final check: warn if large files were created
check_large_files(output_dir, min_mb = 10)
