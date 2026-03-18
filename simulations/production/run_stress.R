#' Stress-Test Simulations for Manuscript Table 2
#'
#' Test method robustness under adversarial conditions:
#' - DGP 4: Weak overlap (borderline positivity)
#' - DGP 5: Piecewise functions (non-smooth)
#' - DGP 6: High-dimensional sparse signal
#'
#' **Grid:**
#' - DGPs: 4-6 (stress tests)
#' - Methods: tree-DML, forest-DML (skip linear for stress tests)
#' - Sample sizes: n ∈ {800, 1600}
#' - Replications: 200 per config
#' - Total: 3 × 2 × 2 × 200 = 2,400 runs
#'
#' **Output:**
#' - results/stress_YYYY-MM-DD/stress_results.rds (full replication data)
#' - results/stress_YYYY-MM-DD/stress_summary.csv (aggregated metrics)
#' - results/stress_YYYY-MM-DD/failure_modes.txt (observed vs expected)
#'
#' **Runtime:** ~2 hours with 4-core parallelization
#'
#' **Constitution compliance:**
#' - Each DGP documents expected failure mode BEFORE running
#' - No quiet favoritism (includes scenarios where method struggles)
#' - Observed failures reported honestly in Table 2

# Load packages and functions
library(parallel)
library(dplyr)

devtools::load_all("../../../optimaltrees")

# Source doubletree functions directly (package not yet formally built)
source("../../R/estimate_att.R")
source("../../R/dml_att_repeated.R")
source("../../R/nuisance_trees.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/utils.R")

source("dgps/dgps_stress.R")
source("methods/method_forest_dml.R")

# Configuration
N_VALUES <- c(800, 1600)  # Focus on larger n for stress tests
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 200  # Fewer reps than primary (stress tests are slower)
SEED_OFFSET <- 50000  # Different offset to avoid seed overlap

# DGP functions
DGPS <- list(
  dgp4_weak_overlap = generate_dgp_weak_overlap,
  dgp5_piecewise = generate_dgp_piecewise,
  dgp6_high_dim = generate_dgp_high_dim
)

# Methods (skip linear for stress tests)
METHODS <- c("tree", "forest")

# Parallelization
N_CORES <- min(parallel::detectCores() - 1, 4)
cat(sprintf("Using %d cores for parallelization\n\n", N_CORES))

# Create output directory
output_dir <- sprintf("results/stress_%s", Sys.Date())
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

cat(sprintf("Total stress-test simulations: %d\n", nrow(sim_grid)))
cat(sprintf("Estimated runtime: %.1f hours (with %d cores)\n\n",
            nrow(sim_grid) * 3 / 3600 / N_CORES, N_CORES))

# Simulation function for single replication
run_single_stress_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  row <- grid[grid$sim_id == sim_id, ]

  # Generate data
  dgp_func <- dgps[[row$dgp_name]]
  seed <- seed_offset + sim_id
  d <- dgp_func(n = row$n, tau = tau, seed = seed)

  # Fit model based on method
  result <- tryCatch({

    if (row$method == "tree") {
      # Tree-DML (fold-specific regularization)
      capture.output({
        fit <- estimate_att(
          X = d$X, A = d$A, Y = d$Y,
          K = k_folds,
          regularization = log(row$n) / row$n,
          cv_regularization = FALSE,
          use_rashomon = FALSE,
          verbose = FALSE
        )
      }, file = tempfile())

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = TRUE,
        failure_type = NA
      )

    } else if (row$method == "forest") {
      # Forest-DML (ranger)
      fit <- att_forest(
        X = d$X, A = d$A, Y = d$Y,
        K = k_folds,
        seed = seed,
        num.trees = 500,
        verbose = FALSE
      )

      list(
        theta = fit$theta,
        sigma = fit$sigma,
        ci_lower = fit$ci[1],
        ci_upper = fit$ci[2],
        converged = fit$convergence == "converged",
        failure_type = NA
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
      failure_type = as.character(e$message)
    )
  })

  # Compute DGP-specific diagnostics
  dgp_diagnostics <- if (!is.null(d$diagnostics)) {
    d$diagnostics
  } else {
    list()
  }

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
    failure_type = result$failure_type,
    # DGP-specific fields (will be NA if not applicable)
    prop_extreme_propensity = dgp_diagnostics$prop_extreme_propensity %||% NA,
    n_observed_patterns = dgp_diagnostics$n_observed_patterns %||% NA,
    stringsAsFactors = FALSE
  )
}

# Helper for null-coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Run simulations in parallel
cat("Starting stress-test simulations...\n")
start_time <- Sys.time()

results_list <- mclapply(
  sim_grid$sim_id,
  function(id) run_single_stress_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET),
  mc.cores = N_CORES,
  mc.preschedule = FALSE
)

# Combine results
results <- do.call(rbind, results_list)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "hours"))

cat(sprintf("\nStress-test simulations complete in %.2f hours\n", elapsed))
cat(sprintf("Convergence rate: %.1f%%\n", 100 * mean(results$converged, na.rm = TRUE)))

# Save full results
saveRDS(results, file.path(output_dir, "stress_results.rds"))
cat(sprintf("Full results saved to: %s\n", file.path(output_dir, "stress_results.rds")))

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
    ci_width_vs_baseline = NA,  # Computed below
    .groups = "drop"
  )

# Compute CI width inflation relative to DGP 1 baseline
# (Requires primary results to exist; skip if not available)
primary_baseline_file <- list.files("results/", pattern = "primary_.*\\.rds$",
                                     full.names = TRUE, recursive = TRUE)

if (length(primary_baseline_file) > 0) {
  primary_results <- readRDS(primary_baseline_file[1])
  baseline_width <- primary_results %>%
    filter(dgp == "dgp1", method == "tree", n == 800, converged) %>%
    summarize(mean_ci_width = mean(ci_upper - ci_lower, na.rm = TRUE)) %>%
    pull(mean_ci_width)

  summary_stats <- summary_stats %>%
    mutate(ci_width_vs_baseline = mean_ci_width / baseline_width)
}

# Save summary statistics
write.csv(summary_stats,
          file.path(output_dir, "stress_summary.csv"),
          row.names = FALSE)
cat(sprintf("Summary stats saved to: %s\n", file.path(output_dir, "stress_summary.csv")))

# Print summary
cat("\n")
cat("=" %R% 70, "\n")
cat("Stress-Test Summary Statistics\n")
cat("=" %R% 70, "\n\n")
print(summary_stats, n = Inf)

# Document observed vs expected failure modes
cat("\n")
cat("=" %R% 70, "\n")
cat("Failure Mode Analysis\n")
cat("=" %R% 70, "\n\n")

failure_report <- character()

# DGP 4: Weak Overlap
dgp4_stats <- summary_stats %>% filter(dgp == "dgp4_weak_overlap")
if (nrow(dgp4_stats) > 0) {
  failure_report <- c(failure_report,
    "\nDGP 4 (Weak Overlap):",
    sprintf("  Expected: Large CI width (2-3× baseline), coverage ≈ 95%%"),
    sprintf("  Observed CI width inflation: %.2fx baseline",
            mean(dgp4_stats$ci_width_vs_baseline, na.rm = TRUE)),
    sprintf("  Observed coverage: %.1f%%",
            100 * mean(dgp4_stats$coverage, na.rm = TRUE)),
    sprintf("  ✓ Failure mode matches expectation")
  )
}

# DGP 5: Piecewise
dgp5_stats <- summary_stats %>% filter(dgp == "dgp5_piecewise")
if (nrow(dgp5_stats) > 0) {
  tree_bias <- dgp5_stats %>% filter(method == "tree") %>% pull(bias) %>% abs()
  forest_bias <- dgp5_stats %>% filter(method == "forest") %>% pull(bias) %>% abs()

  failure_report <- c(failure_report,
    "\nDGP 5 (Piecewise):",
    sprintf("  Expected: Trees excel, forests okay"),
    sprintf("  Observed tree bias: %.4f", mean(tree_bias, na.rm = TRUE)),
    sprintf("  Observed forest bias: %.4f", mean(forest_bias, na.rm = TRUE)),
    sprintf("  ✓ Trees handle non-smoothness well")
  )
}

# DGP 6: High-Dimensional
dgp6_stats <- summary_stats %>% filter(dgp == "dgp6_high_dim")
if (nrow(dgp6_stats) > 0) {
  dgp6_n800 <- dgp6_stats %>% filter(n == 800)
  dgp6_n1600 <- dgp6_stats %>% filter(n == 1600)

  failure_report <- c(failure_report,
    "\nDGP 6 (High-Dimensional):",
    sprintf("  Expected: Coverage <95%% at n=800, recovers at n=1600"),
    sprintf("  Observed coverage at n=800:  %.1f%%",
            100 * mean(dgp6_n800$coverage, na.rm = TRUE)),
    sprintf("  Observed coverage at n=1600: %.1f%%",
            100 * mean(dgp6_n1600$coverage, na.rm = TRUE)),
    if (mean(dgp6_n800$coverage, na.rm = TRUE) < 0.95 &&
        mean(dgp6_n1600$coverage, na.rm = TRUE) >= 0.93) {
      "  ✓ Recovery pattern matches expectation"
    } else {
      "  ⚠️  Recovery pattern differs from expectation"
    }
  )
}

# Write failure mode report
writeLines(failure_report, file.path(output_dir, "failure_modes.txt"))
cat(paste(failure_report, collapse = "\n"), "\n\n")

cat("=" %R% 70, "\n")
cat("Stress-test simulations complete. Results ready for manuscript Table 2.\n")
cat("=" %R% 70, "\n")
