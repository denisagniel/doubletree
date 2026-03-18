# Small test: 3 DGPs × 1 method × 1 size × 5 reps = 15 runs
library(parallel)
library(dplyr)

devtools::load_all("../../../optimaltrees")

# Source doubletree functions directly
source("../../R/estimate_att.R")
source("../../R/nuisance_trees.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/utils.R")

source("dgps/dgps_beta_regimes.R")

# Small grid
N_VALUES <- c(800)
TAU <- 0.10
K_FOLDS <- 5
N_REPS <- 5
SEED_OFFSET <- 20000

DGPS <- list(
  beta_high = generate_dgp_beta_high,
  beta_boundary = generate_dgp_beta_boundary,
  beta_low = generate_dgp_beta_low
)

METHODS <- c("tree")

output_dir <- sprintf("results/beta_test_%s", Sys.Date())
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sim_grid <- expand.grid(
  dgp_name = names(DGPS),
  method = METHODS,
  n = N_VALUES,
  rep = 1:N_REPS,
  stringsAsFactors = FALSE
)
sim_grid$sim_id <- 1:nrow(sim_grid)

cat(sprintf("Test grid: %d simulations\n\n", nrow(sim_grid)))

run_single_sim <- function(sim_id, grid, dgps, tau, k_folds, seed_offset) {
  row <- grid[grid$sim_id == sim_id, ]
  
  dgp_func <- dgps[[row$dgp_name]]
  seed <- seed_offset + sim_id
  d <- dgp_func(n = row$n, tau = tau, seed = seed)
  
  result <- tryCatch({
    fit <- estimate_att(
      X = d$X, A = d$A, Y = d$Y,
      K = k_folds,
      regularization = log(row$n) / row$n,
      cv_regularization = FALSE,
      use_rashomon = FALSE,
      verbose = FALSE  # Suppress verbose output
    )
    
    list(
      theta = fit$theta,
      sigma = fit$sigma,
      ci_lower = fit$ci[1],
      ci_upper = fit$ci[2],
      converged = TRUE
    )
  }, error = function(e) {
    list(
      theta = NA,
      sigma = NA,
      ci_lower = NA,
      ci_upper = NA,
      converged = FALSE,
      error = as.character(e)
    )
  })
  
  data.frame(
    sim_id = sim_id,
    dgp = row$dgp_name,
    method = row$method,
    n = row$n,
    rep = row$rep,
    true_att = d$true_att,
    beta = d$diagnostics$beta,
    rate_regime = d$diagnostics$rate_regime,
    theta = result$theta,
    sigma = result$sigma,
    ci_lower = result$ci_lower,
    ci_upper = result$ci_upper,
    converged = result$converged,
    stringsAsFactors = FALSE
  )
}

cat("Running test simulations...\n")
start_time <- Sys.time()

results_list <- lapply(
  sim_grid$sim_id,
  function(id) {
    cat(sprintf("  Sim %d/%d\r", id, nrow(sim_grid)))
    run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET)
  }
)

results <- do.call(rbind, results_list)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat(sprintf("\n\nTest complete in %.1f seconds\n", elapsed))
cat(sprintf("Convergence rate: %.1f%%\n", 100 * mean(results$converged, na.rm = TRUE)))

saveRDS(results, file.path(output_dir, "test_results.rds"))

# Summary
summary_stats <- results %>%
  filter(converged) %>%
  group_by(dgp, beta, n) %>%
  summarize(
    n_valid = n(),
    coverage = mean(ci_lower <= true_att & ci_upper >= true_att),
    rmse = sqrt(mean((theta - true_att)^2)),
    .groups = "drop"
  )

cat("\nTest Results:\n")
print(summary_stats)

cat(sprintf("\n✓ Test successful. Results in: %s\n", output_dir))
