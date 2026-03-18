# Sequential version of run_primary.R (mclapply removed due to macOS fork issues)
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

# DGP functions
DGPS <- list(
  dgp1 = generate_dgp_binary_att,
  dgp2 = generate_dgp_continuous_att,
  dgp3 = generate_dgp_moderate_att
)

# Method functions
METHODS <- c("tree", "rashomon", "forest", "linear")

cat("Running in SEQUENTIAL mode (mclapply issue on macOS)\\n\\n")

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
sim_grid\$sim_id <- 1:nrow(sim_grid)

cat(sprintf("Total simulations: %d\\n", nrow(sim_grid)))
cat(sprintf("Estimated runtime: %.1f hours (sequential)\\n\\n", nrow(sim_grid) * 3 / 3600))

# Load run_single_sim function from original script
source("run_primary.R", local = TRUE, echo = FALSE)

# Run simulations sequentially
cat("Starting simulations...\\n")
start_time <- Sys.time()

results_list <- lapply(
  sim_grid\$sim_id,
  function(id) run_single_sim(id, sim_grid, DGPS, TAU, K_FOLDS, SEED_OFFSET)
)

# Rest same as original...
results <- dplyr::bind_rows(results_list)
end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "hours"))

cat(sprintf("\\nSimulations complete in %.2f hours\\n", elapsed))
safe_save(results, file.path(output_dir, "simulation_results.rds"))
cat("Results saved\\n")
