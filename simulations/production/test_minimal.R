# Minimal test to diagnose hang
cat("Loading packages...\\n")
suppressMessages({
  library(parallel)
  library(dplyr)
  library(optimaltrees)
})

cat("Loading helpers...\\n")
source("../simulation_helpers.R")

cat("Loading dmltree functions...\\n")
invisible(sapply(c(
  "../../R/estimate_att.R",
  "../../R/nuisance_trees.R",
  "../../R/score_att.R",
  "../../R/inference.R",
  "../../R/utils.R",
  "dgps/dgps_smooth.R"
), safe_source))

cat("Creating minimal grid...\\n")
sim_grid <- expand.grid(
  dgp_name = "dgp1",
  method = "tree",
  n = 400,
  rep = 1:2,
  stringsAsFactors = FALSE
)
sim_grid$sim_id <- 1:nrow(sim_grid)

cat("Grid created:", nrow(sim_grid), "simulations\\n")

# Simplified run function
run_sim <- function(sim_id) {
  cat("  Starting sim", sim_id, "...\\n")

  row <- sim_grid[sim_grid$sim_id == sim_id, ]
  dgp_func <- generate_dgp_binary_att
  seed <- 10000 + sim_id
  d <- dgp_func(n = row$n, tau = 0.1, seed = seed)

  cat("    Data generated, fitting model...\\n")

  fit <- estimate_att(
    X = d$X, A = d$A, Y = d$Y,
    K = 2,
    regularization = 0.01,
    use_rashomon = FALSE,
    verbose = FALSE
  )

  cat("    Fit complete: theta =", fit$theta, "\\n")

  list(
    sim_id = sim_id,
    theta = fit$theta,
    converged = TRUE
  )
}

cat("Running simulations sequentially...\\n")
results <- list()
for (i in 1:nrow(sim_grid)) {
  results[[i]] <- run_sim(i)
}

cat("\\n✓ All simulations complete\\n")
cat("Results:", sapply(results, function(x) x$theta), "\\n")
