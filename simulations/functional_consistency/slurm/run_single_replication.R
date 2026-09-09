#!/usr/bin/env Rscript
# Single replication for functional consistency simulation
# Command-line interface for O2 cluster execution

library(optparse)

# Parse command-line arguments
option_list <- list(
  make_option(c("--n"), type = "integer", dest = "n",
              help = "Sample size"),
  make_option(c("--dgp"), type = "character", dest = "dgp",
              help = "DGP type: simple, complex, or sparse"),
  make_option(c("--method"), type = "character", dest = "method",
              help = "Method: standard_msplit, averaged_tree, or pattern_aggregation"),
  make_option(c("--K"), type = "integer", dest = "K",
              help = "Number of folds"),
  make_option(c("--M"), type = "integer", dest = "M",
              help = "Number of splits"),
  make_option(c("--seed"), type = "integer", dest = "seed",
              help = "Random seed"),
  make_option(c("--output"), type = "character", dest = "output",
              help = "Output RDS file path")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
required <- c("n", "dgp", "method", "K", "M", "seed", "output")
missing <- required[!required %in% names(opt) | sapply(opt[required], is.null)]
if (length(missing) > 0) {
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

# Load packages (try devtools for local development, library for O2)
# TWO PACKAGE-LOADING PATHS, ONE GATE (added 2026-09-09 for SLURM deployment):
#
#   DOUBLETREE_USE_INSTALLED unset/"0"  DEV path (default, unchanged behaviour):
#       devtools::load_all() on local source trees. Correct on the dev box.
#
#   DOUBLETREE_USE_INSTALLED="1"        CLUSTER path: library() against packages
#       that were R CMD INSTALLed on the cluster. Module R does not carry a working
#       pkgload/devtools dev-load.
suppressPackageStartupMessages({
  USE_INSTALLED_PKGS <- Sys.getenv("DOUBLETREE_USE_INSTALLED", "0") == "1"
  if (USE_INSTALLED_PKGS) {
    library(optimaltrees)
    library(doubletree)
    cat("Loaded via library (installed packages)\n")
  } else {
    # Local development: use devtools
    library(devtools)
    local_opt <- file.path(Sys.getenv("HOME"), "RAND/rprojects/global-scholars/optimaltrees")
    local_dbl <- file.path(Sys.getenv("HOME"), "RAND/rprojects/global-scholars/doubletree")
    load_all(local_opt, quiet = TRUE)
    load_all(local_dbl, quiet = TRUE)
    cat("Loaded via devtools (local dev)\n")
  }
})

# Source simulation function
source("run_fc_simulation.R")

# Run simulation
cat("Running simulation with:\n")
cat("  n =", opt$n, "\n")
cat("  dgp =", opt$dgp, "\n")
cat("  method =", opt$method, "\n")
cat("  K =", opt$K, "\n")
cat("  M =", opt$M, "\n")
cat("  seed =", opt$seed, "\n")

result <- run_fc_simulation(
  n = opt$n,
  dgp = opt$dgp,
  method = opt$method,
  K = opt$K,
  M = opt$M,
  seed = opt$seed
)

# Save result
saveRDS(result, file = opt$output)
cat("Results saved to:", opt$output, "\n")
