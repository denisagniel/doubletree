#!/usr/bin/env Rscript

#' Run a single DML-ATT simulation replication
#'
#' Designed for SLURM array jobs where each task runs one replication.
#' This script is optimized for O2 cluster execution with minimal memory footprint.

library(optparse)
library(dplyr)

# Parse command-line arguments
option_list <- list(
  make_option(c("-d", "--dgp"), type = "character", default = NULL,
              help = "DGP name: dgp1, dgp2, or dgp3", metavar = "character"),
  make_option(c("-n", "--sample-size"), type = "integer", default = 400,
              help = "Sample size [default %default]", metavar = "number"),
  make_option(c("-m", "--method"), type = "character", default = "tree",
              help = "Method: tree, rashomon, forest, or linear", metavar = "character"),
  make_option(c("-r", "--replication"), type = "integer", default = 1,
              help = "Replication number [default %default]", metavar = "number"),
  make_option(c("-o", "--output-dir"), type = "character",
              default = "results/o2_primary",
              help = "Output directory for replication results", metavar = "path"),
  make_option(c("--tau"), type = "double", default = 0.10,
              help = "True treatment effect [default %default]", metavar = "number"),
  make_option(c("--k-folds"), type = "integer", default = 5,
              help = "Number of CV folds [default %default]", metavar = "number"),
  make_option(c("--seed-offset"), type = "integer", default = 10000,
              help = "Seed offset [default %default]", metavar = "number"),
  make_option(c("--worker-limit"), type = "integer", default = 1,
              help = "Worker threads (set to 1 for SLURM) [default %default]", metavar = "number"),
  make_option(c("--package-path"), type = "character", default = NULL,
              help = "Path to package root (for devtools::load_all)", metavar = "path")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$dgp)) {
  stop("Must specify --dgp (dgp1, dgp2, or dgp3)")
}

if (!opt$dgp %in% c("dgp1", "dgp2", "dgp3")) {
  stop("Invalid DGP: ", opt$dgp, " (must be dgp1, dgp2, or dgp3)")
}

if (!opt$method %in% c("tree", "rashomon", "forest", "linear")) {
  stop("Invalid method: ", opt$method, " (must be tree, rashomon, forest, or linear)")
}

# Load packages
suppressMessages({
  library(optimaltrees)
})

# Load doubletree from source
# If package-path provided, use devtools::load_all() for development
if (!is.null(opt$`package-path`)) {
  suppressMessages({
    library(devtools)
    devtools::load_all(opt$`package-path`, quiet = TRUE)
  })
} else {
  # Load from relative paths (assumes running from production/ directory)
  doubletree_files <- c(
    "../../R/utils.R",
    "../../R/score_att.R",
    "../../R/inference.R",
    "../../R/nuisance_trees.R",
    "../../R/estimate_att.R"
  )

  for (f in doubletree_files) {
    if (file.exists(f)) {
      source(f)
    } else {
      stop("Cannot find doubletree source file: ", f)
    }
  }
}

# Source DGPs and baseline methods
source("dgps/dgps_smooth.R")
source("methods/method_forest.R")
source("methods/method_linear.R")

# Set seed
seed <- opt$`seed-offset` + opt$replication

cat(sprintf("Replication %d: DGP=%s, n=%d, method=%s, seed=%d\n",
            opt$replication, opt$dgp, opt$`sample-size`, opt$method, seed))
flush.console()

# Generate data based on DGP
if (opt$dgp == "dgp1") {
  d <- generate_dgp_binary_att(n = opt$`sample-size`, tau = opt$tau, seed = seed)
} else if (opt$dgp == "dgp2") {
  d <- generate_dgp_continuous_att(n = opt$`sample-size`, tau = opt$tau, seed = seed)
} else if (opt$dgp == "dgp3") {
  d <- generate_dgp_moderate_att(n = opt$`sample-size`, tau = opt$tau, seed = seed)
}

# Fit model based on method
result <- tryCatch({

  if (opt$method == "tree") {
    # Tree-DML (fold-specific regularization)
    fit <- estimate_att(
      X = d$X, A = d$A, Y = d$Y,
      K = opt$`k-folds`,
      regularization = log(opt$`sample-size`) / opt$`sample-size`,
      cv_regularization = FALSE,
      use_rashomon = FALSE,
      worker_limit = opt$`worker-limit`,
      verbose = FALSE
    )

    list(
      theta = fit$theta,
      sigma = fit$sigma,
      ci_lower = fit$ci[1],
      ci_upper = fit$ci[2],
      converged = TRUE,
      epsilon_n = NA,
      error = NA
    )

  } else if (opt$method == "rashomon") {
    # Rashomon-DML (structure intersection)
    epsilon_n <- 2 * sqrt(log(opt$`sample-size`) / opt$`sample-size`)

    fit <- estimate_att(
      X = d$X, A = d$A, Y = d$Y,
      K = opt$`k-folds`,
      regularization = log(opt$`sample-size`) / opt$`sample-size`,
      cv_regularization = FALSE,
      use_rashomon = TRUE,
      rashomon_bound_multiplier = epsilon_n,
      worker_limit = opt$`worker-limit`,
      verbose = FALSE
    )

    list(
      theta = fit$theta,
      sigma = fit$sigma,
      ci_lower = fit$ci[1],
      ci_upper = fit$ci[2],
      converged = TRUE,
      epsilon_n = epsilon_n,
      error = NA
    )

  } else if (opt$method == "forest") {
    # Forest-DML (ranger)
    fit <- att_forest(
      X = d$X, A = d$A, Y = d$Y,
      K = opt$`k-folds`,
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
      epsilon_n = NA,
      error = NA
    )

  } else if (opt$method == "linear") {
    # Linear-DML (logistic regression)
    fit <- att_linear(
      X = d$X, A = d$A, Y = d$Y,
      K = opt$`k-folds`,
      seed = seed,
      interactions = FALSE,
      verbose = FALSE
    )

    list(
      theta = fit$theta,
      sigma = fit$sigma,
      ci_lower = fit$ci[1],
      ci_upper = fit$ci[2],
      converged = fit$convergence == "converged",
      epsilon_n = NA,
      error = NA
    )

  } else {
    stop("Unknown method: ", opt$method)
  }

}, error = function(e) {
  list(
    theta = NA,
    sigma = NA,
    ci_lower = NA,
    ci_upper = NA,
    converged = FALSE,
    epsilon_n = NA,
    error = conditionMessage(e)
  )
})

# Compile results
result_data <- data.frame(
  dgp = opt$dgp,
  n = opt$`sample-size`,
  method = opt$method,
  replication = opt$replication,
  true_att = d$true_att,
  theta = result$theta,
  sigma = result$sigma,
  ci_lower = result$ci_lower,
  ci_upper = result$ci_upper,
  converged = result$converged,
  epsilon_n = result$epsilon_n,
  error = as.character(result$error),
  stringsAsFactors = FALSE
)

# Create output directory if needed
if (!dir.exists(opt$`output-dir`)) {
  dir.create(opt$`output-dir`, recursive = TRUE)
}

# Save result
output_file <- file.path(
  opt$`output-dir`,
  sprintf("%s_n%d_%s_rep%04d.rds", opt$dgp, opt$`sample-size`, opt$method, opt$replication)
)

saveRDS(result_data, output_file)

cat(sprintf("Saved to: %s\n", output_file))
cat(sprintf("TRUE ATT = %.4f, θ̂ = %.4f [%.4f, %.4f], Converged: %s\n",
            d$true_att, result$theta, result$ci_lower, result$ci_upper,
            ifelse(result$converged, "YES", "NO")))

if (!result$converged) {
  cat(sprintf("ERROR: %s\n", result$error))
}

# Return exit code
if (result$converged) {
  quit(status = 0)
} else {
  quit(status = 1)
}
