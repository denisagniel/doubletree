#!/usr/bin/env Rscript

#' Run a batch of DML-ATT simulation replications
#'
#' Designed for efficient SLURM execution by running multiple replications
#' in a single task (targeting 5 minutes per task).
#'
#' Usage:
#'   Rscript run_batch_replications.R --dgp dgp1 --sample-size 400 --method tree \
#'           --batch-start 1 --batch-size 100 --output-dir /path/to/results

library(optparse)

# Parse command-line arguments
option_list <- list(
  make_option(c("-d", "--dgp"), type = "character", default = NULL,
              help = "DGP name: dgp1 through dgp9", metavar = "character"),
  make_option(c("-n", "--sample-size"), type = "integer", default = 400,
              help = "Sample size [default %default]", metavar = "number"),
  make_option(c("-m", "--method"), type = "character", default = "tree",
              help = "Method: tree, rashomon, forest, or linear", metavar = "character"),
  make_option(c("-s", "--batch-start"), type = "integer", default = 1,
              help = "First replication number in this batch [default %default]", metavar = "number"),
  make_option(c("-b", "--batch-size"), type = "integer", default = 100,
              help = "Number of replications in this batch [default %default]", metavar = "number"),
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
              help = "Worker threads (set to 1 for SLURM) [default %default]", metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$dgp)) {
  stop("Must specify --dgp (dgp1 through dgp9)")
}

if (!opt$dgp %in% c("dgp1", "dgp2", "dgp3", "dgp4", "dgp5", "dgp6", "dgp7", "dgp8", "dgp9")) {
  stop("Invalid DGP: ", opt$dgp, " (must be dgp1 through dgp9)")
}

if (!opt$method %in% c("tree", "rashomon", "forest", "linear")) {
  stop("Invalid method: ", opt$method, " (must be tree, rashomon, forest, or linear)")
}

# Load packages
suppressMessages({
  library(optimaltrees)
  library(doubletree)
})

# Source DGPs and baseline methods
source("dgps/dgps_smooth.R")
source("dgps/dgps_continuous.R")
source("dgps/dgps_phase2.R")
source("methods/method_forest.R")
source("methods/method_linear.R")

# Create output directory if needed
if (!dir.exists(opt$`output-dir`)) {
  dir.create(opt$`output-dir`, recursive = TRUE)
}

# Calculate batch end
batch_end <- opt$`batch-start` + opt$`batch-size` - 1

cat(sprintf("Starting batch: DGP=%s, n=%d, method=%s, reps %d-%d (%d total)\n",
            opt$dgp, opt$`sample-size`, opt$method,
            opt$`batch-start`, batch_end, opt$`batch-size`))
cat(sprintf("Output: %s\n", opt$`output-dir`))
cat(sprintf("Started at: %s\n\n", Sys.time()))

# Track success/failure
n_success <- 0
n_failure <- 0
failed_reps <- integer(0)

# Loop through replications in this batch
for (rep_num in opt$`batch-start`:batch_end) {

  seed <- opt$`seed-offset` + rep_num

  # Generate data based on DGP
  if (opt$dgp == "dgp1") {
    d <- generate_dgp_binary_att(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "binary"
  } else if (opt$dgp == "dgp2") {
    d <- generate_dgp_continuous_att(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "binary"
  } else if (opt$dgp == "dgp3") {
    d <- generate_dgp_moderate_att(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "binary"
  } else if (opt$dgp == "dgp4") {
    d <- generate_dgp_continuous_binary(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "binary"
  } else if (opt$dgp == "dgp5") {
    d <- generate_dgp_continuous_continuous(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "continuous"
  } else if (opt$dgp == "dgp6") {
    d <- generate_dgp_mixed(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "binary"
  } else if (opt$dgp == "dgp7") {
    d <- generate_dgp7(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "binary"
  } else if (opt$dgp == "dgp8") {
    d <- generate_dgp8(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "continuous"
  } else if (opt$dgp == "dgp9") {
    d <- generate_dgp9(n = opt$`sample-size`, tau = opt$tau, seed = seed)
    outcome_type <- "binary"
  } else {
    stop(sprintf("Unknown DGP: %s", opt$dgp))
  }

  # Fit model based on method
  result <- tryCatch({

    if (opt$method == "tree") {
      # Tree-DML (fold-specific regularization)
      fit <- estimate_att(
        X = d$X, A = d$A, Y = d$Y,
        K = opt$`k-folds`,
        outcome_type = outcome_type,
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
        outcome_type = outcome_type,
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
    replication = rep_num,
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

  # Save result
  output_file <- file.path(
    opt$`output-dir`,
    sprintf("%s_n%d_%s_rep%04d.rds", opt$dgp, opt$`sample-size`, opt$method, rep_num)
  )

  saveRDS(result_data, output_file)

  # Track success/failure
  if (result$converged) {
    n_success <- n_success + 1
  } else {
    n_failure <- n_failure + 1
    failed_reps <- c(failed_reps, rep_num)
  }

  # Progress message (every 10 reps or at end)
  if (rep_num %% 10 == 0 || rep_num == batch_end) {
    progress <- rep_num - opt$`batch-start` + 1
    pct <- round(100 * progress / opt$`batch-size`, 1)
    cat(sprintf("[%s] Completed %d/%d (%.1f%%) - Success: %d, Failed: %d\n",
                format(Sys.time(), "%H:%M:%S"), progress, opt$`batch-size`, pct,
                n_success, n_failure))
  }
}

# Final summary
cat(sprintf("\n===========================================\n"))
cat(sprintf("Batch complete at: %s\n", Sys.time()))
cat(sprintf("Configuration: %s, n=%d, %s\n", opt$dgp, opt$`sample-size`, opt$method))
cat(sprintf("Replications: %d-%d (%d total)\n", opt$`batch-start`, batch_end, opt$`batch-size`))
cat(sprintf("Success: %d (%.1f%%)\n", n_success, 100 * n_success / opt$`batch-size`))
cat(sprintf("Failed: %d (%.1f%%)\n", n_failure, 100 * n_failure / opt$`batch-size`))
if (n_failure > 0) {
  cat(sprintf("Failed reps: %s\n", paste(failed_reps, collapse = ", ")))
}
cat(sprintf("===========================================\n"))

# Exit with success
quit(save = "no", status = 0)
