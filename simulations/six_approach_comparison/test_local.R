#!/usr/bin/env Rscript

# Local testing script - loads packages from source instead of requiring install

cat("Loading packages from source...\n")
suppressPackageStartupMessages({
  devtools::load_all("../../../optimaltrees", quiet = TRUE)
  devtools::load_all("../..", quiet = TRUE)  # doubletree
})

cat("Loading simulation code...\n")
source("code/dgps.R")
source("code/estimators.R")
source("code/metrics.R")

# Test parameters
approaches <- 1:6
approach_names <- c(
  "full_sample", "crossfit_separate", "doubletree",
  "doubletree_singlefit", "msplit", "msplit_singlefit"
)

dgps <- c(1, 4)  # Test simple and continuous
dgp_names <- c("simple", "moderate", "complex", "continuous")

# Map approach number to function
approach_map <- list(
  `1` = estimate_att_fullsample,
  `2` = estimate_att_crossfit,
  `3` = estimate_att_doubletree,
  `4` = estimate_att_doubletree_singlefit,
  `5` = estimate_att_msplit,
  `6` = estimate_att_msplit_singlefit
)

# Map DGP number to function
dgp_map <- list(
  `1` = generate_dgp_simple,
  `2` = generate_dgp_moderate,
  `3` = generate_dgp_complex,
  `4` = generate_dgp_continuous
)

cat("\n===========================================\n")
cat("Testing Six-Approach Comparison (Local)\n")
cat("===========================================\n\n")

results_summary <- list()
idx <- 1

# Test each approach with DGP 1 (simple)
for (i in approaches) {
  cat(sprintf("\n--- Testing Approach %d (%s) with DGP 1 (simple) ---\n",
              i, approach_names[i]))

  set.seed(12345)
  data <- dgp_map[[1]](n = 500)

  start_time <- Sys.time()
  result <- tryCatch({
    approach_map[[i]](X = data$X, A = data$A, Y = data$Y)
  }, error = function(e) {
    list(theta = NA, se = NA, error = as.character(e))
  })
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  if (!is.null(result$error)) {
    cat(sprintf("  ERROR: %s\n", result$error))
  } else {
    cat(sprintf("  theta_hat = %.4f\n", result$theta))
    cat(sprintf("  se = %.4f\n", result$se))
    cat(sprintf("  true ATT = %.4f\n", data$true_att))
    cat(sprintf("  bias = %.4f\n", result$theta - data$true_att))
    cat(sprintf("  time = %.2f sec\n", elapsed))
  }

  results_summary[[idx]] <- data.frame(
    approach = i,
    approach_name = approach_names[i],
    dgp = 1,
    dgp_name = "simple",
    theta_hat = result$theta,
    se = result$se,
    true_att = data$true_att,
    elapsed = elapsed,
    error = if(is.null(result$error)) NA_character_ else result$error
  )
  idx <- idx + 1
}

# Test each approach with DGP 4 (continuous)
for (i in approaches) {
  cat(sprintf("\n--- Testing Approach %d (%s) with DGP 4 (continuous) ---\n",
              i, approach_names[i]))

  set.seed(67890)
  data <- dgp_map[[4]](n = 500)

  start_time <- Sys.time()
  result <- tryCatch({
    approach_map[[i]](X = data$X, A = data$A, Y = data$Y)
  }, error = function(e) {
    list(theta = NA, se = NA, error = as.character(e))
  })
  elapsed <- as.numeric(Sys.time() - start_time, units = "secs")

  if (!is.null(result$error)) {
    cat(sprintf("  ERROR: %s\n", result$error))
  } else {
    cat(sprintf("  theta_hat = %.4f\n", result$theta))
    cat(sprintf("  se = %.4f\n", result$se))
    cat(sprintf("  true ATT = %.4f\n", data$true_att))
    cat(sprintf("  bias = %.4f\n", result$theta - data$true_att))
    cat(sprintf("  time = %.2f sec\n", elapsed))
  }

  results_summary[[idx]] <- data.frame(
    approach = i,
    approach_name = approach_names[i],
    dgp = 4,
    dgp_name = "continuous",
    theta_hat = result$theta,
    se = result$se,
    true_att = data$true_att,
    elapsed = elapsed,
    error = if(is.null(result$error)) NA_character_ else result$error
  )
  idx <- idx + 1
}

# Combine results
results_df <- do.call(rbind, results_summary)

cat("\n\n===========================================\n")
cat("Summary of Local Tests\n")
cat("===========================================\n\n")

cat(sprintf("Total tests: %d\n", nrow(results_df)))
cat(sprintf("Successful: %d\n", sum(!is.na(results_df$theta_hat))))
cat(sprintf("Errors: %d\n", sum(is.na(results_df$theta_hat))))

if (sum(is.na(results_df$theta_hat)) > 0) {
  cat("\nErrors by approach:\n")
  error_summary <- results_df[is.na(results_df$theta_hat),
                               c("approach_name", "dgp_name", "error")]
  print(error_summary)
}

cat("\n\nResults by approach:\n")
print(results_df[, c("approach_name", "dgp_name", "theta_hat", "se", "elapsed")])

cat("\n✓ Local testing complete!\n")
