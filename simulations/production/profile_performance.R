# Performance Profiling: Identify bottlenecks in DML-ATT estimation
# Uses Rprof to profile a single DML-ATT run

library(optimaltrees)

cat("\n=== PERFORMANCE PROFILING ===\n\n")

# Source doubletree
source("../../R/utils.R")
source("../../R/score_att.R")
source("../../R/inference.R")
source("../../R/nuisance_trees.R")
source("../../R/estimate_att.R")

# Source DGPs
source("dgps/dgps_smooth.R")

cat("Generating test data (n=400, binary DGP)...\n")
set.seed(12345)
d <- generate_dgp_binary_att(n = 400, tau = 0.10, seed = 12345)

cat("Starting profiler...\n\n")

# Profile tree-DML
cat("=== PROFILING: Tree-DML ===\n")
Rprof("profile_tree.out", interval = 0.01)
fit_tree <- estimate_att(
  X = d$X, A = d$A, Y = d$Y,
  K = 5,
  outcome_type = "binary",
  regularization = log(400) / 400,
  cv_regularization = FALSE,
  use_rashomon = FALSE,
  worker_limit = 4,
  verbose = FALSE
)
Rprof(NULL)
cat("Theta:", fit_tree$theta, "\n\n")

# Profile rashomon-DML
cat("=== PROFILING: Rashomon-DML ===\n")
Rprof("profile_rashomon.out", interval = 0.01)
fit_rashomon <- estimate_att(
  X = d$X, A = d$A, Y = d$Y,
  K = 5,
  outcome_type = "binary",
  regularization = log(400) / 400,
  cv_regularization = FALSE,
  use_rashomon = TRUE,
  worker_limit = 4,
  verbose = FALSE
)
Rprof(NULL)
cat("Theta:", fit_rashomon$theta, "\n\n")

# Analyze profiles
cat("=== PROFILE SUMMARIES ===\n\n")

cat("--- Tree-DML Profile ---\n")
profile_tree <- summaryRprof("profile_tree.out")
cat("\nTop functions by total time:\n")
print(head(profile_tree$by.total, 10))
cat("\nTop functions by self time:\n")
print(head(profile_tree$by.self, 10))

cat("\n--- Rashomon-DML Profile ---\n")
profile_rashomon <- summaryRprof("profile_rashomon.out")
cat("\nTop functions by total time:\n")
print(head(profile_rashomon$by.total, 10))
cat("\nTop functions by self time:\n")
print(head(profile_rashomon$by.self, 10))

# Check for performance issues
cat("\n=== PERFORMANCE DIAGNOSTICS ===\n\n")

# Check if C++ functions dominate
tree_cpp_time <- sum(profile_tree$by.self$self.time[grepl("Rcpp|C\\+\\+", rownames(profile_tree$by.self))])
tree_r_time <- sum(profile_tree$by.self$self.time[!grepl("Rcpp|C\\+\\+", rownames(profile_tree$by.self))])

cat("Tree-DML time breakdown:\n")
cat(sprintf("  C++ time: %.2f sec\n", tree_cpp_time))
cat(sprintf("  R time: %.2f sec\n", tree_r_time))
cat(sprintf("  Ratio: %.1f%% C++ / %.1f%% R\n",
            100*tree_cpp_time/(tree_cpp_time+tree_r_time),
            100*tree_r_time/(tree_cpp_time+tree_r_time)))

cat("\n=== RECOMMENDATIONS ===\n\n")

total_tree_time <- sum(profile_tree$by.total$total.time)
total_rashomon_time <- sum(profile_rashomon$by.total$total.time)

cat(sprintf("Tree-DML total profiled time: %.2f sec\n", total_tree_time))
cat(sprintf("Rashomon-DML total profiled time: %.2f sec\n", total_rashomon_time))

if (tree_r_time > tree_cpp_time) {
  cat("\n⚠️  R code dominates runtime - potential for optimization\n")
} else {
  cat("\n✓ C++ code dominates runtime - expected for tree optimization\n")
}

# Clean up
unlink("profile_tree.out")
unlink("profile_rashomon.out")

cat("\n=== DONE ===\n")
