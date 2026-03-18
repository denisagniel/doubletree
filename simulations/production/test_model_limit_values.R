# Test what model_limit values are being computed in different scenarios

library(optimaltrees)

cat("\n=== TESTING MODEL_LIMIT COMPUTATION ===\n\n")

source("dgps/dgps_smooth.R")
source("../../R/utils.R")

# Generate exact same data
seed <- 100000 + 4*10000 + 1
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = seed)

# Create exact same folds
K <- 5
set.seed(42)
fold_indices <- sample(rep(1:K, length.out = nrow(d$X)))

# Get Fold 1 training data
fold_id <- 1
train_idx <- which(fold_indices != fold_id)

X_tr <- d$X[train_idx, , drop = FALSE]
A_tr <- d$A[train_idx]

cat(sprintf("Training data: n=%d, p=%d\n\n", nrow(X_tr), ncol(X_tr)))

# Compute what fit_tree() would compute for model_limit
# (from fit_tree.R lines 85-113)
p_original <- ncol(X_tr)
n <- nrow(X_tr)

# When discretize_bins = "adaptive":
n_bins_adaptive <- max(2, ceiling(log(n) / 3))
cat(sprintf("n_bins (adaptive): %d\n", n_bins_adaptive))

# Estimate binary feature count
p_estimated <- p_original * (n_bins_adaptive - 1)
cat(sprintf("p_original: %d\n", p_original))
cat(sprintf("p_estimated: %d\n", p_estimated))

# Compute model_limit
if (p_estimated > 100) {
  model_limit <- 1000000
} else if (p_estimated > 50) {
  model_limit <- 100000
} else {
  model_limit <- 0  # 0 = unlimited
}
cat(sprintf("model_limit computed by fit_tree(): %d (0 = unlimited)\n\n", model_limit))

# Test 1: fit_tree() with verbose=TRUE to see what's happening
cat("Test 1: fit_tree() with verbose=TRUE and discretization...\n")
cat("  (Checking if it hangs or if we can see the issue)\n")
flush.console()

t1 <- Sys.time()
result1 <- optimaltrees::fit_tree(
  X_tr, A_tr,
  loss_function = "log_loss",
  discretize_method = "quantiles",
  discretize_bins = "adaptive",
  regularization = log(length(A_tr)) / length(A_tr),
  worker_limit = 4,
  verbose = TRUE
)
t2 <- Sys.time()
cat(sprintf("  Time: %.3fs, n_trees=%d\n\n", as.numeric(difftime(t2, t1, units = "secs")), result1$n_trees))

cat("Test completed!\n")
