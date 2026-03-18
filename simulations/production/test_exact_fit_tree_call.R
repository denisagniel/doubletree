# Test the EXACT fit_tree call that's hanging

library(optimaltrees)

cat("\n=== TESTING EXACT fit_tree CALL ===\n\n")

source("dgps/dgps_smooth.R")
source("../../R/utils.R")

# Generate exact same data
seed <- 100000 + 4*10000 + 1
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = seed)

# Create exact same folds
K <- 5
set.seed(42)  # dmltree default
fold_indices <- sample(rep(1:K, length.out = nrow(d$X)))

# Get Fold 1 training data (EXACT indices that hang)
fold_id <- 1
train_idx <- which(fold_indices != fold_id)

X_tr <- d$X[train_idx, , drop = FALSE]
A_tr <- d$A[train_idx]

cat(sprintf("Training data: n=%d\n", nrow(X_tr)))
cat(sprintf("Outcome: %d zeros, %d ones\n\n", sum(A_tr == 0), sum(A_tr == 1)))

# Test WITHOUT discretization parameters (like our earlier test)
cat("Test 1: fit_tree WITHOUT discretize_* parameters...\n")
t1 <- Sys.time()
fit1 <- optimaltrees::fit_tree(
  X_tr, A_tr,
  loss_function = "log_loss",
  regularization = log(length(A_tr)) / length(A_tr),
  worker_limit = 4,
  verbose = FALSE
)
t2 <- Sys.time()
cat(sprintf("  Time: %.3fs, n_trees=%d\n\n", as.numeric(difftime(t2, t1, units = "secs")), fit1$n_trees))

# Test WITH discretization parameters (like dml_att uses)
cat("Test 2: fit_tree WITH discretize_* parameters...\n")
cat("  (This should hang if discretization is the issue)\n")
flush.console()

t1 <- Sys.time()
fit2 <- optimaltrees::fit_tree(
  X_tr, A_tr,
  loss_function = "log_loss",
  discretize_method = "quantiles",
  discretize_bins = "adaptive",
  regularization = log(length(A_tr)) / length(A_tr),
  worker_limit = 4,
  verbose = FALSE
)
t2 <- Sys.time()
cat(sprintf("  Time: %.3fs, n_trees=%d\n\n", as.numeric(difftime(t2, t1, units = "secs")), fit2$n_trees))

cat("Both tests completed!\n")
