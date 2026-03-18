# Test if hang is in fit_tree() wrapper or in optimaltrees() itself

library(optimaltrees)

cat("\n=== TESTING DIRECT optimaltrees() CALL ===\n\n")

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

# Test 1: Call optimaltrees() directly WITHOUT discretization params
cat("Test 1: optimaltrees() WITHOUT discretization params...\n")
t1 <- Sys.time()
result1 <- optimaltrees::optimaltrees(
  X_tr, A_tr,
  loss_function = "log_loss",
  regularization = log(length(A_tr)) / length(A_tr),
  worker_limit = 4,
  single_tree = TRUE,
  verbose = FALSE
)
t2 <- Sys.time()
cat(sprintf("  Time: %.3fs, n_trees=%d\n\n", as.numeric(difftime(t2, t1, units = "secs")), result1$n_trees))

# Test 2: Call optimaltrees() directly WITH discretization params
cat("Test 2: optimaltrees() WITH discretization params...\n")
cat("  (If this hangs, bug is in optimaltrees() not fit_tree())\n")
flush.console()

t1 <- Sys.time()
result2 <- optimaltrees::optimaltrees(
  X_tr, A_tr,
  loss_function = "log_loss",
  discretize_method = "quantiles",
  discretize_bins = "adaptive",
  regularization = log(length(A_tr)) / length(A_tr),
  worker_limit = 4,
  single_tree = TRUE,
  verbose = FALSE
)
t2 <- Sys.time()
cat(sprintf("  Time: %.3fs, n_trees=%d\n\n", as.numeric(difftime(t2, t1, units = "secs")), result2$n_trees))

cat("Both tests completed!\n")
