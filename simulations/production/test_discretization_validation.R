# Test if discretization + validation is the hang

library(optimaltrees)

cat("\n=== TESTING DISCRETIZATION + VALIDATION ===\n\n")

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

# Test discretization
cat("Step 1: Discretize features...\n")
t1 <- Sys.time()
discret_result <- optimaltrees:::discretize_features(
  X = X_tr,
  method = "quantiles",
  n_bins = "adaptive",
  thresholds = NULL
)
t2 <- Sys.time()
cat(sprintf("  Time: %.3fs\n", as.numeric(difftime(t2, t1, units = "secs"))))
cat(sprintf("  Result: %d rows, %d cols\n\n", nrow(discret_result$X_binary), ncol(discret_result$X_binary)))

X_binary <- discret_result$X_binary

# Test validation (from treefarms.R line 515)
cat("Step 2: Validate binary features...\n")
t1 <- Sys.time()
m <- as.matrix(X_binary)
bad <- !m %in% c(0L, 1L) & !is.na(m)
if (any(bad)) {
  cat("  ERROR: Non-binary values found!\n")
} else {
  cat("  All values are binary\n")
}
t2 <- Sys.time()
cat(sprintf("  Time: %.3fs\n\n", as.numeric(difftime(t2, t1, units = "secs"))))

# Test creating data frame for C++
cat("Step 3: Create data frame with class column...\n")
t1 <- Sys.time()
data_df <- cbind(X_binary, class = A_tr)
t2 <- Sys.time()
cat(sprintf("  Time: %.3fs\n", as.numeric(difftime(t2, t1, units = "secs"))))
cat(sprintf("  Data: %d rows, %d cols\n\n", nrow(data_df), ncol(data_df)))

# Test calling C++ directly
cat("Step 4: Call treefarms_fit_with_config_cpp...\n")
cat("  (This is where it might hang)\n")
flush.console()

# Create minimal config
config <- list(
  loss_function = "log_loss",
  regularization = log(nrow(X_tr)) / nrow(X_tr),
  verbose = FALSE,
  worker_limit = 4
)

config_json <- jsonlite::toJSON(config, auto_unbox = TRUE, pretty = FALSE)

# Convert to CSV string
csv_string <- paste(capture.output(write.table(
  data_df, sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE
)), collapse = "\n")

cat(sprintf("  CSV string length: %d bytes\n", nchar(csv_string)))
cat("  Calling C++...\n")
flush.console()

t1 <- Sys.time()
result <- treefarms_fit_with_config_cpp(csv_string, config_json)
t2 <- Sys.time()

cat(sprintf("  Time: %.3fs\n", as.numeric(difftime(t2, t1, units = "secs"))))
cat("  SUCCESS!\n\n")

cat("All steps completed without hanging!\n")
