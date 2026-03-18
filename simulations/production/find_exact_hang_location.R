# Find EXACT line where optimaltrees hangs after C++ completes

library(optimaltrees)

cat("\n=== INSTRUMENTING OPTIMALTREES TO FIND HANG ===\n\n")

# Generate same problematic data
source("dgps/dgps_smooth.R")
seed <- 100000 + 4*10000 + 1
d <- generate_dgp_continuous_att(n = 800, tau = 0.10, seed = seed)

# Get fold 1 training data (this is what hangs)
K <- 5
fold_ids <- sample(rep(1:K, length.out = nrow(d$X)))
train_idx <- which(fold_ids != 1)

X_train <- d$X[train_idx, ]
y_train <- d$A[train_idx]

cat(sprintf("Training data: n=%d, p=%d\n", nrow(X_train), ncol(X_train)))
cat(sprintf("Outcome: %d zeros, %d ones\n\n", sum(y_train == 0), sum(y_train == 1)))

# Try calling optimaltrees with extreme verbosity
cat("Calling optimaltrees...\n\n")
flush.console()

# Wrap with timeout
start_time <- Sys.time()

result <- system.time({
  fit <- tryCatch({
    optimaltrees(
      X = X_train,
      y = y_train,
      loss_function = "log_loss",
      regularization = log(length(y_train)) / length(y_train),
      worker_limit = 4,
      verbose = TRUE,
      store_training_data = FALSE,  # Try disabling to see if that helps
      compute_probabilities = FALSE  # Try disabling this too
    )
  }, error = function(e) {
    cat(sprintf("\nERROR: %s\n", conditionMessage(e)))
    NULL
  })
})

cat(sprintf("\n=== TIMING ===\n"))
print(result)

if (!is.null(fit)) {
  cat("\nSUCCESS!\n")
  cat(sprintf("n_trees: %d\n", fit$n_trees))
} else {
  cat("\nFAILED or HUNG\n")
}

cat("\n")
