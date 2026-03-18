# test_traversal.R
# Debug tree traversal to see why predictions are constant

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

source("simulations/dgps_realistic.R")

n <- 50  # Smaller dataset for easier debugging
d <- generate_dgp_simple(n, tau = 0.15, seed = 123)

message("Fitting tree...")
tree <- optimaltrees::fit_tree(
  d$X, d$A,
  loss_function = "log_loss",
  regularization = 0.01,
  verbose = FALSE
)

tree_json <- tree$model$tree_json

message("\n=== Manual Tree Traversal Test ===")
message("Root node:")
message("  $feature: ", tree_json$feature, " (", typeof(tree_json$feature), ")")
message("  $feature numeric: ", as.numeric(tree_json$feature))
message("  Feature index (0-based): ", tree_json$feature)
message("  Feature index (1-based): ", as.integer(as.numeric(tree_json$feature) + 1))
message("  Feature name: ", tree_json$name)

message("\nFirst observation in d$X:")
print(d$X[1, ])

# Manually traverse for first observation
obs_1 <- d$X[1, ]
message("\nManual traversal for obs 1:")

feature_idx <- as.integer(as.numeric(tree_json$feature) + 1)
message("  Root split on feature ", feature_idx, " (", names(d$X)[feature_idx], ")")
message("  Obs 1 value: ", obs_1[[feature_idx]])

if (obs_1[[feature_idx]] == 1) {
  message("  → Going to TRUE branch")
  next_node <- tree_json$true
} else {
  message("  → Going to FALSE branch")
  next_node <- tree_json$false
}

if (!is.null(next_node$prediction)) {
  message("  Reached leaf with prediction: ", next_node$prediction)
  if (!is.null(next_node$probabilities)) {
    message("  Leaf probabilities: ", paste(next_node$probabilities, collapse=", "))
  }
} else if (!is.null(next_node$feature)) {
  feature_idx_2 <- as.integer(as.numeric(next_node$feature) + 1)
  message("  Split on feature ", feature_idx_2, " (", names(d$X)[feature_idx_2], ")")
  message("  Obs 1 value: ", obs_1[[feature_idx_2]])
}

message("\n=== Testing get_probabilities_from_tree ===")
get_probabilities_from_tree <- get("get_probabilities_from_tree", envir = asNamespace("optimaltrees"))
probs <- get_probabilities_from_tree(tree_json, d$X)

message("First 5 predictions:")
print(head(probs, 5))
message("\nPrediction SD: ", sd(probs[,2]))

# Check if all rows are getting the same prediction
if (length(unique(probs[,2])) == 1) {
  message("\n✗ ALL PREDICTIONS IDENTICAL: ", unique(probs[,2]))
  message("This suggests traverse_batch is not working correctly")
} else {
  message("\n✓ Predictions vary (", length(unique(probs[,2])), " unique values)")
}
