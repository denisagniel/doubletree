# test_validation.R
# Debug why tree validation is failing

if (!requireNamespace("devtools", quietly = TRUE)) stop("devtools required")
devtools::load_all()

source("simulations/dgps_realistic.R")

n <- 400
d <- generate_dgp_simple(n, tau = 0.15, seed = 123)

message("Fitting tree...")
tree <- optimaltrees::fit_tree(
  d$X, d$A,
  loss_function = "log_loss",
  regularization = 0.01,
  verbose = FALSE
)

message("\n=== Testing tree validation ===")
tree_json <- tree$model$tree_json

# Get validation function
validate_tree_structure <- get("validate_tree_structure", envir = asNamespace("optimaltrees"))

message("Tree JSON structure:")
message("  Is list: ", is.list(tree_json))
message("  Length: ", length(tree_json))
message("  Names: ", paste(names(tree_json), collapse=", "))
message("  Has $feature: ", !is.null(tree_json$feature))
message("  Has $prediction: ", !is.null(tree_json$prediction))
message("  Has $true: ", !is.null(tree_json$true))
message("  Has $false: ", !is.null(tree_json$false))

message("\nValidation result: ", validate_tree_structure(tree_json))

# If validation fails, dig deeper
if (!validate_tree_structure(tree_json)) {
  message("\n✗ VALIDATION FAILED! Debugging...")

  # Check root node details
  message("\nRoot node $feature value: ", tree_json$feature)
  message("Root node $feature type: ", typeof(tree_json$feature))
  message("Root node $feature class: ", class(tree_json$feature))

  # Check children
  if (!is.null(tree_json$true)) {
    message("\n$true branch exists")
    message("  $true is list: ", is.list(tree_json$true))
    message("  $true has $feature: ", !is.null(tree_json$true$feature))
    message("  $true has $prediction: ", !is.null(tree_json$true$prediction))
  }

  if (!is.null(tree_json$false)) {
    message("\n$false branch exists")
    message("  $false is list: ", is.list(tree_json$false))
    message("  $false has $feature: ", !is.null(tree_json$false$feature))
    message("  $false has $prediction: ", !is.null(tree_json$false$prediction))
  }

  # Try validating children separately
  message("\n\nValidating children separately:")
  if (!is.null(tree_json$true)) {
    message("  $true valid: ", validate_tree_structure(tree_json$true))
  }
  if (!is.null(tree_json$false)) {
    message("  $false valid: ", validate_tree_structure(tree_json$false))
  }
} else {
  message("✓ Validation passed!")

  # But predictions still fail, so test prediction
  message("\n=== Testing prediction ===")
  get_probabilities_from_tree <- get("get_probabilities_from_tree", envir = asNamespace("optimaltrees"))
  probs <- get_probabilities_from_tree(tree_json, d$X)
  message("Predictions SD: ", sd(probs[,2]))

  if (sd(probs[,2]) < 0.01) {
    message("✗ Predictions still constant even though validation passed!")
  }
}
