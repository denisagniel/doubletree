#' Apply Discretization Metadata to New Data
#'
#' Helper function to discretize new data using metadata from a treefarmr model.
#' This is needed because treefarmr's predict() doesn't auto-discretize.
#'
#' @param X_new Data.frame of continuous features
#' @param metadata Discretization metadata from treefarmr model
#' @return Data.frame with binary features matching training discretization
#' @noRd
apply_discretization_metadata <- function(X_new, metadata) {
  if (is.null(metadata) || is.null(metadata$features)) {
    # No discretization needed (already binary)
    return(X_new)
  }

  binary_cols_list <- list()

  for (col_name in names(X_new)) {
    x <- X_new[[col_name]]
    feat_meta <- metadata$features[[col_name]]

    if (is.null(feat_meta)) {
      # Feature not in metadata - pass through
      binary_cols_list[[col_name]] <- data.frame(x)
      names(binary_cols_list[[col_name]]) <- col_name
      next
    }

    # Handle different feature types
    if (feat_meta$type == "binary") {
      binary_cols_list[[col_name]] <- data.frame(x)
      names(binary_cols_list[[col_name]]) <- col_name

    } else if (feat_meta$type == "binary_converted") {
      # Map to {0, 1} using original values
      x_binary <- as.numeric(x == max(feat_meta$original_values))
      binary_cols_list[[col_name]] <- data.frame(x_binary)
      names(binary_cols_list[[col_name]]) <- col_name

    } else if (feat_meta$type == "constant") {
      # All zeros (constant feature)
      x_binary <- rep(0, length(x))
      binary_cols_list[[col_name]] <- data.frame(x_binary)
      names(binary_cols_list[[col_name]]) <- feat_meta$new_names

    } else if (feat_meta$type == "continuous") {
      # Apply thresholds to create binary indicators
      thresholds <- feat_meta$thresholds
      new_names <- feat_meta$new_names

      # Create threshold indicators: I(x <= t_i) for each threshold
      threshold_cols <- lapply(seq_along(thresholds), function(i) {
        as.integer(x <= thresholds[i])
      })

      binary_cols <- as.data.frame(threshold_cols)
      names(binary_cols) <- new_names
      binary_cols_list[[col_name]] <- binary_cols

    } else {
      stop("Unknown feature type in metadata: ", feat_meta$type)
    }
  }

  # Combine all binary columns
  X_binary <- do.call(cbind, binary_cols_list)

  return(X_binary)
}
