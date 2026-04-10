#' Structure Utilities for M-Split
#'
#' @description
#' Helper functions for M-split algorithm: modal structure selection,
#' structure frequency analysis.
#'
#' @name structure_utils
NULL

#' Select Modal Structure from List
#'
#' @description
#' Given a list of TreeStructure objects, identify the modal (most frequent)
#' structure using hash-based comparison.
#'
#' @param structures List of TreeStructure objects (from extract_tree_structure)
#' @return List with:
#'   \item{structure}{TreeStructure object (modal structure)}
#'   \item{frequency}{Numeric: proportion of structures matching modal (in [0, 1])}
#'   \item{counts}{Named integer vector: frequency of each unique structure}
#'   \item{hash}{Character: hash of modal structure}
#'
#' @details
#' Uses optimaltrees::structure_hash() for O(1) comparison via hash table
#' instead of O(M^2) pairwise structure comparisons.
#'
#' If there are ties (multiple structures with same frequency), returns the
#' first one encountered.
#'
#' @examples
#' \dontrun{
#' # Fit M models on different splits
#' structures <- vector("list", M)
#' for (m in 1:M) {
#'   model <- fit_tree(X[[m]], y[[m]])
#'   structures[[m]] <- optimaltrees::extract_tree_structure(model)
#' }
#'
#' # Select modal structure
#' modal <- select_structure_modal(structures)
#' cat(sprintf("Modal frequency: %.1f%%\n", modal$frequency * 100))
#' }
#'
#' @export
select_structure_modal <- function(structures) {
  if (!is.list(structures) || length(structures) == 0) {
    stop("structures must be a non-empty list", call. = FALSE)
  }

  # Check all elements are TreeStructure
  is_tree_struct <- vapply(structures, function(s) {
    S7::S7_inherits(s, optimaltrees::TreeStructure)
  }, logical(1))

  if (!all(is_tree_struct)) {
    stop("All elements of structures must be TreeStructure objects", call. = FALSE)
  }

  # Compute hashes
  hashes <- vapply(structures, optimaltrees::structure_hash, character(1))

  # Count frequencies
  counts <- table(hashes)

  # Find modal hash
  modal_hash <- names(counts)[which.max(counts)]
  modal_idx <- which(hashes == modal_hash)[1]  # First occurrence

  list(
    structure = structures[[modal_idx]],
    frequency = as.numeric(max(counts)) / length(structures),
    counts = counts,
    hash = modal_hash
  )
}

#' Analyze Structure Diversity
#'
#' @description
#' Compute summary statistics for structure diversity across M splits.
#'
#' @param structures List of TreeStructure objects
#' @return List with:
#'   \item{n_unique}{Integer: number of unique structures}
#'   \item{shannon_entropy}{Numeric: Shannon entropy of structure distribution}
#'   \item{simpson_index}{Numeric: Simpson diversity index}
#'   \item{modal_frequency}{Numeric: frequency of modal structure}
#'
#' @details
#' Diversity metrics:
#' - n_unique: Higher = more diverse
#' - shannon_entropy: Higher = more diverse (0 = all same, log(M) = uniform)
#' - simpson_index: Higher = more diverse (0 = all same, 1-1/M = uniform)
#' - modal_frequency: Lower = more diverse
#'
#' @examples
#' \dontrun{
#' diversity <- analyze_structure_diversity(structures)
#' cat(sprintf("Unique structures: %d / %d\n",
#'             diversity$n_unique, length(structures)))
#' }
#'
#' @export
analyze_structure_diversity <- function(structures) {
  hashes <- vapply(structures, optimaltrees::structure_hash, character(1))
  counts <- table(hashes)
  M <- length(structures)

  # Number of unique structures
  n_unique <- length(counts)

  # Shannon entropy: H = -sum(p_i * log(p_i))
  probs <- as.numeric(counts) / M
  shannon_entropy <- -sum(probs * log(probs))

  # Simpson index: D = 1 - sum(p_i^2)
  simpson_index <- 1 - sum(probs^2)

  # Modal frequency
  modal_frequency <- max(counts) / M

  list(
    n_unique = n_unique,
    shannon_entropy = shannon_entropy,
    simpson_index = simpson_index,
    modal_frequency = modal_frequency
  )
}
