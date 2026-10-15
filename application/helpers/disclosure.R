## ============================================================================
## application/helpers/disclosure.R
##
## Small-cell helpers shared by the decision report (08) and the Tier-1 exports
## (90). smidata::smi_suppress() protects each exported cell but NOT differences
## between cells or tables, so any count that is the difference of two published
## counts must be handled here: it is hidden by removing a count, not masked.
## ============================================================================

config_report_min_cell <- 11L

#' Which counts of a nested chain must be hidden so no small drop is recoverable.
#'
#' A chain `n[1] >= n[2] >= ...` (attrition steps, a total and its subset) publishes
#' every drop `n[i - 1] - n[i]` implicitly. A drop of 1..(min_cell - 1) patients is
#' a small cell, and masking the drop itself does nothing while both neighbours
#' are shown. The inner count is hidden instead (the last count is kept unless the
#' chain has only two visible counts), and the gaps between the remaining neighbours
#' are re-checked until none is small.
#'
#' @param n Integer counts in chain order (`NA` entries are ignored).
#' @param min_cell Smallest publishable cell.
#' @return A logical vector, `TRUE` where the count must not be published.
hide_small_gaps <- function(n, min_cell = config_report_min_cell) {
  hidden <- rep(FALSE, length(n))
  repeat {
    vis <- which(!hidden & !is.na(n))
    if (length(vis) < 2L) break
    gaps <- n[vis[-length(vis)]] - n[vis[-1L]]
    bad <- which(gaps > 0 & gaps < min_cell)
    if (length(bad) == 0L) break
    k <- bad[[1L]]
    hidden[if (k + 1L < length(vis)) vis[k + 1L] else if (k > 1L) vis[k] else vis[k + 1L]] <- TRUE
  }
  hidden
}

#' Attrition steps with the hidden counts merged into the next published step.
#'
#' @param flow Data frame with `step`, `n`, `reason`.
#' @param min_cell Smallest publishable cell.
#' @return `flow` without the rows `hide_small_gaps()` hides; a removed row's step
#'   and reason are appended to the step that absorbs it.
coarsen_flow <- function(flow, min_cell = config_report_min_cell) {
  hidden <- hide_small_gaps(flow$n, min_cell)
  out <- flow[!hidden, , drop = FALSE]
  kept <- which(!hidden)
  ids <- which(hidden)
  ## Rows with a later kept row are absorbed back to front, so the merged
  ## reasons read in step order; trailing rows are appended front to back.
  for (i in c(rev(ids[ids < max(kept)]), ids[ids > max(kept)])) {
    later <- kept[kept > i]
    to <- if (length(later) > 0L) later[[1L]] else max(kept[kept < i])
    row <- match(to, kept)
    for (col in c("step", "reason")) {
      parts <- if (to > i) c(flow[[col]][i], out[[col]][row]) else c(out[[col]][row], flow[[col]][i])
      out[[col]][row] <- paste(parts[nzchar(parts)], collapse = "; ")
    }
  }
  out
}


#' Mask a partition whose total is published.
#'
#' If any cell is a small positive count, it is masked together with the
#' smallest other positive cell, so the published total leaves two unknowns.
#'
#' @param cells Integer counts that sum to a published total.
#' @param min_cell Smallest publishable cell.
#' @return `cells` with masked entries `NA`.
mask_partition_cells <- function(cells, min_cell = config_report_min_cell) {
  small <- !is.na(cells) & cells > 0 & cells < min_cell
  if (!any(small)) return(cells)
  hide <- small
  if (sum(small) == 1L) {
    other <- which(!is.na(cells) & !small & cells > 0)
    if (length(other) > 0L) hide[other[[which.min(cells[other])]]] <- TRUE
  }
  cells[hide] <- NA
  cells
}


## Prevalence counts with every (column, arm) partition masked so its published
## arm size cannot give a masked cell back.
tier1_disclosable_prevalence <- function(prevalence) {
  out <- prevalence
  for (arm in c("control", "treated")) {
    cols <- paste0(c("n_yes_", "n_no_", "n_na_"), arm)
    masked <- t(apply(as.matrix(out[cols]), 1L, mask_partition_cells))
    out[cols] <- as.data.frame(masked)
  }
  out$prevalence_control <- out$n_yes_control / (out$n_yes_control + out$n_no_control)
  out$prevalence_treated <- out$n_yes_treated / (out$n_yes_treated + out$n_no_treated)
  out
}


## The exact-start sizes in a form whose small cells are not recoverable: a small
## step of the with_rows >= matched >= elig chain hides the inner count, and if
## elig is hidden (or any arm cell is small) the arm cells that would give it back
## are masked together (`mask_partition_cells()`).
tier1_disclosable_sizes <- function(by_band) {
  out <- by_band
  chain <- c("n_with_saa_rows", "n_matched", "n_elig")
  hid <- hide_small_gaps(unlist(out[chain], use.names = FALSE))
  out[chain[hid]] <- NA_integer_
  arms <- c("n_treated", "n_control", "n_achieved_na")
  cells <- unlist(out[arms], use.names = FALSE)
  if (is.na(out$n_elig)) {
    out[arms] <- NA_integer_            # their sum would give elig back
  } else {
    masked <- mask_partition_cells(cells)
    out[arms] <- as.list(masked)
    if (anyNA(masked)) out$n_elig <- NA_integer_   # elig - visible arms = masked arm
  }
  out
}


