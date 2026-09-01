# ============================================================
# stage_b.R
# Study: threshold_superconsistency  (doubletree)
#
# STAGE B: off-grid threshold refinement.
#
# Stage A (the existing optimaltrees solver) can only ever return a cut that is
# one of the quantile grid points it discretized on, so its threshold error is
# floored at the grid mesh -- O(n^{-rho}) with rho = 1/3 at the coarse default,
# no matter how big n gets. Stage B removes that floor WITHOUT refining the grid:
# for each split, it identifies a bracket bounded ONLY by the nearest OTHER
# same-coordinate split in that node's own subtree (+/-Inf if none -- see
# `node_bracket()`), scans every OBSERVED value of that coordinate strictly
# inside the bracket, and keeps whichever exact value minimizes empirical risk.
# The bracket's job is IDENTIFICATION (isolating exactly one true boundary from
# any other), not localization; localization is the scan's job, and the
# attainable resolution is the spacing between order statistics, ~1/n -- which
# is the whole point of the check. The bracket is emphatically NOT "the two
# neighbouring GLOBAL GRID points": that bracket has width O(mesh) -> 0 as the
# grid refines with n, a shrinking target that the empirical/finite-sample
# comparison driving Stage A's grid-point choice does not track, and is exactly
# the mechanism behind the 2/24-reps "shrinking-margin" failure fixed on
# 2026-09-01 (see git history and `quality_reports/plans/
# 2026-09-01_two-stage-theory-session.md`, section 3).
#
# A true boundary interior to one grid atom forces a grid-optimal Stage-A tree
# to spend an EXTRA leaf splitting the SAME coordinate at BOTH grid cutpoints
# bracketing that atom (`prop:transition-leaves` in theory.tex) -- a thin
# spurious "transition leaf" sandwiched in between. `collapse_transitions()`
# MUST run on the coordinate-space tree before `refine_tree()` to merge each
# such pair into one node spanning both original grid cutpoints; otherwise the
# two cuts mutually bound each other at exactly one grid atom apart, silently
# re-creating an O(mesh) bracket inside the new rule.
#
# Nothing here refits or re-searches the tree TOPOLOGY. Structure comes from
# Stage A (modulo the transition-leaf collapse, which is a correction to a
# grid artifact, not a topology search); Stage B only moves thresholds. Leaf
# values are re-fit (they are just within-leaf means for squared-error loss),
# because a cut cannot be scored without them.
#
# The scan is exact, not approximate: for a fixed structure, sweeping the
# candidate cut in increasing order moves observations from the right subtree to
# the left subtree one distinct value at a time, so per-leaf (count, sum, sumsq)
# can be updated incrementally and total SSE recomputed in O(#leaves). Cost is
# O(|S| log|S| + #candidates) per node instead of O(|S| * #candidates); this is
# what keeps n = 32000 cheap and is verified against a naive implementation in
# check_stage_b.R.
#
# Sourced by run_cell.R. Not run directly.
# ============================================================

# ------------------------------------------------------------------
# 1. Read the fitted binary-feature tree back into coordinate space
# ------------------------------------------------------------------

#' Positional lookup from binary column index to (original coordinate, threshold)
#'
#' optimaltrees discretizes each continuous coordinate into indicator columns
#' 1{Xj <= tau_jk} and the fitted tree refers to them by a ZERO-BASED integer
#' index into that binary design matrix. This rebuilds the mapping.
#'
#' The rebuild is POSITIONAL because it mirrors exactly how the binary matrix is
#' assembled upstream (a loop over names(X), and within each coordinate a loop
#' over its thresholds, combined with cbind). It is then cross-checked against
#' metadata$binary_names, so a change in that assembly order fails loudly here
#' rather than silently mislabelling which coordinate a split is on.
#'
#' @param md A model's `@discretization_metadata`.
#' @return Data frame with one row per binary column: coord, k, cut.
bin_lookup <- function(md) {
  if (isTRUE(md$all_binary)) {
    cli::cli_abort("This study needs continuous covariates; metadata says all features are binary.")
  }
  rows <- list()
  for (cn in names(md$features)) {
    fm <- md$features[[cn]]
    if (identical(fm$type, "continuous")) {
      for (k in seq_along(fm$thresholds)) {
        rows[[length(rows) + 1L]] <- data.frame(
          coord = cn, k = k, cut = fm$thresholds[[k]],
          new_name = fm$new_names[[k]], stringsAsFactors = FALSE
        )
      }
    } else {
      rows[[length(rows) + 1L]] <- data.frame(
        coord = cn, k = NA_integer_, cut = NA_real_,
        new_name = fm$new_names[[1L]], stringsAsFactors = FALSE
      )
    }
  }
  lk <- do.call(rbind, rows)

  # Assert the positional reconstruction lines up with the recorded column names.
  if (nrow(lk) != length(md$binary_names)) {
    cli::cli_abort(paste0(
      "Rebuilt ", nrow(lk), " binary columns but metadata records ",
      length(md$binary_names), ". Binary-matrix assembly order has changed."
    ))
  }
  ok <- mapply(function(nm, suffix) endsWith(nm, suffix),
               md$binary_names, lk$new_name)
  if (!all(ok)) {
    bad <- which(!ok)[1L]
    cli::cli_abort(paste0(
      "Binary column ", bad, " is named '", md$binary_names[[bad]],
      "' but was reconstructed as '", lk$new_name[[bad]],
      "'. Binary-matrix assembly order has changed."
    ))
  }
  lk
}

#' Convert a fitted optimaltrees tree into a coordinate-space tree
#'
#' Nodes become list(kind, coord, k, cut, grid_cut, grid_cuts, collapsed, left,
#' right) where `left` is the {Xj <= cut} child. In the fitted object the split
#' is always `binary_feature == 1`, and the binary feature is the indicator
#' 1{Xj <= cut}, so the fitted `$true` branch is the {Xj <= cut} side.
#'
#' No bracket is stored here. The identifying bracket for Stage B is a
#' function of the CURRENT tree structure (nearest other same-coordinate cut),
#' not a static property of one node -- see `node_bracket()`. It also depends
#' on whether a node is later found to be one half of a transition-leaf pair
#' -- see `collapse_transitions()`, which must run before `refine_tree()`.
#'
#' @param node A node of `model@trees[[1]]`.
#' @param lk Output of `bin_lookup()`.
#' @param md A model's `@discretization_metadata`.
#' @return Nested list; leaves are `list(kind = "leaf")`.
as_coord_tree <- function(node, lk, md) {
  if (!is.null(node$prediction)) return(list(kind = "leaf"))
  if (is.null(node$feature)) {
    cli::cli_abort("Tree node has neither a prediction nor a feature.")
  }
  # These two are what make the branch orientation below correct; assert, do not assume.
  if (!identical(node$relation, "==")) {
    cli::cli_abort("Expected relation '==' on a binary feature, got '{node$relation}'.")
  }
  if (!identical(as.numeric(node$reference), 1)) {
    cli::cli_abort("Expected reference 1 on a binary feature, got {node$reference}.")
  }

  j <- as.integer(node$feature) + 1L   # fitted index is ZERO-based
  if (j < 1L || j > nrow(lk)) {
    cli::cli_abort("Split feature index {node$feature} is outside the binary design matrix.")
  }
  coord <- lk$coord[[j]]
  k     <- lk$k[[j]]
  if (is.na(k)) {
    cli::cli_abort("Split lands on a non-continuous column ('{coord}'); nothing to refine.")
  }
  taus <- md$features[[coord]]$thresholds

  list(
    kind      = "split",
    coord     = coord,
    k         = k,
    cut       = taus[[k]],
    grid_cut  = taus[[k]],
    grid_cuts = taus[[k]],
    collapsed = FALSE,
    left  = as_coord_tree(node$true,  lk, md),
    right = as_coord_tree(node$false, lk, md)
  )
}

# ------------------------------------------------------------------
# 2. Structure utilities
# ------------------------------------------------------------------

#' Number of leaves in a coordinate-space subtree
count_leaves <- function(node) {
  if (identical(node$kind, "leaf")) return(1L)
  count_leaves(node$left) + count_leaves(node$right)
}

#' Leaf index (1..count_leaves) for every row of X under a subtree
#'
#' Both children are evaluated on all supplied rows and then selected between.
#' Depth here is at most 2, so this costs a handful of full-length vector ops --
#' cheaper than the bookkeeping of recursive row subsetting.
leaf_index <- function(node, X) {
  if (identical(node$kind, "leaf")) return(rep(1L, nrow(X)))
  n_left <- count_leaves(node$left)
  ifelse(X[[node$coord]] <= node$cut,
         leaf_index(node$left, X),
         n_left + leaf_index(node$right, X))
}

#' Paths to every internal (split) node, root first
#'
#' A path is a character vector of "left"/"right" steps. Root is character(0).
split_paths <- function(node, path = character(0)) {
  if (identical(node$kind, "leaf")) return(list())
  c(list(path),
    split_paths(node$left,  c(path, "left")),
    split_paths(node$right, c(path, "right")))
}

#' Get / set a node at a path
node_at <- function(tree, path) {
  for (step in path) tree <- tree[[step]]
  tree
}

set_node_at <- function(tree, path, value) {
  if (length(path) == 0L) return(value)
  tree[[path[[1L]]]] <- set_node_at(tree[[path[[1L]]]], path[-1L], value)
  tree
}

#' Row indices reaching a node, given the CURRENT cuts
rows_at <- function(tree, X, path) {
  idx <- seq_len(nrow(X))
  node <- tree
  for (step in path) {
    keep <- if (identical(step, "left")) {
      X[[node$coord]][idx] <= node$cut
    } else {
      X[[node$coord]][idx] > node$cut
    }
    idx <- idx[keep]
    node <- node[[step]]
  }
  idx
}

#' All cuts on `coord` anywhere within a coordinate-space subtree
#'
#' Used by `node_bracket()` to find the nearest OTHER split on the same
#' coordinate inside a node's own children -- the only explicit bound Stage B
#' needs. Ancestor constraints on the same coordinate are already enforced by
#' `rows_at()` filtering the row set itself, so they must not be duplicated
#' here (see `node_bracket()`).
#'
#' @param node A coordinate-space (sub)tree.
#' @param coord Coordinate name to match.
#' @return Numeric vector of cut values (possibly length 0).
subtree_cuts_on <- function(node, coord) {
  if (identical(node$kind, "leaf")) return(numeric(0))
  c(if (identical(node$coord, coord)) node$cut else numeric(0),
    subtree_cuts_on(node$left,  coord),
    subtree_cuts_on(node$right, coord))
}

#' Identifying bracket for the split at `path`
#'
#' The nearest other same-coordinate cut in this node's OWN subtree, open at
#' both ends (+/-Inf if none). This bounds IDENTIFICATION -- isolating exactly
#' one true boundary from any other -- not LOCALIZATION, which is Stage B's
#' exact scan's job. Ancestor bounds are not computed here because `xj`
#' (passed to `bracket_candidates()` by `refine_tree()`, via `rows_at()`) is
#' already restricted to the ancestor-consistent range; adding an explicit
#' ancestor bound here would be redundant, not merely harmless, because
#' `rows_at()` uses whatever cut the ancestor CURRENTLY holds (grid before its
#' own refinement, refined after) and duplicating that logic risks disagreeing
#' with it.
#'
#' @param tree Full coordinate-space tree (current -- possibly partially
#'   refined by earlier, shallower calls in the same `refine_tree()` pass).
#' @param path Path to the split node.
#' @return List(lo, hi): an open interval, +/-Inf at either end if no
#'   bounding descendant cut exists on that side.
node_bracket <- function(tree, path) {
  node  <- node_at(tree, path)
  below <- subtree_cuts_on(node$left,  node$coord)
  above <- subtree_cuts_on(node$right, node$coord)
  list(lo = if (length(below)) max(below) else -Inf,
       hi = if (length(above)) min(above) else  Inf)
}

#' Detect and collapse transition-leaf pairs (`prop:transition-leaves`)
#'
#' A true boundary that lands interior to one grid atom forces a grid-optimal
#' Stage-A tree to split the SAME coordinate at both grid cutpoints
#' bracketing that atom, sandwiching a thin spurious leaf (the "transition
#' leaf") in between -- an extra leaf that is not estimating a second true
#' boundary, but half of the machinery for locating one. This function
#' detects that pattern and merges each pair into a single node spanning both
#' original grid cutpoints, so `node_bracket()` naturally treats the pair as
#' one identification problem with no special-casing anywhere else in Stage B.
#'
#' Detection is by GRID ADJACENCY (consecutive `k` indices on the same
#' coordinate, immediate parent/child, matching side) -- NOT by leaf mass,
#' which at ~n/n_bins observations is not "very low" and would both miss real
#' pairs and flag ordinary small leaves. Adjacency is exactly what
#' Assumption (Sep)'s separation floor plus mesh = o(separation) guarantees:
#' two same-coordinate cuts one grid atom apart cannot be two DISTINCT true
#' boundaries. If a DGP ever has two true boundaries on one coordinate closer
#' than one grid atom, this function is wrong by construction -- the
#' chained-pair abort below is the only tripwire for that.
#'
#' @param tree Coordinate-space tree from `as_coord_tree()`, UNREFINED (cuts
#'   are still exactly the grid values `as_coord_tree()` assigned).
#' @param max_iter Safety cap on collapse rounds (one collapse per round).
#' @return List(tree, n_collapses).
collapse_transitions <- function(tree, max_iter = 8L) {
  n_collapses <- 0L
  for (iter in seq_len(max_iter)) {
    paths <- split_paths(tree)
    did_collapse <- FALSE

    for (p in paths) {
      key <- paste(p, collapse = "/")
      P <- node_at(tree, p)
      j <- P$coord

      for (side in c("right", "left")) {
        C <- P[[side]]
        if (!identical(C$kind, "split") || !identical(C$coord, j)) next

        if (isTRUE(P$collapsed) || isTRUE(C$collapsed)) {
          cli::cli_abort(
            "Node '{key}' (or its '{side}' child) is already a collapsed transition pair and is itself adjacent to another same-coordinate split -- a chained (three-way) collapse, which is not supported. This is impossible at max_depth = 2; the depth cap must have changed."
          )
        }
        d <- C$k - P$k
        if (side == "right" && d == -1L) {
          cli::cli_abort(
            "Degenerate split: right child of '{key}' cuts BELOW its parent on coordinate '{j}' (k_parent = {P$k}, k_child = {C$k}). This means the {{Xj <= cut}} child is empty -- a bin_lookup()/branch-orientation bug, not a transition leaf."
          )
        }
        if (side == "left" && d == 1L) {
          cli::cli_abort(
            "Degenerate split: left child of '{key}' cuts ABOVE its parent on coordinate '{j}' (k_parent = {P$k}, k_child = {C$k}). This means the {{Xj > cut}} child is empty -- a bin_lookup()/branch-orientation bug, not a transition leaf."
          )
        }
        is_pair <- (side == "right" && d == 1L) || (side == "left" && d == -1L)
        if (!is_pair) next   # different, genuinely distinct boundary (Case 3 handles it)

        # The transition leaf sits between P's cut and C's cut: it is C's
        # LEFT child if C is P's right child (interval (tau_k, tau_{k+1}]),
        # or C's RIGHT child if C is P's left child (interval (tau_{k-1}, tau_k]).
        trans_side <- if (side == "right") "left" else "right"
        L <- C[[trans_side]]
        if (!identical(L$kind, "leaf")) {
          cli::cli_abort(
            "Transition child at '{key}/{side}/{trans_side}' is itself a split, not a leaf -- a chained (three-way) collapse, which is not supported. This is impossible at max_depth = 2; the depth cap must have changed."
          )
        }

        k_lo   <- min(P$k, C$k);   k_hi   <- max(P$k, C$k)
        tau_lo <- min(P$cut, C$cut); tau_hi <- max(P$cut, C$cut)
        collapsed_node <- list(
          kind      = "split", coord = j,
          k         = NA_integer_, k_lo = k_lo, k_hi = k_hi,
          grid_cuts = c(tau_lo, tau_hi),
          grid_cut  = mean(c(tau_lo, tau_hi)),   # Stage A's implied point estimate
          cut       = mean(c(tau_lo, tau_hi)),   # starting value only; Stage B moves it
          collapsed = TRUE,
          left  = if (side == "right") P$left  else C$left,
          right = if (side == "right") C$right else P$right
        )
        tree <- set_node_at(tree, p, collapsed_node)
        n_collapses  <- n_collapses + 1L
        did_collapse <- TRUE
        break
      }
      if (did_collapse) break   # re-derive split_paths() fresh before the next collapse
    }
    if (!did_collapse) return(list(tree = tree, n_collapses = n_collapses))
  }
  cli::cli_abort("collapse_transitions() did not reach a fixed point within {max_iter} iterations.")
}

#' Total training SSE of a coordinate-space tree with leaf means re-fit
#'
#' Used only for reporting and for the correctness self-check; the Stage B scan
#' itself uses the incremental form below.
tree_sse <- function(tree, X, y) {
  lid <- leaf_index(tree, X)
  K <- count_leaves(tree)
  cnt <- tabulate(lid, nbins = K)
  s <- numeric(K); ss <- numeric(K)
  for (k in which(cnt > 0L)) {
    yk <- y[lid == k]
    s[[k]]  <- sum(yk)
    ss[[k]] <- sum(yk * yk)
  }
  sum(ifelse(cnt > 0L, ss - s^2 / pmax(cnt, 1L), 0))
}

# ------------------------------------------------------------------
# 3. The Stage B scan
# ------------------------------------------------------------------

#' Exact empirical-risk-minimizing cutoff over observed values inside a bracket
#'
#' Sweeps the candidate cut upward. At each distinct candidate value all rows at
#' or below it belong to the left subtree, so rows migrate from the right
#' subtree's leaves to the left subtree's leaves monotonically and per-leaf
#' sufficient statistics update incrementally. SSE is then a K-term sum.
#'
#' @param xj Coordinate values of the rows reaching this node.
#' @param y Outcomes of those rows.
#' @param lid_left,lid_right Leaf id each row would take under the left / right
#'   subtree, on a common 1..K id space.
#' @param K Total number of leaf ids.
#' @param candidates Ascending candidate cutoffs (observed values in the bracket,
#'   plus the incumbent grid cut).
#' @return List(cut, sse, n_candidates).
scan_cutoff <- function(xj, y, lid_left, lid_right, K, candidates) {
  stopifnot(length(candidates) >= 1L)
  ord <- order(xj)
  xs  <- xj[ord]; ys <- y[ord]
  lL  <- lid_left[ord]; lR <- lid_right[ord]

  # Start below every candidate: all rows sit in the RIGHT subtree.
  cnt <- tabulate(lR, nbins = K)
  s   <- numeric(K); ss <- numeric(K)
  for (k in which(cnt > 0L)) {
    sel <- lR == k
    s[[k]]  <- sum(ys[sel]); ss[[k]] <- sum(ys[sel] * ys[sel])
  }

  sse_now <- function() sum(ifelse(cnt > 0L, ss - s^2 / pmax(cnt, 1L), 0))

  best_cut <- NA_real_; best_sse <- Inf
  pos <- 1L; m <- length(xs)
  for (cand in candidates) {
    # Migrate every row with xs <= cand from its right-leaf to its left-leaf.
    while (pos <= m && xs[[pos]] <= cand) {
      kR <- lR[[pos]]; kL <- lL[[pos]]; yv <- ys[[pos]]
      cnt[[kR]] <- cnt[[kR]] - 1L; s[[kR]] <- s[[kR]] - yv; ss[[kR]] <- ss[[kR]] - yv * yv
      cnt[[kL]] <- cnt[[kL]] + 1L; s[[kL]] <- s[[kL]] + yv; ss[[kL]] <- ss[[kL]] + yv * yv
      pos <- pos + 1L
    }
    v <- sse_now()
    if (v < best_sse) { best_sse <- v; best_cut <- cand }
  }
  list(cut = best_cut, sse = best_sse, n_candidates = length(candidates))
}

#' Candidate cutoffs for a node: observed values strictly inside the bracket
#'
#' `br` (from `node_bracket()`) is the identifying bracket: open at both
#' ends, strict inequalities. A candidate EQUAL to a bounding descendant cut
#' would empty one of that descendant's children, so it is excluded, not
#' merely redundant. The node's own grid cutpoint(s) (`grid_cuts` -- length 1
#' normally, length 2 if `collapsed`) are always added as candidates, clipped
#' to the bracket, so Stage B never scores worse than its Stage-A starting
#' point at an ORDINARY node (this guarantee is intentionally NOT claimed at
#' a `collapsed` node -- see `collapse_transitions()`'s docs).
bracket_candidates <- function(xj, node, br) {
  v  <- xj[xj > br$lo & xj < br$hi]
  gc <- node$grid_cuts
  gc <- gc[gc > br$lo & gc < br$hi]
  sort(unique(c(v, gc)))
}

#' Refine every split threshold of a coordinate-space tree, root first
#'
#' Top-down: the root is refined first, then children are refined on the row sets
#' induced by the already-refined ancestors. One pass -- this is a bracket
#' refinement, not an alternating-minimization search.
#'
#' Each node's bracket is recomputed from the CURRENT tree at the moment it is
#' refined (`node_bracket()`), not stored statically: refining an ancestor
#' first means a later same-coordinate descendant's row set (`rows_at()`) is
#' already filtered by the ancestor's REFINED cut, while a not-yet-refined
#' descendant node still contributes its GRID cut as the bounding value seen
#' by its parent. This asymmetry is monotone-safe and requires no second
#' pass, because after `collapse_transitions()` any two same-coordinate
#' splits that remain distinct nodes are separated by the (Sep) margin, not
#' by one grid atom -- do not "fix" this into an alternating minimization.
#'
#' @param tree Coordinate-space tree from `as_coord_tree()`, already passed
#'   through `collapse_transitions()`.
#' @param X,y Training data.
#' @return List(tree, refined) where `refined` is a data frame with one row
#'   per split: path, coord, grid_cut, refined_cut, lo, hi, cand_lo, cand_hi,
#'   collapsed, n_node, n_candidates.
refine_tree <- function(tree, X, y) {
  paths <- split_paths(tree)
  # Root-first order: shorter paths first.
  paths <- paths[order(vapply(paths, length, integer(1)))]

  log_rows <- list()
  for (p in paths) {
    idx  <- rows_at(tree, X, p)
    node <- node_at(tree, p)
    br   <- node_bracket(tree, p)
    cands <- if (length(idx) < 2L) numeric(0) else {
      bracket_candidates(X[[node$coord]][idx], node, br)
    }

    if (length(cands) == 0L) {
      log_rows[[length(log_rows) + 1L]] <- data.frame(
        path = paste(p, collapse = "/"), coord = node$coord,
        grid_cut = node$grid_cut, refined_cut = node$cut,
        lo = br$lo, hi = br$hi, cand_lo = NA_real_, cand_hi = NA_real_,
        collapsed = isTRUE(node$collapsed),
        n_node = length(idx), n_candidates = 0L, stringsAsFactors = FALSE
      )
      next
    }

    Xn <- X[idx, , drop = FALSE]
    yn <- y[idx]
    xj <- Xn[[node$coord]]

    n_left <- count_leaves(node$left)
    lid_left  <- leaf_index(node$left,  Xn)
    lid_right <- n_left + leaf_index(node$right, Xn)
    K <- n_left + count_leaves(node$right)

    res <- scan_cutoff(xj, yn, lid_left, lid_right, K, cands)

    log_rows[[length(log_rows) + 1L]] <- data.frame(
      path = paste(p, collapse = "/"), coord = node$coord,
      grid_cut = node$grid_cut, refined_cut = res$cut,
      lo = br$lo, hi = br$hi,
      cand_lo = min(cands), cand_hi = max(cands),
      collapsed = isTRUE(node$collapsed),
      n_node = length(idx), n_candidates = res$n_candidates,
      stringsAsFactors = FALSE
    )

    node$cut <- res$cut
    tree <- set_node_at(tree, p, node)
  }

  list(tree = tree, refined = do.call(rbind, log_rows))
}
