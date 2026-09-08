# ---------------------------------------------------------------------------
# DGP for the T_a / T_epsilon remainder-term study
#
# Design goals (see ../SPEC.md):
#   * p = 5 continuous covariates, Z_j ~ iid Uniform(0, 1).
#   * The analyst's row-n feature map X_n = Q_n(Z) is a FIXED axis-aligned grid
#     with `n_thresh` thresholds per coordinate. Everything the analyst can fit
#     is piecewise constant on the resulting grid cells.
#   * The TRUE nuisances e_0 and mu_0 are EXACTLY representable as axis-aligned
#     trees with exactly `L` leaves whose split thresholds lie ON that grid.
#     Hence mu_{0,n} = mu_0 and e_{0,n} = e_0: the mesh bias h_n is exactly 0
#     and we are in the well-specified regime the theory asks about.
#   * Moderate SNR: outcome leaf values are spaced ~0.5-1.5 noise SDs apart.
#
# Representation. A leaf is an axis-aligned box recorded as integer interval
# index ranges over the `n_thresh + 1` grid intervals of each coordinate:
# `lo[j]:hi[j]` inclusive. Boxes are therefore exact unions of grid cells, which
# makes every population integral in terms.R an exact finite sum.
# ---------------------------------------------------------------------------

P_COVARIATES <- as.integer(getOption("rtvr.p", 5L))
# Thresholds per coordinate in the analyst's row-n grid. Settable before
# sourcing via options(rtvr.n_thresh = k) so the feasibility probe can sweep it;
# the study fixes one value (see ../SPEC.md).
N_THRESH <- as.integer(getOption("rtvr.n_thresh", 2L))
GRID_THRESHOLDS <- seq_len(N_THRESH) / (N_THRESH + 1L)   # k = 2 -> 1/3, 2/3
N_INTERVAL <- N_THRESH + 1L                     # 4 intervals per coordinate
N_GRID_CELL <- N_INTERVAL^P_COVARIATES          # 1024 equal-mass cells

MU_LEAF_VALUES <- 0.75 * c(-1, -1 / 3, 1 / 3, 1)     # {-0.75,-0.25,0.25,0.75}
E_LEAF_VALUES <- c(0.30, 0.45, 0.55, 0.70)           # bounded away from 0/1
TAU_ATT <- 1                                          # constant effect => ATT = 1
NOISE_SD <- 1

# Coordinate priority orders: different for the two nuisances so that the true
# propensity partition and the true outcome partition genuinely differ (a shared
# partition would make the within-leaf centring in b_i degenerate).
COORD_ORDER_MU <- c(1L, 2L, 3L, 4L, 5L)
COORD_ORDER_E <- c(3L, 4L, 5L, 1L, 2L)


#' Build an axis-aligned partition of [0,1]^p with exactly `L` leaves
#'
#' Deterministic construction: repeatedly split the largest-volume leaf, on the
#' coordinate with the widest remaining interval span (ties broken by
#' `coord_order`), at the median available grid threshold inside that span. This
#' yields a near-balanced tree of depth ceiling(log2(L)) with every leaf of
#' non-negligible mass, so no leaf has vanishing control mass.
#'
#' @return list of leaves; each leaf is `list(lo = int[p], hi = int[p])` giving
#'   inclusive grid-interval index ranges.
build_partition <- function(L, coord_order) {
  stopifnot(L >= 1L, L <= N_INTERVAL^P_COVARIATES)
  leaves <- list(list(lo = rep(1L, P_COVARIATES), hi = rep(N_INTERVAL, P_COVARIATES)))

  while (length(leaves) < L) {
    spans <- vapply(leaves, function(lf) prod(lf$hi - lf$lo + 1L), numeric(1))
    splittable <- vapply(leaves, function(lf) any(lf$hi > lf$lo), logical(1))
    stopifnot("build_partition: no splittable leaf left" = any(splittable))
    spans[!splittable] <- -Inf
    target <- which.max(spans)
    lf <- leaves[[target]]

    widths <- lf$hi - lf$lo
    widths[widths == 0L] <- -1L
    # rank coordinates by (width desc, coord_order position asc)
    ord_pos <- match(seq_len(P_COVARIATES), coord_order)
    j <- order(-widths, ord_pos)[1L]
    stopifnot(widths[j] > 0L)

    # thresholds strictly inside the leaf's span on coordinate j are indexed
    # lo[j] .. hi[j]-1 (threshold k separates interval k from interval k+1)
    cut_candidates <- seq.int(lf$lo[j], lf$hi[j] - 1L)
    k <- cut_candidates[ceiling(length(cut_candidates) / 2)]

    left <- lf
    left$hi[j] <- k
    right <- lf
    right$lo[j] <- k + 1L
    leaves[[target]] <- left
    leaves[[length(leaves) + 1L]] <- right
  }

  stopifnot(length(leaves) == L)
  leaves
}


#' Enumerate the grid cells: one row per cell, columns = interval index per coord
grid_cell_index <- function() {
  as.matrix(expand.grid(rep(list(seq_len(N_INTERVAL)), P_COVARIATES),
                        KEEP.OUT.ATTRS = FALSE))
}

#' Minimum tree depth needed to represent a partition
#'
#' Each leaf box needs one split per finite face, so the number of finite faces
#' on the deepest leaf is the depth the fitter must be allowed to search.
partition_depth <- function(leaves) {
  max(vapply(leaves, function(lf) sum((lf$lo > 1L) + (lf$hi < N_INTERVAL)),
             numeric(1)))
}

#' A representative Z for each grid cell (interval midpoints)
grid_cell_centres <- function(cell_idx = grid_cell_index()) {
  edges <- c(0, GRID_THRESHOLDS, 1)
  Zc <- matrix(0, nrow = nrow(cell_idx), ncol = P_COVARIATES)
  for (j in seq_len(P_COVARIATES)) {
    Zc[, j] <- (edges[cell_idx[, j]] + edges[cell_idx[, j] + 1L]) / 2
  }
  colnames(Zc) <- paste0("z", seq_len(P_COVARIATES))
  Zc
}

#' Map each grid cell to the index of the leaf containing it
cells_to_leaf <- function(leaves, cell_idx) {
  out <- integer(nrow(cell_idx))
  for (l in seq_along(leaves)) {
    lf <- leaves[[l]]
    inside <- rep(TRUE, nrow(cell_idx))
    for (j in seq_len(P_COVARIATES)) {
      inside <- inside & cell_idx[, j] >= lf$lo[j] & cell_idx[, j] <= lf$hi[j]
    }
    out[inside] <- l
  }
  stopifnot("cells_to_leaf: unassigned grid cell" = all(out > 0L))
  out
}


#' Assemble the full truth for a given leaf budget
#'
#' @return list with the two true partitions, per-grid-cell true e_0 and mu_0,
#'   equal cell masses, and the exact ATT.
make_truth <- function(L) {
  cell_idx <- grid_cell_index()
  leaves_mu <- build_partition(L, COORD_ORDER_MU)
  leaves_e <- build_partition(L, COORD_ORDER_E)

  cell_leaf_mu <- cells_to_leaf(leaves_mu, cell_idx)
  cell_leaf_e <- cells_to_leaf(leaves_e, cell_idx)

  mu_by_leaf <- MU_LEAF_VALUES[((seq_len(L) - 1L) %% length(MU_LEAF_VALUES)) + 1L]
  e_by_leaf <- E_LEAF_VALUES[((seq_len(L) - 1L) %% length(E_LEAF_VALUES)) + 1L]

  cell_mu0 <- mu_by_leaf[cell_leaf_mu]
  cell_e0 <- e_by_leaf[cell_leaf_e]
  cell_mass <- rep(1 / N_GRID_CELL, N_GRID_CELL)

  # constant treatment effect => ATT = TAU_ATT exactly, whatever the propensity
  list(
    L = L,
    cell_idx = cell_idx,
    cell_centres = grid_cell_centres(cell_idx),
    cell_mass = cell_mass,
    cell_mu0 = cell_mu0,
    cell_e0 = cell_e0,
    cell_w0 = cell_e0 / (1 - cell_e0),
    pi_true = sum(cell_mass * cell_e0),
    theta0 = TAU_ATT,
    leaves_mu = leaves_mu,
    leaves_e = leaves_e,
    depth_mu = partition_depth(leaves_mu),
    depth_e = partition_depth(leaves_e)
  )
}


#' Binarise continuous Z into the analyst's row-n features X_n = Q_n(Z)
#'
#' Column (j - 1) * N_THRESH + k is 1{Z_j > GRID_THRESHOLDS[k]}.
binarise <- function(Z) {
  n <- nrow(Z)
  out <- matrix(0L, nrow = n, ncol = P_COVARIATES * N_THRESH)
  nm <- character(ncol(out))
  for (j in seq_len(P_COVARIATES)) {
    for (k in seq_len(N_THRESH)) {
      col <- (j - 1L) * N_THRESH + k
      out[, col] <- as.integer(Z[, j] > GRID_THRESHOLDS[k])
      nm[col] <- sprintf("z%d_gt%02d", j, round(100 * GRID_THRESHOLDS[k]))
    }
  }
  colnames(out) <- nm
  out
}

#' Interval index (1..N_INTERVAL) of each coordinate value
interval_of <- function(Z) {
  idx <- matrix(1L, nrow = nrow(Z), ncol = P_COVARIATES)
  for (j in seq_len(P_COVARIATES)) {
    idx[, j] <- 1L + rowSums(outer(Z[, j], GRID_THRESHOLDS, ">"))
  }
  idx
}

#' Row index into the grid-cell table for each observation
cell_id_of <- function(Z) {
  idx <- interval_of(Z)
  # expand.grid varies the FIRST factor fastest
  mult <- N_INTERVAL^(seq_len(P_COVARIATES) - 1L)
  as.integer(1L + (idx - 1L) %*% mult)
}


#' Draw one data set
#'
#' @return list(Z, Xbin, cell, A, Y, e0, mu0, w0, eps)
simulate_data <- function(n, truth) {
  Z <- matrix(stats::runif(n * P_COVARIATES), nrow = n, ncol = P_COVARIATES)
  colnames(Z) <- paste0("z", seq_len(P_COVARIATES))
  cell <- cell_id_of(Z)

  e0 <- truth$cell_e0[cell]
  mu0 <- truth$cell_mu0[cell]
  A <- stats::rbinom(n, 1L, e0)
  # bounded noise (variance 1) so the theory's B_Y moment bound holds literally
  eps <- stats::runif(n, -sqrt(3), sqrt(3)) * NOISE_SD
  Y <- mu0 + TAU_ATT * A + eps

  list(Z = Z, Xbin = binarise(Z), cell = cell, A = A, Y = Y,
       e0 = e0, mu0 = mu0, w0 = e0 / (1 - e0), eps = eps)
}
