# ============================================================
# enumerate_sufficient_class.R
# Study: partition_recovery_clt  (doubletree, S1)
# Spec:  quality_reports/specs/2026-09-01_partition-recovery-clt.md  (see its §2)
#
# SHARED INFRASTRUCTURE. S2 (2026-08-21_honest-inference-sparsity-failure.md)
# imports this file to define its own delta_e/delta_mu dial rather than
# re-deriving a second, possibly inconsistent notion of "population
# approximation error". Keep it dependency-light and side-effect-free: sourcing
# this file defines functions and nothing else.
#
# WHAT IT PROVIDES, for a binary-covariate DGP (so T_Lbar is finite by
# construction, manuscript.tex ass:finite):
#
#   enumerate_tree_partitions()  every partition of {0,1}^p realizable by a
#                                decision tree with <= Lbar leaves, deduplicated
#                                as PARTITIONS (not as tree topologies)
#   sufficient_class()           exact membership in S_j (def:sufficient) and the
#                                exact margin Delta_j (prop:selection-rate), both
#                                from closed-form population probabilities --
#                                never by simulation
#   fitted_cell_labels()         a FITTED optimaltrees model -> its partition of
#                                the covariate space, evaluated on all 2^p cells
#   partition_recovered()        the currency-correct recovery indicator:
#                                Pi_tauhat gamma_0 == gamma_0 exactly
#
# THREE DESIGN POINTS THAT ARE EASY TO GET WRONG, STATED UP FRONT
#
# (1) Recovery is checked from the DEFINITION, not by lookup. Condition (P1)'s
#     currency is "gamma_0 is constant on every leaf of tauhat_j", i.e.
#     Pi_tauhat gamma_0 = gamma_0. It is NOT "tauhat equals one fixed reference
#     partition" and NOT "|tauhat| equals the reference leaf count". The
#     manuscript's own remark after cond:P is explicit that S_j need not be a
#     singleton: any refinement of a sufficient partition -- including one that
#     spends a leaf on a pure-noise coordinate -- is still sufficient. A
#     same-partition or same-leaf-count check therefore scores genuine
#     recoveries as failures. (This is exactly what
#     propensity_loss_choice/code/analyze.R's `same_partition_rate_e` measures,
#     and why that number is NOT a (P1) recovery rate: it is agreement between
#     two ARMS, not membership in S_j.) partition_recovered() below evaluates
#     the definition directly, so it is correct even for a fitted partition that
#     is (unexpectedly) outside the enumerated class; `in_enumerated_class` is
#     returned separately as a consistency diagnostic.
#
# (2) The check must run on the CELL GRID, not on the sample. tauhat partitions
#     the whole covariate space, and Pi_tauhat is a population projection. At
#     n = 200 with p = 5 several of the 32 cells are typically empty in a given
#     sample; a sample-based check would silently skip them and could call a
#     non-sufficient partition sufficient. Evaluating the fitted tree on all 2^p
#     patterns is both well-defined and the only faithful check.
#
# (3) TWO margins are reported, on purpose. manuscript.tex defines R^{(j)} as a
#     SQUARED-ERROR population risk (see the display after eq:risk-emp:
#     R^{(e)}(tau) = E[{A - Pi_tau e_0}^2],
#     R^{(mu)}(tau) = E[{Y - Pi^nu_tau mu_0}^2 | A=0]), so Delta_j in
#     prop:selection-rate is a squared-error gap: `delta_sq`. But
#     doubletree::estimate_att() SELECTS with log-loss by default (both for e,
#     via propensity_loss, and for mu, via outcome_type = "binary"), so the
#     margin that actually governs empirical selection is the excess-log-loss
#     (KL) gap: `delta_kl`. Both are computed. S_j itself is loss-free -- it is
#     defined by Pi^nu_tau gamma_0 = gamma_0 for a strictly positive weight nu,
#     which does not mention a loss -- so only the MARGIN, not the membership,
#     depends on this choice.
# ============================================================

# ---- 1 the covariate grid ------------------------------------------------

#' All 2^p binary covariate patterns
#'
#' Row order is fixed and reproducible (\code{expand.grid} with the first
#' variable varying fastest); every other function here indexes cells by row
#' number, so this order is the study's canonical cell ordering.
#'
#' @param p Number of binary covariates.
#' @param var_names Column names. Default \code{X1..Xp}.
#' @return A data.frame with \code{2^p} rows and \code{p} integer 0/1 columns.
make_cell_grid <- function(p, var_names = paste0("X", seq_len(p))) {
  stopifnot(length(p) == 1L, p >= 1L, p == as.integer(p))
  stopifnot(length(var_names) == p, !anyDuplicated(var_names))
  args <- rep(list(c(0L, 1L)), p)
  names(args) <- var_names
  grid <- expand.grid(args, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  # expand.grid() already returns integers here; as.data.frame keeps it a plain
  # data.frame (not a tibble) because optimaltrees indexes it with X[i, k].
  as.data.frame(lapply(grid, as.integer), stringsAsFactors = FALSE)
}

# ---- 2 canonical partition keys -----------------------------------------

#' Relabel a leaf-label vector canonically
#'
#' Two label vectors describe the same partition iff their canonical forms are
#' identical. Labels are renumbered 1, 2, ... in order of first appearance.
#'
#' @param labels Integer vector, one label per cell.
#' @return Integer vector of the same length.
canonical_labels <- function(labels) {
  match(labels, unique(labels))
}

#' Canonical string key for a partition
#'
#' @param labels Integer vector, one label per cell.
#' @return A single character string.
partition_key <- function(labels) {
  paste(canonical_labels(labels), collapse = ".")
}

# ---- 3 enumerate T_Lbar --------------------------------------------------

#' Every partition of \{0,1\}^p realizable by a tree with at most Lbar leaves
#'
#' Enumerates tree TOPOLOGIES recursively (at each node: stop, or split on any
#' coordinate not already fixed on the path) and then deduplicates the resulting
#' PARTITIONS. Deduplication is not cosmetic: splitting X1 then X2 in both
#' children and splitting X2 then X1 in both children are different topologies
#' realizing the identical 4-leaf partition, and \code{|T_Lbar|} in
#' \code{prop:selection-rate}'s bound counts partitions, not topologies.
#'
#' A split that would send every cell of the current region to one side is
#' skipped: on binary covariates that only happens when the coordinate is
#' already fixed on the path, and the resulting empty leaf is not a partition
#' cell.
#'
#' @param p Number of binary covariates.
#' @param leaf_budget Lbar, the maximum number of leaves.
#' @param var_names Covariate names, matching \code{make_cell_grid()}.
#' @return A list with
#'   \item{cells}{the \code{2^p x p} cell grid (data.frame)}
#'   \item{labels}{integer matrix, one ROW per partition, one COLUMN per cell}
#'   \item{keys}{character vector of canonical partition keys}
#'   \item{n_leaves}{integer vector of leaf counts}
#'   \item{n_partitions}{\code{|T_Lbar|}}
enumerate_tree_partitions <- function(p, leaf_budget,
                                      var_names = paste0("X", seq_len(p))) {
  stopifnot(length(leaf_budget) == 1L, leaf_budget >= 1L,
            leaf_budget == as.integer(leaf_budget))
  leaf_budget <- as.integer(leaf_budget)
  cells <- make_cell_grid(p, var_names)
  cell_mat <- as.matrix(cells)
  n_cells <- nrow(cell_mat)

  # Memoised recursion. The cache key is (region signature, budget); the same
  # (region, budget) subproblem recurs once per coordinate ordering that reaches
  # it, which without memoisation is ~10^5 redundant subproblems at p = 6.
  cache <- new.env(parent = emptyenv())

  enum <- function(region, free_coords, budget) {
    ck <- paste0(paste(region, collapse = ","), "|", budget)
    hit <- cache[[ck]]
    if (!is.null(hit)) return(hit)

    out <- list(list(region))                      # option: stop, one leaf
    if (budget >= 2L) {
      for (k in free_coords) {
        left <- region[cell_mat[region, k] == 0L]
        right <- region[cell_mat[region, k] == 1L]
        if (length(left) == 0L || length(right) == 0L) next   # coord already fixed
        rem <- free_coords[free_coords != k]
        subs_l <- enum(left, rem, budget - 1L)
        subs_r <- enum(right, rem, budget - 1L)
        for (a in subs_l) {
          for (b in subs_r) {
            if (length(a) + length(b) <= budget) {
              out[[length(out) + 1L]] <- c(a, b)
            }
          }
        }
      }
    }
    cache[[ck]] <- out
    out
  }

  raw <- enum(seq_len(n_cells), seq_len(p), leaf_budget)

  # Topologies -> label vectors -> dedupe by canonical key.
  keys <- character(length(raw))
  labs <- vector("list", length(raw))
  for (i in seq_along(raw)) {
    lab <- integer(n_cells)
    for (l in seq_along(raw[[i]])) lab[raw[[i]][[l]]] <- l
    lab <- canonical_labels(lab)
    labs[[i]] <- lab
    keys[[i]] <- paste(lab, collapse = ".")
  }
  keep <- !duplicated(keys)
  labs <- labs[keep]
  keys <- keys[keep]

  label_mat <- do.call(rbind, labs)
  list(
    cells = cells,
    labels = label_mat,
    keys = keys,
    n_leaves = apply(label_mat, 1L, function(z) length(unique(z))),
    n_partitions = length(keys),
    p = p,
    leaf_budget = leaf_budget,
    var_names = var_names
  )
}

# ---- 4 population risks, exactly ----------------------------------------

#' Weighted leaf means of gamma_0 under a partition
#'
#' These are exactly the projections Pi^nu_tau gamma_0 of manuscript.tex: the
#' nu-weighted leaf mean minimises BOTH squared error and log-loss within a
#' leaf, so one set of leaf means serves both margins.
#'
#' @param labels Integer leaf label per cell.
#' @param gamma0 True nuisance value per cell.
#' @param w Non-negative cell weights (need not sum to 1).
#' @return Numeric vector of length \code{length(gamma0)}: the projected value
#'   at each cell.
project_on_partition <- function(labels, gamma0, w) {
  stopifnot(length(labels) == length(gamma0), length(w) == length(gamma0))
  num <- tapply(w * gamma0, labels, sum)
  den <- tapply(w, labels, sum)
  # A leaf with zero total weight cannot influence any population quantity; its
  # projected value is undefined, so carry the truth through unchanged rather
  # than manufacturing an NaN that would poison the excess-risk sums.
  ratio <- ifelse(den > 0, num / den, NA_real_)
  out <- as.numeric(ratio[as.character(labels)])
  zero_mass <- is.na(out)
  out[zero_mass] <- gamma0[zero_mass]
  out
}

#' Excess population risk of a partition, both currencies
#'
#' Excess risk relative to the unrestricted Bayes risk. For squared error this
#' is \code{E_w[(gamma_0 - Pi_tau gamma_0)^2]}; for log-loss it is
#' \code{E_w[KL(gamma_0 || Pi_tau gamma_0)]}. Both vanish exactly when
#' \code{Pi_tau gamma_0 = gamma_0}, which is why either can be used to define
#' the margin -- and why the margin's NUMERIC VALUE differs between them.
#'
#' @param labels Integer leaf label per cell.
#' @param gamma0 True nuisance value per cell, in (0, 1) for the KL currency.
#' @param w Cell weights; normalised internally.
#' @return Named numeric vector \code{c(sq = , kl = )}.
excess_risk <- function(labels, gamma0, w) {
  wn <- w / sum(w)
  gbar <- project_on_partition(labels, gamma0, w)
  sq <- sum(wn * (gamma0 - gbar)^2)
  kl <- if (all(gamma0 > 0 & gamma0 < 1) && all(gbar > 0 & gbar < 1)) {
    sum(wn * (gamma0 * log(gamma0 / gbar) +
                (1 - gamma0) * log((1 - gamma0) / (1 - gbar))))
  } else {
    NA_real_
  }
  # Both quantities are >= 0 exactly (Jensen; gbar is the w-weighted leaf mean).
  # On a SUFFICIENT partition every term is an exact zero mathematically, but the
  # floating-point sum lands at +/- 1e-18, and a negative KL in a results table
  # reads as a bug rather than as cancellation noise. Clamped at 0, which cannot
  # mask a real positive excess: the smallest genuine margin in this study is
  # ~1e-4, fourteen orders of magnitude above the noise floor.
  c(sq = max(0, sq), kl = if (is.na(kl)) NA_real_ else max(0, kl))
}

#' Is gamma_0 exactly constant on every leaf of this partition?
#'
#' The STRUCTURAL criterion, i.e. def:sufficient read literally. Used in
#' preference to "excess risk == 0" because it is a comparison of DGP constants
#' rather than of a floating-point sum, and because it stays meaningful for
#' leaves of zero population mass.
#'
#' @param labels Integer leaf label per cell.
#' @param gamma0 True nuisance value per cell.
#' @param tol Numerical slack on the within-leaf range of \code{gamma0}. The
#'   caller is expected to check (via \code{sufficient_class()}'s
#'   \code{min_value_gap}) that distinct true values are separated by orders of
#'   magnitude more than this.
#' @return \code{TRUE}/\code{FALSE}.
is_sufficient_partition <- function(labels, gamma0, tol = 1e-12) {
  rng <- tapply(gamma0, labels, function(z) diff(range(z)))
  all(rng <= tol)
}

# ---- 5 the sufficient class and the margin ------------------------------

#' Exact sufficient class S_j and margin Delta_j for one nuisance
#'
#' Everything here is closed-form in the DGP's own population probabilities.
#' No data is generated and no model is fitted.
#'
#' @param gamma0 True nuisance value at each cell of \code{enumeration$cells},
#'   in the canonical cell order.
#' @param w Population cell weights for this nuisance's risk. For the
#'   propensity this is \code{P(X = x)}; for the control outcome it is the
#'   CONTROL-weighted measure \code{nu_x = P(X = x){1 - e_0(x)}}, because
#'   \code{R^{(mu)}} is an expectation given \code{A = 0} and
#'   \code{estimate_att()} fits the mu tree on controls only.
#' @param enumeration Output of \code{enumerate_tree_partitions()}.
#' @param tol Structural tolerance forwarded to \code{is_sufficient_partition()}.
#' @return A list with \code{sufficient} (logical, per partition),
#'   \code{excess_sq}/\code{excess_kl}, \code{delta_sq}/\code{delta_kl},
#'   \code{n_sufficient}, \code{min_sufficient_leaves}, \code{argmin_key_*}
#'   (the tightest non-sufficient competitor), \code{min_value_gap}, and the
#'   \code{enumeration} it was built against.
sufficient_class <- function(gamma0, w, enumeration, tol = 1e-12) {
  n_cells <- nrow(enumeration$cells)
  stopifnot(length(gamma0) == n_cells, length(w) == n_cells)
  stopifnot(all(is.finite(gamma0)), all(is.finite(w)), all(w >= 0), sum(w) > 0)

  labels_mat <- enumeration$labels
  n_part <- nrow(labels_mat)
  suff <- logical(n_part)
  ex_sq <- numeric(n_part)
  ex_kl <- numeric(n_part)
  for (i in seq_len(n_part)) {
    lab <- labels_mat[i, ]
    suff[[i]] <- is_sufficient_partition(lab, gamma0, tol = tol)
    e <- excess_risk(lab, gamma0, w)
    ex_sq[[i]] <- e[["sq"]]
    ex_kl[[i]] <- e[["kl"]]
  }

  if (!any(suff)) {
    stop("No partition in T_Lbar is sufficient: S_j is EMPTY at this ",
         "leaf_budget, so Instantiation 1's precondition fails and this DGP ",
         "cannot be used in this study. Raise leaf_budget or simplify the DGP.",
         call. = FALSE)
  }
  if (all(suff)) {
    stop("EVERY partition in T_Lbar is sufficient, so Delta_j is undefined ",
         "(there is no competitor outside S_j). gamma_0 is probably constant.",
         call. = FALSE)
  }

  # Consistency: the structural criterion and the risk criterion must agree.
  # A mismatch means tol is wrong for this DGP's value spacing.
  max_suff_sq <- max(ex_sq[suff])
  if (max_suff_sq > 1e-10) {
    stop("A structurally sufficient partition has excess squared risk ",
         signif(max_suff_sq, 3), " > 1e-10. tol = ", tol,
         " is inconsistent with this DGP's value spacing.", call. = FALSE)
  }

  vals <- sort(unique(gamma0))
  min_value_gap <- if (length(vals) >= 2L) min(diff(vals)) else Inf

  non <- !suff
  delta_sq <- min(ex_sq[non])
  delta_kl <- min(ex_kl[non])

  list(
    enumeration = enumeration,
    gamma0 = gamma0,
    w = w,
    sufficient = suff,
    excess_sq = ex_sq,
    excess_kl = ex_kl,
    # Delta_j of prop:selection-rate. `delta_sq` is the manuscript's own
    # currency (R^{(j)} is a squared-error risk there); `delta_kl` is the
    # currency estimate_att() actually minimises by default.
    delta_sq = delta_sq,
    delta_kl = delta_kl,
    n_partitions = n_part,
    n_sufficient = sum(suff),
    min_sufficient_leaves = min(enumeration$n_leaves[suff]),
    # The tightest competitor outside S_j -- the partition the margin is
    # measured against, reported so the near-tie ST2 builds is inspectable.
    argmin_key_sq = enumeration$keys[[which(non)[which.min(ex_sq[non])]]],
    argmin_key_kl = enumeration$keys[[which(non)[which.min(ex_kl[non])]]],
    n_distinct_values = length(vals),
    min_value_gap = min_value_gap,
    tol = tol
  )
}

# ---- 6 fitted model -> partition -> recovery indicator ------------------

#' Partition of the covariate space induced by a FITTED optimaltrees model
#'
#' Evaluates the fitted tree on every cell of the grid. Uses
#' \code{optimaltrees:::assign_leaf_ids()}, which is the same internal routine
#' \code{optimaltrees::bisect_lambda_to_budget()} itself uses to verify the
#' \code{m_n} leaf-mass floor (via \code{check_leaf_feasibility()}), so the leaf
#' assignment here is the one the fitting path used. It is unexported, hence the
#' \code{:::}; there is no exported equivalent for \code{OptimalTreesModel}
#' (\code{leaf_assignments()} is for \code{RefinedTreeModel}, a different class,
#' and \code{predict()} returns values, not leaf ids).
#'
#' Relies on the fitted split feature indices addressing the ORIGINAL columns of
#' X, which holds exactly because X is binary: \code{discretize_features()} takes
#' its all-binary fast path and returns X unchanged
#' (\code{discretization_metadata$all_binary == TRUE}). Asserted below rather
#' than assumed.
#'
#' @param model An \code{OptimalTreesModel} (e.g.
#'   \code{estimate_att()$nuisance_fits$e_model}).
#' @param cells The cell grid, with columns in the SAME ORDER as the X the
#'   model was fitted on.
#' @return Integer leaf label per cell, canonicalised.
fitted_cell_labels <- function(model, cells) {
  md <- model@discretization_metadata
  if (!isTRUE(md$all_binary)) {
    stop("Fitted model's discretization_metadata does not report ",
         "all_binary = TRUE, so its split feature indices address ",
         "threshold-indicator features rather than the original columns of X. ",
         "This study is binary-covariate only; a non-binary X would make the ",
         "cell-grid leaf assignment silently wrong.", call. = FALSE)
  }
  if (model@n_trees != 1L) {
    stop("Expected exactly 1 fitted tree, got ", model@n_trees,
         ". Rashomon enumeration is not used on this path.", call. = FALSE)
  }
  lab <- optimaltrees:::assign_leaf_ids(model@trees[[1L]], cells)
  if (anyNA(lab) || any(lab < 1L)) {
    stop("assign_leaf_ids() left ", sum(is.na(lab) | lab < 1L),
         " of ", length(lab), " cells unassigned.", call. = FALSE)
  }
  canonical_labels(lab)
}

#' The currency-correct partition-recovery indicator
#'
#' Decides whether \code{Pi_tauhat gamma_0 = gamma_0} exactly, i.e. whether
#' \code{tauhat in S_j} -- Condition (P1)'s actual claim. See design point (1)
#' at the top of this file for why this is NOT a same-partition or
#' same-leaf-count check.
#'
#' @param model A fitted \code{OptimalTreesModel}.
#' @param class_obj Output of \code{sufficient_class()} for the same nuisance.
#' @return A one-row data.frame:
#'   \item{recovered}{\code{TRUE} iff \code{tauhat in S_j} (the headline metric)}
#'   \item{n_leaves_grid}{leaf count of \code{tauhat} as a partition of the grid}
#'   \item{excess_sq, excess_kl}{\code{tauhat}'s own excess population risk --
#'     \code{0} iff recovered, and the size of the approximation error when not}
#'   \item{in_enumerated_class}{consistency check: is \code{tauhat} in the
#'     enumerated \code{T_Lbar}? Expected \code{TRUE}; a \code{FALSE} means the
#'     fit escaped the class the margin was computed over, which would
#'     invalidate the bound comparison and must be reported, not ignored}
#'   \item{partition_key}{canonical key, so within-\code{S_j} VARIATION across
#'     replications is measurable (lem:uniform's empirical signature, spec §7)}
partition_recovered <- function(model, class_obj) {
  cells <- class_obj$enumeration$cells
  lab <- fitted_cell_labels(model, cells)
  ex <- excess_risk(lab, class_obj$gamma0, class_obj$w)
  key <- paste(lab, collapse = ".")
  data.frame(
    recovered = is_sufficient_partition(lab, class_obj$gamma0,
                                        tol = class_obj$tol),
    n_leaves_grid = length(unique(lab)),
    excess_sq = ex[["sq"]],
    excess_kl = ex[["kl"]],
    in_enumerated_class = key %in% class_obj$enumeration$keys,
    partition_key = key,
    stringsAsFactors = FALSE
  )
}

#' prop:selection-rate's bound at a stated c_1
#'
#' \code{Pr(tauhat_j notin S_j) <= |T_Lbar| exp(-c_1 n Delta_j^2)}. The
#' proposition states \code{c_1 > 0} "depending only on the outcome bound, c,
#' and the positive leaf masses" -- it is NOT identified by the proposition, so
#' this function takes it as an explicit argument and the study reports
#' \code{c_1 = 1} as a STATED CONVENTION. Per spec §4, \code{c_1} is
#' deliberately not fitted post hoc to make the curves agree: what is checkable
#' is the decay SHAPE in \code{n * Delta^2}, which \code{implied_c1()} exposes.
#'
#' @param n Sample size(s).
#' @param delta The margin \code{Delta_j}.
#' @param card_T \code{|T_Lbar|}.
#' @param c1 The bound's constant. Default 1.
#' @return Numeric vector, capped at 1 (it is a probability bound).
selection_bound <- function(n, delta, card_T, c1 = 1) {
  pmin(1, card_T * exp(-c1 * n * delta^2))
}

#' The c_1 an observed failure rate implies, for the shape check
#'
#' Solves \code{rate = |T_Lbar| exp(-c_1 n Delta^2)} for \code{c_1}. If
#' prop:selection-rate's shape is right, this is positive and roughly stable
#' across \code{n} within a DGP; a value that drifts toward 0 as \code{n} grows
#' is evidence the empirical decay is SLOWER than exponential in \code{n}.
#' \code{NA} when the observed rate is 0 (bound vacuously satisfied, no
#' information about the constant).
#'
#' @param rate Empirical \code{Pr(tauhat notin S_j)}.
#' @param n Sample size.
#' @param delta The margin.
#' @param card_T \code{|T_Lbar|}.
implied_c1 <- function(rate, n, delta, card_T) {
  ifelse(rate <= 0 | rate >= card_T, NA_real_,
         log(card_T / rate) / (n * delta^2))
}
