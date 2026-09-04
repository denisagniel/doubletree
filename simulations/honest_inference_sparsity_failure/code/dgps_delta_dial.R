# ============================================================
# dgps_delta_dial.R
# Study: honest_inference_sparsity_failure  (doubletree, S2)
# Spec:  quality_reports/specs/2026-08-21_honest-inference-sparsity-failure.md
#        (see its "DGPs -- a single family with a controlled sparsity-violation
#        dial"; the 2026-09-01 changelog at the top of that file is authoritative)
#
# THE ONE GENUINELY NEW PIECE OF INFRASTRUCTURE THIS STUDY NEEDS. Every existing
# DGP in simulations/dgps/ either satisfies ass:sparsity exactly by construction
# or fails it by an UNMEASURED amount; this study's whole point is a KNOWN,
# CONTROLLABLE magnitude of sparsity violation.
#
# WHAT IS REUSED AND WHAT IS NEW
#
#   REUSED, unmodified, by source()-ing S1's file:
#     ../../partition_recovery_clt/code/enumerate_sufficient_class.R
#   provides make_cell_grid(), enumerate_tree_partitions(), project_on_partition(),
#   excess_risk(), is_sufficient_partition(), fitted_cell_labels(),
#   canonical_labels(). Spec changelog item 5 requires this rather than a second,
#   independent definition of delta_e/delta_mu: "both studies must agree on what
#   these quantities mean for the same underlying DGP family."
#
#   NEW here:
#     dgp_spec_delta_dial()  the eps-dial DGP family itself
#     approximation_error()  delta_j = the EXCESS RISK OF THE BEST <=Lbar-leaf
#                            TREE (spec's own words). NOTE this is a DIFFERENT
#                            quantity from S1's `delta_sq`/`delta_kl`, which are
#                            the selection MARGIN (min excess over partitions
#                            OUTSIDE S_j). S1 needed the margin because its DGPs
#                            satisfy ass:sparsity and the question was whether
#                            selection FINDS a sufficient partition. Here
#                            ass:sparsity FAILS by design, S_j is EMPTY for
#                            eps != 0, and the quantity ass:sparsity is stated in
#                            terms of -- delta_e = delta_mu = 0 -- is the
#                            approximation error. Both are computed and reported;
#                            they coincide (at 0) exactly when eps = 0.
#     population_bias()      the EXACT population bias of the plug-in EIF
#                            estimator at fixed nuisance functions, in closed
#                            form over the 2^p cells. This is what makes
#                            cor:width's plateau height a PREDICTION rather than
#                            a post-hoc description, and what makes the
#                            blind-spot variant's "b = 0" checkable to machine
#                            precision instead of asserted.
#     bilinear_ip()          the {1-e_0}-weighted inner product of prop:bilinear.
#
# ------------------------------------------------------------
# THE CONSTRUCTION, AND WHY EACH CHOICE IS FORCED
# ------------------------------------------------------------
#
# p = 5 binary covariates (M = 32 cells, exact enumeration trivial), X_k iid
# Bernoulli(1/2). Roles:
#   X1  the CONFOUNDING GATE: it drives BOTH e_step and mu_step.
#   X4  the residual coordinate for `g` (propensity) and, in the ALIGNED
#       variants, for `h` (control outcome) as well.
#   X5  the residual coordinate for `h` in the ORTHOGONAL blind-spot variant.
#   X2, X3  pure noise, never entering e_0 or mu_0 at any eps.
#
#   e_0(x)  = e_step(x)  + eps * g(x)
#   mu_0(x) = mu_step(x) + eps * h(x)          mu_1 = mu_0 + tau
#
#   e_step  = e_lo  if X1 = 0, e_hi  if X1 = 1     (2 leaves, depth 1)
#   mu_step = mu_lo if X1 = 0, mu_hi if X1 = 1     (2 leaves, depth 1)
#   g = h = r4 := 2*X4 - 1  in {-1,+1}             (ALIGNED, bias-growing)
#   g = r4, h = r5 := 2*X5 - 1                     (ORTHOGONAL, blind spot)
#
# (1) WHY leaf_budget = 2, AND WHY THE STEP USES EXACTLY 2 LEAVES.
#     delta_j has to be a clean dial, and the pseudo-true partition has to be
#     UNIQUE and eps-invariant, or "the exact bias at the pseudo-true limit" (and
#     hence cor:width's predicted plateau height) is not even well defined. Both
#     properties hold iff the step consumes the ENTIRE budget, leaving no spare
#     leaf that could partially absorb the residual. With a spare leaf the argmin
#     of the population risk over T_Lbar becomes a set of symmetric
#     partial-absorption partitions (8 of them, in the Lbar = 4 version tried
#     first), each pairing differently against mu's own argmin set, and the
#     "predicted plateau" degenerates into a range. At Lbar = 2:
#       T_2 = {trivial} U {single split on X1..X5}, so |T_2| = 6;
#       excess(X1-split) = Var(eps*g) = eps^2   (g has zero mean in each X1 leaf)
#       excess(Xk-split), k != 1 = Var(step) + eps^2 (up to O(eps^3) weighting)
#       excess(trivial)          = Var(step) + eps^2
#     so the argmin is the X1 split, and delta_j = eps^2, PROVIDED the step's own
#     variance dominates the residual's IN THE MEASURE THAT NUISANCE IS FIT UNDER.
#
#     THE TRAP THAT ACTUALLY SPRUNG, recorded because a comment asserting
#     "eps^2 < Var(step)" is exactly what would have hidden it. The two nuisances
#     are fit under DIFFERENT measures: the propensity tree sees P(x) (all rows),
#     the mu tree sees nu_x = P(x){1 - e_0(x)} (controls only). Confounding makes
#     nu strongly asymmetric across the X1 gate -- control mass is depleted
#     exactly where the propensity is high -- so
#         Var_P(e_step)  = ((e_hi - e_lo)/2)^2                       = 0.0506
#         Var_nu(mu_step) = q(1-q)(mu_hi - mu_lo)^2,  q = (1-e_hi)/{(1-e_lo)+(1-e_hi)}
#     and q = 0.286 at the pinned e values, so Var_nu is roughly a THIRD of what
#     a naive p = 1/2 reading would give. A first draft of this file checked only
#     the propensity condition and set (mu_lo, mu_hi, eps_max) = (0.30, 0.60,
#     0.18); the enumeration then reported mu's argmin as the X4 SPLIT, not the
#     X1 split, because eps^2 = 0.0324 > Var_nu(mu_step) = 0.0184. The visible
#     symptom was small and easy to read past -- delta_mu came out at 0.55 * eps^2
#     instead of ~eps^2 -- but the consequence was not: the mu tree's pseudo-true
#     limit stopped tracking mu_step at all, so the realised bias became
#     NON-MONOTONE in eps (0.049 at eps = 0.10 against 0.022 at eps = 0.18),
#     which would have destroyed both cor:width's "width grows with eps then
#     plateaus" reading and prop:spectest's power story, while every individual
#     number still looked plausible.
#
#     The fix is therefore twofold and both halves are load-bearing:
#       (i)  the pinned values satisfy eps_max^2 < Var_nu(mu_step) with margin
#            (0.0225 < 0.0356, factor 1.58) as well as eps_max^2 < Var_P(e_step)
#            (0.0225 < 0.0506, factor 2.25), checked in closed form in
#            dgp_spec_delta_dial(); and
#       (ii) build_dgp_dial() ASSERTS, from the enumeration itself, that BOTH
#            argmin partitions are the X1 split. (i) is a necessary heuristic in
#            the squared-error currency; (ii) is the actual gate, because
#            selection happens in the LOG-LOSS currency, where the two excess
#            risks are reweighted by 1/{m(1-m)} and the ordering is not implied by
#            the variance comparison.
#
# (2) WHY THE RESIDUAL IS ADDITIVE IN ONE COORDINATE AND NOT A PARITY TERM.
#     The spec suggests "e.g. a parity-type interaction term" -- an example, not
#     a requirement; what it REQUIRES is that no <=Lbar-leaf partition absorbs
#     the residual exactly. At Lbar = 2 a single additive +/- eps step on X4
#     already forces e_0 to take 4 distinct values, so 4 leaves are needed and
#     Lbar = 2 cannot represent it: ass:sparsity fails, by a measured amount.
#     Parity(X4, X5) would ALSO work and would make delta exactly eps^2 by the
#     same balance argument -- but it puts the TRUTH at depth 4 with 12-16
#     leaves, and the anchor (estimate_att_crossfit(), max_depth = 4L by default)
#     then has to find a depth-4, 12-leaf tree before Arm 0 can pass. The
#     additive version puts the truth at depth 2 with 4 leaves, i.e. only 2
#     leaves beyond what the constrained estimator may use. That matters
#     specifically because of Arm 0 (spec's new 2026-09-01 arm): the anchor must
#     satisfy ass:rate on THIS DGP, and the anchor is itself tree-based, so the
#     study is only informative about thm:anchor if the anchor can escape the
#     sparsity violation. Its escape hatch is that it carries NO leaf budget --
#     making the truth cheap to reach is therefore a first-class design goal, not
#     a convenience. Whether it actually escapes is measured (Arm 0), not assumed.
#
# (3) WHY X1 IS SHARED (the S1 lesson, applied). S1's own dgps.R header note (c)
#     records a defect that a documented design-audit fix failed to reach the code
#     that ran: e_0 on (X1, X2) and mu_0 on the DISJOINT (X3, X4) makes A
#     independent of Y(0), so NO nuisance misfit can move theta_hat and the
#     study's central metric pairing is vacuous. Here the whole point is that
#     sparsity failure BIASES theta_hat, so confounding is load-bearing twice
#     over. X1 raises both e_0 (0.25 -> 0.70) and mu_0 (0.30 -> 0.60), giving
#     E[mu_0 | A=1] - E[mu_0 | A=0] = +0.135 in closed form.
#     check_variable_roles_dial() makes the shared gate machine-checkable and
#     code/verify_dgps.R additionally computes cor(A, Y0) and the two conditional
#     means EMPIRICALLY from an actual draw -- deliberately not trusting this
#     comment, which is exactly the discipline the S1 defect demands.
#
# (4) WHY THE BLIND SPOT MOVES ONE COORDINATE AND NOTHING ELSE.
#     prop:spectest's stated blind spot is orthogonal misspecification: both
#     nuisances badly wrong, bilinear remainder zero. The ORTHOGONAL variant is
#     the ALIGNED severe variant with h's coordinate changed from X4 to X5 --
#     same eps, same e_0, same |g| and |h|, same delta_e, and (see below) the
#     same delta_mu to machine precision. The ONLY difference is which coordinate
#     h lives on, and the bias goes from ~0.16 to EXACTLY 0. Orthogonality is
#     genuine, not support separation:
#         <g, h>_nu = sum_x P(x){1 - e_0(x)} r4(x) r5(x)
#                   = sum_x P(x){1 - e_step(x) - eps r4(x)} r4(x) r5(x)
#     Every term either carries a bare r5 summed over X5 at fixed everything
#     else, or r4^2 = 1 times a bare r5 -- both vanish, so the inner product is
#     exactly 0. (A disjoint-support construction, g on {X1=1} and h on {X1=0},
#     was the alternative; it also gives b = 0 but HALVES both deltas, weakening
#     the "both nuisances substantially wrong" half of the claim and confounding
#     the contrast with a magnitude change.) bilinear_ip() computes the inner
#     product directly and population_bias() computes b at the pseudo-true
#     limits; both are printed by code/verify_dgps.R.
#
#     Note also that delta_mu = eps^2 EXACTLY in the orthogonal variant but
#     eps^2 * {1 - eps^2/(1-e)^2}-ish in the aligned ones: nu depends on X4
#     through e_0, so h = r4 is very slightly non-centred under nu, while h = r5
#     is exactly centred. The difference is O(eps^4) and is reported, not glossed.
#
# BOUNDED OUTCOME. Y is binary in every variant, so prop:selection-rate's and
# thm:main's boundedness preconditions hold. Asserted in new_dgp_spec_dial().
#
# theta_0 IS EXACT. tau is additive on the probability scale and never clipped,
# so theta_0 = tau; it is nevertheless computed from the closed-form cell
# probabilities (theta_0 = sum_x p_x e_0 (mu_1 - mu_0) / sum_x p_x e_0) and the
# identity is checked, not assumed.
# ============================================================

# ---- 0 the reused enumeration --------------------------------------------

#' Source S1's enumeration helper (spec changelog item 5)
#'
#' Read-only use of another study's code. Sourced through a function so the path
#' is asserted and the failure message names the dependency, rather than an
#' opaque "object 'enumerate_tree_partitions' not found" 200 lines later.
#'
#' @param pkg_root doubletree package root.
#' @return Invisibly, the path sourced.
source_s1_enumeration <- function(pkg_root = getwd()) {
  path <- file.path(pkg_root, "simulations", "partition_recovery_clt", "code",
                    "enumerate_sufficient_class.R")
  if (!file.exists(path)) {
    stop("S1's enumeration helper is missing: ", path, "\n",
         "This study is required by its spec (2026-09-01 changelog item 5) to ",
         "reuse it rather than re-derive delta_e/delta_mu independently.",
         call. = FALSE)
  }
  source(path, local = FALSE)
  needed <- c("make_cell_grid", "enumerate_tree_partitions", "project_on_partition",
              "excess_risk", "is_sufficient_partition", "fitted_cell_labels",
              "canonical_labels")
  absent <- needed[!vapply(needed, function(f) is.function(get0(f)), logical(1))]
  if (length(absent)) {
    stop("S1's enumeration helper defined no ", paste(absent, collapse = ", "),
         "; its contract has changed.", call. = FALSE)
  }
  invisible(path)
}

# ---- 1 the eps dial ------------------------------------------------------

# Pinned population values. Chosen jointly under five constraints, all of which
# code/verify_dgps.R re-checks numerically:
#   (a) e_0 and mu_0, mu_1 stay strictly inside (0, 1) at the LARGEST eps, and
#       e_0 stays clear of doubletree's (0.01, 0.99) clip so ass:construct(c)'s
#       clipping never bites and cannot be confused with a sparsity effect;
#   (b) Var_P(e_step) = ((e_hi - e_lo)/2)^2 = 0.0506 > eps^2 at every eps used;
#   (c) Var_nu(mu_step) = q(1-q)(mu_hi - mu_lo)^2 = 0.0356 > eps^2 at every eps
#       used, with q = (1 - e_hi)/{(1 - e_lo) + (1 - e_hi)} = 0.286 -- the
#       CONTROL-weighted analogue of (b), and the one a first draft omitted (see
#       design note (1)); together (b) and (c) keep the pseudo-true partition of
#       BOTH nuisances at the X1 split;
#   (d) the confounding gap E[mu_0|A=1] - E[mu_0|A=0] is substantial;
#   (e) the implied bias b(eps) ~ 4.9 eps^2 spans BELOW to WELL ABOVE the Wald
#       scale at n in {500, 2000, 8000}, so cor:width's plateau TRANSITION is
#       inside the sweep and not merely its plateaued end (spec's own requirement
#       for the eps grid).
E_LO_DIAL  <- 0.25
E_HI_DIAL  <- 0.70
MU_LO_DIAL <- 0.25
MU_HI_DIAL <- 0.65
TAU_DIAL   <- 0.10
P_DIAL     <- 5L

#' Population spec for one point of the eps dial
#'
#' @param eps Sparsity-violation dial. \code{0} is the exactly-sparse favourable
#'   case (ass:sparsity holds, delta_e = delta_mu = 0).
#' @param residual_mode \code{"aligned"} (\code{g = h = r4}; the bilinear
#'   remainder is maximal, so the bias grows like \code{eps^2}) or
#'   \code{"orthogonal"} (\code{g = r4}, \code{h = r5}; \code{<g,h>_nu = 0}
#'   exactly, so the true bias is exactly 0 while both nuisances are equally
#'   misspecified -- prop:spectest's stated blind spot).
#' @param leaf_budget Lbar. See design note (1) for why 2 is not a free choice.
#' @param e_lo,e_hi,mu_lo,mu_hi,tau Population leaf values / treatment effect.
#' @param p Number of binary covariates (>= 5; X2, X3 and any X6.. are noise).
#' @param id,label Identifiers carried into the results.
#' @return A DGP spec list; pass to \code{build_dgp_dial()}.
dgp_spec_delta_dial <- function(eps,
                                residual_mode = c("aligned", "orthogonal"),
                                leaf_budget = 2L,
                                e_lo = E_LO_DIAL, e_hi = E_HI_DIAL,
                                mu_lo = MU_LO_DIAL, mu_hi = MU_HI_DIAL,
                                tau = TAU_DIAL, p = P_DIAL,
                                id = NULL, label = NULL) {
  residual_mode <- match.arg(residual_mode)
  stopifnot(length(eps) == 1L, is.finite(eps), eps >= 0)
  stopifnot(length(p) == 1L, p >= 5L, p == as.integer(p))
  p <- as.integer(p)
  leaf_budget <- as.integer(leaf_budget)

  if (is.null(id)) id <- sprintf("%s_eps%s", substr(residual_mode, 1, 4),
                                 format(eps, trim = TRUE))
  if (is.null(label)) label <- id

  cells <- make_cell_grid(p)
  # X_k iid Bernoulli(1/2) => every cell has probability 2^-p.
  cell_prob <- rep(1 / nrow(cells), nrow(cells))

  r4 <- 2L * cells$X4 - 1L
  r5 <- 2L * cells$X5 - 1L
  g <- as.numeric(r4)
  h <- if (residual_mode == "aligned") as.numeric(r4) else as.numeric(r5)

  e_step <- ifelse(cells$X1 == 0L, e_lo, e_hi)
  # SAME X1 gate as e_step: this is the confounding (design note (3)).
  mu_step <- ifelse(cells$X1 == 0L, mu_lo, mu_hi)

  e0 <- e_step + eps * g
  mu0 <- mu_step + eps * h
  mu1 <- mu0 + tau

  # Fail at the construction, not at replication 200. A clipped e_0 or a mu_1
  # pushed past 1 would silently make tau(x) non-constant and theta_0 != tau.
  if (any(e0 <= 0.02) || any(e0 >= 0.98)) {
    stop("eps = ", eps, " drives e_0 to [", signif(min(e0), 4), ", ",
         signif(max(e0), 4), "], too close to doubletree's (0.01, 0.99) clip. ",
         "Clipping would confound ass:construct(c) with the sparsity effect.",
         call. = FALSE)
  }
  if (any(mu0 <= 0) || any(mu0 >= 1) || any(mu1 <= 0) || any(mu1 >= 1)) {
    stop("eps = ", eps, " drives mu_0 or mu_1 outside (0, 1): mu_0 in [",
         signif(min(mu0), 4), ", ", signif(max(mu0), 4), "], mu_1 in [",
         signif(min(mu1), 4), ", ", signif(max(mu1), 4), "].", call. = FALSE)
  }
  # Design note (1): uniqueness of the pseudo-true partition needs the step's own
  # variance to dominate the residual's IN EACH NUISANCE'S OWN MEASURE. Both are
  # checked in closed form here; build_dgp_dial() then asserts the actual argmin,
  # which is the statement that matters, because selection happens in the
  # log-loss currency where these variance comparisons are only a heuristic.
  var_p_e_step <- ((e_hi - e_lo) / 2)^2
  q_nu <- (1 - e_hi) / ((1 - e_lo) + (1 - e_hi))
  var_nu_mu_step <- q_nu * (1 - q_nu) * (mu_hi - mu_lo)^2
  if (eps > 0 && eps^2 >= var_p_e_step) {
    stop("eps^2 = ", signif(eps^2, 4), " >= Var_P(e_step) = ",
         signif(var_p_e_step, 4), ": the propensity's population-risk argmin ",
         "over T_Lbar would stop being the X1 split, so delta_e would no longer ",
         "be eps^2 and the pseudo-true partition would move with eps.",
         call. = FALSE)
  }
  if (eps > 0 && eps^2 >= var_nu_mu_step) {
    stop("eps^2 = ", signif(eps^2, 4), " >= Var_nu(mu_step) = ",
         signif(var_nu_mu_step, 4), " (nu-weight share of {X1=1} is q = ",
         signif(q_nu, 4), "): the CONTROL-WEIGHTED measure the mu tree is fit ",
         "under makes mu_step's variance much smaller than the propensity's, so ",
         "mu's argmin would flip to the residual coordinate and the realised ",
         "bias would stop being monotone in eps. See design note (1): this is ",
         "the exact failure a propensity-only check let through.", call. = FALSE)
  }

  new_dgp_spec_dial(
    id = id, label = label, p = p, leaf_budget = leaf_budget,
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, mu1 = mu1, tau = tau,
    eps = eps, residual_mode = residual_mode, g = g, h = h,
    e_step = e_step, mu_step = mu_step,
    draw = function(n) draw_delta_dial(n, p, e_lo, e_hi, mu_lo, mu_hi, tau,
                                       eps, residual_mode),
    # Declared coordinate roles, checked on every build_dgp_dial() call. X1 in
    # BOTH lists is the machine-readable statement that this DGP is confounded
    # (S1's dgps.R made exactly this assertion for exactly this reason).
    e_vars = if (eps > 0) c("X1", "X4") else "X1",
    mu_vars = if (eps > 0) {
      c("X1", if (residual_mode == "aligned") "X4" else "X5")
    } else {
      "X1"
    },
    params = list(eps = eps, residual_mode = residual_mode,
                  e_lo = e_lo, e_hi = e_hi, mu_lo = mu_lo, mu_hi = mu_hi,
                  tau = tau, p = p,
                  var_p_e_step = var_p_e_step, var_nu_mu_step = var_nu_mu_step,
                  q_nu = q_nu)
  )
}

#' One sample from the eps-dial DGP
#'
#' Draw order (X columns, then A, then Y0, Y1) is fixed, so a given seed
#' reproduces a given dataset byte for byte. Returns the DETERMINISTIC nuisance
#' vectors as well, so check_spec_vs_draw_dial() can compare population against
#' sample at tolerance 1e-12 instead of at sampling noise.
draw_delta_dial <- function(n, p, e_lo, e_hi, mu_lo, mu_hi, tau, eps,
                            residual_mode) {
  X <- as.data.frame(
    lapply(seq_len(p), function(k) as.integer(stats::rbinom(n, 1L, 0.5)))
  )
  names(X) <- paste0("X", seq_len(p))
  g <- 2L * X$X4 - 1L
  h <- if (residual_mode == "aligned") 2L * X$X4 - 1L else 2L * X$X5 - 1L
  e <- ifelse(X$X1 == 0L, e_lo, e_hi) + eps * g
  A <- as.integer(stats::rbinom(n, 1L, e))
  m0 <- ifelse(X$X1 == 0L, mu_lo, mu_hi) + eps * h
  m1 <- m0 + tau
  Y0 <- as.integer(stats::rbinom(n, 1L, m0))
  Y1 <- as.integer(stats::rbinom(n, 1L, m1))
  list(X = X, A = A, Y = A * Y1 + (1L - A) * Y0,
       e_true = e, mu0_true = m0,
       # Y0 is returned for the CONFOUNDING check: cor(A, Y0) and
       # E[Y0|A=1] - E[Y0|A=0] are computable only with the counterfactual
       # control outcome for every unit, treated ones included.
       Y0 = Y0, Y1 = Y1)
}

# ---- 2 the shared constructor: exact population quantities ---------------

#' Assemble a spec and its exact population quantities
#'
#' Not called directly by study code; \code{dgp_spec_delta_dial()} calls it.
new_dgp_spec_dial <- function(id, label, p, leaf_budget, cells, cell_prob,
                              e0, mu0, mu1, tau, eps, residual_mode, g, h,
                              e_step, mu_step, draw, e_vars, mu_vars, params) {
  n_cells <- nrow(cells)
  stopifnot(length(e0) == n_cells, length(mu0) == n_cells, length(mu1) == n_cells,
            length(cell_prob) == n_cells, length(g) == n_cells, length(h) == n_cells)
  stopifnot(all(e_vars %in% names(cells)), all(mu_vars %in% names(cells)))
  if (any(e0 <= 0 | e0 >= 1)) {
    stop("e_0 must be strictly inside (0, 1) for positivity (ass:causal).",
         call. = FALSE)
  }
  if (any(mu0 < 0 | mu0 > 1) || any(mu1 < 0 | mu1 > 1)) {
    stop("mu_0 and mu_1 must lie in [0, 1]: Y is binary, hence BOUNDED, which ",
         "is thm:main's and prop:selection-rate's own precondition.",
         call. = FALSE)
  }
  if (abs(sum(cell_prob) - 1) > 1e-12) stop("cell_prob must sum to 1.", call. = FALSE)

  p_treated <- sum(cell_prob * e0)
  theta0 <- sum(cell_prob * e0 * (mu1 - mu0)) / p_treated
  # tau is additive and unclipped, so theta_0 = tau identically. Checked rather
  # than assumed: a clipped mu_1 would break it silently.
  if (abs(theta0 - tau) > 1e-12) {
    stop("theta_0 = ", signif(theta0, 10), " != tau = ", tau,
         "; tau(x) is not constant, so mu_1 must have been clipped.",
         call. = FALSE)
  }

  nu <- cell_prob * (1 - e0)

  list(
    id = id, label = label, p = p, leaf_budget = as.integer(leaf_budget),
    var_names = names(cells),
    cells = cells, cell_prob = cell_prob,
    e0 = e0, mu0 = mu0, mu1 = mu1,
    e_step = e_step, mu_step = mu_step,
    eps = eps, residual_mode = residual_mode, g = g, h = h,
    tau_x = mu1 - mu0, tau = tau,
    theta0 = theta0, p_treated = p_treated,
    # Control-weighted measure nu_x = p_x {1 - e_0(x)}: the weighting under which
    # Pi^nu_tau and R^{(mu)} are defined, and the one the mu tree actually sees,
    # since estimate_att() fits it on X[A == 0, ].
    nu = nu,
    e0_range = range(e0), mu0_range = range(mu0),
    draw = draw, e_vars = e_vars, mu_vars = mu_vars, params = params
  )
}

# ---- 3 delta_j: the approximation error of the BEST <=Lbar-leaf tree -----

#' Exact approximation error, pseudo-true partition, and (when non-empty) S_j
#'
#' \code{delta_j} of ass:sparsity is the excess population risk of the BEST
#' partition in \code{T_Lbar} -- \code{0} iff some <=Lbar-leaf tree represents
#' \code{gamma_0} exactly. Computed here over the whole enumerated class, in both
#' currencies, from closed-form population probabilities: no data, no fitting.
#'
#' S1's selection MARGIN is also returned (\code{margin_sq}/\code{margin_kl}, the
#' min excess over partitions OUTSIDE S_j) whenever \code{S_j} is non-empty, so
#' the two studies' quantities can be read side by side on the shared
#' \code{eps = 0} cell. For \code{eps != 0} \code{S_j} is EMPTY by construction
#' and the margin is \code{NA} -- which is the point of the dial, not a defect.
#'
#' @param gamma0 True nuisance value at each cell, in canonical cell order.
#' @param w Population cell weights: \code{P(X = x)} for the propensity,
#'   \code{nu_x = P(X = x){1 - e_0(x)}} for the control outcome.
#' @param enumeration \code{enumerate_tree_partitions()} output.
#' @param loss Currency the pseudo-true partition is selected in. \code{"kl"}
#'   (default) is what \code{estimate_att()} actually minimises for BOTH
#'   nuisances at \code{propensity_loss = "log_loss"} and
#'   \code{outcome_type = "binary"}; \code{"sq"} is manuscript.tex's own currency
#'   for \code{R^{(j)}}. Both excess-risk vectors are returned either way.
#' @param tol Structural tolerance forwarded to \code{is_sufficient_partition()}.
#' @return A list; see the field comments below.
approximation_error <- function(gamma0, w, enumeration, loss = c("kl", "sq"),
                                tol = 1e-12) {
  loss <- match.arg(loss)
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
  if (anyNA(ex_kl)) {
    stop("Excess KL risk is NA for ", sum(is.na(ex_kl)), " partition(s): ",
         "gamma_0 or one of its projections left (0, 1), so the log-loss ",
         "currency is undefined. This study's DGPs are strictly interior.",
         call. = FALSE)
  }

  ex <- if (loss == "kl") ex_kl else ex_sq
  delta <- min(ex)
  # Every argmin, not just the first: if there is more than one, the pseudo-true
  # limit is NOT unique and "the exact bias at the pseudo-true partition" is a
  # set rather than a number. The caller must see that rather than silently get
  # whichever partition came first out of the enumeration.
  argmin_idx <- which(ex <= delta + 1e-14)

  # Consistency between the structural and the risk criterion. A sufficient
  # partition has exactly zero excess risk mathematically; a mismatch here means
  # `tol` is wrong for this DGP's value spacing.
  if (any(suff) && max(ex_sq[suff]) > 1e-10) {
    stop("A structurally sufficient partition has excess squared risk ",
         signif(max(ex_sq[suff]), 3), " > 1e-10 at tol = ", tol, ".",
         call. = FALSE)
  }
  if (any(suff) && delta > 1e-10) {
    stop("S_j is non-empty but the minimum excess risk is ", signif(delta, 3),
         " > 0, which is contradictory.", call. = FALSE)
  }
  if (!any(suff) && delta <= 1e-10) {
    stop("No partition is structurally sufficient yet the minimum excess risk ",
         "is ", signif(delta, 3), " ~ 0. tol = ", tol, " is inconsistent with ",
         "this DGP's value spacing.", call. = FALSE)
  }

  non <- !suff
  vals <- sort(unique(gamma0))

  list(
    enumeration = enumeration,
    gamma0 = gamma0, w = w, loss = loss,
    sufficient = suff, excess_sq = ex_sq, excess_kl = ex_kl,
    # delta_j of ass:sparsity, in BOTH currencies. 0 iff S_j is non-empty.
    delta_sq = min(ex_sq), delta_kl = min(ex_kl),
    # S1's selection margin, for the eps = 0 cell the two studies share.
    margin_sq = if (any(non) && any(suff)) min(ex_sq[non]) else NA_real_,
    margin_kl = if (any(non) && any(suff)) min(ex_kl[non]) else NA_real_,
    n_partitions = n_part,
    n_sufficient = sum(suff),
    sparsity_holds = any(suff),
    # The pseudo-true partition(s): the population-risk argmin in the currency
    # the estimator actually minimises.
    argmin_idx = argmin_idx,
    n_argmin = length(argmin_idx),
    argmin_key = enumeration$keys[[argmin_idx[[1]]]],
    argmin_labels = enumeration$labels[argmin_idx[[1]], ],
    argmin_leaves = enumeration$n_leaves[[argmin_idx[[1]]]],
    # Pi^w gamma_0 at the pseudo-true partition: the nuisance function the
    # constrained estimator converges to, and hence the input to population_bias().
    pseudo_true = project_on_partition(enumeration$labels[argmin_idx[[1]], ],
                                       gamma0, w),
    n_distinct_values = length(vals),
    min_value_gap = if (length(vals) >= 2L) min(diff(vals)) else Inf,
    tol = tol
  )
}

# ---- 4 the exact population bias, in closed form -------------------------

#' Exact population bias of the plug-in EIF estimator at FIXED nuisances
#'
#' \code{estimate_att()} solves \eqn{sum_i psi_i(theta) = 0} with
#' \code{\link[doubletree]{psi_att}}'s score
#'   psi = (A/pi)(Y - m0 - theta) - (1/pi) {e/(1-e)} (1-A)(Y - m0).
#' Taking expectations at fixed \code{(e_hat, m0_hat)} and solving
#' \code{E[psi(theta)] = 0} gives the estimator's population limit
#'   theta_inf = [ E{e_0 (mu_1 - m0_hat)} - E{(1-e_0) odds(e_hat) (mu_0 - m0_hat)} ] / pi,
#' with \code{pi = E[e_0]}. Every expectation is a finite sum over the 2^p cells,
#' so this is EXACT, not simulated.
#'
#' Equivalently (and returned as \code{bilinear_form}, as an algebraic check on
#' the above), \code{pi * (theta_inf - theta_0)} equals the bilinear remainder
#'   E[ (1 - e_0) {odds(e_0) - odds(e_hat)} (mu_0 - m0_hat) ],
#' which is prop:bilinear's own object: it vanishes when the two error functions
#' are orthogonal in the \code{{1-e_0}}-weighted inner product, HOWEVER LARGE
#' either is individually. That identity is what makes the blind-spot variant's
#' \code{b = 0} a computable fact rather than an assertion.
#'
#' @param spec A \code{dgp_spec_delta_dial()} result (built or unbuilt).
#' @param e_hat Propensity function per cell (the fitted/pseudo-true one).
#' @param m0_hat Control-outcome function per cell.
#' @param clip Propensity clip applied inside \code{estimate_att()}. Applied here
#'   too, so the population limit is the limit of what the package computes.
#' @return A list with \code{theta_inf}, \code{bias}, \code{bilinear_form},
#'   \code{bias_from_bilinear} (must equal \code{bias}), and the two error norms
#'   \code{D_w} (odds-scale propensity error) and \code{D_mu}, which are
#'   cor:width's own \eqn{D_w, D_mu}.
population_bias <- function(spec, e_hat, m0_hat, clip = c(0.01, 0.99)) {
  P <- spec$cell_prob
  e0 <- spec$e0
  mu0 <- spec$mu0
  mu1 <- spec$mu1
  stopifnot(length(e_hat) == length(P), length(m0_hat) == length(P))
  e_hat <- pmax(clip[[1]], pmin(clip[[2]], e_hat))

  pi_pop <- sum(P * e0)
  odds_hat <- e_hat / (1 - e_hat)
  odds_0 <- e0 / (1 - e0)

  theta_inf <- (sum(P * e0 * (mu1 - m0_hat)) -
                  sum(P * (1 - e0) * odds_hat * (mu0 - m0_hat))) / pi_pop
  bias <- theta_inf - spec$theta0

  # prop:bilinear's remainder, weighted by nu_x = P_x (1 - e_0(x)).
  bilinear_form <- sum(P * (1 - e0) * (odds_0 - odds_hat) * (mu0 - m0_hat))

  list(
    theta_inf = theta_inf,
    bias = bias,
    bilinear_form = bilinear_form,
    bias_from_bilinear = bilinear_form / pi_pop,
    # cor:width's realised-error norms: L2(nu) norms of the two error functions
    # actually entering the bilinear form.
    D_w = sqrt(sum(spec$nu * (odds_hat - odds_0)^2) / sum(spec$nu)),
    D_mu = sqrt(sum(spec$nu * (m0_hat - mu0)^2) / sum(spec$nu)),
    pi_pop = pi_pop,
    n_clipped = sum(e_hat <= clip[[1]] + 1e-12 | e_hat >= clip[[2]] - 1e-12)
  )
}

#' The \code{{1-e_0}}-weighted inner product of prop:bilinear
#'
#' \code{<u, v>_nu = sum_x P(x) {1 - e_0(x)} u(x) v(x)}, normalised by
#' \code{sum_x nu_x} so it is comparable across DGPs. Also returns the
#' correlation, which is the scale-free statement "these two functions are
#' orthogonal" that the blind-spot construction actually claims.
bilinear_ip <- function(spec, u, v) {
  nu <- spec$nu
  s <- sum(nu)
  ip <- sum(nu * u * v) / s
  nu_u <- sqrt(sum(nu * u^2) / s)
  nu_v <- sqrt(sum(nu * v^2) / s)
  # A correlation between two functions that are identically zero is undefined,
  # not "whatever the ratio of two floating-point dusts happens to be". At
  # eps = 0 both pseudo-true residuals ARE identically zero, and an unguarded
  # ratio there returned a confident-looking -0.78 in the first run of
  # code/verify_dgps.R. The threshold is absolute because these are probabilities
  # and odds on an O(1) scale.
  degenerate <- nu_u < 1e-12 || nu_v < 1e-12
  list(ip = ip, norm_u = nu_u, norm_v = nu_v,
       cor = if (degenerate) NA_real_ else ip / (nu_u * nu_v))
}

# ---- 5 guards: the spec must be what it claims to be ---------------------

#' Assert each nuisance depends only on the coordinates it declares
#'
#' An INTENT-level check. A nuisance keyed on the wrong covariate is still a
#' perfectly self-consistent DGP, errors nowhere, and shows up only as a delta of
#' quietly the wrong size -- or, in this study, as a bilinear form that is not
#' the one the blind spot was supposed to exercise. X1 appearing in BOTH declared
#' lists is the machine-readable statement that the DGP is confounded.
check_variable_roles_dial <- function(spec) {
  cells <- spec$cells
  range_flipping <- function(gamma0, v) {
    others <- setdiff(names(cells), v)
    key <- do.call(paste, c(lapply(others, function(o) cells[[o]]), sep = "-"))
    tapply(gamma0, key, function(z) diff(range(z)))
  }
  check_one <- function(gamma0, claimed, nm) {
    for (v in setdiff(names(cells), claimed)) {
      if (any(range_flipping(gamma0, v) > 1e-12)) {
        stop(nm, " depends on ", v, ", not in its declared variables (",
             paste(claimed, collapse = ", "), "). Either the formula or the ",
             "declared role is wrong; both change delta_j and the bilinear form.",
             call. = FALSE)
      }
    }
    for (v in claimed) {
      if (all(range_flipping(gamma0, v) <= 1e-12)) {
        stop(nm, " does NOT depend on its declared variable ", v,
             "; the declared structure is stale.", call. = FALSE)
      }
    }
    invisible(TRUE)
  }
  check_one(spec$e0, spec$e_vars, "e_0")
  check_one(spec$mu0, spec$mu_vars, "mu_0")
  if (!("X1" %in% spec$e_vars && "X1" %in% spec$mu_vars)) {
    stop("X1 is not declared for both nuisances: the DGP would be ",
         "UNCONFOUNDED, and no amount of sparsity failure could then bias ",
         "theta_hat. This is exactly the defect S1's design audit caught.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Assert the population spec agrees with the sample draw
#'
#' Compares the generator's own DETERMINISTIC \code{e_true}/\code{mu0_true} to
#' the spec's closed-form values, per cell, at tolerance 1e-12 -- not noisy
#' empirical rates, so a transcription error cannot hide behind sampling noise.
check_spec_vs_draw_dial <- function(spec, n = 40000L, tol = 1e-12,
                                    seed = 20260821L) {
  set.seed(seed)
  d <- spec$draw(n)
  for (fld in c("e_true", "mu0_true", "Y0")) {
    if (is.null(d[[fld]])) {
      stop("DGP '", spec$id, "' draw() returns no ", fld, ".", call. = FALSE)
    }
  }
  key_of <- function(df) {
    do.call(paste, c(lapply(names(spec$cells), function(v) df[[v]]), sep = "-"))
  }
  samp_key <- key_of(d$X)
  cell_key <- key_of(spec$cells)
  if (!setequal(unique(samp_key), cell_key)) {
    stop("Sampled covariate patterns do not exhaust the ", nrow(spec$cells),
         " cells at n = ", n, " for DGP '", spec$id, "'.", call. = FALSE)
  }
  idx <- match(cell_key, samp_key)
  cmp <- data.frame(
    dgp = spec$id, cell = cell_key,
    e0_spec = spec$e0, e0_drawn = d$e_true[idx],
    mu0_spec = spec$mu0, mu0_drawn = d$mu0_true[idx],
    stringsAsFactors = FALSE
  )
  bad <- abs(cmp$e0_spec - cmp$e0_drawn) > tol |
    abs(cmp$mu0_spec - cmp$mu0_drawn) > tol
  if (any(bad)) {
    print(cmp[bad, ])
    stop("Population spec disagrees with the sample draw for '", spec$id,
         "' on ", sum(bad), " cell(s).", call. = FALSE)
  }
  invisible(cmp)
}

#' Empirical confounding diagnostics from an actual draw
#'
#' S1's lesson, applied literally: a documented confounding fix that never
#' reached the running code made that study's central metric pairing vacuous, and
#' it was caught only by computing \code{cor(A, Y0)} and
#' \code{E[Y0|A=1] - E[Y0|A=0]} from data rather than reading the comments. So
#' these are computed here, from a draw, for every DGP, every pilot run.
#'
#' \code{Y0} is the counterfactual control outcome for EVERY unit (treated ones
#' included), which is why \code{draw_delta_dial()} returns it.
confounding_check <- function(spec, n = 200000L, seed = 20260821L) {
  set.seed(seed)
  d <- spec$draw(n)
  A <- d$A
  Y0 <- d$Y0
  # Closed-form counterparts, so "the draw and the population agree" is visible
  # rather than inferred.
  P <- spec$cell_prob
  e0 <- spec$e0
  mu0 <- spec$mu0
  pop_treated <- sum(P * e0 * mu0) / sum(P * e0)
  pop_control <- sum(P * (1 - e0) * mu0) / sum(P * (1 - e0))
  data.frame(
    dgp = spec$id, eps = spec$eps, residual_mode = spec$residual_mode, n = n,
    emp_cor_A_Y0 = stats::cor(A, Y0),
    emp_EY0_treated = mean(Y0[A == 1L]),
    emp_EY0_control = mean(Y0[A == 0L]),
    emp_gap = mean(Y0[A == 1L]) - mean(Y0[A == 0L]),
    pop_Emu0_treated = pop_treated,
    pop_Emu0_control = pop_control,
    pop_gap = pop_treated - pop_control,
    stringsAsFactors = FALSE
  )
}

# ---- 6 assemble: enumeration + delta_j + pseudo-true bias ----------------

#' Attach the exact approximation errors, pseudo-true limits and bias to a spec
#'
#' Runs the (deterministic, cheap at \code{p = 5}, \code{Lbar = 2}) enumeration
#' ONCE per DGP. Every replication then reuses the result.
#'
#' \code{check_variable_roles_dial()} runs first, unconditionally: it is cheap,
#' and a mis-keyed nuisance has no other symptom than quantities of the wrong
#' size.
#'
#' @param spec A \code{dgp_spec_delta_dial()} result.
#' @return \code{spec} with \code{enumeration}, \code{class_e}, \code{class_mu},
#'   \code{pseudo} (the exact bias at the pseudo-true limits) and \code{ip_gh}
#'   (the \code{{1-e_0}}-weighted inner product of \code{g} and \code{h}) added.
build_dgp_dial <- function(spec) {
  check_variable_roles_dial(spec)
  en <- enumerate_tree_partitions(spec$p, spec$leaf_budget, spec$var_names)
  spec$enumeration <- en
  # Weights: P(x) for the propensity tree (fit on all rows); nu_x for the mu tree
  # (fit on controls only).
  spec$class_e <- approximation_error(spec$e0, spec$cell_prob, en, loss = "kl")
  spec$class_mu <- approximation_error(spec$mu0, spec$nu, en, loss = "kl")

  # THE HARD GATE (design note (1)). Both pseudo-true partitions must be the X1
  # split -- the partition e_step and mu_step are exactly representable by. If
  # either argmin is anything else, the pseudo-true limit has stopped tracking the
  # step function, delta_j has stopped being ~eps^2, and the realised bias has
  # stopped being monotone in eps: all three silently, with every printed number
  # still looking plausible. This assertion, not the closed-form variance
  # comparisons in dgp_spec_delta_dial(), is what actually catches that, because
  # selection happens in the log-loss currency where the variance ordering does
  # not imply the risk ordering.
  x1_labels <- canonical_labels(spec$cells$X1)
  x1_key <- paste(x1_labels, collapse = ".")
  for (j in c("e", "mu")) {
    cls <- if (j == "e") spec$class_e else spec$class_mu
    if (!identical(cls$argmin_key, x1_key)) {
      stop("DGP '", spec$id, "' (eps = ", spec$eps, "): the pseudo-true ",
           "partition for ", j, " is NOT the X1 split.\n",
           "  argmin leaf count: ", cls$argmin_leaves,
           ", delta_sq = ", signif(cls$delta_sq, 6),
           " (eps^2 = ", signif(spec$eps^2, 6), ", ratio ",
           signif(cls$delta_sq / max(spec$eps^2, .Machine$double.eps), 4), ")\n",
           "  The residual has out-competed the step function in ", j, "'s own ",
           "measure, so the eps dial is no longer a dial. Reduce eps or widen ",
           "the step gap (for mu, remember the measure is nu_x = P_x(1-e_0(x)), ",
           "not P_x -- see design note (1)).", call. = FALSE)
    }
    if (cls$n_argmin != 1L) {
      stop("DGP '", spec$id, "': ", j, " has ", cls$n_argmin, " population-risk ",
           "argmins, so the pseudo-true limit is a SET and 'the exact bias at ",
           "the pseudo-true partition' is not a number.", call. = FALSE)
    }
  }

  # The exact bias the CONSTRAINED estimator converges to. This is cor:width's
  # predicted plateau half-height and, for the blind-spot variant, the number
  # that must be 0 to machine precision.
  spec$pseudo <- population_bias(spec, spec$class_e$pseudo_true,
                                 spec$class_mu$pseudo_true)
  # prop:bilinear's inner product on the RESIDUAL functions g and h themselves
  # (the object the spec's blind-spot instruction names), alongside the same
  # inner product on the realised pseudo-true errors (the object the bias
  # actually depends on). Both are reported; agreement between them is the
  # substantive content of the construction.
  spec$ip_gh <- bilinear_ip(spec, spec$g, spec$h)
  spec$ip_resid <- bilinear_ip(
    spec,
    spec$class_e$pseudo_true / (1 - spec$class_e$pseudo_true) -
      spec$e0 / (1 - spec$e0),
    spec$class_mu$pseudo_true - spec$mu0
  )
  spec
}

#' One-line-per-nuisance summary of a built DGP's population structure
#'
#' Everything a table needs to read the results in the theorems' own quantities,
#' with no simulation involved.
dgp_population_summary_dial <- function(spec) {
  data.frame(
    dgp = spec$id, label = spec$label,
    eps = spec$eps, residual_mode = spec$residual_mode,
    nuisance = c("e", "mu"),
    theta0 = spec$theta0, p = spec$p, leaf_budget = spec$leaf_budget,
    card_T = spec$enumeration$n_partitions,
    # delta_j of ass:sparsity: 0 iff a <=Lbar-leaf tree is exact.
    delta_sq = c(spec$class_e$delta_sq, spec$class_mu$delta_sq),
    delta_kl = c(spec$class_e$delta_kl, spec$class_mu$delta_kl),
    sparsity_holds = c(spec$class_e$sparsity_holds, spec$class_mu$sparsity_holds),
    n_sufficient = c(spec$class_e$n_sufficient, spec$class_mu$n_sufficient),
    # S1's selection margin, defined only where S_j is non-empty.
    margin_kl = c(spec$class_e$margin_kl, spec$class_mu$margin_kl),
    n_argmin = c(spec$class_e$n_argmin, spec$class_mu$n_argmin),
    argmin_key = c(spec$class_e$argmin_key, spec$class_mu$argmin_key),
    argmin_leaves = c(spec$class_e$argmin_leaves, spec$class_mu$argmin_leaves),
    n_distinct_values = c(spec$class_e$n_distinct_values,
                          spec$class_mu$n_distinct_values),
    # cor:width's realised-error norms and the exact bias at the pseudo-true
    # limits. Constant across the two rows by construction (they describe the
    # PAIR), repeated so any single row is self-contained.
    D_w = spec$pseudo$D_w, D_mu = spec$pseudo$D_mu,
    D_product = spec$pseudo$D_w * spec$pseudo$D_mu,
    bias_pseudo = spec$pseudo$bias,
    bilinear_form = spec$pseudo$bilinear_form,
    ip_gh = spec$ip_gh$ip, cor_gh = spec$ip_gh$cor,
    ip_resid = spec$ip_resid$ip, cor_resid = spec$ip_resid$cor,
    e0_min = spec$e0_range[[1]], e0_max = spec$e0_range[[2]],
    stringsAsFactors = FALSE
  )
}
