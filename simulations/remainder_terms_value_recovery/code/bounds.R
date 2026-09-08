# ---------------------------------------------------------------------------
# The crude theoretical bound that the proof audit found (possibly) vacuous.
#
# Localisation radius (VC-entropy form used in the audit):
#     kappa_n = v_{Lbar,p} * log(A_ENT * e * n / v_{Lbar,p}) + log(4 / delta_n)
#     rho_n   = sqrt(kappa_n / n)
#
# with v_{Lbar,p} the VC-subgraph dimension of the leaf-constant tree class. The
# corrected VC result is V_sub = O(Lbar log(Lbar p)); the audit's reference
# anchors are v_{10,5} ~ 185, v_{20,5} ~ 427, v_{30,5} ~ 686. We recover a
# reproducible v(.) by least-squares calibration of the O(Lbar log(Lbar p)) form
# to those three anchors, which also supplies the Lbar = 5 value the audit did
# not tabulate. The literature form (Bartlett, Jordan & McAuliffe 2006, Lemma 9,
# quoted in doubletree/inst/paper/vc-dimension-hybrid-version.tex line 138),
#     V_sub(T_s) <= 4 * d_tilde * s * log2(3 s),
# is reported alongside as a much larger reference point (d_tilde = number of
# BINARY features actually searched).
#
# The lemma bounds being tested (refining-grid-standalone.tex, Step 2):
#     |T_eps| = O_p{ (D_w + rho_n) * rho_n }
#     |T_a|   = O_p{ D_w D_mu + (D_w + D_mu) rho_n + rho_n^2 }
# so after multiplication by sqrt(n) the quantities that must vanish are
#     sqrt(n) |T_eps| <~ sqrt(n) rho_n (D_w + rho_n)
#     sqrt(n) |T_a|   <~ sqrt(n) { D_w D_mu + (D_w + D_mu) rho_n + rho_n^2 }.
# Both are reported (i) evaluated at the REALISED mean D_w, D_mu of the cell, and
# (ii) in pure worst-case form with D ~ rho_n, which is the version the audit
# called vacuous.
# ---------------------------------------------------------------------------

A_ENT <- 2          # uniform-entropy constant in the Giné-Koltchinskii bound
DELTA_RULE <- function(n) 1 / n            # delta_n; log(4/delta_n) = log(4n)

VC_ANCHORS <- data.frame(L = c(10, 20, 30), v = c(185, 427, 686))


#' Calibrate v(Lbar, p) = c1 * Lbar * log2(Lbar p) + c0 to the audit's anchors
vc_calibration <- function(p = 5) {
  x <- VC_ANCHORS$L * log2(VC_ANCHORS$L * p)
  cf <- stats::coef(stats::lm(VC_ANCHORS$v ~ x))
  list(c0 = unname(cf[1]), c1 = unname(cf[2]), p = p)
}

#' v_{Lbar,p} from the calibrated O(Lbar log(Lbar p)) form
vc_dim_calibrated <- function(L, cal = vc_calibration()) {
  cal$c0 + cal$c1 * L * log2(L * cal$p)
}

#' v_{Lbar} from Bartlett-Jordan-McAuliffe Lemma 9 with d_tilde binary features
vc_dim_bjm <- function(L, d_tilde) {
  4 * d_tilde * L * log2(3 * L)
}


#' kappa_n and rho_n for one (n, Lbar) cell
localisation_radius <- function(n, v) {
  kappa <- v * log(A_ENT * exp(1) * n / v) + log(4 / DELTA_RULE(n))
  list(kappa = kappa, rho = sqrt(kappa / n))
}


#' Exact log-cardinality of the class of trees actually searched
#'
#' The standalone proof's localisation radius is r_n = sqrt(s_n / n) with
#' s_n = 1 + log|T_{Lbar,n}|. Here the searched class is exactly
#'   { binary trees over `n_feat` binary splits, depth <= `depth`, leaves <= L },
#' whose count satisfies the recursion
#'   N(d, 1) = 1,  N(0, k) = 0 for k >= 2,
#'   N(d, k) = n_feat * sum_{j=1}^{k-1} N(d-1, j) N(d-1, k-j).
#' Computed in log space (the counts overflow double precision well before
#' depth 7).
log_class_cardinality <- function(L, depth, n_feat) {
  logsumexp <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) return(-Inf)
    m <- max(x)
    m + log(sum(exp(x - m)))
  }
  logN <- matrix(-Inf, nrow = depth + 1L, ncol = L)   # rows: depth 0..depth
  logN[, 1] <- 0                                      # a single leaf
  for (d in seq_len(depth)) {
    for (k in 2:L) {
      terms <- vapply(seq_len(k - 1L),
                      function(j) logN[d, j] + logN[d, k - j],
                      numeric(1))
      logN[d + 1L, k] <- log(n_feat) + logsumexp(terms)
    }
  }
  logsumexp(logN[depth + 1L, ])
}


#' Assemble the full bound table for the design
bound_table <- function(cells, n_feat, cal = vc_calibration()) {
  out <- lapply(seq_len(nrow(cells)), function(i) {
    n <- cells$n[i]; L <- cells$L[i]; depth <- cells$depth[i]
    v_cal <- vc_dim_calibrated(L, cal)
    v_bjm <- vc_dim_bjm(L, n_feat)
    lr_cal <- localisation_radius(n, v_cal)
    lr_bjm <- localisation_radius(n, v_bjm)
    log_card <- log_class_cardinality(L, depth, n_feat)
    s_n <- 1 + log_card
    data.frame(n = n, L = L, depth = depth,
               v_calibrated = v_cal, kappa_n = lr_cal$kappa, rho_n = lr_cal$rho,
               v_bjm = v_bjm, rho_n_bjm = lr_bjm$rho,
               log_card = log_card, s_n = s_n, r_n = sqrt(s_n / n))
  })
  do.call(rbind, out)
}


#' sqrt(n) * lemma bound for T_eps and T_a, at given discrepancies
#'
#' `d_w`, `d_mu` are the realised (mean) population discrepancies of the cell;
#' pass `d_w = d_mu = rho` for the pure worst-case form.
lemma_bound_scaled <- function(n, rho, d_w, d_mu) {
  list(
    T_eps = sqrt(n) * rho * (d_w + rho),
    T_a = sqrt(n) * (d_w * d_mu + (d_w + d_mu) * rho + rho^2)
  )
}
