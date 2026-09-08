#' Full-Sample ATT with Two Optimal Trees (the paper's primary estimator)
#'
#' @description
#' Estimates the Average Treatment Effect on the Treated (ATT) exactly as
#' specified in \code{inst/paper/theory.tex} \strong{Part II} ("Full-sample
#' inference with two optimal trees", \code{\\label{sec:main}}): both nuisance
#' partitions \emph{and} both sets of leaf values are computed from all \eqn{n}
#' observations. There is \strong{no sample splitting, no cross-fitting, no
#' intersection of candidate sets, and no averaging over fits}. The resulting
#' in-sample nuisance predictions are fed to the unchanged efficient
#' influence function (EIF) solver, giving the plain Wald interval whose
#' asymptotic validity is Theorem \code{thm:main} / Corollary
#' \code{cor:variance}.
#'
#' This is the flagship entry point. For the two secondary paths see
#' \code{\link{estimate_att_crossfit}} (Part I, the \eqn{K}-fold cross-fitting
#' fallback that does \emph{not} require structural sparsity) and
#' \code{\link{estimate_att_rashomon}} (the superseded shared-Rashomon-structure
#' variant, excluded from the current theory).
#'
#' @details
#' \strong{Structure selection (\code{theory.tex} eq.~\code{eq:select}).} For
#' each nuisance \eqn{j \in \{e, \mu\}} the tree structure solves
#' \deqn{\hat\tau_j \in \arg\min_{\tau \in \mathcal{T}^{(j)}_n}
#'       \{R^{(j)}_n(\tau) + \lambda_n |\tau|\},}
#' where \eqn{\mathcal{T}^{(j)}_n} restricts to partitions of at most
#' \eqn{\bar L} = \code{leaf_budget} leaves, each carrying at least \eqn{m_n}
#' observations (\emph{control} observations when \eqn{j = \mu}). This is a
#' \strong{fixed} penalty \eqn{\lambda_n}, not a cross-validated one: the
#' \code{cv_regularization_adaptive()} machinery used by
#' \code{\link{estimate_att_crossfit}} is the wrong mechanism here, because it
#' enforces neither a leaf-count cap nor a minimum-leaf-mass floor. The right
#' mechanism is \code{optimaltrees::bisect_lambda_to_budget()}, which tries
#' \code{lambda_n} first (per \code{prop:parsimony}, "the budget binds
#' asymptotically", this alone already respects the budget with probability
#' tending to one) and only bisects upward on overshoot.
#'
#' Because the \eqn{\mu} tree is fit on the control subset alone
#' (\code{X[A == 0, ]}, \code{Y[A == 0]}), the leaf-mass floor on that fit
#' \emph{is} automatically a control-count floor, as \code{eq:feasible}
#' requires -- no group-specific floor argument is needed.
#'
#' \strong{Identifying assumption.} Part II's central limit theorem holds under
#' \strong{structural sparsity} (\code{ass:sparsity}): a tree of at most
#' \eqn{\bar L} leaves represents \emph{both} true nuisances exactly
#' (\eqn{\delta_e = \delta_\mu = 0}; \code{def:sufficient}). This is this
#' estimator's own identifying assumption and it is \emph{not} checkable from
#' data. \code{theory.tex} notes the cost is asymmetric: sparsity is cheap for
#' \emph{hierarchical} structure (active covariate set varying by region --- as
#' few as \eqn{s + 1} leaves) and expensive for \emph{additive} structure (a
#' leaf per configuration of the active coordinates). When sparsity plausibly
#' fails, use \code{\link{estimate_att_crossfit}}, which does not need it.
#'
#' \strong{No automatic fallback.} \code{estimate_att()},
#' \code{\link{estimate_att_crossfit}} and \code{\link{estimate_att_rashomon}}
#' are separate, explicit entry points. This function reports sparsity-proxy
#' \emph{diagnostics} (\code{certified_e}, \code{certified_m0}, etc.) but never
#' switches estimator on the basis of them: doing so would turn estimator
#' choice into data-dependent post-selection inference that the paper does not
#' analyse.
#'
#' \strong{Variance (\code{cor:variance}, eq.~\code{eq:vhat}).} The variance
#' estimator is the ordinary empirical-influence-function plug-in
#' \eqn{\hat V = n^{-1}\sum_i \hat\psi_i^2}, i.e. exactly what
#' \code{\link{att_se}} / the shared \code{eif_att_solve()} already compute; they
#' are reused unchanged.
#'
#' \strong{Known finite-sample caveat (not hidden).} \code{theory.tex}
#' (immediately after eq.~\code{eq:vhat}, and again in its appendix discussion)
#' observes that because \eqn{\hat\mu}, \eqn{\hat e} and \eqn{\hat\theta}
#' together consume \eqn{2\bar L + 1} fitted parameters, this in-sample plug-in
#' is \strong{downward biased by a relative factor of order}
#' \eqn{(2\bar L + 1)/n}, so the reported \code{sigma} is mildly optimistic and
#' the nominal 95\% interval can undercover in small samples. The theory
#' mentions a degrees-of-freedom divisor \eqn{n - (|\hat\tau_e| +
#' |\hat\tau_\mu| + 1)} as "a reasonable finite-sample convention", but this is
#' explicitly parenthetical --- it is \emph{not} required for the asymptotic
#' result (the bias vanishes precisely because \eqn{\bar L} is fixed in \eqn{n},
#' \code{ass:budget}). It is deliberately \strong{not} implemented here:
#' applying it would mean changing the shared \code{att_se()} /
#' \code{eif_att_solve()} used by the other two estimators, or duplicating
#' their logic. This function therefore reports the plain \eqn{n}-normalised
#' plug-in, matching \code{\link{estimate_att_crossfit}}'s and
#' \code{\link{estimate_att_rashomon}}'s convention exactly, so the three are
#' directly comparable. Callers wanting the convention can rescale:
#' \code{sigma * sqrt(n / (n - (n_leaves_e + n_leaves_m0 + 1)))}, using the
#' returned \code{n}, \code{n_leaves_e} and \code{n_leaves_m0}.
#'
#' \strong{Covariates must be binary (an implementation limitation, not a
#' theoretical one).} \code{inst/paper/manuscript.tex} \S\code{sec:trees} and
#' \code{theory.tex} Assumption~\code{ass:finite} both allow \eqn{\mathcal{X}}
#' to be a general (e.g.\ continuous) covariate space, discretised only via a
#' finite set of analyst-chosen cutpoints fixed before seeing the data;
#' \code{ass:finite} itself is explicit that \eqn{\mathcal{X}} is
#' \emph{not} required to be finite. This function does not yet perform that
#' cutpoint-based discretisation internally, and the restriction to already-
#' binary columns is a consequence of that gap, not of the theory:
#' \code{optimaltrees::bisect_lambda_to_budget()} verifies the \eqn{m_n} floor
#' by mapping rows of the \emph{supplied} \code{X} to leaves, while a tree fit
#' on discretised continuous covariates splits on \emph{threshold-indicator}
#' features, so the two disagree unless \code{X} is supplied already binarised
#' to the analyst's own grid atoms. Non-binary \code{X} is therefore rejected
#' up front with a pointer to \code{\link{estimate_att_crossfit}} (which has
#' no leaf-feasibility step and does accept continuous covariates) rather than
#' allowed to fail obscurely or, worse, silently mis-check feasibility.
#'
#' @param X Data.frame or matrix of covariates. Must be binary (0/1); see
#'   "Covariates must be binary" above.
#' @param A Integer or numeric vector of treatment (0/1).
#' @param Y Numeric vector of outcome. Binary (0/1) when \code{outcome_type} is
#'   \code{"binary"}; any numeric when \code{"continuous"}.
#' @param leaf_budget Integer, \strong{required, no default}. The leaf budget
#'   \eqn{\bar L}. Per \code{theory.tex} Assumption~\code{ass:budget} this is a
#'   fixed, analyst-chosen \emph{structural} parameter, constant in \eqn{n} ---
#'   there is no sensible package-wide default, so the caller must supply it.
#'   Larger \eqn{\bar L} makes \code{ass:sparsity} easier to satisfy but the
#'   fitted trees less interpretable and the finite-sample variance bias
#'   (order \eqn{(2\bar L + 1)/n}) larger.
#' @param outcome_type Character. \code{"binary"} (default) or
#'   \code{"continuous"}. Determines the loss used for the \eqn{\mu} tree
#'   (\code{"log_loss"} vs \code{"squared_error"}).
#' @param lambda_n Numeric penalty per leaf in \code{eq:select}. Default
#'   \code{NULL} resolves to \eqn{\log(n)/n}. That rate is what
#'   \code{prop:parsimony} needs: it requires \eqn{\lambda_n \to 0}
#'   (Assumption~\code{ass:global}) \emph{and} \eqn{n\lambda_n \to \infty}
#'   simultaneously, and \eqn{\log(n)/n} satisfies both. (The formula coincides
#'   numerically with \code{optimaltrees::select_epsilon_n()}, but that function
#'   returns a Rashomon \emph{tolerance} --- a different quantity from a
#'   selection penalty --- and is deliberately not called here.)
#' @param m_n Integer minimum observations per leaf, the \eqn{m_n} of
#'   \code{eq:feasible}. Default \code{1L}, matching
#'   \code{optimaltrees::bisect_lambda_to_budget()}'s own default: a
#'   \emph{trivial} floor that forbids empty leaves but imposes no binding
#'   minimum-leaf-mass constraint. That is sufficient, because
#'   Assumption~\code{ass:construct}(d) asks only for \eqn{m_n/n \to 0} and
#'   \eqn{\bar L\, m_n \lesssim n}, both satisfied trivially at \eqn{m_n = 1}.
#'   Raise it to stabilise leaf values at the cost of shrinking the feasible
#'   candidate set.
#' @param verbose Logical. Forwarded to \code{optimaltrees}. Default
#'   \code{FALSE}.
#' @param propensity_loss Character. Loss used to fit and select the
#'   propensity (\eqn{e}) tree: \code{"log_loss"} (default) or
#'   \code{"squared_error"}. Both are proper scoring rules, so for a
#'   \strong{fixed} partition they give the identical leaf estimate --
#'   \code{propensity_loss} can only change which partition
#'   \code{bisect_lambda_to_budget()} selects (log-loss's near-boundary split
#'   weighting, \eqn{\sim 1/[q(1-q)]}, versus squared-error's uniform
#'   weighting), not the leaf values themselves. This is an empirical,
#'   decision-relevant choice, not a validity question: the loss-norm link for
#'   log-loss holds via Pinsker's inequality with a universal constant and no
#'   boundedness assumption on the truth (uniform loss-boundedness for
#'   estimation comes from clipping the \emph{fitted} propensity, orthogonal
#'   to the fitting loss). See
#'   \code{quality_reports/specs/2026-08-20_propensity-loss-choice.md} in the
#'   \code{global-scholars} project for the ablation that motivates the
#'   default. The \eqn{\mu}/outcome tree's loss is unaffected -- it is always
#'   determined by \code{outcome_type}.
#' @param discretize_method Character. Passed through to \code{optimaltrees}.
#'   Default \code{"quantiles"}. With binary \code{X} (the only supported case)
#'   discretisation is a no-op, so this argument exists for signature
#'   consistency with \code{\link{estimate_att_crossfit}}.
#' @param discretize_bins Integer or \code{"adaptive"}. Passed through to
#'   \code{optimaltrees}. Default \code{"adaptive"}. See
#'   \code{discretize_method}.
#' @param ... Additional arguments forwarded to
#'   \code{optimaltrees::bisect_lambda_to_budget()} (and hence to
#'   \code{optimaltrees::fit_tree()}) for \emph{both} the \eqn{e} and \eqn{\mu}
#'   fits.
#'
#' @return A list with elements:
#'   \item{theta}{Point estimate \eqn{\hat\theta} of the ATT.}
#'   \item{sigma}{Wald standard error, the plain \eqn{n}-normalised EIF plug-in
#'     \eqn{\sqrt{n^{-1}\sum_i \hat\psi_i^2 / n}}. \strong{Mildly downward
#'     biased in finite samples} by a relative factor of order
#'     \eqn{(2\bar L + 1)/n}; see "Known finite-sample caveat" in Details.}
#'   \item{ci_95}{Wald 95\% confidence interval \eqn{\hat\theta \pm 1.96\,
#'     \sigma}. Asymptotically valid under \code{ass:sparsity}; inherits the
#'     same finite-sample downward variance bias, so slight undercoverage at
#'     small \eqn{n} is expected.}
#'   \item{ci_95_wald}{Identical to \code{ci_95}, retained so callers can
#'     compare display and Wald intervals uniformly across all three
#'     \code{estimate_att*} entry points.}
#'   \item{score_values}{Influence-function values \eqn{\hat\psi_i} at
#'     \eqn{\hat\theta}.}
#'   \item{nuisance_fits}{List with \code{e_model}, \code{m0_model} (the fitted
#'     \code{OptimalTreesModel} objects) and the in-sample \code{propensity}
#'     (clamped) and \code{outcome_control} predictions, mirroring
#'     \code{\link{estimate_att_crossfit}}'s diagnostic fields.}
#'   \item{n}{Sample size.}
#'   \item{leaf_budget, lambda_n, m_n}{The resolved tuning values actually used
#'     (\code{lambda_n} is the resolved \eqn{\log(n)/n} when the argument was
#'     \code{NULL}).}
#'   \item{certified_e, certified_m0}{Logical. \code{TRUE} iff the returned fit
#'     is \emph{proven} to solve \code{eq:select} exactly (see
#'     \code{optimaltrees::bisect_lambda_to_budget()}).}
#'   \item{used_search_e, used_search_m0}{Logical. \code{FALSE} iff
#'     \code{lambda_n} alone already respected \code{leaf_budget}, so no
#'     bisection was needed --- the asymptotically typical case under
#'     \code{prop:parsimony}.}
#'   \item{n_leaves_e, n_leaves_m0}{Realised leaf counts \eqn{|\hat\tau_e|},
#'     \eqn{|\hat\tau_\mu|}.}
#'   \item{gap_e, gap_m0}{Computable suboptimality slack when the corresponding
#'     \code{certified_*} is \code{FALSE} but the fit is feasible and within
#'     budget; \code{0} when certified, \code{NA} when infeasible or the budget
#'     could not be met.}
#'
#' \strong{How to read the diagnostics.} \code{certified_* == TRUE} together
#' with \code{used_search_* == FALSE} and a small, \eqn{n}-stable
#' \code{n_leaves_*} is the operational signature that \code{prop:parsimony}'s
#' conclusion has kicked in, which is \emph{consistent with}
#' \code{ass:sparsity} holding at this \eqn{\bar L}. It is \strong{not} a proof
#' or a test of it: \code{ass:sparsity} is an assumption about \eqn{P}, is not
#' identified from a single sample, and a certified fit can occur while
#' \eqn{\delta_e, \delta_\mu > 0}. Treat these fields as a signal to reconsider
#' \eqn{\bar L} or to switch deliberately to
#' \code{\link{estimate_att_crossfit}}, never as a licence to auto-switch.
#'
#' @references
#' \code{inst/paper/theory.tex}: Section \code{sec:main} (full-sample
#' construction), eq.~\code{eq:select} (structure selection),
#' eq.~\code{eq:feasible} (leaf-mass feasible set), Assumption
#' \code{ass:budget} (fixed leaf budget), Assumption \code{ass:sparsity}
#' (structural sparsity), Proposition \code{prop:parsimony} (the budget binds
#' asymptotically), Theorem \code{thm:main} (CLT), Corollary
#' \code{cor:variance} / eq.~\code{eq:vhat} (variance estimation).
#'
#' @seealso \code{\link{estimate_att_crossfit}} for the Part I cross-fitting
#'   fallback that does not require structural sparsity;
#'   \code{\link{estimate_att_rashomon}} for the superseded shared-structure
#'   variant; \code{\link{att_se}} and \code{\link{att_ci}} for the shared
#'   variance/interval kernel.
#'
#' @examples
#' \dontrun{
#' library(optimaltrees)
#' set.seed(42)
#' n <- 400
#' # Hierarchically sparse nuisances: e depends on X1 only, m0 on X2 only,
#' # so a 2-leaf tree represents each exactly (ass:sparsity holds at Lbar = 4).
#' X <- data.frame(
#'   X1 = rbinom(n, 1, 0.5), X2 = rbinom(n, 1, 0.5), X3 = rbinom(n, 1, 0.5)
#' )
#' A <- rbinom(n, 1, ifelse(X$X1 == 1, 0.65, 0.35))
#' Y <- rbinom(n, 1, 0.25 + 0.30 * X$X2 + 0.15 * A)
#'
#' fit <- estimate_att(X, A, Y, leaf_budget = 4L)
#' fit$theta
#' fit$ci_95
#'
#' # Sparsity-proxy diagnostics (a signal, never an auto-switch):
#' c(fit$certified_e, fit$certified_m0)
#' c(fit$n_leaves_e, fit$n_leaves_m0)
#' }
#' @export
estimate_att <- function(X, A, Y, leaf_budget,
                         outcome_type = c("binary", "continuous"),
                         lambda_n = NULL, m_n = 1L, verbose = FALSE,
                         propensity_loss = c("log_loss", "squared_error"),
                         discretize_method = "quantiles",
                         discretize_bins = "adaptive",
                         ...) {
  outcome_type <- match.arg(outcome_type)
  propensity_loss <- match.arg(propensity_loss)
  check_att_data(X, A, Y, outcome_type = outcome_type)
  if (is.matrix(X)) X <- as.data.frame(X)
  n <- nrow(X)

  # -- leaf_budget: Lbar of ass:budget. Fixed, analyst-chosen, no default. -----
  if (missing(leaf_budget) || is.null(leaf_budget)) {
    stop("`leaf_budget` is required and has no default. It is the leaf budget ",
         "Lbar of theory.tex Assumption ass:budget: a FIXED, analyst-chosen ",
         "structural parameter (constant in n), so no package-wide default is ",
         "meaningful. Supply e.g. leaf_budget = 4L.", call. = FALSE)
  }
  if (!is.numeric(leaf_budget) || length(leaf_budget) != 1 ||
      is.na(leaf_budget) || leaf_budget < 1 ||
      leaf_budget != as.integer(leaf_budget)) {
    stop("leaf_budget must be a single positive integer, got: ", leaf_budget,
         call. = FALSE)
  }
  leaf_budget <- as.integer(leaf_budget)

  # -- m_n: minimum leaf mass of eq:feasible. ---------------------------------
  if (!is.numeric(m_n) || length(m_n) != 1 || is.na(m_n) || m_n < 1 ||
      m_n != as.integer(m_n)) {
    stop("m_n must be a single positive integer, got: ", m_n, call. = FALSE)
  }
  m_n <- as.integer(m_n)

  # -- lambda_n: FIXED penalty of eq:select, NOT CV-selected. -----------------
  # Default log(n)/n is the prop:parsimony-compatible rate: it needs
  # lambda_n -> 0 (ass:global) AND n*lambda_n -> infinity at the same time.
  # Computed inline on purpose: optimaltrees::select_epsilon_n() happens to use
  # the same formula but returns a Rashomon TOLERANCE, a different quantity.
  if (is.null(lambda_n)) {
    lambda_n <- log(n) / n
  } else if (!is.numeric(lambda_n) || length(lambda_n) != 1 ||
             !is.finite(lambda_n) || lambda_n <= 0) {
    stop("lambda_n must be NULL (default log(n)/n) or a single positive ",
         "finite number, got: ", lambda_n, call. = FALSE)
  }

  # -- discretize_method: same allowed set as estimate_att_crossfit(). --------
  valid_methods <- c("quantiles", "median")
  if (!discretize_method %in% valid_methods) {
    stop("discretize_method must be one of: ",
         paste(valid_methods, collapse = ", "), ", got: ", discretize_method,
         call. = FALSE)
  }

  # -- Covariates must be binary (see roxygen "Covariates must be binary"). ---
  # Not cosmetic: bisect_lambda_to_budget() checks the m_n floor by mapping rows
  # of the SUPPLIED X to leaves, but a tree fit on discretised continuous
  # covariates splits on threshold-indicator features, so the check either errors
  # ("split references feature index k but X has p columns") or, for non-binary
  # integers, silently mis-assigns leaves and mis-verifies eq:feasible.
  non_binary <- !vapply(X, function(col) {
    (is.numeric(col) || is.logical(col)) && all(col %in% c(0, 1))
  }, logical(1))
  if (any(non_binary)) {
    stop("estimate_att() requires binary (0/1) covariates; column(s) ",
         paste(names(X)[non_binary], collapse = ", "), " are not. ",
         "This is a current implementation limitation, not a requirement of ",
         "the theory: theory.tex Assumption ass:finite explicitly allows the ",
         "covariate space to be general (e.g. continuous), discretised via a ",
         "finite set of analyst-chosen cutpoints; this function does not yet ",
         "perform that discretisation internally. Operationally, the leaf-mass ",
         "check inside optimaltrees::bisect_lambda_to_budget() maps rows of ",
         "the supplied X to leaves, which is incompatible with a tree fit on ",
         "discretised threshold indicators. Either discretise X to 0/1 ",
         "indicators yourself, or use estimate_att_crossfit(), which has no ",
         "leaf-feasibility step and accepts continuous covariates.",
         call. = FALSE)
  }

  # -- Sufficient treated/control mass. ---------------------------------------
  # No K-fold framing here: both trees are fit on all n. The e tree needs both
  # arms present; the mu tree is fit on controls only, so eq:feasible's
  # Lbar * m_n <~ n becomes a CONTROL-count requirement for that fit.
  n_treated <- sum(A == 1)
  n_control <- sum(A == 0)
  if (n_treated < 1) {
    stop("No treated units (A = 1); the ATT is not defined. ", call. = FALSE)
  }
  if (n_control < 1) {
    stop("No control units (A = 0); the control-outcome tree cannot be fit.",
         call. = FALSE)
  }
  if (n_control < leaf_budget * m_n) {
    stop("Insufficient control units for leaf_budget = ", leaf_budget,
         " at m_n = ", m_n, ". The control-outcome tree is fit on controls ",
         "only, so eq:feasible requires at least leaf_budget * m_n = ",
         leaf_budget * m_n, " control units, got: ", n_control, ". ",
         "Either reduce leaf_budget or m_n, or collect more data.",
         call. = FALSE)
  }
  if (n < leaf_budget * m_n) {
    stop("Insufficient observations for leaf_budget = ", leaf_budget,
         " at m_n = ", m_n, ": eq:feasible requires at least ",
         leaf_budget * m_n, " observations, got: ", n, ".", call. = FALSE)
  }

  # Treatment proportion must be interior (redundant with the counts above, kept
  # as explicit validation exactly as on the cross-fit path).
  pi_hat <- mean(A)
  if (pi_hat <= 0 || pi_hat >= 1) {
    stop("Invalid treatment proportion: pi_hat = ", pi_hat,
         ". This should not happen after sample size validation.",
         call. = FALSE)
  }

  # == Nuisance fits: BOTH on all n observations (sec:main). ==================
  # Propensity e(X): all rows, treated and control. loss_function is
  # `propensity_loss` (default log_loss); both are proper scoring rules, so
  # for a FIXED partition they give the identical leaf estimate -- the only
  # way they can differ is in which splits/partition get selected. See
  # quality_reports/specs/2026-08-20_propensity-loss-choice.md.
  e_fit <- optimaltrees::bisect_lambda_to_budget(
    X, A,
    leaf_budget = leaf_budget,
    lambda_n = lambda_n,
    m_n = m_n,
    loss_function = propensity_loss,
    discretize_method = discretize_method,
    discretize_bins = discretize_bins,
    verbose = verbose,
    ...
  )

  # Control outcome m0(X) = mu: CONTROL SUBSET ONLY. This restriction is what
  # makes the m_n floor on this fit automatically a control-count floor, as
  # eq:feasible requires for j = mu -- hence no group=/group_value= needed.
  X0 <- X[A == 0, , drop = FALSE]
  Y0 <- Y[A == 0]
  loss_outcome <- if (outcome_type == "continuous") "squared_error" else "log_loss"
  m0_fit <- optimaltrees::bisect_lambda_to_budget(
    X0, Y0,
    leaf_budget = leaf_budget,
    lambda_n = lambda_n,
    m_n = m_n,
    loss_function = loss_outcome,
    discretize_method = discretize_method,
    discretize_bins = discretize_bins,
    verbose = verbose,
    ...
  )

  # == In-sample predictions for ALL n rows. ==================================
  # Both models are queried against the FULL X. There is no held-out/held-in
  # distinction on this path: everything is in-sample by construction, and that
  # in-sample-ness is exactly what Part II's theory analyses and licenses.
  models <- list(e_model = e_fit$fit, m0_model = m0_fit$fit,
                 outcome_type = outcome_type)
  pred <- predict_nuisances_fold(models, X, fold_rows = seq_len(n))

  # predict_nuisances_fold() does NOT clamp (on the cross-fit path the clamp
  # lives one level up, in get_fold_specific_eta()). Apply the same shared
  # bounds here so psi_att()'s propensity precondition holds; ass:construct(c)
  # requires clipping anyway, and in-sample leaf refits can hit 0 or 1 exactly
  # under perfect within-leaf separation.
  e_hat <- pmax(.PROPENSITY_LOWER_BOUND, pmin(.PROPENSITY_UPPER_BOUND, pred$e))
  m0_hat <- pred$m0

  # == Shared EIF solve (cor:variance / eq:vhat), reused UNCHANGED. ===========
  .att <- eif_att_solve(Y, A, e_hat, m0_hat, n)

  nuisance_fits <- list(
    e_model = e_fit$fit,
    m0_model = m0_fit$fit,
    propensity = e_hat,
    outcome_control = m0_hat
  )

  list(
    theta = .att$theta,
    sigma = .att$sigma,                # plain n-normalised EIF plug-in
    ci_95 = .att$ci_95,                # Wald interval
    ci_95_wald = .att$ci_95,           # same, named for cross-estimator parity
    score_values = .att$score_values,
    nuisance_fits = nuisance_fits,
    n = n,
    leaf_budget = leaf_budget,
    lambda_n = lambda_n,
    m_n = m_n,
    # Sparsity-proxy diagnostics: reported, NEVER acted on automatically.
    certified_e = e_fit$certified,
    certified_m0 = m0_fit$certified,
    used_search_e = e_fit$used_search,
    used_search_m0 = m0_fit$used_search,
    n_leaves_e = e_fit$n_leaves,
    n_leaves_m0 = m0_fit$n_leaves,
    gap_e = e_fit$gap,
    gap_m0 = m0_fit$gap
  )
}
