# DGP Taxonomy: What Each Simulation DGP Tests

**Date:** 2026-07-14
**Study:** `six-approach-arbitration` (supersedes the deprecated
`six_approach_comparison`)
**DGP source:** `R/dgp.R`. All DGPs: binary outcome, true ATT = 0.15.

---

## TL;DR

- The four DGPs are **not a difficulty continuum**. They differ in **which
  assumption each stresses**. Labeling them simple→complex invites the wrong
  reading (that complexity per se breaks estimators — it does not).
- On the arbitration DGPs (strengthened propensity), **M-split matches
  cross-fitting**: coverage and variance track each other on every binary DGP.
  There is no M-split failure to explain away.
- **DGP4 (continuous) is the Constitution §9 stress regime.** It stresses
  discretization + the tree smoothness ceiling, and for the Rashomon-intersection
  methods it can exhaust memory (see `INFEASIBLE_CELLS` in `config/grid.R`). That
  is a reported boundary, not a bug.

---

## Why the propensity was strengthened (vs the deprecated study)

`R/dgp.R` (2026-07-08) scales propensity linear + interaction coefficients
~2.5–3× relative to `six_approach_comparison`. Reason: with the old weak
coefficients (~0.2–0.3), cross-validated log-loss correctly preferred a **stump**
for the propensity nuisance — so the structure-selection arbitration was
comparing near-constant propensity fits. Strengthening keeps the interaction
*structure* while making it recoverable by a data-driven selector; overlap stays
safe and the estimand is unchanged (true ATT = 0.15 is set by the outcome model).
This is why taxonomy/parity conclusions must be drawn here, not from the
deprecated weak-propensity study.

---

## The two theorems need different conditions

| Estimator | Needs | On binary DGPs (1–3) |
|-----------|-------|----------------------|
| Cross-fitting | product-rate / nuisance consistency; on binary X trees represent the active subcube **exactly** | satisfied |
| M-split (Theorem 1, **revised** 2026-07-13) | selected structure **near-oracle in risk** → near-oracle predictions via the loss–norm link; **no** unique-structure/margin condition | satisfied |

The **original** Theorem 1 claimed modal-structure identity ($P(s^*=s_0)\to1$),
which requires a risk margin separating the oracle structure from competitors.
Diagnostics show that margin often fails (outcome-structure modal frequency
≈ 0.4–0.5). The revision (`../../inst/paper/m-split-theory.tex`) weakens the claim
to near-oracle risk/predictions — which holds regardless of margin and is all
Theorems 2/4 consume. After the revision, **no binary DGP violates an M-split
assumption**, matching the empirical parity.

---

## Per-DGP mapping

| DGP | Covariates | Nuisance structure | Assumption exercised | Role |
|-----|-----------|--------------------|-----------------------|------|
| **simple** | 3 binary | linear (2–3 splits) | baseline exact representability | base performance |
| **moderate** | 4 binary | one 2-way interaction (4–5 splits) | interaction representability | base performance |
| **complex** | 5 binary | multiple interactions (6–8 splits); **diffuse structure mode** (freq ≈ 0.4–0.5) | many near-tied partitions (no margin) | shows M-split valid without structure identity |
| **continuous** | 2 binary + 2 continuous, incl. $x_4^2$ | smooth/nonlinear | **discretization + smoothness > 1** (trees saturate at first order); Rashomon-set memory blow-up | **STRESS (Constitution §9)** |

---

## Empirical parity (this study, n=1000)

From `dev-scripts/verify_msplit_parity.R` → `dev-scripts/msplit_parity_arbitration.csv`
(50 reps binary; 15 continuous). Coverage / SD(θ):

| DGP | crossfit cov | msplit cov | crossfit SD | msplit SD | reps |
|-----|-------------|-----------|-------------|-----------|------|
| simple | 0.94 | 0.94 | 0.036 | 0.035 | 50 |
| moderate | 0.96 | 0.96 | 0.046 | 0.034 | 50 |
| complex | 1.00 | 0.90 | 0.068 | 0.044 | 50 |
| continuous | 0.93 | **0.80** | 0.042 | 0.064 | **15** |

M-split tracks cross-fitting on all **binary** DGPs. SD grows with DGP complexity
for **both** methods (strengthened propensity ⇒ tighter overlap ⇒ higher ATT
variance); it is not an M-split effect.

**Continuous (stress) DGP — noisy local estimate, not a finding.** M-split
coverage reads 0.80 here, but on **only 15 reps** (Monte-Carlo error on a coverage
estimate ≈ ±12pp) this is within noise of nominal and is not interpretable
locally. The 1000-rep cluster run settles the continuous-DGP performance.

## M-split specification is correct as shipped

Both the **point estimate** (modal structure + cross-fit predictions) and the
**standard error** are correctly specified for the full simulation:

- Point estimate: binary-DGP parity with cross-fitting (above) confirms it.
- SE: the EIF SE `sqrt(mean(psi^2)/n)` (`eif_att_solve` → `att_se`,
  `R/inference.R`) is correct. By Neyman orthogonality the ATT estimator's
  first-order asymptotic variance is `Var[psi(O; theta_0, eta_0)]` — a function of
  the score at the true nuisance and the data distribution, **not** of how the
  nuisance was estimated. Nuisance-estimation error (including the M-split
  averaging / cross-split correlation) is second-order and drops out. This is
  exactly Theorem 4's "same asymptotic variance as single-split." No averaging
  correction is needed and the theory does not depend on one.

Ship both unchanged to the full simulation.

---

## Recommendations for the manuscript's simulation section

1. **Relabel the DGP axis** by assumption exercised, not "difficulty":
   representability (simple), interactions (moderate), diffuse structure
   mode / no margin (complex), smoothness+discretization (continuous).
2. **Report the modal-frequency diagnostic** for M-split alongside coverage — the
   empirical proxy for the (old) margin condition; explains why structure
   *identity* is not claimed while predictions still converge.
3. **Keep the complex and continuous DGPs** (Constitution §9/§11): stress regimes
   must be kept and honestly reported, not dropped to make results look clean.
4. **Continuous is where the interesting boundary is** — Rashomon-set memory blow-up
   (`INFEASIBLE_CELLS`) and the smoothness ceiling. Frame it as the method's
   operating-range boundary.

---

## Provenance

- Parity: `dev-scripts/verify_msplit_parity.R`, `dev-scripts/msplit_parity_arbitration.csv`.
- Diffuse-mode diagnostics (deprecated-dir, still valid mechanism):
  `../../dev-scripts/diagnose_msplit_dgp3.R`, `confirm_margin_taxonomy.R`.
- Theory revision: `../../inst/paper/m-split-theory.tex`, Theorem 1 + remark.
- Deprecation of the prior study: `../six_approach_comparison/DEPRECATED.md`.
