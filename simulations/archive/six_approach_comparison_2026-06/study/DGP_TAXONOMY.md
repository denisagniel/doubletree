# DGP Taxonomy: What Each Simulation DGP Actually Tests

**Date:** 2026-07-13
**Purpose:** Map each DGP in the six-approach comparison to the specific paper
assumption it exercises, and record what the current-code simulations actually
show (as opposed to the stale archived results that motivated a "M-split is
broken" narrative).

---

## TL;DR

- On **current code**, M-split (Approach 5) matches cross-fitting (Approach 2)
  on all binary DGPs (1–3): coverage 0.92–0.98, |bias| < 0.008, no replication
  loss. The apparent DGP3 failure (coverage 0.90, 295/500 reps) lived only in the
  **stale archived `summary_inference.csv`** (July 8) and does not reproduce.
- The DGPs were implicitly treated as a **difficulty continuum** (simple →
  complex). That framing is misleading: complexity is not what breaks any
  estimator here. What varies across DGPs is **which assumption is stressed**,
  and — importantly — none of the binary DGPs actually violate an assumption that
  M-split needs (after the Theorem 1 revision; see below).

---

## The two theorems need different conditions

| Estimator | Needs | On binary covariates |
|-----------|-------|----------------------|
| Cross-fitting (main theorem) | product-rate / nuisance consistency; trees represent functions on the active subcube **exactly** (manuscript §, `manuscript.tex:130`) | satisfied on DGP1–3 |
| M-split (Theorem 1, **revised** 2026-07-13) | selected structure is **near-oracle in risk** → near-oracle in predictions via the loss–norm link; **no** unique-oracle-structure/margin condition | satisfied on DGP1–3 |

The **original** Theorem 1 claimed modal-structure identity ($P(s^*=s_0)\to1$),
which *does* require a risk margin separating the oracle structure from
competitors. That margin is violated whenever several partitions have near-equal
risk — which the diagnostics show is common for the outcome nuisance (modal
frequency ≈ 0.4–0.5). The revision (`m-split-theory.tex`, Theorem 1) weakens the
claim to near-oracle risk/predictions, which holds regardless of margin and is all
Theorems 2/4 consume. **After the revision, DGP3 no longer violates any M-split
assumption** — consistent with the empirical parity.

---

## Per-DGP mapping

DGP definitions: `code/dgps.R`. True ATT = 0.15 throughout.

| DGP | Covariates | Nuisance structure | Assumption exercised | Margin (unique oracle structure)? | Current-code result (n=1000) |
|-----|-----------|--------------------|-----------------------|-----------------------------------|------------------------------|
| **1 Simple** | 3 binary | linear (2–3 splits) | baseline: exact tree representability, few active coords | Moderate — small structure, still multiple equivalent orderings | crossfit cov 1.00 / msplit cov 0.96 |
| **2 Moderate** | 4 binary | one 2-way interaction (4–5 splits) | interaction representability | Weaker — more competing partitions | crossfit 0.98 / msplit 0.92 |
| **3 Complex** | 5 binary | multiple 2-way interactions (6–8 splits) | many active coords; **diffuse structure mode** (freq_m0 ≈ 0.4–0.5) | **Absent** — many near-tied partitions | crossfit 0.96 / msplit 0.98 |
| **4 Continuous** | 2 binary + 2 continuous, incl. $x_4^2$ | smooth/nonlinear | **discretization + smoothness $>1$** (trees saturate at first-order); tree-depth control | N/A (separate issue) | **NOT run** — sim wrapper hangs (depth blow-up); deferred to cluster |

### Notes

- **DGP3 is a stress test of the (old) margin condition, mislabeled as "harder."**
  Its distinguishing feature is not that it is harder to estimate, but that the
  outcome nuisance admits many near-equivalent partitions. Under the revised
  theory this is benign for inference (near-tied-in-risk ⇒ near-equal
  predictions), which is exactly what the current sims show.
- **DGP4 is a genuinely different regime**: continuous covariates + a quadratic
  term stress discretization and the smoothness ceiling (trees saturate at
  first-order smoothness). It also triggers a tree-depth blow-up in the sim's
  `estimate_att_crossfit`/related wrappers (no `max_depth` cap), causing hangs.
  This is a **software** issue (see project memory: `estimate-att-depth-cap-asymmetry`,
  `continuous-dgp-depth-cap`), not an M-split issue. Fix the depth cap in the
  wrappers before drawing DGP4 conclusions.

---

## Recommendations for the manuscript's simulation section

1. **Relabel the DGP axis.** Not "simple → complex difficulty," but "which
   assumption is exercised": representability (DGP1), interactions (DGP2),
   diffuse structure mode / no margin (DGP3), smoothness+discretization (DGP4).
2. **Report the modal-frequency diagnostic** for M-split alongside coverage. It
   is the empirical proxy for the (old) margin condition and explains *why*
   structure identity is not claimed — while predictions still converge.
3. **Keep DGP3.** Per RESEARCH_CONSTITUTION §9/§11, stress regimes must be kept
   and honestly reported; do not drop it to make results look clean. Its role is
   to demonstrate that M-split remains valid even when the structure mode is
   diffuse — a strength of the revised theory, not a failure.
4. **Fix DGP4 tooling** (depth cap) before including continuous-covariate
   results; currently it cannot complete locally.
5. **Regenerate the archived summary.** `results/combined/summary_inference.csv`
   (July 8) is stale and shows failures that no longer reproduce. Re-run the full
   grid on the cluster and replace it before any numbers go into the paper.

---

## Provenance

- Current-code parity: `ground_truth_current.csv` (this dir), via
  `ground_truth_diff.R` (approaches 2 & 5, n=1000, 50 reps, DGP1–3).
- Diffuse-mode diagnostics: `../../dev-scripts/diagnose_msplit_dgp3.R`,
  `confirm_margin_taxonomy.R`.
- Theory revision: `../../inst/paper/m-split-theory.tex`, Theorem 1 + remark.
