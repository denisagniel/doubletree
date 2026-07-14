# DEPRECATED — superseded by `six-approach-arbitration/`

**Status:** Deprecated 2026-07-14. Do not build on this directory.
**Live study:** `../six-approach-arbitration/`

This directory (`six_approach_comparison/`) is the original six-approach study.
It has been superseded and should not be used for new work or for any numbers
that go into the manuscript.

## Why it was superseded

`six-approach-arbitration/` replaces it with a cleaner, correct harness:

1. **Real package calls, not reimplementations.** Arbitration's
   `R/estimators.R` calls `doubletree::estimate_att`,
   `doubletree::estimate_att_msplit`, etc. directly. This directory reimplemented
   the estimators in `code/estimators.R`, which could (and did) drift from the
   package.
2. **Depth caps present.** Arbitration threads `max_depth = 4L` through all
   paths, so the continuous DGP (DGP4) runs. In this directory the crossfit
   wrapper hangs on continuous covariates (unbounded-depth GOSDT blow-up).
3. **Strengthened propensity DGPs.** Arbitration's `R/dgp.R` scales propensity
   coefficients ~2.5–3× (2026-07-08) so CV recovers real propensity *structure*
   instead of a stump. This directory's weak-propensity DGPs made structure
   selection operate on near-constant propensity fits.
4. **Corrected Rashomon tolerance.** Arbitration uses the fixed theory tolerance
   `epsilon_n = log(n)/n` with `auto_tune_intersecting = FALSE`; this directory's
   `eps_n = 2*sqrt(log n/n)` + auto-tune config is invalidated by the
   structural-margin resolution (manuscript Cor.).
5. **Adds `single_tree` (Alt. A) + fidelity diagnostics** (`delta`,
   `intersection_nonempty`, cross-fit twin coverage).

## Do not trust the archived results here

`results/combined/summary_inference.csv` (July 8) is **stale**. It shows an
M-split failure on the complex DGP (coverage 0.90, 295/500 reps) that does **not
reproduce** on current code — a fresh scoped re-run gave coverage 0.98, 50/50
reps (see `ground_truth_current.csv` in this directory, retained for provenance).
Any conclusions about "M-split not working" traced to this stale file, not to the
estimator.

## What was retained from the investigation

- The `estimate_att_msplit.R` fix (dead `structure_selection` parameter) and the
  Theorem 1 revision in `inst/paper/m-split-theory.tex` are DGP-independent and
  apply regardless of which sim directory is used.
- `DGP_TAXONOMY.md` here is being re-homed to `six-approach-arbitration/` against
  the strengthened-propensity DGPs.
