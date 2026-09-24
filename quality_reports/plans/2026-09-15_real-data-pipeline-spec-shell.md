# Plan: Real-Data Application Pipeline Spec (smidata-shaped, uniform template)

**SUPERSEDED 2026-09-16.** The PI decided doubletree reuses dual-bounds' real-data design
wholesale. See `2026-09-15_real-data-pipeline-spec.md` (same directory) for the current spec.
This file is retained unmodified as the historical record of the pre-decision blank state
described below.

**Date:** 2026-09-15
**Status:** SPEC-SHELL ONLY — §1-4 below are genuinely undetermined, not just unverified. This
is NOT a design discussion outcome; it is a record of the current blank state, written so the
eventual real discussion has a starting checklist rather than a blank page. Do not treat
anything in §1-4 as decided.
**Scope:** doubletree's real-data application (the package's methods/simulation content, which
is mature and well-documented, is out of scope — see `README.md` for that).
**Why this doc exists:** Third of five pipeline specs written this session against one uniform
8-point template. Unlike the first two (vector-incremental-effects applied paper, dual-bounds),
which had existing design material to consolidate, this one has none — confirmed by exhausting
every place a design might live (README, DESCRIPTION, vignettes, `inst/paper/`, `simulations/`,
`data/`) and finding no mention of SMI, OMH, Medicaid, or any real dataset anywhere in the repo.

---

## 1. Population / cohort

**UNKNOWN.** No population has ever been specified for doubletree's real-data application. No
file in this repo names a cohort, a data source, or even confirms this is the same NYS OMH SMI
population the other four papers use, rather than a different setting entirely.

## 2. Exposure/treatment

**UNKNOWN.** doubletree's estimator (`estimate_att()`) targets the Average Treatment Effect on
the Treated — so *some* binary treatment/exposure is implied by the method itself, but nothing
in this repo says what it would be in a real application.

## 3. Outcome

**UNKNOWN.** Same — the method supports binary or continuous outcomes (`outcome_type`
parameter), but no real outcome has been named.

## 4. Covariates

**UNKNOWN.** No covariate set has been specified.

## 5. Dataset/measure inventory

**None exists.** Zero dataset_keys, zero table names, zero column names anywhere in this repo
related to a real application. `data/` is an empty placeholder directory (`.gitkeep` only).

## 6. Relationship to other papers

doubletree is the estimator that a "doubletree paper" (per the handoff's list of planned papers
this session works from) would presumably apply to some real dataset. Whether that's the same
NYS OMH SMI population as vector-incremental-effects/health-affairs-case-mix-paper/dual-bounds,
and if so what treatment/outcome pair, is entirely open. If it IS the SMI population, the same
questions raised in the other four specs (which base table, which index-date rule, raw claims vs.
extract) would all apply here too, from scratch.

## 7. Open questions — this is effectively the entire spec

1. **What is the research question?** Is there a target application area at all yet (SMI
   Medicaid, or something else), or is this genuinely undecided?
2. **If SMI:** what population, what treatment, what outcome? (Everything in §1-4 depends on
   this.)
3. **Does this reuse the same base cohort question (D001, `larger_smi_covariates` vs.
   `aim1_smi_charac`) the other SMI papers are navigating, or is a different population/data
   source intended?**
4. **What does the tree-based ATT estimator need from the data that the other estimators don't?**
   (`estimate_att()`'s no-sample-splitting flagship variant requires "grid-exact sparsity" —
   worth understanding what covariate structure that implies before picking a real dataset,
   since it may constrain covariate dimensionality/discretization in a way the other papers'
   estimators don't.)

## 8. Current code state

**Pure methods package. No real-data application code exists.** Confirmed:
- R package structure (DESCRIPTION/NAMESPACE/R/tests), Imports `optimaltrees (>= 0.4.1)` as a
  hard dependency (not Suggests).
- Simulation infrastructure is mature: 3 DGPs × 4 methods (tree-DML, Rashomon-DML, forest-DML,
  linear-DML) × 3 sample sizes × 500 replications, with O2/SLURM distributed-computing setup for
  18,000+ replications in 30-60 minutes.
- `inst/paper/` holds the methods manuscript (`manuscript.tex`) — theory paper, no empirical
  application section drafted yet (per the one relevant grep hit found: a note about
  "rescoping `sec:discussion`'s 'left to future work' sentence to *real-data*" — i.e. the paper
  itself currently defers the real-data application to future work, consistent with nothing
  existing yet).
- House style already documented, shared with optimaltrees
  (`.claude/rules/r-code-conventions.md`).

**This spec cannot progress past this point without a design conversation** — unlike the other
four pipelines, there is no existing sketch, however rough, to consolidate or correct. The next
step is the conversation itself, not further investigation of this repo.
