# Illustration-cell pre-registration: searching for a pedagogical example

**Status:** pre-registered 2026-10-28, before any headroom number in this
document has been computed. Frozen until amended by a dated entry below.

**Purpose:** find ONE small, honest, pedagogical illustration of the
flexibility-vs-legibility tension `doubletree` is built to resolve (flexibility
as insurance against misspecification; a small leaf budget for legibility).
This is NOT a search for "the best" application and does NOT revise the
primary application.

**Explicitly separate from, and does not revise, the primary application**
(SAA -> winsorized service cost, `07_estimate_att.R`/`07b_compare_ols.R`). That
result stands as reported. Nothing here changes `config_leaf_budget`,
`design_matrix_columns()`, or any `OPEN_DECISIONS` entry that gates the
primary application's fingerprint.

**Hard cap: at most 10 cells are ever evaluated under this pre-registration.**
The count below (Tier 0: 2, Tier 1: 2, Tier 2: 5, reference/diagnostic: 2
uncounted) totals 9 candidates. Any amendment that would exceed 10 requires a
new dated entry in [Amendments](#amendments), not a silent edit above it.

**Why this document exists instead of a `92`-style automated screen:** the
candidate space was widened three times during planning (outcomes only ->
all 16 `smidata` measures in either role -> variable covariate sets), and
each widening was evaluated and mostly rejected on the same principle: a
screen that searches a space this large and reports whichever cell looks best
would be the post-selection-inference problem `07_estimate_att.R`'s header
already refuses, one level up. See [Rejected scope](#rejected-scope-and-why)
for what was proposed and declined.

---

## 1. Reference and free diagnostic (not candidates)

| ID | Cell | Role | Note |
|---|---|---|---|
| R | (SAA, winsorized service cost) | Reference, not a candidate | Already run. `xgb - additive_logit` headroom on `m0` ~= 0.2% of outcome MSE (`92_outcome_screen.R`). Calibrates every candidate's number against a measured null. Does not count against the cap. |
| R' | (SAA, log(1 + service cost)) | Free diagnostic | Pure transform of R, zero new plumbing, reuses the committed cohort. Resolves whether R's null was substantive additivity or a tail-masked signal (costs may be multiplicative in risk factors, in which case log-cost shows *less* headroom -- that outcome is itself informative). Illustration-only: log-cost is not a dollar-scale ATT. Does not count against the cap. |

Run R' first. It costs nothing and the paper needs the answer regardless of
anything else in this document.

## 2. Tier 0 -- free, already queued

Both read from `92_outcome_screen.R`'s existing output under the committed SAA
cohort and covariates. **If C1 clears (real headroom, `loworder_share`
comfortably above the floor in [5](#5-the-ranking-metric)), stop: the
illustration is found at zero additional data cost and no exposure switch.**

| ID | Cell | Covariate set | Justification |
|---|---|---|---|
| C1 | (SAA, `hosp_any_fuh30`) | Confirmed 15 **+ `VAR5D_BD_YN`, `VAR5D_SUD_YN`** (p=17) | Whole-cohort, no post-exposure conditioning under the missing-row convention. Bipolar and SUD are named confounders for psychiatric admission specifically, absent from the cost-outcome's X because the cost outcome never needed them. |
| C2 | (SAA, `aap_access`) | Confirmed 15, unchanged | Whole-cohort, contingent on recomputed `DEN >= 1` coverage being near-universal. Kept at the confirmed 15 deliberately, as a comparability anchor against R. Note the semantics caveat in [7](#7-aap-semantics--draft-escalation-note-not-yet-committed): this outcome reads AAP as a care-access achievement indicator, which is probably the correct reading, but see that section before reporting it. |

## 3. Tier 1 -- the pre-registered primary pick and its fallback

Chosen on substantive/clinical reasoning about where the nuisance functions
(propensity and outcome regression -- not treatment-effect heterogeneity) are
plausibly non-additive **in the covariate blocks we actually have**, not from
any headroom number. Full reasoning in
[6](#6-exposure-ranking-reasoned-not-measured).

| ID | Cell | Covariate set | Justification |
|---|---|---|---|
| P1 (primary) | (AAP-as-exposure, corrected semantics, x claims-derived binary "any psychiatric inpatient admission in [s+12mo, s+24mo)") | **Open -- see below** | Best promise-to-cost ratio of any exposure in the catalog: largest denominator (all adults), weakest population restriction (interaction axes survive, unlike SAA's homogeneous schizophrenia+2-fills group), baseline-condition eligibility (not event-anchored), partial `02` plumbing already in git history pre-2026-10-09. Outcome is claims-derived specifically to avoid a second measure-denominator gate (see [8](#8-why-claims-derived-not-a-second-quality-measure)). |
| P2 (fallback, only if P1 is flat) | (FUH30-as-exposure x same claims-derived outcome) | **Open -- see below** | Strongest non-additivity story in the whole catalog on both nuisance surfaces (follow-up receipt after discharge is conjunctive on provider access; non-receipt risk is compounding across SUD x social instability). Not circular with the outcome: distinct events, distinct windows (discharge precedes the exposure window; admission is ascertained afterward). Deferred behind P1 because it requires a patient-specific event-anchored window -- genuinely new machinery, not more of the same `02` plumbing AAP needs. |

**Both P1 and P2 have open items before they can be run, listed in
[9](#9-open-engineering-items-before-p1p2-can-run).** Neither is Tier-0 free.

## 4. Tier 2 -- lower priority, evaluate only if Tiers 0-1 return empty

These were pre-registered in an earlier planning round, before the
substantive exposure-ranking exercise in [6](#6-exposure-ranking-reasoned-not-measured)
identified P1/P2 as stronger candidates. They remain pre-registered (not
dropped) but are now explicitly lower priority.

| ID | Cell | Covariate set | Justification / contingency |
|---|---|---|---|
| C3 | (SAA, `msr_ccs`) | Confirmed 15, unchanged | Women-only restriction is immutable and pre-treatment (unlike a PQI-style comorbidity restriction). Note: this is CCS as an *outcome* under SAA. CCS was separately considered and rejected as a candidate *exposure* ([6](#6-exposure-ranking-reasoned-not-measured), rank 5) for an unrelated reason (unmeasured trauma-history confounder) -- the two roles are independent and this cell is not affected by that rejection. |
| C4 | (SAA, `pqi01` -- diabetes short-term complications admission) | 15 **+ `VAR5G_META_RISK_YN`, `VAR5G_OBES_YN`, `VAR9A_ANTIDM_ALL_YN`** (p=18) | Contingent on two free checks clearing first (see [10](#10-contingency-checks-for-c4-c7-free)): PQI-table outcome observability, and ascertainment timing of the diabetes restriction. Drop `VAR5G_PRE_DM_YN` (near-exclusive of a diabetic population). If restriction is defined by a medication flag, that flag must leave X (constant column). |
| C5 | (SAA, `pqi14` -- uncontrolled diabetes admission) | Identical to C4 | Same population and confounding argument as C4; run only if C4 clears. Nearly free once C4's plumbing exists. |
| C6 | (SAA, `pqi08` -- heart failure admission) | 15 **+ `VAR5G_OBES_YN`, `VAR5G_META_RISK_YN`, `VAR5G_TOBACCO_YN`**, with `VAR5G_HF_YN` used to define the baseline-eligible population and then **excluded from X** (p=18) | `VAR5G_HF_YN` makes the restriction pre-treatment by construction -- cleanest PQI rescue on the population side. Gate B (n per leaf) is the live risk: HF prevalence ~5-8% in SMI Medicaid x 24-month enrollment is thin. |
| C7 | (SAA, `pqi05` -- COPD/asthma admission, age 40-64) | 15 **+ `VAR5G_TOBACCO_YN`** + one declared age split (40-49 / 50-64) (p=17-18) | `VAR5G_TOBACCO_YN` materially improves (does not fully close -- ever/never is a crude proxy for pack-years) the smoking-confounding gap. COPD is common enough in this age band that Gate B is plausibly survivable. |
| ~~C8~~ | ~~(SAA, `pqi15` -- asthma, age 18-39)~~ | **Excluded, not pre-registered** | Smoking confounder fix doesn't rescue it: 18-39 is the age band with the worst 24-month continuous-enrollment retention in this cohort, and it dies at Gate B even at L=4. |

**Excluded from the entire catalog, not pre-registered at any tier:**
`msr_col`, `msr_bcs` (dropped by PI 2026-07-08, reason undocumented -- re-opening
needs its own PI decision, not a screen); `msr_smc`, `msr_smd` (absent from
both the measures registry and the HEDIS/CMIT crosswalk -- no pre-registrable
estimand); `msr_fuh07_v2`/`msr_fuh30_v2`-as-outcomes under a prior exposure
(condition on a post-exposure discharge); `msr_amm_acute`/`msr_amm_cont`-as-outcomes
(condition on a new treatment episode in the measurement window); `(AAP, aap_access)`
read as a lagged-persistence design (AAP achieved in year 1 predicting AAP
achieved in year 2) -- technically valid but a substantively vacuous
illustration, since the estimand is dominated by measure autocorrelation.

## 5. The ranking metric

Rank on **`loworder_share = (pairwise_logit_skill - additive_logit_skill) /
(xgboost_skill - additive_logit_skill)`**, not on raw headroom. Raw headroom
(`xgboost - additive`) is not comparable across cells whose covariate set
differs in size: adding columns raises raw headroom while simultaneously
growing the grid and reducing a small tree's ability to reach it, so ranking
on raw headroom would systematically reward the richest covariate set --
exactly the forking path this document exists to close off. Raw headroom
is retained only as a **feasibility floor**: non-trivial relative to its
paired fold-wise SE, below which a cell is not worth a real `optimaltrees`
fit, never a ranking key.

Report `p` and `2^p` on every row for reference. Use `L = 4` or `L = 8` for
every illustration cell (never 16): it is where Gate B's n-per-leaf
requirement is easiest to clear, where the displayed tree is most legible,
where the package's own simulations validated the method, and where
`optimaltrees`'s advantage over a greedy search is largest.

## 6. Exposure ranking, reasoned not measured

Decided from clinical/behavioral reasoning about where the **nuisance
surfaces** (not treatment-effect heterogeneity -- `estimate_att()` fits trees
to `e(x)` and `mu0(x)`) are plausibly non-additive in a covariate block we
actually have. A real interaction driven by an unmeasured variable produces
no headroom at all; this criterion does most of the sorting below.

| Rank | Exposure | Mechanism | Verdict |
|---|---|---|---|
| 1 | FUH30-as-exposure | Follow-up receipt after discharge is conjunctive (near-automatic with an established outpatient provider, depends on housing/transport/SUD jointly without one); non-receipt risk is compounding (recent discharge AND active SUD AND social instability -> high risk; any one alone -> modest). Both nuisances plausibly non-additive, in blocks we have. | Best story, most expensive (event-anchored window). P2. |
| 2 | AAP (corrected: general care access) | "Barriers bind only for the disengaged": social risk/transport matter little with a usual source of care, compound without one. Largest denominator, weakest population restriction. | Best promise-to-cost ratio. P1. |
| 3 | AMM-continuation | Antidepressant discontinuation interacts sharply with comorbidity (bipolar branch is qualitatively different: manic-switch risk on discontinuation). | Good substance, killed on feasibility: episode-anchored denominator, n collapses in an SMI cohort where most are already on antidepressants. |
| 4 | AMM-acute | Same mechanism, weaker (84-day window leaves less discontinuation variation). | Not pursued. |
| 5 | CCS-as-exposure | Clean immutable denominator (women 30+), but the real effect modifier in SMI women -- trauma history making exams aversive -- is not in any available covariate block. | Killed: real but unmeasured mechanism. |
| 6 | SAA | Already measured: additive, R headroom ~= 0.2%. Homogeneous denominator population (schizophrenia + 2 fills) is the population-level analogue of the covariate-curation problem that produced the null. | Already run; this is why P1/P2 exist. |

## 7. AAP semantics -- draft escalation note, not yet committed

**Not entered into `OPEN_DECISIONS` by this document.** Adding an entry to
`_config.R`'s `OPEN_DECISIONS` changes `config_fingerprint()` (it hashes
`lapply(OPEN_DECISIONS, ...)`), which would force a re-run of `06` for the
*primary* pipeline even though nothing about the SAA-cost application
actually changes. Per Oracle's explicit advice, this is a correctness
question about already-committed work and should reach the PI as its own
item, not be folded into this illustration search. Drafted here so the text
exists; committing it to `_config.R` is a separate, deliberate action.

> **Draft entry, for the PI, not yet in `_config.R`:**
> *Question:* `msr_aap`'s own data dictionary (2025-11-19, authoritative) names
> it HEDIS "Adults' Access to Preventive/Ambulatory Health Services" -- a
> general care-access measure, not antipsychotic-specific. `doubletree`'s
> pre-2026-10-09 exposure framing and `92_outcome_screen.R`'s `aap_access`
> candidate both read it as adherence-adjacent. The sibling `dual-bounds`
> paper inherited the same design and may carry the same reading.
> *Ask:* confirm which reading is correct, and if the dictionary is right,
> correct the framing in both papers' documentation (not necessarily the
> computed quantity, which may be unaffected -- this is a labeling question
> until shown otherwise).

## 8. Why claims-derived, not a second quality measure

Any outcome whose eligibility/denominator is ascertained in the outcome
window `[s+12mo, s+24mo)` is a candidate for exposure-dependent ascertainment:
claims-based coding (and, for some measures, the denominator condition
itself) can be affected by how much care the exposure caused the patient to
receive. This generalizes the diabetes-mediator concern raised for the PQI
cells ([4](#4-tier-2----lower-priority-evaluate-only-if-tiers-0-1-return-empty))
to *every* quality-measure outcome, not just PQI's. Using a **claims-derived
binary event** (any psychiatric inpatient admission in the outcome window,
read directly from the service-claims files, not from a measure's `MSR_DEN`)
eliminates the entire class of validity failure rather than mitigating it
per-cell, and is strictly better than C1 (`hosp_any_fuh30`) on this count:
same underlying quantity, but retires the "a missing measure row means no
discharge" assumption C1 is contingent on.

A continuous, category-disaggregated cost outcome (e.g., psychiatric-inpatient-only
dollars) was considered and rejected as the fix for R's tail-dominance
problem: it has the right intuition (summing heterogeneous cost components
averages toward additivity) but makes the actual power problem worse --
zero-inflation (a binary signal wearing a continuous costume) and worse
relative tails (one admission dwarfs a category mean far more than it
dwarfs the total-cost mean) than R already has.

## 9. Open engineering items before P1/P2 can run

Neither pairing is Tier-0 free. Before either is run:

- **P1 (AAP-as-exposure):** resurrect the pre-2026-10-09 `AAP_MSR_CODE`
  exposure-eligibility path in `02_population_and_eligibility.R` (partial
  history exists in git); confirm AAP's eligibility window is baseline-condition
  type, not event-anchored (per [smidata's measure catalog](#smidata-measure-catalog-reference),
  it should be); construct the claims-derived psychiatric-admission outcome
  (new helper, reading the service-claims files directly rather than a
  measure table); **assign and justify P1's covariate set** following the
  rule in [5](#5-the-ranking-metric)/[6](#6-exposure-ranking-reasoned-not-measured)
  -- a starting hypothesis is the confirmed 15 unchanged, since AAP's named
  mechanism (prior-utilization x social-risk) is already covered by two of
  the four existing blocks, but this needs to be stated and justified before
  any headroom is computed, not left implicit.
- **P2 (FUH30-as-exposure):** a patient-specific, event-anchored window (the
  discharge date, not `INDEX_DT`) -- genuinely new machinery relative to
  `DOUBLETREE_WINDOWS`/`smi_patient_windows()`'s cohort-level anchor. Same
  claims-derived outcome as P1. Same open covariate-set assignment. Buy this
  only if P1 comes back flat, per [3](#3-tier-1----the-pre-registered-primary-pick-and-its-fallback).

## 10. Contingency checks for C4-C7 (free)

Both gate the entire PQI tier and must be checked before any PQI headroom is
computed, regardless of covariate rescue:

1. **Outcome observability:** do the `msr_pqi*` tables carry rows for all
   patients, or only for measure-denominator-eligible patients? If only
   eligibles, restricting on a baseline covariate flag fixes the *population*
   but not the *outcome ascertainment* -- `MSR_NUM` is still only observed
   for patients the measure deemed eligible in the measurement year, which
   is the post-treatment ascertainment the rescue was supposed to escape. If
   so, every PQI cell needs a claims-derived admission (as in P1/P2) instead,
   and moves effectively to the same tier as P1/P2.
2. **Ascertainment timing of the comorbidity restriction:** is it determined
   from a baseline lookback (clean) or from claims within the measurement
   year (contaminated by exposure-driven utilization, for every PQI measure,
   not gradient by clinical plausibility -- see the full reasoning in prior
   planning rounds, available on request, not reproduced here).

## Evaluation order

1. R' (log-cost diagnostic). Free, run now.
2. C1, C2 (Tier 0) against R. **Stop if C1 clears.**
3. PQI contingency checks (#10), free, metadata-only.
4. Resolve P1's open items (#9); run P1. **Stop if P1 clears.**
5. C3 (Tier 2, no open items).
6. C4 -> C5 -> C7 -> C6, ordered by expected n, only if #10 clears and only
   if 1-5 are empty.
7. P2, only if P1 is flat.
8. Tier 3 (a new exposure beyond AAP/FUH30, chosen by the actionability test
   in [Rejected scope](#rejected-scope-and-why)) only if everything above
   returns empty -- and still capped by the 10-cell limit in this document's
   header.

## Rejected scope, and why

Two widenings were proposed and declined during planning; recorded here so
they are not re-proposed without reading this first.

**All 16 `smidata` measures as candidate outcomes, searched quantitatively
and automated.** Rejected: the enumerable-and-affordable space does not scale
with the catalog -- cost scales with distinct *exposures* plumbed (1-3,
ever), not with (exposure, outcome) pairs, and the catalog growing to 16
measures doesn't change that budget.

**All 16 measures in either role (240 nominal ordered pairs), searched and
ranked by headroom.** Rejected on a stronger ground than cost: searching a
wider space makes a credible illustration *less* likely, not more. The
expected maximum of ~240 noisy headroom statistics under a global null is
higher than the expected maximum of a handful, so a wide search can
manufacture an apparently-promising cell out of pure noise that a narrow,
substantively-motivated search cannot -- and any honest multiplicity
correction across 240 cells raises the effect size needed to count as real,
making the search self-defeating relative to its own goal.

**What was kept from that round:** the *actionability test* for whether a
measure is substantively plausible as an exposure at all ("does an
intervention exist that could set this measure to 1 without that
intervention being `cure their disease`" -- process/quality measures pass,
PQI-style adverse-event measures fail) and the observation that
`valid_as_outcome`/`valid_as_exposure` are independent properties of a
measure, not mirror images (`fuh30_followup` is invalid as an outcome under a
prior exposure but valid as an exposure in its own right). A reusable
per-measure validity catalog using these two properties was recommended as a
**separate artifact for `smidata`'s registry**, not for this document, and
is not built here.

## Smidata measure catalog reference

The full 16-measure catalog (adherence/medication-management, preventive
screening, follow-up-after-hospitalization, AHRQ PQI adverse-event measures,
two undocumented measures) that this pre-registration draws candidates from
was compiled during planning from `~/RAND/tools/smidata/inst/dictionary/smi_measures.yml`
and `inst/docs/dictionary/QUALITY_MEASURES.md`. Not reproduced here in full;
consult those files directly for any measure not already named above.

## Amendments

*(None yet. Any change to the frozen cell list, the hard cap, or the ranking
metric after this document's status date goes here, dated, with a reason --
never as a silent edit to the sections above.)*
