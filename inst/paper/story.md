# Paper Story — `interpretable-causal-inference-high-stakes` (doubletree)

**Instance path:** `inst/paper/story.md`
**Governed by:** `.claude/rules/paper-protocol.md`
**Status:** draft as of 2026-09-18; reconciled against the abstract and Introduction 2026-09-18 (see below).

## Reconciliation log (2026-09-18)

Checked every section above against the current abstract (`manuscript.tex:33`) and Introduction
(`:38,40,42,44`) after a whole-document audit found the abstract had drifted from the body during
the §3 restructuring. Fixed to match this file's framing:

- Abstract's "the oracle property" → "partition recovery" (this file's term; `oracle` collided
  with `xuStatisticalOptimalityOptimal2026`'s "oracle inequalities," a different object).
- Abstract's unsupported "sharper... than generic double machine learning" (false — both attain
  the same bound $V$) and "admits richer... structure than a GLM" (contradicted by the
  representational-trade paragraph now in `sec:discussion`, which declines to rank the two) →
  replaced with the true comparative claims: matches DML's efficiency without needing
  cross-fitting; represents structure a GLM must encode explicitly, without an unsupported
  superiority ranking.
- **S3 was entirely absent as a stated claim** from the abstract and contribution paragraph —
  only gestured at via a roadmap section-pointer. Added to both (abstract: one clause; ¶3 of the
  Introduction: one sentence) — this was the largest gap this reconciliation found.
- **S4** ("the anchor itself is a non-interpretable cross-fit object, so the audit trail covers
  the point estimate, not the width") was missing from both the abstract and ¶3 — both now state
  explicitly that the reported point estimate stays the two displayed trees regardless of the
  anchor.
- The scope-boundary disclaimer was previously stated twice, in three different terms, at full
  length (abstract + roadmap). Now stated once, in full, in the abstract; the roadmap points
  forward to §`sec:discussion` without re-litigating it, and only promises what that section
  currently contains (asymptotic + pointwise) — the uniformity and GLM/ML-ensemble comparisons in
  this file's Scope boundary section live in the abstract now, not (yet) restated in the
  Discussion body itself.
- ¶2's "not aware of a construction that delivers both" (vague — cross-fitted interpretable trees
  arguably deliver "both" in that narrow sense) → adopted this file's more precise "without
  splitting, at the efficiency bound," which names the actual missing ingredient.

## Framing

The modern semiparametric stack — AIPW, double machine learning, TMLE — earns its guarantees by
being indifferent to how the propensity score and outcome regression are fit: any learner meeting a
product-rate condition works, provided the own-observation empirical-process term is controlled by
cross-fitting or a Donsker restriction. That indifference is exactly what licenses pairing these
estimators with random forests or boosting, and treating the nuisance model as a means to an end.
In health care payment adjustment, treatment triage, and sentencing, the means *is* the object under
review: a regulator, judge, or auditor needs the calculation traceable step by step, and a black-box
nuisance fit does not supply that trail no matter how well it predicts.

The gap is not merely that black-box learners are opaque. Cross-fitting itself defeats auditability
even with an interpretable base learner: it injects randomness into the nuisance functions and
produces one fit per fold, so five-fold cross-fitting asks stakeholders to read five different
propensity models. Existing constructions deliver flexible nuisance estimation with asymptotic
guarantees, or an interpretable nuisance model — this paper is not aware of one delivering both,
without splitting, at the efficiency bound.

## Contribution

Doubletree estimates the ATT using two globally optimal, size-budgeted decision trees — one for the
propensity score, one for the outcome — each fit on every observation with no splitting or
cross-fitting, and shows that when the truth is exactly tree-representable within the analyst's
pre-specified leaf budget and cutpoint grid, this attains the semiparametric efficiency bound with
the ordinary plug-in variance; when it is not, an anchor interval keeps inference honest instead of
assuming the approximation error away.

## Scope boundary

No real or simulated payment, triage, or sentencing demonstration (a Medicaid serious-mental-illness
application is planned, not run); the recovery guarantee is pointwise in the data-generating process,
not uniform as the structure-defining signal shrinks; misspecification is treated as an
exact-versus-not dichotomy, not a smooth approximation regime; and no claim is made of matching an
unrestricted ML ensemble — the defended comparison is against a GLM working model.

## Major claims

- **S1:** Full-sample fitting is licensed not by the product-rate condition but by the candidate
  class being fixed and finite in advance, which is why no split is needed for either displayed tree.
- **S2:** Under structural sparsity, tree selection lands inside the sufficient partition class with
  probability approaching one, at an exponential rate, giving a root-n CLT with nominal coverage and
  no correction for having chosen the structure adaptively from the same data.
- **S3:** The estimator is regular and efficient in the *full* nonparametric model; it does not
  exploit the smaller bound available in the sparse submodel.
- **S4:** When sparsity fails, the anchor interval retains coverage regardless — but its width
  plateaus at the product of the two realized nuisance errors, and the anchor itself is a
  non-interpretable cross-fit object, so the audit trail covers the point estimate, not the width.
- **S5:** Even outside the class, the point estimate still targets something fixed: when each
  nuisance has a unique best tree-representable approximation within the budget, it is consistent for
  the corresponding best-sparse-tree effect.
- **S6:** The single-outcome-tree plug-in is algebraically doubletree with tied trees; under
  coarsening plus homoskedasticity it can be super-efficient, heteroskedasticity can reverse that,
  and its naive standard error under-covers.
- **S7:** The leaf budget is a single interpretability knob whose price is quantified exactly rather
  than asserted — smaller budget, more displayable tree, wider approximation-aware interval.
- **S8:** The underlying tree representation is not itself a limitation: at the saturated budget it
  represents any grid-measurable function, containing what *any* model restricted to the same
  grid-measurable inputs can compute — a generalized linear working model, but by the identical
  measurability argument also a random forest, boosting ensemble, or kernel smoother so restricted.
  Against the GLM specifically the containment is strict (a tree can represent a cross-coordinate
  departure from additivity on the working model's link scale that no coefficient choice can);
  against a flexible learner grown to purity on the same inputs it is generically an equality, not
  a strict improvement. Doubletree's small leaf budget is a *chosen* restriction for legibility, not
  a representational ceiling — the earlier framing ("neither representation is uniformly
  preferable") understated this and has been corrected.
  The cost of the small budget is bounded and already priced by S4/S7's machinery: a nuisance too
  additive-and-distributed to fit the budget is exactly a nuisance for which structural sparsity
  fails, absorbed into the anchor interval's width, not into validity.

## Reconciliation log, continued (2026-09-18, second pass)

Corrected a real error the author caught: the manuscript's original "trees vs. linear working
models is a two-sided trade, neither uniformly preferable" framing (`sec:lasso-comparison`, then
relocated to `sec:discussion`) understates trees' representational generality. General
(unconstrained) decision trees strictly dominate GLMs representationally — a saturated tree
represents any function constant on the analyst's grid atoms, which includes any GLM's fit,
strictly more so once the truth has any departure from additivity on the GLM's link scale.
Doubletree's small leaf budget is a *deliberate restriction for legibility*, not evidence that
trees are inherently limited. Added as `rem:saturated-general` in §2.2 (the mathematical core) plus
a shortened, corrected Discussion paragraph (the positioning argument, plus a properly-hedged
citation to `rudin2019` for the practical case that legibility need not cost accuracy). A
`domain-reviewer` audit caught and fixed three real issues in the first draft: an overclaim ("any
function on the grid" — contradicted by §2.2's own "$\cX$ need not be finite"), a link-scale
subtlety (log-link GLMs can represent what looks like an interaction on the identity scale), and
an unsupported empirical claim about "real" nuisance functions.

## Reconciliation log, continued (2026-09-25, third pass)

S8 broadened and corrected. `rem:saturated-general`'s containment claim was generalized from
GLM-specific to any model restricted to the same grid-measurable inputs (sound — a pure
measurability argument, no hidden GLM-specific step, per an independent `domain-reviewer` audit,
`quality_reports/reviews/2026-09-25_sec2.2-saturated-general-broadening.md`). That same audit
found the *strictness* claim does not generalize: against a flexible grid-restricted learner grown
to purity, the containment is generically an equality, not strict, since such a learner reproduces
exactly the atom-wise fit the saturated tree already delivers — a real error that had already
propagated into `claims.md`. Also fixed same day, both live only once the leaf-cost formula was
generalized from the binary-only $2^k$ to $\prod_i c_i$: a false coefficient-tie-robustness claim
(within-coordinate ties between adjacent categories *can* collapse leaves, unlike cross-coordinate
ties — vacuous under the old binary-only framing), and a self-contradictory stated reason for it
("not unions of grid atoms," when every level set of an atom-constant function is one by the
remark's own opening premise — corrected to "not *axis-aligned* unions"). `manuscript.tex`,
`claims.md`, and `notation.md` all updated to state the same scope consistently.

Separately, the same session found and fixed an unrelated, pre-existing ambiguity one paragraph
away: the control-weighted projection $\Pi^\nu_\tau$ never stated its reference measure, and a
newly-added explanatory sentence made the ambiguity load-bearing (under the wrong reading, a
result the paper already relies on — "these projections are exactly the squared-error
minimizers," `manuscript.tex:214` — would have been false). Resolved by author confirmation: $\nu$
is a measure, matching `theory.tex:206–207`'s $d\nu=\{1-e_0\}dP$ exactly, with no atom-mass
density factor. The atom-mass symbol $p_x$ this had relied on is removed from `manuscript.tex`
entirely as a result — it had no other use in the document.
