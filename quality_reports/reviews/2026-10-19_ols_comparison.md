# R Code Review: 07b_compare_ols.R + helpers/ols_comparison.R (second pass)
**Date:** 2026-10-19 | **Reviewer:** r-reviewer agent | **Score: 86/100** (first pass 72)
Critical 0, High 0, Medium 2, Low 8.

First-pass fixes confirmed: provenance check (md5/learner/seed types agree), parametric fits before tree refits, repetition-1 folds reused and later repetitions stratified like the anchor, flagship m_n/lambda_n passed to refits, se_method/n_parameters_m0 columns, completion guard, aliasing caveat corrected, input validation, diff_rep1, sound leakage and one-sided-clip tests.

## Open issues (follow-ups)
- M1 tree refits use the current config's leaf_budget/max_depth, which 07's provenance does not record: add them to 07's prov and assert in 07b.
- M2 print() of four data frames to the console (conventions): drop or replace with one cli summary.
- L3 logit aliasing zeroed silently while OLS aliasing is counted. L4 the paired SE ignores training-set overlap (Nadeau-Bengio correction). L5 clip test may not exercise the upper clip (set A[1] <- 0). L6 leakage test omits m0_inter. L7 nuisance_metrics() has no empty-fold guard. L8 missing provenance field gives "subscript out of bounds". L9 mixed clipping when the anchor is the tree crossfit. L10 accepted: doubletree::: internals (guarded), two statements on a line, lm formula path (lm.fit would be faster), README row for 07b.

## Disposition (author)
Committed at 86/100. Not yet run on the server.
