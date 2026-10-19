# R Code Review: estimate_att_crossfit_xgb.R (second pass)
**Date:** 2026-10-19 | **Reviewer:** r-reviewer agent | **Score: 88/100** (first pass 71)
Files: R/estimate_att_crossfit_xgb.R, tests/testthat/test-estimate-att-crossfit-xgb.R, DESCRIPTION, man/estimate_att_crossfit_xgb.Rd
Reviewer read files only (no shell). Critical 0, High 0, Medium 2, Low 7.

First-pass fixes confirmed correct: RNG restore (both seed-present and seed-absent cases, registered before first set.seed), version gate (< 3.0), NULL best_iteration error, explicit eval_metric, reserved params rejected, seed required and range-checked, clipping docs corrected, new leakage/RNG/reserved-params tests valid.

## Open issues (follow-ups)
- M1 Tests fail instead of skipping on xgboost >= 3.0 (CRAN ships 3.x); add a skip helper, a test of the refusal message, `xgboost (< 3.0.0)` in Suggests, and pinning instructions.
- M2 `hit_nrounds_max` misses best iterations in (nrounds_max - early_stopping_rounds, nrounds_max]; return the realised round count and flag on it.
- L3 `clip_c` < 0.01 passes validation, then psi_att() aborts after all K fits (also estimate_att(), pre-existing); validate `clip_c >= 0.01`.
- L4 `@return` omits `hit_nrounds_max`. L5 "validated against 1.7 and 2.x" is unsupported; say 1.7.x. L6 `params` not checked for being a named list. L7 RNG snapshot should be the first lines of the function. L8 untested branches (no-.Random.seed restore, stratified = FALSE, NA seed, valid_frac limits, hit_nrounds_max warning, K = NA). L9 nthread/nrounds_max not validated.

## Disposition (author)
Committed at 88/100 (gate 80). Follow-ups above are not yet done.
