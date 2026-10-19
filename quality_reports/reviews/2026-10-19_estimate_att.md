# R Code Review: estimate_att.R (2026-10-19 delta, second pass)
**Date:** 2026-10-19 | **Reviewer:** r-reviewer agent | **Score: 93/100** (first pass 88)
Files: R/estimate_att.R (feasible_*/certificate_case_* return fields and docs), tests/testthat/test-estimate-att.R, man/estimate_att.Rd. Critical 0, High 0, Medium 0, Low 4.

The feasible_* doc is correct in substance and man/estimate_att.Rd matches. No regressions.

## Open issues (follow-ups, all Low)
- L1 feasible_* doc gives the wrong reason (the up-front check only shows a feasible tree exists; the returned tree is feasible because of the native solver floor).
- L2 "all three estimate_att* entry points" is now four. L3 seealso calls crossfit a "fallback" (the pipeline forbids that reading). L4 hidden `clip_c >= 0.01` floor (shared with the xgb review).

## Disposition (author)
Committed at 93/100.
