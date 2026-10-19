# R Code Review: application pipeline (06, 07, _config.R, winsorize_outcome.R, decision_report.R; second pass)
**Date:** 2026-10-19 | **Reviewer:** r-reviewer agent | **Score: 84/100** (first pass 73)
Critical 0, High 0, Medium 4, Low 8.

First-pass fixes confirmed: winsorization after the structural checks; cap attribute stripped; cutpoints attribute survives; manifest cap round-trips (including Inf for q = 1); sample-max statement removed from _config.R; fingerprint covers the new decision; cf_name defined before use; disclosure.R sourced before use; se_method order correct; provenance seed type matches 07b; decision-report q-mismatch abort correct when the record exists.

## Open issues (follow-ups)
- M1 Per-arm masking of capped counts can be undone by subtraction (total capped is about (1 - q)(n - 1) and n is printed): suppress both arms if either is small; use a helper in disclosure.R.
- M2 07's xgboost pre-check does not test the version, so on xgboost 3.x 07 dies after the multi-hour flagship fit: check < 3.0 up front.
- M3 decision_report accepts a 06 record with no winsor entry (pre-winsorization intermediates): abort on is.null(p06$winsor); add tests for mismatch and missing.
- M4 06's m_n message is stale (estimate_att() default is ceiling(sqrt(n)), so the real control requirement is about 16 * ceiling(sqrt(n))).
- L5 full-precision cap copied into att_provenance.csv (record signif(cap, 2)). L6 estimand string wrong when q = 1. L7 all.equal tolerance too loose for the Y == pmin(Y_raw, cap) check. L8 psi share column names inaccurate when k_min binds. L9 README stale (winsorizing row, no 07b row). L10 winsorize_upper docs promise non-negative input. L11 sensitivity runs overwrite outputs (output names do not carry q). L12 polish (long lines; extend test-script-sources).

## Disposition (author)
Committed at 84/100. Not yet run end to end on the server (06/07/07b are server-only).
