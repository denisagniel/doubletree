## ============================================================================
## 02_population_and_eligibility.R
##
## Base cohort -> AAP eligibility (elig_aap) -> exposure A (aap_achieved).
##
## EVERY DESIGN ELEMENT IN THIS SCRIPT IS INHERITED FROM dual-bounds, VERBATIM
## (smidata inst/analyses/doubletree__application.yml sections 2-4, all
## `status: inherited`). This paper adds no eligibility criterion of its own:
##   population.time_zero     = INDEX_DT
##   population.eligibility   = elig_aap == 1 iff MSR_DEN >= 1 in msr_aap
##   base_cohort.source_table = larger_smi_covariates (381,018 rows)
##   exposure.window          = [INDEX_DT, INDEX_DT + 12mo]
##   exposure.definition      = A = 1 iff MSR_NUM/MSR_DEN >= 0.5
##
## Inheriting the design does NOT mean calling dual-bounds' code. The shared
## implementation is smidata::smi_select_aap_row(), which both projects call
## directly; nothing here reads anything from ~/RAND/rprojects/smi/dual-bounds.
##
## WHAT RUNS HERE TODAY: the msr_code gate below, and nothing else. Everything
## that touches data is written as the intended call in comments, because the
## reads cannot execute (Tier-0 blocked -- see README.md).
##
## Run time when unblocked: msr_aap is 19,012,819 rows. Filtering it to one MSR
## code before the join is what keeps smidata::smi_select_aap_row() tractable;
## expect minutes, not seconds.
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("02 -- population and eligibility")

## ---- gate: msr_aap cannot be read without knowing which MSR code is AAP ----
##
## This runs for real. msr_aap holds ~50-70 rows per patient spanning several
## distinct measures at roughly 12 periods each. Passing an unfiltered table
## into smidata::smi_select_aap_row() makes "the single nearest-month row"
## ambiguous, and it will abort on the tie -- correctly, but 19M rows later.
## Abort here instead, at the point the value is actually missing.
##
## SHARED open decision: identical question to dual-bounds', same files, same
## consequence. See OPEN_DECISIONS$msr_code$registry_ref.

msr_code <- open_value("msr_code")

if (is.null(msr_code)) {
  cli::cli_abort(c(
    "Cannot read {.val {SMI_KEYS$aap}}: the AAP {.field MSR} code is unknown.",
    "x" = OPEN_DECISIONS$msr_code$question,
    "i" = OPEN_DECISIONS$msr_code$note,
    "i" = "Resolve it ONCE, in {.val {OPEN_DECISIONS$msr_code$registry_ref}}, then
           mirror it into {.file application/_config.R}. Do NOT pass {.code NULL}
           through to {.fn smidata::smi_select_aap_row} -- an unfiltered table
           produces a nondeterministic exposure."
  ))
}

## ---- INTENDED PIPELINE (needs real data; not reachable today) --------------
##
## Everything below this line is the call sequence the gate above protects. It
## is written out rather than left vague so the server run is a copy, not a
## design exercise. It mirrors dual-bounds' 02 call-for-call, because the
## design is the same design.
##
## base_cohort <- smidata::smi_read(SMI_KEYS$covariates, columns = SMI_COLS$covariates)
##
## ## Assert immediately after reading, once, near the top:
## stopifnot(nrow(base_cohort) > 0L)
## if (anyDuplicated(base_cohort$ID) > 0L) {
##   cli::cli_abort("{.val {SMI_KEYS$covariates}} is not one row per patient.")
## }
## if (any(is.na(base_cohort$INDEX_DT))) {
##   cli::cli_abort("{.field INDEX_DT} is the time zero for every window; NA is not
##                   recoverable and must not be dropped silently.")
## }
## ## Confirmed row count, smidata snapshot 2026-09-15_2824080e. A mismatch means
## ## the source was replaced in place (it has no upstream archive), which is a
## ## finding, not a nuisance.
## if (nrow(base_cohort) != 381018L) {
##   cli::cli_warn("Expected 381,018 rows; got {nrow(base_cohort)}.")
## }
##
## ## Grain: one row per (patient x measurement period), filtered to ONE measure.
## aap_raw <- smidata::smi_read(SMI_KEYS$aap, columns = SMI_COLS$aap) |>
##   dplyr::filter(.data$MSR == msr_code)
## if (nrow(aap_raw) == 0L) {
##   cli::cli_abort("No {.val {SMI_KEYS$aap}} rows with MSR == {.val {msr_code}}.")
## }
##
## ## Grain change: one row per patient.
## ##
## ## tie_break is left at smidata's own default, "earlier", and is NOT an open
## ## decision in this project: smidata::smi_select_aap_row() aborts on any other
## ## value, so there is no choice being silently made here to flag. (dual-bounds
## ## registers msr_tie_break because its registry flags the default itself; the
## ## tie FREQUENCY is still worth knowing, and its 90_checks_tier1.R CHECK 2
## ## measures it -- for both papers at once, since the input is identical.)
## aap_selected <- smidata::smi_select_aap_row(
##   aap_raw,
##   index_dt = dplyr::select(base_cohort, "ID", "INDEX_DT"),
##   target_offset_months = 12L,
##   achievement_threshold = config_aap_achievement_threshold,
##   max_month_gap = Inf   # see 90_checks_tier1.R for the calibrating check
## )
##
## population <- base_cohort |>
##   dplyr::left_join(aap_selected, by = "ID") |>
##   dplyr::mutate(A = dplyr::if_else(.data$aap_achieved, 1L, 0L))
##
## ## The analysis cohort is the AAP-eligible subset. Report the attrition
## ## explicitly -- an eligibility filter that silently halves the cohort is the
## ## kind of thing that gets discovered in review.
## analysis_cohort <- dplyr::filter(population, .data$elig_aap)
## cli::cli_alert_info(
##   "Eligible: {nrow(analysis_cohort)} of {nrow(population)};
##    unmatched to msr_aap: {sum(!population$matched)}"
## )
##
## ## Both ATT entry points need interior treatment proportion and enough mass in
## ## BOTH arms; estimate_att() additionally needs at least
## ## leaf_budget * m_n CONTROL units, because its mu tree is fit on controls
## ## alone. Report the two arm counts here so a downstream estimator failure is
## ## traceable to the cohort rather than to the estimator.
## cli::cli_alert_info(
##   "A = 1 (AAP achieved): {sum(analysis_cohort$A == 1L)};
##    A = 0: {sum(analysis_cohort$A == 0L)}"
## )

if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("02_population_and_eligibility.R complete.")
}
