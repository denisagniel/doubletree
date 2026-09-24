## ============================================================================
## 01_declare_requirements.R
##
## Declare, against smidata's contract, every dataset and column this pipeline
## reads. Runs FIRST so a missing column is reported here, once, rather than
## 200 lines into stage 03.
##
## This paper requires NO datasets beyond the ones dual-bounds already requires:
## larger_smi_covariates, msr_aap, and all 36 cost-claims keys. Its only
## additional construction -- prior_cost's quartile discretization -- is computed
## from data already being read. The declaration below is therefore the same 38
## requirements dual-bounds declares, and that is a fact worth checking rather
## than assuming: if this list ever grows, the "reuses dual-bounds' design
## wholesale" claim in the registry has quietly stopped being true.
##
## Runs locally with a clear message when smidata is absent -- it must not
## abort, because its whole job is to report the state of the contract.
##
## Run time: seconds (contract metadata only; reads no data).
## ============================================================================

## Resolve this directory whether the script is run from the repo root or from
## application/ itself, so sourcing it interactively always works.
app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "covariate_blocks.R"))
if (config_echo_decisions) echo_open_decisions()

cli::cli_h1("01 -- declare requirements")

if (!config_has_smidata) {
  cli::cli_alert_info(
    "{.pkg smidata} is not installed on this machine; skipping contract verification."
  )
  cli::cli_alert_info(
    "Install from {.path ~/RAND/tools/smidata} to check the
     {length(cost_dataset_keys()) + 2L} dataset requirements."
  )
} else {
  cli::cli_alert_info("Environment: {.val {smidata::smi_env()}}")

  contract <- smidata::smi_contract()
  cli::cli_alert_info("Contract snapshot: {.val {contract$snapshot_id}}")

  ## ---- required datasets, one declaration per dataset ---------------------

  requirements <- c(
    stats::setNames(list(SMI_COLS$covariates), SMI_KEYS$covariates),
    stats::setNames(list(SMI_COLS$aap), SMI_KEYS$aap),
    stats::setNames(
      rep(list(SMI_COLS$cost_claims), length(cost_dataset_keys())),
      cost_dataset_keys()
    )
  )

  cli::cli_alert_info("Declaring {length(requirements)} dataset requirement{?s}.")

  ## Which required datasets does the resolved contract actually contain?
  ## Checked BEFORE the smi_require() loop, and for all 38 at once, because
  ## smi_require() aborts on the first absent dataset -- which on a laptop
  ## resolving smidata's bundled example contract means reporting one missing
  ## dataset instead of the real answer, "this is the wrong contract". This is
  ## a precondition check, not error suppression: the datasets that ARE present
  ## are still declared for real below, and the absent ones are named in full.
  needed_keys <- names(requirements)
  present_keys <- intersect(needed_keys, contract$dataset_keys)
  absent_keys <- setdiff(needed_keys, contract$dataset_keys)

  if (length(absent_keys) > 0L) {
    cli::cli_alert_warning(
      "{length(absent_keys)} of {length(needed_keys)} required datasets are NOT in
       snapshot {.val {contract$snapshot_id}}."
    )
    cli::cli_alert_warning("Absent: {.val {absent_keys}}")
    cli::cli_alert_warning(
      "Snapshot contains only: {.val {contract$dataset_keys}}. On a laptop this is
       smidata's bundled EXAMPLE contract -- the real one is
       {.val 2026-09-15_2824080e}, resolvable on the server."
    )
  }

  ## Same reasoning one level down: a dataset can be present under the right
  ## KEY while carrying an entirely different SCHEMA. smidata's example
  ## contract's msr_aap is exactly that case (MSR_ID/MSR_DENOM, not
  ## MSR/MSR_YR/MSR_DEN), so declaring against it aborts on columns and again
  ## reports the wrong headline. Check first, report the schema mismatch, and
  ## declare only where the contract can actually answer.
  declarable <- character()
  for (key in present_keys) {
    have_cols <- vapply(contract$datasets[[key]]$columns, function(col) col$name, character(1))
    want_cols <- requirements[[key]]
    missing_cols <- setdiff(want_cols, have_cols)
    if (length(missing_cols) > 0L) {
      cli::cli_alert_warning(
        "{.val {key}} is present but its SCHEMA does not match: missing
         {.field {missing_cols}} (has {.field {have_cols}})."
      )
    } else {
      declarable <- c(declarable, key)
    }
  }

  for (key in declarable) {
    smidata::smi_require(key, columns = requirements[[key]])
  }

  cli::cli_alert_success(
    "Declared {length(declarable)} of {length(needed_keys)} requirement{?s} against
     snapshot {.val {contract$snapshot_id}}."
  )
  if (length(declarable) < length(needed_keys)) {
    cli::cli_alert_warning(
      "This run verified NOTHING about the real data. Re-run on the server, where
       {.val 2026-09-15_2824080e} resolves."
    )
  }

  ## ---- the 12 covariate columns are a SEPARATE declaration ----------------
  ##
  ## SMI_COLS$covariates is (ID, INDEX_DT) -- the join keys. The 12 confirmed
  ## _YN covariates are declared here rather than folded into that vector so
  ## that a failure names which kind of column is missing: a missing INDEX_DT
  ## breaks every window in the pipeline, whereas a missing _YN column breaks
  ## the covariate selection specifically and is a different conversation.

  if (SMI_KEYS$covariates %in% declarable) {
    smidata::smi_require(
      SMI_KEYS$covariates,
      columns = unlist(covariate_blocks(), use.names = FALSE)
    )
    cli::cli_alert_success(
      "Declared the 12 confirmed {.field _YN} covariate columns on
       {.val {SMI_KEYS$covariates}}."
    )
  } else {
    cli::cli_alert_warning(
      "Skipped declaring the 12 confirmed {.field _YN} covariates:
       {.val {SMI_KEYS$covariates}} is not declarable against this contract."
    )
  }

  ## ---- Tier-0 blockage: no census means no fixture ------------------------
  ##
  ## smi_fixture() builds from a census. Probe which datasets have one, and
  ## report the gap rather than working around it. Identical blockage to
  ## dual-bounds': the two papers read the same tables, so the same census gap
  ## blocks both, and it is resolved once on the server, not twice here.

  census_keys <- contract$census_keys
  without_census <- setdiff(needed_keys, census_keys)

  if (length(without_census) > 0L) {
    cli::cli_alert_warning(
      "Tier-0 is BLOCKED: {length(without_census)} of {length(needed_keys)} required
       datasets have no {.fn smi_census} record."
    )
    cli::cli_alert_warning(
      "Without a census, {.fn smidata::smi_fixture} cannot build a structurally valid
       fixture for {.val {SMI_KEYS$covariates}}, {.val {SMI_KEYS$aap}}, or the
       cost_claims family."
    )
    cli::cli_alert_warning(
      "The unblocker is a server-side {.fn smi_census} run -- not anything in this
       repo. Do not hand-write a stand-in fixture."
    )
    cli::cli_alert_info(
      "What DOES run at Tier 0 today: {.file application/tests/}, which exercises the
       15-column design matrix and the estimator call at the real SHAPE on
       hand-built toy data. Shape, not distribution."
    )
  } else {
    cli::cli_alert_success(
      "Tier-0 is UNBLOCKED: all {length(needed_keys)} required datasets have a
       {.fn smi_census} record in snapshot {.val {contract$snapshot_id}}
       (promoted 2026-09-23). Same census gap as dual-bounds', resolved once for
       both."
    )
    cli::cli_alert_info(
      "{.fn smidata::smi_fixture} can now build structurally valid fixtures for
       {.val {SMI_KEYS$covariates}}, {.val {SMI_KEYS$aap}}, and the cost_claims
       family, in addition to what already ran at Tier 0
       ({.file application/tests/}'s 15-column design matrix on hand-built toy
       data). Distributional shape from real census, still not real data."
    )
  }
}

## ---- script-only section ---------------------------------------------------
## Literal guard, written inline: inside a function sys.nframe() is >= 1, so
## wrapping this test in a helper silently inverts it.
if (!interactive() && sys.nframe() == 0L) {
  cli::cli_alert_info("01_declare_requirements.R complete.")
}
