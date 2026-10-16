## ============================================================================
## application/helpers/disclosure.R
##
## Small-cell plumbing shared by the decision report (08) and the Tier-1 exports
## (90). The disclosure LOGIC lives in smidata (>= the identity-aware release,
## commit d4d6e76): smi_suppress() masks cells; smi_report() resolves and audits
## the declared identities between cells (smi_identity(), smi_gap()); chains are
## coarsened by smi_coarsen_flow(). This project's own copies of
## hide_small_gaps(), coarsen_flow(), mask_partition_cells(), mask_binary_split()
## and tier1_disclosable_*() were DELETED in favour of those.
##
## What stays here is what smidata does not provide:
##   * mask_counts_primary()    primary suppression ONLY, for counts that are not
##                              a partition (smi_suppress() always adds a pooled
##                              complement, which for independent counts protects
##                              nothing and hides an unrelated cell). Cross-cell
##                              relations are then declared as identities and
##                              repaired by smi_report().
##   * mask_counts_partition()  smi_suppress() + the table tags below.
##   * tag_report_table()       thin wrapper over smidata::smi_report_tag() (the
##                              exported tagger smi_report() requires on every
##                              table: counts, derived columns, labels).
##   * tier1_saa_report()       the Tier-1 tables that share counts, assembled into
##                              ONE smi_report so a single audit covers them.
##
## WHAT THE AUDIT GUARANTEES (and nothing more): no masked value is determined by
## the published cells through the DECLARED identities, no declared sum of masked
## parts and no declared difference of published counts is a count in 1..10.
## Identities nobody declared, dollar dominance and a masked cell narrowed to a
## range inside 1..10 are not covered. See smidata's ?smi_disclosure_audit.
## ============================================================================

config_report_min_cell <- 11L

#' Declare a table's disclosure roles (wrapper over `smidata::smi_report_tag()`).
#'
#' @param tb A data frame (already masked).
#' @param n_cols Count columns (suppressed by smi_export()).
#' @param derived_cols Columns derived from counts (masked with them).
#' @param label_cols Numeric columns that are labels, not counts (written as text).
#' @param identities Identities local to this table (rarely needed; project
#'   identities are passed to `smi_report(identities = )`).
#' @param derived_from Optional `list(derived_col = c(count cols it depends on))`:
#'   a derived column is then masked only when one of ITS sources is masked.
tag_report_table <- function(tb, n_cols, derived_cols = character(),
                             label_cols = character(), identities = list(),
                             derived_from = NULL) {
  smidata::smi_report_tag(tb, n_cols = n_cols, derived_cols = derived_cols,
                          label_cols = label_cols, identities = identities,
                          derived_from = derived_from)
}

#' Primary suppression only: 1..(min_cell - 1) become NA, zero stays zero.
#'
#' A row's `derived_cols` are masked with any masked count in it (a rate or mean
#' returns its numerator otherwise); with `derived_from`, only with a masked
#' count it depends on.
mask_counts_primary <- function(tb, n_cols, derived_cols = character(),
                                label_cols = character(),
                                min_cell = config_report_min_cell,
                                derived_from = NULL) {
  hit_by <- list()
  for (cl in n_cols) {
    v <- tb[[cl]]
    sm <- !is.na(v) & v > 0 & v < min_cell
    tb[[cl]][sm] <- NA
    hit_by[[cl]] <- sm
  }
  any_hit <- Reduce(`|`, hit_by, rep(FALSE, nrow(tb)))
  for (cl in derived_cols) {
    src <- derived_from[[cl]]
    hit <- if (is.null(src)) any_hit else {
      Reduce(`|`, hit_by[intersect(src, n_cols)], rep(FALSE, nrow(tb)))
    }
    tb[[cl]][hit] <- NA
  }
  tag_report_table(tb, n_cols, derived_cols, label_cols, derived_from = derived_from)
}

#' smi_suppress() (primary + pooled complement) for a column that partitions a
#' total, then tagged. Identities linking it to other tables are repaired later.
mask_counts_partition <- function(tb, n_cols, derived_cols = NULL,
                                  label_cols = character(), group = NULL,
                                  min_cell = config_report_min_cell) {
  if (nrow(tb) > 0L) {
    tb <- smidata::smi_suppress(tb, n_cols = n_cols, derived_cols = derived_cols,
                                threshold = min_cell, group = group)
  }
  tag_report_table(tb, n_cols, derived_cols %||% character(), label_cols)
}

`%||%` <- function(x, y) if (is.null(x)) y else x


## ---- Tier-1 (90_checks_tier1.R) -------------------------------------------

#' The exact-start sizes as a long table (`quantity`, `n`), so that the chain and
#' the arm partition can be declared as identities over named cells.
tier1_sizes_table <- function(by_band) {
  q <- c("n_with_saa_rows", "n_matched", "n_elig", "n_treated", "n_control",
         "n_achieved_na")
  tibble::tibble(quantity = q, n = as.integer(unlist(by_band[q], use.names = FALSE)))
}

#' One audited report for the Tier-1 tables that publish related counts.
#'
#' CHECK 6 (exact-start sizes, MSR_YR offsets, gap by offset) and CHECK 3
#' (`_YN` prevalence by arm) publish counts that are linked: an arm size is the
#' sum of a column's yes / no / NA cells, `n_elig` is the sum of the arm sizes,
#' the chain with-rows >= matched >= eligible publishes every drop, and the
#' nearest-period match counts bracket the exact-start ones. Exporting them from
#' separate reports would let a cell masked in one be read off the other, so they
#' are built, resolved and audited together.
#'
#' @param by_band One-row tibble of the six exact-start counts (raw).
#' @param rel_bands `band`, `n_rows` (raw).
#' @param by_offset `target_offset_months`, `n_matched`, `median_gap`, `p90_gap`,
#'   `share_abs_gap_le_1` (raw).
#' @param prevalence `column` and the six `n_{yes,no,na}_{control,treated}`
#'   counts with the two prevalence columns (raw), or `NULL`.
#' @return An audited `smi_report`.
tier1_saa_report <- function(by_band, rel_bands, by_offset, prevalence = NULL,
                             min_cell = config_report_min_cell) {
  sizes <- mask_counts_primary(tier1_sizes_table(by_band), "n", min_cell = min_cell)
  bands <- mask_counts_partition(rel_bands, "n_rows", min_cell = min_cell)
  offs <- mask_counts_primary(
    by_offset, "n_matched",
    derived_cols = c("median_gap", "p90_gap", "share_abs_gap_le_1"),
    min_cell = min_cell
  )

  sz <- function(q) smidata::smi_cell("exact_start_sizes", "n", quantity = q)
  ids <- list(
    smidata::smi_identity(sz("n_elig"),
                          list(sz("n_treated"), sz("n_control"), sz("n_achieved_na")),
                          label = "eligible = treated + control + NA"),
    smidata::smi_gap(sz("n_with_saa_rows"), sz("n_matched"), label = "rows -> matched"),
    smidata::smi_gap(sz("n_matched"), sz("n_elig"), label = "matched -> eligible"),
    smidata::smi_gap(sz("n_with_saa_rows"), sz("n_elig"), label = "rows -> eligible")
  )
  ## The nearest-period match counts (any gap) bracket the exact-start ones, so
  ## their difference is a count nobody should be able to read. smi_gap() needs
  ## big >= small on the published values, so the direction is taken from the RAW
  ## counts (equal counts differ by 0, which is not a small count: nothing to
  ## declare).
  for (k in seq_len(nrow(by_offset))) {
    off <- by_offset$target_offset_months[k]
    raw_off <- by_offset$n_matched[k]
    oc <- smidata::smi_cell("gap_by_offset", "n_matched", target_offset_months = off)
    for (q in c("n_with_saa_rows", "n_matched")) {
      raw_q <- as.numeric(by_band[[q]])
      if (is.na(raw_off) || is.na(raw_q) || raw_off == raw_q) next
      ids <- c(ids, list(
        if (raw_off > raw_q) {
          smidata::smi_gap(oc, sz(q), label = paste0("offset ", off, " vs ", q))
        } else {
          smidata::smi_gap(sz(q), oc, label = paste0(q, " vs offset ", off))
        }
      ))
    }
  }

  tables <- list(exact_start_sizes = sizes, msr_yr_offsets = bands,
                 gap_by_offset = offs)
  if (!is.null(prevalence) && nrow(prevalence) > 0L) {
    cnt <- grep("^n_", names(prevalence), value = TRUE)
    der <- c("prevalence_control", "prevalence_treated")
    ## prevalence_<arm> = yes / (yes + no): masked only with ITS arm's yes / no.
    dfrom <- list(prevalence_control = c("n_yes_control", "n_no_control"),
                  prevalence_treated = c("n_yes_treated", "n_no_treated"))
    tables$yn_prevalence <- mask_counts_primary(prevalence, cnt, der, min_cell = min_cell,
                                                derived_from = dfrom)
    for (col in prevalence$column) {
      for (arm in c("control", "treated")) {
        pv <- function(k) smidata::smi_cell("yn_prevalence", paste0(k, "_", arm), column = col)
        ids <- c(ids, list(smidata::smi_identity(
          sz(paste0("n_", arm)), list(pv("n_yes"), pv("n_no"), pv("n_na")),
          label = paste0("arm size = yes + no + NA, ", col, ", ", arm)
        )))
      }
    }
  }
  do.call(smidata::smi_report, c(
    list("doubletree Tier-1: msr_saa exact-start sizes and design checks",
         levers = list(anchor = config_period_anchor,
                       achievement_threshold = config_saa_achievement_threshold)),
    tables, list(identities = ids, min_cell = min_cell)
  ))
}
