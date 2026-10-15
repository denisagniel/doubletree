## ============================================================================
## application/helpers/period_start.R
##
## The exposure start s (called T0 throughout this directory) and the helpers
## that locate it. PROJECT-WIDE RULE (PI, 2026-10-15):
##
##   s = the first day of the month ON OR AFTER INDEX_DT (INDEX_DT itself when
##   it is the 1st). The msr_saa row (MSR_YR = the period START) must start
##   EXACTLY at s: no nearest period, no gap tolerance, no ties. Patients with
##   no such row, or with MSR_DEN < 1, are ineligible and DROPPED.
##
## The period selection itself is smidata::smi_select_measure_row(anchor =
## "ceiling_month"), called directly from 02_population_and_eligibility.R. The
## local select_saa_after_index() ("start > INDEX_DT, target offset 1, gap cap")
## that used to live here is DELETED, not shimmed: it implemented the superseded
## 2026-10-15-morning rule. What remains is two pure date helpers.
## ============================================================================

#' First day of the month a measurement period starts in.
#'
#' @param msr_yr Character vector of `MSR_YR` labels, e.g. `"RY2019-06"`.
#'   Anything `smidata::smi_parse_msr_yr()` rejects aborts there; no value is
#'   coerced to `NA`.
#' @return A `Date` vector, the same length as `msr_yr`.
msr_yr_period_start <- function(msr_yr) {
  parsed <- smidata::smi_parse_msr_yr(as.character(msr_yr))
  as.Date(sprintf("%04d-%02d-01", as.integer(parsed$yyyy), as.integer(parsed$mm)))
}

#' The exposure start `s`: first day of the month on or after `INDEX_DT`.
#'
#' Mirrors the target `smidata::smi_select_measure_row(anchor = "ceiling_month")`
#' derives internally (this directory uses it to cross-check the selected
#' period and to build windows for patients without needing the msr table).
#' `INDEX_DT` on the 1st is its own `s`; any other day rolls to the 1st of the
#' next month, including 31 December -> 1 January of the next year.
#'
#' @param index_dt `Date` (or coercible) vector.
#' @return A `Date` vector, the same length, always the 1st of a month.
exposure_start <- function(index_dt) {
  d <- as.Date(index_dt)
  first <- lubridate::floor_date(d, "month")
  dplyr::if_else(
    lubridate::day(d) == 1L, first,
    lubridate::add_with_rollback(first, lubridate::period(months = 1L))
  )
}
