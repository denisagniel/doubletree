## ============================================================================
## application/helpers/period_start.R
##
## msr_yr_period_start(): T0, the time zero for every window in this pipeline.
##
## T0 is the first day of the month an msr_* measurement period STARTS in
## (smidata's MSR_YR convention: the period runs from the first day of the
## MSR_YR month through the last day of the month eleven months later). See
## DOUBLETREE_WINDOWS in application/_config.R for why the windows are anchored
## here and not on INDEX_DT.
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
