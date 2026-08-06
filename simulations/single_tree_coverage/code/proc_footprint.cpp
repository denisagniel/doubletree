// ============================================================
// proc_footprint.cpp
// Study: 2026-08-04_single-tree-coverage  (doubletree)
//
// Reads a process's macOS phys_footprint. Compiled at driver startup via
// Rcpp::sourceCpp(); used by the driver watchdog to cap worker memory.
//
// WHY NOT RSS: on macOS, RSS (what ps / ps::ps_memory_info report) EXCLUDES
// compressed pages. Under memory pressure the compressor swallows a leaking
// process's cold pages, so its RSS can FALL while its real memory charge keeps
// climbing -- an RSS watchdog never trips while the host thrashes to death.
// That is the observed 2026-08-04 failure (compressor at 7.5 GB at the time of
// the crash). phys_footprint = resident anonymous + compressed + IOKit mappings;
// it is what jetsam bands on and what Activity Monitor's "Memory" column shows.
//
// proc_pid_rusage() works on same-uid processes without root and costs
// microseconds, so it is safe to poll at 1 Hz.
//
// Not run directly. Sourced by pilot_complex_driver.R.
// ============================================================

#include <Rcpp.h>
#include <libproc.h>
#include <sys/resource.h>

//' Read phys_footprint (and swap precursors) for a pid
//'
//' @param pid Process id (same uid; no root needed).
//' @return Named numeric vector. ok = 0 means the pid is gone or unreadable --
//'   callers MUST branch on ok before trusting the other fields.
//'   footprint    : current phys_footprint, bytes (THE cap metric)
//'   lifetime_max : peak phys_footprint over the process lifetime, bytes
//'                  (post-hoc chunk sizing -- survives the poll-interval blind spot)
//'   pageins      : cumulative pageins; rising => child is already swapping,
//'                  i.e. the watchdog is late
// [[Rcpp::export]]
Rcpp::NumericVector proc_footprint(int pid) {
  rusage_info_current ri;
  if (proc_pid_rusage(pid, RUSAGE_INFO_CURRENT, (rusage_info_t *)&ri) != 0) {
    return Rcpp::NumericVector::create(Rcpp::_["ok"] = 0);
  }
  return Rcpp::NumericVector::create(
    Rcpp::_["ok"]           = 1,
    Rcpp::_["footprint"]    = (double) ri.ri_phys_footprint,
    Rcpp::_["lifetime_max"] = (double) ri.ri_lifetime_max_phys_footprint,
    Rcpp::_["pageins"]      = (double) ri.ri_pageins);
}
