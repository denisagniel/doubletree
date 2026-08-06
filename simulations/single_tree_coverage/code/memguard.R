# ============================================================
# memguard.R
# Study: 2026-08-04_single-tree-coverage  (doubletree)
#
# Memory + crash safety for the chunked-subprocess simulation drivers. Provides:
#   mg_init()        -- compile the footprint probe, open the event log
#   mg_host()        -- host-wide memory state (compressor / swap / free)
#   mg_wait_for_host() -- launch gate: refuse to start a worker onto a stressed host
#   mg_run_guarded() -- spawn ONE worker under a footprint watchdog
#   mg_chunk_done()  -- crash-safe "is this checkpoint really complete?" test
#   mg_event()       -- one flushed line per lifecycle event
#
# WHY A WATCHDOG AND NOT A KERNEL LIMIT (2026-08-05, verified on this host):
# Darwin gives an unprivileged process NO enforceable per-process memory cap.
#   * `ulimit -v` / RLIMIT_AS  -> setrlimit fails outright ("invalid argument").
#   * RLIMIT_DATA              -> setrlimit SUCCEEDS but constrains nothing, because
#                                 macOS malloc uses mach_vm_allocate/mmap, not brk.
#                                 Silent false safety -- deliberately NOT used.
#   * mem.maxVSize()/R_MAX_VSIZE -> bounds only R's VECTOR HEAP, enforced inside
#                                 allocVector. The leak here is native C++ (`new`/
#                                 std::vector) inside TreeFARMS, which never touches
#                                 R's allocator: the R heap stayed flat at 66 MB
#                                 while process memory grew to a SIGKILL. Useless here.
#   * launchd ResidentSetSize, jetsam APIs -> no-ops or entitlement-gated.
# So the driver polls its child and kills it. Note this is a HEURISTIC, not a
# guarantee: the only hard guarantee is running workers in a memory-capped Linux VM.
#
# Sourced by pilot_complex_driver.R. Not run directly.
# ============================================================

suppressPackageStartupMessages({ library(fs) })

GB <- 1024^3

# --- Budget for a 16 GB host --------------------------------------------------
# The per-child footprint cap is anchored on the 2026-08-04 incident (compressor
# reached ~7.7 GB before the host died): trip a child well below that.
#
# HOST-WIDE gating deliberately uses the KERNEL's own pressure signals, not absolute
# byte thresholds. Measured on this host 2026-08-05: idle-ish state was already
# compressor = 7.7 GB / free = 0.06 GB / 18 GB logical committed, so byte thresholds
# like "compressor < 2 GB" never clear and would abort every chunk forever. The
# kernel signals are self-normalising:
#   kern.memorystatus_vm_pressure_level : 1 = normal, 2 = warn, 4 = critical
#   kern.memorystatus_level             : 0-100, % memory available (higher = better)
# CRITICAL is the kernel saying jetsam is imminent -- that is precisely when we must
# kill OUR child first, so the OS never gets to pick our DRIVER as the victim (which
# is what happened on 2026-08-04).
MG <- list(
  hard_cap        = 5.0 * GB,   # child phys_footprint -> SIGKILL the child
  lookahead_s     = 3,          # predictive trip horizon (footprint can jump fast)
  poll_s          = 1,          # child footprint poll (proc_pid_rusage ~ microseconds)
  host_every      = 5,          # poll host state every Nth child poll (shells out)
  launch_avail_pct = 20,        # launch gate: need >= this % memory available
  launch_max_level = 4,         # launch gate: refuse to launch at CRITICAL
  trip_crit_hits  = 2,          # kill child after N consecutive CRITICAL host polls
  trip_avail_pct  = 8,          # ...or if available memory falls below this %
  bail_aborts     = 3,          # N consecutive host aborts -> give up; host > run
  max_chunk_s     = 45 * 60,    # hung native call guard
  kill_grace_s    = 10,         # SIGKILL teardown of a multi-GB task is not instant
  cooldown_s      = 20          # never launch back-to-back after a trip
)

.mg <- new.env(parent = emptyenv())
.mg$event_log <- NULL
.mg$degraded  <- FALSE          # TRUE => footprint probe unavailable, using RSS at half cap
.mg$hard_cap  <- MG$hard_cap    # EFFECTIVE cap; halved when degraded (set in mg_init)

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Initialise memguard: compile the footprint probe and open the event log
#'
#' Compiles proc_footprint.cpp (phys_footprint via proc_pid_rusage). If that fails
#' we fall back to `ps` RSS with the cap HALVED and flag the run DEGRADED, because
#' RSS under-reports on macOS (it excludes compressed pages) and would otherwise
#' trip far too late.
#'
#' @param study_dir Study root.
#' @param event_log Path for the machine-readable event log.
#' @return Invisibly TRUE if the native probe is active, FALSE if degraded.
mg_init <- function(study_dir, event_log) {
  .mg$event_log <- event_log
  fs::dir_create(fs::path_dir(event_log))
  ok <- tryCatch({
    Rcpp::sourceCpp(fs::path(study_dir, "code", "proc_footprint.cpp"))
    TRUE
  }, error = function(e) { .mg$err <- conditionMessage(e); FALSE })
  .mg$degraded <- !ok
  .mg$hard_cap <- if (ok) MG$hard_cap else MG$hard_cap / 2
  if (!ok) {
    mg_event("DEGRADED_METRIC", reason = "sourceCpp_failed",
             detail = gsub("[\t\n]", " ", .mg$err %||% "NA"),
             new_cap_gb = round(.mg$hard_cap / GB, 2))
  }
  invisible(ok)
}

#' Current memory charge of a process, in bytes
#'
#' Uses phys_footprint (resident + COMPRESSED + IOKit) when available. RSS is the
#' degraded fallback only -- see mg_init().
#'
#' @param pid Process id.
#' @return List(ok, bytes, lifetime_max, pageins). ok = FALSE => pid gone/unreadable.
mg_mem <- function(pid) {
  if (!.mg$degraded) {
    fi <- proc_footprint(as.integer(pid))
    if (fi[["ok"]] != 1) return(list(ok = FALSE))
    return(list(ok = TRUE, bytes = fi[["footprint"]],
                lifetime_max = fi[["lifetime_max"]], pageins = fi[["pageins"]]))
  }
  out <- tryCatch({
    mi <- ps::ps_memory_info(ps::ps_handle(as.integer(pid)))
    list(ok = TRUE, bytes = mi[["rss"]], lifetime_max = NA_real_, pageins = mi[["pageins"]])
  }, error = function(e) list(ok = FALSE))
  out
}

#' Host-wide memory state
#'
#' Primary signals are the kernel's own: pressure level and % available. Byte
#' counters are recorded for the event log / post-mortem but are NOT gated on --
#' see the MG comment for why absolute byte thresholds are unusable on this host.
#'
#' Page size is PARSED from the vm_stat header, never assumed: it is 16384 on
#' Apple Silicon and 4096 on Intel, so a hardcoded 4096 under-reports by 4x.
#'
#' @return List: ok, level (1 normal / 2 warn / 4 critical), avail_pct (0-100), and
#'   bytes for compressor / compressor_stored / free / purgeable / swap_used.
mg_host <- function() {
  sysctl_num <- function(key) {
    v <- tryCatch(system2("sysctl", c("-n", key), stdout = TRUE, stderr = FALSE),
                  error = function(e) NA_character_)
    suppressWarnings(as.numeric(v[1]))
  }
  level <- sysctl_num("kern.memorystatus_vm_pressure_level")
  avail <- sysctl_num("kern.memorystatus_level")

  vs <- tryCatch(system2("vm_stat", stdout = TRUE, stderr = FALSE),
                 error = function(e) character())
  psz <- if (length(vs)) {
    suppressWarnings(as.numeric(sub(".*page size of ([0-9]+) bytes.*", "\\1", vs[1])))
  } else NA_real_
  pages <- function(key) {
    if (!length(vs) || is.na(psz)) return(NA_real_)
    ln <- grep(key, vs, fixed = TRUE, value = TRUE)[1]
    if (is.na(ln)) return(NA_real_)
    suppressWarnings(as.numeric(gsub("[^0-9]", "", sub("^[^:]*:", "", ln)))) * psz
  }
  su <- tryCatch(system2("sysctl", c("-n", "vm.swapusage"), stdout = TRUE, stderr = FALSE),
                 error = function(e) NA_character_)
  swap <- suppressWarnings(as.numeric(sub(".*used = ([0-9.]+)M.*", "\\1", su[1]))) * 1024^2

  list(ok = !is.na(level) || !is.na(avail),
       level = level, avail_pct = avail, page = psz,
       # NB: "occupied by compressor" (physical cost) and "stored in compressor"
       # (logical bytes held) are DIFFERENT vm_stat lines, differing by the ~2.3x
       # compression ratio. Grepping the wrong one badly misstates the cost.
       compressor        = pages("Pages occupied by compressor"),
       compressor_stored = pages("Pages stored in compressor"),
       free      = pages("Pages free"),
       purgeable = pages("Pages purgeable"),
       swap_used = if (length(swap) == 1 && !is.na(swap)) swap else NA_real_)
}

#' Append one flushed event line (open-write-close, so a kill cannot lose it)
#'
#' Buffered logging is why the 2026-08-04 driver death left NO final log line; every
#' event here is flushed immediately so a post-mortem always has the last state.
#'
#' @param event Event name (see the taxonomy in the driver header).
#' @param ... Named scalar fields.
#' @return The formatted line, invisibly.
mg_event <- function(event, ...) {
  f <- list(...)
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\t", event,
                 if (length(f)) paste0("\t", paste(names(f), unlist(f), sep = "=", collapse = " ")) else "")
  if (!is.null(.mg$event_log)) cat(line, "\n", sep = "", file = .mg$event_log, append = TRUE)
  invisible(line)
}

#' Launch gate: block until the host can afford another worker
#'
#' Gates on the kernel's pressure level and % available, NOT on absolute bytes.
#' Deliberately permissive (it blocks only CRITICAL or genuinely low availability),
#' because this host's normal state is already WARN with a multi-GB compressor -- a
#' stricter gate would never clear and would stall the run forever.
#'
#' @param max_tries,sleep_s Retry budget.
#' @return TRUE if clear to launch, FALSE if the host stayed stressed (caller aborts).
mg_wait_for_host <- function(max_tries = 20, sleep_s = 30) {
  for (k in seq_len(max_tries)) {
    h <- mg_host()
    if (!isTRUE(h$ok)) return(TRUE)          # cannot read host state -> do not block the run
    clear <- (is.na(h$level)     || h$level     <  MG$launch_max_level) &&
             (is.na(h$avail_pct) || h$avail_pct >= MG$launch_avail_pct)
    if (clear) return(TRUE)
    mg_event("HOST_PRESSURE_WAIT", try = k, level = h$level, avail_pct = h$avail_pct,
             comp_gb = round(h$compressor / GB, 2),
             swap_gb = round(h$swap_used / GB, 2))
    Sys.sleep(sleep_s)
  }
  FALSE
}

#' Has this chunk been permanently blocked by the memory cap?
#'
#' TRUE => already attempted and exceeded the cap; retrying just burns time. Clear the
#' .oom markers (or set GS_RETRY_OOM=1) to force a retry, e.g. after raising the cap
#' or moving to a bigger host.
#'
#' @param out Checkpoint path.
#' @return Logical scalar.
mg_chunk_blocked <- function(out) fs::file_exists(paste0(out, ".oom"))

#' Is this checkpoint genuinely complete?
#'
#' Sentinel first (worker writes <out>.done via atomic rename AFTER the .rds is
#' closed). Legacy fallback: a pre-sentinel .rds counts as done only if it actually
#' DESERIALISES -- which is itself the integrity test, since a write_rds truncated by
#' a SIGKILL fails to read. A bare file.exists() would treat a truncated checkpoint
#' as complete and silently poison the merge.
#'
#' @param out Checkpoint path.
#' @return Logical scalar.
mg_chunk_done <- function(out) {
  if (fs::file_exists(paste0(out, ".done"))) return(TRUE)
  if (!fs::file_exists(out)) return(FALSE)
  ok <- tryCatch({ readr::read_rds(out); TRUE }, error = function(e) FALSE)
  if (!ok) mg_event("CHECKPOINT_CORRUPT", file = fs::path_file(out))
  ok
}

#' Run ONE worker subprocess under a footprint watchdog
#'
#' Output goes to FILES, never a pipe: an undrained processx pipe deadlocks the child
#' at 64 KB, and piping into R would grow the driver's own footprint (the driver must
#' stay small so it keeps getting scheduled under pressure).
#'
#' Threads are pinned to 1 -- per-thread malloc arenas across 10 cores inflate peak
#' footprint substantially, and this workload is memory-bound, not CPU-bound.
#'
#' supervise = TRUE kills the child if the DRIVER dies, so a driver death can never
#' leave a multi-GB orphan behind (exactly what 2026-08-04 risked).
#'
#' @param worker,args Worker script and its positional args.
#' @param out_log,err_log Child stdout/stderr paths.
#' @param sentinel Completion sentinel the worker writes on success.
#' @param label Chunk label for events.
#' @return List(status, exit, peak, lifetime_max, elapsed_s). status is one of
#'   CHUNK_OK / ESTIMATOR_ERROR / NO_SENTINEL_EXIT_0 / OOM_CAP_CHILD /
#'   OOM_CAP_CHILD_PREDICT / HOST_PRESSURE_ABORT / OOM_KILLED_BY_OS / TIMEOUT.
mg_run_guarded <- function(worker, args, out_log, err_log, sentinel, label) {
  if (fs::file_exists(sentinel)) fs::file_delete(sentinel)   # never trust a stale sentinel

  rscript <- c("Rscript", worker, args)
  nice <- unname(Sys.which("nice"))
  cmd  <- if (nzchar(nice)) list(nice, c("-n", "10", rscript)) else list("Rscript", c(worker, args))

  proc <- processx::process$new(
    cmd[[1]], cmd[[2]],
    stdout = out_log, stderr = err_log,
    env = c("current", OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
            VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1"),
    cleanup_tree = TRUE, supervise = TRUE
  )
  pid <- proc$get_pid()
  mg_event("CHUNK_START", chunk = label, pid = pid)

  t_start <- Sys.time()
  peak <- 0; lmax <- NA_real_; prev_b <- NA_real_; prev_t <- NA_real_
  trip <- NA_character_; iter <- 0L; crit_hits <- 0L

  repeat {
    proc$wait(timeout = MG$poll_s * 1000)      # wakes immediately on child exit
    if (!proc$is_alive()) break
    iter <- iter + 1L
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

    m <- mg_mem(pid)
    if (isTRUE(m$ok)) {
      peak <- max(peak, m$bytes)
      if (!is.na(m$lifetime_max)) lmax <- m$lifetime_max
      now <- as.numeric(Sys.time())
      rate <- if (!is.na(prev_b) && now > prev_t) (m$bytes - prev_b) / (now - prev_t) else 0
      prev_b <- m$bytes; prev_t <- now
      if (m$bytes > .mg$hard_cap) {
        trip <- "OOM_CAP_CHILD"
      } else if (rate > 0 && m$bytes + MG$lookahead_s * rate > .mg$hard_cap) {
        # Footprint can jump GBs inside one native call; trip on trajectory, not just level.
        trip <- "OOM_CAP_CHILD_PREDICT"
      }
      if (!is.na(trip)) {
        mg_event(trip, chunk = label, pid = pid,
                 footprint_gb = round(m$bytes / GB, 2),
                 rate_mb_s = round(rate / 1024^2, 1),
                 pageins = m$pageins, elapsed_s = round(elapsed))
      }
    }

    if (is.na(trip) && iter %% MG$host_every == 0) {
      h <- mg_host()
      if (isTRUE(h$ok)) {
        # CRITICAL means jetsam is imminent. Require consecutive hits so a transient
        # spike does not needlessly discard a chunk, but act BEFORE the OS picks a
        # victim -- on 2026-08-04 it picked our driver.
        crit <- (!is.na(h$level) && h$level >= 4) ||
                (!is.na(h$avail_pct) && h$avail_pct < MG$trip_avail_pct)
        crit_hits <- if (crit) crit_hits + 1L else 0L
        if (crit_hits >= MG$trip_crit_hits) {
          trip <- "HOST_PRESSURE_ABORT"
          mg_event(trip, chunk = label, level = h$level, avail_pct = h$avail_pct,
                   hits = crit_hits, comp_gb = round(h$compressor / GB, 2),
                   swap_gb = round(h$swap_used / GB, 2))
        }
      }
    }

    if (is.na(trip) && elapsed > MG$max_chunk_s) {
      trip <- "TIMEOUT"
      mg_event(trip, chunk = label, elapsed_s = round(elapsed))
    }

    if (!is.na(trip)) {
      proc$kill_tree()
      proc$wait(timeout = MG$kill_grace_s * 1000)
      if (proc$is_alive()) {
        mg_event("KILL_SLOW", chunk = label, pid = pid)
        proc$wait(timeout = 30000)
      }
      break
    }
  }

  exit <- tryCatch(proc$get_exit_status(), error = function(e) NA_integer_)
  elapsed_s <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

  status <- if (!is.na(trip)) {
    trip
  } else if (fs::file_exists(sentinel) && !is.na(exit) && exit == 0) {
    "CHUNK_OK"
  } else if (!is.na(exit) && exit == 0) {
    # Clean exit, no sentinel: worker logic bug, NOT an OOM. Must not look like success.
    "NO_SENTINEL_EXIT_0"
  } else if (!is.na(exit) && (exit == 137 || exit < 0)) {
    # Signal death with no watchdog trip => we were too slow. Tuning alarm.
    "OOM_KILLED_BY_OS"
  } else {
    "ESTIMATOR_ERROR"
  }

  # A rep that blew the cap will blow it again on every resume, so record a PERMANENT
  # marker or the run can never finish. It stays VISIBLE, not silent: the driver emits
  # an explicit status = "skipped_oom_local" row for each blocked rep, because these
  # are the reps with the LARGEST Rashomon sets and dropping them is not
  # missing-at-random -- it would bias the coverage / non-recovery endpoints.
  # HOST_PRESSURE_ABORT deliberately gets NO marker: that was the host's fault, not
  # this rep's, so it must be retried.
  if (status %in% c("OOM_CAP_CHILD", "OOM_CAP_CHILD_PREDICT", "OOM_KILLED_BY_OS")) {
    cat(sprintf("status=%s peak_gb=%.3f cap_gb=%.2f when=%s\n", status, peak / GB,
                .mg$hard_cap / GB, format(Sys.time())),
        file = sub("\\.done$", ".oom", sentinel))
  }

  mg_event(status, chunk = label, exit = exit, peak_gb = round(peak / GB, 3),
           lifetime_max_gb = if (is.na(lmax)) NA else round(lmax / GB, 3),
           elapsed_s = round(elapsed_s, 1))

  if (!is.na(trip)) Sys.sleep(MG$cooldown_s)   # let the OS reclaim before the next launch

  list(status = status, exit = exit, peak = peak, lifetime_max = lmax,
       elapsed_s = elapsed_s, fatal = identical(trip, "HOST_PRESSURE_ABORT"))
}
