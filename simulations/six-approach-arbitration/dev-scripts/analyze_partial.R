# =============================================================================
# analyze_partial.R -- completeness-first analysis of a (possibly partial) run
# =============================================================================
# The 20260722 run OOM-killed the expensive Rashomon cells (complex/continuous,
# escalate-complex) at a 4 GB cap, so combine.R wrote a PARTIAL result. Averaging
# blindly over half-empty stress cells would misreport them. This script therefore
# leads with per-cell COMPLETENESS and only summarizes cells against their realized
# rep count, flagging any cell below a reporting threshold.
#
# Base-R only (runs under module R). Reproduces run_one.R conventions: truth=0.15,
# covered/covered_crossfit already computed per unit.
#
#   Rscript dev-scripts/analyze_partial.R [path/to/run.rds]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
rds  <- if (length(args) >= 1) args[1] else {
  cand <- c("20260722-114141_fffe596.rds",
            "results/20260722-114141_fffe596.rds")
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) stop("No run rds found; pass a path.", call. = FALSE)
  hit[1]
}
cat(sprintf("Reading: %s\n", rds))
r <- readRDS(rds)

# Study identity: rebuild the grid to know which cells are EXPECTED and which are
# INFEASIBLE (legitimately empty, not failures).
study_root <- {
  fa <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(fa) == 1) dirname(dirname(normalizePath(sub("^--file=", "", fa)))) else "."
}
source(file.path(study_root, "config", "grid.R"))
TOTAL_REPS <- get0("TOTAL_REPS", ifnotfound = 1000L)

# --- Expected reps per cell (0 if infeasible) --------------------------------
cells <- unique(GRID[, c("method", "n", "dgp", "escalate")])
cells$feasible <- is_feasible(cells)
cells$expected <- ifelse(cells$feasible, TOTAL_REPS, 0L)

# --- Realized reps per cell ---------------------------------------------------
key <- function(d) paste(d$method, d$n, d$dgp, d$escalate, sep = "\r")
got <- as.data.frame(table(key(r)), stringsAsFactors = FALSE)
names(got) <- c("key", "got")
cells$key <- key(cells)
cells <- merge(cells, got, by = "key", all.x = TRUE)
cells$got[is.na(cells$got)] <- 0L
cells$pct <- ifelse(cells$expected > 0, round(100 * cells$got / cells$expected, 1), NA)

ord <- order(factor(cells$method, levels = unique(GRID$method)),
             cells$escalate, cells$dgp, cells$n)
cells <- cells[ord, ]

cat("\n===================== COMPLETENESS BY CELL =====================\n")
show <- cells[, c("method", "n", "dgp", "escalate", "expected", "got", "pct")]
print(show, row.names = FALSE)

cat("\n--- Rollup ---\n")
full_cells <- sum(cells$expected > 0 & cells$got >= cells$expected)
part_cells <- sum(cells$expected > 0 & cells$got > 0 & cells$got < cells$expected)
empty_feas <- sum(cells$expected > 0 & cells$got == 0)
infeas     <- sum(cells$expected == 0)
cat(sprintf("  feasible cells complete (>=%d reps): %d\n", TOTAL_REPS, full_cells))
cat(sprintf("  feasible cells PARTIAL (1..%d reps) : %d\n", TOTAL_REPS - 1L, part_cells))
cat(sprintf("  feasible cells EMPTY (0 reps, OOM)  : %d\n", empty_feas))
cat(sprintf("  infeasible cells (expected empty)   : %d\n", infeas))

if (part_cells > 0) {
  cat("\n  PARTIAL cells (report with caution):\n")
  pc <- cells[cells$expected > 0 & cells$got > 0 & cells$got < cells$expected,
              c("method", "n", "dgp", "escalate", "got", "pct")]
  print(pc, row.names = FALSE)
}
if (empty_feas > 0) {
  cat("\n  EMPTY feasible cells (OOM-killed, NOT run -> exclude from all summaries):\n")
  ec <- cells[cells$expected > 0 & cells$got == 0,
              c("method", "n", "dgp", "escalate")]
  print(ec, row.names = FALSE)
}

# --- Summaries, only where reps are adequate ---------------------------------
# Threshold: report a cell only if it has >= MIN_REPS. Coverage MC error at 1000
# reps ~ +/-1.4pp; at 200 ~ +/-3pp. Use 200 as a soft floor, flag 200..999.
MIN_REPS <- 200L
r$cellkey <- key(r)
adequate  <- cells$key[cells$got >= MIN_REPS & cells$expected > 0]
ra <- r[r$cellkey %in% adequate, ]

cov_se <- function(x) { p <- mean(x, na.rm = TRUE); sqrt(p * (1 - p) / sum(is.finite(x))) }
summ_cell <- function(d) {
  data.frame(
    method   = d$method[1], n = d$n[1], dgp = d$dgp[1], escalate = d$escalate[1],
    reps     = sum(is.finite(d$estimate)),
    bias     = round(mean(d$error, na.rm = TRUE), 4),
    emp_sd   = round(sd(d$estimate, na.rm = TRUE), 4),
    mean_se  = round(mean(d$std_error, na.rm = TRUE), 4),
    se_ratio = round(mean(d$std_error, na.rm = TRUE) / sd(d$estimate, na.rm = TRUE), 3),
    cover    = round(mean(d$covered, na.rm = TRUE), 3),
    cover_mc = round(cov_se(d$covered), 3),
    cover_cf = round(mean(d$covered_crossfit, na.rm = TRUE), 3),
    conv     = round(mean(d$converged, na.rm = TRUE), 3),
    int_ne   = round(mean(d$intersection_nonempty, na.rm = TRUE), 3),
    row.names = NULL, stringsAsFactors = FALSE)
}
if (nrow(ra) > 0) {
  parts <- split(ra, ra$cellkey)
  tab <- do.call(rbind, lapply(parts, summ_cell))
  tab <- tab[order(factor(tab$method, levels = unique(GRID$method)),
                   tab$escalate, tab$dgp, tab$n), ]
  cat(sprintf("\n============ PER-CELL SUMMARIES (reps >= %d) ============\n", MIN_REPS))
  cat("bias=mean(est-0.15); se_ratio=mean_se/emp_sd (want ~1); cover=CI coverage (target 0.95);\n")
  cat("cover_cf=cross-fit twin coverage; int_ne=Rashomon intersection nonempty rate.\n\n")
  print(tab, row.names = FALSE)
} else {
  cat("\n(no cells meet the reporting threshold)\n")
}

cat("\nDONE. Cells below threshold or empty are excluded above; see completeness table.\n")
