## ============================================================================
## application/tests/test-disclosure.R
##
## helpers/disclosure.R and the Tier-1 report. The disclosure logic itself
## (identities, audit, repair, chain coarsening) is smidata's and is tested there;
## this file tests what THIS project adds: the primary-only masker, the table
## tags smi_report_export() reads, and tier1_saa_report(), which declares the
## Tier-1 tables' identities. The property under test is that NO masked value is
## determined by the published cells (smi_report_audit) AND that the fixture's
## small value is not a difference of two published counts, checked by hand.
## Every expected value is worked out in a comment.
## ============================================================================

skip_if_not_installed("smidata")
skip_if_not("smi_report_audit" %in% getNamespaceExports("smidata"))

audit_rows_d <- function(rep) nrow(smidata::smi_report_audit(rep, abort = FALSE))

test_that("mask_counts_primary() masks 1..10, keeps zero, and masks derived columns of the row", {
  tb <- tibble::tibble(k = c("a", "b", "c", "d"), n = c(5L, 0L, 11L, 10L),
                       rate = c(.1, 0, .2, .3))
  m <- mask_counts_primary(tb, "n", "rate", min_cell = 11L)
  ## 5 and 10 are small; 0 is a degeneracy flag, 11 is publishable.
  expect_equal(m$n, c(NA, 0L, 11L, NA))
  expect_equal(is.na(m$rate), c(TRUE, FALSE, FALSE, TRUE))
  ## The tags smi_report_export() reads.
  expect_equal(attr(m, "smi_n_cols"), "n")
  expect_equal(attr(m, "smi_derived_cols"), "rate")
})

test_that("blank_cell() blanks one row's cell and its derived columns, nested tables included", {
  tabs <- list(
    a = mask_counts_primary(tibble::tibble(k = c("x", "y"), n = c(50L, 60L), r = c(.5, .6)),
                            "n", "r"),
    grp = list(inner = mask_counts_primary(tibble::tibble(k = c("x", "y"), n = c(70L, 80L)), "n"))
  )
  out <- blank_cell(tabs, smidata::smi_cell("a", "n", k = "y"))
  expect_equal(out$a$n, c(50L, NA))
  expect_equal(out$a$r, c(.5, NA))
  out <- blank_cell(out, smidata::smi_cell("grp.inner", "n", k = "x"))
  expect_equal(out$grp$inner$n, c(NA, 80L))
})

mk_sizes <- function(rows, matched, elig, tr, ct, na) {
  tibble::tibble(n_with_saa_rows = rows, n_matched = matched, n_elig = elig,
                 n_treated = tr, n_control = ct, n_achieved_na = na)
}
bands_ok <- tibble::tibble(band = c("-12", "0", "12"), n_rows = c(30L, 400L, 70L))
offsets_ok <- tibble::tibble(target_offset_months = c("0", "12"), n_matched = c(480L, 470L),
                             median_gap = c(0, 12), p90_gap = c(1, 13),
                             share_abs_gap_le_1 = c(.9, .8))
sz_of <- function(rep) {
  t <- rep$tables$exact_start_sizes
  stats::setNames(t$n, t$quantity)
}

test_that("tier1_saa_report(): nothing small -> everything published, no repair", {
  ## 500 rows >= 400 matched >= 300 eligible = 100 + 150 + 50.
  r <- tier1_saa_report(mk_sizes(500L, 400L, 300L, 100L, 150L, 50L), bands_ok, offsets_ok)
  expect_equal(r$resolved, 0L)
  expect_equal(unname(sz_of(r)), c(500L, 400L, 300L, 100L, 150L, 50L))
  expect_equal(audit_rows_d(r), 0L)
})

test_that("tier1_saa_report(): a matched -> eligible step of 5 hides the larger count", {
  ## 400 -> 395 is 5 patients. Eligible = 200 + 150 + 45 = 395. Both 400 and 395 are
  ## published only if 400 - 395 may be read, so the declared gap masks the larger
  ## (matched); 500 - 395 = 105 and 500 - 400 are fine.
  r <- tier1_saa_report(mk_sizes(500L, 400L, 395L, 200L, 150L, 45L), bands_ok, offsets_ok)
  s <- sz_of(r)
  expect_true(is.na(s[["n_matched"]]))
  expect_equal(s[["n_elig"]], 395L)
  expect_equal(s[["n_with_saa_rows"]], 500L)
  expect_equal(audit_rows_d(r), 0L)
  ## by hand: no published pair differs by 1..10
  pub <- s[!is.na(s)]
  d <- outer(pub, pub, "-")
  expect_false(any(d >= 1 & d <= 10))
})

test_that("tier1_saa_report(): two consecutive small steps leave no small difference", {
  ## 400 -> 398 -> 395: steps of 2 and 3 (rows 400, matched 398, eligible 395 =
  ## 200 + 150 + 45). Each declared gap masks its larger count: rows, then matched.
  ## Only eligible (395) remains of the chain; 400 - 395 = 5 is not readable.
  r <- tier1_saa_report(mk_sizes(400L, 398L, 395L, 200L, 150L, 45L), bands_ok,
                        tibble::tibble(target_offset_months = c("0", "12"),
                                       n_matched = c(399L, 399L), median_gap = 0,
                                       p90_gap = 1, share_abs_gap_le_1 = .9))
  s <- sz_of(r)
  expect_true(all(is.na(s[c("n_with_saa_rows", "n_matched")])))
  expect_equal(s[["n_elig"]], 395L)
  expect_equal(audit_rows_d(r), 0L)
})

test_that("tier1_saa_report(): a small arm cell is masked with a complement inside the partition", {
  ## Eligible 300 = treated 3 + control 150 + NA 147. The 3 is small; with the
  ## eligible count published it would be 300 - 150 - 147, so the declared
  ## partition masks the smallest other part, NA (147 < 150): 3 + 147 = 150, not a
  ## small sum. Eligible stays published.
  r <- tier1_saa_report(mk_sizes(500L, 400L, 300L, 3L, 150L, 147L), bands_ok, offsets_ok)
  s <- sz_of(r)
  expect_true(is.na(s[["n_treated"]]) && is.na(s[["n_achieved_na"]]))
  expect_equal(c(s[["n_control"]], s[["n_elig"]]), c(150L, 300L))
  expect_equal(audit_rows_d(r), 0L)
})

test_that("tier1_saa_report(): _YN cells are tied to the published arm sizes", {
  prev <- tibble::tibble(
    column = c("a_YN", "b_YN"),
    n_yes_control = c(100L, 4L), n_no_control = c(50L, 146L), n_na_control = c(0L, 0L),
    n_yes_treated = c(30L, 60L), n_no_treated = c(70L, 40L), n_na_treated = c(0L, 0L)
  ) |>
    dplyr::mutate(
      prevalence_control = .data$n_yes_control / (.data$n_yes_control + .data$n_no_control),
      prevalence_treated = .data$n_yes_treated / (.data$n_yes_treated + .data$n_no_treated)
    )
  ## Arm sizes: control 150 (= 100 + 50 = 4 + 146), treated 100 (= 30 + 70 = 60 + 40).
  r <- tier1_saa_report(mk_sizes(500L, 400L, 300L, 100L, 150L, 50L), bands_ok, offsets_ok, prev)
  p <- r$tables$yn_prevalence
  ## Column a: nothing small, unchanged.
  expect_equal(unlist(p[1, c("n_yes_control", "n_no_control", "n_yes_treated")]),
               c(n_yes_control = 100L, n_no_control = 50L, n_yes_treated = 30L))
  ## Column b, control arm: yes = 4 is small. The control arm size (150) is
  ## published in exact_start_sizes, so 150 - 146 would return the 4: the 146 is
  ## masked as well, and so is the row's derived prevalence.
  expect_true(is.na(p$n_yes_control[2]) && is.na(p$n_no_control[2]))
  expect_true(is.na(p$prevalence_control[2]))
  ## The treated cells of column b are still published.
  expect_equal(c(p$n_yes_treated[2], p$n_no_treated[2]), c(60L, 40L))
  expect_equal(audit_rows_d(r), 0L)
})

test_that("tier1_saa_report(): nearest-period match counts do not reveal a small exact-start gap", {
  ## nearest-period matches 405 (offset 0) vs exact-start matched 400: 5 patients
  ## whose period does not start at s. The declared gaps mask the larger count.
  off <- tibble::tibble(target_offset_months = c("0", "12"), n_matched = c(405L, 470L),
                        median_gap = c(0, 12), p90_gap = c(1, 13), share_abs_gap_le_1 = c(.9, .8))
  r <- tier1_saa_report(mk_sizes(500L, 400L, 300L, 100L, 150L, 50L), bands_ok, off)
  expect_true(is.na(r$tables$gap_by_offset$n_matched[1]))
  expect_true(is.na(r$tables$gap_by_offset$median_gap[1]))   # derived columns go with it
  expect_equal(r$tables$gap_by_offset$n_matched[2], 470L)
  expect_equal(audit_rows_d(r), 0L)
})

test_that("property: random small arm partitions always pass the audit, with no small cell published", {
  vals <- c(0L, 3L, 10L, 11L, 40L)
  grid <- expand.grid(tr = vals, ct = vals, na = vals)
  grid <- grid[rowSums(grid) > 0L, ]
  for (i in seq_len(nrow(grid))) {
    g <- unlist(grid[i, ])
    elig <- sum(g)
    r <- tryCatch(
      tier1_saa_report(mk_sizes(elig + 200L, elig + 100L, elig, g[["tr"]], g[["ct"]], g[["na"]]),
                       bands_ok, offsets_ok),
      error = function(e) e
    )
    lab <- paste("treated, control, NA =", paste(g, collapse = ", "))
    expect_false(inherits(r, "error"), label = lab)
    if (!inherits(r, "error")) {
      expect_equal(audit_rows_d(r), 0L, label = lab)
      v <- r$tables$exact_start_sizes$n
      expect_false(any(!is.na(v) & v > 0L & v < 11L), label = lab)
    }
  }
})

test_that("the Tier-1 report exports through smi_report_export() into a temporary outbox", {
  skip_if(smidata::smi_is_server(), "never export test tables on the server")
  withr::local_envvar(c(SMI_OUTBOX = file.path(withr::local_tempdir(), "outbox")))
  r <- tier1_saa_report(mk_sizes(500L, 400L, 395L, 200L, 150L, 45L), bands_ok, offsets_ok)
  recs <- suppressMessages(suppressWarnings(smidata::smi_report_export(r, "doubletree_tier1")))
  expect_equal(nrow(recs), length(r$tables))     # sizes, msr_yr_offsets, gap_by_offset
  ## the matched count the report masked is NA in the file that leaves.
  csv <- list.files(Sys.getenv("SMI_OUTBOX"), pattern = "exact_start_sizes[.]csv$", full.names = TRUE)
  expect_length(csv, 1L)
  out <- utils::read.csv(csv, comment.char = "#")
  expect_true(is.na(out$n[out$quantity == "n_matched"]))
  expect_equal(out$n[out$quantity == "n_elig"], 395L)
})
