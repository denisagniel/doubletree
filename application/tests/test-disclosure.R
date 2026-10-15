## ============================================================================
## application/tests/test-disclosure.R
##
## helpers/disclosure.R: the property under test is that NO masked cell can be
## recovered from the published ones by a single subtraction, not that a given
## literal appears. Every expected value is worked out by hand in a comment.
## ============================================================================

## Published cells that must sum to a published total leave exactly one unknown
## if one is masked; the helpers must leave zero or at least two.
n_unknown <- function(x) sum(is.na(x))

test_that("mask_partition_cells() never leaves a single masked cell next to a total", {
  ## 40, 5, 60: the 5 is small; the smallest other positive cell (40) goes too.
  expect_equal(mask_partition_cells(c(40L, 5L, 60L), 11L), c(NA, NA, 60L))
  ## Two small cells already leave two unknowns: nothing else is masked.
  expect_equal(mask_partition_cells(c(40L, 5L, 7L), 11L), c(40L, NA, NA))
  ## A zero is not a small cell and is not a complement.
  expect_equal(mask_partition_cells(c(0L, 5L, 60L), 11L), c(0L, NA, NA))
  ## Nothing small: unchanged.
  expect_equal(mask_partition_cells(c(40L, 12L, 60L), 11L), c(40L, 12L, 60L))
  ## Over a grid of small partitions the number masked is never exactly 1, except
  ## when the only positive cell is the small one (nothing to pair it with; the
  ## caller must then mask the total, as tier1_disclosable_sizes() does).
  grid <- expand.grid(a = 0:30, b = 0:30, c = 0:30)
  has_partner <- rowSums(grid > 0L) >= 2L
  nm <- apply(grid[has_partner, ], 1L,
              function(r) n_unknown(mask_partition_cells(as.integer(r), 11L)))
  expect_false(any(nm == 1L))
})

test_that("tier1_disclosable_sizes() keeps chain and partition identities unrecoverable", {
  mk <- function(rows, matched, elig, tr, ct, na) {
    tibble::tibble(n_with_saa_rows = rows, n_matched = matched, n_elig = elig,
                   n_treated = tr, n_control = ct, n_achieved_na = na)
  }
  ## Nothing small: published unchanged.
  ok <- mk(500L, 400L, 300L, 100L, 150L, 50L)
  expect_equal(tier1_disclosable_sizes(ok), ok)
  ## matched 400 -> elig 395: a drop of 5. The inner count (matched) is hidden;
  ## 500 - 395 = 105 is a fine cell, and nothing shows 400 - 395.
  gap <- tier1_disclosable_sizes(mk(500L, 400L, 395L, 200L, 150L, 45L))
  expect_true(is.na(gap$n_matched))
  expect_equal(gap$n_elig, 395L)
  ## 400 -> 398 -> 395: drops of 2 and 3. Hiding matched leaves 400 -> 395, a
  ## drop of 5; with only two visible counts the last (elig) goes, and with it the
  ## arms whose sum would give elig back.
  tiny <- tier1_disclosable_sizes(mk(400L, 398L, 395L, 200L, 150L, 45L))
  expect_true(all(is.na(c(tiny$n_matched, tiny$n_elig))))
  expect_true(all(is.na(c(tiny$n_treated, tiny$n_control, tiny$n_achieved_na))))
  ## A small arm cell (3 of 300): masked with its complement, and elig hidden so
  ## elig - treated - control cannot give the masked cell back.
  arm <- tier1_disclosable_sizes(mk(500L, 400L, 300L, 3L, 150L, 147L))
  expect_true(is.na(arm$n_treated))
  expect_true(is.na(arm$n_elig))
  expect_gte(n_unknown(c(arm$n_treated, arm$n_control, arm$n_achieved_na, arm$n_elig)), 2L)
})

test_that("tier1_disclosable_prevalence() masks per (column, arm) partition", {
  prev <- tibble::tibble(
    column = c("a_YN", "b_YN"),
    n_yes_control = c(100L, 4L), n_no_control = c(50L, 146L), n_na_control = c(0L, 0L),
    n_yes_treated = c(30L, 60L), n_no_treated = c(70L, 40L), n_na_treated = c(0L, 0L)
  ) |>
    dplyr::mutate(
      prevalence_control = .data$n_yes_control / (.data$n_yes_control + .data$n_no_control),
      prevalence_treated = .data$n_yes_treated / (.data$n_yes_treated + .data$n_no_treated)
    )
  out <- tier1_disclosable_prevalence(prev)
  ## Column a: nothing small, unchanged.
  expect_equal(unlist(out[1, c("n_yes_control", "n_no_control", "n_yes_treated")]),
               c(n_yes_control = 100L, n_no_control = 50L, n_yes_treated = 30L))
  ## Column b, control arm: yes = 4 is small. The arm has 150 patients, which is
  ## published, so 150 - 146 would give the 4 back: the 146 is masked as well.
  expect_true(is.na(out$n_yes_control[2]) && is.na(out$n_no_control[2]))
  expect_true(is.na(out$prevalence_control[2]))
  ## The treated arm of column b is untouched.
  expect_equal(out$n_yes_treated[2], 60L)
  expect_equal(out$prevalence_treated[2], 0.6)
})
