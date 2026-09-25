## ============================================================================
## application/tests/test-no-unpinned-reads.R
##
## Enforcement: every read in this pipeline must go through smi_read_pinned()
## (application/_config.R), never smidata::smi_read() / smi_read() directly.
## A direct call would silently read the newest ingest instead of the one
## config_ingest_id pins, and (independently, see test-smi-read-pinned.R)
## would ABORT locally the moment anyone set config_ingest_id to a real
## value -- exactly the footgun smi_read_pinned() exists to route around.
##
## Static source-text check, not a behavioral one: it greps the pipeline's
## own numbered scripts (and run_pipeline.R/run_tests.R, harmlessly -- they
## contain no smi_read() call at all) for the literal call. _config.R is
## excluded -- it is where smi_read_pinned() is DEFINED, so its own body
## legitimately contains the one real smidata::smi_read() call this test
## enforces against everywhere else. application/helpers/ and
## application/tests/ are excluded by construction: the glob below is
## application/'s own top-level *.R files, non-recursive, so it never
## descends into either subdirectory.
## ============================================================================

test_that("no pipeline script calls smi_read() directly instead of smi_read_pinned()", {
  pipeline_files <- list.files(app_dir_for_tests, pattern = "\\.R$", full.names = TRUE)
  pipeline_files <- pipeline_files[!grepl("_config\\.R$", pipeline_files)]
  ## Guards the guard: an empty file list would make the assertion below
  ## vacuously pass, which would hide a broken glob rather than a clean repo.
  expect_gt(length(pipeline_files), 0L)

  src <- unlist(lapply(pipeline_files, readLines, warn = FALSE))
  src <- src[!grepl("^\\s*#", src)]   # drop full-line comments (this file's own `##` convention)

  hits <- grep("smi(data::)?_read\\(", src, value = TRUE)
  hits <- hits[!grepl("smi_read_pinned", hits)]   # belt-and-suspenders on top of the regex itself

  expect_equal(hits, character(0))
})
