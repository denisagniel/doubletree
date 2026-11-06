## ============================================================================
## 93_ss_desc_census.R
##
## CENSUS of the service-setting field (SS_DESC) in the service cost-claims
## files, so the match rules for the claims-derived outcomes in 92b A1 (any
## psychiatric inpatient admission, any all-cause inpatient admission, any ED
## visit) can be written against the real values instead of guessed.
## Reads no outcome and estimates nothing. STANDALONE: not part of the 01-07
## cascade.
##
## Per SS_DESC value it reports: claim lines, distinct patients, and the share of
## lines whose ADMIT_DT differs from SRV_DT (smidata's COST.md: an inpatient stay
## has ADMIT_DT != SRV_DT). Two families (aim3_svc, new_svc) are read for ONE
## year (default 22; override on the command line), enough to see the values.
##
## Written to output/export/ (aggregates). Any row with fewer than
## config_report_min_cell lines OR patients is suppressed whole (its numbers are
## NA; the SS_DESC label is kept so the value is still visible). Look at the file
## before pulling it back.
##   ss_desc_census.csv
##
## SERVER-ONLY: reads patient-level files.
##   Rscript application/93_ss_desc_census.R          # year 22
##   Rscript application/93_ss_desc_census.R 19       # another suffix year
## ============================================================================

app_dir <- if (dir.exists("application")) "application" else "."
source(file.path(app_dir, "_config.R"))
source(file.path(app_dir, "helpers", "disclosure.R"))   # config_report_min_cell

cli::cli_h1("93 -- SS_DESC census (setting values in the service claims)")

if (!config_has_smidata || identical(smidata::smi_env(), "local")) {
  cli::cli_abort(c(
    "This census reads the real claims files and must run on the server.",
    "i" = "Locally {.fn smidata::smi_read} would return a synthetic fixture."
  ))
}

year <- suppressWarnings(as.integer(commandArgs(trailingOnly = TRUE)[1]))
if (is.na(year)) year <- 22L
if (!year %in% config_cost_years) {
  cli::cli_abort("Year {.val {year}} is not in config_cost_years ({.val {config_cost_years}}).")
}
families <- unlist(SMI_KEYS$cost_claims, use.names = FALSE)
keys <- paste0(families, "_", year)

census_one <- function(key) {
  cli::cli_alert_info("Reading {.val {key}} (SS_DESC, SRV_DT, ADMIT_DT, ID).")
  d <- smi_read_pinned(key, columns = c("ID", config_cost_setting_col, "SRV_DT", "ADMIT_DT"))
  d <- as.data.frame(d)
  setting <- as.character(d[[config_cost_setting_col]])
  setting[is.na(setting)] <- "<NA>"
  admit_differs <- as.Date(d$ADMIT_DT) != as.Date(d$SRV_DT)
  by_val <- split(seq_len(nrow(d)), setting)
  out <- do.call(rbind, lapply(names(by_val), function(v) {
    i <- by_val[[v]]
    data.frame(
      file = key,
      ss_desc = v,
      n_lines = length(i),
      n_patients = length(unique(d$ID[i])),
      share_admit_ne_srv = mean(admit_differs[i], na.rm = TRUE),
      share_admit_missing = mean(is.na(d$ADMIT_DT[i])),
      stringsAsFactors = FALSE
    )
  }))
  rm(d); invisible(gc())
  out
}

census <- do.call(rbind, lapply(keys, census_one))
census <- census[order(census$file, -census$n_lines), ]

small <- census$n_lines < config_report_min_cell | census$n_patients < config_report_min_cell
num_cols <- c("n_lines", "n_patients", "share_admit_ne_srv", "share_admit_missing")
census[small, num_cols] <- NA
census$suppressed <- small

print(census, row.names = FALSE, digits = 3)

out_path <- output_file("export", "ss_desc_census.csv", app_dir = app_dir)
utils::write.csv(census, out_path, row.names = FALSE)
cli::cli_alert_success("Wrote {.file {out_path}} ({nrow(census)} rows; {sum(small)} suppressed).")
cli::cli_alert_warning("Read the file before pulling it back: suppression is simple and not audited by smi_report().")
