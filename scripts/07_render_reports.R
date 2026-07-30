#!/usr/bin/env Rscript
# scripts/07_render_reports.R
#
# Purpose
#   Render the Quarto reports, but only once the bibliography has been verified.
#   The audit check at the top is the point of this wrapper: it makes an unverified
#   citation a build failure rather than something a reader might notice after
#   publication.
#
# Preconditions
#   scripts/00_verify_references.R must have run and left a clean audit CSV. The
#   check re-reads that CSV rather than trusting its existence, so an audit that
#   ran and failed also blocks the render. The blocking flags are the same three
#   used in R/00_reference_audit.R; WARN rows do not block.
#
#   Aggregation should also have run, since the reports read the summary CSVs, but
#   that is not enforced here: a missing input surfaces as a Quarto error naming the
#   file, which is clear enough.
#
# Inputs
#   --report all|preliminary   Which report to render. "all" renders every entry in
#                              the `reports` list below.
#   CLAPS_OUTPUTS_ROOT (environment variable) to locate the audit CSV.
#
# Usage
#   Rscript scripts/00_verify_references.R && Rscript scripts/07_render_reports.R
#
# Behaviour on an unknown or missing report
#   Warns and moves on rather than stopping, so that rendering several reports is
#   not aborted by one absent file. Read the warnings: the final "All reports
#   complete" message is printed even when a report was skipped.

suppressPackageStartupMessages({
  library(quarto)
  library(optparse)
})

option_list <- list(
  optparse::make_option("--report", default = "all",
    help = "Which report to render: all, preliminary")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Verify reference audit output exists and is clean
audit_csv <- file.path(Sys.getenv("CLAPS_OUTPUTS_ROOT", "outputs"),
                       "reference_audit", "reference_audit.csv")
if (!file.exists(audit_csv)) {
  stop("[render] Reference audit not found. Run scripts/00_verify_references.R first.")
}
audit_df <- readr::read_csv(audit_csv, show_col_types = FALSE)
n_errors <- sum(audit_df$flag %in% c("ERROR", "TITLE_MISMATCH", "YEAR_MISMATCH"), na.rm = TRUE)
if (n_errors > 0) {
  stop("[render] Reference audit has ", n_errors, " failures. Fix references.bib first.")
}

reports <- list(
  preliminary = "reports/preliminary_sample_size_analysis.qmd"
)

to_render <- if (opt$report == "all") names(reports) else opt$report

for (nm in to_render) {
  qmd_path <- reports[[nm]]
  if (is.null(qmd_path)) {
    warning("[render] Unknown report name: ", nm)
    next
  }
  if (!file.exists(qmd_path)) {
    warning("[render] QMD not found: ", qmd_path)
    next
  }
  message("[render] Rendering: ", qmd_path)
  quarto::quarto_render(qmd_path)
  message("[render] Done: ", qmd_path)
}

message("[render] All reports complete.")
