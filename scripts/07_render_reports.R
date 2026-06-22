#!/usr/bin/env Rscript
# scripts/07_render_reports.R
# Render all QMD reports after reference audit and aggregation are complete.
# Fails immediately if the reference audit has not been run or has failed.
# Usage: Rscript scripts/07_render_reports.R [--report all|preliminary]

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
