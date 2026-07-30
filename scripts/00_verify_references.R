#!/usr/bin/env Rscript
# scripts/00_verify_references.R
#
# Purpose
#   Command-line wrapper around R/00_reference_audit.R, which verifies every
#   citation in the bibliography against Crossref. Numbered 00 because it is a
#   gate: scripts/07_render_reports.R refuses to render until this has run and
#   passed.
#
# Inputs
#   --bib PATH   Bibliography to check. Defaults to references.bib.
#   CLAPS_OUTPUTS_ROOT (environment variable) relocates the output directory,
#   which is how the CI job writes its audit somewhere other than outputs/.
#
# Output
#   <outputs root>/reference_audit/reference_audit.csv
#
# Usage
#   Rscript scripts/00_verify_references.R
#
# Exit codes
#   0 on a clean audit, 1 on any failure. The explicit quit() calls matter: this
#   script exists to be a gate in a shell pipeline and in CI, both of which decide
#   whether to continue from the exit status. Requires network access to Crossref.

source("R/00_reference_audit.R")

# Argument parsing is done by hand rather than with optparse, so that the gate has
# no package dependency beyond what the audit itself needs.
args <- commandArgs(trailingOnly = TRUE)
bib_path <- "references.bib"
for (i in seq_along(args)) {
  # The i < length(args) guard means a trailing "--bib" with no value is ignored
  # rather than reading past the end of the vector.
  if (args[i] == "--bib" && i < length(args)) bib_path <- args[i + 1L]
}

tryCatch({
  audit_out <- file.path(Sys.getenv("CLAPS_OUTPUTS_ROOT", "outputs"), "reference_audit")
  run_reference_audit(bib_path = bib_path, out_dir = audit_out)
  message("[00_verify_references] Reference audit passed.")
  quit(status = 0)
}, error = function(e) {
  message("[00_verify_references] AUDIT FAILED: ", conditionMessage(e))
  quit(status = 1)
})
