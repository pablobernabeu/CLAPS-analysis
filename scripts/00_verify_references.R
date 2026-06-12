#!/usr/bin/env Rscript
# scripts/00_verify_references.R
# Run the reference audit before any report rendering.
# Fails with a non-zero exit code if any citation fails verification.
# Usage: Rscript scripts/00_verify_references.R [--bib PATH]

source("R/00_reference_audit.R")

args <- commandArgs(trailingOnly = TRUE)
bib_path <- "references.bib"
for (i in seq_along(args)) {
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
