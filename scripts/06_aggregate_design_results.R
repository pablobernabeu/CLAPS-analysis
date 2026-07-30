#!/usr/bin/env Rscript
# scripts/06_aggregate_design_results.R
#
# Purpose
#   Collect every per-cell .rds written by the design-analysis array into the
#   summary CSVs the report reads, and record the run's provenance.
#
# When to run
#   After the array has finished. hpc/submit_aggregate_afterok.sh submits this with
#   a SLURM afterok dependency so it starts automatically once the array completes.
#   Running it early is not harmful but yields a summary over whatever subset had
#   finished, which is misleading unless the cell counts are checked.
#
# Inputs
#   --out_dir  Directory of per-cell .rds files. Default outputs/design_analysis.
#   --sum_dir  Destination for the summary CSVs. Default outputs/design_summary.
#
# Outputs
#   The seven CSVs documented in R/08_summarise_design.R, plus outputs/manifest.csv
#   recording the git SHA, package versions and the cell counts behind the summary.
#
# Usage
#   Rscript scripts/06_aggregate_design_results.R
#
# Note on interpreting the counts
#   n_design_cells in the manifest counts ROWS, not cells: a cell contributes one
#   row per hypothesis tested. Compare the cell count in the console message
#   against the grid's row count to detect cells that never wrote an output at all,
#   which is the one failure mode invisible to the summaries themselves.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
})

# Null-coalescing operator, defined before anything below can use it. Base R gained
# %||% only in 4.4, and config/arc_modules.yaml sets min_version 4.3.0 while CI pins
# 4.3.3, so it cannot be assumed present. Neither module sourced below defines it,
# so this definition is the only one in force for this script.
`%||%` <- function(a, b) if (!is.null(a)) a else b

source("R/08_summarise_design.R")
source("R/10_job_status.R")

option_list <- list(
  optparse::make_option("--out_dir", default = "outputs/design_analysis",
    help = "Directory containing .rds cell outputs [default: outputs/design_analysis]"),
  optparse::make_option("--sum_dir", default = "outputs/design_summary",
    help = "Directory to write summary tables [default: outputs/design_summary]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

out_dir <- opt$out_dir
sum_dir <- opt$sum_dir

message("[aggregate] Loading design cells from: ", out_dir)
df <- load_design_cells(out_dir)
message("[aggregate] ", nrow(df), " rows loaded across ",
        dplyr::n_distinct(df$cell_id %||% seq_len(nrow(df))), " cells.")

# Write summary tables
write_design_summary(df, out_dir = sum_dir)

# Write manifest
write_manifest(
  out_path        = "outputs/manifest.csv",
  additional_cols = list(
    n_design_cells  = nrow(df),
    n_success       = sum(df$status == "success", na.rm = TRUE),
    aggregated_at   = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
)

message("[aggregate] Done.")
