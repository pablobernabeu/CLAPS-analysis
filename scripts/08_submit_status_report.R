#!/usr/bin/env Rscript
# scripts/08_submit_status_report.R
#
# Purpose
#   Answer "how far has the design analysis got?" by comparing the design grid
#   against the outputs actually on disk, and cross-checking against SLURM where
#   its records are available.
#
# Inputs (hard-coded, not options)
#   config/design_grid.csv         The grid whose completion is being measured.
#   outputs/design_analysis/       The per-cell .rds files.
#   outputs/slurm_job_ids.csv      Optional. If present, each job_id is queried via
#                                  sacct and its state printed.
#
# Output
#   outputs/job_status_report.csv, one row per grid row with a done/pending status.
#
# Usage
#   Rscript scripts/08_submit_status_report.R
#
# One caveat before trusting the numbers
#   "done" means the cell wrote an output, not that its model converged. A cell that
#   errored still writes an .rds and is counted as done here. Read the status column
#   of the aggregated summary (R/08_summarise_design.R) for convergence.
#
#   The filename-suffix problem that used to make this script understate progress for
#   the gender grids was fixed on 2026-07-30; check_grid_completion() now reproduces
#   run_design_cell()'s naming in full.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/10_job_status.R")

grid_path <- "config/design_grid.csv"
out_dir   <- "outputs/design_analysis"

if (!file.exists(grid_path)) {
  stop("[status_report] Design grid not found: ", grid_path)
}

grid      <- readr::read_csv(grid_path, show_col_types = FALSE)
completed <- list_completed_cells(out_dir)
status    <- check_grid_completion(grid, completed)

write_status_report(status, out_path = "outputs/job_status_report.csv")

# Optional: query SLURM if job IDs are available
slurm_jobs_path <- "outputs/slurm_job_ids.csv"
if (file.exists(slurm_jobs_path)) {
  jobs <- readr::read_csv(slurm_jobs_path, show_col_types = FALSE)
  purrr::walk(jobs$job_id, function(jid) {
    sq <- query_slurm_status(jid)
    if (!is.null(sq)) {
      message("[status] Job ", jid, ": ",
              paste(unique(sq$State), collapse = ", "))
    }
  })
}

message("[status_report] Complete.")
