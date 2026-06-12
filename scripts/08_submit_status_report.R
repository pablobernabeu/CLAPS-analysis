#!/usr/bin/env Rscript
# scripts/08_submit_status_report.R
# Generate a status report comparing the design grid against completed outputs,
# query SLURM for job statuses if available, and write a summary CSV.
# Usage: Rscript scripts/08_submit_status_report.R

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
