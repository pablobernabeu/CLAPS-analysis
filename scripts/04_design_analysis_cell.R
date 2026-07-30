#!/usr/bin/env Rscript
# scripts/04_design_analysis_cell.R
#
# Purpose
#   Run exactly one design-analysis cell. This is the unit of work a SLURM array
#   task executes: the array index selects a row of the design grid, and that row
#   fully determines what is simulated and fitted.
#
# Why one cell per process
#   Cells are independent, so an array of single-cell jobs parallelises across the
#   cluster with no coordination, and a cell that exhausts its memory or walltime
#   takes down only itself. It also means a partial run is resumable: rerunning the
#   array skips cells whose .rds already exists (see `overwrite` below).
#
# Inputs
#   --row_index N   1-based row of the grid. Defaults to SLURM_ARRAY_TASK_ID, which
#                   is why an array must be submitted with 1-based indices
#                   (--array=1-N), not 0-based.
#   --grid FILE     Design grid CSV. Default config/design_grid.csv.
#   --config FILE   config/analysis_config.yaml, consulted only for per-language
#                   settings the grid does not carry.
#   --outdir DIR    Where the per-cell .rds is written.
#   --overwrite     Recompute a cell even if its output exists.
#
# Output
#   One .rds in --outdir, named from the cell's parameters by run_design_cell().
#   Written atomically, so an interrupted job leaves no partial file.
#
# Usage
#   Rscript scripts/04_design_analysis_cell.R --row_index 1
#   sbatch hpc/submit_design_analysis_array.sh      # the normal route
#
# Exit behaviour
#   An out-of-range row index is an error, which surfaces a mismatch between the
#   --array range in the submission script and the grid's row count. A cell whose
#   model fails to fit is NOT an error: it writes an .rds with status = "error" and
#   exits 0, because a failure to converge is a result the aggregation must see
#   rather than a job fault. Judge a run by the statuses in its outputs, not by the
#   SLURM exit codes.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(yaml)
  library(optparse)
})

source("R/03_define_priors.R")
source("R/04_model_formulas.R")
source("R/05_hypothesis_tests.R")
source("R/06_simulate_design.R")
source("R/07_extract_diagnostics.R")

option_list <- list(
  optparse::make_option("--row_index", default = NULL, type = "integer",
    help = "Row in design_grid.csv (1-based). Defaults to SLURM_ARRAY_TASK_ID."),
  optparse::make_option("--grid",      default = "config/design_grid.csv"),
  optparse::make_option("--config",    default = "config/analysis_config.yaml"),
  optparse::make_option("--outdir",    default = "outputs/design_analysis"),
  optparse::make_option("--overwrite", action  = "store_true", default = FALSE)
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# Resolve the row index: an explicit --row_index wins, so a single cell can be
# rerun by hand for debugging; otherwise fall back to the SLURM array index. The
# error rather than a default of 1 is deliberate, since silently running row 1
# would look like success while doing the wrong work.
row_idx <- opt$row_index
if (is.null(row_idx)) {
  slurm_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  if (nchar(slurm_id) == 0) stop("[cell] --row_index or SLURM_ARRAY_TASK_ID required.")
  row_idx <- as.integer(slurm_id)
}

grid <- readr::read_csv(opt$grid, show_col_types = FALSE)
if (row_idx < 1 || row_idx > nrow(grid)) {
  stop("[cell] row_index ", row_idx, " out of bounds (grid has ", nrow(grid), " rows).")
}

cell <- grid[row_idx, ]
cfg  <- yaml::read_yaml(opt$config)

# The grid may carry has_pseudo_passive per row; when it does not, fall back to the
# per-language setting in the config. isTRUE() means an unknown language yields
# FALSE rather than NULL, which would fail later inside the ladder builder.
if (!"has_pseudo_passive" %in% names(cell)) {
  cell$has_pseudo_passive <- isTRUE(
    cfg$languages[[cell$language]]$has_pseudo_passive
  )
}

# Provenance header. These lines exist so that the SLURM log alone is enough to
# identify what ran: the job and array IDs tie the log to the scheduler's records,
# and the git SHA ties it to the code. Without the SHA in the log, a result found
# months later cannot be matched to the version that produced it.
message(replicate(60, "-") |> paste(collapse = ""))
message("[cell] Job: ", Sys.getenv("SLURM_JOB_ID", "local"),
        " | Array task: ", Sys.getenv("SLURM_ARRAY_TASK_ID", row_idx),
        " | Row: ", row_idx)
message("[cell] Git SHA: ", tryCatch(
  system2("git", c("rev-parse","--short","HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(e) "unknown"))
message("[cell] Language: ",    cell$language,
        " | Model: ",          cell$model_level,
        " | Prior: ",          cell$prior_regime,
        " | Threshold: ",      cell$threshold_mode,
        " | N participants: ", cell$n_participants,
        " | N verbs: ",        cell$n_verbs,
        " | Seed: ",           cell$seed)

t_total_start <- proc.time()[["elapsed"]]

result <- run_design_cell(cell, out_dir = opt$outdir, overwrite = opt$overwrite)

t_total <- proc.time()[["elapsed"]] - t_total_start
message("[cell] Wall time: ", round(t_total, 1), "s")
message("[cell] Status: ", if (is.list(result)) result$summary$status else result$status)
