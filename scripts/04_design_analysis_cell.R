#!/usr/bin/env Rscript
# scripts/04_design_analysis_cell.R
# Run a single design-analysis cell identified by a row index into config/design_grid.csv.
# Called as a SLURM array task: SLURM_ARRAY_TASK_ID maps to the row index.
# Usage: Rscript scripts/04_design_analysis_cell.R [--row_index N] [--overwrite]

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

# Resolve row index
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

# Merge config defaults into cell if columns missing
if (!"has_pseudo_passive" %in% names(cell)) {
  cell$has_pseudo_passive <- isTRUE(
    cfg$languages[[cell$language]]$has_pseudo_passive
  )
}

# Log header
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
