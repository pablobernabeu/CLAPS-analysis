#!/usr/bin/env Rscript
# scripts/06_databased_cell_v2.R
# Run one v2 data-grounded cell (assurance / safeguard / point) from a grid row.
# Loads the v2 DGP (with posterior draws) for the row's (prior_source, language),
# then simulates -> refits -> Bayes factors via run_databased_cell_v2().
# Usage: Rscript scripts/06_databased_cell_v2.R --row_index N --grid G --dgpdir D --outdir O
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr); library(yaml) })
source("R/03_define_priors.R"); source("R/04_model_formulas.R"); source("R/05_hypothesis_tests.R")
source("R/07_extract_diagnostics.R"); source("R/10_simulate_from_pilot_v2.R")

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--row_index", default = NULL, type = "integer"),
  optparse::make_option("--grid",      default = "config/design_grid_databased_v2.csv"),
  optparse::make_option("--dgpdir",    default = "outputs/pilot_models"),
  optparse::make_option("--outdir",    default = "outputs/design_databased_v2"),
  optparse::make_option("--overwrite", action = "store_true", default = FALSE)
)))
row_idx <- opt$row_index
if (is.null(row_idx)) {
  sid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  if (nchar(sid) == 0) stop("[cell2] --row_index or SLURM_ARRAY_TASK_ID required.")
  row_idx <- as.integer(sid)
}
grid <- readr::read_csv(opt$grid, show_col_types = FALSE)
if (row_idx < 1 || row_idx > nrow(grid)) stop("[cell2] row ", row_idx, " out of bounds (", nrow(grid), ")")
cell <- grid[row_idx, ]

dgp_file <- file.path(opt$dgpdir, paste0("pilot_dgp_v2_", cell$prior_source, "_", cell$language, ".rds"))
if (!file.exists(dgp_file)) stop("[cell2] missing DGP: ", dgp_file)
dgp <- readRDS(dgp_file)

message(strrep("-", 60))
message("[cell2] Row ", row_idx, " | ", cell$language, " | source ", cell$prior_source,
        " | mode ", cell$mode, " | N ", cell$n_participants,
        " | draw ", cell$draw_index, " | seed ", cell$seed)

t0 <- proc.time()[["elapsed"]]
result <- run_databased_cell_v2(cell, dgp, out_dir = opt$outdir, overwrite = opt$overwrite)
message("[cell2] wall ", round(proc.time()[["elapsed"]] - t0, 1), "s | status ",
        if (is.list(result)) result$summary$status else "?")
