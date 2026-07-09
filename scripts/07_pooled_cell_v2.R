#!/usr/bin/env Rscript
# scripts/07_pooled_cell_v2.R
# Run one pooled cross-language design-analysis cell from a grid row.
# Loads the three per-language v2 pilot DGPs (with posterior draws), simulates
# each language independently, binds them, and fits the cross-language analysis
# model via run_pooled_cell_v2().
# Usage: Rscript scripts/07_pooled_cell_v2.R --row_index N --grid G --dgpdir D --outdir O
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr) })
source("R/03_define_priors.R"); source("R/04_model_formulas.R"); source("R/05_hypothesis_tests.R")
source("R/07_extract_diagnostics.R"); source("R/11_simulate_pooled_v2.R")

LANGS <- c("English", "Turkish", "Norwegian")

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--row_index", default = NULL, type = "integer"),
  optparse::make_option("--grid",      default = "config/design_grid_pooled_v2.csv"),
  optparse::make_option("--dgpdir",    default = "outputs/pilot_models"),
  optparse::make_option("--outdir",    default = "outputs/design_pooled_v2"),
  optparse::make_option("--overwrite", action = "store_true", default = FALSE)
)))
row_idx <- opt$row_index
if (is.null(row_idx)) {
  sid <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
  if (nchar(sid) == 0) stop("[pooled2] --row_index or SLURM_ARRAY_TASK_ID required.")
  row_idx <- as.integer(sid)
}
grid <- readr::read_csv(opt$grid, show_col_types = FALSE)
if (row_idx < 1 || row_idx > nrow(grid)) stop("[pooled2] row ", row_idx, " out of bounds (", nrow(grid), ")")
cell <- grid[row_idx, ]

dgps <- setNames(lapply(LANGS, function(lang) {
  f <- file.path(opt$dgpdir, paste0("pilot_dgp_v2_", cell$prior_source, "_", lang, ".rds"))
  if (!file.exists(f)) stop("[pooled2] missing DGP: ", f)
  readRDS(f)
}), LANGS)

message(strrep("-", 60))
message("[pooled2] Row ", row_idx, " | N/lang ", cell$n_participants,
        " | mode ", cell$mode, " | draws E", cell$draw_index_english,
        "/T", cell$draw_index_turkish, "/N", cell$draw_index_norwegian,
        " | seed ", cell$seed)

t0 <- proc.time()[["elapsed"]]
result <- run_pooled_cell_v2(cell, dgps, out_dir = opt$outdir, overwrite = opt$overwrite)
message("[pooled2] wall ", round(proc.time()[["elapsed"]] - t0, 1), "s | status ",
        if (is.list(result)) result$summary$status else "?")
