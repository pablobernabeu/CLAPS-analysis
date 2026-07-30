#!/usr/bin/env Rscript
# scripts/07_pooled_cell_v2.R
#
# Purpose
#   One cell of the pooled cross-language design analysis. Each language is
#   simulated from its OWN pilot data-generating spec, the three data sets are
#   bound together, and a single cross-language model is fitted to the result.
#
# Why simulate per language and pool afterwards
#   The three languages differ in their pilot-estimated effects, thresholds and
#   variance components, and Norwegian additionally lacks pseudo-passives.
#   Generating each from its own spec preserves that heterogeneity, so the pooled
#   fit faces the same between-language variation the real analysis will. Simulating
#   from one averaged spec would understate it and overstate the power of the
#   cross-language model.
#
# Inputs
#   --row_index N  Row of the grid; defaults to SLURM_ARRAY_TASK_ID.
#   --grid FILE    Default config/design_grid_pooled_v2.csv.
#   --dgpdir DIR   Holds pilot_dgp_v2_<prior_source>_<language>.rds for all three
#                  languages. All three must be present; a missing one is an error
#                  rather than a silently smaller pool.
#   --outdir DIR, --overwrite
#
# Note on the grid columns
#   The grid carries a separate draw index per language
#   (draw_index_english/turkish/norwegian). They are independent because the three
#   pilot posteriors are independent: there is no joint posterior to draw a single
#   index from, and tying them to one index would impose a correspondence between
#   unrelated draws.
#
# Usage
#   Rscript scripts/07_pooled_cell_v2.R --row_index 1
#   sbatch hpc/submit_pooled_v2_array.sh
#
# Cost
#   The most expensive cell type in the repository. A cross-language fit takes
#   roughly 8-12 hours, which is why config/analysis_config.yaml sets a lower
#   replicate count for cross-language design points than for single-language ones.

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr) })
source("R/03_define_priors.R"); source("R/04_model_formulas.R"); source("R/05_hypothesis_tests.R")
source("R/07_extract_diagnostics.R"); source("R/11_simulate_pooled_v2.R")

# The three CLAPS pilot languages, hard-coded because the pooled grids are built
# around exactly these three and the grid's per-language draw-index columns are
# named after them. Adding a language means updating both this vector and the grid
# generator scripts/generate_pooled_v2_grid.R.
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
