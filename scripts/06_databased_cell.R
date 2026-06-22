#!/usr/bin/env Rscript
# scripts/06_databased_cell.R
# Run one cell of the data-grounded BAYESIAN power analysis: load the pilot DGP
# for the cell's language, simulate from it, refit the model, compute BFs, save.
suppressPackageStartupMessages({ library(optparse); library(readr); library(dplyr) })
source("R/10_simulate_from_pilot.R")

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--row_index", type = "integer"),
  optparse::make_option("--grid",      default = "config/design_grid_databased.csv"),
  optparse::make_option("--dgpdir",    default = "outputs/pilot_models"),
  optparse::make_option("--outdir",    default = "outputs/design_databased"),
  optparse::make_option("--overwrite", action = "store_true", default = FALSE)
)))

grid <- readr::read_csv(opt$grid, show_col_types = FALSE)
if (opt$row_index < 1 || opt$row_index > nrow(grid)) {
  stop("[databased] row_index ", opt$row_index, " out of range (", nrow(grid), ")")
}
cell <- as.list(grid[opt$row_index, ])

dgp_file <- file.path(opt$dgpdir, paste0("pilot_dgp_", cell$language, ".rds"))
if (!file.exists(dgp_file)) stop("[databased] missing pilot DGP: ", dgp_file)
dgp <- readRDS(dgp_file)

cat("[databased] row", opt$row_index, "|", cell$language, "N", cell$n_participants,
    "mult", cell$effect_mult, "seed", cell$seed, "\n")
run_databased_cell(cell, dgp, opt$outdir, overwrite = opt$overwrite)
