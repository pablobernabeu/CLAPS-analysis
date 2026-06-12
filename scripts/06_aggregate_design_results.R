#!/usr/bin/env Rscript
# scripts/06_aggregate_design_results.R
# Aggregate all design-analysis cells into summary tables and manifests.
# Run after all SLURM array jobs complete (via afterok dependency).
# Usage: Rscript scripts/06_aggregate_design_results.R [--out_dir PATH] [--sum_dir PATH]

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(optparse)
})

source("R/08_summarise_design.R")
source("R/10_job_status.R")

option_list <- list(
  optparse::make_option("--out_dir", default = "outputs/design_analysis",
    help = "Directory containing .rds cell outputs [default: outputs/design_analysis]"),
  optparse::make_option("--sum_dir", default = "outputs/design_summary",
    help = "Directory to write summary tables [default: outputs/design_summary]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

cfg     <- yaml::read_yaml("config/analysis_config.yaml")
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

`%||%` <- function(a, b) if (!is.null(a)) a else b
