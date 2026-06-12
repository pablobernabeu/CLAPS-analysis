#!/usr/bin/env Rscript
# scripts/05_bf_calibration_cell.R
# Bayes-factor calibration: fit a small subset of design cells using both
# Savage-Dickey (primary) and bridge sampling (calibration check).
# Outputs a comparison table for the statistical lead to verify agreement.
# Only run for cells where SLURM_ARRAY_TASK_ID is in the calibration subset.

suppressPackageStartupMessages({
  library(brms)
  library(bridgesampling)
  library(dplyr)
  library(tibble)
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
  optparse::make_option("--row_index", default = NULL, type = "integer"),
  optparse::make_option("--grid",    default = "config/design_grid.csv"),
  optparse::make_option("--config",  default = "config/analysis_config.yaml"),
  optparse::make_option("--outdir",  default = "outputs/bf_calibration"),
  optparse::make_option("--overwrite", action = "store_true", default = FALSE)
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

row_idx <- opt$row_index %||% as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
grid <- readr::read_csv(opt$grid, show_col_types = FALSE)
cell <- grid[row_idx, ]
cfg  <- yaml::read_yaml(opt$config)

if (!"has_pseudo_passive" %in% names(cell)) {
  cell$has_pseudo_passive <- isTRUE(cfg$languages[[cell$language]]$has_pseudo_passive)
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
cell_id  <- paste(cell$language, cell$model_level, cell$prior_regime,
                  cell$threshold_mode, cell$n_participants, cell$n_verbs, cell$seed,
                  sep = "_")
out_path <- file.path(opt$outdir, paste0(cell_id, "_calibration.csv"))

if (file.exists(out_path) && !opt$overwrite) {
  message("[calibration] Skipping existing: ", cell_id)
  quit(status = 0)
}

message("[calibration] Running calibration for: ", cell_id)

has_pp <- isTRUE(cell$has_pseudo_passive)

sim_data <- simulate_claps_data(
  n_participants = cell$n_participants,
  n_verbs        = cell$n_verbs,
  has_pseudo_passive = has_pp,
  seed           = cell$seed
)

prior_obj <- build_brms_prior(cell$prior_regime, cell$threshold_mode,
                               has_pseudo_passive = has_pp)
ladder    <- build_model_ladder(has_pseudo_passive = has_pp)
formula   <- ladder[[cell$model_level]]

samp_args <- production_sampling(seed = cell$seed,
  cores = as.integer(Sys.getenv("STAN_NUM_THREADS", "4")))
ctrl_args <- production_control()

fit <- do.call(brms::brm, c(
  list(formula = formula, data = sim_data, prior = prior_obj,
       backend = "cmdstanr", sample_prior = "yes", control = ctrl_args),
  samp_args
))

# 1. Savage-Dickey BFs (primary)
sd_bfs <- compute_all_bf(fit, has_pseudo_passive = has_pp)

# 2. Bridge sampling BFs (calibration only)
# Requires a null model per hypothesis; we fit the null for the main semantics slope
message("[calibration] Fitting null model for bridge sampling calibration...")
formula_null <- ladder[["L0_random_intercepts_only"]]
fit_null <- do.call(brms::brm, c(
  list(formula = formula_null, data = sim_data, prior = prior_obj,
       backend = "cmdstanr", sample_prior = "yes", control = ctrl_args),
  samp_args
))

bridge_bfs <- tryCatch({
  bf_bs <- brms::bayes_factor(fit, fit_null)
  tibble::tibble(
    hypothesis     = "overall_vs_intercepts_only",
    BF_bs          = bf_bs$bf,
    log_BF_bs      = log(bf_bs$bf),
    method         = "bridge_sampling"
  )
}, error = function(e) {
  tibble::tibble(
    hypothesis = "overall_vs_intercepts_only",
    BF_bs      = NA_real_,
    log_BF_bs  = NA_real_,
    method     = "bridge_sampling",
    error      = conditionMessage(e)
  )
})

calibration <- dplyr::left_join(sd_bfs, bridge_bfs, by = "hypothesis") |>
  dplyr::mutate(cell_id = cell_id)

tmp <- paste0(out_path, ".tmp")
readr::write_csv(calibration, tmp)
file.rename(tmp, out_path)
message("[calibration] Written: ", out_path)

`%||%` <- function(a, b) if (!is.null(a)) a else b
