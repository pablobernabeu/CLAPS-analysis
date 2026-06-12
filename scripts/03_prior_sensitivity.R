#!/usr/bin/env Rscript
# scripts/03_prior_sensitivity.R
# Run prior sensitivity analysis across all prior regimes × threshold modes
# for a given language, using the pilot-fitted model as the base.
# Saves BF and diagnostic results per cell.
# Usage: Rscript scripts/03_prior_sensitivity.R --language English --model_level L4_uncorrelated_maximal

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(yaml)
  library(optparse)
  library(readr)
  library(purrr)
})

source("R/01_read_validate_data.R")
source("R/02_preprocess_factors.R")
source("R/03_define_priors.R")
source("R/04_model_formulas.R")
source("R/05_hypothesis_tests.R")
source("R/07_extract_diagnostics.R")

option_list <- list(
  optparse::make_option("--language",    default = "English"),
  optparse::make_option("--model_level", default = "L4_uncorrelated_maximal"),
  optparse::make_option("--include_gender", action = "store_true", default = FALSE,
    help = "Fit the gender model variation (adds the Gender covariate)."),
  optparse::make_option("--semantics_source", default = NULL,
    help = "Source column for Semantics (e.g. affectedness_scores_agent)."),
  optparse::make_option("--config",      default = "config/analysis_config.yaml"),
  optparse::make_option("--outdir",      default = "outputs/prior_sensitivity"),
  optparse::make_option("--seed",        default = 2025L, type = "integer"),
  optparse::make_option("--overwrite",   action  = "store_true", default = FALSE)
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
variant_tag <- if (isTRUE(opt$include_gender)) "_gender" else ""

# The gender model variation requires the agent/gender-specific affectedness
# (affectedness_scores_agent), not the standard whole-event Semantics.
if (isTRUE(opt$include_gender) && is.null(opt$semantics_source)) {
  opt$semantics_source <- "affectedness_scores_agent"
  message("[sensitivity] Gender variation: sourcing Semantics from 'affectedness_scores_agent'.")
}

cfg    <- yaml::read_yaml(opt$config)
has_pp <- isTRUE(cfg$languages[[opt$language]]$has_pseudo_passive)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Load pilot data
pilot_data_path <- cfg$pilot_data_path
raw   <- read_raw_data(pilot_data_path)
valid <- validate_raw_data(raw, pilot_data_path)
pilot <- split_pilot_confirmatory(valid)$pilot
df    <- dplyr::filter(pilot, Language == opt$language) |>
  preprocess_data(has_pseudo_passive = has_pp,
                  semantics_source = opt$semantics_source,
                  include_gender   = opt$include_gender)

# Build sensitivity grid
sens_grid <- prior_sensitivity_grid()
ladder    <- build_model_ladder(has_pseudo_passive = has_pp,
                                include_gender = opt$include_gender)
formula   <- ladder[[opt$model_level]]
if (is.null(formula)) stop("[sensitivity] Unknown model_level: ", opt$model_level)

samp_args <- production_sampling(seed = opt$seed,
  cores = as.integer(Sys.getenv("STAN_NUM_THREADS", "4")))
ctrl_args <- production_control()

# Run each prior × threshold cell
results <- purrr::pmap_dfr(sens_grid, function(regime_name, threshold_mode) {
  cell_label <- paste0(opt$language, "_", opt$model_level, "_",
                       regime_name, "_", threshold_mode, variant_tag)
  out_rds  <- file.path(opt$outdir, paste0(cell_label, ".rds"))
  diag_csv <- file.path(opt$outdir, paste0(cell_label, "_diag.csv"))
  bf_csv   <- file.path(opt$outdir, paste0(cell_label, "_bf.csv"))

  if (all(file.exists(c(diag_csv, bf_csv))) && !opt$overwrite) {
    message("[sensitivity] Skipping existing: ", cell_label)
    diag <- readr::read_csv(diag_csv, show_col_types = FALSE)
    bf   <- readr::read_csv(bf_csv,   show_col_types = FALSE)
    return(dplyr::bind_cols(
      tibble::tibble(cell_label = cell_label, regime_name = regime_name,
                     threshold_mode = threshold_mode, status = "loaded"),
      dplyr::select(bf, hypothesis, BF_10, bf_category),
      dplyr::select(diag, max_rhat, n_divergent, convergence_ok)
    ))
  }

  message("[sensitivity] Fitting: ", cell_label)

  threshold_params <- NULL
  if (threshold_mode == "ceiling_calibrated") {
    threshold_params <- compute_ceiling_calibrated_thresholds(df, opt$language)
  }
  prior_obj <- build_brms_prior(regime_name, threshold_mode, threshold_params, has_pp)

  fit_result <- tryCatch({
    fit <- do.call(brms::brm, c(
      list(formula = formula, data = df, prior = prior_obj,
           backend = "cmdstanr", sample_prior = "yes", control = ctrl_args),
      samp_args
    ))
    list(fit = fit, error = NULL)
  }, error = function(e) list(fit = NULL, error = conditionMessage(e)))

  if (!is.null(fit_result$error)) {
    return(tibble::tibble(cell_label = cell_label, status = "error",
                          error_msg = fit_result$error))
  }

  fit  <- fit_result$fit
  diag <- extract_convergence_diagnostics(fit)
  bf   <- compute_all_bf(fit, has_pseudo_passive = has_pp)

  # Atomic saves
  tmp_rds <- paste0(out_rds, ".tmp")
  saveRDS(fit, tmp_rds); file.rename(tmp_rds, out_rds)
  save_diagnostics_csv(diag, diag_csv)
  readr::write_csv(bf, bf_csv)

  dplyr::bind_cols(
    tibble::tibble(cell_label = cell_label, regime_name = regime_name,
                   threshold_mode = threshold_mode, status = "success"),
    dplyr::select(bf, hypothesis, BF_10, bf_category),
    dplyr::select(diag, max_rhat, n_divergent, convergence_ok)
  )
})

summary_path <- file.path(opt$outdir,
  paste0(opt$language, "_", opt$model_level, variant_tag, "_sensitivity_summary.csv"))
readr::write_csv(results, summary_path)
message("[sensitivity] Summary written: ", summary_path)
