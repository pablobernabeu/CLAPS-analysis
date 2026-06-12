#!/usr/bin/env Rscript
# scripts/02_fit_pilot_models.R
# Fit pilot data models across the model ladder for a given language and prior.
# This is a CALIBRATION-ONLY analysis. Pilot data must never enter confirmatory inference.
# Usage: Rscript scripts/02_fit_pilot_models.R --language English --prior proposal \
#                --threshold broad [--start_level L5_correlated_maximal] [--overwrite]

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(yaml)
  library(optparse)
  library(readr)
})

source("R/01_read_validate_data.R")
source("R/02_preprocess_factors.R")
source("R/03_define_priors.R")
source("R/04_model_formulas.R")
source("R/07_extract_diagnostics.R")
source("R/09_model_ladder.R")

option_list <- list(
  optparse::make_option("--language",    default = "English"),
  optparse::make_option("--prior",       default = "proposal"),
  optparse::make_option("--threshold",   default = "broad"),
  optparse::make_option("--start_level", default = "L5_correlated_maximal"),
  optparse::make_option("--include_gender", action = "store_true", default = FALSE,
    help = "Fit the gender model variation (adds the Gender covariate)."),
  optparse::make_option("--semantics_source", default = NULL,
    help = "Source column for Semantics (e.g. affectedness_scores_agent). Default keeps the existing column."),
  optparse::make_option("--config",      default = "config/analysis_config.yaml"),
  optparse::make_option("--outdir",      default = "outputs/pilot_models"),
  optparse::make_option("--seed",        default = 2025L, type = "integer"),
  optparse::make_option("--overwrite",   action  = "store_true", default = FALSE)
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# The gender model variation requires the agent/gender-specific affectedness
# (affectedness_scores_agent), NOT the standard whole-event Semantics used in
# the previous studies. Default the source accordingly unless overridden.
if (isTRUE(opt$include_gender) && is.null(opt$semantics_source)) {
  opt$semantics_source <- "affectedness_scores_agent"
  message("[fit_pilot] Gender variation: sourcing Semantics from 'affectedness_scores_agent'.")
}

cfg     <- yaml::read_yaml(opt$config)
has_pp  <- isTRUE(cfg$languages[[opt$language]]$has_pseudo_passive)

# Load and validate pilot data
pilot_data_path <- cfg$pilot_data_path
if (is.null(pilot_data_path) || !file.exists(pilot_data_path)) {
  stop("[fit_pilot] Pilot data not found at: ", pilot_data_path)
}

message("[fit_pilot] Reading pilot data: ", pilot_data_path)
raw    <- read_raw_data(pilot_data_path)
valid  <- validate_raw_data(raw, pilot_data_path)
split  <- split_pilot_confirmatory(valid)
pilot  <- split$pilot

if (nrow(pilot) == 0) stop("[fit_pilot] No pilot rows found after split.")

pilot_lang <- dplyr::filter(pilot, Language == opt$language)
if (nrow(pilot_lang) == 0) {
  stop("[fit_pilot] No pilot rows for language: ", opt$language)
}

df <- preprocess_data(pilot_lang, has_pseudo_passive = has_pp,
                      semantics_source = opt$semantics_source,
                      include_gender   = opt$include_gender)

# Tag outputs so the gender variation does not overwrite baseline fits/logs.
variant_tag <- if (isTRUE(opt$include_gender)) "_gender" else ""

# Threshold calibration
threshold_params <- NULL
if (opt$threshold == "ceiling_calibrated") {
  threshold_params <- compute_ceiling_calibrated_thresholds(df, opt$language)
}

prior_obj <- build_brms_prior(opt$prior, opt$threshold, threshold_params, has_pp)

samp_args <- production_sampling(seed = opt$seed,
  cores = as.integer(Sys.getenv("STAN_NUM_THREADS", "4")))
ctrl_args <- production_control()

message("[fit_pilot] Starting model ladder: ", opt$language, " / ",
        opt$prior, " / ", opt$threshold, " from ", opt$start_level)

result <- run_model_ladder(
  df                = df,
  prior_obj         = prior_obj,
  sampling_args     = samp_args,
  control_args      = ctrl_args,
  has_pseudo_passive = has_pp,
  start_level       = opt$start_level,
  out_dir           = opt$outdir,
  language          = opt$language,
  prior_label       = paste0(opt$prior, "_", opt$threshold, variant_tag),
  include_gender    = opt$include_gender,
  overwrite         = opt$overwrite
)

message("[fit_pilot] Selected model: ", result$selected_model %||% "NONE")

# Write ladder log summary
log_path <- file.path(opt$outdir, paste0(
  opt$language, "_", opt$prior, "_", opt$threshold, variant_tag, "_ladder_log.csv"
))
ladder_log_df <- purrr::imap_dfr(result$ladder_log, function(entry, lv) {
  tibble::tibble(
    model_level        = lv,
    status             = entry$status %||% NA_character_,
    convergence_status = entry$convergence_status %||% NA_character_,
    runtime_sec        = entry$runtime_sec %||% NA_real_,
    error_msg          = substr(entry$error %||% "", 1, 300)
  )
})
readr::write_csv(ladder_log_df, log_path)
message("[fit_pilot] Ladder log: ", log_path)

`%||%` <- function(a, b) if (!is.null(a)) a else b
