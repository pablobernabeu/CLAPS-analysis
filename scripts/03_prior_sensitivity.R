#!/usr/bin/env Rscript
# scripts/03_prior_sensitivity.R
#
# Purpose
#   Refit the same pilot data under every prior regime crossed with every threshold
#   mode, and record how far the Bayes factors move. The question is whether a
#   conclusion is a property of the data or an artefact of the prior. It matters
#   more here than in an estimation analysis, because a Savage-Dickey Bayes factor
#   depends on the prior's width directly, not just on its centre.
#
# What varies and what does not
#   The data, the model formula and the sampler settings are held fixed across all
#   eight cells. Only the prior changes. Any difference in the resulting Bayes
#   factors is therefore attributable to the prior alone.
#
# Inputs
#   --language          Which language's pilot data to use. One language per run.
#   --model_level       Ladder level to fit. Defaults to L4 rather than the maximal
#                       L5, because this script fits eight models and L4 is the
#                       level that reliably converges on pilot-sized data; use the
#                       language's maximal feasible level once that is known.
#   --include_gender    Fit the gender model variation. Implies a different
#                       affectedness source; see below.
#   --semantics_source  Column to take affectedness from.
#   --config            config/analysis_config.yaml.
#   --outdir            Destination for the per-cell files and the summary.
#   --seed, --overwrite
#
# Outputs, per cell (8 cells = 4 regimes x 2 threshold modes)
#   <label>.rds        The fitted model.
#   <label>_diag.csv   Convergence diagnostics.
#   <label>_bf.csv     Bayes factors.
#   plus one <language>_<model_level>[_gender]_sensitivity_summary.csv across cells.
#
# Usage
#   Rscript scripts/03_prior_sensitivity.R --language English \
#     --model_level L4_uncorrelated_maximal
#
# Cost
#   Eight full brms fits, run sequentially in one process. Submit via
#   hpc/submit_prior_sensitivity_array.sh rather than running it on a login node.
#
# Reading the result
#   Compare bf_category across regimes rather than raw BF values. A Bayes factor
#   moving from 40 to 120 does not change any conclusion, whereas one crossing the
#   preregistered threshold of 10 does. The literature_centred regime is expected
#   to give the largest Bayes factors, since it encodes the predicted direction;
#   that is why it is a sensitivity check and never the primary result.

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

  # Resume support: a cell is skipped only if BOTH its CSVs exist, since the two are
  # written separately and a job killed between them would otherwise be treated as
  # complete. The saved summaries are re-read so the returned table is identical
  # whether a cell was fitted now or previously, and `status` records which.
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

  # Threshold priors are derived from data only in ceiling_calibrated mode; broad
  # mode uses the generic Student-t and needs no parameters.
  #
  # CAVEAT: the thresholds are computed from `df`, the same data the model is then
  # fitted to. compute_ceiling_calibrated_thresholds() documents that its input
  # should be an INDEPENDENT pilot sample. Using the analysis data makes this cell
  # mildly optimistic, since the threshold prior is centred on the very response
  # distribution it will be tested against. It is acceptable here because this
  # script is a sensitivity check on the pilot rather than a confirmatory analysis,
  # and because the thresholds are nuisance parameters that no hypothesis concerns.
  # For the confirmatory analysis the pilot and confirmatory samples must be split,
  # which is what split_pilot_confirmatory() is for.
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
