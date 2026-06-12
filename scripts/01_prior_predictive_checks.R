#!/usr/bin/env Rscript
# scripts/01_prior_predictive_checks.R
# Run prior predictive simulations for each prior regime × threshold mode × language.
# Saves plots and summary statistics to outputs/prior_predictive/.
# Usage: Rscript scripts/01_prior_predictive_checks.R --language English --prior primary --threshold broad

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(ggplot2)
  library(yaml)
  library(optparse)
})

source("R/02_preprocess_factors.R")
source("R/03_define_priors.R")
source("R/04_model_formulas.R")

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
option_list <- list(
  optparse::make_option("--language",  default = "English",  help = "Language label"),
  optparse::make_option("--prior",     default = "primary", help = "Prior regime name"),
  optparse::make_option("--threshold", default = "broad",    help = "Threshold mode"),
  optparse::make_option("--config",    default = "config/analysis_config.yaml"),
  optparse::make_option("--outdir",    default = "outputs/prior_predictive"),
  optparse::make_option("--seed",      default = 42L, type = "integer"),
  optparse::make_option("--overwrite", action = "store_true", default = FALSE)
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)
has_pp <- isTRUE(cfg$languages[[opt$language]]$has_pseudo_passive)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

out_base <- file.path(opt$outdir,
  paste0(opt$language, "_", opt$prior, "_", opt$threshold))
out_rds  <- paste0(out_base, "_ppc.rds")

if (file.exists(out_rds) && !opt$overwrite) {
  message("[ppc] Skipping existing: ", out_rds)
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Build a tiny synthetic dataset for prior predictive simulation only
# ---------------------------------------------------------------------------
set.seed(opt$seed)
n_pp <- 10L; n_vb <- 10L
s_types <- if (has_pp) c("Passive","Active","Pseudo_Passive") else c("Passive","Active")

syn_data <- tidyr::crossing(
  Participant    = paste0("P", seq_len(n_pp)),
  Verb           = paste0("V", seq_len(n_vb)),
  S_Type         = s_types,
  Semantics_scaled = c(-0.5, 0, 0.5)
) |>
  dplyr::mutate(Response = 1L) |>   # placeholder; ignored for prior predictive
  preprocess_data(has_pseudo_passive = has_pp)

# ---------------------------------------------------------------------------
# Build model and prior
# ---------------------------------------------------------------------------
ladder  <- build_model_ladder(has_pseudo_passive = has_pp)
formula <- ladder[["L0_random_intercepts_only"]]  # cheapest for ppc

threshold_params <- NULL
if (opt$threshold == "ceiling_calibrated") {
  pilot_path <- cfg$pilot_data_path %||% NULL
  if (!is.null(pilot_path) && file.exists(pilot_path)) {
    pilot_df <- readr::read_csv(pilot_path, show_col_types = FALSE) |>
      dplyr::filter(Language == opt$language)
    threshold_params <- compute_ceiling_calibrated_thresholds(pilot_df, opt$language)
  } else {
    message("[ppc] No pilot data found; falling back to broad thresholds.")
    opt$threshold <- "broad"
  }
}

prior_obj <- build_brms_prior(opt$prior, opt$threshold, threshold_params, has_pp)

# ---------------------------------------------------------------------------
# Fit prior predictive model (sample_prior = "only")
# ---------------------------------------------------------------------------
message("[ppc] Fitting prior predictive for ", opt$language, " / ", opt$prior,
        " / ", opt$threshold)

fit_ppc <- brms::brm(
  formula,
  data         = syn_data,
  prior        = prior_obj,
  backend      = "cmdstanr",
  sample_prior = "only",
  iter         = 2000,
  warmup       = 1000,
  chains       = 2,
  cores        = 2,
  seed         = opt$seed,
  silent       = 2
)

# ---------------------------------------------------------------------------
# Plot and save
# ---------------------------------------------------------------------------
ppc_plot <- brms::pp_check(fit_ppc, ndraws = 100, type = "hist") +
  ggplot2::labs(
    title    = paste0("Prior Predictive Check: ", opt$language,
                      " | ", opt$prior, " | ", opt$threshold),
    subtitle = "Responses must be plausible 1-7 integers; ceiling effects visible if present"
  ) +
  ggplot2::theme_bw()

ggsave_path <- paste0(out_base, "_ppc_plot.pdf")
ggplot2::ggsave(ggsave_path, ppc_plot, width = 8, height = 5)
message("[ppc] Plot saved: ", ggsave_path)

tmp_rds <- paste0(out_rds, ".tmp")
saveRDS(fit_ppc, tmp_rds)
file.rename(tmp_rds, out_rds)
message("[ppc] RDS saved: ", out_rds)

`%||%` <- function(a, b) if (!is.null(a)) a else b
