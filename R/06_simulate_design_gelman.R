# R/06_simulate_design_gelman.R
# ---------------------------------------------------------------------------
# Prior-scale-corrected variant of the literature-anchored design simulator.
#
# WHY THIS EXISTS (audit finding, 2026-06-20)
# The base simulator (R/06_simulate_design.R) stores the focal predictor as the
# raw uniform draw Semantics_scaled = sem ~ U(-0.5, 0.5), whose SD is 1/sqrt(12)
# ~= 0.289. The REAL analysis pipeline (R/02_preprocess_factors.R) instead uses
# Gelman scaling, Semantics_scaled = (x - mean(x)) / (2 * sd(x)), whose SD is 0.5;
# the data-grounded engine (R/10_simulate_from_pilot.R) does the same. The SAME
# weakly informative prior normal(0, 0.5) on the focal slope is applied in all
# cases (R/03_define_priors.R). Because the Savage-Dickey Bayes factor depends on
# the prior's width RELATIVE to the coefficient it constrains, holding the prior
# fixed while the predictor SD differs by a factor of sqrt(3) ~= 1.73 makes the
# base ("0.289") grid overstate the focal-slope Bayes factor by ~sqrt(3) relative
# to the real, Gelman-scaled confirmatory analysis. Per-SD equivalence of the
# TRUE effect is necessary but not sufficient for Bayes-factor power transfer when
# the prior is specified per raw predictor unit.
#
# THE FIX (pure reparametrisation)
# We re-express the ANALYSED predictor on the same Gelman scale the real pipeline
# uses, WITHOUT touching the data-generating linear predictor. We source the base
# simulators and override only the two simulate_* functions: each calls the
# original (identical RNG stream -> byte-identical simulated responses) and then
# overwrites Semantics_scaled with the Gelman transform of the (unchanged)
# Semantics column. Because no new RNG draw is introduced and eta is untouched,
# the simulated data are identical to the base run; only the scale on which the
# fixed prior acts changes, so the recovered focal slope (~0.46 on the Gelman
# scale for the assumed beta_semantics = 0.8) and its Bayes factor now match what
# the real analysis will see. The grid, betas, random-effect SDs, and thresholds
# are deliberately left unchanged so the comparison with the base corrected run
# isolates the prior-scale effect alone.
#
# Row-level scaling (over the assembled rows, grouped by Language for the
# multilanguage path) mirrors R/02_preprocess_factors.R::scale_semantics exactly;
# for the balanced simulated design it equals verb-level scaling up to the
# negligible (n-1) finite-sample correction.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
})

# Base simulators + run_design_cell + helpers. run_design_cell looks up
# simulate_claps_data / simulate_claps_data_multilanguage by name from the global
# environment at call time, so the overrides below take effect for it.
source("R/06_simulate_design.R")

# Gelman scaling, identical to R/02_preprocess_factors.R::scale_semantics.
.gelman_scale_semantics <- function(x) {
  (x - mean(x, na.rm = TRUE)) / (2 * stats::sd(x, na.rm = TRUE))
}

.orig_simulate_claps_data <- simulate_claps_data
simulate_claps_data <- function(...) {
  df <- .orig_simulate_claps_data(...)
  # Single-language data carry no Language column: scale over all rows.
  df$Semantics_scaled <- .gelman_scale_semantics(df$Semantics)
  df
}

.orig_simulate_claps_data_multilanguage <- simulate_claps_data_multilanguage
simulate_claps_data_multilanguage <- function(...) {
  df <- .orig_simulate_claps_data_multilanguage(...)
  # Scale within Language, mirroring the real pipeline's per-language centring.
  df <- df |>
    dplyr::group_by(Language) |>
    dplyr::mutate(Semantics_scaled = .gelman_scale_semantics(Semantics)) |>
    dplyr::ungroup()
  df
}

message("[06_simulate_design_gelman] Gelman-scaled predictor override active ",
        "(Semantics_scaled SD ~= 0.5; eta and responses unchanged).")
