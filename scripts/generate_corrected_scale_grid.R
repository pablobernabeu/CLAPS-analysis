#!/usr/bin/env Rscript
# scripts/generate_corrected_scale_grid.R
# ---------------------------------------------------------------------------
# Prior-SCALE-corrected single-language power BFDA grid.
#
# Same design-generating assumptions as scripts/generate_corrected_power_grid.R
# (per-verb affectedness fix; assumed effects beta_semantics = 0.8, active
# interaction = -0.5, pseudo interaction = 0.2). The ONLY difference at run time
# is that the cells are executed with scripts/04_design_analysis_cell_gelman.R,
# which sources R/06_simulate_design_gelman.R so the ANALYSED focal predictor is
# Gelman-scaled (SD 0.5, matching R/02_preprocess_factors.R and the real
# confirmatory analysis) rather than the raw U(-0.5,0.5) draw (SD 0.289). This
# corrects the prior-scale mismatch surfaced by the 2026-06-20 audit, under which
# the base corrected grid overstates the focal-slope Bayes factor by ~sqrt(3).
#
# Focused confirmation: verb counts {40, 72} only (the design-relevant region the
# report uses; 12/20 are known to be insufficient and are dropped). Distinct seed
# base (7e5) so it never collides with the base corrected grid (6e5). Writes to
# config/design_grid_corrected_scale.csv; cells write to outputs/design_corrected_scale.
#
# Usage (from design_analysis/ root):
#   Rscript scripts/generate_corrected_scale_grid.R [--out config/design_grid_corrected_scale.csv] [--b 50]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(optparse) })

option_list <- list(
  optparse::make_option("--out", default = "config/design_grid_corrected_scale.csv"),
  optparse::make_option("--b",   default = 50L, type = "integer", help = "replicates per cell")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
B <- as.integer(opt$b)

# Assumed data-generating effect sizes (identical to the base corrected grid).
BSEM <- 0.8; BACT <- -0.5; BPSEUDO <- 0.2
ITER <- 3000L; WARMUP <- 1000L; CHAINS <- 4L

langs <- tibble::tibble(
  language           = c("English", "Turkish", "Norwegian"),
  has_pseudo_passive = c(TRUE, TRUE, FALSE),
  beta_pseudo        = c(BPSEUDO, BPSEUDO, 0)      # Norwegian has no pseudo-passive
)
N_SWEEP     <- c(50L, 60L, 70L, 80L, 90L, 100L, 120L)
NVERB_SWEEP <- c(40L, 72L)                          # design-relevant region only

conditions <- langs |>
  tidyr::crossing(n_participants = N_SWEEP, n_verbs = NVERB_SWEEP) |>
  dplyr::mutate(
    model_level             = "L5_correlated_maximal",
    prior_regime            = "primary",
    threshold_mode          = "broad",
    n_items_per_cell        = 1L,
    beta_semantics          = BSEM,
    beta_active_interaction = BACT,
    beta_pseudo_interaction = beta_pseudo,
    iter = ITER, warmup = WARMUP, chains = CHAINS,
    gender_spec             = "none",
    include_gender          = FALSE,
    beta_gender             = 0.3,
    beta_gender_sem_passive = 0.15
  ) |>
  dplyr::select(-beta_pseudo) |>
  dplyr::mutate(.cond = dplyr::row_number())

grid <- tidyr::crossing(conditions, .rep = seq_len(B)) |>
  dplyr::mutate(seed = as.integer(700000L + (.cond - 1L) * B + (.rep - 1L))) |>   # base 7e5: collision-free vs 6e5/8e5
  dplyr::select(-.cond, -.rep)

COL_ORDER <- c(
  "language", "model_level", "prior_regime", "threshold_mode",
  "n_participants", "n_verbs", "n_items_per_cell",
  "beta_semantics", "beta_active_interaction", "beta_pseudo_interaction",
  "has_pseudo_passive", "iter", "warmup", "chains", "seed",
  "gender_spec", "include_gender", "beta_gender", "beta_gender_sem_passive"
)
grid <- dplyr::select(grid, dplyr::all_of(COL_ORDER))
readr::write_csv(grid, opt$out)
message(sprintf("[corrected-scale grid] %d cells (%d langs x %d N x %d n_verbs x %d reps) -> %s",
                nrow(grid), nrow(langs), length(N_SWEEP), length(NVERB_SWEEP), B, opt$out))
