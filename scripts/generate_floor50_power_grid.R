#!/usr/bin/env Rscript
# scripts/generate_floor50_power_grid.R
# ---------------------------------------------------------------------------
# Floor-of-50-verbs power BFDA grid.
#
# The design's minimum viable materials are 50 verbs per language. This grid
# fixes n_verbs = 50 and sweeps the participant count, so we can read off the
# sample size needed for 90% power at the floor verb count. It mirrors
# generate_corrected_power_grid.R in every other respect (same per-verb
# affectedness fix, same L5 single-language model, same assumed effects, same
# primary priors deployed on ARC) so the cells slot straight into the
# design_corrected surface alongside the 12/20/40/72-verb cells.
#
# Cells write to outputs/design_corrected (n_verbs = 50 in the file name, and a
# distinct seed base, so there is no collision with the existing run).
#
# Usage (from design_analysis/ root):
#   Rscript scripts/generate_floor50_power_grid.R [--out config/design_grid_floor50.csv] [--b 50]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(optparse) })

option_list <- list(
  optparse::make_option("--out", default = "config/design_grid_floor50.csv"),
  optparse::make_option("--b",   default = 50L, type = "integer", help = "replicates per cell")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
B <- as.integer(opt$b)

# Assumed data-generating effect sizes (match the corrected power grid).
BSEM <- 0.8; BACT <- -0.5; BPSEUDO <- 0.2
ITER <- 3000L; WARMUP <- 1000L; CHAINS <- 4L

langs <- tibble::tibble(
  language           = c("English", "Turkish", "Norwegian"),
  has_pseudo_passive = c(TRUE, TRUE, FALSE),
  beta_pseudo        = c(BPSEUDO, BPSEUDO, 0)      # Norwegian has no pseudo-passive
)
N_SWEEP     <- c(50L, 60L, 70L, 80L, 90L, 100L, 120L)   # brackets the 90% target at 50 verbs
NVERB_SWEEP <- c(50L)                                    # the floor verb count

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
  dplyr::mutate(seed = as.integer(700000L + (.cond - 1L) * B + (.rep - 1L))) |>   # base 7e5: clear of the 6e5 corrected grid
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
message(sprintf("[floor50 grid] %d cells (%d langs x %d N x 50 verbs x %d reps) -> %s",
                nrow(grid), nrow(langs), length(N_SWEEP), B, opt$out))
