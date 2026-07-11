#!/usr/bin/env Rscript
# scripts/generate_databased_v2_decision_grid.R
# ---------------------------------------------------------------------------
# Decision arm for the collaborator-imposed design cap (2026-07-11): data
# collection is limited to 80 participants per language at the current verb
# counts, so the decision-relevant design points are N = 70, 80 and 100
# (80 was never simulated; 70 and 100 bracket it and existed only at 24
# replicates). Assurance mode, single language, same engine and criteria as
# the main v2 grid, with two upgrades that answer the collaborator review:
#
#   * 120 replicates per cell (Monte-Carlo SE about 3.7 points near 80%
#     power, against about 10 points at the original 24 replicates);
#   * draw_index sampled RANDOMLY over the full pilot posterior (8,000
#     draws) instead of the first-k block, removing the representativeness
#     concern entirely.
#
# Cells write to a dedicated output directory so the decision arm can be
# reported separately from (or pooled with) the original block-draw cells.
#
# Seed base 940000: clear of all earlier grids (6e5, 7e5, 9e5, 91e4,
# 929000-932290).
#
# Usage: Rscript generate_databased_v2_decision_grid.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

OUT    <- "config/design_grid_databased_v2_decision.csv"
REPS   <- 120L
NDRAWS <- 8000L

langs <- tibble::tibble(
  language           = c("English", "Turkish", "Norwegian"),
  has_pseudo_passive = c(TRUE, TRUE, FALSE)
)
N_SWEEP <- c(70L, 80L, 100L)

set.seed(20260711)

grid <- langs |>
  tidyr::crossing(n_participants = N_SWEEP, .rep = seq_len(REPS)) |>
  dplyr::mutate(
    mode           = "assurance",
    prior_source   = "pilot",
    model_level    = "L5_correlated_maximal",
    prior_regime   = "primary",
    threshold_mode = "broad",
    effect_mult    = 1.0,
    draw_index     = sample.int(NDRAWS, dplyr::n(), replace = TRUE),
    iter = 3000L, warmup = 1000L, chains = 4L,
    seed = as.integer(940000L + dplyr::row_number() - 1L)
  ) |>
  dplyr::select(-.rep)

COL_ORDER <- c("language", "n_participants", "mode", "draw_index", "prior_source",
               "model_level", "prior_regime", "threshold_mode", "effect_mult",
               "has_pseudo_passive", "iter", "warmup", "chains", "seed")
grid <- dplyr::select(grid, dplyr::all_of(COL_ORDER))
readr::write_csv(grid, OUT)
message(sprintf("[decision grid] %d cells (3 langs x %d N x assurance x %d reps) -> %s",
                nrow(grid), length(N_SWEEP), REPS, OUT))
