#!/usr/bin/env Rscript
# scripts/generate_pooled_v2_grid.R
# ---------------------------------------------------------------------------
# Grid for the pooled cross-language pilot-grounded design analysis (the
# recommended path in reports/preliminary_sample_size_analysis.qmd). Simulates
# all three languages per replicate from their per-language v2 pilot DGPs and
# fits the pre-existing cross-language analysis model L5_cross_maximal.
#
# N is PER LANGUAGE (total sample = 3N). Assurance mode only for this first
# pass; the safeguard variant can follow once the assurance surface is known.
# 30 replicates per cell gives a Monte-Carlo SE of about 0.07-0.09 near 0.80,
# which is enough to read the trend across the five N values.
#
# Draw indices are sampled RANDOMLY over the full pilot posterior (8,000 draws)
# per replicate and per language, rather than reusing the first-k draws, which
# removes the contiguous-draw-block limitation noted in the report.
#
# Seed base 930000 with spacing 10 per row: rows use seed, seed+1, seed+2 for
# the three language simulations (see simulate_pooled_from_pilots), so spacing
# 10 keeps every RNG stream in the project unique (other grids use bases
# 6e5, 7e5, 9e5 and 91e4, all clear of 930000-931500).
#
# Also writes a 1-row smoke-test grid (tiny N, short chains) used to validate
# the pipeline end to end on the cluster before the real array is submitted.
#
# Usage: Rscript scripts/generate_pooled_v2_grid.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

OUT      <- "config/design_grid_pooled_v2.csv"
OUT_TEST <- "config/design_grid_pooled_v2_TEST.csv"
REPS     <- 30L
N_SWEEP  <- c(50L, 70L, 100L, 130L, 150L)
NDRAWS   <- 8000L   # draws in each pilot_dgp_v2_pilot_<lang>.rds posterior

set.seed(20260709)

grid <- tidyr::crossing(n_participants = N_SWEEP, .rep = seq_len(REPS)) |>
  dplyr::mutate(
    mode                 = "assurance",
    prior_source         = "pilot",
    model_level          = "L5_cross_maximal",
    prior_regime         = "primary",
    threshold_mode       = "broad",
    draw_index_english   = sample.int(NDRAWS, dplyr::n(), replace = TRUE),
    draw_index_turkish   = sample.int(NDRAWS, dplyr::n(), replace = TRUE),
    draw_index_norwegian = sample.int(NDRAWS, dplyr::n(), replace = TRUE),
    iter = 3000L, warmup = 1000L, chains = 4L,
    seed = as.integer(930000L + (dplyr::row_number() - 1L) * 10L)
  ) |>
  dplyr::select(-.rep)

readr::write_csv(grid, OUT)
message(sprintf("[pooled v2 grid] %d cells (%d N x %d reps, assurance) -> %s",
                nrow(grid), length(N_SWEEP), REPS, OUT))

# Smoke test: one tiny, fast cell. Distinct seed range (929000) and its own
# output dir at submit time keep it clear of production cells.
test <- grid[1, ] |>
  dplyr::mutate(n_participants = 12L, iter = 800L, warmup = 400L, chains = 2L,
                seed = 929000L)
readr::write_csv(test, OUT_TEST)
message(sprintf("[pooled v2 grid] smoke-test row -> %s", OUT_TEST))
