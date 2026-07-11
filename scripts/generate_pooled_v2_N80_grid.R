#!/usr/bin/env Rscript
# scripts/generate_pooled_v2_N80_grid.R
# ---------------------------------------------------------------------------
# Pooled cross-language cells at the collaborator-imposed cap of 80
# participants per language (2026-07-11). The main pooled grid swept
# {50, 70, 100, 130, 150}; 80 is now the actual design point, so it gets its
# own cells rather than an interpolation. Same engine and settings as the
# main pooled grid (randomised draw indices, assurance mode).
#
# Seed base 932000 with spacing 10 (languages use +0/1/2): clear of the main
# pooled grid (930000-931490) and every other grid.
#
# Usage: Rscript generate_pooled_v2_N80_grid.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(readr) })

OUT    <- "config/design_grid_pooled_v2_N80.csv"
REPS   <- 30L
NDRAWS <- 8000L

set.seed(20260711)

grid <- tibble::tibble(.rep = seq_len(REPS)) |>
  dplyr::mutate(
    n_participants       = 80L,
    mode                 = "assurance",
    prior_source         = "pilot",
    model_level          = "L5_cross_maximal",
    prior_regime         = "primary",
    threshold_mode       = "broad",
    draw_index_english   = sample.int(NDRAWS, dplyr::n(), replace = TRUE),
    draw_index_turkish   = sample.int(NDRAWS, dplyr::n(), replace = TRUE),
    draw_index_norwegian = sample.int(NDRAWS, dplyr::n(), replace = TRUE),
    iter = 3000L, warmup = 1000L, chains = 4L,
    seed = as.integer(932000L + (.rep - 1L) * 10L)
  ) |>
  dplyr::select(-.rep) |>
  dplyr::relocate(n_participants)

readr::write_csv(grid, OUT)
message(sprintf("[pooled N80 grid] %d cells (N = 80 per language, assurance) -> %s",
                nrow(grid), OUT))
