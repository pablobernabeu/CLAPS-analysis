#!/usr/bin/env Rscript
# scripts/generate_databased_v2_extendedN_grid.R
# ---------------------------------------------------------------------------
# Extended-N follow-up to the assurance result: English and Turkish did not
# reach 80% binding power by N=130 in the original v2 sweep, but the pilot's
# own posterior overwhelmingly supports the theoretically predicted direction
# for both focal effects (asymptotic ceiling 0.94-1.0 for both languages), so
# the shortfall is a matter of how much N is needed, not a structural ceiling.
# Norwegian already clears 80-90% by N=100-130 and is not extended here.
# Assurance mode only (safeguard is being recomputed separately after the
# sign-awareness fix). More replicates than the original sweep (40 vs 24) to
# sharpen the estimate given these cells are now pivotal for the recommendation.
#
# Usage: Rscript generate_databased_v2_extendedN_grid.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr) })

OUT  <- "config/design_grid_databased_v2_extendedN.csv"
REPS <- 40L

langs <- tibble::tibble(
  language           = c("English", "Turkish"),
  has_pseudo_passive  = c(TRUE, TRUE)
)
N_SWEEP <- c(150L, 200L, 250L, 300L, 400L, 500L)

conditions <- langs |>
  tidyr::crossing(n_participants = N_SWEEP, mode = "assurance") |>
  dplyr::mutate(
    prior_source   = "pilot",
    model_level    = "L5_correlated_maximal",
    prior_regime   = "primary",
    threshold_mode = "broad",
    effect_mult    = 1.0,
    iter = 3000L, warmup = 1000L, chains = 4L,
    .cond = dplyr::row_number()
  )

grid <- tidyr::crossing(conditions, .rep = seq_len(REPS)) |>
  dplyr::mutate(
    draw_index = .rep,
    seed       = as.integer(910000L + (.cond - 1L) * REPS + (.rep - 1L))   # base 91e4: clear of the 9e5 base grid
  ) |>
  dplyr::select(-.cond, -.rep)

COL_ORDER <- c("language", "n_participants", "mode", "draw_index", "prior_source",
               "model_level", "prior_regime", "threshold_mode", "effect_mult",
               "has_pseudo_passive", "iter", "warmup", "chains", "seed")
grid <- dplyr::select(grid, dplyr::all_of(COL_ORDER))
readr::write_csv(grid, OUT)
message(sprintf("[extended-N v2 grid] %d cells (%d langs x %d N x assurance x %d reps) -> %s",
                nrow(grid), nrow(langs), length(N_SWEEP), REPS, OUT))
