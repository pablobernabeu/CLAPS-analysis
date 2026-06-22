#!/usr/bin/env Rscript
# scripts/generate_databased_v2_grid.R
# ---------------------------------------------------------------------------
# Grid for the amended (assurance / safeguard) data-grounded design analysis.
# Modes (R/10_simulate_from_pilot_v2.R):
#   * assurance : draw_index = replicate, so each replicate uses a distinct pilot
#                 posterior draw of the focal effects; P(BF>=10) over replicates
#                 integrates the pilot's effect-size uncertainty (the A&L fix).
#   * safeguard : focal slopes at their lower posterior bound (principled
#                 replacement for the ad hoc effect_mult = 0.6).
# prior_source defaults to "pilot" (zero-centred DGP). A "blend" source
# (literature_centred DGP) can be added once the fit-comparison shows it moves
# the DGP at the real verb count.
# Refit/BF still use the zero-centred "primary" prior, so the test is unbiased.
#
# Usage: Rscript scripts/generate_databased_v2_grid.R [--out ...] [--reps 24] [--prior_source pilot]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(optparse) })

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--out",          default = "config/design_grid_databased_v2.csv"),
  optparse::make_option("--reps",         default = 24L, type = "integer"),
  optparse::make_option("--prior_source", default = "pilot")
)))
REPS <- as.integer(opt$reps)

langs <- tibble::tibble(
  language           = c("English", "Turkish", "Norwegian"),
  has_pseudo_passive = c(TRUE, TRUE, FALSE)
)
N_SWEEP <- c(30L, 50L, 70L, 100L, 130L)
MODES   <- c("assurance", "safeguard")

conditions <- langs |>
  tidyr::crossing(n_participants = N_SWEEP, mode = MODES) |>
  dplyr::mutate(
    prior_source   = opt$prior_source,
    model_level    = "L5_correlated_maximal",
    prior_regime   = "primary",      # zero-centred refit prior (unbiased BF)
    threshold_mode = "broad",
    effect_mult    = 1.0,
    iter = 3000L, warmup = 1000L, chains = 4L,
    .cond = dplyr::row_number()
  )

grid <- tidyr::crossing(conditions, .rep = seq_len(REPS)) |>
  dplyr::mutate(
    draw_index = .rep,                                   # assurance: distinct posterior draw per replicate
    seed       = as.integer(900000L + (.cond - 1L) * REPS + (.rep - 1L))
  ) |>
  dplyr::select(-.cond, -.rep)

COL_ORDER <- c("language", "n_participants", "mode", "draw_index", "prior_source",
               "model_level", "prior_regime", "threshold_mode", "effect_mult",
               "has_pseudo_passive", "iter", "warmup", "chains", "seed")
grid <- dplyr::select(grid, dplyr::all_of(COL_ORDER))
readr::write_csv(grid, opt$out)
message(sprintf("[databased v2 grid] %d cells (%d langs x %d N x %d modes x %d reps) source=%s -> %s",
                nrow(grid), nrow(langs), length(N_SWEEP), length(MODES), REPS, opt$prior_source, opt$out))
