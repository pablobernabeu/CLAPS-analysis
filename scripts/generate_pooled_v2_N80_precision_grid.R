#!/usr/bin/env Rscript
# scripts/generate_pooled_v2_N80_precision_grid.R
# ---------------------------------------------------------------------------
# Additional pooled cross-language cells at 80 participants per language, to
# settle whether the confirmatory test reaches a stated power.
#
# Why this grid exists
#   The pooled arm carries the recommendation and is the least replicated part
#   of the analysis: 135 simulated studies against 1,814 in the per-language
#   assurance grid and 2,019 in the literature-anchored cross-check. Combining
#   its cells at or below the 100-participant cap gives 84% [0.76, 0.90] for the
#   interaction at a Bayes factor of 6. The lower end of that interval sits below
#   80%, so the analysis does not establish that the confirmatory test reaches a
#   conventional target; it only estimates it.
#
# Why 80 participants only, rather than spreading across the swept range
#   The report combines cells across sample sizes on the grounds that pooled
#   power is flat in N. That is very likely true for a verb-limited design, but
#   it is established from the same thin cells whose imprecision is the problem,
#   so leaning on it to fix that imprecision is close to circular. Eighty is the
#   recommended design point, and estimating there directly needs no such
#   assumption. It is also the conservative choice: if flatness does hold, these
#   cells can still be combined with the others afterwards.
#
# How many
#   At the observed rate of about 0.84, an exact one-sided interval clears 80%
#   at roughly 380 replicates. Twenty-eight are usable at 80 participants now.
#   REPS defaults to 380 new cells rather than 352, allowing for the roughly one
#   cell in fifteen that has failed or been lost in this arm to date.
#
# Cost, from the runtimes this arm has actually recorded
#   The median completed pooled cell at 80 participants took about 65 core-hours,
#   so 380 cells is on the order of 25,000 core-hours. Submit split across both
#   project accounts, and expect days rather than hours.
#
# Same engine and settings as the existing pooled grids: assurance mode, the
# L5_cross_maximal cross-language model, randomised draw indices over the pilot
# posterior. prior_regime is "primary", which is the key ARC still runs; off ARC
# the same regime is called "primary".
#
# Seed base 1300000 with spacing 10 (languages use +0/1/2). This grid first took
# 940000, which generate_databased_v2_decision_grid.R has held since it entered
# the repository on 2026-07-31 and which the registry assigns to it; 108 of the
# rows below collided with that grid's 940000-941079. The decision arm has all
# 1,080 of its cells computed and this one has none, so this is the grid that
# moved. Do not restate other grids' ranges here. The comment this replaces listed
# only the pooled grids, every range in it correct, and the collision came from
# what it left out. The registry in docs/design_power_analysis_pipeline.md is the
# authority.
#
# Usage: Rscript scripts/generate_pooled_v2_N80_precision_grid.R [REPS]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(readr) })

args   <- commandArgs(trailingOnly = TRUE)
REPS   <- if (length(args) >= 1) as.integer(args[[1]]) else 380L
OUT    <- "config/design_grid_pooled_v2_N80_precision.csv"
NDRAWS <- 8000L

stopifnot(is.finite(REPS), REPS > 0L)

# Fixed seed so the grid is reproducible. This is the generator's own RNG, used
# only to draw the pilot-posterior draw indices below, and it is not a cell seed.
# generate_pooled_model_comparison_grid.R uses the same value, and both generators
# draw three times from that one stream, so the two grids reuse the same posterior
# draws in shifted roles: the 120 English indices of the model-comparison grid are
# the first 120 English indices here, and its Turkish and Norwegian indices recur
# here as English indices 121-360. That correlates the two arms' simulated data,
# so treat any comparison between them as paired rather than independent. Cell
# seeds are a separate matter and are allocated from the registry.
set.seed(20260819)

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
    seed = as.integer(1300000L + (.rep - 1L) * 10L)
  ) |>
  dplyr::select(-.rep) |>
  dplyr::relocate(n_participants)

readr::write_csv(grid, OUT)
message(sprintf(
  "[pooled N80 precision grid] %d cells (N = 80 per language, assurance) -> %s",
  nrow(grid), OUT))
message(sprintf(
  "[pooled N80 precision grid] seeds %d to %d; expect roughly %.0f core-hours at the arm's median",
  1300000L, 1300000L + (REPS - 1L) * 10L, REPS * 65))
