#!/usr/bin/env Rscript
# scripts/generate_databased_grid.R
# Grid for the data-grounded (parametric-bootstrap) BAYESIAN power analysis.
# Each cell draws data from the pilot-calibrated DGP (R/10) at a given N and a
# focal-effect multiplier, refits the maximal model, and computes Bayes factors.
#   effect_mult = 1.00 : the language's own pilot effect (primary).
#   effect_mult = 0.60 : a conservative / safeguard effect (winner's-curse guard).
# The optimistic (cross-language literature pool) case is the existing parametric
# run, so it is not repeated here.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(optparse) })

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--out", default = "config/design_grid_databased.csv"),
  optparse::make_option("--b",   default = 24L, type = "integer", help = "replicates per cell")
)))
B <- as.integer(opt$b)

langs   <- tibble::tibble(language = c("English", "Turkish", "Norwegian"),
                          has_pseudo_passive = c(TRUE, TRUE, FALSE))
N_SWEEP <- c(30L, 50L, 70L, 100L, 130L)
MULT    <- c(0.6, 1.0)

cond <- tidyr::crossing(langs, n_participants = N_SWEEP, effect_mult = MULT) |>
  dplyr::mutate(model_level = "L5_correlated_maximal", prior_regime = "primary",
                threshold_mode = "broad", iter = 3000L, warmup = 1000L, chains = 4L,
                .cond = dplyr::row_number())

grid <- tidyr::crossing(cond, .rep = seq_len(B)) |>
  dplyr::mutate(seed = as.integer(800000L + (.cond - 1L) * B + (.rep - 1L))) |>   # base 8e5: collision-free
  dplyr::select(-.cond, -.rep)

COL <- c("language", "has_pseudo_passive", "n_participants", "effect_mult",
         "model_level", "prior_regime", "threshold_mode", "iter", "warmup", "chains", "seed")
grid <- dplyr::select(grid, dplyr::all_of(COL))
readr::write_csv(grid, opt$out)
message(sprintf("[databased grid] %d cells (%d langs x %d N x %d mult x %d reps) -> %s",
                nrow(grid), nrow(langs), length(N_SWEEP), length(MULT), B, opt$out))
