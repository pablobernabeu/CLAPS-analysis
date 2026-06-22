#!/usr/bin/env Rscript
# scripts/generate_safeguard_grid.R
# ---------------------------------------------------------------------------
# Effect-size sensitivity ("safeguard") arm + low-N localisation, at the real
# 72-verb design (English, L5 maximal correlated model, primary priors).
#
# Two purposes:
#   (1) discount = 1.00 at N below 50 localises the minimum adequate N at the
#       assumed effect (the main grid only swept N >= 50, where power is already
#       saturated at 72 verbs);
#   (2) discount = 0.75 and 0.60 quantify the N required under deliberately
#       conservative ("safeguard") effect sizes, addressing the winner's-curse /
#       Type-M qualification: published anchors may overstate the true effect.
#
# The discount multiplies the FOCAL data-generating effects only; priors are
# unchanged (they encode pre-data belief, not the assumed truth). Seeds use base
# 700000 so they never collide with the main corrected run (base 600000).
#
# Usage (from design_analysis/ root):
#   Rscript scripts/generate_safeguard_grid.R [--out config/design_grid_safeguard.csv] [--b 40]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(optparse) })

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--out", default = "config/design_grid_safeguard.csv"),
  optparse::make_option("--b",   default = 40L, type = "integer", help = "replicates per cell")
)))
B <- as.integer(opt$b)

# Assumed (full-strength) data-generating effect sizes, matched to the main grid.
BSEM <- 0.8; BACT <- -0.5; BPSEUDO <- 0.2
ITER <- 3000L; WARMUP <- 1000L; CHAINS <- 4L

# (discount, N) conditions — all English, 72 verbs.
conds <- dplyr::bind_rows(
  tidyr::crossing(discount = 1.00, n_participants = c(20L, 30L, 40L)),         # localise min N
  tidyr::crossing(discount = 0.75, n_participants = c(40L, 60L, 80L, 100L)),   # safeguard 75%
  tidyr::crossing(discount = 0.60, n_participants = c(60L, 90L, 120L, 150L))   # safeguard 60%
) |>
  dplyr::mutate(
    language                = "English",
    has_pseudo_passive      = TRUE,
    n_verbs                 = 72L,
    model_level             = "L5_correlated_maximal",
    prior_regime            = "primary",
    threshold_mode          = "broad",
    n_items_per_cell        = 1L,
    beta_semantics          = BSEM    * discount,   # discounted focal effects
    beta_active_interaction = BACT    * discount,
    beta_pseudo_interaction = BPSEUDO * discount,
    iter = ITER, warmup = WARMUP, chains = CHAINS,
    gender_spec             = "none",
    include_gender          = FALSE,
    beta_gender             = 0.3,
    beta_gender_sem_passive = 0.15
  ) |>
  dplyr::mutate(.cond = dplyr::row_number())

grid <- tidyr::crossing(conds, .rep = seq_len(B)) |>
  dplyr::mutate(seed = as.integer(700000L + (.cond - 1L) * B + (.rep - 1L))) |>   # base 7e5: collision-free
  dplyr::select(-.cond, -.rep, -discount)

COL_ORDER <- c(
  "language", "model_level", "prior_regime", "threshold_mode",
  "n_participants", "n_verbs", "n_items_per_cell",
  "beta_semantics", "beta_active_interaction", "beta_pseudo_interaction",
  "has_pseudo_passive", "iter", "warmup", "chains", "seed",
  "gender_spec", "include_gender", "beta_gender", "beta_gender_sem_passive"
)
grid <- dplyr::select(grid, dplyr::all_of(COL_ORDER))
readr::write_csv(grid, opt$out)
message(sprintf("[safeguard grid] %d cells (%d conditions x %d reps) -> %s",
                nrow(grid), nrow(conds), B, opt$out))
