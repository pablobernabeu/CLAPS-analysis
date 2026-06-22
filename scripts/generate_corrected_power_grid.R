#!/usr/bin/env Rscript
# scripts/generate_corrected_power_grid.R
# ---------------------------------------------------------------------------
# Corrected SINGLE-LANGUAGE power BFDA grid.
#
# Two corrections relative to the original power grid:
#   (1) the per-verb affectedness fix in R/06_simulate_design.R (one fixed
#       Semantics value per Verb_ID, drawn before the participant/S_Type loops),
#       which re-specifies affectedness as the verb-level predictor it really is;
#   (2) n_verbs is SWEPT as a design dimension alongside N, so the recommendation
#       is an N x n_verbs surface rather than a single N at n_verbs = 20.
#
# Single-language L5 only (the per-language recommendation, and the model that
# correctly omits a by-verb Semantics slope, so it is unaffected by the
# cross-language identifiability issue the fix exposes). prior_regime = "primary"
# to match the priors deployed on ARC.
#
# Usage (from design_analysis/ root):
#   Rscript scripts/generate_corrected_power_grid.R [--out config/design_grid_corrected.csv] [--b 50]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(optparse) })

option_list <- list(
  optparse::make_option("--out", default = "config/design_grid_corrected.csv"),
  optparse::make_option("--b",   default = 50L, type = "integer", help = "replicates per cell")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
B <- as.integer(opt$b)

# Assumed data-generating effect sizes (match the original power grid).
BSEM <- 0.8; BACT <- -0.5; BPSEUDO <- 0.2
ITER <- 3000L; WARMUP <- 1000L; CHAINS <- 4L

langs <- tibble::tibble(
  language           = c("English", "Turkish", "Norwegian"),
  has_pseudo_passive = c(TRUE, TRUE, FALSE),
  beta_pseudo        = c(BPSEUDO, BPSEUDO, 0)      # Norwegian has no pseudo-passive
)
N_SWEEP     <- c(50L, 60L, 70L, 80L, 90L, 100L, 120L)   # extended up: the fix may raise required N
NVERB_SWEEP <- c(12L, 20L, 40L, 72L)                    # 20 = original; 72 ~ the pilot's verb count

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
  dplyr::mutate(seed = as.integer(600000L + (.cond - 1L) * B + (.rep - 1L))) |>   # base 6e5: collision-free
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
message(sprintf("[corrected grid] %d cells (%d langs x %d N x %d n_verbs x %d reps) -> %s",
                nrow(grid), nrow(langs), length(N_SWEEP), length(NVERB_SWEEP), B, opt$out))
