#!/usr/bin/env Rscript
# scripts/generate_feasibility_grid.R
# ---------------------------------------------------------------------------
# Generate the cross-language FEASIBILITY / convergence-timing grid.
#
# This is NOT a power run. It sweeps the full cross-language model ladder
# (L0_cross .. L5_cross) across simulated language counts (3, 10, 20) at the
# full per-language sample size (N = 80), with agent gender entered as the full
# 3-way S_Type x Semantics x Gender FIXED-effect interaction (no gender random
# slopes), to map where convergence and/or the ARC walltime break BEFORE
# data collection. A few seeds per cell bound run-to-run timing variance;
# convergence is read per .rds (runtime_sec + diagnostics), not as a power
# estimate.
#
# Standalone by design: it does NOT touch the power-analysis grids or the main
# generate_design_grid.R. It writes ONE new file, config/design_grid_feasibility.csv.
# prior_regime is written with whichever regime key the copy of R/03_define_priors.R
# deployed on the cluster defines. That deployment lags the name used here, so the value
# below is deliberately the cluster's key rather than this repository's.
#
# Usage (from design_analysis/ root):
#   Rscript scripts/generate_feasibility_grid.R \
#     [--out config/design_grid_feasibility.csv] [--b_feas 3]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)
})

option_list <- list(
  optparse::make_option("--out",    default = "config/design_grid_feasibility.csv"),
  optparse::make_option("--b_feas", default = 3L, type = "integer",
    help = "Seeds (replicates) per feasibility cell [default 3]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

# --- Assumed data-generating effect sizes (match the power-analysis grids) ---
BSEM    <- 0.8     # semantics (affectedness) main effect
BACT    <- -0.5    # Active x Semantics interaction
BPSEUDO <- 0.2     # pseudo-passive interaction
BGENDER             <- 0.3     # agent-gender main (nuisance)
BGENDER_SEM_PASSIVE <- 0.15    # small Gender x Semantics, passive-only (3-way nuisance)

# --- Replication sampler (so timings transfer to the production replicates) ---
ITER   <- 3000L
WARMUP <- 1000L
CHAINS <- 4L

# --- Feasibility design ------------------------------------------------------
N_FEAS       <- 80L                  # full per-language sample size
N_LANG_FEAS  <- c(3L, 10L, 20L)      # simulated language counts to map
B_FEAS       <- as.integer(opt$b_feas)
CROSS_LEVELS <- c(                   # full cross ladder, cheap -> heavy
  "L0_cross_intercepts_only",
  "L1_cross_intercepts_only_ppt_verb",
  "L2_cross_stype_participant_verb",
  "L3_cross_no_participant_interaction",
  "L4_cross_uncorrelated",
  "L5_cross_maximal",
  "L6_cross_gender_maximal",          # + gender random main effect (intermediate)
  "L7_cross_truly_maximal"            # full 3-way random slopes (hypothesised to fail)
)

# expand_grid keeps order: model_level varies slowest (L0 block first), then
# n_languages (3, 10, 20). Cheap cells therefore come first in the array.
conditions <- tidyr::expand_grid(
  model_level = CROSS_LEVELS,
  n_languages = N_LANG_FEAS
) |>
  dplyr::mutate(
    language                = "AllLanguages",
    prior_regime            = "primary",   # cluster regime key; see header
    threshold_mode          = "broad",
    n_participants          = N_FEAS,
    n_verbs                 = 20L,
    n_items_per_cell        = 1L,
    beta_semantics          = BSEM,
    beta_active_interaction = BACT,
    beta_pseudo_interaction = BPSEUDO,
    has_pseudo_passive      = TRUE,         # explicit: AllLanguages is not a config key
    iter                    = ITER,
    warmup                  = WARMUP,
    chains                  = CHAINS,
    gender_spec             = "interaction", # full 3-way S_Type*Semantics*Gender, FIXED-only
    include_gender          = TRUE,          # derived; lets a stale runner still get gender (main)
    beta_gender             = BGENDER,
    beta_gender_sem_passive = BGENDER_SEM_PASSIVE,
    .cond                   = dplyr::row_number()
  )

# Expand each condition into B_FEAS seeded replicates (deterministic seeds,
# base 500000 to avoid collision with the power grids' 1e5/2e5/3e5/4e5/9e5 bases).
grid <- tidyr::expand_grid(conditions, .rep = seq_len(B_FEAS)) |>
  dplyr::mutate(seed = as.integer(500000L + (.cond - 1L) * B_FEAS + (.rep - 1L))) |>
  dplyr::select(-.cond, -.rep)

COL_ORDER <- c(
  "language", "model_level", "prior_regime", "threshold_mode",
  "n_participants", "n_verbs", "n_items_per_cell",
  "beta_semantics", "beta_active_interaction", "beta_pseudo_interaction",
  "has_pseudo_passive", "iter", "warmup", "chains", "seed",
  "gender_spec", "include_gender", "beta_gender", "beta_gender_sem_passive", "n_languages"
)
grid <- dplyr::select(grid, dplyr::all_of(COL_ORDER))

readr::write_csv(grid, opt$out)
message(sprintf(
  "[feasibility grid] %d cells (%d cross levels x %d language counts x %d seeds) -> %s",
  nrow(grid), length(CROSS_LEVELS), length(N_LANG_FEAS), B_FEAS, opt$out))
message("[feasibility grid] cheap-first order: ", paste(CROSS_LEVELS, collapse = " -> "))
