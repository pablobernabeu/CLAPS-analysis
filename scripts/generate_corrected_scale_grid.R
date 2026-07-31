#!/usr/bin/env Rscript
# scripts/generate_corrected_scale_grid.R
# ---------------------------------------------------------------------------
# SUPERSEDED ARM — RETAINED FOR THE RECORD, NOT PART OF THE REPORTED ANALYSIS.
#
# One of the "abandoned literature-anchored refinement arms" of commit 3c5a872
# (2026-07-24), which records that the literature cross-check "stands on its
# completed English and Turkish cells" rather than on the refinements that were
# once promised. The report was reworded at that point to stop promising results
# from this arm.
#
# Status on the cluster. The held array job 12419137 was CANCELLED on 2026-07-30,
# having sat JobHeldUser since 20 June. It had produced NOTHING: no
# outputs/design_corrected_scale directory was ever created. Nothing reads that
# directory in any case — no aggregator takes it as input, and the report cites
# it nowhere. Releasing it would have committed thousands of CPU-hours of shared
# credit to a 2,100-cell grid whose first output would also have been its last
# use.
#
# The prior-scale concern this arm was built to quantify is real and is not
# dismissed by the cancellation: it is documented in
# R/06_simulate_design_gelman.R, which explains the sqrt(3) overstatement and
# implements the Gelman-scaled analysis path. What was abandoned is the plan to
# quantify it with a dedicated 2,100-cell sweep, not the finding itself.
#
# The generator is kept in case the arm is revived. Note its seed base moved from
# 7e5 to 1.0e6 on 2026-07-30; since nothing was ever computed, that re-basing
# orphaned nothing.
#
# ---------------------------------------------------------------------------
# Prior-SCALE-corrected single-language power BFDA grid.
#
# Same design-generating assumptions as scripts/generate_corrected_power_grid.R
# (per-verb affectedness fix; assumed effects beta_semantics = 0.8, active
# interaction = -0.5, pseudo interaction = 0.2). The ONLY difference at run time
# is that the cells are executed with scripts/04_design_analysis_cell_gelman.R,
# which sources R/06_simulate_design_gelman.R so the ANALYSED focal predictor is
# Gelman-scaled (SD 0.5, matching R/02_preprocess_factors.R and the real
# confirmatory analysis) rather than the raw U(-0.5,0.5) draw (SD 0.289). This
# corrects the prior-scale mismatch surfaced by the 2026-06-20 audit, under which
# the base corrected grid overstates the focal-slope Bayes factor by ~sqrt(3).
#
# Focused confirmation: verb counts {40, 72} only (the design-relevant region the
# report uses; 12/20 are known to be insufficient and are dropped). Distinct seed
# base (7e5) so it never collides with the base corrected grid (6e5). Writes to
# config/design_grid_corrected_scale.csv; cells write to outputs/design_corrected_scale.
#
# Usage (from design_analysis/ root):
#   Rscript scripts/generate_corrected_scale_grid.R [--out config/design_grid_corrected_scale.csv] [--b 50]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(optparse) })

option_list <- list(
  optparse::make_option("--out", default = "config/design_grid_corrected_scale.csv"),
  optparse::make_option("--b",   default = 50L, type = "integer", help = "replicates per cell")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
B <- as.integer(opt$b)

# Assumed data-generating effect sizes (identical to the base corrected grid).
BSEM <- 0.8; BACT <- -0.5; BPSEUDO <- 0.2
ITER <- 3000L; WARMUP <- 1000L; CHAINS <- 4L

langs <- tibble::tibble(
  language           = c("English", "Turkish", "Norwegian"),
  has_pseudo_passive = c(TRUE, TRUE, FALSE),
  beta_pseudo        = c(BPSEUDO, BPSEUDO, 0)      # Norwegian has no pseudo-passive
)
N_SWEEP     <- c(50L, 60L, 70L, 80L, 90L, 100L, 120L)
NVERB_SWEEP <- c(40L, 72L)                          # design-relevant region only

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
  # Seed base 1.0e6, exclusive to this grid: 1000000-1002099 (42 conditions x B = 50).
  # Moved off 7e5 on 2026-07-30, where it had overlapped generate_floor50_power_grid.R
  # and generate_safeguard_grid.R. This grid was chosen to move because it has never
  # been run (no outputs/design_corrected_scale directory exists on ARC), so re-basing
  # orphaned nothing. See the seed registry in
  # docs/design_power_analysis_pipeline.md.
  dplyr::mutate(seed = as.integer(1000000L + (.cond - 1L) * B + (.rep - 1L))) |>
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
message(sprintf("[corrected-scale grid] %d cells (%d langs x %d N x %d n_verbs x %d reps) -> %s",
                nrow(grid), nrow(langs), length(N_SWEEP), length(NVERB_SWEEP), B, opt$out))
