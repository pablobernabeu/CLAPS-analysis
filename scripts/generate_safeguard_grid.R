#!/usr/bin/env Rscript
# scripts/generate_safeguard_grid.R
# ---------------------------------------------------------------------------
# SUPERSEDED ARM — RETAINED FOR THE RECORD, NOT PART OF THE REPORTED ANALYSIS.
#
# Both purposes described below were abandoned on 2026-07-24 (commit 3c5a872,
# "Report: de-stale the abandoned literature-anchored refinement arms"), which
# states that "the low-N (20-40 participant) arm ... and the 75/60%
# discounted-effect sensitivity arm were all promised as forthcoming, but their
# jobs have been held (admin/user) and idle for a month and are superseded: ...
# the safeguard is now built into the data-grounded engine, and the discounted
# scenario is realised by conditioning on the pilot effects directly."
#
# Those are precisely the two arms this file generates: discount = 1.00 at
# N = 20/30/40 is the low-N arm, and discount = 0.75/0.60 is the
# discounted-effect arm. What replaced them is the safeguard MODE of the
# pilot-grounded engine (R/10_simulate_from_pilot_v2.R, mode = "safeguard"),
# whose results reach the report through outputs/design_databased_v2 and
# outputs/design_summary_pilot/joint_power_pilot.csv. Do not confuse the two: the
# report's "safeguard" figures come from that mode, never from this grid.
#
# Status on the cluster. The held array job 12409496 (tasks 123-440) was
# CANCELLED on 2026-07-30, having sat JobHeldUser since 19 June. Its 122
# completed outputs remain in outputs/design_safeguard and are read by nothing:
# no aggregator takes that directory as input, and the report cites it nowhere.
# Finishing the array would have cost roughly 3,000 CPU-hours of shared credit
# for outputs no analysis consumes.
#
# The generator is kept because the arm may be wanted again for a revision or a
# reviewer request, and because deleting it would erase the record of what was
# tried. If it is ever revived, note that its seed base moved from 7e5 to 1.1e6
# on 2026-07-30, so a regenerated grid will not match the 122 outputs on disk.
#
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
  # Seed base 1.1e6, exclusive to this grid: 1100000-1100439 (11 conditions x B = 40).
  # Moved off 7e5 on 2026-07-30, where it had overlapped
  # generate_floor50_power_grid.R (700000-701049) and
  # generate_corrected_scale_grid.R. One of the three had to move regardless, and
  # this grid had the fewest computed cells to lose.
  #
  # CONSEQUENCE OF THE MOVE: the 122 cells already computed under
  # outputs/design_safeguard carry the old seeds (700280-700401) in their filenames,
  # so they no longer match what this grid asks for and the resume-by-existing-output
  # check will not find them. They are not deleted; regenerating this grid and
  # resubmitting recomputes them under the new seeds. See the seed registry in
  # docs/design_power_analysis_pipeline.md.
  dplyr::mutate(seed = as.integer(1100000L + (.cond - 1L) * B + (.rep - 1L))) |>
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
