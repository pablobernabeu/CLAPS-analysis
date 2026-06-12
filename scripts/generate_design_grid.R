#!/usr/bin/env Rscript
# scripts/generate_design_grid.R
# ---------------------------------------------------------------------------
# Generate the replicated design grids for the Bayes-factor *power* analysis.
#
# The earlier feasibility study used one seed per design point (n_sims = 1),
# which yields a 0/1 BF-exceedance indicator rather than a power estimate.
# This generator expands every design CONDITION into `n_simulations_per_cell`
# independently seeded replicates — one SLURM array task per replicate — so the
# exceedance proportion across replicates is a genuine Monte-Carlo power estimate.
#
# Output is split into separate grids per "study" so that (a) each stays within
# the cluster MaxArraySize (5000) and (b) each can be submitted with tailored
# resources:
#   design_grid_single.csv  single-language power curves + prior sensitivity
#   design_grid_gender.csv  gender-covariate variation power curves
#   design_grid_cross.csv   cross-language (pooled) reduced power analysis
#
# Sampler: replicates use the lighter `replication_sampler` (4 chains x 3000
# iter, 1000 warmup, adapt_delta 0.99) — adequate for a single BF per fit and
# convergence-checked per replicate. The heavy 16x5000 sampler in `model:` is
# reserved for the one-off maximal-model convergence demonstration.
#
# Usage:
#   Rscript scripts/generate_design_grid.R \
#     [--config config/analysis_config.yaml] [--out_dir config] \
#     [--b 200] [--b_cross 50]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(yaml)
  library(optparse)
})

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0) return(b)
  if (length(a) == 1 && is.na(a)) return(b)
  a
}

option_list <- list(
  optparse::make_option("--config",  default = "config/analysis_config.yaml"),
  optparse::make_option("--out_dir", default = "config"),
  optparse::make_option("--b",       default = NA_integer_, type = "integer",
    help = "Replicates per single-language/gender design point [default: config]"),
  optparse::make_option("--b_cross", default = NA_integer_, type = "integer",
    help = "Replicates per cross-language design point [default: config]")
)
opt <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

cfg <- yaml::read_yaml(opt$config)

B       <- opt$b       %||% cfg$design_analysis$n_simulations_per_cell       %||% 200L
B_cross <- opt$b_cross %||% cfg$design_analysis$n_simulations_per_cell_cross %||% 50L
B       <- as.integer(B)
B_cross <- as.integer(B_cross)

rs     <- cfg$replication_sampler %||% list()
ITER   <- as.integer(rs$iter   %||% 3000L)
WARMUP <- as.integer(rs$warmup %||% 1000L)
CHAINS <- as.integer(rs$chains %||% 4L)

# --- Assumed true effect sizes (the data-generating values the power analysis
# conditions on). These match the pilot/literature anchors used throughout. ---
BSEM    <- 0.8    # semantics (affectedness) main effect (H1a)
BACT    <- -0.5   # Active x Semantics interaction (H1b)
BPSEUDO <- 0.2    # pseudo-passive interaction (languages with a pseudo-passive)
BGENDER <- 0.3    # referent-gender covariate (gender variation)

# Languages and their pseudo-passive availability.
langs <- tibble::tibble(
  language           = c("English", "Turkish", "Norwegian"),
  has_pseudo_passive = c(TRUE, TRUE, FALSE),
  beta_pseudo        = c(BPSEUDO, BPSEUDO, 0)   # Norwegian has no pseudo-passive level
)

N_SWEEP <- c(30L, 40L, 50L, 60L)   # participant counts for the power curves

COL_ORDER <- c(
  "language", "model_level", "prior_regime", "threshold_mode",
  "n_participants", "n_verbs", "n_items_per_cell",
  "beta_semantics", "beta_active_interaction", "beta_pseudo_interaction",
  "has_pseudo_passive", "iter", "warmup", "chains", "seed",
  "include_gender", "beta_gender"
)

# Expand a table of design CONDITIONS into B seeded replicates. Seeds are unique
# and deterministic within a grid: seed_base + (condition - 1) * B + (rep - 1).
expand_reps <- function(conditions, B, seed_base) {
  conditions |>
    dplyr::mutate(.cond = dplyr::row_number()) |>
    tidyr::crossing(.rep = seq_len(B)) |>
    dplyr::mutate(seed = as.integer(seed_base + (.cond - 1L) * B + (.rep - 1L))) |>
    dplyr::select(-.cond, -.rep)
}

mk_conditions <- function(df) {
  df |>
    dplyr::mutate(
      n_verbs                 = 20L,
      n_items_per_cell        = 1L,
      beta_semantics          = BSEM,
      beta_active_interaction = BACT,
      beta_pseudo_interaction = beta_pseudo,
      iter = ITER, warmup = WARMUP, chains = CHAINS
    ) |>
    dplyr::select(-beta_pseudo)
}

# === Study A: single-language baseline power curves + prior sensitivity ======
power_curves <- langs |>
  tidyr::crossing(n_participants = N_SWEEP) |>
  dplyr::mutate(model_level = "L5_correlated_maximal", prior_regime = "primary",
                threshold_mode = "broad", include_gender = FALSE, beta_gender = 0) |>
  mk_conditions()

prior_sens <- langs |>
  tidyr::crossing(prior_regime = c("weak", "literature_centred", "heavy_tailed")) |>
  dplyr::mutate(model_level = "L5_correlated_maximal", threshold_mode = "broad",
                n_participants = 50L, include_gender = FALSE, beta_gender = 0) |>
  mk_conditions()

grid_single <- dplyr::bind_rows(
  expand_reps(power_curves, B, seed_base = 100000L),
  expand_reps(prior_sens,  B, seed_base = 200000L)
)

# === Study B: gender-variation power curves ==================================
gender_curves <- langs |>
  tidyr::crossing(n_participants = N_SWEEP) |>
  dplyr::mutate(model_level = "L5_correlated_maximal", prior_regime = "primary",
                threshold_mode = "broad", include_gender = TRUE, beta_gender = BGENDER) |>
  mk_conditions()

grid_gender <- expand_reps(gender_curves, B, seed_base = 300000L)

# === Study C: cross-language (pooled) reduced power analysis =================
# Uses the faster, better-converging L4 cross-uncorrelated model. The L5
# cross-maximal convergence/feasibility is documented by the earlier heavy run.
cross_curves <- tibble::tibble(
  language = "AllLanguages", model_level = "L4_cross_uncorrelated",
  prior_regime = "primary", threshold_mode = "broad",
  n_participants = N_SWEEP, has_pseudo_passive = TRUE, beta_pseudo = BPSEUDO,
  include_gender = FALSE, beta_gender = 0
) |>
  mk_conditions()

grid_cross <- expand_reps(cross_curves, B_cross, seed_base = 900000L)

# === Study D: extended single-language N-sweep (70-100) =====================
# The 30-60 sweep showed the H1b Active x affectedness interaction does not reach
# 80% power by N=60 in any language; this extension locates the required sample
# size per language. Baseline (non-gender), primary prior, maximal model.
N_EXTEND <- c(70L, 80L, 90L, 100L)
power_curves_ext <- langs |>
  tidyr::crossing(n_participants = N_EXTEND) |>
  dplyr::mutate(model_level = "L5_correlated_maximal", prior_regime = "primary",
                threshold_mode = "broad", include_gender = FALSE, beta_gender = 0) |>
  mk_conditions()
grid_extend <- expand_reps(power_curves_ext, B, seed_base = 400000L)

# --- Write ------------------------------------------------------------------
write_grid <- function(g, name) {
  g <- dplyr::select(g, dplyr::all_of(COL_ORDER))
  path <- file.path(opt$out_dir, name)
  readr::write_csv(g, path)
  message(sprintf("[grid] %-26s %5d rows (%d conditions x B)", name, nrow(g),
                  nrow(dplyr::distinct(dplyr::select(g, -seed)))))
  if (nrow(g) > 5000) warning("[grid] ", name, " exceeds MaxArraySize 5000!")
  invisible(path)
}

message("[grid] B = ", B, " (single/gender), B_cross = ", B_cross,
        " | sampler: ", CHAINS, " chains x ", ITER, " iter (", WARMUP, " warmup)")
write_grid(grid_single, "design_grid_single.csv")
write_grid(grid_gender, "design_grid_gender.csv")
write_grid(grid_cross,  "design_grid_cross.csv")
write_grid(grid_extend, "design_grid_extend.csv")
message("[grid] Total replicate fits: ",
        nrow(grid_single) + nrow(grid_gender) + nrow(grid_cross) + nrow(grid_extend))
