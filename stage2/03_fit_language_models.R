#!/usr/bin/env Rscript
# stage2/03_fit_language_models.R
# ---------------------------------------------------------------------------
# Stage 2, step 3 of 6: the individual-language models.
#
# Fits each variety's cumulative-logit model by the registered ladder (R/09, starting
# at the maximal model of the Stage 1 report), with the primary priors, the
# registered sampler and the registered convergence criteria, and extracts the two
# interactions per SD of affectedness from the most complex level that converged. A
# variety whose models all fail the criteria is recorded as such and left out of the
# meta-analysis.
#
# One variety per call, so that the cluster can run them as an array
# (hpc/submit_stage2_language_fits.sh). Each call writes
#   fits/<variety>_primary_<level>.rds       the selected fit (by R/09)
#   fits/diagnostics/...                     diagnostics of every level tried
#   results/estimates/<variety>.csv          the extracted interactions, or a row
#                                            stating that no level converged
#
# Usage (from design_analysis/):
#   Rscript stage2/03_fit_language_models.R --variety "Spanish (Chilean)" [--dry-run]
#   Rscript stage2/03_fit_language_models.R --index 7 [--dry-run]   (row 7 of the analysed varieties)
#   Rscript stage2/03_fit_language_models.R --all [--dry-run]
# A level killed by the cluster's wall time is a permitted reason to descend
# (R/09_model_ladder.R, FALLBACK_TRIGGERS). The ladder is then resumed at the next
# level with --start-level <level>, and the timeout is recorded in the run log.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("R/02_preprocess_factors.R")
source("R/03_define_priors.R")
source("R/09_model_ladder.R")
source("R/14_stage2.R")

args <- commandArgs(trailingOnly = TRUE)
aval <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) NULL else args[i + 1]
}
cfg <- stage2_settings(stage2_config(), dry_run = "--dry-run" %in% args)
P <- cfg$paths
lang <- read_language_table(cfg$languages_table)

acc <- readr::read_csv(file.path(P$processed, "acceptability_analysed.csv"), show_col_types = FALSE)
scores <- readr::read_csv(file.path(P$processed, "verb_affectedness.csv"), show_col_types = FALSE)
check <- readr::read_csv(file.path(P$results, "outcome_neutral.csv"), show_col_types = FALSE)
analysed <- check$Variety[check$passes]
# Varieties in the order of the language table, so that --index and the seeds are stable.
analysed <- lang$variety[lang$variety %in% analysed]

todo <- if ("--all" %in% args) {
  analysed
} else if (!is.null(aval("--variety"))) {
  aval("--variety")
} else if (!is.null(aval("--index"))) {
  analysed[as.integer(aval("--index"))]
} else {
  stop("give --variety, --index or --all")
}
if (anyNA(todo) || !all(todo %in% analysed)) {
  stop("[stage2:03] not an analysed variety: ", paste(setdiff(todo, analysed), collapse = ", "))
}

label_of <- function(v) gsub("[^A-Za-z0-9]+", "_", v)
dir.create(file.path(P$results, "estimates"), recursive = TRUE, showWarnings = FALSE)
for (v in todo) {
  d <- build_variety_data(dplyr::filter(acc, Variety == v), dplyr::filter(scores, Variety == v))
  seed <- cfg$seed_base + match(v, lang$variety)
  message(sprintf("[stage2:03] %s: %d ratings, %d participants, %d verbs, pseudo-passive %s, seed %d",
                  v, nrow(d), nlevels(d$Participant), nlevels(d$Verb),
                  attr(d, "has_pseudo_passive"), seed))
  res <- fit_variety_model(d, cfg, out_dir = P$fits, label = label_of(v), seed = seed,
                           start_level = aval("--start-level"))
  out <- if (is.null(res$fit)) {
    tibble::tibble(Variety = v, hypothesis = NA_character_, model_level = NA_character_,
                   status = "no_converged_model")
  } else {
    dplyr::mutate(extract_interactions(res$fit), Variety = v, model_level = res$selected_model,
                  status = "converged", .before = 1)
  }
  out$has_pseudo_passive <- attr(d, "has_pseudo_passive")
  out$levels_tried <- paste(names(res$ladder_log), collapse = " > ")
  readr::write_csv(out, file.path(P$results, "estimates", paste0(label_of(v), ".csv")))
}
