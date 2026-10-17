#!/usr/bin/env Rscript
# stage2/02_affectedness.R
# ---------------------------------------------------------------------------
# Stage 2, step 2 of 6: the affectedness predictor and the outcome-neutral check.
#
# For each variety, takes the mean affectedness rating of each verb over the
# analysed raters and both AGENT versions, checks that the verbs' means have an SD of
# at least the registered minimum (the outcome-neutral check, which reads only the
# ratings), and centres the means and divides them by two of their SDs (Gelman,
# 2008). A variety that fails the check is left out of the confirmatory analyses and
# reported.
#
# Writes:
#   processed/verb_affectedness.csv   per variety and verb: mean rating, scaled score
#   results/outcome_neutral.csv       per variety: number of verbs, SD, passes
#
# Usage (from design_analysis/):
#   Rscript stage2/02_affectedness.R [--dry-run]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("R/03_define_priors.R")
source("R/14_stage2.R")

args <- commandArgs(trailingOnly = TRUE)
cfg <- stage2_settings(stage2_config(), dry_run = "--dry-run" %in% args)
P <- cfg$paths

rat <- readr::read_csv(file.path(P$processed, "ratings_analysed.csv"), show_col_types = FALSE)
means <- verb_affectedness(rat)
check <- outcome_neutral_check(means, cfg$outcome_neutral$min_affectedness_sd)
scaled <- scale_verb_affectedness(means)

readr::write_csv(scaled, file.path(P$processed, "verb_affectedness.csv"))
readr::write_csv(check, file.path(P$results, "outcome_neutral.csv"))

message(sprintf("[stage2:02] outcome-neutral check passed by %d of %d varieties",
                sum(check$passes), nrow(check)))
if (any(!check$passes)) {
  message("[stage2:02] left out of the confirmatory analyses: ",
          paste(sprintf("%s (SD %.2f)", check$Variety[!check$passes],
                        check$sd_verb_means[!check$passes]), collapse = ", "))
}
