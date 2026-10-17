#!/usr/bin/env Rscript
# stage2/05_robustness_single_model.R
# ---------------------------------------------------------------------------
# Stage 2, step 5 of 6: the first planned robustness check, a single model across
# languages.
#
# Stacks the analysed varieties, with participants and verbs nested within variety
# and Language as the language, and fits the single cumulative-logit model by the
# registered ladder (R/14_stage2.R, build_single_model_ladder: from the maximal model
# of the Stage 1 report down to random intercepts), with the primary priors, the
# registered sampler and convergence criteria. The two interactions are reported with
# their directional Bayes factors from the most complex level that converged. The
# check is reported beside the confirmatory tests and decides neither hypothesis.
#
# Writes:
#   fits/single_model_primary_<level>.rds and fits/diagnostics/...
#   results/robustness_single_model.csv
#
# Usage (from design_analysis/):
#   Rscript stage2/05_robustness_single_model.R [--dry-run] [--start-level R2_no_participant_interaction]
#   (on the cluster: hpc/submit_stage2_single_model.sh)
# --start-level resumes the ladder after a level was killed by the wall time, a
# permitted reason to descend (R/09_model_ladder.R, FALLBACK_TRIGGERS).
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
start_level <- if ("--start-level" %in% args) args[match("--start-level", args) + 1] else NULL
cfg <- stage2_settings(stage2_config(), dry_run = "--dry-run" %in% args)
P <- cfg$paths
lang <- read_language_table(cfg$languages_table)

acc <- readr::read_csv(file.path(P$processed, "acceptability_analysed.csv"), show_col_types = FALSE)
scores <- readr::read_csv(file.path(P$processed, "verb_affectedness.csv"), show_col_types = FALSE)
check <- readr::read_csv(file.path(P$results, "outcome_neutral.csv"), show_col_types = FALSE)
analysed <- lang$variety[lang$variety %in% check$Variety[check$passes]]

pooled <- build_pooled_data(lapply(analysed, function(v) {
  build_variety_data(dplyr::filter(acc, Variety == v), dplyr::filter(scores, Variety == v))
}), lang)
has_pp <- "Pseudo_Passive" %in% levels(pooled$S_Type)
message(sprintf("[stage2:05] %d ratings, %d languages, %d varieties, %d participants, %d verbs",
                nrow(pooled), nlevels(pooled$Language), length(analysed),
                nlevels(pooled$Participant), nlevels(pooled$Verb)))

s <- cfg$sampler
res <- run_model_ladder(
  df            = pooled,
  prior_obj     = build_brms_prior("primary", "broad", has_pseudo_passive = has_pp),
  sampling_args = list(iter = s$iter, warmup = s$warmup, chains = s$chains, cores = s$cores,
                       seed = cfg$seed_base),
  control_args  = list(adapt_delta = s$adapt_delta, max_treedepth = s$max_treedepth),
  has_pseudo_passive = has_pp,
  start_level   = start_level,
  out_dir       = P$fits,
  language      = "single_model",
  prior_label   = "primary",
  ladder        = build_single_model_ladder()
)

out <- if (is.null(res$fit)) {
  tibble::tibble(hypothesis = c("H1", "H2"), model_level = NA_character_, status = "no_converged_model")
} else {
  terms <- STAGE2_TERMS[STAGE2_TERMS %in% brms::variables(res$fit)]
  dplyr::bind_rows(lapply(names(terms), function(h) {
    x <- as.data.frame(res$fit, variable = terms[[h]])[[1]]
    b <- bf_from_draws(x)
    tibble::tibble(hypothesis = h, term = terms[[h]], model_level = res$selected_model,
                   status = "converged", per_sd_mean = mean(x) / 2, per_sd_sd = stats::sd(x) / 2,
                   per_sd_lo = stats::quantile(x, 0.025) / 2, per_sd_hi = stats::quantile(x, 0.975) / 2,
                   p_negative = b$p_negative, bf = b$bf, bf_bound = b$bound,
                   result = interpret_bf(b$bf, cfg$inference$bf_threshold))
  }))
}
out$levels_tried <- paste(names(res$ladder_log), collapse = " > ")
readr::write_csv(out, file.path(P$results, "robustness_single_model.csv"))
print(out)
