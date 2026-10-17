#!/usr/bin/env Rscript
# stage2/06_robustness_families.R
# ---------------------------------------------------------------------------
# Stage 2, step 6 of 6: the second planned robustness check, languages nested within
# their Glottolog families.
#
# Repeats the meta-analysis of step 4 on the same language estimates, with a family
# level added: the same normal prior on the mean, the same half-normal(0.5) prior on
# the between-language SD and a half-normal(0.5) prior on the between-family SD,
# fitted in brms with the registered sampler and convergence criteria. Reported beside
# the confirmatory tests; it decides neither hypothesis.
#
# Writes:
#   results/robustness_families.csv
#
# Usage (from design_analysis/):
#   Rscript stage2/06_robustness_families.R [--dry-run]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source("R/03_define_priors.R")
source("R/07_extract_diagnostics.R")
source("R/14_stage2.R")

args <- commandArgs(trailingOnly = TRUE)
cfg <- stage2_settings(stage2_config(), dry_run = "--dry-run" %in% args)
P <- cfg$paths

by_lang <- readr::read_csv(file.path(P$results, "language_estimates.csv"), show_col_types = FALSE)
out <- dplyr::bind_rows(lapply(c("H1", "H2"), function(h) {
  d <- dplyr::filter(by_lang, hypothesis == h)
  if (nrow(d) < 2) return(tibble::tibble(hypothesis = h, k = nrow(d), status = "not testable"))
  r <- fit_family_meta(d, prior_mu_sd = cfg$meta_prior_sd[[h]],
                       tau_scale = cfg$inference$tau_prior_halfnormal_scale,
                       sampler = cfg$sampler, seed = cfg$seed_base + 1000L)
  sds <- brms::VarCorr(r$fit)
  tibble::tibble(hypothesis = h, k = nrow(d), families = dplyr::n_distinct(d$family),
                 status = r$status, mu_mean = r$mu_mean, mu_sd = r$mu_sd,
                 sd_family = sds$family$sd[1, "Estimate"],
                 sd_language = sds$`family:language`$sd[1, "Estimate"],
                 p_negative = r$p_negative, bf = r$bf, bf_bound = r$bound,
                 result = interpret_bf(r$bf, cfg$inference$bf_threshold))
}))
readr::write_csv(out, file.path(P$results, "robustness_families.csv"))
print(out)
