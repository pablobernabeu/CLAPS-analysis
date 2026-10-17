#!/usr/bin/env Rscript
# stage2/04_meta_analysis.R
# ---------------------------------------------------------------------------
# Stage 2, step 4 of 6: the confirmatory tests of H1 and H2.
#
# Collects the per-variety estimates from step 3, averages each language's varieties
# into one estimate (the mean of their posterior means, with the mean of their
# posterior SDs as its standard error), and runs the registered Bayesian
# random-effects meta-analysis with bayesmeta: H1 over every language, H2 over the
# languages whose data include a pseudo-passive. Each hypothesis is decided by its
# directional Bayes factor: 10 or more is strong evidence for it, 1/10 or less strong
# evidence that the interaction runs the other way, anything in between inconclusive.
#
# Writes:
#   results/language_estimates.csv   one row per language and hypothesis
#   results/meta_analysis.csv        one row per hypothesis: the test and its result
#   results/left_out.csv             varieties with no converged model
#
# Usage (from design_analysis/):
#   Rscript stage2/04_meta_analysis.R [--dry-run]
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
lang <- read_language_table(cfg$languages_table)

files <- list.files(file.path(P$results, "estimates"), pattern = "[.]csv$", full.names = TRUE)
if (!length(files)) stop("[stage2:04] no estimates in ", file.path(P$results, "estimates"))
est <- dplyr::bind_rows(lapply(files, readr::read_csv, show_col_types = FALSE))

left_out <- dplyr::distinct(dplyr::filter(est, status != "converged"), Variety, status, levels_tried)
readr::write_csv(left_out, file.path(P$results, "left_out.csv"))
est <- dplyr::filter(est, status == "converged")

by_lang <- average_varieties(est, lang)
readr::write_csv(by_lang, file.path(P$results, "language_estimates.csv"))

tests <- dplyr::bind_rows(lapply(c("H1", "H2"), function(h) {
  d <- dplyr::filter(by_lang, hypothesis == h)
  if (nrow(d) < 2) {
    return(tibble::tibble(hypothesis = h, k = nrow(d), result = "not testable (fewer than 2 languages)"))
  }
  r <- run_meta_test(d$y, d$se, prior_mu_sd = cfg$meta_prior_sd[[h]],
                     tau_scale = cfg$inference$tau_prior_halfnormal_scale)
  dplyr::mutate(r, hypothesis = h, prior_mu_sd = cfg$meta_prior_sd[[h]],
                result = interpret_bf(bf, cfg$inference$bf_threshold), .before = 1)
}))
n_lang <- dplyr::n_distinct(by_lang$language)
tests$languages_total <- n_lang
tests$minimum_reached <- n_lang >= 20
readr::write_csv(tests, file.path(P$results, "meta_analysis.csv"))

for (i in seq_len(nrow(tests))) {
  t <- tests[i, ]
  if (is.null(t$bf) || is.na(t$bf)) {
    message(sprintf("[stage2:04] %s: %s", t$hypothesis, t$result))
  } else {
    message(sprintf("[stage2:04] %s: %d languages, mean %.3f per SD [%.3f, %.3f], tau %.3f, BF %.3g: %s",
                    t$hypothesis, t$k, t$mu_mean, t$mu_lo, t$mu_hi, t$tau_median, t$bf, t$result))
  }
}
if (nrow(left_out)) message("[stage2:04] left out, no converged model: ", paste(left_out$Variety, collapse = ", "))
if (n_lang < 20) message("[stage2:04] ", n_lang, " languages, fewer than the registered minimum of 20")
