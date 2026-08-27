#!/usr/bin/env Rscript
# scripts/aggregate_pooled.R
# ---------------------------------------------------------------------------
# Aggregate the pooled cross-language design analysis (design_pooled_v2) into
# a tidy summary CSV the report can read, mirroring aggregate_pilot.R.
# One row per per-language sample size, with detection rates for the two focal
# predictions and their joint at Bayes-factor thresholds 10, 6 and 3, the
# interaction alone at 19, and the rates over the converged subset. The
# language label on every cell is "AllLanguages"; n_participants is per
# language (total sample is three times that).
#
#   Rscript scripts/aggregate_pooled.R \
#     --cells  <data>/outputs/design_pooled_v2 \
#     --outdir <data>/outputs/design_summary_pilot
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr) })

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--cells",  default = "outputs/design_pooled_v2"),
  optparse::make_option("--outdir", default = "outputs/design_summary_pilot")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Unusable cells are counted by reason rather than dropped silently; see the
# equivalent block in aggregate_pilot.R for why (a failed fit still exits
# 0 at the scheduler level, so this tally is the first place a cluster-wide
# fault becomes visible). Pooled cells are the expensive ones, at roughly one to
# four days each, so losing them unnoticed is correspondingly costlier.
skipped <- c(unreadable = 0L, no_summary = 0L, fit_error = 0L,
             no_bf_column = 0L, bf_not_extractable = 0L)
bump <- function(reason) skipped[[reason]] <<- skipped[[reason]] + 1L

files <- list.files(opt$cells, pattern = "\\.rds$", full.names = TRUE)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- tryCatch(readRDS(files[[i]]), error = function(e) NULL)
  if (is.null(d)) { bump("unreadable"); next }
  s <- d$summary; bf <- d$bf_results
  if (is.null(s) || is.null(bf) || is.null(s$status)) { bump("no_summary"); next }
  if (s$status != "success") { bump("fit_error"); next }
  if (!("BF_10" %in% names(bf))) { bump("no_bf_column"); next }
  g <- function(h) suppressWarnings(as.numeric(bf$BF_10[bf$hypothesis == h]))
  h1a <- g("H1a_semantics_positive"); h1b <- g("H1b_active_interaction_negative")
  if (length(h1a) != 1 || length(h1b) != 1 || !is.finite(h1a) || !is.finite(h1b)) {
    bump("bf_not_extractable"); next
  }
  dg <- d$diagnostics
  ok <- if (!is.null(dg) && !is.null(dg$convergence_ok)) isTRUE(dg$convergence_ok[[1]]) else NA
  rows[[i]] <- data.frame(n_participants = as.integer(s$n_participants),
                          bf_h1a = h1a, bf_h1b = h1b, converged = ok)
}
cells <- dplyr::bind_rows(rows)

if (sum(skipped) > 0) {
  message(sprintf("[aggregate pooled] skipped %d of %d files (%s)", sum(skipped), length(files),
                  paste(sprintf("%s=%d", names(skipped)[skipped > 0], skipped[skipped > 0]),
                        collapse = ", ")))
}
if (nrow(cells) == 0) { message("[aggregate pooled] no usable cells yet"); quit(save = "no") }

pooled <- cells |>
  dplyr::group_by(n_participants) |>
  dplyr::summarise(
    reps           = dplyr::n(),
    p_h1a          = mean(bf_h1a >= 10), p_h1b = mean(bf_h1b >= 10),
    p_joint        = mean(bf_h1a >= 10 & bf_h1b >= 10),
    # The MARGINAL rates at the relaxed thresholds, added 2026-08-12. Only the
    # joint rate was carried at 6 and 3, which left the single most important
    # quantity for the collaborators uncomputable from this file: the chance of
    # detecting the interaction ALONE, without also requiring the passive-only
    # affectedness slope. Powering for H1b by itself is the design question now
    # on the table, and the stored per-cell Bayes factors already answer it, so
    # this is a re-scoring rather than a re-run.
    p_h1a_bf6      = mean(bf_h1a >= 6),
    p_h1b_bf6      = mean(bf_h1b >= 6),
    p_joint_bf6    = mean(bf_h1a >= 6  & bf_h1b >= 6),
    p_h1a_bf3      = mean(bf_h1a >= 3),
    p_h1b_bf3      = mean(bf_h1b >= 3),
    p_joint_bf3    = mean(bf_h1a >= 3  & bf_h1b >= 3),
    median_bf_h1a  = median(bf_h1a),  median_bf_h1b = median(bf_h1b),
    # Reported, not applied: see aggregate_pilot.R. Filtering the pooled
    # cells on convergence moves the headline rate by only a few points, so the
    # low joint power at BF >= 10 is a property of the design rather than an
    # artefact of imperfect sampling.
    n_converged       = sum(converged %in% TRUE),
    p_joint_converged = if (any(converged %in% TRUE)) {
      mean(bf_h1a[converged %in% TRUE] >= 10 & bf_h1b[converged %in% TRUE] >= 10)
    } else NA_real_,
    # The interaction alone over the converged subset, and at a Bayes factor of
    # 19 (posterior sign probability 0.95). Added 2026-08-27; see aggregate_pilot.R.
    p_h1b_converged   = if (any(converged %in% TRUE)) {
      mean(bf_h1b[converged %in% TRUE] >= 10)
    } else NA_real_,
    p_h1b_bf19        = mean(bf_h1b >= 19),
    # Monte-Carlo standard error of the headline joint rate, so that readers of
    # the CSV are not left to infer the precision from the replicate count.
    mcse_p_joint      = sqrt(p_joint * (1 - p_joint) / dplyr::n()),
    # The same for the marginal H1b rate. The pooled cells carry far fewer
    # replicates than the per-language ones (tens rather than hundreds), so any
    # pooled rate quoted without this is quoted more precisely than it is known.
    mcse_p_h1b        = sqrt(p_h1b * (1 - p_h1b) / dplyr::n()),
    .groups = "drop") |>
  dplyr::arrange(n_participants)
readr::write_csv(pooled, file.path(opt$outdir, "pooled_power_pilot.csv"))
message(sprintf("[aggregate pooled] %d usable cells -> pooled_power_pilot.csv (%d N rows)",
                nrow(cells), nrow(pooled)))
