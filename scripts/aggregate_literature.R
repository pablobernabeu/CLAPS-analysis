#!/usr/bin/env Rscript
# scripts/aggregate_literature.R
# ---------------------------------------------------------------------------
# Aggregate the literature-anchored design analysis (design_corrected) into the
# summary CSV the report reads, mirroring aggregate_pilot.R and
# aggregate_pooled.R. One row per language x participant count x verb count x
# hypothesis.
#
# This exists because the report's literature summary previously had no
# committed producer. It was written once by hand in June, while one language
# was still computing, and could not be regenerated when those cells landed, so
# the report went on reporting a two-language cross-check for six weeks. The
# repo's general aggregator, scripts/06_aggregate_design_results.R, writes a
# different and wider set of columns for the whole design_analysis pipeline; it
# omits the convergence-subset columns the report needs, which is why this
# narrower script exists alongside it rather than replacing it.
#
#   Rscript scripts/aggregate_literature.R \
#     --cells  <data>/outputs/design_corrected \
#     --outdir <data>/outputs/design_summary_literature
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr) })

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--cells",  default = "outputs/design_corrected"),
  optparse::make_option("--outdir", default = "outputs/design_summary_literature")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# Unusable cells are counted by reason rather than dropped in silence; see the
# equivalent block in aggregate_pilot.R for why. A failed fit still exits 0 at
# the scheduler level, so this tally is the first place a fault that destroyed
# part of a grid becomes visible.
skipped <- c(unreadable = 0L, no_summary = 0L, fit_error = 0L, no_bf_column = 0L)
bump <- function(reason) skipped[[reason]] <<- skipped[[reason]] + 1L

files <- list.files(opt$cells, pattern = "\\.rds$", full.names = TRUE)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- tryCatch(readRDS(files[[i]]), error = function(e) NULL)
  if (is.null(d)) { bump("unreadable"); next }
  s <- d$summary; bf <- d$bf_results; g <- d$diagnostics
  if (is.null(s) || is.null(bf) || is.null(s$status)) { bump("no_summary"); next }
  if (!identical(s$status, "success")) { bump("fit_error"); next }
  if (!("BF_10" %in% names(bf))) { bump("no_bf_column"); next }
  rows[[i]] <- data.frame(
    language       = s$language,
    n_participants = as.integer(s$n_participants),
    n_verbs        = as.integer(s$n_verbs),
    hypothesis     = bf$hypothesis,
    bf             = suppressWarnings(as.numeric(bf$BF_10)),
    converged      = if (!is.null(g$convergence_ok)) isTRUE(g$convergence_ok[[1]]) else NA
  )
}
cells <- dplyr::bind_rows(rows)

if (sum(skipped) > 0) {
  message(sprintf("[aggregate literature] skipped %d of %d files (%s)",
                  sum(skipped), length(files),
                  paste(sprintf("%s=%d", names(skipped)[skipped > 0], skipped[skipped > 0]),
                        collapse = ", ")))
}
if (nrow(cells) == 0) stop("[aggregate literature] no usable cells in ", opt$cells)

# p_bf_primary is the detection rate at the preregistered threshold of 10 and
# p_bf_secondary the rate at the moderate threshold of 3. The *_conv columns
# repeat the primary rate over the converged fits alone, so the report can show
# that the recommendation does not rest on marginal fits. Convergence is
# reported rather than applied, as in the other two aggregators.
exc <- cells |>
  dplyr::filter(is.finite(bf)) |>
  dplyr::group_by(language, n_participants, n_verbs, hypothesis) |>
  dplyr::summarise(
    reps              = dplyr::n(),
    n_conv            = sum(converged %in% TRUE),
    p_bf_primary      = mean(bf >= 10),
    p_bf_primary_conv = if (any(converged %in% TRUE)) {
      mean(bf[converged %in% TRUE] >= 10)
    } else NA_real_,
    p_bf_secondary    = mean(bf >= 3),
    p_convergence_ok  = mean(converged %in% TRUE),
    median_bf         = median(bf),
    .groups = "drop") |>
  dplyr::arrange(language, n_verbs, n_participants, hypothesis)

readr::write_csv(exc, file.path(opt$outdir, "bf_exceedance_literature.csv"))
message(sprintf("[aggregate literature] %d usable cells -> bf_exceedance_literature.csv (%d rows, %s)",
                nrow(cells), nrow(exc), paste(sort(unique(exc$language)), collapse = ", ")))

thin <- dplyr::filter(exc, reps < 10L)
if (nrow(thin) > 0) {
  message(sprintf("[aggregate literature] %d cell(s) below 10 replicates, treat as provisional",
                  nrow(thin)))
}
