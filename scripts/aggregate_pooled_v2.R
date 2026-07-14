#!/usr/bin/env Rscript
# scripts/aggregate_pooled_v2.R
# ---------------------------------------------------------------------------
# Aggregate the pooled cross-language design analysis (design_pooled_v2) into
# a tidy summary CSV the report can read, mirroring aggregate_databased_v2.R.
# One row per per-language sample size, with detection rates for the two focal
# predictions and their joint at Bayes-factor thresholds 10, 6 and 3. The
# language label on every cell is "AllLanguages"; n_participants is per
# language (total sample is three times that).
#
#   Rscript scripts/aggregate_pooled_v2.R \
#     --cells  <data>/outputs/design_pooled_v2 \
#     --outdir <data>/outputs/design_summary_databased_v2
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr) })

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--cells",  default = "outputs/design_pooled_v2"),
  optparse::make_option("--outdir", default = "outputs/design_summary_databased_v2")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(opt$cells, pattern = "\\.rds$", full.names = TRUE)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- tryCatch(readRDS(files[[i]]), error = function(e) NULL)
  if (is.null(d)) next
  s <- d$summary; bf <- d$bf_results
  if (is.null(s) || is.null(bf) || is.null(s$status) || s$status != "success") next
  if (!("BF_10" %in% names(bf))) next
  g <- function(h) suppressWarnings(as.numeric(bf$BF_10[bf$hypothesis == h]))
  h1a <- g("H1a_semantics_positive"); h1b <- g("H1b_active_interaction_negative")
  if (length(h1a) != 1 || length(h1b) != 1 || is.na(h1a) || is.na(h1b)) next
  rows[[i]] <- data.frame(n_participants = as.integer(s$n_participants),
                          bf_h1a = h1a, bf_h1b = h1b)
}
cells <- dplyr::bind_rows(rows)
if (nrow(cells) == 0) { message("[aggregate pooled] no usable cells yet"); quit(save = "no") }

pooled <- cells |>
  dplyr::group_by(n_participants) |>
  dplyr::summarise(
    reps           = dplyr::n(),
    p_h1a          = mean(bf_h1a >= 10), p_h1b = mean(bf_h1b >= 10),
    p_joint        = mean(bf_h1a >= 10 & bf_h1b >= 10),
    p_joint_bf6    = mean(bf_h1a >= 6  & bf_h1b >= 6),
    p_joint_bf3    = mean(bf_h1a >= 3  & bf_h1b >= 3),
    median_bf_h1a  = median(bf_h1a),  median_bf_h1b = median(bf_h1b),
    .groups = "drop") |>
  dplyr::arrange(n_participants)
readr::write_csv(pooled, file.path(opt$outdir, "pooled_power_databased_v2.csv"))
message(sprintf("[aggregate pooled] %d usable cells -> pooled_power_databased_v2.csv (%d N rows)",
                nrow(cells), nrow(pooled)))
