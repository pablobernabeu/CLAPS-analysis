#!/usr/bin/env Rscript
# scripts/aggregate_databased_v2.R
# ---------------------------------------------------------------------------
# Aggregate the pilot-grounded (v2) design analysis into tidy summary CSVs.
# The report (reports/preliminary_sample_size_analysis.qmd) reads these files,
# so no simulation result needs to be hard-coded in the report source.
#
# Outputs (to --outdir):
#   joint_power_databased_v2.csv   One row per language x mode x N: replicate
#     count, marginal detection rates for H1a and H1b, and the joint rate
#     (both focal predictions clear BF >= 10 in the same simulated study).
#   pilot_params_ceilings.csv      One row per language: fitted pilot
#     quantities (verb count, affectedness spread, by-verb SDs, focal
#     posterior means and per-SD values, the slope's posterior SD, the joint
#     sign probability) and the analytic detection ceilings at BF >= 10, 6
#     and 3, with the verb counts needed for 80% and 90% joint ceilings.
#
# The ceiling formula: a one-sided Savage-Dickey criterion BF >= k under a
# symmetric zero-centred prior is met exactly when the posterior sign
# probability reaches k/(k+1), i.e. when the estimate sits qnorm(k/(k+1))
# standard errors from zero. With unlimited participants the standard error
# of a verb-level effect falls to tau/sqrt(Sxx) (by-verb SD over the
# affectedness spread; Clark 1973, doi:10.1016/S0022-5371(73)80014-3;
# Westfall, Kenny & Judd 2014, doi:10.1037/xge0000014), so the ceiling is
# the detection rate at that floor, averaged over the pilot posterior draws.
#
# Run on the cluster, next to the per-cell outputs and pilot DGP files:
#   Rscript scripts/aggregate_databased_v2.R \
#     --cells  <data>/outputs/design_databased_v2 \
#     --dgpdir <data>/outputs/pilot_models \
#     --outdir <data>/outputs/design_summary_databased_v2
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr) })

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--cells",  default = "outputs/design_databased_v2"),
  optparse::make_option("--dgpdir", default = "outputs/pilot_models"),
  optparse::make_option("--outdir", default = "outputs/design_summary_databased_v2")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Joint power from the per-cell results -------------------------------
files <- list.files(opt$cells, pattern = "\\.rds$", full.names = TRUE)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  d <- tryCatch(readRDS(files[[i]]), error = function(e) NULL)
  if (is.null(d)) next
  s <- d$summary; bf <- d$bf_results
  if (is.null(s) || is.null(bf) || is.null(s$status) || s$status != "success") next
  if (!("BF_10" %in% names(bf))) next
  h1a <- suppressWarnings(as.numeric(bf$BF_10[bf$hypothesis == "H1a_semantics_positive"]))
  h1b <- suppressWarnings(as.numeric(bf$BF_10[bf$hypothesis == "H1b_active_interaction_negative"]))
  if (length(h1a) != 1 || length(h1b) != 1 || is.na(h1a) || is.na(h1b)) next
  rows[[i]] <- data.frame(language = s$language, mode = s$mode,
                          n_participants = as.integer(s$n_participants),
                          h1a = h1a >= 10, h1b = h1b >= 10)
}
cells <- dplyr::bind_rows(rows)
joint <- cells |>
  dplyr::group_by(language, mode, n_participants) |>
  dplyr::summarise(reps = dplyr::n(), p_h1a = mean(h1a), p_h1b = mean(h1b),
                   p_joint = mean(h1a & h1b), .groups = "drop") |>
  dplyr::arrange(mode, language, n_participants)
readr::write_csv(joint, file.path(opt$outdir, "joint_power_databased_v2.csv"))
message(sprintf("[aggregate v2] %d usable cells -> joint_power_databased_v2.csv (%d rows)",
                nrow(cells), nrow(joint)))

# ---- 2. Pilot parameters and analytic ceilings ------------------------------
ceiling_at <- function(ba, bb, s_a, s_b, bf) {
  z <- qnorm(bf / (bf + 1))
  pa <- pnorm(ba / s_a - z)
  pb <- pnorm(-bb / s_b - z)
  c(h1a = mean(pa), h1b = mean(pb), joint = mean(pa * pb))
}
verbs_for <- function(ba, bb, s_a, s_b, n_verbs, target, bf = 10, cap = 1200L) {
  z <- qnorm(bf / (bf + 1))
  for (V in seq(n_verbs, cap, by = 2L)) {
    sc <- sqrt(V / n_verbs)
    j <- mean(pnorm(ba / (s_a / sc) - z) * pnorm(-bb / (s_b / sc) - z))
    if (j >= target) return(as.integer(V))
  }
  NA_integer_
}

LANGS <- c("English", "Turkish", "Norwegian")
params <- dplyr::bind_rows(lapply(LANGS, function(L) {
  d <- readRDS(file.path(opt$dgpdir, paste0("pilot_dgp_v2_pilot_", L, ".rds")))
  vs  <- d$verb_affectedness
  sxx <- sum((vs - mean(vs))^2)
  tau_int <- sqrt(d$Sigma_verb["Intercept", "Intercept"])
  tau_act <- sqrt(d$Sigma_verb["S_TypeActive", "S_TypeActive"])
  s_a <- tau_int / sqrt(sxx)
  s_b <- tau_act / sqrt(sxx)
  ba <- d$fixef_draws[, "Semantics_scaled"]
  bb <- d$fixef_draws[, "S_TypeActive:Semantics_scaled"]
  sd_x <- sd(vs)
  c10 <- ceiling_at(ba, bb, s_a, s_b, 10)
  c6  <- ceiling_at(ba, bb, s_a, s_b, 6)
  c3  <- ceiling_at(ba, bb, s_a, s_b, 3)
  data.frame(
    language = L, n_verbs = length(vs), sxx = sxx, affectedness_sd = sd_x,
    ndraws = nrow(d$fixef_draws),
    tau_verb_intercept = tau_int, tau_verb_active = tau_act,
    se_floor_h1a = s_a, se_floor_h1b = s_b,
    slope_mean = mean(ba), interaction_mean = mean(bb),
    slope_per_sd = mean(ba) * sd_x, interaction_per_sd = mean(bb) * sd_x,
    slope_post_sd = sd(ba),
    sign_prob_joint = mean(ba > 0 & bb < 0),
    ceiling_h1a_bf10 = c10[["h1a"]], ceiling_h1b_bf10 = c10[["h1b"]],
    ceiling_joint_bf10 = c10[["joint"]],
    ceiling_joint_bf6 = c6[["joint"]], ceiling_joint_bf3 = c3[["joint"]],
    verbs_for_80 = verbs_for(ba, bb, s_a, s_b, length(vs), 0.80),
    verbs_for_90 = verbs_for(ba, bb, s_a, s_b, length(vs), 0.90)
  )
}))
readr::write_csv(params, file.path(opt$outdir, "pilot_params_ceilings.csv"))
message("[aggregate v2] pilot_params_ceilings.csv written")
