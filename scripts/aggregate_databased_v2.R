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
# Every cell that cannot be used is counted against a reason rather than being
# dropped in silence. This matters operationally: a SLURM array task whose model
# fails still exits 0, because the failure is caught and recorded inside the
# result file, so a fault that destroys most of a grid is invisible in the
# scheduler's own accounting and shows up here only as a smaller cell count.
# The July 2026 node-image fault (Stan failing to compile on part of the
# cluster) destroyed 901 of 1,080 decision-arm cells and was diagnosed only by
# reading the stored status field by hand. The tally below surfaces exactly that
# class of failure on every run.
skipped <- c(unreadable = 0L, no_summary = 0L, fit_error = 0L,
             no_bf_column = 0L, bf_not_extractable = 0L)
bump <- function(reason) skipped[[reason]] <<- skipped[[reason]] + 1L

files <- list.files(opt$cells, pattern = "\\.rds$", full.names = TRUE)
rows <- vector("list", length(files))
for (i in seq_along(files)) {
  # A truncated or partially written .rds (job killed mid-save) fails to
  # deserialise; treat it as missing rather than aborting the whole aggregation.
  d <- tryCatch(readRDS(files[[i]]), error = function(e) NULL)
  if (is.null(d)) { bump("unreadable"); next }
  s <- d$summary; bf <- d$bf_results
  if (is.null(s) || is.null(bf) || is.null(s$status)) { bump("no_summary"); next }
  # status is written by the cell runner: "success" only if the fit completed
  # and the Bayes factors were computed. Anything else carries an error message.
  if (s$status != "success") { bump("fit_error"); next }
  if (!("BF_10" %in% names(bf))) { bump("no_bf_column"); next }
  h1a <- suppressWarnings(as.numeric(bf$BF_10[bf$hypothesis == "H1a_semantics_positive"]))
  h1b <- suppressWarnings(as.numeric(bf$BF_10[bf$hypothesis == "H1b_active_interaction_negative"]))
  # Guard against a hypothesis label appearing zero or twice, and against a
  # non-finite Bayes factor (a Savage-Dickey ratio can overflow when the
  # posterior density at the null is numerically zero).
  if (length(h1a) != 1 || length(h1b) != 1 || !is.finite(h1a) || !is.finite(h1b)) {
    bump("bf_not_extractable"); next
  }
  # Convergence is recorded per cell but deliberately does NOT filter the main
  # power columns, because excluding non-converged fits would silently condition
  # the estimate on the fit having gone well, which is not how the real analysis
  # would behave. It is carried through as a separate column so the sensitivity
  # can be reported without re-running anything (see the *_converged columns).
  g <- d$diagnostics
  ok <- if (!is.null(g) && !is.null(g$convergence_ok)) isTRUE(g$convergence_ok[[1]]) else NA
  rows[[i]] <- data.frame(language = s$language, mode = s$mode,
                          n_participants = as.integer(s$n_participants),
                          bf_h1a = h1a, bf_h1b = h1b, converged = ok)
}
cells <- dplyr::bind_rows(rows)

if (sum(skipped) > 0) {
  message(sprintf("[aggregate v2] skipped %d of %d files (%s)", sum(skipped), length(files),
                  paste(sprintf("%s=%d", names(skipped)[skipped > 0], skipped[skipped > 0]),
                        collapse = ", ")))
}
if (nrow(cells) == 0) stop("[aggregate v2] no usable cells in ", opt$cells)
# The continuous Bayes factors are stored per cell, so detection rates at any
# threshold are a re-scoring, not a re-run. BF >= 10 keeps the original
# column names (the report reads them); 6 and 3 quantify the criterion
# relaxations under discussion with the collaborators.
joint <- cells |>
  dplyr::group_by(language, mode, n_participants) |>
  dplyr::summarise(
    reps         = dplyr::n(),
    p_h1a        = mean(bf_h1a >= 10),
    p_h1b        = mean(bf_h1b >= 10),
    p_joint      = mean(bf_h1a >= 10 & bf_h1b >= 10),
    p_h1a_bf6    = mean(bf_h1a >= 6),
    p_h1b_bf6    = mean(bf_h1b >= 6),
    p_joint_bf6  = mean(bf_h1a >= 6 & bf_h1b >= 6),
    p_h1a_bf3    = mean(bf_h1a >= 3),
    p_h1b_bf3    = mean(bf_h1b >= 3),
    p_joint_bf3  = mean(bf_h1a >= 3 & bf_h1b >= 3),
    # Convergence sensitivity, reported alongside rather than applied. n_converged
    # counts cells meeting all of R-hat < 1.01, bulk and tail ESS >= 400, zero
    # divergences and zero max-treedepth saturations (Vehtari et al. 2021,
    # doi:10.1214/20-BA1221); p_joint_converged recomputes the headline rate on
    # that subset, and is NA when no cell in the group qualifies.
    n_converged  = sum(converged %in% TRUE),
    p_joint_converged = if (any(converged %in% TRUE)) {
      mean(bf_h1a[converged %in% TRUE] >= 10 & bf_h1b[converged %in% TRUE] >= 10)
    } else NA_real_,
    .groups = "drop") |>
  dplyr::arrange(mode, language, n_participants)
readr::write_csv(joint, file.path(opt$outdir, "joint_power_databased_v2.csv"))
message(sprintf("[aggregate v2] %d usable cells -> joint_power_databased_v2.csv (%d rows)",
                nrow(cells), nrow(joint)))

# A cell holding only a handful of replicates yields a power estimate that can
# read as 0% or 100% on chance alone; the report suppresses these via its own
# minimum-replicate guard, but the CSV is also read directly, so name them here.
thin <- dplyr::filter(joint, reps < 10L)
if (nrow(thin) > 0) {
  message(sprintf("[aggregate v2] %d cell(s) below 10 replicates, treat as provisional: %s",
                  nrow(thin), paste(sprintf("%s/%s N=%d (%d)", thin$language, thin$mode,
                                            thin$n_participants, thin$reps), collapse = "; ")))
}

# ---- 2. Pilot parameters and analytic ceilings ------------------------------
# Detection ceiling at threshold `bf`: the probability of clearing it with an
# unlimited number of participants, i.e. once the standard error has fallen to
# the verb-level floor. Averaging over the pilot posterior draws (`ba`, `bb`)
# makes this an assurance rather than a point calculation, so it carries the
# uncertainty in the pilot effects rather than conditioning on a point estimate
# (Albers & Lakens 2018, doi:10.1016/j.jesp.2017.09.004).
#
# The joint term multiplies the two probabilities within each draw. That is
# exact in the draw-to-draw correlation between the effects, which is what the
# assurance framing needs, but it treats the two test statistics as independent
# given the draw. They in fact share the by-verb random effects, so the joint
# ceiling is a mild over-estimate; the simulations, which have no such
# assumption, come in slightly below it, as expected.
ceiling_at <- function(ba, bb, s_a, s_b, bf) {
  z <- qnorm(bf / (bf + 1))          # BF >= k  <=>  estimate >= qnorm(k/(k+1)) SEs from 0
  pa <- pnorm(ba / s_a - z)          # H1a is directional-positive
  pb <- pnorm(-bb / s_b - z)         # H1b is directional-negative, hence the sign flip
  c(h1a = mean(pa), h1b = mean(pb), joint = mean(pa * pb))
}

# Smallest verb count reaching `target` joint ceiling. The standard error of a
# verb-level effect scales as 1/sqrt(V), so adding verbs shrinks the floor while
# adding participants cannot. Steps in twos because the design pairs verbs
# across agent gender, and returns NA when `cap` is reached, which is a real
# answer rather than a failure: for Turkish the requirement runs into the
# hundreds of verbs and the remedy has to be pooling instead.
verbs_for <- function(ba, bb, s_a, s_b, n_verbs, target, bf = 10, cap = 1200L) {
  z <- qnorm(bf / (bf + 1))
  for (V in seq(n_verbs, cap, by = 2L)) {
    sc <- sqrt(V / n_verbs)          # SE shrinks by sqrt(V / V_current)
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
