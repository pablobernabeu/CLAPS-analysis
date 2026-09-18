#!/usr/bin/env Rscript
# scripts/gender_versions_power.R
# ---------------------------------------------------------------------------
# What the second agent-gender version of every item is worth in power.
#
# THE QUESTION
#   The PI asked on 16 September 2026 whether the reversed male/female agent
#   version of each item ("The woman pushed the man" against "The man pushed the
#   woman", in every sentence type) could be dropped, halving a session from 40 to
#   20 minutes. Every pilot participant rated both versions, so the question is
#   what the second rating of each participant-by-verb-by-sentence-type cell buys
#   in precision on the verb-level focal terms, and what that precision is worth.
#
# WHAT THIS SCRIPT DOES, AND WHAT IT LEAVES TO OTHERS
#   It does not estimate the precision ratio. It takes any
#
#       r = SE(two versions) / SE(one version)
#
#   and turns it into percentage points of power, so that whichever estimate of r
#   the pilot supports can be read straight off the table.
#
#   The baseline is already the one-version design. simulate_from_pilot_v2() in
#   R/10_simulate_from_pilot_v2.R builds one row per participant, verb and
#   sentence type, so the decision-arm figures sent to the PI (English 0.90 for
#   H1b at N = 80) describe the 20-minute session. Dropping the second version
#   lowers none of them; keeping it could only raise them, by what this computes.
#
# THE REDUCTION
#   A directional Savage-Dickey Bayes factor under a symmetric prior clears B
#   exactly when the posterior sign probability reaches B/(B+1), which under a
#   normal posterior is |estimate| / SE >= z* = qnorm(B/(B+1)), 1.335 at B = 10
#   (R/05_hypothesis_tests.R). Assurance power at an effective SE s is then
#
#       H1b (predicted negative):  mean over pilot draws of Phi(-theta_d / s - z*)
#       H1a (predicted positive):  mean over pilot draws of Phi( theta_d / s - z*)
#
#   with theta_d the pilot posterior draws that the decision arm itself drew its
#   generating values from. One number, s, is solved per language, hypothesis and
#   threshold so that this reproduces the brms decision arm at N = 80. s is an
#   EFFECTIVE standard error: it absorbs what the normal approximation leaves out
#   (prior shrinkage, a posterior SD that differs from the sampling SD, fits that
#   went badly), which is why it is calibrated and not computed. The same
#   reduction gives the ceilings in scripts/aggregate_pilot.R and carries the
#   K-language engine in R/13_simulate_klanguage_pooled.R.
#
#   Scaling s by r assumes the second version changes only the spread of the
#   estimate. That is right for a precision question: the verb-level estimand and
#   the verb affectedness predictor are the same in both versions (the
#   gender-version affectedness varies within verb by about 0.6% of its spread
#   between verbs).
#
# CHECKS THE SCRIPT RUNS ON ITSELF
#   1. s is also solved at N = 70 and N = 100 and set against the verb-information
#      floor tau_verb / sqrt(Sxx), below which no participant count can take it.
#   2. A model s(N)^2 = a + b/N is tried three ways: freely from the three
#      calibrations (with a parametric bootstrap of their Monte-Carlo error), with
#      a pinned at the floor and b taken from the N = 80 cell, and with a pinned
#      and b fitted to all three cells by binomial maximum likelihood. Only an
#      identified model converts a gain in precision into extra participants.
#   3. An independent route to the same model from the DGP variance components,
#      using the Fisher information of the cumulative-logit likelihood, which also
#      gives the value r would take if the two versions were independent replicates.
#
# OUTPUT
#   outputs/design_summary_pilot/gender_versions_power.csv, one long table whose
#   `block` column separates:
#     calibration    s at N = 70, 80, 100 per language, hypothesis and threshold
#     identification the a + b/N model, free and floor-anchored
#     decomposition  the DGP-component route and the independent-replicate r
#     grid           power at r = 0.70 ... 1.00 and the gain over one version
#     pooled         K-language pooled power under the same r, and the analytic
#                    pooled SE ratio
#
# USAGE (from design_analysis/)
#   "C:/Program Files/R/R-4.6.1/bin/Rscript.exe" scripts/gender_versions_power.R
#   Options: --boot=1000 --pooled_reps=4000
#   About 22 minutes on one core at the defaults, 15 of them in the pooled
#   section; --pooled_reps=1000 brings the total to about 10.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(stats)
  library(dplyr)
  library(MASS)
})

# The other scripts assume design_analysis/ as the working directory; accept the
# repository root too, since that is where an interactive session usually starts.
if (!dir.exists("R") && dir.exists("design_analysis/R")) setwd("design_analysis")
stopifnot(file.exists("R/13_simulate_klanguage_pooled.R"))

args <- commandArgs(trailingOnly = TRUE)
aval <- function(f, d) {
  h <- grep(paste0("^", f, "="), args, value = TRUE)
  if (!length(h)) d else sub(paste0("^", f, "="), "", h[1])
}
N_BOOT      <- as.integer(aval("--boot", 1000L))
POOLED_REPS <- as.integer(aval("--pooled_reps", 4000L))

OUT     <- "outputs/design_summary_pilot"
DGP_DIR <- "outputs/pilot_models"
LANGS   <- c("English", "Turkish", "Norwegian")
BFS     <- c(10, 6, 3)
N_REF   <- 80L
N_ALL   <- c(70L, 80L, 100L)
R_GRID  <- round(seq(0.70, 1.00, by = 0.01), 2)
R_KEY   <- c(0.85, 0.90, 0.95, 0.98)
TERM    <- c(h1a = "Semantics_scaled", h1b = "S_TypeActive:Semantics_scaled")
# Sign of the predicted effect. H1b is the primary per-language test since the
# PI's decision of 25 August; H1a is carried for the joint figure only.
DIR     <- c(h1a = 1, h1b = -1)
# Minutes per participant with and without the second version, from the PI's email.
MIN_TWO <- 40
MIN_ONE <- 20

zstar <- function(B) qnorm(B / (B + 1))
rule  <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

# ---------------------------------------------------------------------------
# 1. Inputs: the decision arm, the pilot draws, the floors
# ---------------------------------------------------------------------------
arm <- read.csv(file.path(OUT, "joint_power_pilot.csv"))
arm <- arm[arm$mode == "assurance" & arm$n_participants %in% N_ALL, ]
PP  <- read.csv(file.path(OUT, "pilot_params_ceilings.csv"))

arm_col <- function(h, B) paste0("p_", h, if (B == 10) "" else paste0("_bf", B))
arm_get <- function(L, N, col) {
  row <- arm[arm$language == L & arm$n_participants == N, ]
  if (nrow(row) != 1L) stop("decision arm: expected one row for ", L, " N = ", N)
  c(p = row[[col]], reps = row$reps)
}

DGP <- setNames(lapply(LANGS, function(L)
  readRDS(file.path(DGP_DIR, paste0("pilot_dgp_v2_pilot_", L, ".rds")))), LANGS)
THETA <- lapply(DGP, function(d) {
  fd <- as.matrix(d$fixef_draws)
  list(h1a = as.numeric(fd[, TERM[["h1a"]]]), h1b = as.numeric(fd[, TERM[["h1b"]]]))
})
floor_of <- function(L, h) PP[PP$language == L, paste0("se_floor_", h)]

rule("INPUTS")
print(arm[, c("language", "n_participants", "reps", "p_h1a", "p_h1b", "p_joint", "mcse_p_h1b")],
      row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# 2. The normal approximation and its inversion
# ---------------------------------------------------------------------------
power_at <- function(s, theta, dir, z) {
  vapply(s, function(si) mean(pnorm(dir * theta / si - z)), numeric(1))
}

# As s -> 0 power tends to the posterior probability of the predicted sign, and
# as s -> Inf to Phi(-z*) = 1/(B+1). A target outside that range has no s, and
# returning the boundary (0 or Inf) is the honest answer: the approximation
# cannot reproduce that figure at any precision.
calibrate_s <- function(p, theta, dir, z) {
  hi <- mean(dir * theta > 0)
  lo <- pnorm(-z)
  if (p >= hi) return(0)
  if (p <= lo) return(Inf)
  exp(uniroot(function(ls) power_at(exp(ls), theta, dir, z) - p,
              c(log(1e-4), log(1e3)), tol = 1e-10)$root)
}

# Power is not guaranteed monotone in s: a wrong-sign draw gains power as s grows.
# The inversion assumes monotonicity, so it is checked rather than assumed.
for (L in LANGS) for (h in names(TERM)) for (B in BFS) {
  sg <- exp(seq(log(0.01), log(3), length.out = 400))
  pw <- power_at(sg, THETA[[L]][[h]], DIR[[h]], zstar(B))
  if (any(diff(pw) > 1e-9)) warning("power not monotone in s: ", L, " ", h, " BF", B)
}

# A fast inverse for the bootstrap, from a fine grid on log s. Exact inversion by
# uniroot at every bootstrap draw would cost minutes for no gain in accuracy.
make_inverse <- function(theta, dir, z) {
  sg <- exp(seq(log(1e-3), log(20), length.out = 1500))
  pw <- power_at(sg, theta, dir, z)
  hi <- mean(dir * theta > 0); lo <- pnorm(-z)
  function(p) {
    out <- suppressWarnings(approx(rev(pw), rev(log(sg)), xout = p, ties = "ordered")$y)
    out <- exp(out)
    out[p >= max(pw)] <- ifelse(p[p >= max(pw)] >= hi, 0, min(sg))
    out[p <= min(pw)] <- ifelse(p[p <= min(pw)] <= lo, Inf, max(sg))
    out
  }
}

# ---------------------------------------------------------------------------
# 3. Calibration at N = 70, 80, 100
# ---------------------------------------------------------------------------
rule("CALIBRATED EFFECTIVE SE, s, BY LANGUAGE, HYPOTHESIS, THRESHOLD AND N")

calib <- do.call(rbind, lapply(LANGS, function(L) do.call(rbind, lapply(names(TERM), function(h)
  do.call(rbind, lapply(BFS, function(B) do.call(rbind, lapply(N_ALL, function(N) {
    a  <- arm_get(L, N, arm_col(h, B))
    k  <- round(a[["p"]] * a[["reps"]])
    ci <- binom.test(k, a[["reps"]])$conf.int
    th <- THETA[[L]][[h]]; z <- zstar(B)
    fl <- floor_of(L, h)
    s  <- calibrate_s(a[["p"]], th, DIR[[h]], z)
    data.frame(block = "calibration", language = L, hypothesis = h, bf = B,
               n_participants = N, reps = a[["reps"]], p_arm = a[["p"]],
               mcse_arm = sqrt(a[["p"]] * (1 - a[["p"]]) / a[["reps"]]),
               p_arm_lo = ci[1], p_arm_hi = ci[2],
               s = s,
               # The CI on s is the image of the Clopper-Pearson interval on p, so it
               # is exact to the approximation and carries only Monte-Carlo error.
               s_lo = calibrate_s(ci[2], th, DIR[[h]], z),
               s_hi = calibrate_s(ci[1], th, DIR[[h]], z),
               se_floor = fl, s_over_floor = s / fl,
               ceiling_at_floor = power_at(fl, th, DIR[[h]], z),
               sign_prob = mean(DIR[[h]] * th > 0))
  }))))))))

# A calibration is flagged, not silently used, when its figure sits at or above
# what the verb count allows. Above the ceiling means the simulations beat the
# DGP's own infinite-participant limit, which Monte-Carlo error can do near it.
calib$flag <- with(calib, ifelse(
  s == 0 | is.infinite(s), "outside the approximation's range",
  ifelse(s < se_floor, "below verb floor: arm exceeds the infinite-N ceiling",
  ifelse(s_lo <= se_floor, "near ceiling: MC interval reaches the verb floor", "ok"))))
print(calib[calib$bf == 10, c("language", "hypothesis", "n_participants", "reps", "p_arm",
                              "mcse_arm", "s", "s_lo", "s_hi", "se_floor", "s_over_floor",
                              "ceiling_at_floor", "flag")],
      row.names = FALSE, digits = 3)
cat("\nBF 6 and BF 3 at N = 80 (under the approximation s should not depend on the threshold):\n")
print(calib[calib$n_participants == N_REF, c("language", "hypothesis", "bf", "p_arm", "s", "s_lo",
                                              "s_hi", "flag")], row.names = FALSE, digits = 3)

s_ref <- function(L, h, B) calib$s[calib$language == L & calib$hypothesis == h &
                                     calib$bf == B & calib$n_participants == N_REF]

# ---------------------------------------------------------------------------
# 4. Identification of s(N)^2 = a + b/N
# ---------------------------------------------------------------------------
rule("IS s(N)^2 = a + b/N IDENTIFIED BY THE DECISION ARM?")

set.seed(20260916L)
ident <- do.call(rbind, lapply(LANGS, function(L) do.call(rbind, lapply(names(TERM), function(h)
  do.call(rbind, lapply(BFS, function(B) {
    z   <- zstar(B); th <- THETA[[L]][[h]]
    inv <- make_inverse(th, DIR[[h]], z)
    cc  <- calib[calib$language == L & calib$hypothesis == h & calib$bf == B, ]
    cc  <- cc[order(cc$n_participants), ]
    fl  <- floor_of(L, h)
    free_fit <- function(s) {
      if (any(!is.finite(s))) return(c(a = NA, b = NA))
      co <- coef(lm(y ~ x, data.frame(y = s^2, x = 1 / cc$n_participants)))
      c(a = unname(co[1]), b = unname(co[2]))
    }
    pt <- free_fit(cc$s)
    # Parametric bootstrap of the three binomial rates: whatever the point fit says,
    # this shows whether 120 to 144 replicates a cell can pin down a slope in 1/N at
    # all over a range of 70 to 100.
    bs <- t(replicate(N_BOOT, {
      p <- rbinom(3, cc$reps, cc$p_arm) / cc$reps
      free_fit(inv(p))
    }))
    ok <- stats::complete.cases(bs)
    b80 <- N_REF * (s_ref(L, h, B)^2 - fl^2)
    # The floor-anchored slope from all three cells at once, by binomial maximum
    # likelihood. A single cell can land high or low by Monte-Carlo error alone,
    # and b80 inherits that error in full; pooling the three cells does not. The
    # deviance against the saturated fit (2 df) says whether one slope in 1/N can
    # account for all three cells.
    k   <- round(cc$p_arm * cc$reps)
    nll <- function(lb) {
      p <- power_at(sqrt(fl^2 + exp(lb) / cc$n_participants), th, DIR[[h]], z)
      -sum(dbinom(k, cc$reps, pmin(pmax(p, 1e-12), 1 - 1e-12), log = TRUE))
    }
    rng <- c(log(1e-4), log(1e3))
    opt <- optimize(nll, rng)
    prof <- function(side) {
      lim <- if (side < 0) rng[1] else rng[2]
      g <- function(lb) nll(lb) - opt$objective - qchisq(0.95, 1) / 2
      if (g(lim) < 0) return(exp(lim))
      exp(uniroot(g, sort(c(opt$minimum, lim)))$root)
    }
    b_ml <- exp(opt$minimum)
    ll_sat <- sum(dbinom(k, cc$reps, k / cc$reps, log = TRUE))
    g2 <- 2 * (opt$objective + ll_sat)
    data.frame(block = "identification", language = L, hypothesis = h, bf = B,
               a_free = pt[["a"]], b_free = pt[["b"]],
               b_free_lo = if (any(ok)) quantile(bs[ok, "b"], 0.025, names = FALSE) else NA,
               b_free_hi = if (any(ok)) quantile(bs[ok, "b"], 0.975, names = FALSE) else NA,
               p_boot_b_nonpositive = if (any(ok)) mean(bs[ok, "b"] <= 0) else NA,
               p_boot_a_below_zero = if (any(ok)) mean(bs[ok, "a"] < 0) else NA,
               boot_unusable = mean(!ok),
               a_anchored = fl^2, b_anchored = b80,
               anchored_identified = is.finite(b80) && b80 > 0,
               b_ml = b_ml, b_ml_lo = prof(-1), b_ml_hi = prof(1),
               s80_ml = sqrt(fl^2 + b_ml / N_REF),
               p80_ml = power_at(sqrt(fl^2 + b_ml / N_REF), th, DIR[[h]], z),
               ml_deviance = g2, ml_gof_p = pchisq(g2, df = 2, lower.tail = FALSE),
               # What the anchored model predicts at the other two sample sizes, set
               # against the arm with its own Monte-Carlo interval.
               p_pred_n70 = if (b80 > 0) power_at(sqrt(fl^2 + b80 / 70), th, DIR[[h]], z) else NA,
               p_arm_n70 = cc$p_arm[cc$n_participants == 70],
               p_arm_n70_lo = cc$p_arm_lo[cc$n_participants == 70],
               p_arm_n70_hi = cc$p_arm_hi[cc$n_participants == 70],
               p_pred_n100 = if (b80 > 0) power_at(sqrt(fl^2 + b80 / 100), th, DIR[[h]], z) else NA,
               p_arm_n100 = cc$p_arm[cc$n_participants == 100],
               p_arm_n100_lo = cc$p_arm_lo[cc$n_participants == 100],
               p_arm_n100_hi = cc$p_arm_hi[cc$n_participants == 100])
  }))))))
ident$pred_n70_inside_ci  <- with(ident, p_pred_n70 >= p_arm_n70_lo & p_pred_n70 <= p_arm_n70_hi)
ident$pred_n100_inside_ci <- with(ident, p_pred_n100 >= p_arm_n100_lo & p_pred_n100 <= p_arm_n100_hi)
print(ident[ident$bf == 10, c("language", "hypothesis", "a_free", "b_free", "b_free_lo", "b_free_hi",
                              "p_boot_b_nonpositive", "a_anchored", "b_anchored",
                              "p_pred_n70", "p_arm_n70", "pred_n70_inside_ci",
                              "p_pred_n100", "p_arm_n100", "pred_n100_inside_ci")],
      row.names = FALSE, digits = 3)
cat("\nFloor-anchored slope fitted to all three cells (maximum likelihood):\n")
print(ident[, c("language", "hypothesis", "bf", "b_anchored", "b_ml", "b_ml_lo", "b_ml_hi",
                "s80_ml", "p80_ml", "ml_deviance", "ml_gof_p")], row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# 5. The same model from the DGP variance components
# ---------------------------------------------------------------------------
# An independent route to b that does not touch the simulations. Per participant,
# the interaction is a difference of two affectedness slopes (Active and Passive
# rows), each estimated from 72 verbs with a per-row information I(eta) set by the
# cumulative-logit likelihood; averaging over N participants adds the participant
# slope variance. So
#
#     s(N)^2 = tau_verb^2 / Sxx  +  b / N,
#     b = 1 / E[ 1 / (Sigma_part_focal + v_p) ]
#
# with v_p = 1/Sxx_w(Active) + 1/Sxx_w(Passive) for H1b and 1/Sxx_w(Passive) for
# H1a, Sxx_w the information-weighted spread of affectedness for participant p.
# The harmonic mean is what GLS weighting over participants delivers: a participant
# pinned at the top of the scale carries almost no information and should count
# for almost nothing. Two INDEPENDENT versions double every I, halving v_p, which
# gives the most favourable r the fitted model allows. The pilot model has no
# participant-by-item term, so whatever the two versions share in a participant's
# reaction to a verb is treated there as independent noise, which pushes the real
# r towards 1. The one thing that could push the other way is also absent from the
# model: a verb-by-version item variance, which a second version would average
# down at the verb level. It matters only if a one-version session shows every
# participant the same version of a verb, and not if versions are counterbalanced.
rule("DGP VARIANCE COMPONENTS, AND r IF THE TWO VERSIONS WERE INDEPENDENT")

ordinal_info <- function(eta, thr) {
  tt <- c(-Inf, thr, Inf)
  Fm <- vapply(tt, function(t) plogis(t - eta), numeric(length(eta)))
  fm <- vapply(tt, function(t) dlogis(t - eta), numeric(length(eta)))
  P  <- Fm[, -1, drop = FALSE] - Fm[, -ncol(Fm), drop = FALSE]
  D  <- fm[, -ncol(fm), drop = FALSE] - fm[, -1, drop = FALSE]
  rowSums(D^2 / pmax(P, 1e-300))
}

set.seed(20260917L)
N_SIM <- 4000L
decomp <- do.call(rbind, lapply(LANGS, function(L) {
  d  <- DGP[[L]]; x <- unname(d$verb_affectedness); V <- length(x)
  ff <- d$fixef; SP <- d$Sigma_part; SV <- d$Sigma_verb
  bP <- MASS::mvrnorm(N_SIM, rep(0, ncol(SP)), SP); colnames(bP) <- colnames(SP)
  # Verb effects redrawn per simulated participant: this averages the information
  # over verb samples, which is what the decision arm does across replicates.
  bV <- MASS::mvrnorm(N_SIM * V, rep(0, ncol(SV)), SV); colnames(bV) <- colnames(SV)
  X  <- matrix(x, N_SIM, V, byrow = TRUE)
  eta_P <- ff[["Semantics_scaled"]] * X + bP[, "Intercept"] + bP[, "Semantics_scaled"] * X +
    matrix(bV[, "Intercept"], N_SIM, V)
  eta_A <- eta_P + ff[["S_TypeActive"]] + ff[["S_TypeActive:Semantics_scaled"]] * X +
    bP[, "S_TypeActive"] + bP[, "S_TypeActive:Semantics_scaled"] * X +
    matrix(bV[, "S_TypeActive"], N_SIM, V)
  info <- function(eta) matrix(ordinal_info(as.vector(eta), d$thresholds), N_SIM, V)
  sxx_w <- function(W) { xb <- rowSums(W * X) / rowSums(W); rowSums(W * (X - xb)^2) }
  vP <- 1 / sxx_w(info(eta_P)); vA <- 1 / sxx_w(info(eta_A))
  sxx <- sum((x - mean(x))^2)
  do.call(rbind, lapply(names(TERM), function(h) {
    sp_f <- SP[TERM[[h]], TERM[[h]]]
    v    <- if (h == "h1b") vA + vP else vP
    tau  <- sqrt(SV[if (h == "h1b") "S_TypeActive" else "Intercept",
                    if (h == "h1b") "S_TypeActive" else "Intercept"])
    a    <- tau^2 / sxx
    b1   <- 1 / mean(1 / (sp_f + v))
    b2   <- 1 / mean(1 / (sp_f + v / 2))
    th   <- THETA[[L]][[h]]
    s80c <- s_ref(L, h, 10)
    b80  <- N_REF * (s80c^2 - a)
    # Pilot sample sizes (claps_pilot_harmonised.csv: 70 English, 70 Turkish,
    # 57 Norwegian participants).
    n_pilot <- if (L == "Norwegian") 57L else 70L
    data.frame(block = "decomposition", language = L, hypothesis = h, bf = 10,
               tau_verb = tau, sxx = sxx, a_dgp = a,
               sigma_part_focal = sp_f, resid_part_one = b1 - sp_f,
               b_dgp_one = b1, b_dgp_two_indep = b2,
               share_resid_in_b = (b1 - sp_f) / b1,
               s80_dgp_one = sqrt(a + b1 / N_REF),
               s80_dgp_two_indep = sqrt(a + b2 / N_REF),
               s80_calibrated = s80c,
               r_indep_dgp = sqrt((a + b2 / N_REF) / (a + b1 / N_REF)),
               # The same bound with the calibrated participant term, scaled by the
               # DGP's own two-to-one ratio. It uses the calibration for the size of
               # the participant term and the DGP only for how much of it is residual.
               r_indep_anchored = if (b80 > 0) sqrt((a + b80 * (b2 / b1) / N_REF) / s80c^2) else NA,
               # The most the second version can add under the pilot model, on the
               # calibrated N = 80 baseline.
               gain_pp_at_r_indep_dgp = 100 * (power_at(sqrt((a + b2 / N_REF) / (a + b1 / N_REF)) * s80c,
                                                        th, DIR[[h]], zstar(10)) -
                                                 power_at(s80c, th, DIR[[h]], zstar(10))),
               power_one_dgp = power_at(sqrt(a + b1 / N_REF), th, DIR[[h]], zstar(10)),
               # A diagnostic, not an estimate of r: the posterior SD from the real
               # pilot (both versions, N = 70 or 57), next to what the DGP route
               # predicts for one and for two independent versions at that N.
               n_pilot = n_pilot, pilot_post_sd_two_versions = sd(th),
               s_pilot_n_dgp_one = sqrt(a + b1 / n_pilot),
               s_pilot_n_dgp_two_indep = sqrt(a + b2 / n_pilot),
               s_pilot_n_anchored_one = if (b80 > 0) sqrt(a + b80 / n_pilot) else NA)
  }))
}))
print(decomp[, c("language", "hypothesis", "a_dgp", "sigma_part_focal", "resid_part_one",
                 "share_resid_in_b", "s80_dgp_one", "s80_calibrated", "r_indep_dgp",
                 "r_indep_anchored", "gain_pp_at_r_indep_dgp", "pilot_post_sd_two_versions", "s_pilot_n_dgp_one",
                 "s_pilot_n_dgp_two_indep", "s_pilot_n_anchored_one")],
      row.names = FALSE, digits = 3)

# Correlation of the two focal estimates under the DGP, for the joint figure. The
# Passive slope is shared by both terms (entering H1b with a minus sign), and the
# verb intercepts and Active deviations are correlated in Sigma_verb.
focal_corr <- function(L, n, versions = 1) {
  d  <- DGP[[L]]; x <- unname(d$verb_affectedness); sxx <- sum((x - mean(x))^2)
  dc <- decomp[decomp$language == L, ]
  ra <- dc$resid_part_one[dc$hypothesis == "h1a"] / versions
  rb <- dc$resid_part_one[dc$hypothesis == "h1b"] / versions
  va <- d$Sigma_verb["Intercept", "Intercept"] / sxx + (d$Sigma_part[TERM[["h1a"]], TERM[["h1a"]]] + ra) / n
  vb <- d$Sigma_verb["S_TypeActive", "S_TypeActive"] / sxx + (d$Sigma_part[TERM[["h1b"]], TERM[["h1b"]]] + rb) / n
  cv <- d$Sigma_verb["Intercept", "S_TypeActive"] / sxx + (d$Sigma_part[TERM[["h1a"]], TERM[["h1b"]]] - ra) / n
  cv / sqrt(va * vb)
}

# Bivariate normal CDF by a one-dimensional integral on the probability scale,
# which needs nothing beyond base R and is exact to the rule's resolution here.
pbvn <- function(u, v, c, m = 100L) {
  if (abs(c) < 1e-12) return(pnorm(u) * pnorm(v))
  w  <- (seq_len(m) - 0.5) / m
  Pu <- pmax(pnorm(u), 1e-300)
  xs <- qnorm(outer(Pu, w))
  Pu * rowMeans(pnorm((v - c * xs) / sqrt(1 - c^2)))
}

# The correlated variant costs a hundred times the independent one, so it is
# computed at the preregistered threshold only.
joint_at <- function(L, sa, sb, B, corr) {
  z <- zstar(B)
  # With estimate = theta + s Z, H1a clears when -Z_a <= u and H1b when Z_b <= v,
  # so the pair (-Z_a, Z_b) carries correlation -corr.
  u <- THETA[[L]]$h1a / sa - z
  v <- -THETA[[L]]$h1b / sb - z
  c(indep = mean(pnorm(u) * pnorm(v)),
    dgp_corr = if (B == 10) mean(pbvn(u, v, -corr)) else NA_real_)
}

# ---------------------------------------------------------------------------
# 6. The r grid
# ---------------------------------------------------------------------------
rule("POWER OF THE TWO-VERSION DESIGN IF ITS SE IS r TIMES THE ONE-VERSION SE (N = 80)")

grid <- do.call(rbind, lapply(LANGS, function(L) do.call(rbind, lapply(BFS, function(B) {
  z <- zstar(B)
  single <- do.call(rbind, lapply(names(TERM), function(h) {
    th <- THETA[[L]][[h]]; s1 <- s_ref(L, h, B); fl <- floor_of(L, h)
    cr <- calib[calib$language == L & calib$hypothesis == h & calib$bf == B &
                  calib$n_participants == N_REF, ]
    p1 <- power_at(s1, th, DIR[[h]], z)
    p2 <- power_at(R_GRID * s1, th, DIR[[h]], z)
    idr <- ident[ident$language == L & ident$hypothesis == h & ident$bf == B, ]
    dcr <- decomp[decomp$language == L & decomp$hypothesis == h, ]
    # N of the ONE-version design with the same precision: a + b/N = r^2 s80^2.
    # No finite N exists once r s80 is at or below the verb floor; Inf says so.
    n_eq <- function(a, b, s80) {
      den <- R_GRID^2 * s80^2 - a
      ifelse(b > 0 & den > 0, b / den, Inf)
    }
    tidy0 <- function(v) ifelse(is.finite(v) & abs(v) < 1e-6, 0, v)
    ne_anch <- if (isTRUE(idr$anchored_identified)) n_eq(fl^2, idr$b_anchored, s1) else NA
    s80_dgp <- sqrt(dcr$a_dgp + dcr$b_dgp_one / N_REF)
    ne_dgp  <- n_eq(dcr$a_dgp, dcr$b_dgp_one, s80_dgp)
    s_ml    <- idr$s80_ml
    ne_ml   <- n_eq(fl^2, idr$b_ml, s_ml)
    data.frame(block = "grid", language = L, hypothesis = h, bf = B,
               n_participants = N_REF, r = R_GRID,
               p_arm = cr$p_arm, mcse_arm = cr$mcse_arm, s_one = s1, s_two = R_GRID * s1,
               se_floor = fl, calibration_flag = cr$flag,
               power_one = p1, power_two = p2, gain_pp = 100 * (p2 - p1),
               mcse_arm_pp = 100 * cr$mcse_arm,
               # if/else, not ifelse(): a scalar test would return one value and
               # recycle it down the whole r grid.
               gain_over_mcse = if (cr$mcse_arm > 0) (p2 - p1) / cr$mcse_arm else NA_real_,
               n_one_version_equiv_anchored = ne_anch,
               extra_participants_anchored = tidy0(ne_anch - N_REF),
               n_one_version_equiv_dgp = ne_dgp,
               extra_participants_dgp = tidy0(ne_dgp - N_REF),
               # Sensitivity: the same grid from the three-cell fit, whose baseline
               # does not lean on the N = 80 cell alone.
               s_one_ml = s_ml, power_one_ml = power_at(s_ml, th, DIR[[h]], z),
               power_two_ml = power_at(R_GRID * s_ml, th, DIR[[h]], z),
               gain_pp_ml = 100 * (power_at(R_GRID * s_ml, th, DIR[[h]], z) -
                                     power_at(s_ml, th, DIR[[h]], z)),
               extra_participants_ml = tidy0(ne_ml - N_REF),
               # Participant-minutes: 80 people at 40 minutes against the one-version
               # N that matches their precision at 20 minutes each.
               minutes_two_version = N_REF * MIN_TWO,
               minutes_one_version_equiv_anchored = ne_anch * MIN_ONE)
  }))
  sa <- s_ref(L, "h1a", B); sb <- s_ref(L, "h1b", B)
  j1 <- joint_at(L, sa, sb, B, focal_corr(L, N_REF, 1))
  joint <- do.call(rbind, lapply(R_GRID, function(r) {
    # The same r is applied to both terms and the one-version correlation is kept.
    # The log prints the correlation for two independent versions beside it; the
    # two differ little, so holding it fixed moves the joint gain very little.
    j2 <- joint_at(L, r * sa, r * sb, B, focal_corr(L, N_REF, 1))
    data.frame(r = r, power_one = j1[["indep"]], power_two = j2[["indep"]],
               power_one_corr = j1[["dgp_corr"]], power_two_corr = j2[["dgp_corr"]])
  }))
  jp <- arm_get(L, N_REF, paste0("p_joint", if (B == 10) "" else paste0("_bf", B)))
  joint <- data.frame(block = "grid", language = L, hypothesis = "joint", bf = B,
                      n_participants = N_REF, joint,
                      p_arm = jp[["p"]], mcse_arm = sqrt(jp[["p"]] * (1 - jp[["p"]]) / jp[["reps"]]),
                      gain_pp = 100 * (joint$power_two - joint$power_one),
                      gain_pp_corr = 100 * (joint$power_two_corr - joint$power_one_corr),
                      mcse_arm_pp = 100 * sqrt(jp[["p"]] * (1 - jp[["p"]]) / jp[["reps"]]),
                      joint_dependence = "indep = independent given the pilot draw; corr = DGP-implied correlation of the two estimates")
  bind_rows(single, joint)
}))))

cat("\nFocal-estimate correlation implied by the DGP at N = 80:",
    paste(sprintf("%s %.2f (one version) %.2f (two independent)", LANGS,
                  vapply(LANGS, focal_corr, numeric(1), n = N_REF, versions = 1),
                  vapply(LANGS, focal_corr, numeric(1), n = N_REF, versions = 2)), collapse = "; "), "\n")

show <- grid[grid$bf == 10 & grid$r %in% c(R_KEY, 1) & grid$hypothesis != "joint",
             c("language", "hypothesis", "r", "p_arm", "mcse_arm", "s_one", "power_two",
               "gain_pp", "gain_over_mcse", "extra_participants_anchored",
               "extra_participants_dgp", "calibration_flag")]
print(show, row.names = FALSE, digits = 3)
cat("\nJoint H1a-and-H1b, BF 10:\n")
print(grid[grid$bf == 10 & grid$r %in% c(R_KEY, 1) & grid$hypothesis == "joint",
           c("language", "r", "p_arm", "power_one", "power_two", "gain_pp",
             "power_one_corr", "power_two_corr", "gain_pp_corr")],
      row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# 7. The pooled K-language test
# ---------------------------------------------------------------------------
# The K-language engine is reused, not re-implemented: its meta-analytic posterior
# and Bayes factor are called directly, and the per-language SE split is rebuilt
# from the same inputs scripts/run_klanguage_design_analysis.R uses, with a check
# that it lands on the SE that script recorded. The replicate loop is written out
# here only so that the population draw, the true effects and the standard-normal
# noise are SHARED across r. The difference in power is then a paired contrast,
# whose Monte-Carlo error is far smaller than that of two independent cells.
rule("POOLED K-LANGUAGE TEST UNDER THE SAME r")

source("R/13_simulate_klanguage_pooled.R")
H <- read.csv(file.path(OUT, "cross_language_effects_harmonised.csv"))
SE_VERB <- mean(PP$se_floor_h1b  * PP$affectedness_sd)
SE_POST <- mean(PP$slope_post_sd * PP$affectedness_sd)
SE_PPT  <- sqrt(max(SE_POST^2 - SE_VERB^2, 0))
SE_LANG <- sqrt(SE_VERB^2 + SE_PPT^2)
RHO_EFF <- 0.164
SE_SHARED <- SE_VERB * sqrt(RHO_EFF)
SE_INDEP  <- sqrt(SE_VERB^2 * (1 - RHO_EFF) + SE_PPT^2)
kp <- read.csv(file.path(OUT, "klanguage_power.csv"))
if (abs(SE_INDEP - kp$se_lang[1]) > 1e-9 || abs(SE_SHARED - kp$shared_verb_sd[1]) > 1e-9)
  stop("per-language SE split has drifted from klanguage_power.csv; re-read run_klanguage_design_analysis.R")

Hc <- H[H$se_kind %in% c("slope_sd", "posterior_sd"), ]
sdx <- setNames(PP$affectedness_sd, PP$language)
hit <- Hc$source == "CLAPS_pilot" & Hc$language %in% LANGS
Hc$se[hit] <- vapply(Hc$language[hit], function(L) sd(THETA[[L]]$h1b) * sdx[[L]], numeric(1))
Hc <- Hc[!duplicated(paste(Hc$source, Hc$language)), ]
one_per_language <- function(df) {
  df <- df[order(df$source != "CLAPS_pilot"), ]
  df[!duplicated(df$language), ]
}
POPS <- list(claps_protocol  = subset(Hc, source == "CLAPS_pilot"),
             seven           = one_per_language(Hc),
             glossa_protocol = subset(Hc, source == "Glossa2023" & language != "English"))
POPS <- lapply(POPS, function(d) data.frame(language = d$language, y = d$y, s = d$se))

R_POOL <- c(1, 0.95, 0.90, 0.85)
# Two ways the within-language change can reach the pooled test. "all" scales the
# whole within-language SE, shared verb part included, and so overstates the gain:
# no change to the session can alter which 72 verbs the teams translate. "floor"
# holds the verb part fixed and asks the participant part to absorb all of r, and
# is undefined once r * SE_LANG falls below the verb part.
split_for <- function(r, variant) {
  if (variant == "all") return(c(indep = r * SE_INDEP, shared = r * SE_SHARED))
  ppt2 <- r^2 * SE_LANG^2 - SE_VERB^2
  if (ppt2 < 0) return(c(indep = NA, shared = NA))
  c(indep = sqrt(SE_VERB^2 * (1 - RHO_EFF) + ppt2), shared = SE_SHARED)
}

pooled <- do.call(rbind, lapply(names(POPS), function(pn) do.call(rbind, lapply(c(3L, 6L, 12L), function(K) {
  set.seed(5100L + K + 10L * match(pn, names(POPS)))
  pop  <- POPS[[pn]]
  pars <- draw_population_params(POOLED_REPS, pop$y, pop$s, prior_tau_scale = DESIGN_TAU_PRIOR_SCALE)
  combos <- expand.grid(r = R_POOL, variant = c("all", "floor"), stringsAsFactors = FALSE)
  combos <- combos[!(combos$r == 1 & combos$variant == "floor"), ]
  hitm <- matrix(NA, POOLED_REPS, nrow(combos))
  for (i in seq_len(POOLED_REPS)) {
    theta <- rnorm(K, pars$mu[i], pars$tau[i])
    e0 <- rnorm(1); e <- rnorm(K)
    for (j in seq_len(nrow(combos))) {
      sp <- split_for(combos$r[j], combos$variant[j])
      if (anyNA(sp)) next
      y  <- theta + sp[["shared"]] * e0 + sp[["indep"]] * e
      post <- ra_meta_posterior(y, rep(sp[["indep"]], K), prior_tau_scale = CALIBRATED_TAU_PRIOR_SCALE)
      hitm[i, j] <- directional_bf_meta(post) >= 10
    }
  }
  base <- hitm[, which(combos$r == 1)]
  do.call(rbind, lapply(seq_len(nrow(combos)), function(j) {
    sp <- split_for(combos$r[j], combos$variant[j])
    d  <- hitm[, j] - base
    data.frame(block = "pooled", regime = pn, K = K, bf = 10, r = combos$r[j],
               pooled_variant = if (combos$r[j] == 1) "baseline" else combos$variant[j],
               reps = POOLED_REPS, se_indep = sp[["indep"]], se_shared = sp[["shared"]],
               power_one = mean(base), power_two = mean(hitm[, j]),
               gain_pp = 100 * mean(d), mcse_gain_pp = 100 * sd(d) / sqrt(POOLED_REPS),
               mcse_power_pp = 100 * sqrt(mean(base) * (1 - mean(base)) / POOLED_REPS))
  }))
}))))
print(pooled[, c("regime", "K", "r", "pooled_variant", "se_indep", "power_one", "power_two",
                 "gain_pp", "mcse_gain_pp")], row.names = FALSE, digits = 3)

# The mechanism, in closed form. With known tau the pooled SE is
# sqrt((tau^2 + s_indep^2)/K + s_shared^2); the ratio below is how far r moves it.
# The tau values are each regime's posterior median under the design prior, and
# a second within-language SE, from the calibrated s80 for H1b on the per-SD
# scale, shows the bound under the larger SE the brms arm implies.
tau_med <- vapply(POPS, function(p)
  ra_meta_posterior(p$y, p$s, prior_tau_scale = DESIGN_TAU_PRIOR_SCALE)$tau_median, numeric(1))
s80_sd <- mean(vapply(LANGS, function(L) s_ref(L, "h1b", 10) * sdx[[L]], numeric(1)))
analytic <- do.call(rbind, lapply(names(tau_med), function(pn) do.call(rbind, lapply(c(3L, 6L, 12L), function(K)
  do.call(rbind, lapply(c("engine", "calibrated_s80"), function(src) {
    si <- if (src == "engine") SE_INDEP else sqrt(max(s80_sd^2 - SE_SHARED^2, 0))
    ss <- SE_SHARED
    se0 <- sqrt((tau_med[[pn]]^2 + si^2) / K + ss^2)
    data.frame(block = "pooled", regime = pn, K = K, bf = 10, r = R_KEY,
               pooled_variant = paste0("analytic_all_", src),
               tau = tau_med[[pn]], se_indep = si, se_shared = ss,
               pooled_se_ratio = sqrt((tau_med[[pn]]^2 + R_KEY^2 * si^2) / K + R_KEY^2 * ss^2) / se0)
  }))))))
cat("\nAnalytic pooled-SE ratio when the whole within-language SE shrinks by r (an upper bound):\n")
print(analytic[analytic$K %in% c(6L, 12L), c("regime", "K", "r", "pooled_variant", "tau", "se_indep",
                                             "pooled_se_ratio")], row.names = FALSE, digits = 3)

# ---------------------------------------------------------------------------
# 8. Write
# ---------------------------------------------------------------------------
res <- bind_rows(calib, ident, decomp, grid, pooled, analytic)
lead <- c("block", "language", "hypothesis", "regime", "K", "bf", "n_participants", "r")
res <- res[, c(lead, setdiff(names(res), lead))]
write.csv(res, file.path(OUT, "gender_versions_power.csv"), row.names = FALSE)
cat("\nwrote", file.path(OUT, "gender_versions_power.csv"), "(", nrow(res), "rows )\n")

rule("HEADLINE: H1b AT BF 10, N = 80")
hl <- grid[grid$bf == 10 & grid$hypothesis == "h1b" & grid$r %in% c(R_KEY, 1),
           c("language", "r", "power_two", "gain_pp", "mcse_arm_pp", "extra_participants_anchored",
             "power_one_ml", "gain_pp_ml", "extra_participants_ml")]
print(hl, row.names = FALSE, digits = 3)
