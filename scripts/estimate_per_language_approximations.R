#!/usr/bin/env Rscript
# scripts/estimate_per_language_approximations.R
# ---------------------------------------------------------------------------
# Two per-language quantities the simulations do not give directly, by a normal
# approximation to each language's directional Bayes factor:
#
#   1. The expected detection rate in a language WITHOUT a pilot, at 80 participants
#      and 72 verbs, when its effect is drawn from the population of languages that
#      the 20-language synthesis assumes. Thirteen of the 20 planned languages have
#      no pilot, so the decision arm, which simulates only English, Turkish and
#      Norwegian, cannot speak for them.
#   2. The rate of misleading evidence: how often a language reaches a Bayes factor
#      of 10 in the predicted direction when the true effect is exactly zero.
#
# THE APPROXIMATION
#   A language's estimate is normal around its true effect theta with standard
#   error s, the prior on the effect is normal around zero with SD P (the priors of
#   the Design-Analysis Method, per SD of affectedness), so the posterior is normal
#   with mean w * est and SD sqrt(w) * s, where w = P^2 / (P^2 + s^2). The
#   directional Bayes factor reaches 10 exactly when the posterior probability of
#   the predicted sign reaches 10/11 (symmetric prior), that is when
#   w * est / (sqrt(w) * s) <= -qnorm(10/11) for a negative prediction.
#
# CALIBRATION
#   The full model is wider than this normal-normal sketch, so s is first scaled by
#   one factor k, chosen so that the approximation reproduces the decision arm's
#   own rates at 80 participants for the pilot languages, with each language's
#   true effect drawn from its pilot posterior as the decision arm draws it. The
#   same k then scales the standard error a language without a pilot is given in
#   the 20-language engine.
#
# Outputs: outputs/design_summary_pilot/per_language_approximations.csv
# Run from design_analysis/:  Rscript scripts/estimate_per_language_approximations.R
# ---------------------------------------------------------------------------

source("R/13_simulate_klanguage_pooled.R")
IN <- "outputs/design_summary_pilot"
PP  <- read.csv(file.path(IN, "pilot_params_ceilings.csv"))
PPD <- read.csv(file.path(IN, "pseudo_passive_decision_arm.csv"))
H   <- read.csv(file.path(IN, "cross_language_effects_harmonised.csv"))
N_PLAN <- 80
Z10 <- qnorm(10 / 11)
# Prior SDs on the model's scale (R/03_define_priors.R, primary regime); per SD of
# affectedness they are halved with the 2-SD scaling of Semantics_scaled.
PRIOR_MODEL <- c(h1 = 0.5, h2 = 0.6)

pilot <- function(L) readRDS(file.path("outputs/pilot_models", paste0("pilot_dgp_v2_pilot_", L, ".rds")))
N_PILOT <- c(English = 70, Turkish = 70, Norwegian = 57)
TERM <- c(h1 = "S_TypeActive:Semantics_scaled", h2 = "S_TypePseudo_Passive:Semantics_scaled")
VERB_IX <- c(h1 = 2L, h2 = 3L)  # Sigma_verb rows: intercept, Active, Pseudo_Passive

# One language and hypothesis: draws of the true effect and the SE at 80 participants,
# per SD of affectedness, with the participant part rescaled from the pilot sample.
lang_inputs <- function(L, h) {
  p <- pilot(L)
  sdx <- PP$affectedness_sd[PP$language == L]
  sxx <- PP$sxx[PP$language == L]
  d <- as.numeric(as.matrix(p$fixef_draws)[, TERM[[h]]]) * sdx
  s_verb <- sqrt(p$Sigma_verb[VERB_IX[[h]], VERB_IX[[h]]] / sxx) * sdx
  s_ppt  <- sqrt(max(sd(d)^2 - s_verb^2, 0)) * sqrt(N_PILOT[[L]] / N_PLAN)
  list(theta = d, s = sqrt(s_verb^2 + s_ppt^2), P = PRIOR_MODEL[[h]] * sdx)
}

# P(BF >= 10 for a negative effect), averaged over draws of the true effect.
p_detect <- function(theta, s, P) {
  w <- P^2 / (P^2 + s^2)
  mean(pnorm((-Z10 * s / sqrt(w) - theta) / s))
}

cells <- rbind(data.frame(language = c("English", "Turkish", "Norwegian"), h = "h1"),
               data.frame(language = c("English"), h = "h2"))
cells$target <- mapply(function(L, h) {
  r <- PPD[PPD$language == L & PPD$n_participants == N_PLAN, ]
  r[[paste0("n_", if (h == "h1") "h1b" else "h2b")]] / r$reps
}, cells$language, cells$h)
inp <- Map(lang_inputs, cells$language, cells$h)

# One k for all cells, by least squares on the logit scale of the rates.
loss <- function(k) {
  pr <- vapply(inp, function(x) p_detect(x$theta, k * x$s, x$P), numeric(1))
  sum((qlogis(pmin(pmax(pr, 1e-4), 1 - 1e-4)) - qlogis(cells$target))^2)
}
k <- optimize(loss, c(0.5, 4))$minimum
cells$s <- vapply(inp, `[[`, numeric(1), "s")
cells$approx_k1 <- vapply(inp, function(x) p_detect(x$theta, x$s, x$P), numeric(1))
cells$approx_k  <- vapply(inp, function(x) p_detect(x$theta, k * x$s, x$P), numeric(1))
cat(sprintf("calibration factor k = %.3f\n", k))
print(cells, row.names = FALSE, digits = 3)

# --- 1. A language without a pilot ---------------------------------------------
# Its SE is the engine's per-language SE at 80 participants and 72 verbs, the same
# as in the 20-language runs, scaled by k. Its true effect comes from the
# population posterior of each set of measured languages, as for a new place there.
RHO_EFF <- 0.164
SE_VERB <- mean(PP$se_floor_h1b * PP$affectedness_sd)
SE_POST <- mean(PP$slope_post_sd * PP$affectedness_sd)
SE_H1 <- sqrt(SE_VERB^2 + max(SE_POST^2 - SE_VERB^2, 0))
pp_in <- lapply(c("English", "Turkish"), lang_inputs, h = "h2")
SE_H2 <- {
  fl <- vapply(c("English", "Turkish"), function(L) {
    p <- pilot(L); sdx <- PP$affectedness_sd[PP$language == L]
    sqrt(p$Sigma_verb[3, 3] / PP$sxx[PP$language == L]) * sdx
  }, numeric(1))
  po <- vapply(pp_in, function(x) sd(x$theta), numeric(1))
  sqrt(mean(fl)^2 + max(mean(po)^2 - mean(fl)^2, 0) * mean(N_PILOT[c("English", "Turkish")]) / N_PLAN)
}
claps_se <- setNames(vapply(c("English", "Turkish", "Norwegian"), function(L) sd(lang_inputs(L, "h1")$theta),
                            numeric(1)), c("English", "Turkish", "Norwegian"))
P1 <- build_populations(prepare_language_table(H, claps_se), "canonical")
g <- H[H$source == "Glossa2023", ]
known5 <- rbind(
  data.frame(language = c("English", "Turkish"), source = "CLAPS_pilot",
             y = vapply(pp_in, function(x) mean(x$theta), numeric(1)),
             s = vapply(pp_in, function(x) sd(x$theta), numeric(1))),
  data.frame(language = c("Balinese", "Indonesian", "Mandarin"), source = "Glossa2023",
             y = c(-0.53, -0.34, -1.09) * g$sd_x[match(c("Balinese", "Indonesian", "Mandarin"), g$language)],
             s = c(0.13, 0.10, 0.18) * g$sd_x[match(c("Balinese", "Indonesian", "Mandarin"), g$language)]))
P2 <- list(seven = known5, glossa_protocol = known5[known5$source == "Glossa2023", ],
           claps_protocol = known5[known5$source == "CLAPS_pilot", ])
P_H1 <- PRIOR_MODEL[["h1"]] * mean(PP$affectedness_sd)
P_H2 <- PRIOR_MODEL[["h2"]] * mean(PP$affectedness_sd)

set.seed(20261007L)
NDRAW <- 200000L
no_pilot <- do.call(rbind, lapply(c("claps_protocol", "seven", "glossa_protocol"), function(pn) {
  d1 <- draw_population_params(NDRAW, P1[[pn]]$y, P1[[pn]]$se %||% P1[[pn]]$s)
  d2 <- draw_population_params(NDRAW, P2[[pn]]$y, P2[[pn]]$s)
  th1 <- rnorm(NDRAW, d1$mu, d1$tau)
  th2 <- rnorm(NDRAW, d2$mu, d2$tau)
  data.frame(quantity = "detection, language without a pilot", population = pn,
             h1 = p_detect(th1, k * SE_H1, P_H1), h2 = p_detect(th2, k * SE_H2, P_H2))
}))

# --- 2. Misleading evidence at a true effect of zero -----------------------------
zero <- rbind(
  data.frame(quantity = "BF >= 10 at a true effect of zero", population = c("English", "Turkish", "Norwegian"),
             h1 = vapply(c("English", "Turkish", "Norwegian"), function(L) {
               x <- lang_inputs(L, "h1"); p_detect(0, k * x$s, x$P) }, numeric(1)),
             h2 = c(vapply(c("English", "Turkish"), function(L) {
               x <- lang_inputs(L, "h2"); p_detect(0, k * x$s, x$P) }, numeric(1)), NA)),
  data.frame(quantity = "BF >= 10 at a true effect of zero", population = "language without a pilot",
             h1 = p_detect(0, k * SE_H1, P_H1), h2 = p_detect(0, k * SE_H2, P_H2)))

out <- rbind(no_pilot, zero)
out$k <- k
cat("\n"); print(out, row.names = FALSE, digits = 3)
write.csv(out, file.path(IN, "per_language_approximations.csv"), row.names = FALSE)
