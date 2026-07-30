# R/03_define_priors.R
# Prespecified prior regimes for the CLAPS Bayesian cumulative-logit
# mixed-effects model.
#
# Literature basis (all sources verified in references.bib):
#
#   * Ordinal cumulative-logit model and brms thresholds:
#       Bürkner & Vuorre (2019); Liddell & Kruschke (2018); Bürkner (2017, 2018);
#       Veríssimo (2021); Taylor et al. (2023).
#
#   * General principles for proper, weakly-to-moderately informative priors
#     in the context of the ordinal likelihood and Bayes factors:
#       Gelman et al. (2008, 2017); Schad et al. (2021, 2023);
#       Gabry et al. (2019).
#
#   * Variance-component priors in hierarchical models:
#       Gelman (2006); Chung et al. (2015); Simpson et al. (2017).
#
#   * LKJ correlation prior:
#       Lewandowski et al. (2009).
#
#   * Substantive direction and magnitude of the affectedness slope and the
#     active vs passive interaction, used to scale the focal-slope priors:
#       Ambridge et al. (2016, 2023); Aryawibawa & Ambridge (2018);
#       Darmasetiyawan & Ambridge (2022); Liu & Ambridge (2021);
#       Bidgood et al. (2020); Paolazzi et al. (2022).
#
#   * Acquisition evidence supporting a semantic constraint on the passive:
#       Maratsos et al. (1985); Pinker, Lebeaux & Frost (1987);
#       Nguyen & Pearl (2021); Agostinho, Gavarró & Santos (2025).
#
# Four prior regimes are defined:
#
#   1. primary             - prespecified, weakly-to-moderately informative,
#                              focal slopes zero-centred but scaled to
#                              previously observed magnitudes
#   2. weak                  - weak-prior sensitivity check
#   3. literature_centred    - sensitivity-only prior centred on pooled
#                              previously observed values (NOT used for
#                              primary Bayes factors because it encodes the
#                              direction)
#   4. heavy_tailed          - Student-t robustness check
#
# Two threshold modes:
#
#   * broad                  - generic Student-t(3, 0, 2.5) on each
#                              ordinal intercept (Bürkner & Vuorre, 2019)
#   * ceiling_calibrated     - per-threshold normal priors centred on
#                              logit cumulative proportions from
#                              independent pilot data (Schad et al., 2021)

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(purrr)
})

# ---------------------------------------------------------------------------
# Empirical anchors from previously published passive-affectedness studies.
# The passive Semantics slope and the active-vs-passive interaction are
# cross-language pooled posterior means from the Bayesian meta-analytic
# synthesis (Ambridge, Arnon & Bekman, 2023), which draws together the
# antecedent single-language studies (Ambridge et al., 2016; Aryawibawa &
# Ambridge, 2018; Darmasetiyawan & Ambridge, 2022; Liu & Ambridge, 2021).
#
# NOTE: the pseudo-passive interaction is NOT a cross-language pooled value.
# Pseudo-passives are attested in only a subset of the contributing languages,
# so the figure below is a directional anchor reflecting the meta-analytic
# conclusion that affectedness does not raise pseudo-passive acceptability
# (i.e. a negative pseudo-vs-passive interaction, estimated chiefly from the
# Hebrew data). It is used only by the direction-encoding literature_centred
# sensitivity regime, never for the primary Bayes factor. All magnitudes are
# reported by the antecedent studies on their original rating scales and are
# used here only as order-of-magnitude anchors for the log-odds-scale priors.
# ---------------------------------------------------------------------------

EMPIRICAL_ANCHORS <- list(
  semantics_pooled                  = 0.47,   # pooled (cross-language) passive Semantics slope
  s_type_active_interaction         = -0.31,  # pooled (cross-language) active - passive
  s_type_pseudo_passive_interaction = -0.36,  # directional anchor (NOT pooled); see note above
  # Single-language Semantics range reported in the antecedent studies:
  semantics_min = 0.27,
  semantics_max = 0.80
)

# ---------------------------------------------------------------------------
# Prior regimes
# ---------------------------------------------------------------------------
# Each regime is a list of brms prior strings indexed by:
#   * b_default          : non-focal fixed effects (sentence-type main effects
#                          and any other regression coefficients)
#   * b_semantics        : passive affectedness slope (focal)
#   * b_active_int       : S_TypeActive : Semantics interaction (primary test)
#   * b_pseudo_int       : S_TypePseudo_Passive : Semantics interaction
#                          (secondary test, less constrained direction)
#   * Intercept          : ordinal thresholds in broad mode
#   * sd                 : group-level standard deviations
#   * cor                : group-level correlation matrices
#
# Coefficient names assume treatment coding with Passive as the reference
# level. They should be confirmed with brms::get_prior() before fitting.

PRIOR_REGIMES <- list(

  # ---- Primary prior set ------------------------------------------------
  # Zero-centred focal-slope priors avoid pre-loading the directional
  # hypothesis tests; scales are calibrated to cover the previously
  # observed range of single-language Semantics slopes (~0.27 to 0.80) and
  # interactions (Ambridge et al., 2016; Ambridge, Arnon & Bekman, 2023;
  # Aryawibawa & Ambridge, 2018; Darmasetiyawan & Ambridge, 2022;
  # Liu & Ambridge, 2021). Non-focal fixed effects use the
  # weakly-informative logistic-regression prior tradition
  # (Gelman et al., 2008, 2017). Threshold, sd and cor priors follow
  # Bürkner & Vuorre (2019), Gelman (2006), Chung et al. (2015),
  # Simpson et al. (2017) and Lewandowski et al. (2009).
  #
  # How the scales were chosen. Every focal prior is centred at zero, so it does
  # not presuppose the direction under test, and the scale is set so that the
  # previously observed magnitudes sit comfortably inside the prior rather than in
  # its tail. With SD 0.5, roughly 95% of the prior mass lies within +/- 1.0 on
  # the log-odds scale, which spans the published single-language Semantics slopes
  # of about 0.27 to 0.80 (EMPIRICAL_ANCHORS) with room to spare. The
  # pseudo-passive interaction gets a slightly wider 0.6 because it is attested in
  # fewer languages and its magnitude is correspondingly less well established.
  #
  # A caution specific to Bayes factors: unlike a posterior mean, a Savage-Dickey
  # Bayes factor depends on the prior scale directly, and an arbitrarily wide prior
  # inflates evidence for the null (the Jeffreys-Lindley effect). That is why these
  # scales are argued from prior evidence rather than made vague for safety, and
  # why the `weak` regime below exists as a declared sensitivity check instead of
  # as the default.
  primary = list(
    b_default     = "normal(0, 1.5)",
    b_semantics   = "normal(0, 0.5)",
    b_active_int  = "normal(0, 0.5)",
    b_pseudo_int  = "normal(0, 0.6)",
    Intercept     = "student_t(3, 0, 2.5)",
    # Student-t rather than normal on the group-level SDs: the heavier tail lets a
    # variance component be larger than expected without the prior fighting the
    # data, while still shrinking it away from implausibly large values
    # (Gelman, 2006; Chung et al., 2015).
    sd            = "student_t(3, 0, 1)",
    # lkj(2) places mild mass away from perfect correlation, which regularises the
    # random-effects correlation matrix that L5 must estimate; lkj(1) would be
    # uniform over all valid matrices (Lewandowski et al., 2009).
    cor           = "lkj(2)"
  ),

  # ---- Weak-prior sensitivity ------------------------------------------
  # Replaces focal slope SDs with normal(0, 1) and non-focal SDs with
  # normal(0, 2). Conclusions should not depend on the more concentrated
  # primary prior (Schad et al., 2021, 2023).
  # Roughly a doubling of every focal scale, which is the magnitude that matters:
  # a Bayes factor that survives a halving of prior precision is not an artefact
  # of the prior. lkj(1) removes the mild correlation regularisation as part of the
  # same check.
  weak = list(
    b_default     = "normal(0, 2)",
    b_semantics   = "normal(0, 1)",
    b_active_int  = "normal(0, 1)",
    b_pseudo_int  = "normal(0, 1)",
    Intercept     = "student_t(3, 0, 2.5)",
    sd            = "student_t(3, 0, 2)",
    cor           = "lkj(1)"
  ),

  # ---- Literature-centred sensitivity ----------------------------------
  # Centres focal slopes on the previously observed values: pooled
  # cross-language means for the passive slope and the active interaction,
  # and a directional anchor for the pseudo-passive interaction (see
  # EMPIRICAL_ANCHORS). This regime is appropriate for estimation
  # sensitivity but is reported as a sensitivity-only check for Bayes
  # factors because it encodes the predicted direction (Schad et al., 2023).
  literature_centred = list(
    b_default     = "normal(0, 1.5)",
    b_semantics   = "normal(0.47, 0.35)",
    b_active_int  = "normal(-0.31, 0.4)",
    b_pseudo_int  = "normal(-0.36, 0.5)",
    Intercept     = "student_t(3, 0, 2.5)",
    sd            = "student_t(3, 0, 1)",
    cor           = "lkj(2)"
  ),

  # ---- Heavy-tailed robustness check ------------------------------------
  # Student-t priors on the focal slopes accommodate larger language-
  # specific effects without inflating prior mass at exactly zero
  # (Gelman et al., 2008).
  heavy_tailed = list(
    b_default     = "student_t(3, 0, 1.5)",
    b_semantics   = "student_t(3, 0, 0.5)",
    b_active_int  = "student_t(3, 0, 0.5)",
    b_pseudo_int  = "student_t(3, 0, 0.6)",
    Intercept     = "student_t(3, 0, 2.5)",
    sd            = "student_t(3, 0, 1)",
    cor           = "lkj(2)"
  )
)

# ---------------------------------------------------------------------------
# Threshold prior modes
# ---------------------------------------------------------------------------
# Generic mode places a single Student-t(3, 0, 2.5) prior on every ordinal
# intercept; this is wide enough to allow mass at either end of a 1-7
# scale while remaining proper (required for valid Bayes factors;
# Bürkner & Vuorre, 2019; Schad et al., 2023).
# Ceiling-calibrated mode places per-threshold normal priors centred on
# logit cumulative proportions from an independent pilot sample. This is
# a preregistered sensitivity analysis intended for languages with
# documented ceiling effects in the independent pilot data (Schad et al.,
# 2021).

THRESHOLD_MODES <- c("broad", "ceiling_calibrated")

#' Compute ceiling-calibrated threshold priors from independent pilot data.
#'
#' Uses smoothed cumulative category proportions transformed to the logit
#' scale. Smoothing prevents 0/1 cumulative probabilities at the
#' boundaries (Bürkner & Vuorre, 2019). The pilot sample must be
#' independent of the confirmatory sample; the split is enforced by
#' R/01_read_validate_data.R::split_pilot_confirmatory().
#'
#' @param pilot_df Pilot data frame containing a 1-7 Response column.
#' @param language Character; language label for logging.
#' @param smooth_alpha Dirichlet smoothing parameter.
#' @return Named list with threshold_means and threshold_sds vectors
#'   (length 6, one per cumulative threshold of a 7-point scale).
compute_ceiling_calibrated_thresholds <- function(pilot_df,
                                                  language     = "unknown",
                                                  smooth_alpha = 0.5) {
  stopifnot("Response" %in% names(pilot_df))
  n_cats <- 7L
  # tabulate() with nbins guarantees a length-7 vector, including zeros for
  # categories no participant used. table() would silently omit them and shift
  # every subsequent threshold.
  counts <- tabulate(pilot_df$Response, nbins = n_cats)
  # Additive (Dirichlet / Laplace) smoothing. An unused category would otherwise
  # give a cumulative proportion equal to its neighbour's, and a category at the
  # top of the scale would give exactly 1, whose logit is infinite. alpha = 0.5 is
  # the Jeffreys prior for a multinomial, and is set in
  # config/analysis_config.yaml under ceiling_calibration$smooth_alpha.
  smoothed <- (counts + smooth_alpha) / (sum(counts) + n_cats * smooth_alpha)
  # A 7-category ordinal model has 6 thresholds, so the final cumulative
  # probability (which is 1 by construction) is dropped.
  cum_prob <- cumsum(smoothed)[seq_len(n_cats - 1)]
  # Belt-and-braces bound before the logit. Smoothing already prevents exact 0 and
  # 1, but with a large pilot sample a genuine ceiling can still push a cumulative
  # proportion close enough to 1 that its logit becomes an extreme threshold mean.
  # Clipping at 0.01/0.99 caps the resulting prior mean at about +/- 4.6 log-odds.
  cum_prob_clipped <- pmin(pmax(cum_prob, 0.01), 0.99)
  threshold_means  <- qlogis(cum_prob_clipped)
  # Per-threshold SDs widen as the threshold moves away from zero, because extreme
  # thresholds sit where the pilot has fewest observations and so are the least
  # precisely located. The particular constants (a 1.5 baseline shrinking by 0.10
  # per log-odds, floored at 0.8) are a pragmatic choice rather than a derived
  # result: they keep every threshold prior clearly weaker than the pilot
  # likelihood while remaining proper, which is what Bayes-factor validity
  # requires. Given the 0.01/0.99 clipping above, a threshold mean cannot exceed
  # about 4.6 in absolute value, so the realised SDs span roughly 1.04 to 1.5 and
  # the 0.8 floor never actually binds; it is a guard against a future change to
  # the clipping bounds rather than an active constraint.
  threshold_sds <- pmax(0.8, 1.5 - abs(threshold_means) * 0.10)

  message("[thresholds] Language: ", language,
          " | Cum probs: ", paste(round(cum_prob_clipped, 3), collapse = ", "))
  message("[thresholds] Means (logit): ",
          paste(round(threshold_means, 2), collapse = ", "))

  list(
    threshold_means = threshold_means,
    threshold_sds   = threshold_sds
  )
}

#' Build a brms prior object for a given regime and threshold mode.
#'
#' Coefficient names assume treatment coding with Passive as the reference
#' level and a Semantics_scaled predictor. If a language lacks Pseudo_Passive,
#' the pseudo_passive interaction prior is omitted.
#'
#' @param regime_name Character; one of names(PRIOR_REGIMES).
#' @param threshold_mode Character; one of THRESHOLD_MODES.
#' @param threshold_params Named list from compute_ceiling_calibrated_thresholds()
#'   or NULL for broad mode.
#' @param has_pseudo_passive Logical; whether the language has a
#'   Pseudo_Passive level.
#' @return A brmsprior object suitable for passing to brms::brm().
build_brms_prior <- function(regime_name        = "primary",
                             threshold_mode     = "broad",
                             threshold_params   = NULL,
                             has_pseudo_passive = TRUE) {
  stopifnot(regime_name %in% names(PRIOR_REGIMES))
  stopifnot(threshold_mode %in% THRESHOLD_MODES)
  r <- PRIOR_REGIMES[[regime_name]]

  prior_list <- c(
    brms::prior_string(r$b_default,   class = "b"),
    brms::prior_string(r$b_semantics, class = "b", coef = "Semantics_scaled"),
    brms::prior_string(r$b_active_int, class = "b",
                       coef = "S_TypeActive:Semantics_scaled"),
    brms::prior_string(r$sd,           class = "sd"),
    brms::prior_string(r$cor,          class = "cor")
  )

  if (isTRUE(has_pseudo_passive)) {
    prior_list <- c(
      prior_list,
      brms::prior_string(r$b_pseudo_int, class = "b",
                         coef = "S_TypePseudo_Passive:Semantics_scaled")
    )
  }

  if (threshold_mode == "ceiling_calibrated" && !is.null(threshold_params)) {
    for (k in seq_along(threshold_params$threshold_means)) {
      prior_list <- c(
        prior_list,
        brms::prior_string(
          paste0("normal(",
                 round(threshold_params$threshold_means[k], 3), ", ",
                 round(threshold_params$threshold_sds[k],   3), ")"),
          class = "Intercept",
          coef  = as.character(k)
        )
      )
    }
  } else {
    prior_list <- c(
      prior_list,
      brms::prior_string(r$Intercept, class = "Intercept")
    )
  }

  prior_list
}

#' Drop prior rows that do not correspond to any parameter of the given model.
#' brms errors if a prior references a non-existent parameter — e.g. an `lkj`
#' correlation prior for a model with only uncorrelated (`||`) or
#' intercept-only random effects (no correlation parameters), or a coefficient
#' absent under the data's factor coding. Keeping only matching prior rows lets
#' the same prior object be reused safely across the whole model ladder.
#' @param prior_obj A brmsprior object.
#' @param formula A brms formula / brmsformula (with family).
#' @param data The data the model will be fit to.
#' @return The subset of prior_obj whose (class, coef, group, dpar) exist in the model.
align_prior_to_model <- function(prior_obj, formula, data) {
  # Ask brms which priors this model actually admits. default_prior() is the
  # current name; get_prior() is the deprecated alias, kept as a fallback because
  # ARC's brms may be older than the local one.
  valid <- tryCatch(
    brms::default_prior(formula, data = data),
    error = function(e) tryCatch(brms::get_prior(formula, data = data),
                                 error = function(e2) NULL)
  )
  # If brms could not be asked at all, pass the prior through untouched rather than
  # filtering against an empty table, which would discard every prior and silently
  # fit the model with brms defaults.
  if (is.null(valid) || nrow(valid) == 0) return(prior_obj)

  # In a brms prior table, an empty coef/group/dpar means "every parameter of this
  # class", so a blank field must be treated as a wildcard that matches rather than
  # as a value that must be equal. Without this, the class-wide priors (b, sd, cor)
  # would match nothing and be dropped, leaving the model with default priors.
  blank <- function(x) is.na(x) | x == ""
  keep <- vapply(seq_len(nrow(prior_obj)), function(i) {
    p <- prior_obj[i, , drop = FALSE]
    # Keep the row if at least one real parameter of the model matches it.
    any(valid$class == p$class &
        (blank(p$coef)  | valid$coef  == p$coef) &
        (blank(p$group) | valid$group == p$group) &
        (blank(p$dpar)  | valid$dpar  == p$dpar))
  }, logical(1))

  dropped <- prior_obj[!keep, , drop = FALSE]
  if (nrow(dropped) > 0) {
    message("[prior] Dropped ", nrow(dropped), " prior(s) absent from this model: ",
            paste(unique(paste0(dropped$class,
                                ifelse(blank(dropped$coef), "", paste0(":", dropped$coef)))),
                  collapse = ", "))
  }
  kept <- prior_obj[keep, , drop = FALSE]
  class(kept) <- class(prior_obj)   # preserve brmsprior class after subsetting
  kept
}

#' Every prior regime crossed with every threshold mode, for the sensitivity grid.
#'
#' @return A tibble with one row per (regime_name, threshold_mode) combination:
#'   four regimes by two threshold modes, so eight rows.
#' @details Used by scripts/03_prior_sensitivity.R to enumerate the sensitivity
#'   runs. Deriving the grid from names(PRIOR_REGIMES) and THRESHOLD_MODES rather
#'   than listing it means a regime added above is picked up automatically and
#'   cannot be left out of the sensitivity analysis by oversight.
prior_sensitivity_grid <- function() {
  tidyr::crossing(
    regime_name    = names(PRIOR_REGIMES),
    threshold_mode = THRESHOLD_MODES
  )
}
