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
  primary = list(
    b_default     = "normal(0, 1.5)",
    b_semantics   = "normal(0, 0.5)",
    b_active_int  = "normal(0, 0.5)",
    b_pseudo_int  = "normal(0, 0.6)",
    Intercept     = "student_t(3, 0, 2.5)",
    sd            = "student_t(3, 0, 1)",
    cor           = "lkj(2)"
  ),

  # ---- Weak-prior sensitivity ------------------------------------------
  # Replaces focal slope SDs with normal(0, 1) and non-focal SDs with
  # normal(0, 2). Conclusions should not depend on the more concentrated
  # primary prior (Schad et al., 2021, 2023).
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
  counts <- tabulate(pilot_df$Response, nbins = n_cats)
  smoothed <- (counts + smooth_alpha) / (sum(counts) + n_cats * smooth_alpha)
  cum_prob <- cumsum(smoothed)[seq_len(n_cats - 1)]
  cum_prob_clipped <- pmin(pmax(cum_prob, 0.01), 0.99)
  threshold_means  <- qlogis(cum_prob_clipped)
  # Per-threshold SDs are wider when the threshold is far from zero, to
  # reflect that extreme thresholds are estimated from fewer observations.
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
  valid <- tryCatch(
    brms::default_prior(formula, data = data),
    error = function(e) tryCatch(brms::get_prior(formula, data = data),
                                 error = function(e2) NULL)
  )
  if (is.null(valid) || nrow(valid) == 0) return(prior_obj)

  blank <- function(x) is.na(x) | x == ""
  keep <- vapply(seq_len(nrow(prior_obj)), function(i) {
    p <- prior_obj[i, , drop = FALSE]
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

#' Return a data frame of all prior x threshold combinations for the
#' sensitivity grid.
prior_sensitivity_grid <- function() {
  tidyr::crossing(
    regime_name    = names(PRIOR_REGIMES),
    threshold_mode = THRESHOLD_MODES
  )
}
