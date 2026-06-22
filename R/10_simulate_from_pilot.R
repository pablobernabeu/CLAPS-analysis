# R/10_simulate_from_pilot.R
# ---------------------------------------------------------------------------
# Data-grounded ("parametric bootstrap") engine for the BAYESIAN design analysis.
#
# Instead of simulating from ASSUMED effect sizes and an assumed predictor
# distribution, this generates data from parameters CALIBRATED to a model fit to
# the real pilot data:
#   - the real per-verb affectedness ratings (Gelman-scaled, as in the analysis),
#   - the pilot-estimated fixed effects (focal effects optionally discounted to
#     probe a conservative / safeguard scenario),
#   - the pilot-estimated by-participant and by-verb random-effect covariance,
#   - the pilot-estimated ordinal thresholds.
# The refit and the Savage-Dickey Bayes factor downstream are unchanged, so the
# whole design analysis stays Bayesian end-to-end (no frequentist tooling).
#
# extract_dgp_params()   : pull a DGP spec out of a fitted brms pilot model.
# simulate_from_pilot()  : simulate one ordinal dataset from a DGP spec.
# run_databased_cell()   : simulate -> refit (brms) -> Bayes factors -> save.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(brms); library(MASS); library(dplyr); library(tibble); library(purrr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Design value of a named population/RE term (treatment coding, Passive reference)
# for a given sentence type and (scaled) affectedness. Vectorised over a data frame.
.term_columns <- function(s_type, sem) {
  Active <- as.numeric(s_type == "Active")
  Pseudo <- as.numeric(s_type == "Pseudo_Passive")
  list(
    "Intercept"                               = rep(1, length(sem)),
    "S_TypeActive"                            = Active,
    "S_TypePseudo_Passive"                    = Pseudo,
    "Semantics_scaled"                        = sem,
    "S_TypeActive:Semantics_scaled"           = Active * sem,
    "S_TypePseudo_Passive:Semantics_scaled"   = Pseudo * sem
  )
}

#' Extract a data-generating spec from a fitted brms cumulative pilot model.
#' @param fit brmsfit of Response ~ S_Type * Semantics_scaled + (..|Participant) + (..|Verb).
#' @param verb_affectedness named numeric: Verb -> Gelman-scaled affectedness (one per verb).
#' @param s_types ordered character of sentence-type levels (Passive first).
#' @return list DGP spec consumed by simulate_from_pilot().
extract_dgp_params <- function(fit, verb_affectedness, s_types, n_cats = 7L) {
  fe        <- brms::fixef(fit)[, "Estimate"]
  is_thr    <- grepl("^Intercept\\[", names(fe))
  thresholds <- unname(fe[is_thr])
  fixef_pop  <- fe[!is_thr]                      # S_TypeActive, Semantics_scaled, interactions, ...

  vc <- brms::VarCorr(fit)
  build_Sigma <- function(grp) {
    sd <- vc[[grp]]$sd[, "Estimate"]
    terms <- names(sd)
    R <- tryCatch(vc[[grp]]$cor[, "Estimate", ], error = function(e) NULL)
    if (is.null(R) || is.null(dim(R))) R <- diag(length(sd))
    dimnames(R) <- list(terms, terms)
    S <- diag(sd, length(sd)) %*% R %*% diag(sd, length(sd))
    dimnames(S) <- list(terms, terms)
    S
  }
  list(
    fixef             = fixef_pop,
    thresholds        = thresholds,
    Sigma_part        = build_Sigma("Participant"),
    Sigma_verb        = build_Sigma("Verb"),
    verb_affectedness = verb_affectedness,
    s_types           = s_types,
    n_cats            = length(thresholds) + 1L   # categories implied by the fitted cutpoints
  )
}

#' Simulate one CLAPS dataset from a pilot-calibrated DGP spec.
#' @param dgp spec from extract_dgp_params().
#' @param n_participants integer.
#' @param effect_mult multiplier on the FOCAL fixed effects (Semantics slope and its
#'   S_Type interactions) only; 1 = pilot estimate, <1 = conservative/safeguard.
#' @param seed integer.
#' @return tibble (Participant, Verb, S_Type, Semantics_scaled, Response).
simulate_from_pilot <- function(dgp, n_participants, effect_mult = 1.0, seed = 1L) {
  set.seed(seed)
  verbs   <- names(dgp$verb_affectedness)
  s_types <- dgp$s_types
  n_cats  <- dgp$n_cats %||% 7L

  grid <- expand.grid(pi = seq_len(n_participants), Verb = verbs, S_Type = s_types,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  sem  <- unname(dgp$verb_affectedness[grid$Verb])
  tc   <- .term_columns(grid$S_Type, sem)

  # Fixed effects, focal terms discounted by effect_mult.
  ff <- dgp$fixef
  focal <- c("Semantics_scaled", "S_TypeActive:Semantics_scaled",
             "S_TypePseudo_Passive:Semantics_scaled")
  ff[names(ff) %in% focal] <- ff[names(ff) %in% focal] * effect_mult
  eta <- Reduce(`+`, lapply(names(ff), function(t) ff[[t]] * tc[[t]]))

  # By-participant random effects.
  pterms <- colnames(dgp$Sigma_part)
  bP <- MASS::mvrnorm(n_participants, rep(0, length(pterms)), dgp$Sigma_part)
  bP <- matrix(bP, nrow = n_participants, dimnames = list(NULL, pterms))
  eta <- eta + Reduce(`+`, lapply(pterms, function(t) bP[grid$pi, t] * tc[[t]]))

  # By-verb random effects.
  vterms <- colnames(dgp$Sigma_verb)
  bV <- MASS::mvrnorm(length(verbs), rep(0, length(vterms)), dgp$Sigma_verb)
  bV <- matrix(bV, nrow = length(verbs), dimnames = list(verbs, vterms))
  vidx <- match(grid$Verb, verbs)
  eta <- eta + Reduce(`+`, lapply(vterms, function(t) bV[vidx, t] * tc[[t]]))

  # Cumulative-logit -> sample an ordinal response per row.
  cum   <- vapply(dgp$thresholds, function(tk) plogis(tk - eta), numeric(length(eta)))
  cum   <- cbind(0, cum, 1)
  probs <- t(apply(cum, 1, diff))
  resp  <- apply(probs, 1, function(p) { p <- pmax(p, 0); sample.int(n_cats, 1, prob = p / sum(p)) })

  tibble::tibble(
    Participant      = paste0("P", grid$pi),
    Verb             = grid$Verb,
    S_Type           = factor(grid$S_Type, levels = s_types),
    Semantics_scaled = sem,
    Response         = as.integer(resp)
  )
}

#' Run one data-grounded design cell: simulate from pilot DGP, refit, compute BF.
#' Mirrors run_design_cell() but with the pilot-calibrated simulator.
run_databased_cell <- function(cell, dgp, out_dir, overwrite = FALSE) {
  source("R/03_define_priors.R");  source("R/04_model_formulas.R");  source("R/05_hypothesis_tests.R")
  source("R/07_extract_diagnostics.R")
  cell_id <- paste("databased", cell$language, sprintf("N%03d", as.integer(cell$n_participants)),
                   sprintf("mult%03d", as.integer(round(100 * cell$effect_mult))), cell$seed, sep = "_")
  out_file <- file.path(out_dir, paste0(cell_id, ".rds"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(out_file) && !overwrite) { message("[databased] skip ", cell_id); return(invisible(readRDS(out_file))) }

  has_pp <- isTRUE(cell$has_pseudo_passive)
  t0 <- proc.time()[["elapsed"]]
  sim_data <- simulate_from_pilot(dgp, as.integer(cell$n_participants),
                                  effect_mult = as.numeric(cell$effect_mult), seed = as.integer(cell$seed))
  formula   <- build_model_ladder(has_pseudo_passive = has_pp)[[cell$model_level]]
  prior_obj <- align_prior_to_model(
    build_brms_prior(regime_name = cell$prior_regime, threshold_mode = cell$threshold_mode,
                     has_pseudo_passive = has_pp), formula, sim_data)
  samp <- production_sampling(iter = as.integer(cell$iter %||% 3000),
                              warmup = as.integer(cell$warmup %||% 1000),
                              chains = as.integer(cell$chains %||% 4), seed = as.integer(cell$seed))
  ctrl <- production_control()

  fitres <- tryCatch({
    fit <- brms::brm(formula, data = sim_data, prior = prior_obj, backend = "cmdstanr",
                     sample_prior = "yes", iter = samp$iter, warmup = samp$warmup,
                     chains = samp$chains, cores = samp$cores, seed = samp$seed,
                     control = ctrl, silent = 2)
    list(fit = fit, error = NULL)
  }, error = function(e) list(fit = NULL, error = conditionMessage(e)))
  rt <- proc.time()[["elapsed"]] - t0

  if (!is.null(fitres$error)) {
    result <- list(summary = tibble::tibble(cell_id = cell_id, status = "error",
                   error_message = fitres$error, runtime_sec = rt,
                   language = cell$language, n_participants = cell$n_participants,
                   n_verbs = length(dgp$verb_affectedness), effect_mult = cell$effect_mult,
                   prior_regime = cell$prior_regime, seed = cell$seed))
  } else {
    result <- list(
      summary = tibble::tibble(cell_id = cell_id, status = "success", runtime_sec = rt,
                   language = cell$language, n_participants = cell$n_participants,
                   n_verbs = length(dgp$verb_affectedness), effect_mult = cell$effect_mult,
                   model_level = cell$model_level, prior_regime = cell$prior_regime,
                   threshold_mode = cell$threshold_mode, seed = cell$seed,
                   iter = samp$iter, warmup = samp$warmup, chains = samp$chains),
      bf_results  = tryCatch(compute_all_bf(fitres$fit, has_pseudo_passive = has_pp),
                             error = function(e) tibble::tibble(error = conditionMessage(e))),
      diagnostics = extract_convergence_diagnostics(fitres$fit))
  }
  tmp <- paste0(out_file, ".tmp"); saveRDS(result, tmp); file.rename(tmp, out_file)
  message("[databased] done ", cell_id, " (", round(rt, 1), "s)")
  result
}
