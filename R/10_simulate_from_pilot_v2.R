# R/10_simulate_from_pilot_v2.R
# ---------------------------------------------------------------------------
# Amended data-grounded engine addressing the Albers & Lakens (2018) pilot-power
# critique of the point-estimate plug-in in R/10_simulate_from_pilot.R.
#
# The base engine plugs the pilot posterior MEANS of the focal slopes in as the
# true effect (a conditional-power / fixed-point design analysis). This variant
# adds two data-generating MODES on top of the same machinery:
#
#   * "point"      : base behaviour (focal = posterior mean, optional effect_mult).
#                    Kept as the labelled reference.
#   * "assurance"  : on each replicate, draw the WHOLE fixed-effect + threshold
#                    vector from the pilot POSTERIOR (one posterior draw per
#                    replicate). Averaging P(BF>=10) over replicates then
#                    integrates the pilot's estimation uncertainty in the focal
#                    effects -> assurance, not point power (the A&L-recommended
#                    move).
#   * "safeguard"  : focal slopes set to a LOWER posterior credible bound
#                    (default 10th percentile), a principled conservative effect
#                    in place of the ad hoc effect_mult = 0.6.
#
# The DGP spec additionally carries the posterior draws of the population
# coefficients and thresholds (extract_dgp_params_v2). Random-effect covariance
# is still taken at the posterior mean (its uncertainty is second-order for the
# focal-effect question; documented as a residual limitation).
#
# The refit + Savage-Dickey Bayes factor downstream are unchanged (zero-centred
# analysis prior), so the test itself is never pre-loaded.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(brms); library(posterior); library(MASS); library(dplyr); library(tibble); library(purrr)
})

# Base engine: .term_columns(), extract_dgp_params(), simulate_from_pilot(),
# run_databased_cell(), and `%||%`.
source("R/10_simulate_from_pilot.R")

.FOCAL_TERMS <- c("Semantics_scaled",
                  "S_TypeActive:Semantics_scaled",
                  "S_TypePseudo_Passive:Semantics_scaled")

#' Extract a DGP spec that ALSO carries posterior draws (for assurance) and a
#' lower credible bound on the focal slopes (for safeguard).
#' @inheritParams extract_dgp_params
#' @param lwr_q numeric in (0,1); lower-tail quantile for the safeguard effect.
extract_dgp_params_v2 <- function(fit, verb_affectedness, s_types, n_cats = 7L, lwr_q = 0.10) {
  dgp <- extract_dgp_params(fit, verb_affectedness, s_types, n_cats)

  dr  <- posterior::as_draws_matrix(fit)
  cn  <- colnames(dr)
  bc  <- cn[grepl("^b_", cn)]
  B   <- dr[, bc, drop = FALSE]
  colnames(B) <- sub("^b_", "", bc)

  is_thr   <- grepl("^Intercept\\[", colnames(B))
  thr_draws <- B[, is_thr, drop = FALSE]
  pop_draws <- B[, !is_thr, drop = FALSE]

  # Align population-draw columns to the point-estimate fixef ordering, and
  # thresholds to Intercept[1..K] order, so downstream indexing is unambiguous.
  pop_draws <- pop_draws[, names(dgp$fixef), drop = FALSE]
  thr_idx   <- as.integer(gsub(".*\\[(\\d+)\\].*", "\\1", colnames(thr_draws)))
  thr_draws <- thr_draws[, order(thr_idx), drop = FALSE]

  focal     <- intersect(.FOCAL_TERMS, colnames(pop_draws))
  focal_lwr <- vapply(focal, function(f) as.numeric(stats::quantile(pop_draws[, f], probs = lwr_q)), numeric(1))
  names(focal_lwr) <- focal

  c(dgp, list(
    fixef_draws      = pop_draws,
    thresholds_draws = thr_draws,
    focal_lwr        = focal_lwr,
    lwr_q            = lwr_q,
    ndraws           = nrow(pop_draws)
  ))
}

#' Simulate one CLAPS dataset from a pilot DGP under a given MODE.
#' @param dgp spec from extract_dgp_params_v2() (point mode also works with
#'   extract_dgp_params()).
#' @param mode "point" | "assurance" | "safeguard".
#' @param draw_index integer; which posterior draw to use in "assurance" mode
#'   (wrapped modulo ndraws). Ignored otherwise.
#' @param effect_mult multiplier on the focal terms in "point" mode only.
simulate_from_pilot_v2 <- function(dgp, n_participants, mode = "assurance",
                                   draw_index = 1L, effect_mult = 1.0, seed = 1L) {
  set.seed(seed)
  verbs   <- names(dgp$verb_affectedness)
  s_types <- dgp$s_types
  n_cats  <- dgp$n_cats %||% 7L

  # --- choose the data-generating fixed effects + thresholds for this mode ---
  if (identical(mode, "point")) {
    ff  <- dgp$fixef
    thr <- dgp$thresholds
    ff[names(ff) %in% .FOCAL_TERMS] <- ff[names(ff) %in% .FOCAL_TERMS] * effect_mult
  } else if (identical(mode, "assurance")) {
    if (is.null(dgp$fixef_draws)) stop("[v2] assurance mode needs a v2 DGP (posterior draws).")
    di  <- ((as.integer(draw_index) - 1L) %% dgp$ndraws) + 1L
    ff  <- dgp$fixef_draws[di, ]
    thr <- dgp$thresholds_draws[di, ]
    names(ff) <- colnames(dgp$fixef_draws)
  } else if (identical(mode, "safeguard")) {
    if (is.null(dgp$focal_lwr)) stop("[v2] safeguard mode needs a v2 DGP (focal_lwr).")
    ff  <- dgp$fixef
    thr <- dgp$thresholds
    for (f in names(dgp$focal_lwr)) ff[f] <- dgp$focal_lwr[[f]]
  } else stop("[v2] unknown mode: ", mode)

  grid <- expand.grid(pi = seq_len(n_participants), Verb = verbs, S_Type = s_types,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  sem  <- unname(dgp$verb_affectedness[grid$Verb])
  tc   <- .term_columns(grid$S_Type, sem)

  eta <- Reduce(`+`, lapply(names(ff), function(t) if (!is.null(tc[[t]])) ff[[t]] * tc[[t]] else 0))

  pterms <- colnames(dgp$Sigma_part)
  bP <- MASS::mvrnorm(n_participants, rep(0, length(pterms)), dgp$Sigma_part)
  bP <- matrix(bP, nrow = n_participants, dimnames = list(NULL, pterms))
  eta <- eta + Reduce(`+`, lapply(pterms, function(t) bP[grid$pi, t] * tc[[t]]))

  vterms <- colnames(dgp$Sigma_verb)
  bV <- MASS::mvrnorm(length(verbs), rep(0, length(vterms)), dgp$Sigma_verb)
  bV <- matrix(bV, nrow = length(verbs), dimnames = list(verbs, vterms))
  vidx <- match(grid$Verb, verbs)
  eta <- eta + Reduce(`+`, lapply(vterms, function(t) bV[vidx, t] * tc[[t]]))

  cum   <- vapply(thr, function(tk) plogis(tk - eta), numeric(length(eta)))
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

#' Run one v2 data-grounded cell (assurance / safeguard / point), refit, BF, save.
#' cell columns: language, n_participants, mode, draw_index, prior_source,
#'   model_level, prior_regime, threshold_mode, has_pseudo_passive, iter, warmup,
#'   chains, seed, effect_mult (point only).
run_databased_cell_v2 <- function(cell, dgp, out_dir, overwrite = FALSE) {
  source("R/03_define_priors.R");  source("R/04_model_formulas.R");  source("R/05_hypothesis_tests.R")
  source("R/07_extract_diagnostics.R")
  mode <- as.character(cell$mode %||% "assurance")
  psrc <- as.character(cell$prior_source %||% "pilot")
  cell_id <- paste("databased2", psrc, mode, cell$language,
                   sprintf("N%03d", as.integer(cell$n_participants)),
                   cell$seed, sep = "_")
  out_file <- file.path(out_dir, paste0(cell_id, ".rds"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(out_file) && !overwrite) { message("[databased2] skip ", cell_id); return(invisible(readRDS(out_file))) }

  has_pp <- isTRUE(cell$has_pseudo_passive)
  t0 <- proc.time()[["elapsed"]]
  sim_data <- simulate_from_pilot_v2(dgp, as.integer(cell$n_participants), mode = mode,
                                     draw_index = as.integer(cell$draw_index %||% 1L),
                                     effect_mult = as.numeric(cell$effect_mult %||% 1.0),
                                     seed = as.integer(cell$seed))
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

  base_cols <- tibble::tibble(cell_id = cell_id, language = cell$language,
                  n_participants = cell$n_participants, mode = mode, prior_source = psrc,
                  draw_index = as.integer(cell$draw_index %||% 1L),
                  n_verbs = length(dgp$verb_affectedness), prior_regime = cell$prior_regime,
                  seed = cell$seed, runtime_sec = rt)
  if (!is.null(fitres$error)) {
    result <- list(summary = dplyr::bind_cols(tibble::tibble(status = "error",
                   error_message = fitres$error), base_cols))
  } else {
    result <- list(
      summary     = dplyr::bind_cols(tibble::tibble(status = "success"), base_cols,
                       tibble::tibble(model_level = cell$model_level,
                                      threshold_mode = cell$threshold_mode,
                                      iter = samp$iter, warmup = samp$warmup, chains = samp$chains)),
      bf_results  = tryCatch(compute_all_bf(fitres$fit, has_pseudo_passive = has_pp),
                             error = function(e) tibble::tibble(error = conditionMessage(e))),
      diagnostics = extract_convergence_diagnostics(fitres$fit))
  }
  tmp <- paste0(out_file, ".tmp"); saveRDS(result, tmp); file.rename(tmp, out_file)
  message("[databased2] done ", cell_id, " (", round(rt, 1), "s)")
  result
}
