# R/11_simulate_pooled_v2.R
# ---------------------------------------------------------------------------
# Pooled cross-language design analysis, grounded in the per-language pilots.
#
# Rationale (see reports/preliminary_sample_size_analysis.qmd, Paths Forward):
# the affectedness effect is verb-level, so per-language power is capped by the
# per-language verb count. Pooling the three languages multiplies the verb-level
# information (about 210 verb-by-language units instead of 72) and tests the
# cross-linguistic claim directly, which is where the theory lives
# (Ambridge, Arnon & Bekman, 2023, doi:10.5070/G6011177).
#
# DGP design choice: each replicate simulates the three languages INDEPENDENTLY
# from their already-validated per-language v2 DGPs (R/10_simulate_from_pilot_v2.R),
# then binds the datasets. No new estimation machinery is introduced on the
# data-generating side; the only new component is the pooled ANALYSIS model,
# which is the pre-existing cross-language ladder rung L5_cross_maximal
# (R/04_model_formulas.R, the OSF reference model). Under assurance, each
# language's focal effects are drawn from its own pilot posterior; draws are
# independent across languages because the pilots were fitted independently.
# The simulated languages are the three CLAPS languages themselves, not draws
# from a language population, so the by-Language random effects are estimated
# by the analysis model but are not part of the DGP.
# ---------------------------------------------------------------------------

# Engine: simulate_from_pilot_v2(), `%||%`, .FOCAL_TERMS (and, transitively,
# the base engine in R/10_simulate_from_pilot.R).
source("R/10_simulate_from_pilot_v2.R")

.POOLED_S_LEVELS <- c("Passive", "Active", "Pseudo_Passive")

#' Simulate one pooled cross-language dataset from the per-language pilot DGPs.
#'
#' @param dgps Named list of v2 DGP specs (from extract_dgp_params_v2), one per
#'   language, e.g. list(English = ..., Turkish = ..., Norwegian = ...).
#' @param n_per_language Integer; participants PER LANGUAGE.
#' @param mode "point" | "assurance" | "safeguard", passed through per language.
#' @param draw_indices Named integer vector (same names as dgps); which pilot
#'   posterior draw each language uses in assurance mode. Supplying these from
#'   the grid (randomly sampled over the full posterior) avoids the
#'   contiguous-first-draws limitation of draw_index = replicate.
#' @param seed Integer; language i uses seed + (i - 1), so grid seeds must be
#'   spaced at least length(dgps) apart (the grid generator uses spacing 10).
#' @return Tibble with Participant, Verb (both language-prefixed to prevent
#'   accidental ID sharing across languages), S_Type (Passive-reference factor
#'   over the union of levels), Semantics_scaled, Response, Language.
simulate_pooled_from_pilots <- function(dgps, n_per_language, mode = "assurance",
                                        draw_indices = NULL, seed = 1L) {
  langs <- names(dgps)
  stopifnot(length(langs) >= 2, !is.null(langs))
  if (is.null(draw_indices)) draw_indices <- setNames(rep(1L, length(langs)), langs)

  sims <- lapply(seq_along(langs), function(i) {
    lang <- langs[[i]]
    sim  <- simulate_from_pilot_v2(dgps[[lang]], n_per_language, mode = mode,
                                   draw_index = as.integer(draw_indices[[lang]]),
                                   seed = as.integer(seed) + (i - 1L))
    # Prefix IDs by language unless already prefixed (Verb_ID in the harmonised
    # pilot, and hence the DGP verb names, already carry the language prefix).
    pfx <- paste0(lang, "_")
    sim$Participant <- ifelse(startsWith(sim$Participant, pfx),
                              sim$Participant, paste0(pfx, sim$Participant))
    sim$Verb        <- ifelse(startsWith(sim$Verb, pfx),
                              sim$Verb, paste0(pfx, sim$Verb))
    sim$S_Type      <- as.character(sim$S_Type)
    sim$Language    <- lang
    sim
  })

  out <- dplyr::bind_rows(sims)
  out$S_Type   <- factor(out$S_Type, levels = intersect(.POOLED_S_LEVELS, unique(out$S_Type)))
  out$Language <- factor(out$Language, levels = langs)
  out
}

#' Run one pooled cross-language cell: simulate -> fit L5_cross_maximal -> BF.
#'
#' Mirrors run_databased_cell_v2() so downstream aggregation works unchanged:
#' the saved rds is list(summary, bf_results, diagnostics) with the same column
#' names (language = "AllLanguages"; n_participants = PER-LANGUAGE N).
#'
#' cell columns: n_participants, mode, draw_index_english, draw_index_turkish,
#'   draw_index_norwegian, prior_source, model_level, prior_regime,
#'   threshold_mode, iter, warmup, chains, seed.
run_pooled_cell_v2 <- function(cell, dgps, out_dir, overwrite = FALSE) {
  source("R/03_define_priors.R");  source("R/04_model_formulas.R")
  source("R/05_hypothesis_tests.R"); source("R/07_extract_diagnostics.R")
  mode <- as.character(cell$mode %||% "assurance")
  psrc <- as.character(cell$prior_source %||% "pilot")
  cell_id <- paste("pooled2", psrc, mode,
                   sprintf("N%03d", as.integer(cell$n_participants)),
                   cell$seed, sep = "_")
  out_file <- file.path(out_dir, paste0(cell_id, ".rds"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(out_file) && !overwrite) { message("[pooled2] skip ", cell_id); return(invisible(readRDS(out_file))) }

  draw_idx <- c(English   = as.integer(cell$draw_index_english   %||% 1L),
                Turkish   = as.integer(cell$draw_index_turkish   %||% 1L),
                Norwegian = as.integer(cell$draw_index_norwegian %||% 1L))
  draw_idx <- draw_idx[names(dgps)]

  t0 <- proc.time()[["elapsed"]]
  sim_data <- simulate_pooled_from_pilots(dgps, as.integer(cell$n_participants),
                                          mode = mode, draw_indices = draw_idx,
                                          seed = as.integer(cell$seed))

  # Pooled data always contains pseudo-passive rows (from English + Turkish),
  # so the pseudo interaction prior and hypothesis apply. NB the pooled H2
  # (pseudo-by-affectedness) coefficient is identified from English and Turkish
  # only, since Norwegian has no pseudo-passive; despite the "AllLanguages"
  # label it is a two-language estimand, mirroring the real planned analysis.
  formula   <- build_multilanguage_ladder()[[cell$model_level]]
  prior_obj <- align_prior_to_model(
    build_brms_prior(regime_name = cell$prior_regime, threshold_mode = cell$threshold_mode,
                     has_pseudo_passive = TRUE), formula, sim_data)
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

  base_cols <- tibble::tibble(cell_id = cell_id, language = "AllLanguages",
                  n_participants = cell$n_participants,
                  n_total = as.integer(cell$n_participants) * length(dgps),
                  mode = mode, prior_source = psrc,
                  draw_index_english = draw_idx[["English"]],
                  draw_index_turkish = draw_idx[["Turkish"]],
                  draw_index_norwegian = draw_idx[["Norwegian"]],
                  n_verbs = sum(vapply(dgps, function(d) length(d$verb_affectedness), integer(1))),
                  prior_regime = cell$prior_regime, seed = cell$seed, runtime_sec = rt)
  if (!is.null(fitres$error)) {
    result <- list(summary = dplyr::bind_cols(tibble::tibble(status = "error",
                   error_message = fitres$error), base_cols))
  } else {
    result <- list(
      summary     = dplyr::bind_cols(tibble::tibble(status = "success"), base_cols,
                       tibble::tibble(model_level = cell$model_level,
                                      threshold_mode = cell$threshold_mode,
                                      iter = samp$iter, warmup = samp$warmup, chains = samp$chains)),
      bf_results  = tryCatch(compute_all_bf(fitres$fit, has_pseudo_passive = TRUE),
                             error = function(e) tibble::tibble(error = conditionMessage(e))),
      # Protected like the fit and BF steps (see run_databased_cell_v2).
      diagnostics = tryCatch(extract_convergence_diagnostics(fitres$fit),
                             error = function(e) tibble::tibble(error = conditionMessage(e))))
  }
  tmp <- paste0(out_file, ".tmp"); saveRDS(result, tmp); file.rename(tmp, out_file)
  message("[pooled2] done ", cell_id, " (", round(rt, 1), "s)")
  result
}
