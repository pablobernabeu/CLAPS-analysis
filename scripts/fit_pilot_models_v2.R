#!/usr/bin/env Rscript
# scripts/fit_pilot_models_v2.R
# ---------------------------------------------------------------------------
# Re-fit ONE language's REAL pilot and save a v2 DGP spec that CARRIES THE
# POSTERIOR DRAWS of the fixed effects + thresholds (needed for the assurance
# and safeguard modes in R/10_simulate_from_pilot_v2.R). The base pilot fit
# (scripts/fit_pilot_models.R) saved only point estimates, so assurance cannot
# be added downstream without this re-fit.
#
# --regime selects the prior the DGP is extracted under:
#   * "primary" (zero-centred; = the local "primary") -> prior_source "pilot":
#       data-dominated DGP, NOT blended with the literature.
#   * "literature_centred" (centred on the pooled anchors 0.47/-0.31/-0.36)
#       -> prior_source "blend": the DGP blends pilot + published evidence
#       (addresses Albers & Lakens "use all available evidence"). Whether this
#       actually MOVES the DGP at the real verb count is an empirical question
#       this script answers by printing the focal-slope comparison.
#
# Usage:
#   Rscript scripts/fit_pilot_models_v2.R --language English --regime primary --prior_source pilot
#   Rscript scripts/fit_pilot_models_v2.R --language English --regime literature_centred --prior_source blend
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr); library(brms) })
source("R/02_preprocess_factors.R"); source("R/03_define_priors.R")
source("R/04_model_formulas.R");      source("R/05_hypothesis_tests.R")
source("R/07_extract_diagnostics.R"); source("R/10_simulate_from_pilot_v2.R")

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--language",     default = "English"),
  optparse::make_option("--regime",       default = "primary"),
  optparse::make_option("--prior_source", default = "pilot"),
  optparse::make_option("--pilot",        default = "data/pilot/claps_pilot_harmonised.csv"),
  optparse::make_option("--outdir",       default = "outputs/pilot_models"),
  optparse::make_option("--iter",         default = 3000L, type = "integer"),
  optparse::make_option("--warmup",       default = 1000L, type = "integer"),
  optparse::make_option("--chains",       default = 4L,    type = "integer")
)))
lang <- opt$language
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

d <- readr::read_csv(opt$pilot, show_col_types = FALSE) |> dplyr::filter(Language == lang)
keep <- c("Passive", "Active", "Pseudo_Passive")
keep <- keep[keep %in% unique(d$S_Type)]
d    <- dplyr::filter(d, S_Type %in% keep)
d    <- scale_semantics(d, centre_by = "Language")
d$S_Type   <- factor(d$S_Type, levels = keep)
d$Verb     <- d$Verb_ID
d$Response <- as.integer(d$Response)
has_pp     <- "Pseudo_Passive" %in% keep

va_df    <- dplyr::distinct(d, Verb, Semantics_scaled)
verb_aff <- setNames(va_df$Semantics_scaled, va_df$Verb)

formula   <- build_model_ladder(has_pseudo_passive = has_pp)[["L5_correlated_maximal"]]
prior_obj <- align_prior_to_model(build_brms_prior(opt$regime, "broad", has_pp), formula, d)
samp      <- production_sampling(iter = opt$iter, warmup = opt$warmup, chains = opt$chains, seed = 2025)

cat("[pilot fit v2]", lang, "| regime", opt$regime, "| source", opt$prior_source,
    "| verbs", length(verb_aff), "| obs", nrow(d), "| ppts", dplyr::n_distinct(d$Participant), "\n")

fit <- brms::brm(formula, data = d, prior = prior_obj, backend = "cmdstanr",
                 sample_prior = "no", iter = samp$iter, warmup = samp$warmup,
                 chains = samp$chains, cores = samp$cores, seed = samp$seed,
                 control = production_control(), silent = 2)

dgp <- extract_dgp_params_v2(fit, verb_affectedness = verb_aff, s_types = keep)
out <- file.path(opt$outdir, paste0("pilot_dgp_v2_", opt$prior_source, "_", lang, ".rds"))
saveRDS(dgp, out)

# Compact summary + focal-slope comparison anchors (for the blend-vs-pilot check).
saveRDS(list(language = lang, regime = opt$regime, prior_source = opt$prior_source,
             fixef = dgp$fixef, focal_lwr = dgp$focal_lwr, ndraws = dgp$ndraws,
             thresholds = dgp$thresholds, n_verbs = length(verb_aff), n_obs = nrow(d),
             diagnostics = tryCatch(extract_convergence_diagnostics(fit), error = function(e) NULL)),
        file.path(opt$outdir, paste0("pilot_dgp_v2_summary_", opt$prior_source, "_", lang, ".rds")))

cat("[pilot fit v2] saved", out, "(", dgp$ndraws, "draws )\n")
cat("focal posterior MEANS:\n"); print(round(dgp$fixef[intersect(.FOCAL_TERMS, names(dgp$fixef))], 3))
cat("focal lower", dgp$lwr_q, "quantile (safeguard):\n"); print(round(dgp$focal_lwr, 3))
