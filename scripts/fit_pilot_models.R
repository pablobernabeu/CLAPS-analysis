#!/usr/bin/env Rscript
# scripts/fit_pilot_models.R
# ---------------------------------------------------------------------------
# Fit the maximal Bayesian cumulative model to ONE language's REAL pilot data
# and save a data-generating spec (DGP) for the data-grounded design analysis.
# This is the "fit to the real data" step the PI intends; everything downstream
# (simulate from this fit -> refit -> Bayes factor) stays Bayesian.
#
# Usage:  Rscript scripts/fit_pilot_models.R --language English [--outdir ...]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(brms)
})
source("R/02_preprocess_factors.R"); source("R/03_define_priors.R")
source("R/04_model_formulas.R");      source("R/05_hypothesis_tests.R")
source("R/07_extract_diagnostics.R"); source("R/10_simulate_from_pilot.R")

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--language", default = "English"),
  optparse::make_option("--pilot",    default = "data/pilot/claps_pilot_harmonised.csv"),
  optparse::make_option("--outdir",   default = "outputs/pilot_models"),
  optparse::make_option("--iter",     default = 2000L, type = "integer"),
  optparse::make_option("--warmup",   default = 1000L, type = "integer"),
  optparse::make_option("--chains",   default = 4L,    type = "integer")
)))
lang <- opt$language
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

d <- readr::read_csv(opt$pilot, show_col_types = FALSE) |> dplyr::filter(Language == lang)
present <- unique(d$S_Type)
keep    <- c("Passive", "Active", "Pseudo_Passive")
keep    <- keep[keep %in% present]                       # drops Pseudo for Norwegian, excludes Synthetic_Passive
d       <- dplyr::filter(d, S_Type %in% keep)
d       <- scale_semantics(d, centre_by = "Language")    # Gelman scaling: (x - mean)/(2 sd)
d$S_Type   <- factor(d$S_Type, levels = keep)            # Passive reference, treatment coding
d$Verb     <- d$Verb_ID                                  # grouping factor name expected by the ladder
d$Response <- as.integer(d$Response)
has_pp     <- "Pseudo_Passive" %in% keep

va_df    <- dplyr::distinct(d, Verb, Semantics_scaled)
verb_aff <- setNames(va_df$Semantics_scaled, va_df$Verb)

formula   <- build_model_ladder(has_pseudo_passive = has_pp)[["L5_correlated_maximal"]]
prior_obj <- align_prior_to_model(build_brms_prior("primary", "broad", has_pp), formula, d)
samp      <- production_sampling(iter = opt$iter, warmup = opt$warmup, chains = opt$chains, seed = 2025)

cat("[pilot fit]", lang, "| S_Types:", paste(keep, collapse = ","),
    "| verbs", length(verb_aff), "| obs", nrow(d),
    "| ppts", dplyr::n_distinct(d$Participant), "\n")

fit <- brms::brm(formula, data = d, prior = prior_obj, backend = "cmdstanr",
                 sample_prior = "no", iter = samp$iter, warmup = samp$warmup,
                 chains = samp$chains, cores = samp$cores, seed = samp$seed,
                 control = production_control(), silent = 2)

dgp <- extract_dgp_params(fit, verb_affectedness = verb_aff, s_types = keep)
saveRDS(dgp, file.path(opt$outdir, paste0("pilot_dgp_", lang, ".rds")))
saveRDS(list(language = lang, fixef = dgp$fixef, thresholds = dgp$thresholds,
             sd_part = sqrt(diag(dgp$Sigma_part)), sd_verb = sqrt(diag(dgp$Sigma_verb)),
             n_verbs = length(verb_aff), n_obs = nrow(d),
             n_ppts = dplyr::n_distinct(d$Participant),
             diagnostics = tryCatch(extract_convergence_diagnostics(fit), error = function(e) NULL)),
        file.path(opt$outdir, paste0("pilot_dgp_summary_", lang, ".rds")))

cat("[pilot fit] saved DGP for", lang, "\n")
print(round(dgp$fixef, 3))
cat("thresholds:", paste(round(dgp$thresholds, 2), collapse = ", "), "\n")
cat("participant SDs:", paste(round(sqrt(diag(dgp$Sigma_part)), 2), collapse = ", "), "\n")
cat("verb SDs:",        paste(round(sqrt(diag(dgp$Sigma_verb)), 2), collapse = ", "), "\n")
