#!/usr/bin/env Rscript
# scripts/fit_pilot_models_pooled.R
# ---------------------------------------------------------------------------
# Fit the pooled cross-language analysis model (L5_cross_maximal, the OSF
# reference) to the COMBINED real pilot data of all three languages, and save a
# light extract: population-coefficient posterior draws, random-effect
# summaries, and per-language verb affectedness.
#
# This fit is NOT consumed by the pooled design-analysis array, whose DGP is
# the three per-language v2 pilot DGPs. It serves two purposes: (a) it shows
# that the pooled analysis model converges on real CLAPS data before thousands
# of simulated fits assume it does, and (b) its draws and by-verb variability
# feed the analytic ceiling calculation for the pooled test, mirroring the
# per-language ceilings in the report.
#
# Usage: Rscript scripts/fit_pilot_models_pooled.R [--regime primary]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr); library(brms); library(posterior) })
source("R/02_preprocess_factors.R"); source("R/03_define_priors.R")
source("R/04_model_formulas.R");      source("R/07_extract_diagnostics.R")

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--regime", default = "primary"),
  optparse::make_option("--pilot",  default = "data/pilot/claps_pilot_harmonised.csv"),
  optparse::make_option("--outdir", default = "outputs/pilot_models"),
  optparse::make_option("--iter",   default = 3000L, type = "integer"),
  optparse::make_option("--warmup", default = 1000L, type = "integer"),
  optparse::make_option("--chains", default = 4L,    type = "integer")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

LANGS <- c("English", "Turkish", "Norwegian")
KEEP  <- c("Passive", "Active", "Pseudo_Passive")

d <- readr::read_csv(opt$pilot, show_col_types = FALSE) |>
  dplyr::filter(Language %in% LANGS, S_Type %in% KEEP) |>
  scale_semantics(centre_by = "Language")
# Language-prefixed IDs prevent accidental sharing of a Participant or Verb
# level across languages (verb inventories and participant pools are disjoint).
d$Participant <- paste(d$Language, d$Participant, sep = "_")
d$Verb        <- paste(d$Language, d$Verb_ID,     sep = "_")
d$S_Type      <- factor(d$S_Type, levels = KEEP)
d$Language    <- factor(d$Language, levels = LANGS)
d$Response    <- as.integer(d$Response)

formula   <- build_multilanguage_ladder()[["L5_cross_maximal"]]
prior_obj <- align_prior_to_model(build_brms_prior(opt$regime, "broad", TRUE), formula, d)
samp      <- production_sampling(iter = opt$iter, warmup = opt$warmup, chains = opt$chains, seed = 2026)

cat("[pooled pilot fit] obs", nrow(d), "| ppts", dplyr::n_distinct(d$Participant),
    "| verbs", dplyr::n_distinct(d$Verb), "| languages", nlevels(d$Language), "\n")

fit <- brms::brm(formula, data = d, prior = prior_obj, backend = "cmdstanr",
                 sample_prior = "no", iter = samp$iter, warmup = samp$warmup,
                 chains = samp$chains, cores = samp$cores, seed = samp$seed,
                 control = production_control(), silent = 2)

# Light extract: everything the analytic-ceiling calculation and reporting
# need, without persisting the multi-gigabyte brmsfit.
dr <- posterior::as_draws_matrix(fit)
bc <- colnames(dr)[grepl("^b_", colnames(dr)) & !grepl("^b_Intercept\\[", colnames(dr))]
B  <- dr[, bc, drop = FALSE]; colnames(B) <- sub("^b_", "", colnames(B))

va <- d |> dplyr::distinct(Language, Verb, Semantics_scaled)
vc <- brms::VarCorr(fit)

out <- list(
  fixef_draws       = B,
  fixef             = colMeans(B),
  varcorr           = vc,
  verb_affectedness = setNames(va$Semantics_scaled, va$Verb),
  verb_language     = setNames(as.character(va$Language), va$Verb),
  n_obs             = nrow(d),
  regime            = opt$regime,
  diagnostics       = tryCatch(extract_convergence_diagnostics(fit), error = function(e) NULL)
)
saveRDS(out, file.path(opt$outdir, "pilot_fit_pooled_extract.rds"))

cat("[pooled pilot fit] saved pilot_fit_pooled_extract.rds |",
    nrow(B), "draws x", ncol(B), "coefficients\n")
cat("pooled focal posterior means:\n")
print(round(out$fixef[grep("Semantics_scaled", names(out$fixef))], 3))
