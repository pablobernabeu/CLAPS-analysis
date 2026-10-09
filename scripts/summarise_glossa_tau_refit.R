#!/usr/bin/env Rscript
# scripts/summarise_glossa_tau_refit.R
# ---------------------------------------------------------------------------
# Summarise the three arms of the Glossa tau refit into one tracked table.
#
# scripts/refit_glossa_pooled_tau.R writes one .rds per arm on the cluster, and
# those files are not tracked. The design note and the K-language driver quote
# the refit's between-language SD of the interaction, so without a tracked
# summary those figures could not be traced to anything in the repository. This
# script reads the three arms and writes the figures they quote, with the
# diagnostics that say how far each arm can be trusted.
#
# The refit's tau is a one-stage SD over the five Glossa languages. It is NOT the
# meta-analytic Glossa tau the K-language engine uses (0.53 over four languages),
# so this table documents the refit and is not read by the engine.
#
# Output (to --outdir):
#   glossa_tau_refit_summary.csv   One row per arm: the specification, the
#     number of observations and runtime; the posterior of tau (the by-Language SD
#     of the active-by-affectedness interaction) and of the interaction itself;
#     the pseudo-passive-by-affectedness interaction and its by-Language SD; the
#     passive slope's by-Language SD; the sampler diagnostics; and, for the arms
#     that split the verb term, the share of each verb term's variance common to all
#     languages (columns added on 9 October 2026, the earlier columns unchanged).
#
# Run from design_analysis/, with the arms copied from the cluster:
#   Rscript scripts/summarise_glossa_tau_refit.R --indir <copy of glossa_tau_refit>
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(optparse))

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  # No default. The arms are not tracked, so a default path would only ever find
  # nothing, and a summary of nothing must not be written.
  optparse::make_option("--indir"),
  optparse::make_option("--outdir", default = "outputs/design_summary_pilot")
)))
tag <- "[summarise glossa tau]"
if (is.null(opt$indir) || !dir.exists(opt$indir)) {
  stop(tag, " --indir must name the directory holding the refit's arm .rds files")
}

# The arms in the order the refit script defines them, with what each changes.
ARMS <- c(
  M0_published  = "published specification, reproduced",
  M1_verbsplit  = "verb term split into shared and language-specific parts",
  M2_thresholds = "as M1, plus per-language response thresholds"
)

# One posterior summary row, or NAs when a term is absent from an arm.
fixed_row <- function(fx, term) {
  if (!term %in% rownames(fx)) return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  unname(as.numeric(fx[term, c("Estimate", "Est.Error", "Q2.5", "Q97.5")]))
}
sd_row <- function(lg, term) {
  key <- paste0("sd(", term, ")")
  if (!key %in% rownames(lg)) return(c(NA_real_, NA_real_, NA_real_))
  unname(as.numeric(lg[key, 1:4][c(1, 3, 4)]))
}

# The verb terms whose variance M1 and M2 split into a part common to all languages
# (group "verb") and a language-specific part (group "VerbLang").
VERB_TERMS <- c(intercept = "Intercept", active = "S_TypeActive",
                pseudo_passive = "S_TypePseudo_Passive")

#' Share of a verb term's variance common to all languages, or NA without the split.
#'
#' With a shared and a language-specific deviation for each verb, the share
#' sd_shared^2 / (sd_shared^2 + sd_specific^2) is the correlation the model implies
#' between two languages' deviations for the same verb. It is computed from the
#' posterior means of the two SDs.
verb_share <- function(sr, term) {
  key <- paste0("sd(", term, ")")
  if (!all(c("verb", "VerbLang") %in% names(sr)) ||
        !key %in% rownames(sr$verb) || !key %in% rownames(sr$VerbLang)) {
    return(NA_real_)
  }
  shared <- sr$verb[key, "Estimate"]^2
  shared / (shared + sr$VerbLang[key, "Estimate"]^2)
}

rows <- lapply(names(ARMS), function(arm) {
  f <- file.path(opt$indir, paste0("glossa_tau_", arm, ".rds"))
  # Every arm must be present: the point of the table is the comparison between
  # them, and a table missing one would read as though that arm were not run.
  if (!file.exists(f)) stop(tag, " missing arm: ", basename(f))
  r <- readRDS(f)
  if (!identical(r$arm, arm)) stop(tag, " ", basename(f), " records arm '", r$arm, "'")
  tau  <- r$tau_language_interaction
  beta <- r$beta_interaction
  fx   <- r$summary_fixed
  lg   <- r$summary_random$Language
  pp   <- fixed_row(fx, "S_TypePseudo_Passive:Semantics")
  sdpp <- sd_row(lg, "S_TypePseudo_Passive:Semantics")
  sdse <- sd_row(lg, "Semantics")
  vs   <- vapply(VERB_TERMS, verb_share, numeric(1), sr = r$summary_random)
  dg   <- r$diagnostics
  data.frame(
    arm = arm, specification = ARMS[[arm]],
    n_obs = r$n_obs, runtime_h = round(r$runtime_h, 1),
    tau_median = tau[["median"]], tau_mean = tau[["mean"]], tau_sd = tau[["sd"]],
    tau_q025 = tau[["q025"]], tau_q975 = tau[["q975"]],
    interaction_median = beta[["median"]], interaction_mean = beta[["mean"]],
    interaction_q025 = beta[["q025"]], interaction_q975 = beta[["q975"]],
    pseudo_passive_estimate = pp[1], pseudo_passive_se = pp[2],
    pseudo_passive_q025 = pp[3], pseudo_passive_q975 = pp[4],
    sd_language_pseudo_passive = sdpp[1], sd_language_pseudo_passive_q025 = sdpp[2],
    sd_language_pseudo_passive_q975 = sdpp[3],
    sd_language_semantics = sdse[1], sd_language_semantics_q025 = sdse[2],
    sd_language_semantics_q975 = sdse[3],
    rhat_max = dg$rhat_max, ess_bulk_min = dg$ess_bulk_min,
    divergent = dg$divergent, treedepth_saturated = dg$treedepth_sat,
    verb_share_intercept = vs[["intercept"]], verb_share_active = vs[["active"]],
    verb_share_pseudo_passive = vs[["pseudo_passive"]],
    stringsAsFactors = FALSE
  )
})
out <- do.call(rbind, rows)

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
path <- file.path(opt$outdir, "glossa_tau_refit_summary.csv")
utils::write.csv(out, path, row.names = FALSE)
message(tag, " wrote ", path)
