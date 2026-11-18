#!/usr/bin/env Rscript
# scripts/summarise_pilot_fits.R
# ---------------------------------------------------------------------------
# The pilot models in full, for the Supplementary Information of the Stage 1 report:
# every fixed effect and rating threshold of each pilot language's model (posterior mean,
# SD, 95% quantile interval, R-hat and bulk and tail effective sample sizes, computed from
# the saved draws of the four chains), and the convergence summary over all parameters
# that the fitting run recorded.
#
# Inputs (gitignored, as the draws are not public):
#   outputs/pilot_models/pilot_dgp_v2_pilot_<language>.rds          draws, from
#     scripts/fit_pilot_models_v2.R (hpc/submit_pilot_fits_v2.sh, tasks 1-3)
#   outputs/pilot_models/pilot_dgp_v2_summary_pilot_<language>.rds  the run's summary,
#     whose `diagnostics` covers every parameter, random-effect SDs and correlations
#     included (R/07_extract_diagnostics.R)
#
# Writes:
#   outputs/design_summary_pilot/pilot_fit_parameters.csv   one row per language and parameter
#   outputs/design_summary_pilot/pilot_fit_diagnostics.csv  one row per language
#
# Usage (from design_analysis/): Rscript scripts/summarise_pilot_fits.R
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(posterior))

pm <- "outputs/pilot_models"
langs <- c("English", "Norwegian", "Turkish")
par_rows <- list()
diag_rows <- list()
for (L in langs) {
  p <- readRDS(file.path(pm, paste0("pilot_dgp_v2_pilot_", L, ".rds")))
  s <- readRDS(file.path(pm, paste0("pilot_dgp_v2_summary_pilot_", L, ".rds")))
  stopifnot(identical(s$language, L), identical(s$prior_source, "pilot"))
  for (nm in c("fixef_draws", "thresholds_draws")) {
    x <- p[[nm]]
    # The R-hat and ESS below are only meaningful with the chains kept apart.
    stopifnot(attr(x, "nchains") == s$diagnostics$n_chains)
    sm <- summarise_draws(x, mean, sd, ~quantile(.x, probs = c(0.025, 0.975)), rhat,
                          ess_bulk, ess_tail)
    names(sm) <- c("parameter", "mean", "sd", "lo95", "hi95", "rhat", "ess_bulk", "ess_tail")
    sm$language <- L
    sm$block <- if (nm == "fixef_draws") "fixed effect" else "threshold"
    sm$n_draws <- ndraws(x)
    par_rows[[length(par_rows) + 1L]] <- as.data.frame(sm)
  }
  dg <- as.data.frame(s$diagnostics)
  dg$language <- L
  dg$n_verbs <- s$n_verbs
  dg$n_obs <- s$n_obs
  diag_rows[[length(diag_rows) + 1L]] <- dg
}
pars <- do.call(rbind, par_rows)
pars <- pars[, c("language", "block", "parameter", "mean", "sd", "lo95", "hi95", "rhat",
                 "ess_bulk", "ess_tail", "n_draws")]
diags <- do.call(rbind, diag_rows)
diags <- diags[, c("language", "n_obs", "n_verbs", "n_params", "n_chains", "n_iter", "n_warmup",
                   "max_rhat", "min_ess_bulk", "min_ess_tail", "n_divergent", "n_max_treedepth",
                   "convergence_ok")]
write.csv(pars, "outputs/design_summary_pilot/pilot_fit_parameters.csv", row.names = FALSE)
write.csv(diags, "outputs/design_summary_pilot/pilot_fit_diagnostics.csv", row.names = FALSE)
print(diags, digits = 5)
