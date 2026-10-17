#!/usr/bin/env Rscript
# scripts/summarise_pilot_interactions.R
# ---------------------------------------------------------------------------
# The two focal interactions in each pilot language, per SD of affectedness: posterior
# mean, 95% quantile interval and the posterior probability of the predicted (negative)
# sign, from the posterior draws of the pilot models that the design analysis uses
# (outputs/pilot_models/pilot_dgp_v2_pilot_<language>.rds, scripts/fit_pilot_models_v2.R).
# Each coefficient is multiplied by the SD of the scaled predictor in that pilot, as the
# design analysis does throughout.
#
# The report's table of pilot estimates and the Stage 1 report's Pilot data section read
# this summary, so that neither needs the pilot draws, which are not public.
#
# Writes outputs/design_summary_pilot/pilot_interactions.csv.
# Usage (from design_analysis/): Rscript scripts/summarise_pilot_interactions.R
# ---------------------------------------------------------------------------
pm <- "outputs/pilot_models"
pp <- read.csv("outputs/design_summary_pilot/pilot_params_ceilings.csv")
terms <- c(H1 = "S_TypeActive:Semantics_scaled", H2 = "S_TypePseudo_Passive:Semantics_scaled")
rows <- list()
for (L in c("English", "Turkish", "Norwegian")) {
  p <- readRDS(file.path(pm, paste0("pilot_dgp_v2_pilot_", L, ".rds")))
  d <- as.matrix(p$fixef_draws)
  sdx <- pp$affectedness_sd[pp$language == L]
  stopifnot(length(sdx) == 1L)
  for (h in names(terms)) {
    if (!terms[[h]] %in% colnames(d)) next   # Norwegian has no pseudo-passive
    x <- d[, terms[[h]]] * sdx
    row <- data.frame(language = L, hypothesis = h, term = terms[[h]],
                      per_sd_mean = mean(x), per_sd_lo95 = unname(stats::quantile(x, 0.025)),
                      per_sd_hi95 = unname(stats::quantile(x, 0.975)), p_negative = mean(x < 0),
                      n_draws = length(x), affectedness_sd = sdx)
    rows[[length(rows) + 1L]] <- row
  }
}
out <- do.call(rbind, rows)
write.csv(out, "outputs/design_summary_pilot/pilot_interactions.csv", row.names = FALSE)
print(out, digits = 3)
