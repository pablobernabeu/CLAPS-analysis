#!/usr/bin/env bash
# Build the clean BASELINE-ONLY design summary (single-language power curve,
# N = 30-100, primary prior, maximal model), excluding the gender variation and
# the cross-language analysis BY SEED RANGE. The primary/maximal cells are
# otherwise indistinguishable between the baseline, gender, and extension studies
# (same language/model_level/prior_regime) -- only the simulation seed separates
# them:
#   baseline single power curve : seeds 100000-102399 (N 30,40,50,60)
#   single prior-sensitivity    : seeds 200000-201799 (N 50, non-primary regimes)
#   gender variation            : seeds 300000-302399 (EXCLUDED)
#   extension                   : seeds 400000-402399 (N 70,80,90,100)
#   cross-language              : seeds 900000+        (EXCLUDED)
# We take seeds 100000-201799 (baseline + prior-sens; prior-sens drops out under
# the primary filter) plus 400000-402399 (extension). The focal table then
# filters prior_regime==primary, model_level==L5_correlated_maximal, n_verbs==20.
set -euo pipefail

module load R/4.4.2-gfbf-2024a 2>/dev/null || true
export R_LIBS_USER=/data/PROJECT_GROUP/PROJECT_ID/PROJECT_GROUP/R/library_4.4

export SRC=/data/PROJECT_GROUP/PROJECT_ID/PROJECT_GROUP/outputs/design_analysis
export TMP=/data/PROJECT_GROUP/PROJECT_ID/PROJECT_GROUP/outputs/_baseline_tmp_rds
export SUM=/data/PROJECT_GROUP/PROJECT_ID/PROJECT_GROUP/outputs/design_summary_baseline

rm -rf "$TMP" "$SUM"
mkdir -p "$TMP" "$SUM"

# Select baseline + extension by seed (awk coerces sci-notation 1e+05 -> 100000).
ls "$SRC"/*.rds | awk -F/ '{
  fn=$NF; n=fn; sub(/\.rds$/,"",n); m=split(n,a,"_"); v=a[m]+0;
  if ((v>=100000 && v<=201799) || (v>=400000 && v<=402399)) print $0
}' | while read -r f; do ln -s "$f" "$TMP"/; done

echo "[select] linked $(ls "$TMP" | wc -l) cells into TMP (expect 6600)"

cd "$HOME/design_analysis"
mkdir -p outputs
Rscript scripts/06_aggregate_design_results.R --out_dir "$TMP" --sum_dir "$SUM" 2>&1 | tail -15

echo "== baseline focal table (primary / maximal / 20 verbs) =="
Rscript -e '
suppressMessages({library(dplyr); library(readr); library(tidyr)})
exc <- read_csv(file.path(Sys.getenv("SUM"), "bf_exceedance.csv"), show_col_types=FALSE)
focal <- c(H1a="H1a_semantics_positive", H1b="H1b_active_interaction_negative")
pc <- exc %>% filter(model_level=="L5_correlated_maximal", prior_regime=="primary",
                     n_verbs==20, hypothesis %in% focal) %>%
  mutate(H = names(focal)[match(hypothesis, focal)])
w <- pc %>% select(language, n_participants, H, p_bf_primary, n_sims) %>%
  tidyr::pivot_wider(names_from = H, values_from = p_bf_primary) %>%
  arrange(language, n_participants)
print(as.data.frame(w), row.names = FALSE)
cat("\n-- smallest N with H1b P(BF>10) >= 0.80 (recommended sample size) --\n")
rec <- pc %>% filter(H == "H1b", p_bf_primary >= 0.80) %>%
  group_by(language) %>% summarise(N_recommended = min(n_participants), .groups = "drop")
print(as.data.frame(rec), row.names = FALSE)
'
