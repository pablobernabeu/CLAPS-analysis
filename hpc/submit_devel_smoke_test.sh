#!/bin/bash
# hpc/submit_devel_smoke_test.sh
# Interactive/devel smoke test: compile CmdStan, restore renv, run a minimal brms fit.
# Run this on an ARC interactive node ONLY (never on a login node).
# Usage: bash hpc/submit_devel_smoke_test.sh
#
# Request an interactive node first:
#   srun --partition=devel --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G \
#        --time=01:00:00 --pty bash

set -euo pipefail

# Anchor to the submit directory. Run interactively via `bash`, but stay robust
# if ever launched under sbatch (where BASH_SOURCE points to a spool copy).
cd "${SUBMIT_DIR:-$HOME/design_analysis}"

# ---------------------------------------------------------------------------
# Load modules
# ---------------------------------------------------------------------------
module purge
ARC_R_MODULE="${ARC_R_MODULE:-R}"
module load "$ARC_R_MODULE"

# Load CmdStan module if available; otherwise use renv-installed version
if module load cmdstan 2>/dev/null; then
  echo "[smoke] CmdStan loaded via module."
else
  echo "[smoke] CmdStan module not available; relying on cmdstanr installation."
fi

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export CMDSTANR_OUTPUT_DIR="${TMPDIR:-outputs/cmdstan_tmp}"
mkdir -p "$CMDSTANR_OUTPUT_DIR"
mkdir -p outputs/logs

# ---------------------------------------------------------------------------
# Restore renv
# ---------------------------------------------------------------------------
echo "[smoke] Restoring renv library..."
Rscript -e "renv::restore(prompt = FALSE)"

# ---------------------------------------------------------------------------
# Verify CmdStan
# ---------------------------------------------------------------------------
echo "[smoke] Verifying CmdStan installation..."
Rscript -e "
  library(cmdstanr)
  check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
  cat('[smoke] CmdStan path:', cmdstan_path(), '\n')
  cat('[smoke] CmdStan version:', cmdstan_version(), '\n')
"

# ---------------------------------------------------------------------------
# Minimal brms smoke test (L0 intercepts-only, 2 chains, 500 iter)
# ---------------------------------------------------------------------------
echo "[smoke] Running minimal brms smoke test..."
Rscript - <<'REOF'
suppressPackageStartupMessages({
  library(brms)
})

set.seed(1L)
smoke_data <- data.frame(
  Response    = sample(1:7, 60, replace = TRUE),
  S_Type      = factor(rep(c("Passive","Active"), 30), levels = c("Passive","Active")),
  Semantics_scaled = rnorm(60),
  Participant = factor(paste0("P", rep(1:10, each = 6))),
  Verb        = factor(paste0("V", rep(1:5, times = 12)))
)

fit_smoke <- brm(
  bf(Response ~ S_Type * Semantics_scaled + (1|Participant) + (1|Verb),
     family = cumulative(link = "logit", threshold = "flexible")),
  data    = smoke_data,
  prior   = c(prior(normal(0,1), class = b),
              prior(normal(0,3), class = Intercept),
              prior(exponential(1), class = sd)),
  backend = "cmdstanr",
  iter    = 500, warmup = 250, chains = 2, cores = 2,
  seed    = 1L, silent = 2
)

cat("[smoke] Rhat max:", max(rhat(fit_smoke), na.rm = TRUE), "\n")
cat("[smoke] Smoke test PASSED.\n")
REOF

echo "[smoke] All smoke tests passed."
