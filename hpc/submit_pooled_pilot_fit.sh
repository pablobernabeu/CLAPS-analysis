#!/bin/bash
#SBATCH --job-name=claps_pooledfit
# long, not medium: the first fit of this model to the real pilot data took
# 7 days 23 hours, far beyond medium's 48-hour cap and beyond the 42 hours this
# script used to request. The declared limit now reflects the measured runtime.
#SBATCH --partition=long
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# 4 CPUs, not the 8 this script used to ask for: brms runs one process per chain and
# this fit uses 4 chains with no within-chain threading, so the extra cores sat idle
# for all 7d23h of job 8377363, which came in at 45% CPU efficiency. 8G, not 32G:
# that same run peaked at 2.23 GiB. Together these free ~838 core-hours per run and
# leave a job the scheduler can backfill.
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/pooledfit_%j.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/pooledfit_%j.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_pooled_pilot_fit.sh
# Fit the pooled cross-language model (L5_cross_maximal) to the combined REAL
# pilot data and save a light extract (scripts/fit_pilot_models_pooled.R).
# Independent of the pooled design-analysis array; used for the analytic
# ceiling of the pooled test and as a convergence check on real data.
#
#   sbatch --clusters=htc --account=PROJECT_GROUP hpc/submit_pooled_pilot_fit.sh

set -euo pipefail
SUBMIT_DIR="$HOME/design_analysis"; cd "$SUBMIT_DIR"
echo "Job $SLURM_JOB_ID | acct ${SLURM_JOB_ACCOUNT:-?} | $(hostname) | $(date -Iseconds)"

module purge
module load "${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"

if [[ -z "${DATA:-}" ]]; then echo "ERROR: \$DATA is not set." >&2; exit 1; fi
PROJECT_DATA="${DATA}/PROJECT_GROUP"
export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export CMDSTANR_OUTPUT_DIR="${TMPDIR:-${PROJECT_DATA}/cmdstan_tmp}"
export CMDSTAN="$(ls -d "${PROJECT_DATA}/cmdstan/cmdstan-"* 2>/dev/null | sort -V | tail -1)"
if [[ -z "${CMDSTAN}" || ! -x "${CMDSTAN}/bin/stanc" ]]; then
  echo "ERROR: CmdStan not found under ${PROJECT_DATA}/cmdstan." >&2; exit 1
fi
mkdir -p "$CMDSTANR_OUTPUT_DIR" outputs/logs

# MODEL_LEVEL and OUTFILE may be set on the sbatch line, exported with
# --export=ALL,MODEL_LEVEL,OUTFILE. Naming them is required: SLURM does not pass
# the submitting shell's variables unless it is told to.
MODEL_LEVEL="${MODEL_LEVEL:-L5_cross_maximal}"
OUTFILE="${OUTFILE:-pilot_fit_pooled_extract.rds}"
echo "Model level: ${MODEL_LEVEL} | output: ${OUTFILE}"
Rscript scripts/fit_pilot_models_pooled.R --regime primary \
  --model_level "${MODEL_LEVEL}" \
  --outfile "${OUTFILE}" \
  --outdir "${PROJECT_DATA}/outputs/pilot_models"

echo "End $(date -Iseconds) | exit $?"
