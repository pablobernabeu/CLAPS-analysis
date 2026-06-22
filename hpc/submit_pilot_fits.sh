#!/bin/bash
#SBATCH --job-name=claps_pilotfit
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/pilotfit_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/pilotfit_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_pilot_fits.sh
# Fit the maximal Bayesian model to each language's REAL pilot data (one array
# task per language) and save the data-generating spec used by the data-grounded
# power analysis. Submit:  sbatch --account=PROJECT_GROUP --array=1-3 hpc/submit_pilot_fits.sh

set -euo pipefail
SUBMIT_DIR="$HOME/design_analysis"; cd "$SUBMIT_DIR"
echo "Job $SLURM_JOB_ID | task $SLURM_ARRAY_TASK_ID | acct ${SLURM_JOB_ACCOUNT:-?} | $(hostname) | $(date -Iseconds)"

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
OUTDIR="${PROJECT_DATA}/outputs/pilot_models"
mkdir -p "$OUTDIR" "$CMDSTANR_OUTPUT_DIR" outputs/logs

LANGS=(English Turkish Norwegian)
LANGUAGE="${LANGS[$((SLURM_ARRAY_TASK_ID - 1))]}"
echo "Pilot fit language: $LANGUAGE"

Rscript scripts/fit_pilot_models.R --language "$LANGUAGE" --outdir "$OUTDIR"
echo "End $(date -Iseconds) | exit $?"
