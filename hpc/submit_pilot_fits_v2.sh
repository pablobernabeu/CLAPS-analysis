#!/bin/bash
#SBATCH --job-name=claps_pilotv2
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/pilotv2_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/pilotv2_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_pilot_fits_v2.sh
# Re-fit the pilots saving POSTERIOR DRAWS (for the assurance/safeguard engine),
# under the zero-centred ("primary" -> source "pilot") AND literature_centred
# ("blend") priors. 6 array tasks = 3 languages x 2 prior sources.
# Writes pilot_dgp_v2_{pilot,blend}_<lang>.rds alongside the existing point DGPs
# (distinct names; nothing existing is overwritten).
#   sbatch --account=PROJECT_GROUP --array=1-6 hpc/submit_pilot_fits_v2.sh

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
DGP_DIR="${PROJECT_DATA}/outputs/pilot_models"
mkdir -p "$DGP_DIR" "$CMDSTANR_OUTPUT_DIR" outputs/logs

case "$SLURM_ARRAY_TASK_ID" in
  1) LANG=English;   REGIME=primary;           SRC=pilot ;;
  2) LANG=Turkish;   REGIME=primary;           SRC=pilot ;;
  3) LANG=Norwegian; REGIME=primary;           SRC=pilot ;;
  4) LANG=English;   REGIME=literature_centred; SRC=blend ;;
  5) LANG=Turkish;   REGIME=literature_centred; SRC=blend ;;
  6) LANG=Norwegian; REGIME=literature_centred; SRC=blend ;;
  *) echo "[pilotv2] task $SLURM_ARRAY_TASK_ID has no mapping; exiting cleanly."; exit 0 ;;
esac
echo "[pilotv2] LANG=$LANG REGIME=$REGIME SRC=$SRC"

Rscript scripts/fit_pilot_models_v2.R \
  --language     "$LANG" \
  --regime       "$REGIME" \
  --prior_source "$SRC" \
  --outdir       "$DGP_DIR"

echo "End $(date -Iseconds) | exit $?"
