#!/bin/bash
#SBATCH --job-name=claps_bf_calibration
#SBATCH --partition=long
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-20%5
#SBATCH --output=/home/%u/design_analysis/outputs/logs/bf_calibration_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/bf_calibration_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_bf_calibration_array.sh
# Bayes-factor calibration: compare Savage-Dickey vs bridge sampling for a
# small subset of design cells (defined in config/bf_calibration_subset.csv).
# CPU-only. Requires prior design-analysis array to have completed.

set -euo pipefail

# Under sbatch, BASH_SOURCE can point to a spool copy in a non-writable path, so
# anchor to the submit directory explicitly (mirrors submit_design_analysis_array.sh).
if [[ -z "${DATA:-}" ]]; then
  echo "ERROR: \$DATA is not set. Cannot write to project storage." >&2
  exit 1
fi
PROJECT_DATA="${DATA}/PROJECT_GROUP"
SUBMIT_DIR="$HOME/design_analysis"
cd "$SUBMIT_DIR"

echo "=========================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Array task:   $SLURM_ARRAY_TASK_ID"
echo "Host:         $(hostname)"
echo "Start time:   $(date -Iseconds)"
echo "Git SHA:      $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "=========================================="

module purge
ARC_R_MODULE="${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"
module load "$ARC_R_MODULE"

export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export CMDSTANR_OUTPUT_DIR="${TMPDIR:-${PROJECT_DATA}/cmdstan_tmp}"
CALIB_OUT="${PROJECT_DATA}/outputs/bf_calibration"
mkdir -p "$CMDSTANR_OUTPUT_DIR" "$CALIB_OUT" outputs/logs

# The calibration subset reuses design_grid.csv rows; indices listed here
# or in a separate calibration_subset.txt file (one row index per line).
SUBSET_FILE="config/bf_calibration_subset.txt"
if [[ ! -f "$SUBSET_FILE" ]]; then
  echo "[calibration] ERROR: Subset file not found: $SUBSET_FILE" >&2
  exit 1
fi

ROW_INDEX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SUBSET_FILE")
if [[ -z "$ROW_INDEX" ]]; then
  echo "[calibration] Array task $SLURM_ARRAY_TASK_ID has no matching row; exiting."
  exit 0
fi

echo "Design grid row: $ROW_INDEX"

Rscript scripts/05_bf_calibration_cell.R \
  --row_index "$ROW_INDEX" \
  --grid      "${GRID:-config/design_grid.csv}" \
  --config    config/analysis_config.yaml \
  --outdir    "$CALIB_OUT" \
  ${OVERWRITE:+--overwrite}

echo "End time: $(date -Iseconds)"
