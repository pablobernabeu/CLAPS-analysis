#!/bin/bash
#SBATCH --job-name=claps_databased2
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/databased2_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/databased2_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_databased_v2_array.sh
# Amended data-grounded BAYESIAN power analysis (assurance + safeguard modes).
# Depends on the v2 pilot fits (submit_pilot_fits_v2.sh) having written
# pilot_dgp_v2_<source>_<lang>.rds. Uses a NEW grid, cell runner, simulator, and
# output dir, so the running point-estimate run (design_databased) is untouched.
#   sbatch --account=PROJECT_GROUP --dependency=afterok:<pilotv2_jobid> --array=1-720%20 hpc/submit_databased_v2_array.sh

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
OUTPUT_DIR="${PROJECT_DATA}/outputs/design_databased_v2"
DGP_DIR="${PROJECT_DATA}/outputs/pilot_models"
mkdir -p "$OUTPUT_DIR" "$CMDSTANR_OUTPUT_DIR" outputs/logs

GRID="${GRID:-config/design_grid_databased_v2.csv}"
ROW="$SLURM_ARRAY_TASK_ID"
N_ROWS=$(tail -n +2 "$GRID" | wc -l)
if [[ "$ROW" -gt "$N_ROWS" ]]; then
  echo "[databased2] task $ROW exceeds grid size $N_ROWS; exiting cleanly."; exit 0
fi
echo "Grid row $ROW / $N_ROWS | $(awk -F',' -v r="$((ROW+1))" 'NR==r{print $1,"N="$2,"mode="$3,"src="$5}' "$GRID")"

Rscript scripts/06_databased_cell_v2.R \
  --row_index "$ROW" \
  --grid      "$GRID" \
  --dgpdir    "$DGP_DIR" \
  --outdir    "$OUTPUT_DIR" \
  ${OVERWRITE:+--overwrite}

echo "End $(date -Iseconds) | exit $?"
