#!/bin/bash
#SBATCH --job-name=claps_corrected
#SBATCH --partition=short
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/corrected_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/corrected_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_corrected_array.sh
# Corrected single-language power BFDA (per-verb affectedness fix + n_verbs sweep).
# Single-language L5 cells are small, so this is deliberately backfill-friendly:
# short partition, 8 CPU, 16 GB, 8 h walltime. The ACCOUNT and ARRAY range are
# passed at submit so the run can be split across PROJECT_GROUP and PROJECT_GROUP, e.g.:
#   sbatch --account=PROJECT_GROUP    --array=1-2940%40   hpc/submit_corrected_array.sh
#   sbatch --account=PROJECT_GROUP --array=2941-4200%30 hpc/submit_corrected_array.sh
# Writes to a NEW dir, outputs/design_corrected (the original power cells untouched).

set -euo pipefail
SUBMIT_DIR="$HOME/design_analysis"
cd "$SUBMIT_DIR"
echo "Job $SLURM_JOB_ID | task $SLURM_ARRAY_TASK_ID | acct ${SLURM_JOB_ACCOUNT:-?} | $(hostname) | $(date -Iseconds)"

module purge
ARC_R_MODULE="${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"
module load "$ARC_R_MODULE"

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
OUTPUT_DIR="${PROJECT_DATA}/outputs/design_corrected"
mkdir -p "$R_LIBS_USER" "$RENV_PATHS_CACHE" "$CMDSTANR_OUTPUT_DIR" "$OUTPUT_DIR" outputs/logs

GRID="${GRID:-config/design_grid_corrected.csv}"
ROW_INDEX="$SLURM_ARRAY_TASK_ID"
N_ROWS=$(tail -n +2 "$GRID" | wc -l)
if [[ "$ROW_INDEX" -gt "$N_ROWS" ]]; then
  echo "[corrected] task $ROW_INDEX exceeds grid size $N_ROWS; exiting cleanly."; exit 0
fi

SEED_COL=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="seed"){print i; exit}}' "$GRID")
SEED=$(awk -F',' -v row="$((ROW_INDEX+1))" -v c="${SEED_COL:-15}" 'NR==row {print $c}' "$GRID")
echo "Grid row $ROW_INDEX / $N_ROWS | $(awk -F',' -v row="$((ROW_INDEX+1))" 'NR==row{print $1,$2,"N="$5,"nverb="$6}' "$GRID") | seed $SEED"

Rscript scripts/04_design_analysis_cell.R \
  --row_index "$ROW_INDEX" \
  --grid      "$GRID" \
  --config    config/analysis_config.yaml \
  --outdir    "$OUTPUT_DIR" \
  ${OVERWRITE:+--overwrite}

EXIT_CODE=$?
echo "End $(date -Iseconds) | exit $EXIT_CODE"
exit $EXIT_CODE
