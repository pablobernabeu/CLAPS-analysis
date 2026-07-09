#!/bin/bash
#SBATCH --job-name=claps_pooled2
#SBATCH --partition=medium
#SBATCH --time=1-18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/pooled2_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/pooled2_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_pooled_v2_array.sh
# Pooled cross-language pilot-grounded design analysis (assurance mode).
# Depends on the three per-language v2 pilot DGPs
# (pilot_dgp_v2_pilot_{English,Turkish,Norwegian}.rds) already present under
# outputs/pilot_models. Own grid, runner, and output dir, so nothing shared
# with the running single-language sweeps.
#
# The medium partition (42 h here) is required: a pooled fit at N=150 per
# language holds about 85k rows plus a third grouping factor, comparable to the
# single-language N=400-500 cells that exceeded short's 12 h.
#
#   sbatch --clusters=htc --account=PROJECT_GROUP --array=1-75%20  hpc/submit_pooled_v2_array.sh
#   sbatch --clusters=htc --account=PROJECT_GROUP    --array=76-150%20 hpc/submit_pooled_v2_array.sh
#
# Smoke test (isolated grid AND output dir, short partition):
#   GRID=config/design_grid_pooled_v2_TEST.csv \
#   OUTPUT_DIR=$DATA/PROJECT_GROUP/outputs/design_pooled_test \
#   sbatch --clusters=htc --account=PROJECT_GROUP --partition=short --time=2:00:00 \
#     --export=ALL,GRID,OUTPUT_DIR --array=1 hpc/submit_pooled_v2_array.sh

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
OUTPUT_DIR="${OUTPUT_DIR:-${PROJECT_DATA}/outputs/design_pooled_v2}"
DGP_DIR="${PROJECT_DATA}/outputs/pilot_models"
mkdir -p "$OUTPUT_DIR" "$CMDSTANR_OUTPUT_DIR" outputs/logs

GRID="${GRID:-config/design_grid_pooled_v2.csv}"
ROW="$SLURM_ARRAY_TASK_ID"
N_ROWS=$(tail -n +2 "$GRID" | wc -l)
if [[ "$ROW" -gt "$N_ROWS" ]]; then
  echo "[pooled2] task $ROW exceeds grid size $N_ROWS; exiting cleanly."; exit 0
fi
echo "Grid row $ROW / $N_ROWS | $(awk -F',' -v r="$((ROW+1))" 'NR==r{print "N/lang="$1,"mode="$2}' "$GRID")"

Rscript scripts/07_pooled_cell_v2.R \
  --row_index "$ROW" \
  --grid      "$GRID" \
  --dgpdir    "$DGP_DIR" \
  --outdir    "$OUTPUT_DIR" \
  ${OVERWRITE:+--overwrite}

echo "End $(date -Iseconds) | exit $?"
