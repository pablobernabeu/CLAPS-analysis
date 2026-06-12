#!/bin/bash
#SBATCH --job-name=claps_highmem_rescue
#SBATCH --partition=long
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/highmem_rescue_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/highmem_rescue_%A_%a.err
#SBATCH --mail-type=FAIL,END
#SBATCH --array=1-50%5

# hpc/submit_highmem_rescue_htc_cpu.sh
# CPU high-memory rescue path on ARC HTC cluster.
# Used ONLY when ordinary ARC arc jobs fail with OOM.
# NO GPU resources are requested. This is a CPU-only job on HTC.
# Retry the same prespecified cell with more memory before considering model simplification.
#
# To submit to HTC instead of arc, change the partition:
#   --partition=long   (if HTC uses the same partition names)
# or add cluster flag:
#   #SBATCH --clusters=htc
#
# Failure log: outputs/failed_cells.txt should list row indices of OOM failures.
# Usage: sbatch --clusters=htc hpc/submit_highmem_rescue_htc_cpu.sh

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
echo "Cluster:      HTC CPU HIGH-MEMORY RESCUE"
echo "Start time:   $(date -Iseconds)"
echo "Git SHA:      $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "NOTE: CPU-only. No GPU resources requested."
echo "=========================================="

module purge
ARC_R_MODULE="${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"
module load "$ARC_R_MODULE"

export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export CMDSTANR_OUTPUT_DIR="${TMPDIR:-${PROJECT_DATA}/cmdstan_tmp}"
DESIGN_OUT="${PROJECT_DATA}/outputs/design_analysis"
mkdir -p "$CMDSTANR_OUTPUT_DIR" "$DESIGN_OUT" outputs/logs

# Read failed row index from failed_cells list
FAILED_LIST="${FAILED_LIST:-${PROJECT_DATA}/outputs/failed_cells.txt}"
if [[ ! -f "$FAILED_LIST" ]]; then
  echo "[rescue] ERROR: Failed cells list not found: $FAILED_LIST" >&2
  exit 1
fi

ROW_INDEX=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FAILED_LIST")
if [[ -z "$ROW_INDEX" ]]; then
  echo "[rescue] Array task $SLURM_ARRAY_TASK_ID has no matching failed row; exiting."
  exit 0
fi

echo "[rescue] Retrying design grid row: $ROW_INDEX on HTC high-memory CPU node"
echo "[rescue] Memory available: ${SLURM_MEM_PER_NODE:-256G} MB"

# Force overwrite so the rescue attempt replaces the OOM record
Rscript scripts/04_design_analysis_cell.R \
  --row_index "$ROW_INDEX" \
  --grid      "${GRID:-config/design_grid.csv}" \
  --config    config/analysis_config.yaml \
  --outdir    "$DESIGN_OUT" \
  --overwrite

EXIT_CODE=$?
if [[ $EXIT_CODE -eq 0 ]]; then
  echo "[rescue] SUCCESS: Row $ROW_INDEX completed on HTC."
  # Remove from failed list
  grep -v "^${ROW_INDEX}$" "$FAILED_LIST" > "${FAILED_LIST}.tmp" \
    && mv "${FAILED_LIST}.tmp" "$FAILED_LIST"
else
  echo "[rescue] FAILED again: Row $ROW_INDEX. Manual intervention required." >&2
  echo "RESCUE_FAILED:${ROW_INDEX}" >> outputs/logs/rescue_failures.log
fi

echo "End time: $(date -Iseconds)"
exit $EXIT_CODE
