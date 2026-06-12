#!/bin/bash
#SBATCH --job-name=claps_aggregate
#SBATCH --partition=short
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/aggregate_%j.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/aggregate_%j.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_aggregate_afterok.sh
# Aggregate design results after all design array jobs complete.
# Submit with afterok dependency on the design-analysis and calibration job arrays:
#   DESIGN_JOB=$(sbatch hpc/submit_design_analysis_array.sh | awk '{print $4}')
#   CALIB_JOB=$(sbatch hpc/submit_bf_calibration_array.sh | awk '{print $4}')
#   sbatch --dependency=afterok:${DESIGN_JOB}:${CALIB_JOB} hpc/submit_aggregate_afterok.sh

set -euo pipefail

# Under sbatch, BASH_SOURCE can point to a spool copy in a non-writable path, so
# anchor to the submit directory explicitly (mirrors submit_design_analysis_array.sh).
if [[ -z "${DATA:-}" ]]; then
  echo "ERROR: \$DATA is not set. Cannot locate project storage." >&2
  exit 1
fi
PROJECT_DATA="${DATA}/PROJECT_GROUP"
SUBMIT_DIR="$HOME/design_analysis"
cd "$SUBMIT_DIR"

echo "=========================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Host:       $(hostname)"
echo "Start time: $(date -Iseconds)"
echo "=========================================="

module purge
ARC_R_MODULE="${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"
module load "$ARC_R_MODULE"

# Match the array job's environment so the same R library and output tree are used.
export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Cells write to (and the summary is read from) the project-storage output tree,
# not the home directory; point the aggregator at it explicitly.
OUTPUT_DIR="${PROJECT_DATA}/outputs/design_analysis"
SUM_DIR="${PROJECT_DATA}/outputs/design_summary"
mkdir -p outputs/logs "$SUM_DIR"

echo "[aggregate] Aggregating design results from $OUTPUT_DIR ..."
Rscript scripts/06_aggregate_design_results.R --out_dir "$OUTPUT_DIR" --sum_dir "$SUM_DIR"

echo "[aggregate] Running status report..."
Rscript scripts/08_submit_status_report.R || \
  echo "[aggregate] status report step failed (non-fatal)"

echo "[aggregate] Done."
echo "End time: $(date -Iseconds)"
