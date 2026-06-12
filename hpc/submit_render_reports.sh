#!/bin/bash
#SBATCH --job-name=claps_render_reports
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/render_reports_%j.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/render_reports_%j.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_render_reports.sh
# Render all QMD reports after aggregation.
# Submit with afterok dependency on submit_aggregate_afterok.sh.
# References must be verified before rendering; this script will fail
# if the reference audit has not passed.

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

# Quarto is not an ARC module; use the standalone install on PATH (~/bin/quarto ->
# $DATA/PROJECT_GROUP/quarto-*). Its bundled deno SIGTRAPs on login nodes but renders
# fine on compute nodes, where this batch job runs. TMPDIR must be writable/exec.
export PATH="$HOME/bin:$PATH"
export QUARTO_PATH="$HOME/bin/quarto"
export TMPDIR="${TMPDIR:-$HOME/tmp}"
mkdir -p "$HOME/tmp"
echo "[render] quarto: $(command -v quarto || echo NOT-FOUND)"

# Match the array job's environment, and point the audit + reports at the
# project-storage output tree (where aggregation wrote the summary CSVs).
export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export CLAPS_OUTPUTS_ROOT="${PROJECT_DATA}/outputs"
mkdir -p outputs/logs "${CLAPS_OUTPUTS_ROOT}/reference_audit"

echo "[render] Verifying references..."
Rscript scripts/00_verify_references.R --bib references.bib

echo "[render] Rendering reports..."
Rscript scripts/07_render_reports.R --report all

echo "[render] Done."
echo "End time: $(date -Iseconds)"
