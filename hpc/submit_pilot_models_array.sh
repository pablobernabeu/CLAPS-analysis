#!/bin/bash
#SBATCH --job-name=claps_pilot_models
#SBATCH --partition=long
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-8
#SBATCH --output=/home/%u/design_analysis/outputs/logs/pilot_models_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/pilot_models_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_pilot_models_array.sh
# SLURM array job: fit pilot models across languages × prior regimes.
# Each task runs one language × prior_regime × threshold_mode combination.
# Array index 1..N maps to rows of config/pilot_model_grid.csv.
# CPU-only on ARC arc cluster.

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

# ---------------------------------------------------------------------------
# Modules
# ---------------------------------------------------------------------------
module purge
ARC_R_MODULE="${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"
module load "$ARC_R_MODULE"

export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export CMDSTANR_OUTPUT_DIR="${TMPDIR:-${PROJECT_DATA}/cmdstan_tmp}"
PILOT_OUT="${PROJECT_DATA}/outputs/pilot_models"
mkdir -p "$CMDSTANR_OUTPUT_DIR" "$PILOT_OUT" outputs/logs

# ---------------------------------------------------------------------------
# Read grid row
# ---------------------------------------------------------------------------
GRID="config/pilot_model_grid.csv"
if [[ ! -f "$GRID" ]]; then
  echo "[pilot] ERROR: Grid not found: $GRID" >&2
  exit 1
fi

ROW_INDEX="$SLURM_ARRAY_TASK_ID"
LANGUAGE=$(awk -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $1}' "$GRID")
PRIOR=$(awk    -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $2}' "$GRID")
THRESHOLD=$(awk -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $3}' "$GRID")

if [[ -z "$LANGUAGE" ]]; then
  echo "[pilot] ERROR: Empty row at index $ROW_INDEX" >&2
  exit 1
fi

echo "Language:    $LANGUAGE"
echo "Prior:       $PRIOR"
echo "Threshold:   $THRESHOLD"

# Check for existing output unless OVERWRITE=1
OUT_MARKER="${PILOT_OUT}/${LANGUAGE}_${PRIOR}_${THRESHOLD}_ladder_log.csv"
if [[ -f "$OUT_MARKER" && "${OVERWRITE:-0}" != "1" ]]; then
  echo "[pilot] Output exists, skipping: $OUT_MARKER"
  exit 0
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
Rscript scripts/02_fit_pilot_models.R \
  --language    "$LANGUAGE" \
  --prior       "$PRIOR" \
  --threshold   "$THRESHOLD" \
  --config      config/analysis_config.yaml \
  --outdir      "$PILOT_OUT" \
  --seed        2025

echo "[pilot] Completed: $LANGUAGE / $PRIOR / $THRESHOLD"
echo "End time: $(date -Iseconds)"
