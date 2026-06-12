#!/bin/bash
#SBATCH --job-name=claps_prior_sensitivity
#SBATCH --partition=long
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-16
#SBATCH --output=/home/%u/design_analysis/outputs/logs/prior_sensitivity_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/prior_sensitivity_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_prior_sensitivity_array.sh
# Prior sensitivity analysis across all prior regimes × threshold modes × languages.
# Array index maps to rows of config/prior_sensitivity_grid.csv.
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

module purge
ARC_R_MODULE="${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"
module load "$ARC_R_MODULE"

export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export CMDSTANR_OUTPUT_DIR="${TMPDIR:-${PROJECT_DATA}/cmdstan_tmp}"
SENS_OUT="${PROJECT_DATA}/outputs/prior_sensitivity"
mkdir -p "$CMDSTANR_OUTPUT_DIR" "$SENS_OUT" outputs/logs

GRID="${GRID:-config/prior_sensitivity_grid.csv}"
ROW_INDEX="$SLURM_ARRAY_TASK_ID"
LANGUAGE=$(awk    -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $1}' "$GRID")
MODEL_LEVEL=$(awk -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $2}' "$GRID")

if [[ -z "$LANGUAGE" ]]; then
  echo "[sensitivity] ERROR: Empty row $ROW_INDEX" >&2; exit 1
fi

echo "Language:    $LANGUAGE"
echo "Model level: $MODEL_LEVEL"

OUT_MARKER="${SENS_OUT}/${LANGUAGE}_${MODEL_LEVEL}_sensitivity_summary.csv"
if [[ -f "$OUT_MARKER" && "${OVERWRITE:-0}" != "1" ]]; then
  echo "[sensitivity] Output exists, skipping."
  exit 0
fi

Rscript scripts/03_prior_sensitivity.R \
  --language    "$LANGUAGE" \
  --model_level "$MODEL_LEVEL" \
  --config      config/analysis_config.yaml \
  --outdir      "$SENS_OUT" \
  --seed        2025

echo "[sensitivity] Completed: $LANGUAGE / $MODEL_LEVEL"
echo "End time: $(date -Iseconds)"
