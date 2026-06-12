#!/bin/bash
#SBATCH --job-name=claps_design_analysis
#SBATCH --partition=long
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-39%39   # Default: rows 1-39 (single-language: English 20 +
                          # Turkish 14 + Norwegian 5). Cross-language rows 40-51
                          # (AllLanguages) are now supported via
                          # simulate_claps_data_multilanguage(); run them with
                          # `sbatch --array=40-51 hpc/submit_design_analysis_array.sh`.
#SBATCH --output=/home/%u/design_analysis/outputs/logs/design_analysis_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/design_analysis_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_design_analysis_array.sh
# Bayes-factor design analysis: each array task = one cell in config/design_grid.csv.
# 39 single-language cells run in parallel on ARC arc, long partition.
# AllLanguages cells (rows 40-51) use simulate_claps_data_multilanguage() and the
# cross-language ladder (by-Language random effects); submit with --array=40-51.
# Wall time 14 days: the long partition has unlimited max; previous analysis
# exceeded the prior 7-day cap before all chains completed.
#
# Submit from design_analysis/ root:
#   cd /path/to/repo/design_analysis && sbatch hpc/submit_design_analysis_array.sh

set -euo pipefail

# Under sbatch, BASH_SOURCE can point to a spool copy in a non-writable path.
# Always anchor paths to the submit directory instead.
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

if [[ -z "${DATA:-}" ]]; then
  echo "ERROR: \$DATA is not set. Cannot write to project storage." >&2
  exit 1
fi
PROJECT_DATA="${DATA}/PROJECT_GROUP"

export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export CMDSTANR_OUTPUT_DIR="${TMPDIR:-${PROJECT_DATA}/cmdstan_tmp}"
OUTPUT_DIR="${PROJECT_DATA}/outputs/design_analysis"
mkdir -p "$R_LIBS_USER" "$RENV_PATHS_CACHE" "$CMDSTANR_OUTPUT_DIR" "$OUTPUT_DIR" outputs/logs

# Compute grid row from array task ID. GRID is overridable at submit time, e.g.
#   sbatch --export=ALL,GRID=config/design_grid_single.csv --array=1-4200 ...
GRID="${GRID:-config/design_grid.csv}"
ROW_INDEX="$SLURM_ARRAY_TASK_ID"
N_ROWS=$(tail -n +2 "$GRID" | wc -l)

if [[ "$ROW_INDEX" -gt "$N_ROWS" ]]; then
  echo "[design] Array task $ROW_INDEX exceeds grid size $N_ROWS; exiting cleanly."
  exit 0
fi

echo "Grid row: $ROW_INDEX / $N_ROWS"

# Record seed and model info for log. The seed column is looked up by name so
# that adding grid columns (e.g. include_gender, beta_gender) does not shift it
# (it is no longer the last field).
SEED_COL=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="seed"){print i; exit}}' "$GRID")
SEED=$(awk    -F',' -v row="$((ROW_INDEX+1))" -v c="${SEED_COL:-15}" 'NR==row {print $c}' "$GRID")
LANGUAGE=$(awk -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $1}' "$GRID")
MODEL=$(awk   -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $2}' "$GRID")
PRIOR=$(awk   -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $3}' "$GRID")
THRESH=$(awk  -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $4}' "$GRID")

echo "Language:   $LANGUAGE | Model: $MODEL | Prior: $PRIOR | Threshold: $THRESH | Seed: $SEED"

Rscript scripts/04_design_analysis_cell.R \
  --row_index "$ROW_INDEX" \
  --grid      "$GRID" \
  --config    config/analysis_config.yaml \
  --outdir    "$OUTPUT_DIR" \
  ${OVERWRITE:+--overwrite}

EXIT_CODE=$?
echo "End time: $(date -Iseconds)"
echo "Exit code: $EXIT_CODE"
exit $EXIT_CODE
