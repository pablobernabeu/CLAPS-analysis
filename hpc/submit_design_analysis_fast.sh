#!/bin/bash
#SBATCH --job-name=claps_da_fast
#SBATCH --partition=medium
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-39%39
#SBATCH --output=/home/%u/design_analysis/outputs/logs/da_fast_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/da_fast_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_design_analysis_fast.sh
# Fast-screening tier of the Bayes-factor design analysis.
# Runs all 39 single-language cells (English + Turkish + Norwegian) at
# reduced sampling (iter=1500, warmup=750, chains=4) to obtain initial BF
# operating-characteristic estimates before committing to full production.
#
# Seeds are offset +10000 from the production seeds so outputs are distinct.
# AllLanguages cells (rows 40-51) are excluded until simulate_claps_data_multilanguage()
# is implemented (current gap documented in R/06_simulate_design.R).
#
# Partition: medium (max 2 days) — typical wall time per cell ~30 min – 4 hours.
#
# Submit from design_analysis/ root:
#   cd /path/to/repo/design_analysis && sbatch hpc/submit_design_analysis_fast.sh

set -euo pipefail

# Under sbatch, BASH_SOURCE can point to a spool copy in a non-writable path.
# Always anchor paths to the submit directory instead.
SUBMIT_DIR="$HOME/design_analysis"
cd "$SUBMIT_DIR"
mkdir -p outputs/logs

echo "=========================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Array task:   $SLURM_ARRAY_TASK_ID"
echo "Host:         $(hostname)"
echo "Tier:         FAST SCREENING (iter=1500, chains=4)"
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
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export CMDSTANR_OUTPUT_DIR="${TMPDIR:-${PROJECT_DATA}/cmdstan_tmp}"
OUTPUT_DIR="${PROJECT_DATA}/outputs/design_analysis_fast"
mkdir -p "$R_LIBS_USER" "$RENV_PATHS_CACHE" "$CMDSTANR_OUTPUT_DIR" "$OUTPUT_DIR" outputs/logs

GRID="config/design_grid_fast.csv"
ROW_INDEX="$SLURM_ARRAY_TASK_ID"
N_ROWS=$(tail -n +2 "$GRID" | wc -l)

if [[ "$ROW_INDEX" -gt "$N_ROWS" ]]; then
  echo "[da_fast] Array task $ROW_INDEX exceeds grid size $N_ROWS; exiting."
  exit 0
fi

SEED=$(    awk -F',' -v r="$((ROW_INDEX+1))" 'NR==r {print $NF}' "$GRID")
LANGUAGE=$(awk -F',' -v r="$((ROW_INDEX+1))" 'NR==r {print $1}'  "$GRID")
MODEL=$(   awk -F',' -v r="$((ROW_INDEX+1))" 'NR==r {print $2}'  "$GRID")
PRIOR=$(   awk -F',' -v r="$((ROW_INDEX+1))" 'NR==r {print $3}'  "$GRID")

echo "Language: $LANGUAGE | Model: $MODEL | Prior: $PRIOR | Seed: $SEED"

Rscript scripts/04_design_analysis_cell.R \
  --row_index "$ROW_INDEX" \
  --grid      "$GRID" \
  --config    config/analysis_config.yaml \
  --outdir    "$OUTPUT_DIR" \
  ${OVERWRITE:+--overwrite}

EXIT_CODE=$?
echo "End time: $(date -Iseconds) | Exit: $EXIT_CODE"
exit $EXIT_CODE
