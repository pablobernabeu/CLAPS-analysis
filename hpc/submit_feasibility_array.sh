#!/bin/bash
#SBATCH --job-name=claps_feasibility
#SBATCH --partition=long
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=1-72%10
#SBATCH --output=/home/%u/design_analysis/outputs/logs/feasibility_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/feasibility_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_feasibility_array.sh
# Cross-language convergence/timing FEASIBILITY curve (NOT a power run).
# Each array task = one cell in config/design_grid_feasibility.csv: the full
# cross ladder (L0..L5) x simulated language counts (3, 10, 20) at N=80 with the
# within-subjects gender covariate on. Maps where convergence and/or the ARC
# walltime break before data collection.
#
# Writes to a SEPARATE output directory (design_analysis_feasibility_v2) so the
# 9200 power-analysis cells (outputs/design_analysis) and the earlier
# feasibility_v1 run are left completely untouched.
#
# Submit from design_analysis/ root:
#   cd "$HOME/design_analysis" && sbatch hpc/submit_feasibility_array.sh
# Stage the cheap floor cells first if preferred (L0_cross at 3/10/20 langs x3 seeds):
#   sbatch --array=1-9 hpc/submit_feasibility_array.sh
# 64 GB matches the production cross fits; cells that OOM (no .rds written) are a
# feasibility finding and can be re-run on the 256 GB HTC high-memory rescue tier.

set -euo pipefail

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
# CmdStan is installed in the project data area (by submit_setup_renv.sh), not in
# cmdstanr's default ~/.cmdstan, so point cmdstanr/brms at it explicitly. Fail
# fast with a clear message if it is missing, rather than letting every cell error
# with an opaque "CmdStan path has not been set yet".
export CMDSTAN="$(ls -d "${PROJECT_DATA}/cmdstan/cmdstan-"* 2>/dev/null | sort -V | tail -1)"
if [[ -z "${CMDSTAN}" || ! -x "${CMDSTAN}/bin/stanc" ]]; then
  echo "ERROR: CmdStan not found under ${PROJECT_DATA}/cmdstan. Run 'sbatch hpc/submit_setup_renv.sh' first." >&2
  exit 1
fi
echo "CmdStan: $CMDSTAN"
OUTPUT_DIR="${PROJECT_DATA}/outputs/design_analysis_feasibility_v2"
mkdir -p "$R_LIBS_USER" "$RENV_PATHS_CACHE" "$CMDSTANR_OUTPUT_DIR" "$OUTPUT_DIR" outputs/logs

# GRID overridable at submit time; defaults to the feasibility grid.
GRID="${GRID:-config/design_grid_feasibility.csv}"
ROW_INDEX="$SLURM_ARRAY_TASK_ID"
N_ROWS=$(tail -n +2 "$GRID" | wc -l)

if [[ "$ROW_INDEX" -gt "$N_ROWS" ]]; then
  echo "[feasibility] Array task $ROW_INDEX exceeds grid size $N_ROWS; exiting cleanly."
  exit 0
fi

echo "Grid row: $ROW_INDEX / $N_ROWS"

# Look up columns by name so an added n_languages column does not shift them.
SEED_COL=$(awk  -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="seed"){print i; exit}}' "$GRID")
NLANG_COL=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="n_languages"){print i; exit}}' "$GRID")
SEED=$(awk     -F',' -v row="$((ROW_INDEX+1))" -v c="${SEED_COL:-15}"  'NR==row {print $c}' "$GRID")
NLANG=$(awk    -F',' -v row="$((ROW_INDEX+1))" -v c="${NLANG_COL:-18}" 'NR==row {print $c}' "$GRID")
LANGUAGE=$(awk -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $1}' "$GRID")
MODEL=$(awk    -F',' -v row="$((ROW_INDEX+1))" 'NR==row {print $2}' "$GRID")

echo "Language: $LANGUAGE | Model: $MODEL | n_languages: $NLANG | Seed: $SEED"

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
