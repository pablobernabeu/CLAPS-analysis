#!/bin/bash
#SBATCH --job-name=claps_safeguard
#SBATCH --partition=short
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/safeguard_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/safeguard_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_safeguard_array.sh
# Effect-size sensitivity ("safeguard") arm at the real 72-verb design: discounted
# focal effects (1.00 at low N to localise the minimum; 0.75 and 0.60 to quantify
# the N needed under conservative effects). Cells run to N=150 x 72 verbs, so the
# resource request is larger than the main run (24 GB, 12 h on short). The ACCOUNT
# and ARRAY range are passed at submit, e.g.:
#   sbatch --account=PROJECT_GROUP --array=1-440%40 hpc/submit_safeguard_array.sh
# Writes to a NEW dir, outputs/design_safeguard (everything else untouched).

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
OUTPUT_DIR="${PROJECT_DATA}/outputs/design_safeguard"
mkdir -p "$R_LIBS_USER" "$RENV_PATHS_CACHE" "$CMDSTANR_OUTPUT_DIR" "$OUTPUT_DIR" outputs/logs

GRID="${GRID:-config/design_grid_safeguard.csv}"
ROW_INDEX="$SLURM_ARRAY_TASK_ID"
N_ROWS=$(tail -n +2 "$GRID" | wc -l)
if [[ "$ROW_INDEX" -gt "$N_ROWS" ]]; then
  echo "[safeguard] task $ROW_INDEX exceeds grid size $N_ROWS; exiting cleanly."; exit 0
fi

SEED_COL=$(awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="seed"){print i; exit}}' "$GRID")
SEED=$(awk -F',' -v row="$((ROW_INDEX+1))" -v c="${SEED_COL:-15}" 'NR==row {print $c}' "$GRID")
echo "Grid row $ROW_INDEX / $N_ROWS | $(awk -F',' -v row="$((ROW_INDEX+1))" 'NR==row{print $1,"N="$5,"nverb="$6,"bsem="$8}' "$GRID") | seed $SEED"

Rscript scripts/04_design_analysis_cell.R \
  --row_index "$ROW_INDEX" \
  --grid      "$GRID" \
  --config    config/analysis_config.yaml \
  --outdir    "$OUTPUT_DIR" \
  ${OVERWRITE:+--overwrite}

EXIT_CODE=$?
echo "End $(date -Iseconds) | exit $EXIT_CODE"
exit $EXIT_CODE
