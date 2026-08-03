#!/bin/bash
#SBATCH --job-name=claps_pooled2
#SBATCH --partition=medium
# 5 days. Raised from 1-18:00:00 (42 h) on 2026-07-31 after measuring what these
# fits actually cost: across 331 completed pooled/decision cells on htc the range
# was 1-04:37 to 4-00:00, i.e. 29 to 96 hours. A 42-hour limit therefore killed a
# large share of the expensive cells, and did so silently from the grid's point of
# view — a timed-out task leaves no .rds, so the cell simply looks "not yet run".
# The array currently in flight (8377361) survives only because it was submitted
# with an explicit --time=10-00:00:00 override; the script default would have lost
# its N=130 and N=150 cells. See the runtime note below.
#SBATCH --time=5-00:00:00
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
# MEASURED RUNTIME (2026-07-31). A pooled fit at N=150 per language holds about
# 85k rows plus a third grouping factor. Measured across 331 completed cells:
#
#   minimum   1-04:37  (29 h)
#   maximum   4-00:00  (96 h)
#
# That is three to eight times the "8-12 h per fit" figure quoted elsewhere in the
# repository, which was an early estimate and is corrected in the same commit as
# this note. Budget accordingly: at 30 concurrent tasks the 150-cell grid is a
# multi-week run, not a multi-day one, and the expensive N=130/150 cells dominate.
#
# The medium partition is required for this reason; short's 12 h is nowhere near
# enough.
#
#   sbatch --clusters=htc --account=PROJECT_GROUP --array=1-75%20  hpc/submit_pooled_v2_array.sh
#   sbatch --clusters=htc --account=PROJECT_GROUP    --array=76-150%20 hpc/submit_pooled_v2_array.sh
#
# Smoke test (isolated grid AND output dir, short partition):
#   GRID=config/design_grid_pooled_v2_TEST.csv \
#   OUTPUT_DIR=$DATA/PROJECT_GROUP/outputs/design_pooled_test \
#   sbatch --clusters=htc --account=PROJECT_GROUP --partition=short --time=2:00:00 \
#     --export=ALL,GRID,OUTPUT_DIR --array=1 hpc/submit_pooled_v2_array.sh
#
# Note the --export=ALL,GRID,OUTPUT_DIR in the smoke-test form: SLURM does not pass
# the submitting shell's variables to the job unless they are named, so setting
# GRID= on the sbatch line without exporting it has no effect inside the script.
#
# The shared preamble below (module handling, R_LIBS_USER, thread limits, CmdStan
# discovery, the over-wide-array clean exit) is explained line by line under
# "Anatomy of a submission script" in docs/arc_submission_guide.md.

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
# Keep every scratch file this job writes on project storage rather than on the
# node's /tmp. cmdstanr puts the COMPILED MODEL in R's tempdir(), which follows
# TMPDIR, and a pooled cell can run for days: cell N130 seed 931000 died on
# 2026-08-03 after 121 hours with "File does not exist:
# '/tmp/RtmpW8tYOr/model_5ae97f0f...'", the compiled model having been removed
# from /tmp underneath the running job. Setting TMPDIR here fixes the model, the
# sampler's CSVs and anything else R writes to tempdir(), in one place. It also
# makes the fallback on the next line unconditional, which is deliberate: the
# node-local default is exactly what fails on the longest cells.
export TMPDIR="${PROJECT_DATA}/cmdstan_tmp"
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
