#!/bin/bash
#SBATCH --job-name=claps_stage2_single
#SBATCH --partition=long
#SBATCH --time=14-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# 16 CPUs for the registered 16 chains. The single model holds every language's
# data at once, hence the larger memory request.
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/stage2_single_%j.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/stage2_single_%j.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_stage2_single_model.sh
# Stage 2, step 5: the single model across languages, by the registered ladder
# (stage2/05_robustness_single_model.R). Needs steps 1 and 2; independent of step 3.
#
#   sbatch --account=<project> hpc/submit_stage2_single_model.sh
#
# A level killed by the wall time is a permitted reason to descend; resume at the
# next level with
#   sbatch --account=<project> --export=ALL,START_LEVEL=<level> hpc/submit_stage2_single_model.sh

set -euo pipefail
SUBMIT_DIR="$HOME/design_analysis"; cd "$SUBMIT_DIR"
echo "Job $SLURM_JOB_ID | acct ${SLURM_JOB_ACCOUNT:-?} | $(hostname) | $(date -Iseconds)"

module purge
module load "${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"

if [[ -z "${DATA:-}" ]]; then echo "ERROR: \$DATA is not set." >&2; exit 1; fi
PROJECT_DATA="${DATA}/PROJECT_GROUP"
export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export STAN_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"
export CMDSTAN="$(ls -d "${PROJECT_DATA}/cmdstan/cmdstan-"* 2>/dev/null | sort -V | tail -1)"
if [[ -z "${CMDSTAN}" || ! -x "${CMDSTAN}/bin/stanc" ]]; then
  echo "ERROR: CmdStan not found under ${PROJECT_DATA}/cmdstan." >&2; exit 1
fi
mkdir -p outputs/logs

EXTRA=()
if [[ -n "${START_LEVEL:-}" ]]; then EXTRA=(--start-level "$START_LEVEL"); fi
Rscript stage2/05_robustness_single_model.R ${EXTRA[@]+"${EXTRA[@]}"}

echo "End $(date -Iseconds) | exit $?"
