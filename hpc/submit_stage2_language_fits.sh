#!/bin/bash
#SBATCH --job-name=claps_stage2_lang
#SBATCH --partition=long
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# 16 CPUs for the registered 16 chains, one process per chain.
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/stage2_lang_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/stage2_lang_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_stage2_language_fits.sh
# Stage 2, step 3: one individual-language model per array task, task i fitting the
# i-th analysed variety in the order of config/stage2_languages.csv. Steps 1 and 2
# (stage2/01_check_and_exclude.R, stage2/02_affectedness.R) must have run, and the
# repository deployed to the cluster must be the current one
# (scripts/verify_arc_deployment.sh), since the registered priors are read under the
# regime name "primary".
#
#   sbatch --account=<project> --array=1-<number of analysed varieties> hpc/submit_stage2_language_fits.sh
#
# To resume a variety whose level was killed by the wall time, at the next level:
#   sbatch --account=<project> --array=<i> --export=ALL,START_LEVEL=<level> hpc/submit_stage2_language_fits.sh

set -euo pipefail
SUBMIT_DIR="$HOME/design_analysis"; cd "$SUBMIT_DIR"
echo "Job $SLURM_JOB_ID | task $SLURM_ARRAY_TASK_ID | acct ${SLURM_JOB_ACCOUNT:-?} | $(hostname) | $(date -Iseconds)"

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
Rscript stage2/03_fit_language_models.R --index "$SLURM_ARRAY_TASK_ID" ${EXTRA[@]+"${EXTRA[@]}"}

echo "End $(date -Iseconds) | exit $?"
