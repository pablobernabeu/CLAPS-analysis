#!/bin/bash
#SBATCH --job-name=claps_glossatau
# long, not medium. The published fit of this model took 490,250 s elapsed (136 h)
# with 16 chains under rstan, and the comparable CLAPS pooled fits at a similar
# number of observations run 29-96 h with 4 chains. 10 days leaves headroom for the
# two corrected arms, which carry an extra grouping factor.
#SBATCH --partition=long
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
# One process per chain, no within-chain threading, so 4 CPUs for 4 chains.
# Memory: the first attempt (job 13106297) requested 16G and every arm was
# OOM-killed after sampling, peaking at 15.4 GiB inside brms::loo(). LOO is no
# longer computed inline, but the post-sampling summaries still allocate on the
# scale of the draws array, so the request is doubled rather than trimmed back to
# the sampling peak.
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --array=1-3
#SBATCH --output=/home/%u/design_analysis/outputs/logs/glossatau_%A_%a.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/glossatau_%A_%a.err
#SBATCH --mail-type=FAIL,END

# hpc/submit_glossa_tau_refit.sh
# Re-estimate the between-language SD of the interaction from the Ambridge, Arnon
# & Bekman (2023) five-language data under a correctly specified model.
# One array task per arm: 1 M0_published, 2 M1_verbsplit, 3 M2_thresholds.
#
#   sbatch hpc/submit_glossa_tau_refit.sh
#
# The five per-language CSVs must be at $SUBMIT_DIR/data/glossa2023/ before this
# runs; they are copied by hand like everything else in this project (see
# scripts/verify_arc_deployment.sh).

set -euo pipefail
SUBMIT_DIR="$HOME/design_analysis"; cd "$SUBMIT_DIR"
echo "Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} | acct ${SLURM_JOB_ACCOUNT:-?} | $(hostname) | $(date -Iseconds)"

module purge
module load "${ARC_R_MODULE:-R/4.4.2-gfbf-2024a}"

if [[ -z "${DATA:-}" ]]; then echo "ERROR: \$DATA is not set." >&2; exit 1; fi
PROJECT_DATA="${DATA}/PROJECT_GROUP"
export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
# Persistent, per-arm, on shared storage. NOT $TMPDIR: the auto_tmpdir SPANK
# plugin (see /etc/slurm/plugstack.conf) gives every job a private bind-mounted
# /tmp that is destroyed when the job ends, so draws written there die with the
# job. The first attempt lost three arms of finished sampling that way. The R
# script turns this into options(cmdstanr_output_dir=), which is what cmdstanr
# actually reads -- the exported variable alone has never had any effect.
export CMDSTANR_OUTPUT_DIR="${PROJECT_DATA}/outputs/glossa_tau_refit/cmdstan_csv/arm${SLURM_ARRAY_TASK_ID:-0}"
export CMDSTAN="$(ls -d "${PROJECT_DATA}/cmdstan/cmdstan-"* 2>/dev/null | sort -V | tail -1)"
if [[ -z "${CMDSTAN}" || ! -x "${CMDSTAN}/bin/stanc" ]]; then
  echo "ERROR: CmdStan not found under ${PROJECT_DATA}/cmdstan." >&2; exit 1
fi
mkdir -p "$CMDSTANR_OUTPUT_DIR" outputs/logs

if [[ ! -f data/glossa2023/English_Passives.csv ]]; then
  echo "ERROR: data/glossa2023/ is missing the Ambridge et al. CSVs." >&2; exit 1
fi

Rscript scripts/refit_glossa_pooled_tau.R \
  --datadir data/glossa2023 \
  --outdir "${PROJECT_DATA}/outputs/glossa_tau_refit"

echo "End $(date -Iseconds) | exit $?"
