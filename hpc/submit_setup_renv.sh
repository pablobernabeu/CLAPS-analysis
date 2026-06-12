#!/bin/bash
#SBATCH --job-name=claps_setup
#SBATCH --partition=short
#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/home/%u/design_analysis/outputs/logs/setup_%j.out
#SBATCH --error=/home/%u/design_analysis/outputs/logs/setup_%j.err

# One-shot setup job: installs all renv packages + cmdstan on a compute node.
# Run ONCE before submitting any analysis arrays.
# Re-running is safe (renv::restore is idempotent; cmdstan skips if already present).
#
# Submit from design_analysis/ root:
#   cd ~/design_analysis && sbatch hpc/submit_setup_renv.sh

set -euo pipefail
cd ~/design_analysis
mkdir -p outputs/logs

echo "=========================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Host:       $(hostname)"
echo "Start time: $(date -Iseconds)"
echo "=========================================="

module purge
module load R/4.4.2-gfbf-2024a

if [[ -z "${DATA:-}" ]]; then
  echo "ERROR: \$DATA is not set. Cannot write to project storage." >&2
  exit 1
fi
PROJECT_DATA="${DATA}/PROJECT_GROUP"

export R_LIBS_USER="${PROJECT_DATA}/R/library_4.4"
export RENV_PATHS_CACHE="${PROJECT_DATA}/renv/cache"
export CMDSTAN_INSTALL_DIR="${PROJECT_DATA}/cmdstan"
mkdir -p "$R_LIBS_USER" "$RENV_PATHS_CACHE" "$CMDSTAN_INSTALL_DIR"

echo ""
echo "=== Step 1a: bootstrap renv ==="
Rscript --no-save -e '
  lib <- Sys.getenv("R_LIBS_USER")
  .libPaths(c(lib, .libPaths()))
  if (!requireNamespace("renv", quietly = TRUE)) {
    cat("Installing renv into", lib, "\n")
    install.packages("renv", lib = lib, repos = "https://cloud.r-project.org")
  } else {
    cat("renv already present:", as.character(packageVersion("renv")), "\n")
  }
'

echo ""
echo "=== Step 1b: renv::restore() ==="
Rscript --no-save -e '
  lib <- Sys.getenv("R_LIBS_USER")
  .libPaths(c(lib, .libPaths()))
  library(renv, lib.loc = lib)
  renv::restore(prompt = FALSE, library = lib)
  cat("renv::restore() complete\n")
'

echo ""
echo "=== Step 2: install cmdstan ==="
Rscript --no-save -e '
  lib <- Sys.getenv("R_LIBS_USER")
  .libPaths(c(lib, .libPaths()))
  library(cmdstanr, lib.loc = lib)
  already_installed <- tryCatch({
    p <- cmdstanr::cmdstan_path()
    nchar(p) > 0
  }, error = function(e) FALSE)
  if (!already_installed) {
    cmdstanr::install_cmdstan(cores = 4, overwrite = FALSE,
      dir = Sys.getenv("CMDSTAN_INSTALL_DIR"))
  } else {
    cat("cmdstan already installed at:", cmdstanr::cmdstan_path(), "\n")
  }
  cat("cmdstan version:", cmdstanr::cmdstan_version(), "\n")
'

echo ""
echo "=== Step 3: smoke check ==="
Rscript --no-save -e '
  lib <- Sys.getenv("R_LIBS_USER")
  .libPaths(c(lib, .libPaths()))
  library(brms,      lib.loc = lib)
  library(cmdstanr,  lib.loc = lib)
  cat("brms version:   ", as.character(packageVersion("brms")), "\n")
  cat("cmdstan version:", cmdstanr::cmdstan_version(), "\n")
  cat("Setup complete\n")
'

echo ""
echo "End time: $(date -Iseconds)"
