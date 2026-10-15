#!/bin/bash
#SBATCH --job-name=f4_full
#SBATCH --partition=debug,batch
#SBATCH --account=pi-hohndor
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/09_f4_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/09_f4_%j.err

set -euo pipefail

SCRIPT_DIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345"
GUARD="${SCRIPT_DIR}/require_slurm_worker.sh"
PY_ENV="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"

module load hmmer/3.3.2

echo "=== [09 Full F4 Ectoine Pipeline] Started at $(date -u) ==="
echo "Host: $(hostname -s)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

"${GUARD}" "${PY_ENV}" "${SCRIPT_DIR}/scripts/09_f4_ectc_full.py"

echo "=== [09 Full F4 Ectoine Pipeline] Completed at $(date -u) ==="
