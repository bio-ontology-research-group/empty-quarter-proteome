#!/bin/bash
#SBATCH --job-name=f12_ref_ctrl
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch,debug
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/ref_ctrl_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/ref_ctrl_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
DIAMOND_DIR="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin"
export PATH="${DIAMOND_DIR}:${PATH}"

module load hmmer/3.3.2

echo "=== [Reference Control Audit] Started at $(date -u) ==="
echo "Host: $(hostname)"
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Python: $("$PYTHON" --version)"
echo "HMMER: $(hmmsearch -h | head -n 2 | tail -n 1)"
echo "DIAMOND: $(diamond version)"

./require_slurm_worker.sh "$PYTHON" scripts/test_references_control.py

echo "=== Reference Control Audit Complete at $(date -u) ==="
