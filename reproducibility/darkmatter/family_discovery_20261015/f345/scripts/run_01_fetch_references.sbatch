#!/bin/bash
#SBATCH --job-name=f345_01_refs
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug,batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/01_refs_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/logs/01_refs_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
DIAMOND_DIR="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin"
export PATH="${DIAMOND_DIR}:${PATH}"

module load hmmer/3.3.2

echo "=== [F3/F4/F5 Pilot Stage 1] Fetch References Started at $(date -u) ==="
echo "Host: $(hostname)"
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "Python: $("$PYTHON" --version)"
echo "HMMER: $(hmmsearch -h | head -n 2 | tail -n 1)"
echo "DIAMOND: $(diamond version)"

./require_slurm_worker.sh "$PYTHON" scripts/01_fetch_references.py

echo "=== Stage 1 Complete at $(date -u) ==="
