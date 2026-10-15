#!/bin/bash
#SBATCH --job-name=f1_full_pipeline
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/full_f1_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/full_f1_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
DIAMOND_DIR="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin"
MMSEQS_DIR="/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin"
export PATH="${MMSEQS_DIR}:${DIAMOND_DIR}:${PATH}"

module load hmmer/3.3.2

echo "=== [Full F1 Pipeline Execution] Started at $(date -u) ==="
echo "Host: $(hostname)"
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "CPUs allocated: ${SLURM_CPUS_PER_TASK}"
echo "Python: $("$PYTHON" --version)"
echo "HMMER: $(hmmsearch -h | head -n 2 | tail -n 1)"
echo "DIAMOND: $(diamond version)"
echo "MMseqs2: $(mmseqs version)"

echo "--- Stage 1: Reconcile Pilot & Finalize Track F2 Status ---"
./require_slurm_worker.sh "$PYTHON" scripts/01_reconcile_pilot_and_finalize_f2.py

echo "--- Stage 2: Track F1 Full Census (All 5,438 Representative Proteomes) ---"
./require_slurm_worker.sh "$PYTHON" scripts/02_f1_full_census.py

echo "--- Stage 3: Track F1 Environmental Comparators (SMAG & GEM Pools) ---"
./require_slurm_worker.sh "$PYTHON" scripts/03_f1_environmental_comparators.py

echo "--- Stage 4: Diversity, Novelty, Clustering & ProstT5 Manifest ---"
./require_slurm_worker.sh "$PYTHON" scripts/04_f1_diversity_novelty_and_manifest.py

echo "--- Stage 5: Full Audit, Checksumming & Results Synthesis ---"
./require_slurm_worker.sh "$PYTHON" scripts/05_f1_full_audit_and_report.py

echo "=== [Full F1 Pipeline Execution] Successfully Completed at $(date -u) ==="
