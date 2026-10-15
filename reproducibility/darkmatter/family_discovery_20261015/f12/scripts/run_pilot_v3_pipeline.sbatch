#!/bin/bash
#SBATCH --job-name=f12_pilot_v3
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug,batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/pilot_v3_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/pilot_v3_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12"
cd "${WORKDIR}"

PYTHON="/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python"
DIAMOND_DIR="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin"
export PATH="${DIAMOND_DIR}:${PATH}"

module load hmmer/3.3.2

echo "=== [F1/F2 Feasibility Pilot v3] Started at $(date -u) ==="
echo "Host: $(hostname)"
echo "Slurm Job ID: ${SLURM_JOB_ID}"
echo "CPUs allocated: ${SLURM_CPUS_PER_TASK}"
echo "Python: $("$PYTHON" --version)"
echo "HMMER: $(hmmsearch -h | head -n 2 | tail -n 1)"
echo "DIAMOND: $(diamond version)"

echo "--- Stage 1: Reference Prep & 256-Cohort Deterministic Selection ---"
./require_slurm_worker.sh "$PYTHON" scripts/01b_prep_dedup_and_cohort.py

echo "--- Stage 2: Track F1 Photolyase/Cryptochrome Pilot Screen ---"
./require_slurm_worker.sh "$PYTHON" scripts/02_f1_pilot_v3.py

echo "--- Stage 3: Track F2 Phosphatase Remapping & Competing HAD Test ---"
./require_slurm_worker.sh "$PYTHON" scripts/03_f2_pilot_v3.py

echo "--- Stage 4: Pilot Audit, Checksums, and Verdict Summarization ---"
./require_slurm_worker.sh "$PYTHON" scripts/04_audit_and_summarize_v3.py

echo "=== [F1/F2 Feasibility Pilot v3] Successfully Completed at $(date -u) ==="
