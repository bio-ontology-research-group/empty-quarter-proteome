#!/bin/bash
#SBATCH --job-name=h4_v7_recovery
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/h4_v7_recovery_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/logs/h4_v7_recovery_%j.err

set -euo pipefail

echo "[$(date -Iseconds)] Starting Hypothesis 4 V7 Recovery & Post-Processing on Ibex"
echo "Job ID: ${SLURM_JOB_ID:-LOCAL}"
echo "Node: ${SLURMD_NODENAME:-$(hostname)}"
echo "Allocated CPUs: ${SLURM_CPUS_PER_TASK:-16}"

W=/ibex/scratch/projects/c2014/rob/metagenomics
H4_DIR=$W/06_darkmatter/p7_analyses_clean/h4_structure
mkdir -p "$H4_DIR/logs" "$H4_DIR/results"

PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python
module load foldseek/10-941cd33

echo "[$(date -Iseconds)] Foldseek version:"
foldseek version

echo "[$(date -Iseconds)] Python environment:"
$PY --version

# Previous initial run log from job 53738540
PREV_LOG="$H4_DIR/logs/h4_v7_53738540.out"

# Execute recovery post-processing wrapper
echo "[$(date -Iseconds)] Executing v7_recovery_postprocess.py..."
$PY "$H4_DIR/scripts/v7_recovery_postprocess.py" --threads "${SLURM_CPUS_PER_TASK:-16}" --log "$PREV_LOG"

# Verify all checksums and closure assertions on IBEX in strict mode
echo "[$(date -Iseconds)] Running independent strict closure verification on IBEX..."
$PY "$H4_DIR/scripts/v7_verify_closure.py" --mode strict --dir "$H4_DIR/results"

echo "[$(date -Iseconds)] SUCCESS: Hypothesis 4 V7 recovery and post-processing complete."
