#!/bin/bash
# Wrapper for genecall_reconcile.py (revision item 3b). Resources are set per stage at submission
# by genecall_reconcile_submit.sh; usage: sbatch [...] genecall_reconcile.sbatch prep|shard|merge
#SBATCH --job-name=recon3b
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
set -euo pipefail
export PYTHONNOUSERSITE=1
PY=/ibex/scratch/projects/c2014/rob/miniconda3/envs/eco_probes/bin/python
D=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/genecall_reconcile_20261102
echo "[$(date -Iseconds)] $(hostname) job ${SLURM_JOB_ID} task ${SLURM_ARRAY_TASK_ID:-NA} mode $* root ${RECON_ROOT:-default} test ${RECON_TEST:-0}"
"$PY" -u "$D/scripts/genecall_reconcile.py" "$@"
echo "[$(date -Iseconds)] done"
