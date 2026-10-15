#!/bin/bash
#SBATCH --job-name=h4-gpu-root-audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/orix_gpu_20261014/audit_%j.log
set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean
exec "$B/require_slurm_worker.sh" /ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python -u "$B/h4_structure/orix_gpu_20261014/audit_h4_gpu.py"
