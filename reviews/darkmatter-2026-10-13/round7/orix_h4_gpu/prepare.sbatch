#!/bin/bash
#SBATCH --job-name=h4-orix-prepare
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/h4_structure/orix_gpu_20261014/prepare_%j.log
set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean
exec "$B/require_slurm_worker.sh" /ibex/scratch/projects/c2014/rob/miniconda3/envs/metabolism-env/bin/python "$B/h4_structure/orix_gpu_20261014/prepare.py"
