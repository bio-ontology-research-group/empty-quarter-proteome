#!/bin/bash
#SBATCH --job-name=dm7-review-audit
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/manuscript_review_20261014/verify-%j.log
set -euo pipefail
cd /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p7_analyses_clean/manuscript_review_20261014
bash reproducibility/darkmatter/round7_analyses/require_slurm_worker.sh sha256sum -c reviews/darkmatter-2026-10-13/round7/manuscript/ibex_mirror.sha256
