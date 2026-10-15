#!/bin/bash
#SBATCH --job-name=ale_pfam_cand
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/pfam_candidates_clean/logs/02_hmmsearch_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/pfam_candidates_clean/logs/02_hmmsearch_%j.err
# Thin wrapper over the plain-bash worker.
set -euo pipefail
bash /ibex/scratch/projects/c2014/rob/metagenomics/04_structure/pfam_candidates_clean/scripts/02_hmmsearch_pfam.sh
