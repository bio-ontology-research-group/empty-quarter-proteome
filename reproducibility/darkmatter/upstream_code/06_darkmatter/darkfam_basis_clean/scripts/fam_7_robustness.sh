#!/bin/bash
# 30%-identity recluster robustness one-liner (review note; Methods sentence).
#SBATCH --job-name=ale_fam_robust
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam7_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=batch
set -euo pipefail
R=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate mmseqs2
mkdir -p "$R/robust_tmp"
mmseqs easy-cluster "$R/dark_reps.faa" "$R/robust30" "$R/robust_tmp" \
    --min-seq-id 0.3 -c 0.8 --cov-mode 0 --threads "$SLURM_CPUS_PER_TASK"
cut -f1 "$R/robust30_cluster.tsv" | sort | uniq -c | awk '$1>=3{n++} END{print "clusters >=3 members at 30% id:", n}'
rm -rf "$R/robust_tmp"
