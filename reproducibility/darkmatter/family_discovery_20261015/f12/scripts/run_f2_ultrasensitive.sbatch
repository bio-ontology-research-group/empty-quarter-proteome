#!/bin/bash
#SBATCH --job-name=f2_ultrasensitive
#SBATCH --account=pi-hohndor
#SBATCH --partition=debug,batch
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/f2_ultra_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12/logs/f2_ultra_%j.err

set -euo pipefail

WORKDIR="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12"
cd "${WORKDIR}"

DIAMOND_DIR="/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin"
export PATH="${DIAMOND_DIR}:${PATH}"

./require_slurm_worker.sh diamond version

echo "Running DIAMOND blastp --ultra-sensitive on f2_leads.faa..."
diamond blastp \
    --threads 8 \
    --query results/pilot_v3_results/f2_leads.faa \
    --db references_v2/phosphatase_references_dedup.dmnd \
    --out results/pilot_v3_results/f2_leads_diamond_ultrasensitive.tsv \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
    --evalue 1e-3 \
    --max-target-seqs 0 \
    --ultra-sensitive

echo "DIAMOND --ultra-sensitive finished."
echo "Number of hits:"
wc -l results/pilot_v3_results/f2_leads_diamond_ultrasensitive.tsv || true

