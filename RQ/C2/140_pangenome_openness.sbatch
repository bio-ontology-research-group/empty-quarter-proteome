#!/bin/bash
# RQ-C2: historical exploratory pangenome workflow.
#   Reproduces the archived raw-union/dRep analysis, not the current checksum-closed
#   catalogue. Core/accessory vs UniRef50-no-hit is an association, not a mechanism.
#   Refs: Tettelin 2005 (PNAS); McInerney et al. 2017 (Nat Microbiol).
#
#   Clustering (the real compute) runs here on a batch node -- NEVER a login node.
#   Submit:  sbatch /ibex/scratch/projects/c2014/rob/metagenomics/scripts/ibex/140_pangenome_openness.sbatch
#SBATCH --job-name=hohndor_pangenome
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/pangenome_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/pangenome_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
set -uo pipefail
export M=/ibex/scratch/projects/c2014/rob/metagenomics
mkdir -p "$M/logs" "$M/06_pangenome"
source "$M/../miniconda3/etc/profile.d/conda.sh"; conda activate eggnog
echo "[$(date -Iseconds)] host=$(hostname) cpus=${SLURM_CPUS_PER_TASK} mmseqs=$(mmseqs version 2>/dev/null)"
export OUTDIR="$M/06_pangenome"
export WORK="$M/06_pangenome/work"
export MIN_GENOMES=10 CORE_FRAC=0.95 SOFT_FRAC=0.90 SHELL_FRAC=0.15 PERMS=50
export MIN_SEQ_ID=0.5 COV=0.8
python "$M/scripts/ibex/140_pangenome_openness.py"
echo "[$(date -Iseconds)] pangenome job finished rc=$?"
