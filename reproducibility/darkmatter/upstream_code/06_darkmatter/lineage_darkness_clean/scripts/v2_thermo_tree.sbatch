#!/bin/bash
# v2_thermo_tree - approximate-ML tree of the desert Thermoplasmatota against GTDB references.
#   step 1: python v2_thermo_tree.py      (builds results/thermo_tree.fasta, seconds, login ok)
#   step 2: sbatch scripts/v2_thermo_tree.sbatch
#SBATCH --job-name=ale_thermo_tree
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean/logs/thermo_tree_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean/logs/thermo_tree_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=batch
set -euo pipefail

P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean
FASTTREE=/home/velazqal/miniconda3/envs/phabox2/bin/FastTreeMP
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-16}

test -s "$P/results/thermo_tree.fasta" || { echo "alignment missing; run v2_thermo_tree.py first" >&2; exit 1; }

echo "taxa: $(grep -c '^>' "$P/results/thermo_tree.fasta")"
echo "threads: $OMP_NUM_THREADS"
"$FASTTREE" -lg -gamma "$P/results/thermo_tree.fasta" > "$P/results/thermo_tree.nwk"
echo "done: $(wc -c < "$P/results/thermo_tree.nwk") bytes"
