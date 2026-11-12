#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=shared
#SBATCH --job-name=eq_a4_clean
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis4_structures/logs/%x-%j.out
# Inode hygiene: pack the 1,709 loose colabfold_search a3m files into one tar, drop search intermediates,
# and refresh the checksum list.
set -euo pipefail
A4=/scratch/project/k10018/eq-protein-external-20261110/analysis4_structures
cd $A4/work/cfsearch
n=$(ls *.a3m | wc -l); [ "$n" -eq 1709 ]
tar -cf $A4/results/alignments/colabfold_search_uniref30_2302_envdb_202108_a3m.tar *.a3m
[ "$(tar -tf $A4/results/alignments/colabfold_search_uniref30_2302_envdb_202108_a3m.tar | wc -l)" -eq 1709 ]
cd $A4; rm -rf work/cfsearch; rm -f work/search/res*
find $A4 -type f | wc -l
cd /scratch/project/k10018/eq-protein-external-20261110
ls analysis4_structures/results/alignments/* analysis4_structures/scripts/* | xargs -P 4 -n 1 sha256sum | sort -k2 > $A4/results/sha256_alignments_scripts_final.txt
cat $A4/results/sha256_alignments_scripts_final.txt
echo CLEAN_DONE
