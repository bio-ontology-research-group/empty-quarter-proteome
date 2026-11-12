#!/bin/bash
#SBATCH --account=k10018
#SBATCH --partition=shared
#SBATCH --job-name=eq_a4_sha
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/project/k10018/eq-protein-external-20261110/analysis4_structures/logs/%x-%j.out
# SHA-256 of the inputs and of the alignment outputs, for the receipt.
set -euo pipefail
ROOT=/scratch/project/k10018/eq-protein-external-20261110
A4=$ROOT/analysis4_structures
cd $ROOT
{
ls genecall/harmonized_proteins.faa genecall/harmonized_gene_manifest.tsv.gz metadata/harmonized_genome_manifest.tsv \
   mmseqs/harmonized_uniref50_per_gene.tsv.gz signalp/harmonized_signalp6_predictions.tsv.gz \
   analysis2_families/clusters/families_id30_c80.tsv \
   analysis2_families/results/dark_exported_widespread_families_id30_c80.tsv \
   source-metadata/ani_groups/ani_95_groups_1957.tsv analysis2_families/mmseqs_db/all \
   analysis4_structures/work/step1/*.tsv analysis4_structures/work/step1/queries.faa \
   analysis4_structures/results/alignments/* analysis4_structures/results/queries_pfam.domtbl \
   analysis4_structures/scripts/* | xargs -P 8 -n 1 sha256sum
} | sort -k2 > $A4/results/sha256_inputs_outputs.txt
cat $A4/results/sha256_inputs_outputs.txt
echo SHA_DONE
