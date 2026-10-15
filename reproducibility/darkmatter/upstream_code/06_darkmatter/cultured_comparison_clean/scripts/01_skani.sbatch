#!/bin/bash
#SBATCH --job-name=poc05_skani
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/skani_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean/logs/skani_%j.err
# RAK-PR cultured genomes vs the 5,438 soil MAG species representatives.
# skani 0.3.1 from the IBEX GTDB-Tk 2.5.2 conda env, as used by
# W/07_aquatic_paper/07_endemicity (sketch_gtdb.sh, search_gtdb.sh).
# Deviation from that analysis: no --qi. In skani 0.3.1 --qi means "use individual sequences
# of a multi-fasta as separate queries", i.e. per-contig ANI. Here we need whole-genome ANI
# and a genome-level aligned fraction, so queries are passed as whole files via --ql.
set -uo pipefail
SK=/ibex/sw/rl9c/gtdb-tk/2.5.2/rl9.4_conda3/Miniconda3/envs/gtdbtk_v2.5.2/bin/skani
P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean
DB=$P/skani_mag_db
T=${SLURM_CPUS_PER_TASK:-32}

echo "skani: $($SK --version 2>&1 | head -1)"
echo "queries (cultured): $(wc -l < $P/scripts/cultured_paths.txt)"
echo "refs (MAG reps):    $(wc -l < $P/scripts/mag_paths.txt)"

echo "=== STAGE 1: sketch the 5,438 MAG species representatives ==="; date
rm -rf $DB
$SK sketch -l $P/scripts/mag_paths.txt -o $DB -t $T
echo "sketch exit $?; sketched files: $(ls $DB/*.sketch 2>/dev/null | wc -l)"; du -sh $DB

echo "=== STAGE 2: search each cultured genome against the sketch DB (-s 75) ==="; date
$SK search -d $DB --ql $P/scripts/cultured_paths.txt \
   -o $P/results/skani_search.tsv -t $T -s 75
echo "search exit $?; rows: $(wc -l < $P/results/skani_search.tsv)"

echo "=== STAGE 3: all-vs-all dist, 722 x 5,438, as a cross-check and for the reverse direction ==="; date
$SK dist --ql $P/scripts/cultured_paths.txt --rl $P/scripts/mag_paths.txt \
   -o $P/results/skani_dist.tsv -t $T -s 75
echo "dist exit $?; rows: $(wc -l < $P/results/skani_dist.tsv)"

echo "=== DONE ==="; date
