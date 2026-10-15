#!/bin/bash
# v2_thermo_relatives_gtdbtk - GTDB-Tk 2.7.2 classify_wf (release 232) on the MGnify comparators only,
# to confirm their class at r232 (MGnify lineages come from each catalogue's own GTDB release).
# Invocation copied from rebuild_20260729/coassembly_novelty/legacy_taxonomic_screen_r232 (classify_wf,
# ANI screen on); pplacer on one CPU to keep the ar53 placement memory bounded.
#   after v2_thermo_relatives_fetch.py:  sbatch scripts/v2_thermo_relatives_gtdbtk.sbatch
#SBATCH --job-name=ale_thermo_rel_gtdbtk
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean/logs/thermo_relatives_gtdbtk_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean/logs/thermo_relatives_gtdbtk_%j.err
set -euo pipefail

P=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/lineage_darkness_clean
T=$P/tmp/thermo_relatives
module load gtdb-tk/2.7.2
grep -q "VERSION_DATA=r232" $GTDBTK_DATA_PATH/metadata/metadata.txt \
    || { echo "GTDB-Tk data at $GTDBTK_DATA_PATH is not release 232" >&2; exit 1; }
echo "GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH ($(grep VERSION_DATA $GTDBTK_DATA_PATH/metadata/metadata.txt))"
echo "genomes: $(ls $T/mgnify/*.fna | wc -l)"; date

rm -rf $T/gtdbtk $T/gtdbtk_scratch; mkdir -p $T/gtdbtk_scratch
gtdbtk classify_wf --genome_dir $T/mgnify --extension fna --out_dir $T/gtdbtk \
    --cpus ${SLURM_CPUS_PER_TASK:-32} --pplacer_cpus 1 --scratch_dir $T/gtdbtk_scratch

S=$T/gtdbtk/gtdbtk.ar53.summary.tsv
test -s $S
cp $S $P/results/v2_thermo_relatives_gtdbtk_mgnify.tsv
test -s $T/gtdbtk/gtdbtk.bac120.summary.tsv && echo "WARNING: some MGnify genomes classified as bacteria"
echo "archaeal rows: $(($(wc -l < $S) - 1))"; date
