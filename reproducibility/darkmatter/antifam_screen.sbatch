#!/bin/bash
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=03:00:00
set -euo pipefail
BASE=/ibex/scratch/projects/c2014/rob/metagenomics
ROOT=$BASE/06_darkmatter/p4_review_20261013
OUT=$ROOT/results/antifam
mkdir -p "$OUT"
HMM=$BASE/../../rob/miniconda3/envs/eggnog/bin/hmmsearch
HMM=/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/hmmsearch
"$HMM" -h > "$OUT/hmmer_version.txt"
sha256sum "$HMM" "$ROOT/scripts/AntiFam.hmm" > "$OUT/program_database.sha256"
for ARM in dark ctrl; do
    FAA=$BASE/06_darkmatter/darkfam_basis_clean/results/${ARM}_reps.faa
    sha256sum "$FAA" >> "$OUT/input.sha256"
    "$HMM" --cpu 16 --cut_ga --noali --tblout "$OUT/${ARM}.tbl" --domtblout "$OUT/${ARM}.domtbl" "$ROOT/scripts/AntiFam.hmm" "$FAA" > "$OUT/${ARM}.log"
done
