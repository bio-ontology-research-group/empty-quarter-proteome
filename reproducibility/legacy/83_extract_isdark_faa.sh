#!/bin/bash
# Extract the sequence-dark candidate proteins for the "really dark" triple filter
# (sequence-dark AND structure-dark AND InterProScan-negative). Produces:
#   isdark_candidates.faa        — is_dark==1 (no UniRef50 hit OR pident<30)  [563,619]
#   isdark_nohit_ids.txt         — strict subset: NO UniRef50 hit at all      [487,619]
# Stops (*) are stripped (InterProScan/foldseek dislike them). Chunked for the IPS array.
#
# Submit: sbatch /ibex/scratch/projects/c2014/rob/metagenomics/scripts/83_extract_isdark_faa.sh
#
#SBATCH --job-name=hohndor_isdark_faa
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/isdark_faa_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --partition=batch
set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate pyrodigal   # seqkit
R=/ibex/scratch/projects/c2014/rob/metagenomics
ISD=$R/06_darkmatter/isdark/is_dark.tsv
FAA=$R/06_darkmatter/orthogroups/all_rep_proteins.faa
OUT=$R/06_darkmatter/reallydark
mkdir -p "$OUT" "$OUT/chunks"

# id lists
awk -F'\t' 'NR>1 && $3==1{print $1}'            "$ISD" > "$OUT/isdark_ids.txt"
awk -F'\t' 'NR>1 && $3==1 && $2==""{print $1}'  "$ISD" > "$OUT/isdark_nohit_ids.txt"
echo "[extract] is_dark ids: $(wc -l < "$OUT/isdark_ids.txt")   nohit(strict): $(wc -l < "$OUT/isdark_nohit_ids.txt")"

# pull seqs, strip '*', single-line
seqkit grep -f "$OUT/isdark_ids.txt" "$FAA" \
  | seqkit replace -s -p '\*' -r '' \
  | seqkit seq -w 0 > "$OUT/isdark_candidates.faa"
echo "[extract] seqs: $(grep -c '^>' "$OUT/isdark_candidates.faa")"

# split into 200 chunks for the InterProScan array
seqkit split2 -p 200 -O "$OUT/chunks" -f "$OUT/isdark_candidates.faa" 2>&1 | tail -3
ls "$OUT/chunks" | head
echo "[extract] chunks: $(ls "$OUT/chunks"/*.faa 2>/dev/null | wc -l)"
