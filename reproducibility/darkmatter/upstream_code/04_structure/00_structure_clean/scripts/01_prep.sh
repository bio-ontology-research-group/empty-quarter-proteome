#!/bin/bash
# 01_prep — split the 4.29M genome-resolved dark reps into:
#   reuse_ids   : already have a ProstT5 3Di in tier1/tier2  (~2.9M) -> reuse
#   missing_ids : no 3Di yet                                  (~1.4M) -> 02_prostt5 generates them
# and extract the missing reps' sequences, chunked for the ProstT5 GPU array.
#
#SBATCH --job-name=sc_prep
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc_prep_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/sc_prep_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=batch

set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
S=$M/04_structure
OUT=$S/00_structure_clean/01_prep; mkdir -p "$OUT/missing_chunks"
DARK=$M/04_protein_catalog/darkmatter_clean/05_mag_dark/mag_dark_rep_ids.txt.gz
REPS=$M/04_protein_catalog/all_linclust90/all_reps90.faa.gz

# 1. dark-MAG rep IDs (4.29M) and all 3Di IDs (24.5M = tier1 ∪ tier2)
zcat "$DARK" | sort -u > "$OUT/dark_ids.txt"
cat "$S/tier1_cluster/concat/"*.lookup "$S/tier2_cluster/concat/"*.lookup | cut -f2 | sort -u > "$OUT/threedi_ids.txt"

# 2. split reuse vs missing
comm -12 "$OUT/dark_ids.txt" "$OUT/threedi_ids.txt" > "$OUT/reuse_ids.txt"
comm -23 "$OUT/dark_ids.txt" "$OUT/threedi_ids.txt" > "$OUT/missing_ids.txt"
echo "[$(date -Iseconds)] reuse=$(wc -l < "$OUT/reuse_ids.txt")  missing=$(wc -l < "$OUT/missing_ids.txt")"

# 3. extract the missing reps' sequences from the catalogue rep FASTA (streaming; header 1st token = protein id)
python3 - "$OUT/missing_ids.txt" "$REPS" "$OUT/missing.faa" <<'PY'
import sys, gzip
ids = {l.strip() for l in open(sys.argv[1])}
op = gzip.open if sys.argv[2].endswith(".gz") else open
keep=False; n=0
with op(sys.argv[2],"rt") as fh, open(sys.argv[3],"w") as out:
    for line in fh:
        if line[0]==">":
            keep = line[1:].split()[0] in ids
            if keep: n+=1
        if keep: out.write(line)
print(f"[info] wrote {n} / {len(ids)} missing sequences", flush=True)
PY

# 4. chunk into ~8 files (~175k proteins each) for the A100 array (same shape as 06c_prostt5_prep)
NC=8
N=$(grep -c '^>' "$OUT/missing.faa"); PER=$(( (N + NC - 1) / NC ))
echo "[$(date -Iseconds)] chunking $N missing seqs into $NC files (~$PER each)"
awk -v per="$PER" -v od="$OUT/missing_chunks" -v nc="$NC" '
  /^>/ { c++; nch=int((c-1)/per)+1; if(nch>nc)nch=nc; if(nch!=ch){ if(of)close(of); ch=nch; of=sprintf("%s/chunk_%03d.fa",od,ch) } }
  { print > of }' "$OUT/missing.faa"
ls -l "$OUT/missing_chunks/"
echo "DONE_PREP"
