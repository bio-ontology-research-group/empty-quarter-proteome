#!/bin/bash
# Phase D (heavy) .1 — download the MGnify90 protein cluster catalogue (the largest global
# environmental protein set, AA-native, biome-tagged) + build a DIAMOND DB, and extract the 23,382
# dark family rep proteins as the query. Tests whether EQ desert dark matter is globally shared or
# desert-endemic (these reps are already UniRef50-dark; MGnify covers env sequences absent from UniRef).
#   sbatch scripts/ibex/98a_mgnify_download_makedb.sh
#SBATCH --job-name=hohndor_mgnify_dl
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/mgnify_dl_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/mgnify_dl_%j.err
#SBATCH --cpus-per-task=24
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
DB=$WD/databases/mgnify90
DM=$WD/06_darkmatter
OUT=$DM/darkfam_ecology/endemism
mkdir -p "$DB" "$OUT"
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eggnog
T=${SLURM_CPUS_PER_TASK:-24}
FTP=https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release

echo "[$(date -Iseconds)] download mgy_clusters (74G) + biomes (resumable)"
cd "$DB"
wget -c -q "$FTP/mgy_clusters.fa.gz" -O mgy_clusters.fa.gz
wget -c -q "$FTP/mgy_biomes.tsv.gz"  -O mgy_biomes.tsv.gz
sz=$(stat -c %s mgy_clusters.fa.gz)
echo "  mgy_clusters.fa.gz size=$sz (expect 79301981114)"
[ "$sz" -ge 79000000000 ] || { echo "ERROR: incomplete download"; exit 1; }

echo "[$(date -Iseconds)] diamond makedb (~718M proteins)"
diamond makedb --in mgy_clusters.fa.gz -d mgy_clusters -p "$T" 2>&1 | tail -3
ls -lh mgy_clusters.dmnd | awk '{print "  db:",$5}'

echo "[$(date -Iseconds)] extract 23,382 dark rep proteins as query"
python3 - "$DM/orthogroups/orthogroups_rep_seq.fasta" "$DM/darkfam_ecology/catalog/darkfam_meta.tsv" "$OUT/query.faa" <<'PY'
import sys
rep_fa, meta, outp = sys.argv[1:4]
want = set()
for i, l in enumerate(open(meta)):
    if i == 0: continue
    want.add(l.split("\t")[0])
o = open(outp, "w"); keep = False; n = 0
for line in open(rep_fa):
    if line.startswith(">"):
        keep = line[1:].split()[0] in want
        if keep: o.write(line); n += 1
    elif keep:
        o.write(line)
o.close()
print(f"  query proteins extracted: {n} / {len(want)}", flush=True)
PY
grep -c '^>' "$OUT/query.faa"
date -Iseconds > "$DB/MGNIFY_DB_DONE"
echo "[$(date -Iseconds)] download+makedb+query done"
