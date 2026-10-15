#!/bin/bash
# §3.1 SMAG leg, step a: fetch the SMAG genomes (Ma et al. 2023, Nat Commun 14:7318; Zenodo 8223844,
# CC-BY-4.0). 40,039 soil MAGs (>=50% complete, <10% contaminated), NUCLEOTIDE ONLY, shipped as one
# gzip tarball split into 373 parts mag.tar.gz.aa .. .oi (39.0 GB). The record's own README says
# `tar -xj`; the first bytes are 1f8b, so it is gzip and we use `tar xz`.
#
# What this does, idempotently:
#   1. file manifest (key, size, md5, url) from the Zenodo API, written to tmp/smag/zenodo_files.tsv
#   2. every part downloaded with retries and md5-verified; a part already present with the right md5
#      is skipped, a wrong one is re-fetched
#   3. Supplementary Data 2.xlsx (per-MAG GTDB-Tk taxonomy + CheckM2, one row per MAG) -> results/
#      (Strict OOXML: openpyxl reads it as zero sheets; kept as provenance, not parsed here)
#   4. `cat parts | tar xz` streamed into tmp/smag/genomes/ (no concatenated archive is ever staged),
#      then the genome list = files named <SD2 user_genome>.fa, asserted to be exactly 40,039; the
#      parts are deleted only after both tar and that assertion pass
#   sbatch spec_6a_smag_fetch.sbatch
#SBATCH --job-name=ale_spec_smag_fetch
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_smag_fetch_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean/logs/spec_smag_fetch_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
#SBATCH --account=pi-hohndor
set -uo pipefail
S=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/specificity_clean
T=$S/tmp/smag
PARTS=$T/parts
G=$T/genomes
PY=/home/velazqal/miniconda3/envs/python-utils/bin/python
API=https://zenodo.org/api/records/8223844
N_PARTS_EXPECT=373
N_MAGS_EXPECT=40039
mkdir -p "$PARTS" "$G" "$S/results"
log(){ echo "[$(date -Iseconds)] $*"; }

A=$(df -B1G --output=avail "$S" | tail -1 | tr -d ' ')
log "free on $S: ${A} GB (peak need ~200 GB)"
[ "$A" -lt 400 ] && { log "ABORT: <400 GB free"; exit 1; }

if [ -f "$T/EXTRACTED.ok" ]; then
  log "genomes already extracted ($(wc -l < "$T/genome_list.txt") listed) — nothing to do"; exit 0
fi

# ---------------- 1. manifest ----------------
curl -sfL --retry 5 "$API" -o "$T/record.json" || { log "ABORT: cannot reach Zenodo from $(hostname)"; exit 1; }
"$PY" - "$T/record.json" "$T/zenodo_files.tsv" <<'EOF'
import json, sys
d = json.load(open(sys.argv[1]))
with open(sys.argv[2], "w") as o:
    for f in sorted(d["files"], key=lambda f: f["key"]):
        o.write("%s\t%d\t%s\t%s\n" % (f["key"], f["size"], f["checksum"].split(":", 1)[1], f["links"]["self"]))
EOF
NP=$(awk -F'\t' '$1 ~ /^mag\.tar\.gz\./' "$T/zenodo_files.tsv" | wc -l)
log "parts in record: $NP (expect $N_PARTS_EXPECT)"
[ "$NP" -ne "$N_PARTS_EXPECT" ] && { log "ABORT: part count differs"; exit 1; }

# ---------------- 2. parts, md5-verified ----------------
fetch(){ # key md5 url dir
  local key=$1 md5=$2 url=$3 dir=$4 out="$4/$1" try
  if [ -f "$out" ] && [ "$(md5sum "$out" | cut -d' ' -f1)" = "$md5" ]; then return 0; fi
  for try in 1 2 3 4 5; do
    rm -f "$out.part"
    curl -sfL --retry 10 --retry-all-errors --retry-delay 20 -o "$out.part" "$url"
    if [ -f "$out.part" ] && [ "$(md5sum "$out.part" | cut -d' ' -f1)" = "$md5" ]; then
      mv "$out.part" "$out"; return 0
    fi
    echo "[$key] md5 mismatch or failed download, attempt $try" >&2; sleep 30
  done
  echo "[$key] FAILED after 5 attempts" >&2; return 1
}
export -f fetch
log "downloading parts (4 in parallel)"
awk -F'\t' '$1 ~ /^mag\.tar\.gz\./ {print $1"\t"$3"\t"$4}' "$T/zenodo_files.tsv" \
  | xargs -P 4 -L 1 bash -c 'fetch "$0" "$1" "$2" '"$PARTS"
NOK=0
while IFS=$'\t' read -r key size md5 url; do
  [[ $key == mag.tar.gz.* ]] || continue
  [ -f "$PARTS/$key" ] && [ "$(md5sum "$PARTS/$key" | cut -d' ' -f1)" = "$md5" ] && NOK=$((NOK+1))
done < "$T/zenodo_files.tsv"
log "parts md5-verified: $NOK / $N_PARTS_EXPECT"
[ "$NOK" -ne "$N_PARTS_EXPECT" ] && { log "ABORT: not every part verified — rerun to resume"; exit 1; }

# ---------------- 3. per-MAG taxonomy table ----------------
IFS=$'\t' read -r _ _ SMD5 SURL < <(awk -F'\t' '$1=="Supplementary Data 2.xlsx"' "$T/zenodo_files.tsv")
fetch "Supplementary Data 2.xlsx" "$SMD5" "$SURL" "$T" || { log "ABORT: SD2 download failed"; exit 1; }
cp "$T/Supplementary Data 2.xlsx" "$S/results/smag_supplementary_data_2.xlsx"
log "SD2 (per-MAG GTDB taxonomy) -> results/smag_supplementary_data_2.xlsx md5=$SMD5"

# ---------------- 4. stream-extract, then drop the parts ----------------
if [ -f "$T/TAR.ok" ]; then
  log "tar already completed — skipping extraction"
else
  log "cat parts | tar xz -> $G"
  find "$PARTS" -maxdepth 1 -name 'mag.tar.gz.*' ! -name '*.part' -print0 | LC_ALL=C sort -z \
    | xargs -0 cat | tar xz -C "$G"
  RC=("${PIPESTATUS[@]}")
  log "pipeline rc: ${RC[*]}"
  for r in "${RC[@]}"; do [ "$r" -ne 0 ] && { log "ABORT: extraction failed — parts kept"; exit 1; }; done
  touch "$T/TAR.ok"
fi

# The catalogue is DEFINED by SD2's user_genome column (40,039 MAGs). The tarball holds more files
# than that (checked 2026-09-24: 40,352 in MAG/ = 40,039 catalogue MAGs + 312 timestamped duplicate
# copies such as "BL1bin.10 (2022:5:31 ...).fa", byte-identical to their originals, + a few bins not in
# SD2 + stray k141/nohup.out). Only files whose name is <user_genome>.fa are gene-called.
"$PY" - "$T/Supplementary Data 2.xlsx" "$T/sd2_user_genome.txt" <<'EOF'
import zipfile, re, sys          # Strict OOXML: openpyxl sees no sheets, so read the XML directly
z = zipfile.ZipFile(sys.argv[1])
ss = [''.join(re.findall(r'<t[^>]*>([^<]*)</t>', s))
      for s in re.findall(r'<si>(.*?)</si>', z.read('xl/sharedStrings.xml').decode(), re.S)]
ids = [ss[int(v)] if 't="s"' in a else v for a, v in
       re.findall(r'<c r="B\d+"([^>]*)><v>([^<]*)</v>', z.read('xl/worksheets/sheet1.xml').decode())]
ids = ids[ids.index('user_genome') + 1:]
open(sys.argv[2], 'w').write('\n'.join(ids) + '\n')
EOF
NS=$(sort -u "$T/sd2_user_genome.txt" | wc -l)
NF=$(find "$G" -type f | wc -l)
awk 'NR==FNR {want[$1".fa"]=1; next} { f=$0; sub(/.*\//, "", f); if (f in want) print }' \
  "$T/sd2_user_genome.txt" <(find "$G" -type f) | LC_ALL=C sort > "$T/genome_list.txt"
NG=$(wc -l < "$T/genome_list.txt")
log "SD2 MAG ids: $NS; files extracted: $NF; catalogue genomes matched: $NG (expect $N_MAGS_EXPECT); ignored $((NF - NG)) files; extracted size $(du -sh "$G" | cut -f1)"
[ "$NS" -ne "$N_MAGS_EXPECT" ] || [ "$NG" -ne "$N_MAGS_EXPECT" ] && { log "ABORT: catalogue not fully matched — inspect $G, parts kept"; exit 1; }

rm -rf "$PARTS"
touch "$T/EXTRACTED.ok"
log "parts deleted; DONE_SMAG_FETCH"
