#!/bin/bash
#SBATCH --job-name=vfdb_census
#SBATCH --account=pi-hohndor
#SBATCH --partition=batch
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_mobilome/vfdb/vfdb_%j.out
# Virulence-factor census: mmseqs2 search of the 5,438 dRep-representative MAG proteomes
# vs VFDB setB (full). Uses the cluster's installed mmseqs2 module (no external binaries).
# Headers MAG-prefixed so each hit maps back to its genome. Stringent thresholds
# (id>=40%, cov>=70%, e<=1e-10) to limit dual-use false positives; finer id tiers downstream.
# NB: VF hits in soil bacteria are largely competition/symbiosis genes, not human pathogenesis.
set -euo pipefail
M=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$M/06_mobilome/vfdb
mkdir -p "$OUT"
cd "$OUT"

module load mmseqs2/14.7 2>/dev/null || module load mmseqs2/14.7e284 2>/dev/null || module load mmseqs2/13.4
echo "mmseqs: $(mmseqs version 2>&1 | head -1)"

# VFDB setB full protein set (reference data download)
if [ ! -s VFDB_setB_pro.fas ]; then
  wget -q "http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz" -O VFDB_setB_pro.fas.gz
  gunzip -kf VFDB_setB_pro.fas.gz
fi
echo "VFDB proteins: $(grep -c '^>' VFDB_setB_pro.fas)"

# query = all dRep rep proteomes, headers prefixed with MAG (filename)
REP=$M/06_darkmatter/drep_rep_proteomes
if [ ! -s allreps.faa ]; then
  for f in "$REP"/*.faa.gz; do
    mag=$(basename "$f" .faa.gz)
    zcat "$f" | sed "s/^>/>${mag}|/"
  done > allreps.faa
fi
echo "query proteins: $(grep -c '^>' allreps.faa) from $(ls "$REP"/*.faa.gz | wc -l) MAGs"

rm -rf "$OUT/tmp"
mmseqs easy-search allreps.faa VFDB_setB_pro.fas vfdb_hits.tsv "$OUT/tmp" \
  --min-seq-id 0.4 -c 0.7 --cov-mode 0 -e 1e-10 --threads 32 --max-seqs 3 \
  --format-output "query,target,pident,alnlen,evalue,bits,theader"
echo "raw hits: $(wc -l < vfdb_hits.tsv)"
touch VFDB_DONE
echo "DONE -> $OUT/vfdb_hits.tsv"
