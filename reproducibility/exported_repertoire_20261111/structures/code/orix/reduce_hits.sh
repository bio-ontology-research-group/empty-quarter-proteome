#!/bin/bash
# Login-node text handling (seconds): keep the five best hits per query and the number of hits per query.
set -euo pipefail
cd /mnt/data/u/hohndor/eq-exported-dark-20261111/foldseek/hits; mkdir -p ../top
for a in whole crops; do
  LC_ALL=C sort -t$'\t' -k1,1 -k4,4gr pdb_$a.tsv | awk -F'\t' '$1!=q{q=$1;n=0} ++n<=5' > ../top/pdb_$a.tsv
  cut -f1 pdb_$a.tsv | LC_ALL=C sort | uniq -c | awk '{print $2"\t"$1}' > ../top/pdb_$a.nhits.tsv
  LC_ALL=C sort -t$'\t' -k1,1 -k3,3g afdb_$a.tsv | awk -F'\t' '$1!=q{q=$1;n=0} ++n<=5' > ../top/afdb_$a.tsv
  cut -f1 afdb_$a.tsv | LC_ALL=C sort | uniq -c | awk '{print $2"\t"$1}' > ../top/afdb_$a.nhits.tsv
done
wc -l ../top/*
