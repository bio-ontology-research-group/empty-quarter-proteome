#!/bin/bash
# darkclean_1 — derive the UNCOVERED representatives: all_cluster90 clusters with NO member
# carrying a UniRef50 status in the assembly-namespace delta search. Those reps are the ones to
# re-search. Also writes covered/uncovered counts (independent re-confirmation of the 37.8M).
# Assembly namespace only; the Bakta seq_vs_uniref50 layer is intentionally NOT used.
#
# Submit: sbatch 04_protein_catalog/darkmatter_clean/scripts/darkclean_1_uncovered_ids.sh
#
#SBATCH --job-name=velazqal_darkclean1
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean1_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean1_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --time=08:00:00
#SBATCH --partition=batch

set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate instrain

M=/ibex/scratch/projects/c2014/rob/metagenomics
CAT=$M/04_protein_catalog
OUT=$CAT/darkmatter_clean
mkdir -p "$OUT"

python3 - <<'PY'
import os
CAT="/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog"
OUT=f"{CAT}/darkmatter_clean"
CLUSTER=f"{CAT}/all_linclust90/all_cluster90.tsv"
HIT=f"{CAT}/delta_vs_uniref50/hit_ids.txt"
NOHIT=f"{CAT}/delta_vs_uniref50/nohit_ids.txt"

# 1. load the assembly-namespace delta search results (the "searched" id set)
searched = set()
for p in (HIT, NOHIT):
    n0 = len(searched)
    with open(p) as fh:
        for line in fh:
            v = line.strip()
            if v:
                searched.add(v)
    print(f"[info] +{len(searched)-n0:,} from {os.path.basename(p)} (set={len(searched):,})", flush=True)
print(f"[info] searched (delta hit+nohit) = {len(searched):,}", flush=True)

# 2. stream the cluster tsv (grouped by representative); a cluster is COVERED iff any member
#    is in the searched set. Emit the reps of UNCOVERED clusters.
cur = None
covered_flag = False
nrows = nclust = covered = uncovered = 0
u = open(f"{OUT}/01_uncovered_rep_ids.txt", "w")
with open(CLUSTER) as fh:
    for line in fh:
        nrows += 1
        if nrows % 100_000_000 == 0:
            print(f"[info] {nrows:,} rows; {nclust:,} clusters", flush=True)
        t = line.find("\t")
        if t < 0:
            continue
        rep = line[:t]
        mem = line[t+1:].rstrip("\n")
        if rep != cur:
            if cur is not None:
                nclust += 1
                if covered_flag:
                    covered += 1
                else:
                    uncovered += 1
                    u.write(cur + "\n")
            cur = rep
            covered_flag = False
        if mem in searched:
            covered_flag = True
    # final cluster
    if cur is not None:
        nclust += 1
        if covered_flag:
            covered += 1
        else:
            uncovered += 1
            u.write(cur + "\n")
u.close()

with open(f"{OUT}/01_counts.txt", "w") as c:
    c.write(f"clusters\t{nclust}\ncovered\t{covered}\nuncovered\t{uncovered}\n")
print(f"[info] clusters {nclust:,}; covered {covered:,}; uncovered {uncovered:,}", flush=True)
print("DONE_DARKCLEAN1")
PY
