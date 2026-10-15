#!/bin/bash
# Phase D (heavy) .3 — biome breakdown of the globally-shared EQ dark families: which global habitats do
# their MGnify90 matches come from? Joins hit MGYPs -> mgy_biomes.tsv.gz (MGYP -> biome lineage).
#   sbatch scripts/ibex/98c_mgnify_biomes.sh
#SBATCH --job-name=hohndor_mgnify_biomes
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/mgnify_biomes_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/mgnify_biomes_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
DB=$WD/databases/mgnify90
O=$WD/06_darkmatter/darkfam_ecology/endemism

python3 - "$O/mgnify_hits.m8" "$DB/mgy_biomes.tsv.gz" "$O" <<'PY'
import sys, gzip
from collections import defaultdict, Counter
hits, biomes, O = sys.argv[1], sys.argv[2], sys.argv[3]

# MGYP -> set(EQ reps that hit it)
mgyp_reps = defaultdict(set)
for line in open(hits):
    c = line.split("\t")
    mgyp_reps[c[1]].add(c[0])
print(f"[biomes] unique hit MGYPs: {len(mgyp_reps)}", flush=True)

# stream biomes; for matched MGYP, map its biome (level-2/3 category) to the EQ reps
def cat(lineage):
    p = lineage.strip().split(":")               # root:A:B:C...
    return ":".join(p[1:3]) if len(p) >= 3 else (p[1] if len(p) >= 2 else "root")
rep_biomes = defaultdict(set)                     # EQ rep -> set(biome categories)
n = 0
for line in gzip.open(biomes, "rt"):
    c = line.rstrip("\n").split("\t")
    if len(c) < 3: continue
    reps = mgyp_reps.get(c[0])
    if not reps: continue
    b = cat(c[2])
    for r in reps: rep_biomes[r].add(b)
    n += 1
print(f"[biomes] biome rows matched: {n} ; shared reps with biome: {len(rep_biomes)}", flush=True)

# per-biome: how many EQ shared families occur there
bc = Counter()
for r, bs in rep_biomes.items():
    for b in bs: bc[b] += 1
nshared = len(rep_biomes)
with open(f"{O}/BIOME_BREAKDOWN.txt", "w") as o:
    def w(x): o.write(x + "\n"); print(x)
    w("# Phase D — biomes of the globally-shared EQ dark families (MGnify90 matches)")
    w(f"shared EQ dark families with a biome-annotated hit: {nshared}")
    w("\n# EQ shared families occurring in each global biome category (a family can span several):")
    for b, k in bc.most_common(30):
        w(f"  {k:6d} ({100*k/nshared:.1f}%)  {b}")
    multi = sum(1 for bs in rep_biomes.values() if len(bs) >= 3)
    w(f"\ncosmopolitan (>=3 distinct biomes): {multi} ({100*multi/nshared:.1f}%)")
PY
echo "[$(date -Iseconds)] biome breakdown done -> $O/BIOME_BREAKDOWN.txt"
