#!/bin/bash
# darkclean_4 — combine the delta + uncovered UniRef50 searches into the COMPLETE dark set over all
# 449.3M clusters, in one clean pass (assembly namespace only; no Bakta). A cluster is DARK iff it has
# a searched member and NO hit member. Recomputes from scratch (does not depend on catalog_ladder).
#
# Submit (after darkclean_3c): sbatch .../darkmatter_clean/scripts/darkclean_4_combine.sh
#
#SBATCH --job-name=velazqal_darkclean4
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean4_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkclean4_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G
#SBATCH --time=08:00:00
#SBATCH --partition=batch

set -euo pipefail
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate instrain

python3 - <<'PY'
import os, gzip
CAT="/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog"
OUT=f"{CAT}/darkmatter_clean"
CLUSTER=f"{CAT}/all_linclust90/all_cluster90.tsv"
HIT=[f"{CAT}/delta_vs_uniref50/hit_ids.txt",   f"{OUT}/03_uncovered_vs_uniref50/hit_ids.txt"]
NOHIT=[f"{CAT}/delta_vs_uniref50/nohit_ids.txt", f"{OUT}/03_uncovered_vs_uniref50/nohit_ids.txt"]
os.makedirs(f"{OUT}/04_catalog_dark", exist_ok=True)

def load(paths):
    s=set()
    for p in paths:
        n0=len(s)
        with open(p) as fh:
            for line in fh:
                v=line.strip()
                if v: s.add(v)
        print(f"[info] +{len(s)-n0:,} from {p} (set={len(s):,})", flush=True)
    return s

hits=load(HIT); nohits=load(NOHIT)
print(f"[info] hits={len(hits):,} nohits={len(nohits):,}", flush=True)

cur=None; cur_hit=False; cur_srch=False
nrows=nclust=hitc=darkc=uncov=0
dark=gzip.open(f"{OUT}/04_catalog_dark/dark_rep_ids.txt.gz","wt")
with open(CLUSTER) as fh:
    for line in fh:
        nrows+=1
        if nrows%100_000_000==0: print(f"[info] {nrows:,} rows; {nclust:,} clusters", flush=True)
        t=line.find("\t")
        if t<0: continue
        rep=line[:t]; mem=line[t+1:].rstrip("\n")
        if rep!=cur:
            if cur is not None:
                nclust+=1
                if cur_hit: hitc+=1
                elif cur_srch: darkc+=1; dark.write(cur+"\n")
                else: uncov+=1
            cur=rep; cur_hit=False; cur_srch=False
        if mem in hits: cur_hit=True; cur_srch=True
        elif mem in nohits: cur_srch=True
    if cur is not None:
        nclust+=1
        if cur_hit: hitc+=1
        elif cur_srch: darkc+=1; dark.write(cur+"\n")
        else: uncov+=1
dark.close()

covered=hitc+darkc
with open(f"{OUT}/04_catalog_dark/dark_summary.tsv","w") as s:
    s.write("metric\tvalue\tfraction\n")
    s.write(f"clusters\t{nclust}\t1.0\n")
    s.write(f"hit\t{hitc}\t{hitc/nclust:.6f}\n")
    s.write(f"dark\t{darkc}\t{darkc/nclust:.6f}\n")
    s.write(f"still_uncovered\t{uncov}\t{uncov/nclust:.6f}\n")
    s.write(f"dark_of_covered\t{darkc}\t{darkc/covered:.6f}\n")
print(f"[info] clusters {nclust:,}; hit {hitc:,}; dark {darkc:,}; still_uncovered {uncov:,}", flush=True)
print("DONE_DARKCLEAN4")
PY
