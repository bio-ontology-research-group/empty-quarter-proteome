#!/bin/bash
# REBUILD step 5 (REWRITTEN 2026-08-13 after multi-agent review) — length-matched ANNOTATED control
# for the dN/dS comparison (§1.4 leg 2). Review fixes incorporated:
#  (F1) MAG-mapped gate now replicates darkclean_5_mag.sh EXACTLY: contigs count only if their bin
#       is in the 37,949 HQ+MQ set from results/catalogue_provenance.tsv (not any-bin).
#  (F2) darkness excluded via the FULL 69.9M dark set (belt-and-braces vs mag_dark subset).
#  (F3) control input sampled PER LENGTH BIN (10aa bins) to match the dark length distribution
#       (dark median 69aa vs non-dark pool 211aa — uniform sampling would confound dN/dS).
#  (F4) all outputs atomic (.part+mv) with hard count guards; mmseqs gated on completion marker.
#  (F5) control CDS via the shared parameterized scripts/slice_cds.py (NOT fam_3's hardcoded copy).
# Pipeline: HQ+MQ contig set -> non-dark MAG-mapped pool + lengths + pool FASTA (one pass) ->
# length-stratified sample (seed 20260813) -> mmseqs 50/80 -> size-binned family matching.
#
# Submit (after fam_1; independent of fam_2/3/4): sbatch scripts/fam_5_control.sh
#
#SBATCH --job-name=ale_fam_control
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam5_%j.out
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=24:00:00
#SBATCH --partition=batch

set -euo pipefail
B=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean
CAT=/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog
PROV=/ibex/scratch/projects/c2014/rob/metagenomics/results/catalogue_provenance.tsv
R=$B/results
export LC_ALL=C

# ---- 1+2. non-dark MAG-mapped pool (HQ+MQ gate, full-dark exclusion) + pool FASTA + lengths ----
if [ ! -f "$R/CTRL_POOL_DONE" ]; then
python3 - <<'PY'
import gzip, csv, os
from collections import defaultdict
CAT="/ibex/scratch/projects/c2014/rob/metagenomics/04_protein_catalog"
PROV="/ibex/scratch/projects/c2014/rob/metagenomics/results/catalogue_provenance.tsv"
R="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results"

def strip(g): return g[:-6] if g.endswith(".fa.gz") else g
hqmq = set()
with open(PROV) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        hqmq.add(strip(row["genome"]))
print(f"[ctrl] HQ+MQ bins: {len(hqmq)} (expect 37949)", flush=True)
assert len(hqmq) == 37949

pref = {"single_sample": "ss", "self_supervised": "selfsup"}
binned = set()
with gzip.open(f"{CAT}/contig_to_bin.tsv.gz", "rt") as fh:
    for line in fh:
        c, strat, sb = line.rstrip("\n").split("\t")
        if f"{pref[strat]}__{c.split('::')[0]}__{sb}" in hqmq:
            binned.add(c)
print(f"[ctrl] HQ+MQ-binned contigs: {len(binned)} (darkclean_5 saw 15,366,244)", flush=True)

dark = set()
with gzip.open(f"{CAT}/darkmatter_clean/04_catalog_dark/dark_rep_ids.txt.gz", "rt") as fh:
    for line in fh: dark.add(line.strip())
print(f"[ctrl] full dark set: {len(dark)} (expect 69,938,628)", flush=True)

# pass over cluster table: rep is MAG-mapped if ANY member contig is HQ+MQ-binned
nd = open(f"{R}/nondark_magmapped_ids.txt.part", "w")
cur, mapped, n = None, False, 0
with open(f"{CAT}/all_linclust90/all_cluster90.tsv") as fh:
    for line in fh:
        rep, _, member = line.partition("\t")
        if rep != cur:
            if cur and mapped and cur not in dark: nd.write(cur + "\n"); n += 1
            cur, mapped = rep, False
        if not mapped:
            m = member.rstrip("\n")
            if m[:m.rfind("_")] in binned: mapped = True
    if cur and mapped and cur not in dark: nd.write(cur + "\n"); n += 1
nd.close()
os.replace(f"{R}/nondark_magmapped_ids.txt.part", f"{R}/nondark_magmapped_ids.txt")
print(f"[ctrl] non-dark MAG-mapped reps: {n}", flush=True)

# one pass over reps FASTA: write pool FASTA + lengths for the non-dark MAG-mapped set
keep = set(l.strip() for l in open(f"{R}/nondark_magmapped_ids.txt"))
pool = open(f"{R}/ctrl_pool.faa.part", "w"); lens = open(f"{R}/ctrl_pool_len.tsv.part", "w")
pid, buf, wrote = None, [], 0
def flush():
    global wrote
    if pid and buf:
        s = "".join(buf); pool.write(f">{pid}\n{s}\n"); lens.write(f"{pid}\t{len(s)}\n"); wrote += 1
with gzip.open(f"{CAT}/all_linclust90/all_reps90.faa.gz", "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            flush()
            h = line[1:].split()[0]
            pid = h if h in keep else None; buf = []
        elif pid: buf.append(line.strip())
flush(); pool.close(); lens.close()
assert wrote == len(keep), f"pool mismatch {wrote} vs {len(keep)}"
os.replace(f"{R}/ctrl_pool.faa.part", f"{R}/ctrl_pool.faa")
os.replace(f"{R}/ctrl_pool_len.tsv.part", f"{R}/ctrl_pool_len.tsv")
print(f"[ctrl] pool FASTA + lengths written: {wrote}", flush=True)
PY
touch "$R/CTRL_POOL_DONE"
fi

# ---- 3. length-stratified sample matching the dark length histogram ----
if [ ! -s "$R/ctrl_reps.faa" ]; then
python3 - <<'PY'
import random, os
random.seed(20260813)
R="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results"
BIN=10
want = {}
pid, ln = None, 0
for line in open(f"{R}/dark_reps.faa"):
    if line.startswith(">"):
        if pid: want[ln//BIN] = want.get(ln//BIN, 0) + 1
        pid, ln = line, 0
    else: ln += len(line.strip())
if pid: want[ln//BIN] = want.get(ln//BIN, 0) + 1
total_want = sum(want.values())
print(f"[ctrl] dark histogram: {total_want} reps, {len(want)} length bins", flush=True)
assert total_want == 4292446

from collections import defaultdict
pool = defaultdict(list)
for line in open(f"{R}/ctrl_pool_len.tsv"):
    p, l = line.rstrip("\n").split("\t")
    pool[int(l)//BIN].append(p)
chosen, short = [], 0
for b, cnt in sorted(want.items()):
    have = pool.get(b, [])
    random.shuffle(have)
    take = have[:cnt]
    if len(take) < cnt:
        short += cnt - len(take)
        print(f"[ctrl] bin {b*BIN}-{b*BIN+BIN-1}aa short by {cnt-len(take)}", flush=True)
    chosen += take
print(f"[ctrl] sampled {len(chosen)} / target {total_want} (shortfall {short}, {100*short/total_want:.2f}%)", flush=True)
sel = set(chosen)
out = open(f"{R}/ctrl_reps.faa.part", "w")
pid, wrote = None, 0
for line in open(f"{R}/ctrl_pool.faa"):
    if line.startswith(">"):
        pid = line[1:].strip() if line[1:].strip() in sel else None
        if pid: out.write(line); wrote += 1
    elif pid: out.write(line)
out.close()
assert wrote == len(chosen), f"extract mismatch {wrote} vs {len(chosen)}"
os.replace(f"{R}/ctrl_reps.faa.part", f"{R}/ctrl_reps.faa")
open(f"{R}/ctrl_sample_n.txt", "w").write(str(len(chosen)) + "\n")
print(f"[ctrl] ctrl_reps.faa written: {wrote}", flush=True)
PY
fi

# ---- 4. cluster identically to the dark arm ----
source /ibex/scratch/projects/c2014/rob/miniconda3/etc/profile.d/conda.sh
conda activate mmseqs2
mkdir -p "$R/ctrl_tmp"
if [ ! -f "$R/CTRL_CLUSTER_DONE" ]; then
    mmseqs easy-cluster "$R/ctrl_reps.faa" "$R/ctrl" "$R/ctrl_tmp" \
        --min-seq-id 0.5 -c 0.8 --cov-mode 0 --threads "$SLURM_CPUS_PER_TASK"
    touch "$R/CTRL_CLUSTER_DONE"
fi

# ---- 5. size-binned family matching (3..9 exact, 10-19, 20+) ----
python3 - <<'PY'
import random
random.seed(20260813)
R="/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/results"
def szbin(m): return str(m) if m <= 9 else ("10-19" if m <= 19 else "20+")
def sizes(path):
    d = {}
    for line in open(path):
        f, m = line.rstrip("\n").split("\t"); d[f] = d.get(f, 0) + 1
    return d
dark = sizes(f"{R}/fam_cluster.tsv")
adm = set(l.strip() for l in open(f"{R}/FAMILIES_FINAL.txt"))
from collections import Counter, defaultdict
want = Counter(szbin(dark[f]) for f in adm)
ctrl = sizes(f"{R}/ctrl_cluster.tsv")
pool = defaultdict(list)
for f, m in ctrl.items():
    if m >= 3: pool[szbin(m)].append(f)
chosen, short = [], 0
for b, cnt in sorted(want.items()):
    have = pool.get(b, [])
    random.shuffle(have)
    take = have[:cnt]
    if len(take) < cnt:
        short += cnt - len(take)
        print(f"[ctrl] size-bin {b} short by {cnt-len(take)}", flush=True)
    chosen += take
open(f"{R}/CONTROL_FAMILIES.txt", "w").write("\n".join(chosen) + "\n")
sel = set(chosen)
with open(f"{R}/ctrl_member_fam.tsv", "w") as out:
    for line in open(f"{R}/ctrl_cluster.tsv"):
        f, m = line.rstrip("\n").split("\t")
        if f in sel: out.write(f"{m}\t{f}\n")
print(f"[ctrl] dark families {len(adm)}  matched {len(chosen)}  shortfall {short}")
PY
echo "[$(date -Iseconds)] step 5 done — control CDS next via scripts/slice_cds.py on ctrl_member_fam.tsv -> ctrl_cds_parts/"
