#!/bin/bash
# Phase C.1 — extract member nucleotide CDS for multi-member DARK families (rep is sequence-novel,
# >=3 members) and a size-matched NON-DARK control set (the null), for within-family dN/dS. Members
# come from the 5,438 dRep genomes (orthogroups were built from dRep rep proteomes); CDS sliced from
# the genome FASTA via the pyrodigal GFF, as in 90. Members capped (~40/family, deterministic) to
# bound compute; big families thus give precise dN/dS, small ones contribute to the distribution.
# LEGACY/WITHDRAWN: the control arm is sample-confounded and arm assignment is
# based on UniRef detectability, itself related to evolutionary rate.  This
# script is retained only to reproduce the historical table; it is not an
# approved analysis path for new claims.
#   sbatch scripts/ibex/93_darkfam_member_cds.sh
#SBATCH --job-name=hohndor_darkfam_memcds
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_memcds_%j.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_memcds_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --partition=batch
set -euo pipefail
if [ "${ALLOW_WITHDRAWN_DNDS_PIPELINE:-0}" != "1" ]; then
    echo "withdrawn dN/dS control pipeline; set ALLOW_WITHDRAWN_DNDS_PIPELINE=1 only for explicit historical reproduction" >&2
    exit 64
fi
WD=/ibex/scratch/projects/c2014/rob/metagenomics
DM=$WD/06_darkmatter
OG=$DM/orthogroups
GEN=$WD/05_refine/_drep/out/dereplicated_genomes
DRP=$DM/drep_rep_proteomes
DARK=$DM/reallydark/isdark_ids.txt
OUT=$DM/darkfam_ecology/selection
mkdir -p "$OUT/members/dark" "$OUT/members/control"
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eco_probes
CAP=40; MIN=3

python3 - "$OG" "$GEN" "$DRP" "$DARK" "$OUT" "$CAP" "$MIN" <<'PY'
import os, sys, gzip, zlib
OG, GEN, DRP, DARK, OUT, CAP, MIN = sys.argv[1:8]
CAP = int(CAP); MIN = int(MIN)

dark = set(l.strip() for l in open(DARK) if l.strip())
# classify focal reps from og_recurrence: dark (>=MIN members) + size-matched non-dark control
nmemb = {}
dark_reps = []; ctrl_pool = []
for line in open(os.path.join(OG, "og_recurrence.tsv")):
    c = line.rstrip("\n").split("\t")
    if len(c) < 3: continue
    try: nm = int(c[1])
    except ValueError: continue
    if nm < MIN: continue
    nmemb[c[0]] = nm
    (dark_reps if c[0] in dark else ctrl_pool).append(c[0])
print(f"[memcds] dark families: {len(dark_reps)} ; non-dark pool: {len(ctrl_pool)}", flush=True)

# size-match control to dark by binning on log2(n_members); sample equal count per bin
import math
def bink(n): return int(math.log2(n))
from collections import defaultdict
dbins = defaultdict(int)
for r in dark_reps: dbins[bink(nmemb[r])] += 1
cby = defaultdict(list)
for r in ctrl_pool: cby[bink(nmemb[r])].append(r)
ctrl = []
for b, want in dbins.items():
    pool = sorted(cby.get(b, []))          # deterministic
    take = pool[:want] if len(pool) >= want else pool
    ctrl.extend(take)
focal = {r: "dark" for r in dark_reps}
for r in ctrl: focal[r] = "control"
print(f"[memcds] size-matched control families: {len(ctrl)}", flush=True)

# stream cluster membership; include each member with prob CAP/n_members (deterministic via crc32)
by_mag = defaultdict(list)   # mag -> [(rep, locus)]
fam_n = defaultdict(int)
for line in open(os.path.join(OG, "orthogroups_cluster.tsv")):
    rep, mem = line.rstrip("\n").split("\t")[:2]
    grp = focal.get(rep)
    if grp is None: continue
    nm = nmemb[rep]
    if nm > CAP:
        if (zlib.crc32(mem.encode()) % nm) >= CAP:   # keep ~CAP per family
            continue
    mag, locus = mem.split("::", 1)
    by_mag[mag].append((rep, locus)); fam_n[rep] += 1
print(f"[memcds] members to extract across {len(by_mag)} MAGs", flush=True)

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]
def rd(p):
    n, s = None, []
    op = gzip.open(p, "rt") if p.endswith(".gz") else open(p, errors="ignore")
    for line in op:
        if line.startswith(">"):
            if n is not None: yield n, "".join(s)
            n = line[1:].split()[0]; s = []
        else: s.append(line.strip())
    if n is not None: yield n, "".join(s)

# group target loci per mag, parse GFF once, slice CDS, append to per-family fna
fh = {}   # rep -> open file handle
def fam_path(rep):
    grp = focal[rep]; safe = rep.replace("::", "_").replace("/", "_")
    return os.path.join(OUT, "members", grp, safe + ".fna")
nseq = 0
for mag, items in by_mag.items():
    g = os.path.join(GEN, mag + ".fa"); gff = os.path.join(DRP, mag + ".gff.gz")
    if not (os.path.exists(g) and os.path.exists(gff)): continue
    want = {}
    for rep, locus in items: want.setdefault(locus, []).append(rep)
    contigs = {n: s for n, s in rd(g)}
    coords = {}
    for line in gzip.open(gff, "rt"):
        if line.startswith("#") or "\tCDS\t" not in line: continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9: continue
        i = c[8].find("ID=")
        if i < 0: continue
        loc = c[8][i+3:].split(";")[0]
        if loc in want: coords[loc] = (c[0], int(c[3]), int(c[4]), c[6])
    for locus, reps in want.items():
        cc = coords.get(locus)
        if not cc: continue
        contig, s, e, strand = cc; seq = contigs.get(contig)
        if not seq: continue
        nt = seq[s-1:e]
        if strand == "-": nt = rc(nt)
        if len(nt) < 90: continue           # need a sensible CDS
        for rep in reps:
            if rep not in fh: fh[rep] = open(fam_path(rep), "a")
            fh[rep].write(f">{mag}::{locus}\n{nt}\n"); nseq += 1
for h in fh.values(): h.close()

# write family index (only families that ended with >=3 extracted members are usable)
idx = open(os.path.join(OUT, "family_index.tsv"), "w")
idx.write("family_id\tset\tn_members_total\tn_members_extracted\tfna\n")
usable = 0
for rep, grp in focal.items():
    p = fam_path(rep)
    n = 0
    if os.path.exists(p):
        n = sum(1 for l in open(p) if l.startswith(">"))
    if n >= 3: usable += 1
    idx.write(f"{rep}\t{grp}\t{nmemb[rep]}\t{n}\t{p}\n")
idx.close()
print(f"[memcds] extracted {nseq} member CDS ; usable families (>=3 extracted): {usable}", flush=True)
PY

echo "=== counts ==="
echo "dark family fna:"; ls "$OUT/members/dark" 2>/dev/null | wc -l
echo "control family fna:"; ls "$OUT/members/control" 2>/dev/null | wc -l
tail -n +2 "$OUT/family_index.tsv" | awk -F"\t" '$4>=3{print $2}' | sort | uniq -c
date -Iseconds > "$OUT/MEMCDS_DONE"
echo "[$(date -Iseconds)] member CDS extraction done -> $OUT"
