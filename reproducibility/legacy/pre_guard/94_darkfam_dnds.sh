#!/bin/bash
# Phase C.2 — per-family dN/dS. For each usable family (>=3 member CDS): translate -> mafft protein
# alignment -> back-translate to a codon alignment -> pairwise NG86 dN/dS (biopython cal_dn_ds) over a
# capped set of pairs, with pairwise gap/N stripping and saturation filtering. Batched array.
#   sbatch --array=0-N%30 scripts/ibex/94_darkfam_dnds.sh     (N = ceil(usable/BATCH)-1)
#SBATCH --job-name=hohndor_darkfam_dnds
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_dnds_%A_%a.out
#SBATCH --error=/ibex/scratch/projects/c2014/rob/metagenomics/logs/darkfam_dnds_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$WD/06_darkmatter/darkfam_ecology/selection
MAFFT=$WD/../miniconda3/envs/busco/bin/mafft
BATCH=200
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eco_probes
mkdir -p "$OUT/dnds_parts"

python3 - "$OUT" "$MAFFT" "$BATCH" "${SLURM_ARRAY_TASK_ID:-0}" <<'PY'
import os, sys, subprocess, tempfile, zlib
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
OUT, MAFFT, BATCH, TASK = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])

fams = []
for line in open(os.path.join(OUT, "family_index.tsv")):
    c = line.rstrip("\n").split("\t")
    if c[0] == "family_id" or len(c) < 5: continue
    if int(c[3]) >= 3: fams.append((c[0], c[1], c[4]))   # id, set, fna
fams.sort()
chunk = fams[TASK*BATCH:(TASK+1)*BATCH]
out = open(os.path.join(OUT, "dnds_parts", f"part_{TASK}.tsv"), "w")
out.write("family_id\tset\tn_seq\tn_pairs\tmean_dN\tmean_dS\tmean_dNdS\tmedian_dNdS\tmean_pid\n")

def read_fna(p):
    n, s, d = None, [], []
    for line in open(p):
        if line.startswith(">"):
            if n is not None: d.append((n, "".join(s)))
            n = line[1:].strip(); s = []
        else: s.append(line.strip())
    if n is not None: d.append((n, "".join(s)))
    return d

def codon_align(seqs):
    # translate (table 11, trim to codon), mafft proteins, back-translate to codon alignment
    prot = {}; cds = {}
    for nm, nt in seqs:
        nt = nt[:len(nt)//3*3]
        if len(nt) < 90: continue
        p = str(Seq(nt).translate(table=11)).rstrip("*")
        if "*" in p or len(p) < 25: continue
        prot[nm] = p; cds[nm] = nt
    if len(prot) < 3: return {}
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as tf:
        for nm, p in prot.items(): tf.write(f">{nm}\n{p}\n")
        fa = tf.name
    try:
        al = subprocess.run([MAFFT, "--auto", "--quiet", "--thread", "1", fa],
                            capture_output=True, text=True, timeout=300).stdout
    except Exception:
        os.unlink(fa); return {}
    os.unlink(fa)
    aln = {}; cur = None; buf = []
    for line in al.splitlines():
        if line.startswith(">"):
            if cur is not None: aln[cur] = "".join(buf)
            cur = line[1:].strip(); buf = []
        else: buf.append(line.strip())
    if cur is not None: aln[cur] = "".join(buf)
    codon = {}
    for nm, ap in aln.items():
        nt = cds.get(nm);
        if not nt: continue
        out_s = []; i = 0
        for aa in ap:
            if aa == "-": out_s.append("---")
            else:
                out_s.append(nt[i:i+3] if i+3 <= len(nt) else "---"); i += 3
        codon[nm] = "".join(out_s)
    return codon

def pair_dnds(a, b):
    # strip codons where either has gap/N; require divergence; NG86
    ca, cb = [], []
    L = min(len(a), len(b))
    for k in range(0, L-2, 3):
        x, y = a[k:k+3], b[k:k+3]
        if "-" in x or "-" in y or "N" in x.upper() or "N" in y.upper(): continue
        ca.append(x); cb.append(y)
    if len(ca) < 30: return None
    sa, sb = "".join(ca), "".join(cb)
    if sa == sb: return None
    pid = sum(1 for x, y in zip(sa, sb) if x == y) / len(sa)
    try:
        dN, dS = cal_dn_ds(CodonSeq(sa), CodonSeq(sb), method="NG86")
    except Exception:
        return None
    if dS is None or dS <= 0 or dS > 2 or dN is None or dN < 0: return None
    return dN, dS, pid

import statistics as st
for fid, fset, fna in chunk:
    if not os.path.exists(fna): continue
    seqs = read_fna(fna)
    codon = codon_align(seqs)
    if len(codon) < 3: continue
    ids = sorted(codon)
    pairs = [(ids[i], ids[j]) for i in range(len(ids)) for j in range(i+1, len(ids))]
    if len(pairs) > 80:                                   # deterministic cap
        pairs = [pairs[k] for k in range(0, len(pairs), max(1, len(pairs)//80))][:80]
    dn, ds, rr, pids = [], [], [], []
    for a, b in pairs:
        r = pair_dnds(codon[a], codon[b])
        if r:
            dn.append(r[0]); ds.append(r[1]); rr.append(r[0]/r[1]); pids.append(r[2])
    if len(rr) < 1: continue
    out.write(f"{fid}\t{fset}\t{len(codon)}\t{len(rr)}\t{st.mean(dn):.4f}\t{st.mean(ds):.4f}\t"
              f"{st.mean(rr):.4f}\t{st.median(rr):.4f}\t{st.mean(pids):.4f}\n")
out.close()
print(f"[dnds] task {TASK}: {len(chunk)} families processed")
PY
echo "[$(date -Iseconds)] dnds task ${SLURM_ARRAY_TASK_ID:-0} done"
