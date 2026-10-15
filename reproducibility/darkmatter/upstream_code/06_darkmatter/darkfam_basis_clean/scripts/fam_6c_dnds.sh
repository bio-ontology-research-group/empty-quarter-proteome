#!/bin/bash
# Phase C.2 — per-family dN/dS. For each usable family (>=3 member CDS): translate -> mafft protein
# alignment -> back-translate to a codon alignment -> pairwise NG86 dN/dS (biopython cal_dn_ds) over a
# capped set of pairs, with pairwise gap/N stripping and saturation filtering. Batched array.
#   sbatch --export=ARM=dark --array=0-530 scripts/fam_6c_dnds.sh  (verbatim June algorithm, repointed)
#   sbatch --export=ARM=ctrl --array=0-530 scripts/fam_6c_dnds.sh
#
# CALIBRATION ARMS (added 2026-09-20). The paper compares dark (0.425) with a matched annotated
# control (0.388); both are coding, so nothing anchors the scale at 1. ARM=null_frame1 reads the
# SAME dark member CDS in the +1 frame: drop the first nt, trim to a codon multiple, rewrite every
# codon that translates to a stop under table 11 as NNN (the existing pairwise gap/N stripping then
# removes those codons from the counts). All other filters are identical. ARM=dark_subset is the
# same code path with shift 0 (true frame) on the same families = positive control for the masking.
# ARM=dark and ARM=ctrl are untouched: SHIFT=0 and no masking, same index, same output bytes.
#   MODE=prep  ARM=null_frame1 bash scripts/fam_6c_dnds.sh      # -> family_index_null_frame1.tsv
#       (seeded 20260920 random 10,000 of the evaluable >=3-CDS rows of family_index_dark.tsv)
#   sbatch --export=ARM=null_frame1 --array=0-49 scripts/fam_6c_dnds.sh
#   sbatch --export=ARM=dark_subset --array=0-49 scripts/fam_6c_dnds.sh
#   MODE=merge ARM=null_frame1 bash scripts/fam_6c_dnds.sh      # -> dnds_null_frame1.tsv
#   MODE=merge ARM=dark_subset bash scripts/fam_6c_dnds.sh      # -> dnds_dark_subset.tsv
# MODE=merge is refused for dark|ctrl so results/selection/dnds.tsv can never be rewritten here.
#SBATCH --job-name=ale_fam_dnds
#SBATCH --output=/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/darkfam_basis_clean/logs/fam6c_%A_%a.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --partition=batch
set -euo pipefail
WD=/ibex/scratch/projects/c2014/rob/metagenomics
OUT=$WD/06_darkmatter/darkfam_basis_clean/results/selection
MAFFT=$WD/../miniconda3/envs/busco/bin/mafft
BATCH=200
source $WD/../miniconda3/etc/profile.d/conda.sh
conda activate eco_probes
mkdir -p "$OUT/dnds_parts"

ARM=${ARM:?set ARM}
export ARM
MODE=${MODE:-run}

# --- MODE=prep (calibration arms only): seeded 10,000-family subset index -----------------------
if [ "$MODE" = prep ]; then
  [ "$ARM" = null_frame1 ] || { echo "ERROR: MODE=prep is defined only for ARM=null_frame1"; exit 1; }
  python3 - "$OUT" <<'PY'
import os, random, sys
OUT = sys.argv[1]
SEED, N = 20260920, 10000
src = os.path.join(OUT, "family_index_dark.tsv")
hdr, rows = None, []
for line in open(src):
    c = line.rstrip("\n").split("\t")
    if c[0] == "family_id": hdr = line.rstrip("\n"); continue
    if len(c) < 5: continue
    if int(c[3]) >= 3: rows.append(c)          # evaluable = the worker's own >=3 CDS gate
rows.sort(key=lambda c: c[0])                  # fixed order so the draw is reproducible
pick = sorted(random.Random(SEED).sample(rows, N), key=lambda c: c[0])
dst = os.path.join(OUT, "family_index_null_frame1.tsv")
with open(dst + ".part", "w") as out:
    out.write(hdr + "\n")
    for c in pick:
        out.write("\t".join([c[0], "null_frame1", c[2], c[3], c[4]]) + "\n")
os.replace(dst + ".part", dst)
print(f"[prep] evaluable dark families {len(rows)}, sampled {len(pick)} (seed {SEED}) -> {dst}")
PY
  exit 0
fi

# --- MODE=merge (calibration arms only) ---------------------------------------------------------
if [ "$MODE" = merge ]; then
  case "$ARM" in
    null_frame1|dark_subset) ;;
    *) echo "ERROR: MODE=merge refused for ARM=$ARM (dnds.tsv is not rewritten by this script)"; exit 1;;
  esac
  M=$OUT/dnds_$ARM.tsv
  shopt -s nullglob
  parts=("$OUT"/dnds_parts/part_"$ARM"_[0-9]*.tsv)      # glob, not ls|head: pipefail + SIGPIPE
  [ ${#parts[@]} -gt 0 ] || { echo "ERROR: no dnds_parts for ARM=$ARM"; exit 1; }
  head -1 "${parts[0]}" > "$M.part"
  for f in "${parts[@]}"; do tail -n +2 "$f"; done >> "$M.part"
  mv "$M.part" "$M"
  echo "[$(date -Iseconds)] merged $ARM: $(( $(wc -l < "$M") - 1 )) families -> $M"
  exit 0
fi

python3 - "$OUT" "$MAFFT" "$BATCH" "${SLURM_ARRAY_TASK_ID:-0}" <<'PY'
import os, sys, subprocess, tempfile, zlib
from Bio.Seq import Seq
from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
OUT, MAFFT, BATCH, TASK = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
ARM = os.environ["ARM"]
# calibration arms: both run on the seeded 10,000-family subset index; null_frame1 is read in the
# +1 frame. dark|ctrl keep SHIFT=0 / MASK=False / their own index, i.e. the original code path.
CALIB = ARM in ("null_frame1", "dark_subset")
SHIFT = 1 if ARM == "null_frame1" else 0
MASK = CALIB
IDX = "family_index_null_frame1.tsv" if CALIB else "family_index_" + ARM + ".tsv"

fams = []
for line in open(os.path.join(OUT, IDX)):
    c = line.rstrip("\n").split("\t")
    if c[0] == "family_id" or len(c) < 5: continue
    if int(c[3]) >= 3: fams.append((c[0], ARM if CALIB else c[1], c[4]))   # id, set, fna
fams.sort()
chunk = fams[TASK*BATCH:(TASK+1)*BATCH]
out = open(os.path.join(OUT, "dnds_parts", "part_" + ARM + f"_{TASK}.tsv"), "w")
out.write("family_id\tset\tn_seq\tn_pairs\tmean_dN\tmean_dS\tmean_dNdS\tmedian_dNdS\tmean_pid"
          + ("\tfrac_stop_codons\n" if CALIB else "\n"))

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
    # SHIFT=0 and MASK=False for dark|ctrl, so their path through this function is unchanged.
    # Returns (codon alignment, codons kept, codons masked as stops).
    prot = {}; cds = {}; ncod = 0; nmask = 0
    for nm, nt in seqs:
        nt = nt[SHIFT:]                       # +1 frame for null_frame1
        nt = nt[:len(nt)//3*3]
        if len(nt) < 90: continue
        nmk = 0
        if MASK:                              # stop codons -> NNN, dropped later by pair_dnds
            t = str(Seq(nt).translate(table=11))
            if "*" in t:
                b = list(nt)
                for i, aa in enumerate(t):
                    if aa == "*": b[3*i:3*i+3] = "NNN"
                nt = "".join(b); nmk = t.count("*")
        p = str(Seq(nt).translate(table=11)).rstrip("*")
        if "*" in p or len(p) < 25: continue
        prot[nm] = p; cds[nm] = nt; ncod += len(nt)//3; nmask += nmk
    if len(prot) < 3: return {}, ncod, nmask
    with tempfile.NamedTemporaryFile("w", suffix=".faa", delete=False) as tf:
        for nm, p in prot.items(): tf.write(f">{nm}\n{p}\n")
        fa = tf.name
    try:
        al = subprocess.run([MAFFT, "--auto", "--quiet", "--thread", "1", fa],
                            capture_output=True, text=True, timeout=300).stdout
    except Exception:
        os.unlink(fa); return {}, ncod, nmask
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
    return codon, ncod, nmask

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
n_missing = n_fail_min3 = n_no_pairs = n_written = 0
for fid, fset, fna in chunk:
    if not os.path.exists(fna): n_missing += 1; continue
    seqs = read_fna(fna)
    codon, ncod, nmask = codon_align(seqs)
    if len(codon) < 3: n_fail_min3 += 1; continue
    ids = sorted(codon)
    pairs = [(ids[i], ids[j]) for i in range(len(ids)) for j in range(i+1, len(ids))]
    if len(pairs) > 80:                                   # deterministic cap
        pairs = [pairs[k] for k in range(0, len(pairs), max(1, len(pairs)//80))][:80]
    dn, ds, rr, pids = [], [], [], []
    for a, b in pairs:
        r = pair_dnds(codon[a], codon[b])
        if r:
            dn.append(r[0]); ds.append(r[1]); rr.append(r[0]/r[1]); pids.append(r[2])
    if len(rr) < 1: n_no_pairs += 1; continue
    out.write(f"{fid}\t{fset}\t{len(codon)}\t{len(rr)}\t{st.mean(dn):.4f}\t{st.mean(ds):.4f}\t"
              f"{st.mean(rr):.4f}\t{st.median(rr):.4f}\t{st.mean(pids):.4f}"
              + (f"\t{nmask/ncod:.4f}\n" if CALIB else "\n"))
    n_written += 1
out.close()
print(f"[dnds] task {TASK}: {len(chunk)} families processed")
if CALIB:                                     # per-task accounting for the calibration report
    qd = os.path.join(OUT, "dnds_qc"); os.makedirs(qd, exist_ok=True)
    with open(os.path.join(qd, f"qc_{ARM}_{TASK}.tsv"), "w") as q:
        q.write("arm\ttask\tn_families\tn_missing_fna\tn_fail_min3seq\tn_no_usable_pairs\tn_written\n")
        q.write(f"{ARM}\t{TASK}\t{len(chunk)}\t{n_missing}\t{n_fail_min3}\t{n_no_pairs}\t{n_written}\n")
    print(f"[dnds] {ARM} task {TASK}: missing {n_missing}, fail_min3seq {n_fail_min3}, "
          f"no_pairs {n_no_pairs}, written {n_written}")
PY
echo "[$(date -Iseconds)] dnds task ${SLURM_ARRAY_TASK_ID:-0} done"
