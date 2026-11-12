#!/usr/bin/env python3
"""Stage 2: per-protein mature start, TMbed topology and mature-region sequence features for EQ members.

TMbed labels come from (a) the reused 2026-10-30 predictions matched by exact sequence and (b) the new
TMbed 1.0.2 CPU run of this analysis (work/tmbed_new.pred.gz). A QC sample predicted by both is compared.
Output: work/eq_protein_features.tsv.gz, work/features_info.json
"""
import argparse, gzip, glob, json, math, re, sys, hashlib
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

EQ = "EQ_Thermoplasmatota_237"
EXCLUDED = {"External_Comparators_CrossResource_Duplicate_54", "SPIRE_ProGenomes3_SpecI_Excluded_148"}
CS_RE = re.compile(r"CS pos: (\d+)-(\d+)\. Pr: ([0-9.]+)")
SEQUON = re.compile(r"(?=N[^P][ST])")
# Bjellqvist pK values (as in ExPASy ProtParam / Biopython IsoelectricPoint)
PKP = {"K": 10.0, "R": 12.0, "H": 5.98}
PKN = {"D": 4.05, "E": 4.45, "C": 9.0, "Y": 10.0}
PK_NTERM, PK_CTERM = 7.5, 3.55


def pI(seq):
    if not seq:
        return float("nan")
    cp = {a: seq.count(a) for a in PKP}
    cn = {a: seq.count(a) for a in PKN}
    lo, hi = 0.0, 14.0
    for _ in range(40):
        ph = (lo + hi) / 2
        ch = 1 / (1 + 10 ** (ph - PK_NTERM)) - 1 / (1 + 10 ** (PK_CTERM - ph))
        for a, n in cp.items():
            ch += n / (1 + 10 ** (ph - PKP[a]))
        for a, n in cn.items():
            ch -= n / (1 + 10 ** (PKN[a] - ph))
        if ch > 0:
            lo = ph
        else:
            hi = ph
    return round((lo + hi) / 2, 2)


def lowcomplexity(seq, w=12, thr=2.2):
    """Fraction of residues covered by >=1 window of length w with Shannon entropy <= thr bits (SEG-like)."""
    n = len(seq)
    if n < w:
        return 0.0
    cov = bytearray(n)
    cnt = {}
    for c in seq[:w]:
        cnt[c] = cnt.get(c, 0) + 1
    plogp = [0.0] + [-(k / w) * math.log2(k / w) for k in range(1, w + 1)]
    i = 0
    while True:
        h = sum(plogp[k] for k in cnt.values())
        if h <= thr:
            cov[i:i + w] = b"\x01" * w
        if i + w >= n:
            break
        a, b = seq[i], seq[i + w]
        cnt[a] -= 1
        if not cnt[a]:
            del cnt[a]
        cnt[b] = cnt.get(b, 0) + 1
        i += 1
    return sum(cov) / n


def st_window_max(seq, w=30):
    n = len(seq)
    if n < w:
        return (seq.count("S") + seq.count("T")) / n if n else 0.0
    v = np.frombuffer(seq.encode(), dtype="S1")
    x = ((v == b"S") | (v == b"T")).astype(np.int32)
    c = np.concatenate([[0], np.cumsum(x)])
    return float((c[w:] - c[:-w]).max()) / w


def helix_segments(lab):
    """Runs of H or h (a change between H and h starts a new helix). 0-based [start, end)."""
    out = []
    i, n = 0, len(lab)
    while i < n:
        c = lab[i]
        if c in "Hh":
            j = i
            while j < n and lab[j] == c:
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out


def feat2(args):
    gid, seq, lab, ms, c_complete = args
    L = len(seq)
    mat = seq[ms:]
    r = {"gene_id": gid}
    if lab is not None:
        hs = helix_segments(lab)
        mature_h = [(a, b) for a, b in hs if (a + b) / 2 >= ms]
        pre_h = [(a, b) for a, b in hs if (a + b) / 2 < ms]
        r["tm_helices_total"] = len(hs)
        r["tm_helices_mature"] = len(mature_h)
        r["tm_helix_in_signal_region"] = len(pre_h)
        r["tm_mature_coords"] = ";".join(f"{a + 1}-{b}" for a, b in mature_h)
        last_end = mature_h[-1][1] if mature_h else 0
        r["cterm_tm"] = bool(mature_h) and (L - last_end) <= 30
        r["cterm_tail_len"] = (L - last_end) if mature_h else -1
        r["tm_beta_strands"] = len(re.findall(r"B+|b+", lab))
    else:
        for k in ("tm_helices_total", "tm_helices_mature", "tm_helix_in_signal_region", "tm_beta_strands"):
            r[k] = -1
        r["tm_mature_coords"] = ""
        r["cterm_tm"] = False
        r["cterm_tail_len"] = -1
    ml = len(mat)
    r["mature_length"] = ml
    r["mature_pI"] = pI(mat)
    nseq = len(SEQUON.findall(mat))
    r["n_sequons"] = nseq
    r["sequons_per_100aa"] = round(100 * nseq / ml, 3) if ml else float("nan")
    r["st_fraction"] = round((mat.count("S") + mat.count("T")) / ml, 4) if ml else float("nan")
    r["st_max_window30"] = round(st_window_max(mat), 4) if ml else float("nan")
    r["lowcomplexity_fraction"] = round(lowcomplexity(mat), 4) if ml else float("nan")
    # heuristic C-terminal sorting-signal-like tail: PGF/PEF tripeptide, then >=12 hydrophobic-rich residues,
    # then a short tail containing K/R, all within the last 45 residues (C-complete proteins only).
    tail = seq[-45:]
    m = re.search(r"P[GE]F.{0,6}([AVLIFMWGSTC]{12,25}).{0,3}[KR].{0,8}$", tail)
    r["cterm_PGF_PEF_motif_heuristic"] = bool(m) and bool(c_complete)
    return r




def seqhash(s):
    return hashlib.sha256(s.encode()).hexdigest()[:24]


def read_pred(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip(); lab = fh.readline().strip()
            assert h.startswith(">") and len(s) == len(lab), h
            out[h[1:].split()[0]] = (s, lab)
    return out


def nhel(lab):
    return len(helix_segments(lab))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=64)
    a = ap.parse_args()
    W = a.work
    info = {}
    mem = pd.read_csv(W / "members.tsv.gz", sep="\t", dtype={"gene_id": str})
    eqm = mem[mem.is_eq].set_index("gene_id")
    seqs = {}
    g = None
    for line in open(W / "eq_members_full.faa"):
        if line[0] == ">":
            g = line[1:].strip()
        else:
            seqs[g] = line.strip()
    assert set(seqs) == set(eqm.index)
    reused = dict(pd.read_csv(W / "tmbed_reused.tsv.gz", sep="\t", dtype=str).values)
    new = {}
    for h, (s, lab) in read_pred(W / "tmbed_new.pred.gz").items():
        assert seqhash(s) == h, h
        new[h] = lab
    qc = [l.strip() for l in open(W / "tmbed_qc_hashes.txt") if l.strip()]
    both = [h for h in qc if h in new and h in reused]
    info["tmbed_reused_sequences"] = len(reused); info["tmbed_new_sequences"] = len(new)
    info["tmbed_qc_n"] = len(both)
    if both:
        info["tmbed_qc_identical_label_string"] = sum(new[h] == reused[h] for h in both)
        info["tmbed_qc_same_helix_count"] = sum(nhel(new[h]) == nhel(reused[h]) for h in both)
        info["tmbed_qc_same_signal_presence"] = sum(("S" in new[h]) == ("S" in reused[h]) for h in both)
        info["tmbed_qc_mean_per_residue_identity"] = round(float(np.mean(
            [np.mean([x == y for x, y in zip(new[h], reused[h])]) for h in both])), 5)
    labs = {}; src = {}
    for gid, s in seqs.items():
        h = seqhash(s)
        if h in reused:
            labs[gid] = reused[h]; src[gid] = "reused_20261030"
        elif h in new:
            labs[gid] = new[h]; src[gid] = "new_tmbed_1.0.2_cpu"
    sp_pos = eqm.signalp_class != "OTHER"
    info["tmbed_coverage_eq_members"] = [len(labs), len(seqs)]
    info["tmbed_coverage_sp_pos"] = [int(sum(g_ in labs for g_ in eqm.index[sp_pos])), int(sp_pos.sum())]
    info["tmbed_source_counts_sp_pos"] = pd.Series([src.get(g_, "none") for g_ in eqm.index[sp_pos]]).value_counts().to_dict()

    rows = []
    for gid in eqm.index:
        lab = labs.get(gid); cls = eqm.signalp_class[gid]; cs = int(eqm.cs_pos[gid])
        s_end = 0
        if lab is not None:
            m = re.match(r"S+", lab)
            s_end = m.end() if m else 0
        if cls == "OTHER":
            ms, sc = 0, "no_signal_peptide"
        elif cs > 0:
            ms, sc = cs, "signalp_cs"
        elif s_end > 0:
            ms, sc = s_end, "tmbed_signal_end"
        else:
            ms, sc = -1, "imputed_class_median"
        rows.append((gid, len(seqs[gid]), ms, sc, lab is not None, s_end, src.get(gid, "none")))
    d = pd.DataFrame(rows, columns=["gene_id", "aa_length", "mature_start", "mature_start_source", "tmbed_available",
                                    "tmbed_signal_len", "tmbed_source"]).set_index("gene_id")
    d["signalp_class"] = eqm.signalp_class
    med = d[d.mature_start_source == "signalp_cs"].groupby("signalp_class").mature_start.median().round().astype(int)
    info["class_median_cs"] = med.to_dict()
    imp = d.mature_start < 0
    d.loc[imp, "mature_start"] = d.loc[imp, "signalp_class"].map(med).fillna(25).astype(int)
    d["mature_start"] = np.minimum(d.mature_start, d.aa_length - 1).clip(lower=0)
    info["mature_start_source_counts_sp_pos"] = d[d.signalp_class != "OTHER"].mature_start_source.value_counts().to_dict()
    jobs = [(g_, seqs[g_], labs.get(g_), int(d.mature_start[g_]), bool(eqm.c_complete[g_])) for g_ in d.index]
    with Pool(a.threads) as p:
        f2 = pd.DataFrame(p.map(feat2, jobs, chunksize=500)).set_index("gene_id")
    out = eqm[["family_rep", "genome_id", "dark", "signalp_class", "sp_prob", "cs_pos", "cs_prob", "signalp_status",
               "partial_flag", "n_complete", "c_complete"]].join(d.drop(columns=["signalp_class"])).join(f2)
    out.to_csv(W / "eq_protein_features.tsv.gz", sep="\t")
    info["n_eq_feature_rows"] = int(len(out))
    json.dump(info, open(W / "features_info.json", "w"), indent=1, default=str)
    print("FEATURES_DONE", json.dumps(info, default=str), flush=True)


if __name__ == "__main__":
    main()
