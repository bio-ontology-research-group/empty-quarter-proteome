#!/usr/bin/env python3
"""Step 3: a3m alignments (monomer, paired complexes) and depth statistics.

Input: MMseqs2 result2msa a3m database (queries = mature sequences, targets = all 3,617,011 proteins),
optionally per-query a3m files from a local colabfold_search (UniRef30 + ColabFold envdb).
Rules
- a hit row is kept when it covers >= 50 % of the query match columns;
- redundancy: greedy in hit order, a row is dropped when >= 95 % of its aligned residues are identical
  to an already kept row (the query is always kept);
- Neff: sum over kept rows of 1 / (number of kept rows with >= 80 % identity over all query columns);
- real pairs: one row pair per genome, the two hits must be adjacent genes (rank difference 1, same contig);
- control pairs: one row pair per genome, the best hit of each chain in that genome (different genes);
- control paired depth is capped at the non-redundant paired depth of the matched real pair;
- unpaired blocks: the filtered single-chain alignments, padded with gaps (ColabFold paired a3m convention).
"""
import argparse, json, os, tarfile, io, sys
from pathlib import Path
from multiprocessing import Pool
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--step1", type=Path, required=True)
ap.add_argument("--msa-db", type=Path, required=True)
ap.add_argument("--qdb-lookup", type=Path, required=True)
ap.add_argument("--extra-dir", type=Path, default=None, help="colabfold_search output directory (QID.a3m)")
ap.add_argument("--out", type=Path, required=True)
ap.add_argument("--tag", required=True)
ap.add_argument("--threads", type=int, default=96)
ap.add_argument("--max-mono", type=int, default=5000)
ap.add_argument("--max-unpaired", type=int, default=2000)
ap.add_argument("--max-paired", type=int, default=2000)
a = ap.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
GAP = ord("-")

G = pd.read_pickle(a.step1 / "genes_min.pkl").set_index("gene_id")
Q = pd.read_csv(a.step1 / "queries.tsv", sep="\t").set_index("qid")
REPS = pd.read_csv(a.step1 / "representatives.tsv", sep="\t")
PAIRS = pd.read_csv(a.step1 / "pairs.tsv", sep="\t")
FAMC = pd.read_csv(a.step1 / "family_codes.tsv", sep="\t", header=None, index_col=1)[0]
qseq = {}
with open(a.step1 / "queries.faa") as f:
    for line in f:
        if line[0] == ">": k = line[1:].strip()
        else: qseq[k] = line.strip()
lookup = pd.read_csv(a.qdb_lookup, sep="\t", header=None, usecols=[0, 1], names=["key", "name"])
key_of = dict(zip(lookup.name, lookup.key))
idx = pd.read_csv(str(a.msa_db) + ".index", sep="\t", header=None, names=["key", "off", "len"]).set_index("key")
data_files = [a.msa_db] if a.msa_db.exists() else sorted(a.msa_db.parent.glob(a.msa_db.name + ".[0-9]*"),
                                                         key=lambda p: int(p.suffix[1:]))
sizes = np.cumsum([0] + [p.stat().st_size for p in data_files])

def read_entry(qid):
    r = idx.loc[key_of[qid]]
    k = int(np.searchsorted(sizes, r.off, side="right") - 1)
    with open(data_files[k], "rb") as f:
        f.seek(int(r.off - sizes[k]))
        return f.read(int(r.len)).rstrip(b"\x00").decode()

def parse_a3m(txt):
    names, rows = [], []
    for block in txt.split("\x00"):
        lines = block.splitlines()
        i = 0
        while i < len(lines):
            if lines[i].startswith(">"):
                names.append(lines[i][1:].split()[0] if len(lines[i]) > 1 else "")
                rows.append(lines[i + 1] if i + 1 < len(lines) else "")
                i += 2
            else:
                i += 1
    return names, rows

_tbl = bytes.maketrans(b"", b"")
_lower = bytes(range(ord("a"), ord("z") + 1))
def match_cols(row):
    return np.frombuffer(row.encode().translate(_tbl, _lower), dtype=np.uint8)

def load_hits(qid):
    """Rows of the own-protein search (genome-resolved), query coverage >= 0.5, in hit order."""
    q = qseq[qid]; L = len(q)
    names, rows = parse_a3m(read_entry(qid))
    out_n, out_r, out_m = [], [], []
    seen = set()
    for nm, r in zip(names[1:], rows[1:]):     # first row is the query
        if nm in seen or nm not in G.index: continue
        m = match_cols(r)
        if len(m) != L: continue
        if (m != GAP).sum() < 0.5 * L: continue
        seen.add(nm); out_n.append(nm); out_r.append(r.rstrip("abcdefghijklmnopqrstuvwxyz")); out_m.append(m)
    return out_n, out_r, (np.vstack(out_m) if out_m else np.zeros((0, L), np.uint8)), len(names) - 1

EXTRA = {}
if a.extra_dir is not None:      # map files by their first (query) sequence, whatever the file naming
    for p in a.extra_dir.glob("*.a3m"):
        with open(p, errors="replace") as fh:
            for line in fh:
                line = line.replace("\x00", "").strip()
                if line and line[0] not in ">#":
                    EXTRA[line.translate(str.maketrans("", "", "abcdefghijklmnopqrstuvwxyz-"))] = p
                    break

def load_extra(qid):
    if a.extra_dir is None: return [], [], None
    p = EXTRA.get(qseq[qid])
    if p is None: return [], [], None
    q = qseq[qid]; L = len(q)
    names, rows = parse_a3m(p.read_text(errors="replace"))
    out_n, out_r, out_m, seen = [], [], [], set()
    for nm, r in zip(names, rows):
        m = match_cols(r)
        if len(m) != L or (m != GAP).sum() < 0.5 * L: continue
        key = m.tobytes()
        if key in seen or nm in ("101", qid): continue
        seen.add(key); out_n.append("ext|" + nm); out_r.append(r.rstrip("abcdefghijklmnopqrstuvwxyz")); out_m.append(m)
    return out_n, out_r, (np.vstack(out_m) if out_m else np.zeros((0, L), np.uint8))

def greedy95(qm, M, thr=0.95, cap=None):
    """Indices of rows of M kept; query row qm seeds the kept set."""
    L = M.shape[1]
    kept = np.empty((min(len(M), cap or len(M)) + 1, L), np.uint8); kept[0] = qm; nk = 1
    keep = []
    for i in range(len(M)):
        r = M[i]; ng = r != GAP; d = int(ng.sum())
        if d == 0: continue
        if ((kept[:nk][:, ng] == r[ng]).sum(1) >= thr * d).any(): continue
        if cap is not None and len(keep) >= cap: break
        kept[nk] = r; nk += 1; keep.append(i)
    return np.array(keep, dtype=int)

def neff80(qm, M, cap=5000):
    X = np.vstack([qm[None], M[:cap]])
    N, L = X.shape
    cnt = np.zeros(N)
    for s in range(0, N, 256):
        blk = X[s:s + 256]
        ident = np.zeros((len(blk), N), np.int32)
        for c0 in range(0, L, 64):
            ident += (blk[:, None, c0:c0 + 64] == X[None, :, c0:c0 + 64]).sum(2)
        cnt[s:s + 256] = (ident >= 0.8 * L).sum(1)
    return float((1.0 / cnt).sum())

def chain(qid):
    names, rows, M, n_raw = load_hits(qid)
    qm = np.frombuffer(qseq[qid].encode(), np.uint8)
    eq = G.is_eq.reindex(names).to_numpy(bool) if names else np.zeros(0, bool)
    k_all = greedy95(qm, M)
    k_eq = greedy95(qm, M[eq])
    en, er, EM = load_extra(qid)
    st = {"qid": qid, "query_aa": len(qm), "n_hits_raw": n_raw, "n_hits_qcov50": len(names),
          "n_hits_qcov50_eq": int(eq.sum()), "n_nr95_own": len(k_all) + 1, "n_nr95_eq_only": len(k_eq) + 1,
          "neff80_own": neff80(qm, M[k_all]), "neff80_eq_only": neff80(qm, M[eq][k_eq]),
          "n_genomes_own": int(G.gix.reindex(names).nunique()) if names else 0}
    fam = Q.family[qid]
    members = G.index[G.fc.to_numpy() == FAMC[fam]]
    st["family_members"] = len(members)
    st["family_members_in_alignment"] = int(len(set(members) & set(names)))
    fn = [names[i] for i in k_all]; fr = [rows[i] for i in k_all]; FM = M[k_all]
    st["extra_db_used"] = EM is not None
    if EM is not None and len(EM):
        base = np.vstack([qm[None], FM])
        # filter external rows against query + own rows
        kx = []
        kept = base
        L = len(qm)
        buf = np.empty((len(EM), L), np.uint8); nb = 0
        for i in range(len(EM)):
            r = EM[i]; ng = r != GAP; d = int(ng.sum())
            if ((kept[:, ng] == r[ng]).sum(1) >= 0.95 * d).any(): continue
            if nb and ((buf[:nb][:, ng] == r[ng]).sum(1) >= 0.95 * d).any(): continue
            buf[nb] = r; nb += 1; kx.append(i)
            if nb >= 8000: break
        st["n_extra_qcov50"] = len(EM); st["n_nr95_extra_added"] = len(kx)
        fn += [en[i] for i in kx]; fr += [er[i] for i in kx]; FM = np.vstack([FM, EM[kx]]) if kx else FM
    else:
        st["n_extra_qcov50"] = 0; st["n_nr95_extra_added"] = 0
    st["n_nr95_total"] = len(fn) + 1
    st["neff80_total"] = neff80(qm, FM) if st["n_nr95_extra_added"] else st["neff80_own"]
    return st, (names, rows, M), (fn, fr)

def mono_job(args):
    label, qid = args
    st, _, (fn, fr) = chain(qid)
    q = qseq[qid]
    txt = [f">{label}", q]
    for nm, r in list(zip(fn, fr))[: a.max_mono]:
        txt += [">" + nm, r]
    st["label"] = label; st["n_rows_written"] = len(txt) // 2
    return label, "\n".join(txt) + "\n", st

def pair_job(r):
    qa, qb = r["qid_A"], r["qid_B"]
    sa, (na, ra, MA), (fna, fra) = chain(qa)
    sb, (nb_, rb, MB), (fnb, frb) = chain(qb)
    A, B = qseq[qa], qseq[qb]; LA, LB = len(A), len(B)
    ga = G.loc[na] if na else G.iloc[:0]; gb = G.loc[nb_] if nb_ else G.iloc[:0]
    pairs = []; used = set()
    if r["type"] == "real":
        bpos = {}
        for jb, (c, rk) in enumerate(zip(gb.ctg.to_numpy(), gb["rank"].to_numpy())):
            bpos.setdefault((c, rk), jb)
        for ia, (c, rk, gx) in enumerate(zip(ga.ctg.to_numpy(), ga["rank"].to_numpy(), ga.gix.to_numpy())):
            if gx in used: continue
            cand = [bpos[k] for k in ((c, rk - 1), (c, rk + 1)) if k in bpos]
            if cand:
                pairs.append((ia, min(cand))); used.add(gx)
    else:
        firstb = {}
        for jb, gx in enumerate(gb.gix.to_numpy()): firstb.setdefault(gx, []).append(jb)
        for ia, gx in enumerate(ga.gix.to_numpy()):
            if gx in used or gx not in firstb: continue
            jb = next((j for j in firstb[gx] if nb_[j] != na[ia]), None)
            if jb is None: continue
            pairs.append((ia, jb)); used.add(gx)
    qm = np.frombuffer((A + B).encode(), np.uint8)
    if pairs:
        PM = np.hstack([MA[[p[0] for p in pairs]], MB[[p[1] for p in pairs]]])
        kp = greedy95(qm, PM, cap=int(r.get("_cap", a.max_paired)))
        neffp = neff80(qm, PM[kp])
    else:
        kp = np.array([], int); neffp = 1.0
    eqp = sum(bool(G.is_eq[na[pairs[i][0]]]) for i in kp)
    lines = [f"#{LA},{LB}\t1,1", ">101\t102", A + B]
    for i in kp:
        ia, jb = pairs[i]
        lines += [f">{na[ia]}\t{nb_[jb]}", ra[ia] + rb[jb]]
    lines += [">101", A + "-" * LB]
    for nm, row in list(zip(fna, fra))[: a.max_unpaired]:
        lines += [">" + nm, row + "-" * LB]
    lines += [">102", "-" * LA + B]
    for nm, row in list(zip(fnb, frb))[: a.max_unpaired]:
        lines += [">" + nm, "-" * LA + row]
    st = {"pair_id": r["pair_id"], "type": r["type"], "len_A": LA, "len_B": LB, "n_genomes_paired_raw": len(pairs), "paired_cap": int(r.get("_cap", a.max_paired)),
          "n_paired_nr95": int(len(kp)) + 1, "neff80_paired": neffp, "n_paired_nr95_from_eq": int(eqp),
          "n_unpaired_A": min(len(fna), a.max_unpaired) + 1, "n_unpaired_B": min(len(fnb), a.max_unpaired) + 1,
          "neff80_A": sa["neff80_total"], "neff80_B": sb["neff80_total"],
          "neff80_A_own": sa["neff80_own"], "neff80_B_own": sb["neff80_own"]}
    return r["pair_id"], "\n".join(lines) + "\n", st

def write_tar(path, items):
    with tarfile.open(path, "w") as tf:
        for name, txt in items:
            b = txt.encode(); ti = tarfile.TarInfo(name); ti.size = len(b); ti.mtime = 0
            tf.addfile(ti, io.BytesIO(b))

if __name__ == "__main__":
    reps = REPS[REPS.usable & REPS.qid.notna()].reset_index(drop=True)
    reps["label"] = [f"M{i + 1:04d}" for i in range(len(reps))]
    with Pool(a.threads) as pool:
        mono = pool.map(mono_job, list(zip(reps.label, reps.qid)), chunksize=4)
        prs = pool.map(pair_job, PAIRS[PAIRS.type == "real"].to_dict("records"), chunksize=2)
        # controls: paired depth capped at the paired depth of the matched real pair (depth matching)
        depth = {l: s_["n_paired_nr95"] - 1 for l, _, s_ in prs}
        crec = PAIRS[PAIRS.type == "control"].to_dict("records")
        for r_ in crec:
            r_["_cap"] = max(int(depth.get(r_["matched_real_pair"], a.max_paired)), 0)
        prs += pool.map(pair_job, crec, chunksize=2)
    ms = pd.DataFrame([m[2] for m in mono]).merge(
        reps[["label", "family_rep", "in_981", "rep_gene", "tier", "signal_peptide_cut", "mature_aa"]], on="label")
    ms.to_csv(a.out / f"monomer_alignment_depth_{a.tag}.tsv", sep="\t", index=False)
    write_tar(a.out / f"monomer_a3m_{a.tag}.tar", [(f"monomer/{l}.a3m", t) for l, t, _ in mono])
    ps = pd.DataFrame([p[2] for p in prs]).merge(PAIRS, on=["pair_id", "type"], suffixes=("", "_step1"))
    ps.to_csv(a.out / f"pair_alignment_depth_{a.tag}.tsv", sep="\t", index=False)
    write_tar(a.out / f"complex_real_a3m_{a.tag}.tar", [(f"real/{l}.a3m", t) for l, t, s in prs if s["type"] == "real"])
    write_tar(a.out / f"complex_control_a3m_{a.tag}.tar",
              [(f"control/{l}.a3m", t) for l, t, s in prs if s["type"] == "control"])
    f = ms[ms.in_981]
    summ = {"tag": a.tag, "n_monomer_alignments": int(len(ms)), "n_981": int(len(f)),
            "median_nr95_eq_only": float(f.n_nr95_eq_only.median()), "median_nr95_own": float(f.n_nr95_own.median()),
            "median_nr95_total": float(f.n_nr95_total.median()),
            "median_neff_eq_only": float(f.neff80_eq_only.median()), "median_neff_own": float(f.neff80_own.median()),
            "median_neff_total": float(f.neff80_total.median()),
            "n_981_neff_ge30_eq_only": int((f.neff80_eq_only >= 30).sum()),
            "n_981_neff_ge30_own": int((f.neff80_own >= 30).sum()),
            "n_981_neff_ge30_total": int((f.neff80_total >= 30).sum()),
            "n_981_neff_ge100_eq_only": int((f.neff80_eq_only >= 100).sum()),
            "n_981_neff_ge100_own": int((f.neff80_own >= 100).sum()),
            "n_981_neff_ge100_total": int((f.neff80_total >= 100).sum()),
            "median_family_member_recall": float((f.family_members_in_alignment / f.family_members).median()),
            "pairs": ps.groupby("type").agg(n=("pair_id", "size"), median_paired=("n_paired_nr95", "median"),
                                            median_neff_paired=("neff80_paired", "median"),
                                            median_len_B=("len_B", "median")).to_dict("index")}
    (a.out / f"alignment_summary_{a.tag}.json").write_text(json.dumps(summ, indent=1))
    print("ALIGN_DONE", json.dumps(summ))
