#!/usr/bin/env python3
"""Collect per-model metrics from ColabFold outputs, write confident-region crops, and archive models.
usage: collect_metrics.py SET SHARD   (SET: monomer | complex_real | complex_control)
Runs at the end of each folding task (seconds) and writes one metrics table per shard.
Confidence scaled 0-1. Crops (paper rule): centred 15-residue moving average of pLDDT (window shortened at
chain ends) >= 0.70; gaps <= 5 residues bridged when their mean raw pLDDT >= 0.55; intervals >= 50 residues."""
import sys, json, re, glob, os, tarfile
from pathlib import Path
import numpy as np
D = Path("/mnt/data/u/hohndor/eq-exported-dark-20261111")
SET = sys.argv[1]; SHARD = sys.argv[2]
MONO = SET.startswith("monomer")   # monomer_tail: long monomers offloaded to extra GPUs
P = D / "predictions" / SET
out = D / "results"; out.mkdir(exist_ok=True)

def read_pdb(path):
    ch, rs, xyz, bf, ca, lines = [], [], [], [], [], []
    for l in open(path):
        if l.startswith("ATOM"):
            lines.append(l)
            ch.append(l[21]); rs.append(int(l[22:26])); xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            bf.append(float(l[60:66])); ca.append(l[12:16].strip() == "CA")
    return np.array(ch), np.array(rs), np.array(xyz), np.array(bf), np.array(ca), lines

def regions(pl):
    n = len(pl); sm = np.empty(n)
    for i in range(n):
        sm[i] = pl[max(0, i - 7):min(n, i + 8)].mean()
    ok = sm >= 0.70
    iv = []; i = 0
    while i < n:
        if ok[i]:
            j = i
            while j + 1 < n and ok[j + 1]: j += 1
            iv.append([i, j]); i = j + 1
        else: i += 1
    merged = []
    for s, e in iv:
        if merged and s - merged[-1][1] - 1 <= 5 and pl[merged[-1][1] + 1:s].mean() >= 0.55:
            merged[-1][1] = e
        else: merged.append([s, e])
    return [(s, e) for s, e in merged if e - s + 1 >= 50]

rows = []; crops_dir = D / "foldseek" / "queries" / "monomer_crops"; whole_dir = D / "foldseek" / "queries" / "monomer_whole"
if MONO:
    crops_dir.mkdir(parents=True, exist_ok=True); whole_dir.mkdir(parents=True, exist_ok=True)
crop_rows = []
for done in sorted(P.glob(f"{SHARD}/*.done.txt")):
    name = done.name[:-len(".done.txt")]; sd = done.parent
    sc = {}
    for js in sd.glob(f"{name}_scores_rank_*.json"):
        m = re.search(r"_scores_rank_(\d+)_(.+)_model_(\d)_seed_(\d+)\.json$", js.name)
        j = json.load(open(js))
        sc[int(m.group(1))] = (int(m.group(3)), j, js)
    if 1 not in sc:
        rows.append({"id": name, "status": "no_scores"}); continue
    pdbs = list(sd.glob(f"{name}_unrelaxed_rank_001_*.pdb"))
    ch, rs, xyz, bf, ca, lines = read_pdb(pdbs[0])
    pl = bf[ca] / 100.0
    r = {"id": name, "status": "ok", "shard": sd.name, "n_res": int(ca.sum()), "top_model": sc[1][0],
         "mean_plddt": float(pl.mean()), "frac_res_plddt_ge70": float((pl >= 0.70).mean()),
         "ptm": float(sc[1][1].get("ptm", np.nan)), "n_models": len(sc),
         "mean_plddt_all_models": ";".join(f"{np.mean(sc[k][1]['plddt']) / 100:.4f}" for k in sorted(sc)),
         "ptm_all_models": ";".join(f"{sc[k][1].get('ptm', float('nan')):.4f}" for k in sorted(sc)),
         "model_order_by_rank": ";".join(str(sc[k][0]) for k in sorted(sc)), "pdb": pdbs[0].name}
    m1 = [k for k in sc if sc[k][0] == 1]
    if m1: r["mean_plddt_model_1"] = float(np.mean(sc[m1[0]][1]["plddt"]) / 100)
    if MONO:
        r["qualified"] = bool(r["n_res"] >= 100 and r["mean_plddt"] >= 0.70)
        (whole_dir / f"{name}.pdb").write_text("".join(lines) + "END\n")
        resid = rs[ca]
        regs = regions(pl)
        r["n_confident_regions"] = len(regs); r["confident_region_residues"] = int(sum(e - s + 1 for s, e in regs))
        for k, (s, e) in enumerate(regs, 1):
            keep = set(resid[s:e + 1].tolist())
            cl = [l for l in lines if int(l[22:26]) in keep]
            cid = f"{name}__crop{k}_{resid[s]}_{resid[e]}"
            (crops_dir / f"{cid}.pdb").write_text("".join(cl) + "END\n")
            crop_rows.append({"crop_id": cid, "id": name, "start": int(resid[s]), "end": int(resid[e]),
                              "n_res": e - s + 1, "mean_plddt": float(pl[s:e + 1].mean())})
    else:
        j = sc[1][1]
        chains = sorted(set(ch)); A = ch == chains[0]; B = ch == chains[1]
        caA = ca & A; caB = ca & B
        LA, LB = int(caA.sum()), int(caB.sum())
        pae = np.array(j["pae"])
        # heavy-atom contacts < 5 A between chains -> interface residues
        XA, XB = xyz[A], xyz[B]
        d2 = ((XA[:, None, :] - XB[None, :, :]) ** 2).sum(-1) if len(XA) * len(XB) < 6e7 else None
        if d2 is None:
            close = np.zeros((len(XA), len(XB)), bool)
            for s in range(0, len(XA), 500):
                close[s:s + 500] = ((XA[s:s + 500, None, :] - XB[None, :, :]) ** 2).sum(-1) < 25.0
        else: close = d2 < 25.0
        ra = np.unique(rs[A][close.any(1)]); rb = np.unique(rs[B][close.any(0)])
        residA = rs[caA]; residB = rs[caB]
        ia = np.flatnonzero(np.isin(residA, ra)); ib = np.flatnonzero(np.isin(residB, rb))
        plA, plB = bf[caA] / 100, bf[caB] / 100
        # residue-level contact map on CB-like (CA) < 8 A
        CA_A, CA_B = xyz[caA], xyz[caB]
        dca = np.sqrt(((CA_A[:, None, :] - CA_B[None, :, :]) ** 2).sum(-1))
        cont = dca < 8.0
        inter = np.concatenate([pae[:LA, LA:].ravel(), pae[LA:, :LA].T.ravel()])
        r.update(len_A=LA, len_B=LB, iptm=float(j.get("iptm", np.nan)),
                 ranking_confidence=float(0.8 * j.get("iptm", np.nan) + 0.2 * j.get("ptm", np.nan)),
                 iptm_all_models=";".join(f"{sc[k][1].get('iptm', float('nan')):.4f}" for k in sorted(sc)),
                 iptm_mean_models=float(np.mean([sc[k][1].get("iptm", np.nan) for k in sc])),
                 iptm_min_models=float(np.min([sc[k][1].get("iptm", np.nan) for k in sc])),
                 n_models_iptm_ge05=int(sum(sc[k][1].get("iptm", 0) >= 0.5 for k in sc)),
                 n_interface_res_A=int(len(ia)), n_interface_res_B=int(len(ib)),
                 n_ca_contacts_8A=int(cont.sum()),
                 interface_plddt=float(np.concatenate([plA[ia], plB[ib]]).mean()) if len(ia) + len(ib) else np.nan,
                 interface_pae=float(((pae[:LA, LA:][cont].mean() + pae[LA:, :LA].T[cont].mean()) / 2)) if cont.any() else np.nan,
                 mean_interchain_pae=float(inter.mean()), min_interchain_pae=float(inter.min()),
                 mean_plddt_A=float(plA.mean()), mean_plddt_B=float(plB.mean()))
    rows.append(r)
import csv
keys = []
for r in rows:
    for k in r:
        if k not in keys: keys.append(k)
with open(out / f"{SET}_metrics.{SHARD}.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, keys, delimiter="\t"); w.writeheader(); w.writerows(rows)
if MONO:
    with open(out / f"{SET}_crops.{SHARD}.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, ["crop_id", "id", "start", "end", "n_res", "mean_plddt"], delimiter="\t")
        w.writeheader(); w.writerows(crop_rows)
print(SET, "models:", len(rows), "ok:", sum(r["status"] == "ok" for r in rows),
      "qualified:" if MONO else "", sum(r.get("qualified", False) for r in rows) if MONO else "")
