#!/usr/bin/env python3
"""RQ-B3 — unified gene-environment association map for dark protein families.

Synthesises the per-variable indicator scans (phase_e_indicators.py) into one
gene-environment map. For each geochemical/climatic axis we recompute the FULL
per-family Spearman association from the raw dark-family abundance matrix + env
matrix (the indicators/<VAR>_top.tsv files only store the top-40 hits, so the
full positive/negative split and the q<0.10 tier cannot be read off them) and,
for the desert master variables, the space-partial association controlling for
ranked lat+lon (phase_e_partial_osmolyte.py).

Method is a faithful re-implementation of the two phase_e scripts so this run
both VERIFIES their documented q<0.05 counts and EXTENDS them with the pos/neg
split, the q<0.10 tier and the decoded-EC fraction.

Outputs (under RQ/B3/):
  gene_env_summary.tsv   one row per variable
  fig_gene_env.pdf/.png  panel a = #sig families/variable (pos vs neg)
                         panel b = strongest family per key variable
Run as:  python3 RQ/B3/run.py    (from repo root)
Deps: stdlib + numpy + scipy + matplotlib.
"""
import os, csv
import numpy as np
from scipy.stats import rankdata, t as tdist
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
ECO = os.path.join(REPO, "results", "darkmatter_ecology")
OUT = HERE

MIN_PRESENT = 25                       # family must be present in >= this many subset samples
# XRF elements (gated on has_xrf, n=251) then climate axes (salinity gated, aridity ungated)
XRF_ELEMENTS = ["V", "Mo", "Ni", "Fe", "Cu", "Co", "Mn", "Zn", "S", "Cl", "Na", "Ca", "K"]
CONT_AXES = [(el, "has_xrf") for el in XRF_ELEMENTS] + \
            [("salinity_proxy", "has_xrf"), ("aridity_index", None)]
PARTIAL_AXES = ["salinity_proxy", "S", "Cl", "Ca", "aridity_index"]
# variable display order requested by RQ-B3
ORDER = ["V", "Mo", "Ni", "Na", "Cl", "Ca", "K", "Mn", "Cu", "Co", "Zn", "Fe",
         "S", "salinity_proxy", "aridity_index"]
# curated osmolyte-relevant EC numbers (desiccation/salinity solutes)
OSMO_EC = {"3.1.3.69", "3.1.3.70", "3.1.3.12", "2.4.1.15", "2.4.1.245",
           "2.7.7.71", "3.2.1.28", "2.4.1.36"}


# ---------- IO ----------
def load_env(p):
    return {r["sample"]: r for r in csv.DictReader(open(p), delimiter="\t")}


def load_matrix(p):
    with open(p) as fh:
        samples = fh.readline().rstrip("\n").split("\t")[1:]
        fams, rows = [], []
        for line in fh:
            c = line.rstrip("\n").split("\t")
            fams.append(c[0])
            rows.append([float(x) for x in c[1:]])
    return samples, fams, np.asarray(rows, dtype=np.float64)


def load_kv(p, key, *vals):
    out = {}
    for row in csv.DictReader(open(p), delimiter="\t"):
        out[row[key]] = tuple(row.get(v, "") for v in vals)
    return out


# ---------- stats (mirrors phase_e_indicators.py / phase_e_partial_osmolyte.py) ----------
def bh_fdr(p):
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        prev = min(prev, p[idx] * n / (i + 1))
        q[idx] = prev
    return q


def partial_spearman(M, y, covars):
    """Partial Spearman of each row of M vs y, controlling covars (list of vectors).
    covars=[] gives the marginal Spearman."""
    n = len(y)
    X = np.column_stack([np.ones(n)] + [rankdata(c).astype(float) for c in covars])
    P = X @ np.linalg.pinv(X.T @ X) @ X.T
    ry = rankdata(y).astype(float)
    ey = ry - P @ ry
    rM = rankdata(M, axis=1).astype(float)
    eM = rM - rM @ P.T
    num = eM @ ey
    den = np.sqrt((eM ** 2).sum(axis=1) * (ey ** 2).sum())
    r = np.divide(num, den, out=np.zeros(M.shape[0]), where=den > 0)
    r = np.clip(r, -0.999999, 0.999999)
    df = n - 2 - len(covars)
    tstat = r * np.sqrt(df / (1 - r ** 2))
    p = 2 * tdist.sf(np.abs(tstat), df)
    return r, p


def subset(samples, env, var, gate):
    cols, y, lat, lon = [], [], [], []
    for j, s in enumerate(samples):
        r = env.get(s)
        if not r or (gate and r.get(gate) != "1"):
            continue
        v = r.get(var, "")
        if v in ("", "NA", None):
            continue
        try:
            yv = float(v); la = float(r["lat"]); lo = float(r["lon"])
        except (ValueError, KeyError):
            continue
        cols.append(j); y.append(yv); lat.append(la); lon.append(lo)
    return np.array(cols, int), np.array(y), np.array(lat), np.array(lon)


def main():
    env = load_env(os.path.join(ECO, "env_matrix.tsv"))
    samples, fams, M = load_matrix(os.path.join(ECO, "darkfam_abundance.tsv"))
    fams = np.array(fams)
    ec = load_kv(os.path.join(ECO, "family_ec.tsv"), "family_id", "best_ec", "rxn_name")
    dnds = load_kv(os.path.join(ECO, "dnds.tsv"), "family_id", "median_dNdS")
    print(f"matrix {M.shape[0]} families x {M.shape[1]} samples")

    def label(fam):
        e = ec.get(fam, ("", ""))
        lab = e[1] or e[0] or "(EC-less)"
        return e[0], lab

    # ---- full marginal association per axis ----
    res = {}
    for var, gate in CONT_AXES:
        cols, y, lat, lon = subset(samples, env, var, gate)
        if len(y) < 30:
            continue
        sub = M[:, cols]
        present = (sub > 0).sum(axis=1)
        mask = present >= MIN_PRESENT
        fm = fams[mask]
        r, p = partial_spearman(sub[mask], y, [])
        q = bh_fdr(p)
        sig = q < 0.05
        sig10 = q < 0.10
        pos = sig & (r > 0)
        neg = sig & (r < 0)
        # decoded-EC fraction among significant
        ecflag = np.array([bool(ec.get(f, ("", ""))[0]) for f in fm])
        osmo = np.array([ec.get(f, ("", ""))[0] in OSMO_EC for f in fm])
        # top positive / negative family (most significant within sign)
        def top(selmask):
            idx = np.where(selmask)[0]
            if idx.size == 0:
                return ("", "", "", "")
            best = idx[np.argmin(q[idx])]
            ecnum, lab = label(fm[best])
            return (fm[best], f"{r[best]:+.3f}", f"{q[best]:.2e}", lab)
        res[var] = dict(n_test=int(mask.sum()), n=len(y),
                        q05=int(sig.sum()), q10=int(sig10.sum()),
                        npos=int(pos.sum()), nneg=int(neg.sum()),
                        frac_ec=float(ecflag[sig].mean()) if sig.sum() else 0.0,
                        n_osmo=int(osmo[sig].sum()),
                        toppos=top(pos), topneg=top(neg))
        print(f"{var:16s} n={len(y):4d} tested={int(mask.sum()):6d} "
              f"q05={int(sig.sum()):5d} (+{int(pos.sum())}/-{int(neg.sum())}) "
              f"q10={int(sig10.sum()):5d}")

    # ---- space-partial for the desert master variables ----
    partial = {}
    for var in PARTIAL_AXES:
        gate = None if var == "aridity_index" else "has_xrf"
        cols, y, lat, lon = subset(samples, env, var, gate)
        sub = M[:, cols]
        present = (sub > 0).sum(axis=1)
        mask = present >= MIN_PRESENT
        rm, pm = partial_spearman(sub[mask], y, [])
        qm = bh_fdr(pm)
        rp, pp = partial_spearman(sub[mask], y, [lat, lon])
        qp = bh_fdr(pp)
        sm, sp = qm < .05, qp < .05
        nm, npar = int(sm.sum()), int(sp.sum())
        both = int((sm & sp).sum())                       # same family sig in both
        partial[var] = dict(marg=nm, part=npar, both=both,
                            retain=(100.0 * both / nm) if nm else 0.0)
        print(f"partial {var:16s} marginal={nm:5d} partial={npar:5d} both={both:5d} "
              f"retained={partial[var]['retain']:.0f}%")

    # ---- write summary tsv ----
    sp = os.path.join(OUT, "gene_env_summary.tsv")
    with open(sp, "w") as out:
        out.write("variable\tn_samples\tn_tested\tn_assoc_q05\tn_assoc_q10\tn_pos\tn_neg\t"
                  "frac_with_ec\tn_osmolyte_ec\t"
                  "top_pos_family\ttop_pos_r\ttop_pos_q\ttop_pos_label\t"
                  "top_neg_family\ttop_neg_r\ttop_neg_q\ttop_neg_label\t"
                  "marginal_q05\tpartial_q05\tpartial_retained_pct\n")
        for v in ORDER:
            d = res.get(v)
            if d is None:
                continue
            pa = partial.get(v, {})
            tp, tn = d["toppos"], d["topneg"]
            out.write("\t".join(str(x) for x in [
                v, d["n"], d["n_test"], d["q05"], d["q10"], d["npos"], d["nneg"],
                f"{d['frac_ec']:.3f}", d["n_osmo"],
                tp[0], tp[1], tp[2], tp[3],
                tn[0], tn[1], tn[2], tn[3],
                pa.get("marg", ""), pa.get("part", ""),
                (f"{pa['retain']:.0f}" if pa else "")]) + "\n")
    print("wrote", sp)

    make_figure(res, partial)
    return res, partial


# ---------- figure ----------
def make_figure(res, partial):
    order = [v for v in ORDER if v in res]
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13, 6.2))

    # panel a: stacked pos/neg significant counts per variable, sorted by total
    so = sorted(order, key=lambda v: res[v]["q05"])
    yp = np.arange(len(so))
    pos = [res[v]["npos"] for v in so]
    neg = [res[v]["nneg"] for v in so]
    axA.barh(yp, pos, color="#c0392b", label="positive (r>0)")
    axA.barh(yp, neg, left=pos, color="#2471a3", label="negative (r<0)")
    axA.set_yticks(yp)
    axA.set_yticklabels(so, fontsize=9)
    axA.set_xlabel("dark families significantly associated (BH q<0.05)")
    axA.set_title("a  Gene-environment recruitment per axis", loc="left", fontweight="bold")
    axA.legend(fontsize=8, loc="lower right", frameon=False)
    for i, v in enumerate(so):
        t = res[v]["q05"]
        if t:
            axA.text(t + max(pos[i] + neg[i] for i in [i]) * 0.0 + 40, i, str(t),
                     va="center", fontsize=7.5, color="#444")
    axA.set_xlim(0, max(res[v]["q05"] for v in so) * 1.12 + 50)

    # panel b: strongest single family per key axis (signed r), labelled
    key = [v for v in ["aridity_index", "salinity_proxy", "Ca", "S", "Cl",
                       "Na", "K", "Fe", "Mo"] if v in res]
    rows = []
    for v in key:
        d = res[v]
        # pick the more-extreme of top pos / top neg by |r|
        cands = []
        for tag in ("toppos", "topneg"):
            f, rr, qq, lab = d[tag]
            if f:
                cands.append((float(rr), qq, lab))
        if not cands:
            continue
        rr, qq, lab = max(cands, key=lambda c: abs(c[0]))
        rows.append((v, rr, lab))
    rows.sort(key=lambda t: t[1])
    yb = np.arange(len(rows))
    cols = ["#c0392b" if r > 0 else "#2471a3" for _, r, _ in rows]
    axB.barh(yb, [r for _, r, _ in rows], color=cols)
    axB.set_yticks(yb)
    axB.set_yticklabels([v for v, _, _ in rows], fontsize=9)
    axB.axvline(0, color="k", lw=.8)
    axB.set_xlabel("Spearman r of strongest dark family")
    axB.set_title("b  Strongest association per key axis", loc="left", fontweight="bold")
    # place EC label in the empty half (opposite side of the bar) to avoid collisions
    for i, (v, r, lab) in enumerate(rows):
        if r > 0:
            axB.text(-0.02, i, lab[:36], va="center", ha="right", fontsize=7)
        else:
            axB.text(0.02, i, lab[:36], va="center", ha="left", fontsize=7)
    axB.set_xlim(-0.62, 0.62)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(OUT, f"fig_gene_env.{ext}"), dpi=200, bbox_inches="tight")
    print("wrote", os.path.join(OUT, "fig_gene_env.pdf"))


if __name__ == "__main__":
    main()
