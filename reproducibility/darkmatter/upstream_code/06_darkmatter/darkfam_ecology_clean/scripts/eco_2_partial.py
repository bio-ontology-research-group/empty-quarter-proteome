#!/usr/bin/env python3
"""eco_2 — space-partial per-family associations (§2.5 rerun on the 106,042 universe).

Part (a) of paper-repo code/phase_e_partial_osmolyte.py (June 2026), repointed to the rebuilt
basis. Method verbatim: for the 5 master axes, recompute the per-family associations as
partial Spearman controlling ranked lat & lon (rank-transform, residualize via hat-matrix
projection, correlate residuals), vs the marginal Spearman; report how many BH q<0.05 hits
survive the space control.

Part (b) of the June script (osmolyte-phosphatase vignette) is NOT ported: it is
EC-decode-dependent and per-protein decoding is deferred to the annotator sibling paper.

Added here: the binary plant-root contrast is partialled on the same footing as the five
continuous axes (point-biserial = Spearman on the 0/1 compartment indicator). Its subset is
the 547 samples that carry lat/lon, not the 551 columns of the abundance matrix eco_1 uses
for the marginal axis, so its marginal count is the one recomputed on this subset.
"""
import os, csv
import numpy as np
from scipy.stats import rankdata, t as tdist

B = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter"
REPO = "/ibex/project/c2014/alelopezv/repos/empty-quarter-mags"
ENVF = os.path.join(REPO, "results", "darkmatter_ecology", "env_matrix.tsv")
ABUND = os.path.join(B, "darkfam_basis_clean", "results", "fam_abundance.tsv")
DNDS = os.path.join(B, "darkfam_basis_clean", "results", "selection", "dnds.tsv")
RES = os.path.join(B, "darkfam_ecology_clean", "results")
OUT = os.path.join(RES, "indicators")
os.makedirs(OUT, exist_ok=True)
MIN_PRESENT = 25


def load_env(p):
    return {r["sample"]: r for r in csv.DictReader(open(p), delimiter="\t")}


def load_matrix(p):
    with open(p) as fh:
        samples = fh.readline().rstrip("\n").split("\t")[1:]
        fams, rows = [], []
        for line in fh:
            c = line.rstrip("\n").split("\t")
            fams.append(c[0]); rows.append([float(x) for x in c[1:]])
    return samples, fams, np.asarray(rows)


def load_kv(p, key, *vals):
    out = {}
    for row in csv.DictReader(open(p), delimiter="\t"):
        out[row[key]] = tuple(row.get(v, "") for v in vals)
    return out


def bh_fdr(p):
    p = np.asarray(p, float); n = len(p); order = np.argsort(p); q = np.empty(n); prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]; prev = min(prev, p[idx] * n / (i + 1)); q[idx] = prev
    return q


def ranks(a, axis=None):
    return rankdata(a, axis=axis).astype(float)


def partial_spearman(Msub, y, covars):
    """Partial Spearman of each row of Msub vs y, controlling covars (list of vectors).
    Rank-transform all, residualize on [1,covars] via projection, correlate residuals."""
    n = len(y)
    X = np.column_stack([np.ones(n)] + [ranks(c) for c in covars])
    P = X @ np.linalg.pinv(X.T @ X) @ X.T            # n x n hat matrix
    ry = ranks(y); ey = ry - P @ ry
    rM = ranks(Msub, axis=1); eM = rM - rM @ P.T     # residualize each row
    num = eM @ ey
    den = np.sqrt((eM**2).sum(axis=1) * (ey**2).sum())
    r = np.divide(num, den, out=np.zeros(Msub.shape[0]), where=den > 0)
    r = np.clip(r, -0.999999, 0.999999)
    df = n - 2 - len(covars)
    tstat = r * np.sqrt(df / (1 - r**2))
    p = 2 * tdist.sf(np.abs(tstat), df)
    return r, p


def marginal_spearman(Msub, y):
    return partial_spearman(Msub, y, [])


def axis_value(row, var):
    """Axis value for one env row: the binary plant-root contrast is the 0/1 compartment
    indicator (as in eco_1), every other axis is the numeric column of that name."""
    if var == "plantroot_vs_other":
        return 1.0 if row.get("compartment") == "PlantRoot" else 0.0
    return float(row.get(var, ""))


def subset(samples, env, var, gate):
    cols, y, lat, lon = [], [], [], []
    for j, s in enumerate(samples):
        r = env.get(s)
        if not r or (gate and r.get(gate) != "1"):
            continue
        try:
            v = axis_value(r, var); la = float(r["lat"]); lo = float(r["lon"])
        except (ValueError, KeyError):
            continue
        cols.append(j); y.append(v); lat.append(la); lon.append(lo)
    return np.array(cols), np.array(y), np.array(lat), np.array(lon)


def main():
    env = load_env(ENVF)
    samples, fams, M = load_matrix(ABUND)
    famarr = np.array(fams)
    dnds = load_kv(DNDS, "family_id", "median_dNdS")

    print("=== space-partial (control lat+lon) ===")
    part_lines = ["# Space-partial per-family associations (control lat+lon)\n",
                  "Marginal vs partial-Spearman significant-family counts (BH q<0.05), "
                  "106,042 universe.\n",
                  "| axis | n | marginal q<0.05 | partial q<0.05 | retained % | median |Dr| |",
                  "|---|---|---|---|---|---|"]
    summary_rows = []
    for var, gate in [("salinity_proxy", "has_xrf"), ("S", "has_xrf"), ("Cl", "has_xrf"),
                      ("Ca", "has_xrf"), ("aridity_index", None),
                      ("plantroot_vs_other", None)]:
        cols, y, lat, lon = subset(samples, env, var, gate)
        if len(y) < 30:
            continue
        sub = M[:, cols]; present = (sub > 0).sum(axis=1); mask = present >= MIN_PRESENT
        rm, pm = marginal_spearman(sub[mask], y); qm = bh_fdr(pm)
        rp, pp = partial_spearman(sub[mask], y, [lat, lon]); qp = bh_fdr(pp)
        nm, npar = int((qm < .05).sum()), int((qp < .05).sum())
        both = (qm < .05) & (qp < .05)
        ret = 100 * both.sum() / max(nm, 1)
        dr = float(np.median(np.abs(rm - rp)))
        print(f"  {var:15s} n={len(y):4d} marg={nm:5d} part={npar:5d} retained={ret:4.0f}%")
        part_lines.append(f"| {var} | {len(y)} | {nm} | {npar} | {ret:.0f}% | {dr:.3f} |")
        summary_rows.append((var, len(y), nm, npar, f"{ret:.1f}", f"{dr:.4f}"))
        fm = famarr[mask]
        rows = sorted(zip(fm, rp, qp, rm, qm, present[mask]), key=lambda t: t[2])
        with open(os.path.join(OUT, f"partial_{var}_top.tsv"), "w") as o:
            o.write("family_id\tpartial_r\tpartial_q\tmarginal_r\tmarginal_q\tn_present\tdnds\n")
            for fam, rpr, qpr, rmr, qmr, pres in rows[:60]:
                d = dnds.get(fam, ("",))[0]
                o.write(f"{fam}\t{rpr:.3f}\t{qpr:.2e}\t{rmr:.3f}\t{qmr:.2e}\t{pres}\t{d}\n")
    open(os.path.join(RES, "PARTIAL.md"), "w").write("\n".join(part_lines) + "\n")
    with open(os.path.join(RES, "partial_summary.tsv"), "w") as o:
        o.write("axis\tn_samples\tn_marginal_q05\tn_partial_q05\tretained_pct\tmedian_abs_dr\n")
        for row in summary_rows:
            o.write("\t".join(str(x) for x in row) + "\n")
    print("done ->", OUT)


if __name__ == "__main__":
    main()
