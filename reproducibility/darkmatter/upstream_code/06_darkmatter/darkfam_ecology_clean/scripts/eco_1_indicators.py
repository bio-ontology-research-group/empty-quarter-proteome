#!/usr/bin/env python3
"""eco_1 — per-family indicators vs environment / XRF (§2.5 rerun on the 106,042 universe).

Copy of paper-repo code/phase_e_indicators.py (June 2026, 23,382 object) repointed to the
rebuilt basis. Method verbatim: per family, Spearman(abundance, env_var) across samples that
carry the variable (XRF-gated / climate subsets); point-biserial for the binary plant-root
contrast; BH-FDR per axis; family tested only if nonzero in >= MIN_PRESENT subset samples.

Deviations from the June script (each forced, documented in ../README.md):
  - EC-decode join dropped: no family_ec.tsv exists for the 106,042 basis; per-protein
    decoding is deferred to the sibling annotator paper (locked structure).
  - Plant-host genus axes dropped: plants_by_site_trip.tsv (amplicon cache) is not on IBEX.
    The plant-root compartment contrast (from env_matrix) is retained.
  - meta join = fam_recurrence.tsv (n_members/n_samples); dnds join = selection/dnds.tsv.
  - Adds results/indicators_summary.tsv (per-axis machine-readable counts + direction split);
    the June script only printed these to stdout.
  - Adds results/indicators/<axis>_all.tsv: the FULL per-family table for every axis
    (family_id, spearman_r, p, q, n_present), so downstream steps join to the persisted
    statistics instead of recomputing them. The <axis>_top.tsv files are unchanged.
"""
import os, csv
import numpy as np
from scipy.stats import rankdata, t as tdist

B = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter"
REPO = "/ibex/project/c2014/alelopezv/repos/empty-quarter-mags"
ENVF = os.path.join(REPO, "results", "darkmatter_ecology", "env_matrix.tsv")
ABUND = os.path.join(B, "darkfam_basis_clean", "results", "fam_abundance.tsv")
META = os.path.join(B, "darkfam_basis_clean", "results", "fam_recurrence.tsv")
DNDS = os.path.join(B, "darkfam_basis_clean", "results", "selection", "dnds.tsv")
RES = os.path.join(B, "darkfam_ecology_clean", "results")
OUT = os.path.join(RES, "indicators")
os.makedirs(OUT, exist_ok=True)

MIN_PRESENT = 25      # family must be nonzero in >= this many subset samples to be tested
TOP_N = 40            # rows written per axis

# XRF elements worth testing (V = nitrogenase hypothesis; Mo = its competitor cofactor;
# Ni/Fe/Cu/Co = metalloenzyme cofactors; S/Cl/Na = salinity chemistry)
XRF_ELEMENTS = ["V", "Mo", "Ni", "Fe", "Cu", "Co", "Mn", "Zn", "S", "Cl", "Na", "Ca", "K"]


def load_env(path):
    return {r["sample"]: r for r in csv.DictReader(open(path), delimiter="\t")}


def load_matrix(path):
    with open(path) as fh:
        samples = fh.readline().rstrip("\n").split("\t")[1:]
        fams, rows = [], []
        for line in fh:
            c = line.rstrip("\n").split("\t")
            fams.append(c[0])
            rows.append([float(x) for x in c[1:]])
    return samples, fams, np.asarray(rows, dtype=np.float64)   # families x samples


def load_kv(path, key, *vals):
    out = {}
    with open(path) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            out[row[key]] = tuple(row.get(v, "") for v in vals)
    return out


def spearman_vec(M, y):
    """Spearman r of each row of M (n_fam x n) against vector y (len n). Returns (r, p)."""
    n = len(y)
    ry = rankdata(y)
    ry = ry - ry.mean()
    rM = rankdata(M, axis=1).astype(np.float64)
    rM -= rM.mean(axis=1, keepdims=True)
    num = rM @ ry
    den = np.sqrt((rM**2).sum(axis=1) * (ry**2).sum())
    r = np.divide(num, den, out=np.zeros(M.shape[0]), where=den > 0)
    r = np.clip(r, -0.999999, 0.999999)
    tstat = r * np.sqrt((n - 2) / (1 - r**2))
    p = 2 * tdist.sf(np.abs(tstat), n - 2)
    return r, p


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


def main():
    env = load_env(ENVF)
    samples, fams, M = load_matrix(ABUND)
    print(f"matrix {M.shape[0]} families x {M.shape[1]} samples")

    dnds = load_kv(DNDS, "family_id", "median_dNdS", "mean_dNdS")
    meta = load_kv(META, "family_rep", "n_members", "n_samples")

    def annotate(fam):
        m = meta.get(fam, ("", ""))
        dn = dnds.get(fam, ("",))
        return {"n_members": m[0], "n_samples": m[1], "dnds": dn[0]}

    results = {}

    # ---- continuous axes: XRF elements + salinity (XRF-gated) + aridity (climate) ----
    cont_axes = [(el, "has_xrf") for el in XRF_ELEMENTS] + \
                [("salinity_proxy", "has_xrf"), ("aridity_index", None)]
    for var, gate in cont_axes:
        yvals = []
        for s in samples:
            r = env.get(s)
            if r is None:
                continue
            if gate and r.get(gate) != "1":
                continue
            v = r.get(var, "")
            if v in ("", "NA", None):
                continue
            try:
                yv = float(v)
            except ValueError:
                continue
            yvals.append((s, yv))
        cols = [samples.index(s) for s, _ in yvals]
        y = np.array([v for _, v in yvals])
        if len(y) < 30:
            continue
        sub = M[:, cols]
        present = (sub > 0).sum(axis=1)
        mask = present >= MIN_PRESENT
        r, p = spearman_vec(sub[mask], y)
        q = bh_fdr(p)
        fmask = np.array(fams)[mask]
        rows = sorted(zip(fmask, r, p, q, present[mask]),
                      key=lambda t: (t[3], -abs(t[1])))
        results[var] = (len(y), int(mask.sum()), rows)
        print(f"{var:16s} n={len(y):4d} tested={int(mask.sum()):6d} "
              f"q<0.05: {int((q<0.05).sum()):5d}  q<0.1: {int((q<0.1).sum()):5d}")

    # ---- binary: plant-root compartment vs others (point-biserial = Spearman w/ 0/1) ----
    yc = np.array([1.0 if env.get(s, {}).get("compartment") == "PlantRoot" else 0.0
                   for s in samples])
    present = (M > 0).sum(axis=1)
    mask = present >= MIN_PRESENT
    r, p = spearman_vec(M[mask], yc)
    q = bh_fdr(p)
    fmask = np.array(fams)[mask]
    rows = sorted(zip(fmask, r, p, q, present[mask]), key=lambda t: (t[3], -abs(t[1])))
    results["plantroot_vs_other"] = (int(yc.sum()), int(mask.sum()), rows)
    print(f"{'plantroot':16s} n={int(yc.sum()):4d} tested={int(mask.sum()):6d} "
          f"q<0.05: {int((q<0.05).sum()):5d}")

    # ---- write per-axis tables + machine-readable summary + markdown ----
    for axis, (n, tested, rows) in results.items():
        with open(os.path.join(OUT, f"{axis}_top.tsv"), "w") as out:
            out.write("family_id\tspearman_r\tp\tq\tn_present\tdnds\t"
                      "n_members\tn_samples_catalog\n")
            for fam, rr, pp, qq, pres in rows[:TOP_N]:
                a = annotate(fam)
                out.write(f"{fam}\t{rr:.3f}\t{pp:.2e}\t{qq:.2e}\t{pres}\t{a['dnds']}\t"
                          f"{a['n_members']}\t{a['n_samples']}\n")
        # full table: every tested family, same row order as the top file (q, then |r|)
        with open(os.path.join(OUT, f"{axis}_all.tsv"), "w") as out:
            out.write("family_id\tspearman_r\tp\tq\tn_present\n")
            for fam, rr, pp, qq, pres in rows:
                out.write(f"{fam}\t{rr:.6f}\t{pp:.6e}\t{qq:.6e}\t{pres}\n")

    with open(os.path.join(RES, "indicators_summary.tsv"), "w") as out:
        out.write("axis\tn_samples\tn_tested\tn_q05\tn_q10\tn_q05_negative\tpct_negative_q05\n")
        for axis, (n, tested, rows) in results.items():
            sig = [t for t in rows if t[3] < 0.05]
            sig10 = [t for t in rows if t[3] < 0.10]
            neg = [t for t in sig if t[1] < 0]
            pct = 100.0 * len(neg) / len(sig) if sig else float("nan")
            out.write(f"{axis}\t{n}\t{tested}\t{len(sig)}\t{len(sig10)}\t{len(neg)}\t{pct:.1f}\n")

    write_markdown(results, annotate)
    print("done ->", OUT)


def write_markdown(results, annotate):
    lines = ["# Per-family indicators (§2.5 rerun, 106,042 universe)\n",
             "Spearman of dark-family abundance vs each axis; BH-FDR per axis. "
             "No EC labels: per-protein decoding is deferred to the annotator sibling paper.\n"]
    for axis, (n, tested, rows) in results.items():
        sig = [r for r in rows if r[3] < 0.05]
        lines.append(f"\n## {axis}  (n={n} samples, {tested} families tested, "
                     f"{len(sig)} at q<0.05)\n")
        if not sig:
            lines.append("_No family passes q<0.05 — honest negative for this axis._\n")
            continue
        lines.append("| family | r | q | n_pres | dN/dS | members | samples |")
        lines.append("|---|---|---|---|---|---|---|")
        for fam, rr, pp, qq, pres in sig[:12]:
            a = annotate(fam)
            short = fam.split("::")[-1]
            lines.append(f"| {short} | {rr:+.2f} | {qq:.1e} | {pres} | "
                         f"{a['dnds']} | {a['n_members']} | {a['n_samples']} |")
    open(os.path.join(RES, "INDICATORS.md"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
