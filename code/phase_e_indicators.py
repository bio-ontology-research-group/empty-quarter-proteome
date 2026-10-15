#!/usr/bin/env python3
"""Phase E+ — per-family functional indicators vs environment / XRF / plant host.

Phase E asked whether *categories* of decoded function are enriched among structured
families (answer: no). This asks the finer question Phase E never did: which *individual*
multi-member dark families track a specific environmental axis — vanadium and other XRF
elements, salinity, aridity, plant-root compartment, dominant plant host — and what does the
C.V.T decoder say those families are?

Method: per family, Spearman(abundance, env_var) across samples that carry the variable
(XRF 251, climate 456); point-biserial for the binary compartment / plant-host contrasts.
BH-FDR per axis. Join decoded EC (family_ec), selection (dnds), endemism (n_samples_catalog).

Outputs results/darkmatter_ecology/indicators/{axis}_top.tsv + INDICATORS.md.
Self-contained (numpy/scipy/csv); abundance + env already pulled to results/.
"""
import os, csv, math
import numpy as np
from scipy.stats import rankdata, t as tdist

HERE = os.path.dirname(os.path.abspath(__file__))
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")
OUT = os.path.join(ECO, "indicators")
os.makedirs(OUT, exist_ok=True)
AMP = os.path.join(HERE, "..", "..", "empty-quarter-amplicon", "repository", "cache")

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


def dominant_plant_by_site_trip(path):
    """Nearest iNaturalist observation per (site,trip) -> genus token. Prefer species/genus rank."""
    best = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                site, trip = row["site"], row["trip"]
                dist = float(row["distance_km"])
            except (KeyError, ValueError):
                continue
            name = row["scientific_name"].strip()
            rank = row.get("rank", "")
            if rank in ("kingdom", "phylum"):   # uninformative ('Plantae')
                continue
            genus = name.split()[0] if name else ""
            if not genus:
                continue
            k = (site, trip)
            if k not in best or dist < best[k][0]:
                best[k] = (dist, genus)
    return {k: v[1] for k, v in best.items()}


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
    env = load_env(os.path.join(ECO, "env_matrix.tsv"))
    samples, fams, M = load_matrix(os.path.join(ECO, "darkfam_abundance.tsv"))
    famidx = {f: i for i, f in enumerate(fams)}
    print(f"matrix {M.shape[0]} families x {M.shape[1]} samples")

    ec = load_kv(os.path.join(ECO, "family_ec.tsv"), "family_id",
                 "best_ec", "rxn_name", "pathways", "C_confidence")
    dnds = load_kv(os.path.join(ECO, "dnds.tsv"), "family_id", *_dnds_cols())
    meta = load_kv(os.path.join(ECO, "darkfam_meta.tsv"), "family_id",
                   "n_members", "n_samples_catalog")
    plant = dominant_plant_by_site_trip(os.path.join(AMP, "plants_by_site_trip.tsv"))

    def annotate(fam):
        e = ec.get(fam, ("", "", "", ""))
        m = meta.get(fam, ("", ""))
        dn = dnds.get(fam, ("",))
        return {"best_ec": e[0], "rxn": e[1], "pathways": e[2], "C": e[3],
                "n_members": m[0], "n_samples": m[1], "dnds": dn[0]}

    results = {}

    # ---- continuous axes: XRF elements (251) + salinity + aridity (456) ----
    cont_axes = [(el, "has_xrf") for el in XRF_ELEMENTS] + \
                [("salinity_proxy", "has_xrf"), ("aridity_index", None)]
    for var, gate in cont_axes:
        keep, yvals = [], []
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
            keep.append(samples.index(s) if False else None)  # placeholder
            yvals.append((s, yv))
        # build column index + y
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
    cols = list(range(len(samples)))
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

    # ---- plant host: one-vs-rest among PlantRoot samples, top genera ----
    pr_cols, pr_genus = [], []
    for j, s in enumerate(samples):
        r = env.get(s, {})
        if r.get("compartment") != "PlantRoot":
            continue
        g = plant.get((r.get("site"), r.get("trip_num")))
        if g:
            pr_cols.append(j)
            pr_genus.append(g)
    pr_cols = np.array(pr_cols)
    pr_genus = np.array(pr_genus)
    from collections import Counter
    counts = Counter(pr_genus)
    print("plant-root host genera (mapped samples):", dict(counts.most_common(8)))
    subPR = M[:, pr_cols]
    presentPR = (subPR > 0).sum(axis=1)
    maskPR = presentPR >= MIN_PRESENT
    for genus, c in counts.most_common():
        if c < 12 or (len(pr_genus) - c) < 12:
            continue
        yv = (pr_genus == genus).astype(float)
        r, p = spearman_vec(subPR[maskPR], yv)
        q = bh_fdr(p)
        fmask = np.array(fams)[maskPR]
        rows = sorted(zip(fmask, r, p, q, presentPR[maskPR]),
                      key=lambda t: (t[3], -abs(t[1])))
        nsig = int((q < 0.05).sum())
        results[f"host_{genus}"] = (c, int(maskPR.sum()), rows)
        print(f"host_{genus:12s} n={c:4d} tested={int(maskPR.sum()):6d} q<0.05: {nsig:5d}")

    # ---- write per-axis tables + markdown ----
    for axis, (n, tested, rows) in results.items():
        with open(os.path.join(OUT, f"{axis}_top.tsv"), "w") as out:
            out.write("family_id\tspearman_r\tp\tq\tn_present\tbest_ec\trxn_name\t"
                      "pathways\tC_confidence\tdnds\tn_members\tn_samples_catalog\n")
            for fam, rr, pp, qq, pres in rows[:TOP_N]:
                a = annotate(fam)
                out.write(f"{fam}\t{rr:.3f}\t{pp:.2e}\t{qq:.2e}\t{pres}\t{a['best_ec']}\t"
                          f"{a['rxn']}\t{a['pathways']}\t{a['C']}\t{a['dnds']}\t"
                          f"{a['n_members']}\t{a['n_samples']}\n")

    write_markdown(results, annotate)
    print("done ->", OUT)


def _dnds_cols():
    # dnds.tsv: detect the dN/dS column name lazily
    return ("median_dNdS", "mean_dNdS")


def write_markdown(results, annotate):
    lines = ["# Per-family functional indicators (Phase E+)\n",
             "Spearman of dark-family abundance vs each axis; BH-FDR per axis. "
             "Top decoded families with q<0.05 shown. Decoder P is calibrated 0.22-0.43 "
             "(C.V.T) so EC labels are *plausible*, not certain.\n"]
    for axis, (n, tested, rows) in results.items():
        sig = [r for r in rows if r[3] < 0.05]
        lines.append(f"\n## {axis}  (n={n} samples, {tested} families tested, "
                     f"{len(sig)} at q<0.05)\n")
        if not sig:
            lines.append("_No family passes q<0.05 — honest negative for this axis._\n")
            continue
        lines.append("| family | r | q | n_pres | decoded EC | reaction/pathway | dN/dS |")
        lines.append("|---|---|---|---|---|---|---|")
        for fam, rr, pp, qq, pres in sig[:12]:
            a = annotate(fam)
            label = a["rxn"] or a["pathways"] or "(EC only)"
            short = fam.split("::")[-1]
            lines.append(f"| {short} | {rr:+.2f} | {qq:.1e} | {pres} | "
                         f"{a['best_ec']} | {label[:42]} | {a['dnds']} |")
    open(os.path.join(ECO, "INDICATORS.md"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
