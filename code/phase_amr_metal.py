#!/usr/bin/env python3
"""#1 Metal / metalloid resistance census + geochemistry link.

AMRFinderPlus STRESS rows (filtered out of the promoted amr_long.tsv) = metal/metalloid
resistance. Census by metal class, phylum, compartment; then test per-sample resistance load
vs XRF. Caveat: the XRF panel has NO As or Hg (the two biggest classes), so cognate correlation
is only possible for Cu/Ni/Cr/Co/Zn (few genes); for As/Hg we run an exploratory scan vs every
measured element (As co-varies geochemically with Fe/S sulfides).

In: results/amr_metal/stress_long.tsv, env_matrix.tsv, workbench GTDB summaries.
Out: results/amr_metal/{metal_census.tsv, metal_xrf_corr.tsv, METAL.md}
"""
import os, csv, math
from collections import Counter, defaultdict
import numpy as np
from scipy.stats import rankdata, t as tdist

HERE = os.path.dirname(os.path.abspath(__file__))
AMRM = os.path.join(HERE, "..", "results", "amr_metal")
ECO = os.path.join(HERE, "..", "results", "darkmatter_ecology")
GTDB = [os.path.join(HERE, "..", "workbench", "results", "soil_gtdbtk", x)
        for x in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv")]
XRF_ELEMENTS = ["Cu", "Ni", "Cr", "Co", "Zn", "Fe", "Mn", "As", "Hg", "S", "Cl",
                "salinity_proxy", "aridity_index", "V", "Mo"]
COGNATE = {"COPPER": "Cu", "COPPER/SILVER": "Cu", "NICKEL": "Ni",
           "CADMIUM/COBALT/NICKEL": "Ni", "CHROMATE": "Cr", "ZINC": "Zn"}


def load_gtdb():
    tax = {}
    for p in GTDB:
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                tax[r["user_genome"]] = r["classification"]
    return tax


def rank(lin, pre):
    for t in lin.split(";"):
        if t.startswith(pre):
            return t[3:]
    return ""


def spearman(x, y):
    n = len(x)
    if n < 12:
        return float("nan"), float("nan")
    rx, ry = rankdata(x), rankdata(y)
    rx, ry = rx - rx.mean(), ry - ry.mean()
    d = math.sqrt((rx**2).sum() * (ry**2).sum())
    if d == 0:
        return float("nan"), float("nan")
    r = float((rx * ry).sum() / d)
    r = max(min(r, 0.999999), -0.999999)
    ts = r * math.sqrt((n - 2) / (1 - r**2))
    return r, float(2 * tdist.sf(abs(ts), n - 2))


def main():
    tax = load_gtdb()
    env = {r["sample"]: r for r in csv.DictReader(open(os.path.join(ECO, "env_matrix.tsv")), delimiter="\t")}

    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(AMRM, "stress_long.tsv"))]
    # cols: mag gene class subclass method set
    mag_metals = defaultdict(set)     # mag -> set(metal class)  (MAG set only)
    cls_count, phy_count, comp_count, set_count = Counter(), defaultdict(Counter), defaultdict(Counter), Counter()
    for mag, gene, mcls, sub, method, dset in rows:
        set_count[dset] += 1
        cls_count[mcls] += 1
        if dset == "MAG":
            mag_metals[mag].add(mcls)
            lin = tax.get(mag, "")
            ph = rank(lin, "p__") if lin else ""
            if ph:
                phy_count[mcls][ph] += 1
            samp = mag.split("__")[0]
            comp = env.get(samp, {}).get("compartment", "?")
            comp_count[mcls][comp] += 1

    # ---- census output ----
    with open(os.path.join(AMRM, "metal_census.tsv"), "w") as o:
        o.write("metal_class\tn_hits_all\tn_MAGs\tdom_phylum\tcompartment_PR_Surf_Deep\n")
        for mcls, n in cls_count.most_common():
            nmag = sum(1 for m, s in mag_metals.items() if mcls in s)
            domp = phy_count[mcls].most_common(1)[0][0] if phy_count[mcls] else ""
            c = comp_count[mcls]
            o.write(f"{mcls}\t{n}\t{nmag}\t{domp}\t"
                    f"{c.get('PlantRoot',0)}:{c.get('Surface',0)}:{c.get('Deep',0)}\n")

    # ---- per-sample resistance load for the big classes + cognate, then XRF corr ----
    # build per-sample MAG membership from stress MAGs (a sample's MAGs carrying a class)
    samp_load = defaultdict(lambda: defaultdict(set))   # sample -> class -> set(mag)
    for mag, metals in mag_metals.items():
        samp = mag.split("__")[0]
        for mcls in metals:
            samp_load[samp][mcls].add(mag)

    samples = [s for s in env if env[s].get("has_xrf") == "1"]
    corr_rows = []
    test_classes = ["ARSENIC", "MERCURY", "COPPER", "NICKEL", "CHROMATE"]
    for mcls in test_classes:
        load = np.array([len(samp_load.get(s, {}).get(mcls, ())) for s in samples], float)
        if load.sum() < 10:
            corr_rows.append((mcls, "n/a", "—", float("nan"), float("nan"), int(load.sum())))
            continue
        # cognate element (if measured) + exploratory scan
        scan = []
        for el in XRF_ELEMENTS:
            vals, ld = [], []
            for s, l in zip(samples, load):
                v = env[s].get(el, "")
                if v not in ("", "NA"):
                    try:
                        vals.append(float(v)); ld.append(l)
                    except ValueError:
                        pass
            if len(vals) >= 30:
                r, p = spearman(np.array(ld), np.array(vals))
                scan.append((el, r, p))
        scan = [s for s in scan if not math.isnan(s[1])]
        cog = COGNATE.get(mcls, "")
        cogr = next((f"{r:+.2f}(p={p:.2g})" for el, r, p in scan if el == cog), "not in XRF panel")
        top = sorted(scan, key=lambda t: t[2])[:1]
        tel, tr, tp = top[0] if top else ("", float("nan"), float("nan"))
        corr_rows.append((mcls, cog or "—", cogr, tel, tr, tp, int(load.sum())))

    with open(os.path.join(AMRM, "metal_xrf_corr.tsv"), "w") as o:
        o.write("metal_class\tcognate_element\tcognate_corr\ttop_xrf_correlate\ttop_r\ttop_p\tn_sample_load\n")
        for row in corr_rows:
            if len(row) == 6:
                mcls, cog, cogr, r, p, n = row
                o.write(f"{mcls}\t{cog}\t{cogr}\t—\t\t\t{n}\n")
            else:
                mcls, cog, cogr, tel, tr, tp, n = row
                o.write(f"{mcls}\t{cog}\t{cogr}\t{tel}\t{tr:+.3f}\t{tp:.2e}\t{n}\n")

    # ---- markdown ----
    n_mag_with_metal = len(mag_metals)
    ml = ["# Metal / metalloid resistance (AMRFinderPlus STRESS)\n",
          f"{set_count.get('MAG',0)} STRESS hits in {n_mag_with_metal} MAGs "
          f"(+{set_count.get('ISOLATE_ENRICH',0)} in isolates/enrichment). The metal-resistance "
          f"layer dropped when only Type=AMR was promoted.\n",
          "## Census (MAGs)\n",
          "| metal class | hits | MAGs | dominant phylum | PR:Surf:Deep |",
          "|---|---|---|---|---|"]
    for mcls, n in cls_count.most_common(8):
        nmag = sum(1 for m, s in mag_metals.items() if mcls in s)
        domp = phy_count[mcls].most_common(1)[0][0] if phy_count[mcls] else ""
        c = comp_count[mcls]
        ml.append(f"| {mcls} | {n} | {nmag} | {domp} | "
                  f"{c.get('PlantRoot',0)}:{c.get('Surface',0)}:{c.get('Deep',0)} |")
    ml += ["\n## XRF link (honest)\n",
           "XRF panel has **no As or Hg** — the two biggest classes can't be cognate-correlated. "
           "Cu/Ni have too few genes to power a test. Exploratory scan vs every measured element:\n",
           "| metal class | cognate | cognate corr | top XRF correlate | r | p | sample-load |",
           "|---|---|---|---|---|---|---|"]
    for row in corr_rows:
        if len(row) == 7:
            mcls, cog, cogr, tel, tr, tp, n = row
            ml.append(f"| {mcls} | {cog} | {cogr} | {tel} | {tr:+.2f} | {tp:.1e} | {n} |")
        else:
            mcls, cog, cogr, r, p, n = row
            ml.append(f"| {mcls} | {cog} | {cogr} | — | — | — | {n} |")
    open(os.path.join(AMRM, "METAL.md"), "w").write("\n".join(ml) + "\n")
    print(f"MAGs with metal resistance: {n_mag_with_metal}")
    print("census ->", os.path.join(AMRM, "metal_census.tsv"))
    for row in corr_rows:
        print(row)


if __name__ == "__main__":
    main()
