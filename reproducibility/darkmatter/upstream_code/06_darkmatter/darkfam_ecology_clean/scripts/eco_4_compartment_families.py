#!/usr/bin/env python3
"""eco_4 — who the plant-root dark families are (§2.5 compartment extension).

Takes the persisted plant-root axis of eco_1 (results/indicators/plantroot_vs_other_all.tsv;
nothing is recomputed here) and joins the q<0.05 families to

  - phylum, via darkfam_basis_clean fam_phylum.tsv, counting a family only when its members
    agree on one phylum (n_phyla == 1; 77,069 of 78,728 attributable families), and
  - expedition persistence, via fam_presence.tsv x sample_site_trip.tsv: a family counts for a
    trip if it is detected in at least one sample of that trip, and the reported figure is the
    share detected in all five trips (same rule as the registered 79.4 % baseline).

Writes results/compartment_families_summary.tsv (counts, effect sizes, phylum shares,
over-representation ratios, persistence) and results/compartment_families_by_phylum.tsv
(the full per-phylum cross-tab). Rerunnable: both files are rewritten from scratch.
"""
import os, csv
import numpy as np

B = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter"
BASIS = os.path.join(B, "darkfam_basis_clean", "results")
PHYLUM = os.path.join(BASIS, "fam_phylum.tsv")
PRESENCE = os.path.join(BASIS, "fam_presence.tsv")
SITETRIP = os.path.join(BASIS, "sample_site_trip.tsv")
RES = os.path.join(B, "darkfam_ecology_clean", "results")
AXIS = os.path.join(RES, "indicators", "plantroot_vs_other_all.tsv")

QCUT = 0.05
TRIPS = ["trip1", "trip2", "trip3", "trip4", "trip5"]
SHARE_MIN_PCT = 2.0    # phyla listed in the summary share block
RATIO_MIN_FAM = 20     # over-representation ranking needs this many root-enriched families


def load_axis(path):
    """Persisted per-family table -> (family_id, r, q) for every tested family."""
    out = []
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out.append((row["family_id"], float(row["spearman_r"]), float(row["q"])))
    return out


def load_phylum(path):
    """family_rep -> phylum, single-phylum families only."""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["n_phyla"] == "1":
                out[row["family_rep"]] = row["phyla"]
    return out


def trips_detected(presence, sitetrip):
    """family_rep -> number of expeditions with at least one detection."""
    sample2trip = {}
    with open(sitetrip) as fh:
        for c in csv.reader(fh, delimiter="\t"):
            sample2trip[c[0]] = c[2]
    with open(presence) as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        cols_by_trip = [[j for j, s in enumerate(header) if sample2trip.get(s) == t]
                        for t in TRIPS]
        print(f"presence matrix: {len(header)} sample columns, "
              f"per trip {[len(c) for c in cols_by_trip]}, "
              f"{sum(1 for s in header if s not in sample2trip)} unmapped")
        out = {}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            vals = c[1:]
            out[c[0]] = sum(1 for idx in cols_by_trip
                            if any(vals[j] != "0" for j in idx))
    return out


def pct(num, den):
    return 100.0 * num / den if den else float("nan")


def quartiles(v):
    return (float(np.percentile(v, 25)), float(np.median(v)), float(np.percentile(v, 75))) \
        if len(v) else (float("nan"),) * 3


def main():
    axis = load_axis(AXIS)
    phylum = load_phylum(PHYLUM)
    ndet = trips_detected(PRESENCE, SITETRIP)

    sig = [(f, r) for f, r, q in axis if q < QCUT]
    enr = [(f, r) for f, r in sig if r > 0]
    dep = [(f, r) for f, r in sig if r < 0]
    print(f"tested={len(axis)} q<{QCUT}={len(sig)} enriched={len(enr)} depleted={len(dep)} "
          f"({pct(len(dep), len(sig)):.1f}% negative)")

    # ---- phylum composition: the two significant sets vs the attributable universe ----
    def counts(fams):
        c = {}
        for f, _ in fams:
            p = phylum.get(f)
            if p:
                c[p] = c.get(p, 0) + 1
        return c, sum(c.values())

    c_enr, n_enr_attr = counts(enr)
    c_dep, n_dep_attr = counts(dep)
    c_all = {}
    for p in phylum.values():
        c_all[p] = c_all.get(p, 0) + 1
    n_all_attr = sum(c_all.values())
    print(f"attributable (single-phylum): all={n_all_attr} enriched={n_enr_attr} "
          f"depleted={n_dep_attr}")

    ratio_enr, ratio_dep = {}, {}
    for p, n in c_all.items():
        base = pct(n, n_all_attr)
        ratio_enr[p] = pct(c_enr.get(p, 0), n_enr_attr) / base if base else float("nan")
        ratio_dep[p] = pct(c_dep.get(p, 0), n_dep_attr) / base if base else float("nan")

    # ---- expedition persistence ----
    base5 = pct(sum(1 for v in ndet.values() if v == 5), len(ndet))
    enr5 = pct(sum(1 for f, _ in enr if ndet.get(f) == 5), len(enr))
    dep5 = pct(sum(1 for f, _ in dep if ndet.get(f) == 5), len(dep))
    print(f"all-five-expedition persistence: baseline {base5:.1f}% "
          f"enriched {enr5:.1f}% depleted {dep5:.1f}%")

    # ---- effect sizes ----
    q1e, mede, q3e = quartiles(np.array([r for _, r in enr]))
    q1d, medd, q3d = quartiles(np.array([r for _, r in dep]))

    with open(os.path.join(RES, "compartment_families_summary.tsv"), "w") as out:
        w = lambda s, m, v: out.write(f"{s}\t{m}\t{v}\n")
        out.write("section\tmetric\tvalue\n")
        w("counts", "n_tested", len(axis))
        w("counts", "n_significant_q05", len(sig))
        w("counts", "n_root_enriched", len(enr))
        w("counts", "n_root_depleted", len(dep))
        w("counts", "pct_root_depleted_of_significant", f"{pct(len(dep), len(sig)):.1f}")
        w("counts", "n_root_enriched_phylum_attributable", n_enr_attr)
        w("counts", "n_root_depleted_phylum_attributable", n_dep_attr)
        w("counts", "n_all_phylum_attributable", n_all_attr)
        w("effect_size", "rho_root_enriched_q25", f"{q1e:.4f}")
        w("effect_size", "rho_root_enriched_median", f"{mede:.4f}")
        w("effect_size", "rho_root_enriched_q75", f"{q3e:.4f}")
        w("effect_size", "rho_root_depleted_q25", f"{q1d:.4f}")
        w("effect_size", "rho_root_depleted_median", f"{medd:.4f}")
        w("effect_size", "rho_root_depleted_q75", f"{q3d:.4f}")
        for label, cnt, tot in [("phylum_share_root_enriched_pct", c_enr, n_enr_attr),
                                ("phylum_share_root_depleted_pct", c_dep, n_dep_attr),
                                ("phylum_share_all_attributable_pct", c_all, n_all_attr)]:
            for p, n in sorted(cnt.items(), key=lambda kv: -kv[1]):
                if pct(n, tot) >= SHARE_MIN_PCT:
                    w(label, p, f"{pct(n, tot):.1f}")
        for p, rr in sorted(ratio_enr.items(), key=lambda kv: -kv[1]):
            if c_enr.get(p, 0) >= RATIO_MIN_FAM and rr >= 1.0:
                w("over_representation_root_enriched", p, f"{rr:.2f}")
        w("persistence", "pct_all_five_expeditions_baseline", f"{base5:.1f}")
        w("persistence", "pct_all_five_expeditions_root_enriched", f"{enr5:.1f}")
        w("persistence", "pct_all_five_expeditions_root_depleted", f"{dep5:.1f}")
        w("persistence", "n_families_in_presence_matrix", len(ndet))

    with open(os.path.join(RES, "compartment_families_by_phylum.tsv"), "w") as out:
        out.write("phylum\tn_all_attributable\tpct_all_attributable\tn_root_enriched\t"
                  "pct_root_enriched\tn_root_depleted\tpct_root_depleted\t"
                  "ratio_enriched_vs_all\tratio_depleted_vs_all\n")
        for p in sorted(c_all, key=lambda k: (-c_enr.get(k, 0), -c_all[k])):
            out.write(f"{p}\t{c_all[p]}\t{pct(c_all[p], n_all_attr):.2f}\t"
                      f"{c_enr.get(p, 0)}\t{pct(c_enr.get(p, 0), n_enr_attr):.2f}\t"
                      f"{c_dep.get(p, 0)}\t{pct(c_dep.get(p, 0), n_dep_attr):.2f}\t"
                      f"{ratio_enr[p]:.2f}\t{ratio_dep[p]:.2f}\n")

    print("done ->", RES)


if __name__ == "__main__":
    main()
