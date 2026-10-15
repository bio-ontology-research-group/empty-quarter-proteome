#!/usr/bin/env python3
"""Step 3 -- does the repair repertoire track sun exposure? Site-paired compartment test.

Design. Replicate samples within a site x compartment cell are not independent, and the cells are
badly unbalanced (1 to 10 samples; 159 cells over 61 sites). So we COLLAPSE FIRST: take the median
across replicates within each site x compartment, giving one value per site per compartment, then
compare Surface vs Deep as a PAIRED test across sites that have both. Pairing holds site, climate
and geochemistry constant, leaving sun exposure as what differs.

  primary   Wilcoxon signed-rank on the per-site Surface - Deep difference
  secondary Kruskal-Wallis across all three compartments on site-level values
  control   the near-universal groups (NER, BER, SOS, oxidative) must stay FLAT. A pipeline that
            manufactures compartment signal would move them too.

Multiple testing: Holm-Bonferroni across the traits tested, per metric.

Out: results/compartment_tests.tsv , results/site_level_values.tsv
"""
import csv, os, collections, statistics
from scipy.stats import wilcoxon, kruskal

M   = "/ibex/scratch/projects/c2014/rob/metagenomics"
RES = os.environ.get("RES_OUT", f"{M}/06_metabolism/osmolyte_clean/results")  # v551 rerun (2026-09-24)
os.makedirs(RES, exist_ok=True)

# the focal-trait aliases duplicate their functional groups; keep the group names only
DROP = {"trehalose_strict", "glucosylglycerol_strict", "mannosylglycerate_strict"}
METRICS = ["occurrence_pct", "abund_whole_pct", "abund_annotated_pct"]

rows = list(csv.DictReader(open(f"{RES}/community_by_sample.tsv"), delimiter="\t"))
traits = sorted({r["trait"] for r in rows} - DROP)

# ---- collapse to one value per site x compartment ----------------------------
cell = collections.defaultdict(list)          # (trait, metric, site, compartment) -> [values]
for r in rows:
    for m in METRICS:
        cell[(r["trait"], m, r["site"], r["compartment"])].append(float(r[m]))
site_val = {k: statistics.median(v) for k, v in cell.items()}

with open(f"{RES}/site_level_values.tsv", "w") as o:
    o.write("trait\tmetric\tsite\tcompartment\tn_samples\tmedian_value\n")
    for k, v in sorted(site_val.items()):
        o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]}\t{len(cell[k])}\t{v:.4f}\n")

sites = sorted({r["site"] for r in rows})
print(f"[test] {len(rows)} sample-rows -> {len(sites)} sites; collapsing replicates by median")

def holm(pvals):
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    adj, running = [0.0] * len(pvals), 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(pvals) - rank) * pvals[i])
        adj[i] = min(1.0, running)
    return adj

results = []
for metric in METRICS:
    stash = []
    for t in traits:
        S = [(s, site_val.get((t, metric, s, "Surface")), site_val.get((t, metric, s, "Deep")))
             for s in sites]
        pairs = [(a, b) for _, a, b in S if a is not None and b is not None]
        if len(pairs) < 10:
            continue
        surf = [a for a, _ in pairs]; deep = [b for _, b in pairs]
        diff = [a - b for a, b in pairs]
        try:
            w = wilcoxon(surf, deep)
            p = w.pvalue
        except ValueError:
            p = float("nan")
        # Kruskal across the three compartments at site level
        groups = []
        for c in ("Surface", "Deep", "PlantRoot"):
            g = [site_val[(t, metric, s, c)] for s in sites if (t, metric, s, c) in site_val]
            if g:
                groups.append(g)
        kp = kruskal(*groups).pvalue if len(groups) == 3 else float("nan")
        n_pos = sum(1 for d in diff if d > 0)
        stash.append(dict(trait=t, metric=metric, n_sites=len(pairs),
                          surface=statistics.median(surf), deep=statistics.median(deep),
                          delta=statistics.median(diff), sites_surface_higher=n_pos,
                          p_wilcoxon=p, p_kruskal3=kp))
    adj = holm([s["p_wilcoxon"] for s in stash])
    for s, a in zip(stash, adj):
        s["p_holm"] = a
    results += stash

with open(f"{RES}/compartment_tests.tsv", "w") as o:
    cols = ["trait", "metric", "n_sites", "surface", "deep", "delta", "sites_surface_higher",
            "p_wilcoxon", "p_holm", "p_kruskal3"]
    o.write("\t".join(cols) + "\n")
    for s in results:
        o.write("\t".join(f"{s[c]:.4g}" if isinstance(s[c], float) else str(s[c]) for c in cols) + "\n")

for metric in METRICS:
    print(f"\n=== {metric} -- site-paired Surface vs Deep ===")
    print(f"{'trait':28s} {'sites':>5s} {'Surf':>7s} {'Deep':>7s} {'delta':>7s} {'S>D':>7s} {'p(Holm)':>10s}")
    for s in [r for r in results if r["metric"] == metric]:
        flag = "  <<<" if s["p_holm"] < 0.05 else ""
        print(f"{s['trait']:28s} {s['n_sites']:5d} {s['surface']:7.1f} {s['deep']:7.1f} "
              f"{s['delta']:+7.1f} {s['sites_surface_higher']:4d}/{s['n_sites']:<3d} "
              f"{s['p_holm']:10.2e}{flag}")
print(f"\n[test] wrote {RES}/compartment_tests.tsv and site_level_values.tsv")
