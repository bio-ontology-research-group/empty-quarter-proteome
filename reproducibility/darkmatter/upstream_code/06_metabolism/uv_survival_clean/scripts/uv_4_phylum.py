#!/usr/bin/env python3
"""Step 3b -- is the compartment effect real, or is it taxonomic composition in disguise?

The threat: photolyase and carotenoid genes are concentrated in particular phyla. If surface
communities simply contain more of those phyla, we would see "more photolyase at the surface"
with no compartment-specific biology at all -- the gene just rides along with the taxon.

Two tests:
  (1) does phylum composition itself differ Surface vs Deep?      (site-paired, abundance)
  (2) does the trait difference survive WITHIN a phylum?          (site-paired, occurrence
      recomputed using only species of that phylum)
If (2) holds across several independent phyla, the effect is not a composition artefact.

Out: results/phylum_composition_test.tsv , results/within_phylum_tests.tsv
"""
import csv, os, glob, collections, statistics
from scipy.stats import wilcoxon

M    = "/ibex/scratch/projects/c2014/rob/metagenomics"
IN   = f"{M}/06_metabolism/uv_survival_clean/results"        # per-species inputs (rep_traits.tsv)
RES  = os.environ.get("RES_OUT", IN)                         # v551 rerun: RES_OUT=.../results/v551 (2026-09-24)
os.makedirs(RES, exist_ok=True)
PROF = f"{M}/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"

FOCAL   = ["photoreactivation", "photoprotective_pigment", "dsb_nhej", "mismatch_repair"]
MIN_SPP = 15          # a phylum must have at least this many detected species in a sample to be scored
MIN_SITES = 15        # and appear in at least this many paired sites to be tested

rows = list(csv.DictReader(open(f"{IN}/rep_traits.tsv"), delimiter="\t"))
annotated = {r["rep"] for r in rows if r["annotated"] == "1"}
phylum    = {r["rep"]: r["phylum"] for r in rows}
trait     = {t: {r["rep"] for r in rows if r["annotated"] == "1" and r[t] == "1"} for t in FOCAL}

meta = {r["sample"]: r for r in csv.DictReader(open(f"{M}/metadata/sample_metadata.tsv"), delimiter="\t")}
have = {os.path.basename(d) for d in glob.glob(f"{PROF}/*") if os.path.isdir(d)}
samples = sorted(set(meta) & have)

TOP = [p for p, _ in collections.Counter(
    phylum[r] for r in annotated if phylum.get(r)).most_common(8)]
print(f"[phylum] testing within: {TOP}")

comp_abund = collections.defaultdict(dict)                  # (site,compartment) -> phylum -> [abund]
within     = collections.defaultdict(list)                  # (trait,phylum,site,compartment) -> [occ]
for s in samples:
    ab = {}
    with open(f"{PROF}/{s}/abund.tsv") as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 2 or c[0] == "unmapped":
                continue
            try:
                v = float(c[1])
            except ValueError:
                continue
            if v > 0 and c[0] in annotated:
                ab[c[0]] = v
    if not ab:
        print(f"[phylum] WARNING {s}: no mapped reads to an annotated rep, skipped")
        continue
    m = meta[s]; key = (m["site"], m["compartment"])
    tot = sum(ab.values()) or 1.0
    byp = collections.defaultdict(float)
    for g, v in ab.items():
        byp[phylum.get(g, "")] += v
    for p in TOP:
        comp_abund[key].setdefault(p, []).append(byp[p] / tot * 100)
    for p in TOP:
        det = {g for g in ab if phylum.get(g) == p}
        if len(det) < MIN_SPP:
            continue
        for t in FOCAL:
            within[(t, p, m["site"], m["compartment"])].append(len(det & trait[t]) / len(det) * 100)

# ---- (1) does composition itself shift? --------------------------------------
sites = sorted({k[0] for k in comp_abund})
with open(f"{RES}/phylum_composition_test.tsv", "w") as o:
    o.write("phylum\tn_sites\tsurface_pct\tdeep_pct\tdelta\tp_wilcoxon\n")
    print(f"\n=== (1) phylum composition, site-paired Surface vs Deep (% of annotated community) ===")
    print(f"{'phylum':22s} {'sites':>5s} {'Surf':>7s} {'Deep':>7s} {'delta':>7s} {'p':>10s}")
    for p in TOP:
        pr = [(statistics.median(comp_abund[(s, 'Surface')][p]),
               statistics.median(comp_abund[(s, 'Deep')][p]))
              for s in sites if (s, 'Surface') in comp_abund and (s, 'Deep') in comp_abund]
        if len(pr) < MIN_SITES:
            continue
        a = [x for x, _ in pr]; b = [y for _, y in pr]
        pv = wilcoxon(a, b).pvalue
        o.write(f"{p}\t{len(pr)}\t{statistics.median(a):.3f}\t{statistics.median(b):.3f}\t"
                f"{statistics.median(a)-statistics.median(b):+.3f}\t{pv:.4g}\n")
        print(f"{p:22s} {len(pr):5d} {statistics.median(a):7.2f} {statistics.median(b):7.2f} "
              f"{statistics.median(a)-statistics.median(b):+7.2f} {pv:10.2e}")

# ---- (2) does the trait effect survive within phylum? ------------------------
site_val = {k: statistics.median(v) for k, v in within.items()}
with open(f"{RES}/within_phylum_tests.tsv", "w") as o:
    o.write("trait\tphylum\tn_sites\tsurface_pct\tdeep_pct\tdelta\tsites_surface_higher\tp_wilcoxon\n")
    for t in FOCAL:
        print(f"\n=== (2) {t}: occurrence WITHIN phylum, site-paired Surface vs Deep ===")
        print(f"{'phylum':22s} {'sites':>5s} {'Surf':>7s} {'Deep':>7s} {'delta':>7s} {'S>D':>8s} {'p':>10s}")
        for p in TOP:
            pr = [(site_val[(t, p, s, 'Surface')], site_val[(t, p, s, 'Deep')])
                  for s in sites
                  if (t, p, s, 'Surface') in site_val and (t, p, s, 'Deep') in site_val]
            if len(pr) < MIN_SITES:
                continue
            a = [x for x, _ in pr]; b = [y for _, y in pr]
            d = [x - y for x, y in pr]; npos = sum(1 for x in d if x > 0)
            pv = wilcoxon(a, b).pvalue
            o.write(f"{t}\t{p}\t{len(pr)}\t{statistics.median(a):.3f}\t{statistics.median(b):.3f}\t"
                    f"{statistics.median(d):+.3f}\t{npos}\t{pv:.4g}\n")
            star = "  <<<" if pv < 0.05 else ""
            print(f"{p:22s} {len(pr):5d} {statistics.median(a):7.1f} {statistics.median(b):7.1f} "
                  f"{statistics.median(d):+7.1f} {npos:4d}/{len(pr):<3d} {pv:10.2e}{star}")

print(f"\n[phylum] wrote {RES}/phylum_composition_test.tsv and within_phylum_tests.tsv")
