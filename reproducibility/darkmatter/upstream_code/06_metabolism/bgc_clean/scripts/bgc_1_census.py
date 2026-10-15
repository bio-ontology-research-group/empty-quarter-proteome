#!/usr/bin/env python3
"""Step 1 -- antiSMASH BGC census on the CURRENT 5,438-representative catalogue.

Aggregates results/antismash/<mag>/<mag>.json (job 50413466, antiSMASH 8.0.4 --minimal,
5,438/5,438 complete) into paper-ready census tables. Ports the counting semantics of
scripts/aggregate_mobilome.py (the historical 15,469-MAG aggregator) exactly:

  - product counts walk records[].areas[].products; a multi-product area contributes
    once per product; n_MAGs is deduplicated per MAG.
  - n_BGCs per MAG = number of areas (cross-checked against *.region*.gbk file count).

Two layers the historical aggregator did not have:
  - broad-category rollup via areas[].protoclusters[*].category (antiSMASH 8 puts the
    category only there; area['categories'] is None in 8.0.4 output). Reported
    alongside, NOT replacing, the product table: rolling terpene-precursor into
    terpene would break comparability with the historical census (decided 2026-08-15,
    keep them separate).
  - contig_edge, read from features[type=='region'].qualifiers, so the fragmentation
    caveat is measured on our own data instead of quoting the historical 90%.

Taxonomy = 05_refine/_gtdbtk_r232 (R232, covers exactly the 5,438 reps); the old
02_taxonomy/gtdbtk_out join is the retired basis. Column is named tax_class to avoid
the historical two-columns-named-class collision.

Out (results/):
  bgc_per_mag.tsv          per MAG: metadata, n_BGCs, n_edge_BGCs, product_types (;-set)
  bgc_product_summary.tsv  per product: n_BGCs, n_edge_BGCs, n_MAGs, pct_species
  bgc_category_summary.tsv per protocluster category: n_protoclusters, n_MAGs, pct_species
  bgc_census_summary.tsv   headline scalars (totals, edge fraction, any-BGC prevalence)

Run: sbatch -p batch -c 16 --mem 32G -t 1:00:00 --wrap "python3 scripts/bgc_1_census.py"
"""
import csv, glob, json, os, re, sys
from multiprocessing import Pool

M   = "/ibex/scratch/projects/c2014/rob/metagenomics"
B   = f"{M}/06_metabolism/bgc_clean"
RES = f"{B}/results"
GTDB = [f"{M}/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv",
        f"{M}/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"]

# ---- metadata helpers (ported verbatim from scripts/aggregate_mobilome.py) ----
TRIP = {'T': 2, 'F': 3, 'S': 4, 'V': 5}
comp_re = re.compile(r'(PR|S|D)r\d')
COMP = {'PR': 'PlantRoot', 'S': 'Surface', 'D': 'Deep'}
def sample_of(mag): return mag.split('__', 1)[0]
def trip_of(s): return 1 if s[:1].isdigit() else TRIP.get(s[0], 0)
def comp_of(s):
    m = comp_re.search(s); return COMP.get(m.group(1), 'NA') if m else 'NA'
def site_of(s):
    m = re.match(r'^[TFSV]?(\d+)', s); return m.group(1) if m else 'NA'

RANKS = ['phylum', 'class', 'order', 'family', 'genus']
tax = {}
for p in GTDB:
    for r in csv.DictReader(open(p), delimiter="\t"):
        d = {}
        for tok in r["classification"].split(';'):
            tok = tok.strip()
            for rk, pre in zip(RANKS, ['p__', 'c__', 'o__', 'f__', 'g__']):
                if tok.startswith(pre):
                    d[rk] = tok[len(pre):] or 'unclassified'
        tax[r["user_genome"]] = d
def meta(mag):
    s = sample_of(mag)
    d = tax.get(mag, {})
    return [s, str(trip_of(s)), site_of(s), comp_of(s)] + [d.get(r, 'NA') for r in RANKS]
META_H = ['sample', 'trip', 'site', 'compartment',
          'phylum', 'tax_class', 'order', 'family', 'genus']

mags = [l.split('\t', 1)[0] for l in open(f"{RES}/rep_index.tsv")]
print(f"[census] MAGs: {len(mags)}; refined-GTDB entries: {len(tax)}", flush=True)
assert len(mags) == 5438 and len(tax) == 5438

# ---- per-MAG JSON walk (parallel) --------------------------------------------
def scan(mag):
    d = f"{RES}/antismash/{mag}"
    J = json.load(open(f"{d}/{mag}.json"))
    prods, cats = [], []
    n_areas = 0
    n_edge = 0
    edge_prods = []
    for rec in J.get('records', []):
        for area in rec.get('areas', []):
            n_areas += 1
            prods += area.get('products', [])
            for pc in (area.get('protoclusters') or {}).values():
                cats.append(pc.get('category', 'NA'))
        for f in rec.get('features', []):
            if f.get('type') != 'region':
                continue
            q = f.get('qualifiers', {})
            if q.get('contig_edge', ['False'])[0] == 'True':
                n_edge += 1
                edge_prods += q.get('product', [])
    n_gbk = len(glob.glob(f"{d}/*.region*.gbk"))
    return mag, n_areas, n_gbk, n_edge, prods, cats, edge_prods

with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', 8))) as pool:
    scanned = pool.map(scan, mags, chunksize=25)

# ---- aggregate ---------------------------------------------------------------
from collections import Counter, defaultdict
prod_cnt, prod_edge, cat_cnt = Counter(), Counter(), Counter()
prod_mags, cat_mags = defaultdict(set), defaultdict(set)
per_mag, mismatch = [], 0
tot_areas = tot_edge = n_any = 0
for mag, n_areas, n_gbk, n_edge, prods, cats, edge_prods in scanned:
    if n_areas != n_gbk:
        mismatch += 1
    tot_areas += n_areas; tot_edge += n_edge
    if n_areas > 0:
        n_any += 1
        for p in set(prods): prod_mags[p].add(mag)
        for p in prods:      prod_cnt[p] += 1
        for p in edge_prods: prod_edge[p] += 1
        for c in set(cats):  cat_mags[c].add(mag)
        for c in cats:       cat_cnt[c] += 1
    per_mag.append([mag] + meta(mag) +
                   [str(n_areas), str(n_edge), ';'.join(sorted(set(prods)))])
print(f"[census] areas-vs-gbk mismatches: {mismatch}", flush=True)

N = len(mags)
def w(path, header, rows):
    with open(path, 'w', newline='') as o:
        wr = csv.writer(o, delimiter='\t'); wr.writerow(header)
        for r in rows: wr.writerow(r)
    print(f"  wrote {path} ({len(rows)} rows)", flush=True)

w(f"{RES}/bgc_per_mag.tsv",
  ['mag'] + META_H + ['n_BGCs', 'n_edge_BGCs', 'product_types'], per_mag)
w(f"{RES}/bgc_product_summary.tsv",
  ['product_type', 'n_BGCs', 'n_edge_BGCs', 'n_MAGs', 'pct_species'],
  [[p, prod_cnt[p], prod_edge.get(p, 0), len(prod_mags[p]),
    f"{100*len(prod_mags[p])/N:.2f}"]
   for p in sorted(prod_cnt, key=lambda x: -len(prod_mags[x]))])
w(f"{RES}/bgc_category_summary.tsv",
  ['category', 'n_protoclusters', 'n_MAGs', 'pct_species'],
  [[c, cat_cnt[c], len(cat_mags[c]), f"{100*len(cat_mags[c])/N:.2f}"]
   for c in sorted(cat_cnt, key=lambda x: -len(cat_mags[x]))])
w(f"{RES}/bgc_census_summary.tsv", ['metric', 'value'],
  [['n_species', N],
   ['n_species_with_BGC', n_any],
   ['pct_species_with_BGC', f"{100*n_any/N:.2f}"],
   ['total_regions', tot_areas],
   ['regions_on_contig_edge', tot_edge],
   ['pct_regions_contig_edge', f"{100*tot_edge/max(tot_areas,1):.2f}"],
   ['areas_vs_gbk_mismatches', mismatch]])
print("[census] done", flush=True)
