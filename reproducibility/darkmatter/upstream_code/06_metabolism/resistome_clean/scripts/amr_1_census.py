#!/usr/bin/env python3
"""Step 1 -- AMRFinderPlus census on the CURRENT 5,438-representative catalogue.

Aggregates results/amrfinder/<mag>.amrfinder.tsv (job 50413467, AMRFinderPlus 4.2.7,
DB 2026-05-15.1, --plus; 5,438/5,438 complete, 2,747 with >=1 hit) into census tables.

Ports scripts/aggregate_mobilome.py lines 101-119 (the historical 15,469-MAG
aggregator): the antibiotic-resistome table keeps its exact convention -- rows with
Type == 'AMR' only, grouped on Class, n_MAGs deduplicated per MAG -- so the new
census is like-for-like against the historical percentages (BETA-LACTAM 3,816/15,469
= 24.67% etc.). Scope is recorded but NOT filtered on, matching the old convention
(core + plus pooled; the core/plus split is reported in the summary).

What the historical aggregator dropped -- and 06_resistome_safety/amr_metal/METAL.md
flags as a defect ("The metal-resistance layer dropped when only Type=AMR was
promoted") -- is recovered here as its own layer: Type == 'STRESS' rows (metal /
metalloid / biocide resistance), grouped on Class and Subclass. VIRULENCE-type rows
are counted in the summary only; the virulence analysis is the PI's VFDB census
(06_resistome_safety/vfdb/, current basis), reused by reference.

Taxonomy = 05_refine/_gtdbtk_r232 (R232, exactly the 5,438 reps). The GTDB class
column is named tax_class and the AMR class column amr_class: the historical
amr_long.tsv had two columns literally named 'class', which silently collapse in any
DictReader/pandas load.

Out (results/):
  amr_long.tsv               one row per Type==AMR hit: metadata, gene, amr_class,
                             subclass, scope, method
  amr_class_summary.tsv      per amr_class: n_hits, n_MAGs, pct_species
  stress_long.tsv            one row per Type==STRESS hit (same columns + subtype)
  stress_class_summary.tsv   per STRESS Class: n_hits, n_MAGs, pct_species
  amr_census_summary.tsv     headline scalars (type/scope splits, any-AMR prevalence)
  rep_amr_traits.tsv         per representative: flags for the community layer
                             (any_amr, per major amr_class, ARSENIC, MERCURY)

Run: python3 scripts/amr_1_census.py    (TSV parsing only; seconds)
"""
import csv, os, re
from collections import Counter, defaultdict

M   = "/ibex/scratch/projects/c2014/rob/metagenomics"
B   = f"{M}/06_metabolism/resistome_clean"
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

mags = [l.split('\t', 1)[0] for l in open(f"{RES}/rep_faa_index.tsv")]
print(f"[census] MAGs: {len(mags)}; refined-GTDB entries: {len(tax)}")
assert len(mags) == 5438 and len(tax) == 5438
N = len(mags)

# ---- walk the per-MAG TSVs ---------------------------------------------------
amr_long, stress_long = [], []
amr_cls_hits, str_cls_hits = Counter(), Counter()
amr_cls_mags, str_cls_mags = defaultdict(set), defaultdict(set)
type_cnt, scope_amr = Counter(), Counter()
amr_mags, str_mags = set(), set()
mag_amr_classes = defaultdict(set)   # for rep_amr_traits
n_files = 0
for mag in mags:
    p = f"{RES}/amrfinder/{mag}.amrfinder.tsv"
    n_files += 1
    with open(p) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            typ = row['Type']
            type_cnt[typ] += 1
            gene = row['Element symbol']
            cls, sub = row['Class'], row['Subclass']
            if typ == 'AMR':
                scope_amr[row['Scope']] += 1
                amr_long.append([mag] + meta(mag) +
                                [gene, cls, sub, row['Scope'], row['Method']])
                amr_cls_hits[cls] += 1; amr_cls_mags[cls].add(mag)
                amr_mags.add(mag); mag_amr_classes[mag].add(cls)
            elif typ == 'STRESS':
                stress_long.append([mag] + meta(mag) +
                                   [gene, cls, sub, row['Subtype'],
                                    row['Scope'], row['Method']])
                str_cls_hits[cls] += 1; str_cls_mags[cls].add(mag)
                str_mags.add(mag); mag_amr_classes[mag].add(f"STRESS:{cls}")
print(f"[census] files: {n_files}; type counts: {dict(type_cnt)}")
print(f"[census] AMR scope split: {dict(scope_amr)}")

def w(path, header, rows):
    with open(path, 'w', newline='') as o:
        wr = csv.writer(o, delimiter='\t'); wr.writerow(header)
        for r in rows: wr.writerow(r)
    print(f"  wrote {path} ({len(rows)} rows)")

w(f"{RES}/amr_long.tsv",
  ['mag'] + META_H + ['gene', 'amr_class', 'subclass', 'scope', 'method'], amr_long)
w(f"{RES}/amr_class_summary.tsv",
  ['amr_class', 'n_hits', 'n_MAGs', 'pct_species'],
  [[c, amr_cls_hits[c], len(amr_cls_mags[c]), f"{100*len(amr_cls_mags[c])/N:.2f}"]
   for c in sorted(amr_cls_hits, key=lambda x: -len(amr_cls_mags[x]))])
w(f"{RES}/stress_long.tsv",
  ['mag'] + META_H + ['gene', 'stress_class', 'subclass', 'subtype',
                      'scope', 'method'], stress_long)
w(f"{RES}/stress_class_summary.tsv",
  ['stress_class', 'n_hits', 'n_MAGs', 'pct_species'],
  [[c, str_cls_hits[c], len(str_cls_mags[c]), f"{100*len(str_cls_mags[c])/N:.2f}"]
   for c in sorted(str_cls_hits, key=lambda x: -len(str_cls_mags[x]))])
w(f"{RES}/amr_census_summary.tsv", ['metric', 'value'],
  [['n_species', N],
   ['n_species_any_hit', len({m for m in mags
                              if m in amr_mags or m in str_mags})],
   ['n_species_AMR', len(amr_mags)],
   ['pct_species_AMR', f"{100*len(amr_mags)/N:.2f}"],
   ['n_species_STRESS', len(str_mags)],
   ['pct_species_STRESS', f"{100*len(str_mags)/N:.2f}"],
   ['n_hits_AMR', sum(amr_cls_hits.values())],
   ['n_hits_STRESS', sum(str_cls_hits.values())],
   ['n_hits_VIRULENCE_type', type_cnt.get('VIRULENCE', 0)],
   ['AMR_scope_core', scope_amr.get('core', 0)],
   ['AMR_scope_plus', scope_amr.get('plus', 0)]])

# ---- per-rep trait flags for the community layer (step 2) --------------------
FLOOR = 50   # an amr_class/stress class becomes a community trait at >=50 carrier MAGs
traits = (['any_AMR'] +
          [c for c in sorted(amr_cls_hits, key=lambda x: -len(amr_cls_mags[x]))
           if len(amr_cls_mags[c]) >= FLOOR] +
          [f"STRESS:{c}" for c in sorted(str_cls_hits,
                                         key=lambda x: -len(str_cls_mags[x]))
           if len(str_cls_mags[c]) >= FLOOR])
with open(f"{RES}/rep_amr_traits.tsv", 'w', newline='') as o:
    wr = csv.writer(o, delimiter='\t')
    wr.writerow(['rep', 'phylum'] + traits)
    for mag in mags:
        has = mag_amr_classes.get(mag, set())
        row = [mag, tax.get(mag, {}).get('phylum', 'NA')]
        row.append('1' if mag in amr_mags else '0')
        for t in traits[1:]:
            row.append('1' if t in has else '0')
        wr.writerow(row)
print(f"  wrote {RES}/rep_amr_traits.tsv (traits: {traits})")
print("[census] done")
