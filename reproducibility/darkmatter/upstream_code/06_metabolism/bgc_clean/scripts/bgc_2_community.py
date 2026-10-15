#!/usr/bin/env python3
"""Step 2 -- per-SAMPLE community carriage of BGC traits.

Same abundance basis and dual-denominator design as sections 2.2/2.3
(lifeonair_3_abundance.py / uv_2_community.py): CoverM relative abundance of the
5,438 representatives, minimap2-sr, canonical soil samples with a profile (548).

  occurrence   = share of species DETECTED in the sample that carry the trait
  abund_whole  = % of all reads, incl. the ~51% mapping to no representative
                 (assumes unmapped organisms lack the trait -> hard LOWER BOUND)
  abund_mapped = % of the abundance that maps to any representative
                 (assumes unmapped organisms carry it at the same rate)

No unannotated exclusion: antiSMASH scanned every one of the 5,438 genomes, so a
zero-region representative is a genuine negative (see amr_2_community.py).

Traits = any_BGC + every product type carried by >=50 representatives, derived
from results/bgc_per_mag.tsv (bgc_1_census.py). terpene and terpene-precursor
stay separate (decided 2026-08-15).

Out: results/community_by_sample.tsv, results/community_summary.tsv

Run: python3 scripts/bgc_2_community.py   (after job 50533513 lands)
"""
import csv, os, glob, statistics, collections

M    = "/ibex/scratch/projects/c2014/rob/metagenomics"
B    = f"{M}/06_metabolism/bgc_clean"
RES  = f"{B}/results"
PROF = f"{M}/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"
DETECT = 0.0
FLOOR  = 50

# ---- traits from the census --------------------------------------------------
rows = list(csv.DictReader(open(f"{RES}/bgc_per_mag.tsv"), delimiter="\t"))
allreps = {r["mag"] for r in rows}
carriers = collections.defaultdict(set)
for r in rows:
    if int(r["n_BGCs"]) > 0:
        carriers["any_BGC"].add(r["mag"])
        for p in r["product_types"].split(";"):
            if p:
                carriers[p].add(r["mag"])
TRAITS = ["any_BGC"] + sorted((p for p in carriers if p != "any_BGC"
                               and len(carriers[p]) >= FLOOR),
                              key=lambda p: -len(carriers[p]))
trait = {t: carriers[t] for t in TRAITS}
print(f"[community] reps: {len(allreps)}; traits (>= {FLOOR} carriers): {TRAITS}")

meta = {r["sample"]: r for r in csv.DictReader(open(f"{M}/metadata/sample_metadata.tsv"), delimiter="\t")}
have = {os.path.basename(d) for d in glob.glob(f"{PROF}/*") if os.path.isdir(d)}
samples = sorted(set(meta) & have)
print(f"[community] canonical={len(meta)} profiled={len(have)} usable={len(samples)}")

def read_profile(s):
    ab, unmapped = {}, 0.0
    with open(f"{PROF}/{s}/abund.tsv") as fh:
        next(fh)
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 2:
                continue
            try:
                v = float(c[1])
            except ValueError:
                continue
            if c[0] == "unmapped":
                unmapped = v
            else:
                ab[c[0]] = v
    return ab, unmapped

out = open(f"{RES}/community_by_sample.tsv", "w")
out.write("sample\tsite\tcompartment\ttrait\tunmapped_pct\tn_detected\tn_detected_trait\t"
          "occurrence_pct\tabund_whole_pct\tabund_mapped_pct\n")
per = collections.defaultdict(lambda: collections.defaultdict(list))
overall = collections.defaultdict(list)
nsamp = 0
for s in samples:
    ab, un = read_profile(s)
    m = meta[s]
    det = {g for g, v in ab.items() if v > DETECT and g in allreps}
    mapped_abund = sum(v for g, v in ab.items() if g in allreps)
    if not det or mapped_abund <= 0:
        print(f"[community] WARNING {s}: no detected species, skipped")
        continue
    nsamp += 1
    for t in TRAITS:
        ts = trait[t]
        occ = len(det & ts) / len(det) * 100
        whole = sum(ab[g] for g in ts if g in ab)
        mapp = whole / mapped_abund * 100
        out.write(f"{s}\t{m['site']}\t{m['compartment']}\t{t}\t{un:.3f}\t{len(det)}\t"
                  f"{len(det & ts)}\t{occ:.4f}\t{whole:.4f}\t{mapp:.4f}\n")
        per[t][m["compartment"]].append((occ, whole, mapp))
        overall[t].append((occ, whole, mapp))
out.close()
print(f"[community] wrote {RES}/community_by_sample.tsv  ({nsamp} samples x {len(TRAITS)} traits)")

with open(f"{RES}/community_summary.tsv", "w") as o:
    o.write("trait\tcompartment\tn_samples\tocc_median\tocc_q1\tocc_q3\t"
            "abund_whole_median\tabund_whole_q1\tabund_whole_q3\tabund_mapped_median\n")
    for t in TRAITS:
        for c, v in [("ALL", overall[t])] + [(c, per[t][c]) for c in ("Surface", "Deep", "PlantRoot")]:
            if not v:
                continue
            occ = sorted(x[0] for x in v)
            wh  = sorted(x[1] for x in v)
            qo = statistics.quantiles(occ, n=4); qw = statistics.quantiles(wh, n=4)
            o.write(f"{t}\t{c}\t{len(v)}\t{statistics.median(occ):.2f}\t{qo[0]:.2f}\t{qo[2]:.2f}\t"
                    f"{statistics.median(wh):.2f}\t{qw[0]:.2f}\t{qw[2]:.2f}\t"
                    f"{statistics.median([x[2] for x in v]):.2f}\n")
print(f"[community] wrote {RES}/community_summary.tsv")

print(f"\n{'trait':22s} {'occ_med':>8s} {'whole_med':>10s} {'mapped_med':>10s}   (ALL samples)")
for t in TRAITS:
    v = overall[t]
    print(f"{t:22s} {statistics.median([x[0] for x in v]):8.1f} "
          f"{statistics.median([x[1] for x in v]):10.1f} "
          f"{statistics.median([x[2] for x in v]):10.1f}")
