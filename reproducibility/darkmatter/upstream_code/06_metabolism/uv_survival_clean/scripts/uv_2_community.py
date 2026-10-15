#!/usr/bin/env python3
"""Step 2 -- turn per-species repair traits into per-SAMPLE community measurements.

Same abundance basis as section 2.2 (lifeonair_3_abundance.py): CoverM relative abundance of the
5,438 representatives, minimap2-sr, over the canonical soil samples that have a profile (548).

For every sample and every trait we report two independent views:

  occurrence  = share of the species DETECTED in that sample which carry the trait
                (unweighted; "of the species that live here, how many can do this?")
  abundance   = share of the community BY ABUNDANCE carrying the trait, under two denominators:
                  whole      -- % of all reads incl. the ~51% that map to no representative
                                (assumes unmapped organisms lack the trait -> hard LOWER BOUND)
                  annotated  -- % of the abundance that maps to an ANNOTATED representative
                                (assumes unmapped organisms carry it at the same rate -> extrapolation)

The 124 representatives with no annotation (see uv_0_patch_repko.py) are UNKNOWN, not negative:
they are removed from both numerator and denominator of every metric.

Out: results/community_by_sample.tsv  (one row per sample per trait)
     results/community_summary.tsv    (medians across samples, by compartment)
"""
import csv, os, glob, statistics, collections

M    = "/ibex/scratch/projects/c2014/rob/metagenomics"
U    = f"{M}/06_metabolism/uv_survival_clean"
IN   = f"{U}/results"                              # per-species inputs (rep_traits.tsv) -- cohort-independent
RES  = os.environ.get("RES_OUT", IN)               # v551 rerun: RES_OUT=.../results/v551 (2026-09-24)
os.makedirs(RES, exist_ok=True)
PROF = f"{M}/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"

DETECT = 0.0        # primary detection rule: relative abundance strictly greater than this
DETECT2 = 0.01      # sensitivity detection rule (%), reported alongside

# ---- traits ------------------------------------------------------------------
rows = list(csv.DictReader(open(f"{IN}/rep_traits.tsv"), delimiter="\t"))
TRAITS = [c for c in rows[0] if c not in ("rep", "annotated", "phylum", "n_panel_kos")]
trait = {t: {r["rep"] for r in rows if r["annotated"] == "1" and r[t] == "1"} for t in TRAITS}
annotated = {r["rep"] for r in rows if r["annotated"] == "1"}
phylum = {r["rep"]: r["phylum"] for r in rows}
print(f"[community] traits: {TRAITS}")
print(f"[community] annotated representatives: {len(annotated)} of {len(rows)}")

# ---- samples -----------------------------------------------------------------
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
out.write("sample\tsite\tcompartment\ttrait\tunmapped_pct\tn_detected_ann\tn_detected_trait\t"
          "occurrence_pct\toccurrence_pct_d2\tabund_whole_pct\tabund_annotated_pct\n")
per = collections.defaultdict(lambda: collections.defaultdict(list))
nsamp = 0
for s in samples:
    ab, un = read_profile(s)
    m = meta[s]
    det  = {g for g, v in ab.items() if v > DETECT and g in annotated}
    det2 = {g for g, v in ab.items() if v > DETECT2 and g in annotated}
    ann_abund = sum(v for g, v in ab.items() if g in annotated)
    if not det or ann_abund <= 0:
        print(f"[community] WARNING {s}: no detected annotated species, skipped")
        continue
    nsamp += 1
    for t in TRAITS:
        ts = trait[t]
        occ  = len(det & ts) / len(det) * 100
        occ2 = (len(det2 & ts) / len(det2) * 100) if det2 else float("nan")
        whole = sum(ab[g] for g in ts if g in ab)
        annp  = whole / ann_abund * 100
        out.write(f"{s}\t{m['site']}\t{m['compartment']}\t{t}\t{un:.3f}\t{len(det)}\t{len(det & ts)}\t"
                  f"{occ:.4f}\t{occ2:.4f}\t{whole:.4f}\t{annp:.4f}\n")
        per[t][m["compartment"]].append((occ, whole, annp))
out.close()
print(f"[community] wrote {RES}/community_by_sample.tsv  ({nsamp} samples x {len(TRAITS)} traits)")

with open(f"{RES}/community_summary.tsv", "w") as o:
    o.write("trait\tcompartment\tn_samples\tocc_median\tocc_q1\tocc_q3\t"
            "abund_whole_median\tabund_annotated_median\n")
    for t in TRAITS:
        for c in ("Surface", "Deep", "PlantRoot"):
            v = per[t][c]
            if not v:
                continue
            occ = sorted(x[0] for x in v); q = statistics.quantiles(occ, n=4)
            o.write(f"{t}\t{c}\t{len(v)}\t{statistics.median(occ):.2f}\t{q[0]:.2f}\t{q[2]:.2f}\t"
                    f"{statistics.median([x[1] for x in v]):.2f}\t"
                    f"{statistics.median([x[2] for x in v]):.2f}\n")

print(f"\n{'trait':28s} {'Surface':>18s} {'Deep':>18s} {'PlantRoot':>18s}    (median occurrence %)")
for t in TRAITS:
    cells = []
    for c in ("Surface", "Deep", "PlantRoot"):
        v = per[t][c]
        cells.append(f"{statistics.median([x[0] for x in v]):.1f}" if v else "-")
    print(f"{t:28s} {cells[0]:>18s} {cells[1]:>18s} {cells[2]:>18s}")
print(f"\n[community] wrote {RES}/community_summary.tsv")
