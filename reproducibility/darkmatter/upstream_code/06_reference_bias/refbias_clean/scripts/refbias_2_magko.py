#!/usr/bin/env python3
"""refbias step 2 -- MAG-encoded KO profiles per sample, on the current catalogue.

Builds the metagenome side of the section 2.1 comparison: for every canonical soil sample that has
a CoverM profile, which KOs the MAG community actually encodes and at what abundance.

Join recipe is the one already used and shipped by sections 2.2 / 2.3 / 2.5
(lifeonair_3_abundance.py, uv_2_community.py, osmo_2_community.py) -- generalised from a handful
of marker KOs to the whole KO universe:

  rep_ko_v2.tsv                              rep -> set of KOs   (the CORRECTED table: 0 unannotated,
                                             _sub reps contig-subsetted, 8,134 spurious KOs removed)
  profiles_full_recompute/<sample>/abund.tsv rep -> relative abundance %, plus an `unmapped` row
  metadata/sample_metadata.tsv               the canonical 551 soil samples

  cohort = canonical 551 INTERSECT profiled  ->  548
  (51Sr1, F13PRr1, F5Dr1 have no profile: they are the CoverM array tasks that failed in the PI's
   573-sample run. We do not rerun them -- 548 is exactly the basis sections 2.2/2.3/2.5 ship on.)

Two normalisations per KO, as in 2.2/2.3:
  whole      sum of pct over reps carrying the KO, /100          -- includes unmapped, hard lower bound
  mapped     same sum / (100 - unmapped)                          -- of the MAG-mapped community only
Presence is any detection (pct > 0) of a carrying representative.

Rep-only by construction: CoverM maps against the 5,438 representatives (reps.mmi), so abundance is
representative-level and is NOT redistributed to cluster members -- same as the shipped sections.

  python3 scripts/refbias_2_magko.py
Writes:
  results/mag_ko_relab_by_sample.tsv.gz     KO x 548, mapped-normalised relative abundance
  results/mag_ko_presence_by_sample.tsv.gz  KO x 548, 0/1 presence
  results/mag_sample_stats.tsv              per-sample richness, unmapped fraction, reps detected
"""
import csv, os, glob, gzip

M = "/ibex/scratch/projects/c2014/rob/metagenomics"
W = f"{M}/06_reference_bias/refbias_clean"
REPKO = f"{M}/04_functional/repko_clean/results/rep_ko_v2.tsv"
PROF = f"{M}/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"

# ---- rep -> KO set -----------------------------------------------------------------
rep_ko = {}
for line in open(REPKO):
    c = line.rstrip("\n").split("\t")
    if not c[0]:
        continue
    rep_ko[c[0]] = [k for k in (c[1].split(",") if len(c) > 1 and c[1] else []) if k]
n_empty = sum(1 for v in rep_ko.values() if not v)
print(f"[magko] representatives in rep_ko_v2: {len(rep_ko)} | with 0 KOs: {n_empty}", flush=True)

# ---- cohort ------------------------------------------------------------------------
canon = [r["sample"] for r in csv.DictReader(open(f"{M}/metadata/sample_metadata.tsv"), delimiter="\t")]
profiled = {os.path.basename(d) for d in glob.glob(f"{PROF}/*") if os.path.isdir(d)}
cohort = sorted(set(canon) & profiled)
print(f"[magko] canonical={len(set(canon))} profiled={len(profiled)} usable={len(cohort)}", flush=True)
missing = sorted(set(canon) - profiled)
print(f"[magko] canonical WITHOUT a profile ({len(missing)}): {missing}", flush=True)

# ---- accumulate ---------------------------------------------------------------------
all_kos = sorted({k for v in rep_ko.values() for k in v})
ko_idx = {k: i for i, k in enumerate(all_kos)}
print(f"[magko] distinct KOs in the catalogue: {len(all_kos)}", flush=True)

relab = {k: {} for k in all_kos}   # ko -> sample -> mapped-normalised fraction
pres = {k: set() for k in all_kos}
stats = []

for s in cohort:
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
            elif v > 0:
                ab[c[0]] = v
    mapped = 100.0 - unmapped
    acc = {}
    for rep, pct in ab.items():
        for k in rep_ko.get(rep, ()):
            acc[k] = acc.get(k, 0.0) + pct
    for k, tot in acc.items():
        if k in relab:
            relab[k][s] = tot / mapped if mapped > 0 else 0.0
            pres[k].add(s)
    stats.append({"sample": s, "n_reps_detected": len(ab), "unmapped_pct": f"{unmapped:.4f}",
                  "mapped_pct": f"{mapped:.4f}", "n_ko_present": len(acc)})
    if len(stats) % 100 == 0:
        print(f"[magko] {len(stats)}/{len(cohort)} samples", flush=True)

# ---- write ---------------------------------------------------------------------------
os.makedirs(f"{W}/results", exist_ok=True)
with gzip.open(f"{W}/results/mag_ko_relab_by_sample.tsv.gz", "wt") as o:
    o.write("KO\t" + "\t".join(cohort) + "\n")
    for k in all_kos:
        row = relab[k]
        o.write(k + "\t" + "\t".join(f"{row[s]:.8g}" if s in row else "0" for s in cohort) + "\n")
with gzip.open(f"{W}/results/mag_ko_presence_by_sample.tsv.gz", "wt") as o:
    o.write("KO\t" + "\t".join(cohort) + "\n")
    for k in all_kos:
        p = pres[k]
        o.write(k + "\t" + "\t".join("1" if s in p else "0" for s in cohort) + "\n")
with open(f"{W}/results/mag_sample_stats.tsv", "w") as o:
    cols = ["sample", "n_reps_detected", "unmapped_pct", "mapped_pct", "n_ko_present"]
    o.write("\t".join(cols) + "\n")
    for d in stats:
        o.write("\t".join(str(d[c]) for c in cols) + "\n")

rich = sorted(d["n_ko_present"] for d in stats)
med = rich[len(rich) // 2]
print(f"[magko] median MAG-encoded KO richness per sample: {med}", flush=True)
print(f"[magko] -> results/mag_ko_{{relab,presence}}_by_sample.tsv.gz, mag_sample_stats.tsv", flush=True)
