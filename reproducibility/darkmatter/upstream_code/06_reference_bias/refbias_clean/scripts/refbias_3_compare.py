#!/usr/bin/env python3
"""refbias step 3 -- 16S/PICRUSt2 prediction versus MAG-encoded gene content.

The section 2.1 comparison. Deliberately framed as COVERAGE, not accuracy:

  The companion amplicon paper (empty-quarter-ecology) reports good agreement between PICRUSt2 and
  shotgun data -- per-sample rho 0.743, community-mean rho 0.806 -- computed over the 3,471 KOs the
  two methods SHARE. That result is not in conflict with this one and we do not contest it. It
  measures accuracy CONDITIONAL on a gene family being inside the KO reference. Section 2.1 measures
  the size of the outside: proteins with no KO at all cannot enter that intersection at any accuracy,
  and KOs outside PICRUSt2's own 10,543-family reference cannot be predicted even when the MAG
  encodes them. Both statements are true; they compose.

Inputs (all present on IBEX, none gated on the PI's stalled chain):
  PICRUSt2  rebuild_20260729/picrust2_current_20260803/inputs/ko_metagenome_unstrat.tsv
            10,543 KOs x 1,270 amplicon profiles. sha256 c71c2256... -- byte-identical to the file
            the companion ecology paper pins in its own input_manifest.tsv, so both papers use the
            same predictions. Staged + checksummed by the PI (INPUTS_COMPLETE).
            weighted_nsti.tsv gives per-profile placement confidence.
  MAG       refbias_2_magko.py output (rep_ko_v2 x CoverM profiles_full_recompute, 548 samples).

Cohort: canonical AND CoverM-profiled AND has an amplicon profile -> 448 samples.
Duplicate amplicon runs (4 of the 448 were sequenced twice) are collapsed by
`sum_then_within_sample_normalize`, the rule the companion paper documents in its own crosswalk.

  python3 scripts/refbias_3_compare.py
Writes:
  results/per_sample_comparison.tsv   per-sample recall / jaccard / rho / NSTI / richness
  results/comparison_summary.tsv      the headline numbers for the paper
  results/ko_universe.tsv             where each MAG-encoded KO stands w.r.t. the PICRUSt2 reference
"""
import csv, os, gzip, math

M = "/ibex/scratch/projects/c2014/rob/metagenomics"
W = f"{M}/06_reference_bias/refbias_clean"
IN = f"{M}/rebuild_20260729/picrust2_current_20260803/inputs"


def spearman(x, y):
    """Spearman rho with average ranks for ties; nan if degenerate."""
    n = len(x)
    if n < 3:
        return float("nan")
    def rank(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = rank(x), rank(y)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx > 0 and dy > 0 else float("nan")


def median(v):
    v = sorted(v)
    if not v:
        return float("nan")
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2.0


def quart(v, q):
    v = sorted(v)
    if not v:
        return float("nan")
    i = (len(v) - 1) * q
    lo, hi = int(math.floor(i)), int(math.ceil(i))
    return v[lo] if lo == hi else v[lo] * (hi - i) + v[hi] * (i - lo)


# ---- amplicon metadata: picrust2 column -> sample id -------------------------------
col2s, dup_of = {}, {}
for r in csv.DictReader(open(f"{IN}/sample_metadata.tsv"), delimiter="\t"):
    if r["is_control"] == "True":
        continue
    col2s[r["picrust2_col"]] = r["sample_id"]
    dup_of.setdefault(r["sample_id"], []).append(r["picrust2_col"])
nsti_col = {r["sample"]: float(r["weighted_NSTI"])
            for r in csv.DictReader(open(f"{IN}/weighted_nsti.tsv"), delimiter="\t")}

# ---- MAG side ----------------------------------------------------------------------
magf = f"{W}/results/mag_ko_relab_by_sample.tsv.gz"
with gzip.open(magf, "rt") as fh:
    mag_samples = fh.readline().rstrip("\n").split("\t")[1:]
    mag = {}
    for line in fh:
        c = line.rstrip("\n").split("\t")
        mag[c[0]] = [float(x) for x in c[1:]]
mag_si = {s: i for i, s in enumerate(mag_samples)}
print(f"[cmp] MAG side: {len(mag)} KOs x {len(mag_samples)} samples", flush=True)

cohort = sorted(s for s in mag_samples if s in dup_of)
print(f"[cmp] comparison cohort (MAG-profiled AND amplicon-profiled): {len(cohort)}", flush=True)
ndup = sum(1 for s in cohort if len(dup_of[s]) > 1)
print(f"[cmp] of these, sequenced twice on the amplicon side: {ndup}", flush=True)

# ---- PICRUSt2 side: read only the columns we need ----------------------------------
need_cols = {c for s in cohort for c in dup_of[s]}
with open(f"{IN}/ko_metagenome_unstrat.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    keep = [i for i, h in enumerate(hdr) if h in need_cols]
    kept_names = [hdr[i] for i in keep]
    pic = {}
    for line in fh:
        c = line.rstrip("\n").split("\t")
        pic[c[0]] = [float(c[i]) for i in keep]
pic_ref = set(pic)
print(f"[cmp] PICRUSt2 reference universe: {len(pic_ref)} KOs; columns read: {len(kept_names)}", flush=True)

ci = {c: i for i, c in enumerate(kept_names)}
# collapse duplicate runs: sum, then within-sample normalise (companion paper's rule)
pic_s = {}
for s in cohort:
    cols = [ci[c] for c in dup_of[s] if c in ci]
    tot = {}
    for k, v in pic.items():
        x = sum(v[i] for i in cols)
        if x > 0:
            tot[k] = x
    ssum = sum(tot.values())
    pic_s[s] = {k: v / ssum for k, v in tot.items()} if ssum > 0 else {}

# ---- per-sample comparison ---------------------------------------------------------
rows = []
for s in cohort:
    mi = mag_si[s]
    mset = {k for k, v in mag.items() if v[mi] > 0}
    pset = set(pic_s[s])
    shared = mset & pset
    outside_ref = mset - pic_ref            # MAG KO absent from PICRUSt2's whole reference
    inside_missed = (mset & pic_ref) - pset  # in reference, but not predicted here
    xs = [mag[k][mi] for k in sorted(shared)]
    ys = [pic_s[s][k] for k in sorted(shared)]
    nst = median([nsti_col[c] for c in dup_of[s] if c in nsti_col] or [float("nan")])
    rows.append({
        "sample": s, "n_mag_ko": len(mset), "n_pic_ko": len(pset), "n_shared": len(shared),
        "recall_of_mag": len(shared) / len(mset) if mset else float("nan"),
        "jaccard": len(shared) / len(mset | pset) if (mset | pset) else float("nan"),
        "mag_outside_pic_reference": len(outside_ref),
        "frac_mag_outside_reference": len(outside_ref) / len(mset) if mset else float("nan"),
        "mag_in_ref_not_predicted": len(inside_missed),
        "pic_only": len(pset - mset),
        "spearman_shared": spearman(xs, ys), "weighted_nsti": nst,
    })

cols = ["sample", "n_mag_ko", "n_pic_ko", "n_shared", "recall_of_mag", "jaccard",
        "mag_outside_pic_reference", "frac_mag_outside_reference", "mag_in_ref_not_predicted",
        "pic_only", "spearman_shared", "weighted_nsti"]
os.makedirs(f"{W}/results", exist_ok=True)
with open(f"{W}/results/per_sample_comparison.tsv", "w") as o:
    o.write("\t".join(cols) + "\n")
    for d in rows:
        o.write("\t".join(f"{d[c]:.6g}" if isinstance(d[c], float) else str(d[c]) for c in cols) + "\n")

# ---- KO universe ladder -------------------------------------------------------------
mag_any = {k for k, v in mag.items() if any(x > 0 for x in (v[mag_si[s]] for s in cohort))}
with open(f"{W}/results/ko_universe.tsv", "w") as o:
    o.write("quantity\tvalue\n")
    o.write(f"mag_kos_detected_in_cohort\t{len(mag_any)}\n")
    o.write(f"picrust_reference_kos\t{len(pic_ref)}\n")
    o.write(f"mag_kos_inside_picrust_reference\t{len(mag_any & pic_ref)}\n")
    o.write(f"mag_kos_outside_picrust_reference\t{len(mag_any - pic_ref)}\n")
    o.write(f"picrust_reference_kos_never_in_mags\t{len(pic_ref - mag_any)}\n")

# ---- summary -------------------------------------------------------------------------
def s_(name, v):
    return (f"{name}\t{median(v):.6g}\t{quart(v,0.25):.6g}\t{quart(v,0.75):.6g}\t"
            f"{min(v):.6g}\t{max(v):.6g}\n")

with open(f"{W}/results/comparison_summary.tsv", "w") as o:
    o.write(f"# cohort\t{len(cohort)}\tsamples with BOTH a CoverM MAG profile and a 16S/PICRUSt2 profile\n")
    o.write("metric\tmedian\tq25\tq75\tmin\tmax\n")
    for nm, key in [("weighted_nsti", "weighted_nsti"),
                    ("n_mag_ko", "n_mag_ko"), ("n_pic_ko", "n_pic_ko"), ("n_shared", "n_shared"),
                    ("recall_of_mag", "recall_of_mag"), ("jaccard", "jaccard"),
                    ("frac_mag_outside_reference", "frac_mag_outside_reference"),
                    ("spearman_shared", "spearman_shared")]:
        v = [d[key] for d in rows if not (isinstance(d[key], float) and math.isnan(d[key]))]
        o.write(s_(nm, v))

print("[cmp] --- headline ---", flush=True)
for key in ("weighted_nsti", "recall_of_mag", "spearman_shared", "frac_mag_outside_reference",
            "n_mag_ko", "n_pic_ko"):
    v = [d[key] for d in rows if not (isinstance(d[key], float) and math.isnan(d[key]))]
    print(f"[cmp]   {key:28s} median {median(v):.4g}  (IQR {quart(v,0.25):.4g}-{quart(v,0.75):.4g})", flush=True)
print(f"[cmp]   MAG KOs outside the PICRUSt2 reference entirely: {len(mag_any - pic_ref)} of {len(mag_any)}", flush=True)
print(f"[cmp] -> results/{{per_sample_comparison,comparison_summary,ko_universe}}.tsv", flush=True)
