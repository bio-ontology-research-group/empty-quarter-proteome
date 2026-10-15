#!/usr/bin/env python3
"""Read-only computational audit for Round 4 review of empty-quarter-mags."""
import csv
import gzip
import json
import math
import os
import statistics as st
import sys
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "/home/leechuck/Documents/papers/empty-quarter-mags"
DM = os.path.join(ROOT, "reproducibility/darkmatter")
INP = os.path.join(DM, "inputs")
SNAP = os.path.join(DM, "source_snapshot")
OUT_DIR = os.path.join(ROOT, "reviews/darkmatter-2026-10-13/round4/reviewer_checks")

# Index all files: SNAP first, then DM, then INP last so INP always wins!
FILE_MAP = {}
for base in [SNAP, DM, INP]:
    for dirpath, _, filenames in os.walk(base):
        if dirpath.startswith(SNAP) and base != SNAP:
            continue
        if dirpath.startswith(INP) and base != INP:
            continue
        for fn in filenames:
            FILE_MAP[fn] = os.path.join(dirpath, fn)


def fpath(name):
    if name in FILE_MAP:
        return FILE_MAP[name]
    if name + ".gz" in FILE_MAP:
        return FILE_MAP[name + ".gz"]
    if name.endswith(".gz") and name[:-3] in FILE_MAP:
        return FILE_MAP[name[:-3]]
    raise FileNotFoundError(name)


out_log = open(os.path.join(OUT_DIR, "r4_audit_part1.txt"), "w")


def pr(*args, **kwargs):
    print(*args, **kwargs)
    print(*args, **kwargs, file=out_log)


def section(title):
    pr("\n" + "=" * 78)
    pr(title)
    pr("=" * 78)


# 1. Partial Spearman vs Marginal Spearman in partial_summary.tsv.gz
section("1. ECOLOGY PARTIAL SPEARMAN ARITHMETIC CHECK")
ps_path = fpath("partial_summary.tsv.gz")
pr("Path:", ps_path)
ps = pd.read_csv(ps_path, sep="\t")
pr(ps.to_string(index=False))
for _, r in ps.iterrows():
    nm = int(r["n_marginal_q05"])
    npar = int(r["n_partial_q05"])
    ret_pct = float(r["retained_pct"])
    implied_both = round(nm * ret_pct / 100.0)
    naive_pct = 100.0 * npar / nm if nm else 0.0
    pr(
        f"  {r['axis']:20s}: n_marg={nm:6d}, n_part={npar:6d}, "
        f"n_part/n_marg={naive_pct:5.1f}%, table retained_pct={ret_pct:5.1f}% "
        f"(implies both≈{implied_both}, new_in_partial≈{npar - implied_both})"
    )

# 2. Check PICRUSt2 per_sample_comparison.tsv.gz
section("2. PICRUST2 REFERENCE BIAS AUDIT")
pic_path = fpath("per_sample_comparison.tsv.gz")
pr("Path:", pic_path)
pic = pd.read_csv(pic_path, sep="\t")
pr("Columns:", list(pic.columns))
pr("n rows:", len(pic))
pr(pic.describe().T[["min", "50%", "mean", "max"]].to_string())
pr("Samples with mag_in_ref_not_predicted == 0:", int((pic["mag_in_ref_not_predicted"] == 0).sum()))
for col in pic.columns:
    s = pic[col]
    if pd.api.types.is_numeric_dtype(s):
        pr(f"  {col:30s}: min={s.min():.4f}, q25={s.quantile(0.25):.4f}, med={s.median():.4f}, q75={s.quantile(0.75):.4f}, max={s.max():.4f}")
if "picrust_kos" in pic.columns:
    pr("picrust_kos value counts (top 10):")
    pr(pic["picrust_kos"].value_counts().head(10).to_string())

# 3. Check content_vs_depth.tsv.gz and rarefaction
section("3. SPECIFICITY / CONTENT VS DEPTH / SPIRE ANCHOR CHECK")
cvd = pd.read_csv(fpath("content_vs_depth.tsv.gz"), sep="\t")
pr(cvd.to_string(index=False))
b6 = pd.read_csv(fpath("b6_absence_occupancy_by_length.tsv.gz"), sep="\t")
pr("\nb6_absence_occupancy_by_length.tsv.gz:")
pr(b6.to_string(index=False))

for fn in ["rarefaction_curve.tsv", "b6_absence_vs_occupancy.tsv", "b6_definition_check.tsv", "b6_occupancy_bins.tsv", "b6_length_confound.tsv", "s14_persistence_by_length.tsv"]:
    try:
        p = fpath(fn)
        pr(f"\n--- {fn} ({p}) ---")
        df = pd.read_csv(p, sep="\t")
        pr(df.head(15).to_string(index=False))
    except FileNotFoundError:
        pr(f"Not found: {fn}")

# 4. Check lineage darkness & matched null
section("4. LINEAGE DARKNESS & MATCHED NULL CHECK")
v2_reps = pd.read_csv(fpath("v2_reps.tsv.gz"), sep="\t")
pr("v2_reps shape:", v2_reps.shape)
pr("dark_fraction median:", v2_reps["dark_fraction"].median(), "mean:", v2_reps["dark_fraction"].mean())
pr("n_proteins sum:", v2_reps["n_proteins"].sum(), "n_dark sum:", v2_reps["n_dark"].sum())
pr("pooled dark_proteins / n_proteins:", v2_reps["n_dark"].sum() / v2_reps["n_proteins"].sum())

mn = pd.read_csv(fpath("v2_matched_null.tsv.gz"), sep="\t")
for pool, sub in mn.groupby("pool"):
    diff = sub["target_dark"] - sub["control_dark"]
    w_scipy = stats.wilcoxon(diff)
    n_ties_abs = len(diff) - len(set(abs(diff)))
    pr(
        f"Pool={pool}: n={len(sub)}, target_med={sub['target_dark'].median():.4f}, "
        f"ctrl_med={sub['control_dark'].median():.4f}, diff_med={diff.median():.4f}, "
        f"darker={(diff > 0).sum()}/{len(diff)}, scipy_p={w_scipy.pvalue:.3e}, "
        f"abs_diff_ties={n_ties_abs}"
    )

comp_ph = pd.read_csv(fpath("v2_dark_composition_by_phylum.tsv.gz"), sep="\t")
pr("\nTop phyla by species count in v2_dark_composition_by_phylum:")
pr(comp_ph.head(10).to_string(index=False))

# 5. Check cultured comparison
section("5. CULTURED COMPARISON CHECK")
sp_path = os.path.join(INP, "06_darkmatter/cultured_comparison_clean/results/step3_paired_by_species.tsv.gz")
pr("Path:", sp_path)
sp_pairs = pd.read_csv(sp_path, sep="\t")
pr(sp_pairs.to_string(index=False))
delta = sp_pairs["mag_dark_pct_mean"] - sp_pairs["cult_dark_pct_mean"]
pr("n species:", len(delta))
pr("scipy median(delta):", st.median(delta))
pr("sorted delta[10], delta[11]:", sorted(delta)[10], sorted(delta)[11])
pr("step3_04_summarise upper-middle diffs[len//2]:", sorted(delta)[len(delta) // 2])

pa_path = os.path.join(INP, "06_darkmatter/cultured_comparison_clean/results/pair_anatomy_per_pair.tsv.gz")
pa = pd.read_csv(pa_path, sep="\t")
pr("\npair_anatomy_per_pair.tsv.gz columns:", list(pa.columns))
pr("n pairs:", len(pa))
pr("sum mag_dark_n:", pa["mag_dark_n"].sum(), "sum cult_dark_n:", pa["cult_dark_n"].sum())
pr("sum mag_genes:", pa["mag_genes"].sum(), "sum cult_genes:", pa["cult_genes"].sum())

pr("Unique MAGs in 35 pairs:", pa["mag_id"].nunique() if "mag_id" in pa.columns else pa.iloc[:, 1].nunique())
pr("Unique isolate assemblies in 35 pairs:", pa["cult_id"].nunique() if "cult_id" in pa.columns else pa.iloc[:, 2].nunique())
pr(pa.iloc[:, :8].to_string(index=False))

out_log.close()
