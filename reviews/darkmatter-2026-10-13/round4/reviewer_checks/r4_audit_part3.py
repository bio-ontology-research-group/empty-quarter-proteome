#!/usr/bin/env python3
"""Part 3 of Round 4 computational audit: Lineage darkness, Thermoplasmatota, Structure, Reference bias."""
import csv
import gzip
import json
import math
import os
import statistics as st
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

ROOT = "/home/leechuck/Documents/papers/empty-quarter-mags"
DM = os.path.join(ROOT, "reproducibility/darkmatter")
INP = os.path.join(DM, "inputs")
SNAP = os.path.join(DM, "source_snapshot")
OUT_DIR = os.path.join(ROOT, "reviews/darkmatter-2026-10-13/round4/reviewer_checks")

out_log = open(os.path.join(OUT_DIR, "r4_audit_part3.txt"), "w")


def pr(*args, **kwargs):
    print(*args, **kwargs)
    print(*args, **kwargs, file=out_log)


def section(title):
    pr("\n" + "=" * 78)
    pr(title)
    pr("=" * 78)


# 1. Lineage darkness regression & novelty test
section("1. LINEAGE DARKNESS REGRESSION & NOVELTY TEST")
d = pd.read_csv(
    os.path.join(INP, "06_darkmatter/lineage_darkness_clean/results/v2_reps.tsv.gz"),
    sep="\t",
)
pr("v2_reps n:", len(d))
pr("dark_fraction median:", d["dark_fraction"].median(), "IQR:", d["dark_fraction"].quantile(0.25), d["dark_fraction"].quantile(0.75))

# Corrected genome size regression: length / (completeness / 100)
d["corr_size"] = d["length"] / (d["completeness"] / 100.0)
d["log10_corr_size"] = np.log10(d["corr_size"])
m_sz = smf.ols("dark_fraction ~ log10_corr_size", data=d).fit()
pr(f"Corrected size OLS: slope={m_sz.params['log10_corr_size']:.5f}, R2={m_sz.rsquared:.5f}")

# Check Model A, B, C, L1, L2, L3
Q = "completeness + contamination + log10_N50 + log10_length"
mA = smf.ols("dark_fraction ~ " + Q, data=d).fit()
mB = smf.ols("dark_fraction ~ " + Q + " + C(phylum)", data=d).fit()
mC = smf.ols("dark_fraction ~ " + Q + " + C(phylum) + C(genus_grp)", data=d).fit()
mL1 = smf.ols("dark_fraction ~ C(phylum)", data=d).fit()
mL2 = smf.ols("dark_fraction ~ C(phylum) + C(genus_grp)", data=d).fit()
mL3 = smf.ols("dark_fraction ~ C(phylum) + C(genus_grp) + " + Q, data=d).fit()
pr(f"mA (quality): R2={mA.rsquared:.4f}, adjR2={mA.rsquared_adj:.4f}, df_model={mA.df_model}")
pr(f"mB (A+phylum): R2={mB.rsquared:.4f}, adjR2={mB.rsquared_adj:.4f}, df_model={mB.df_model}")
pr(f"mC (B+genus): R2={mC.rsquared:.4f}, adjR2={mC.rsquared_adj:.4f}, df_model={mC.df_model}")
pr(f"Unique to lineage (mL3 - mA): {mL3.rsquared - mA.rsquared:.4f}, Unique to quality (mL3 - mL2): {mL3.rsquared - mL2.rsquared:.4f}")

# Check novelty ladder
nA = smf.ols("dark_fraction ~ " + Q + " + first_unnamed_rank", data=d).fit()
nB = smf.ols("dark_fraction ~ " + Q + " + first_unnamed_rank + C(phylum)", data=d).fit()
pr(f"Novelty nA R2={nA.rsquared:.4f}, nB R2={nB.rsquared:.4f}, nB-nA={nB.rsquared - nA.rsquared:.4f}, ratio={(nB.rsquared - nA.rsquared)/(mB.rsquared - mA.rsquared):.4f}")

# First unnamed rank medians
for k, sub in d.groupby("first_unnamed_rank"):
    pr(f"  first_unnamed_rank={k}: n={len(sub)}, median_dark={sub['dark_fraction'].median():.4f}")

# Within-phylum named vs placeholder genus
for ph in ["Actinomycetota", "Pseudomonadota", "Bacteroidota", "Thermoplasmatota"]:
    sub = d[d["phylum"] == ph]
    for gc in ["named", "placeholder", "blank"]:
        s2 = sub[sub["genus_class"] == gc]
        if len(s2):
            pr(f"  {ph:18s} {gc:12s}: n={len(s2):4d}, median={s2['dark_fraction'].median():.4f}")

# Phylum residuals
ph_dark = pd.read_csv(
    os.path.join(INP, "06_darkmatter/lineage_darkness_clean/results/v2_phylum_dark.tsv.gz"),
    sep="\t",
)
pr("\nPhyla with n_reps >= 20 sorted by median_resid_A:")
pr(ph_dark[ph_dark["n_reps"] >= 20].sort_values("median_resid_A", ascending=False).to_string(index=False))

# Check Thermoplasmatota vs other archaea and bacteria
th = d[d["phylum"] == "Thermoplasmatota"]
ar_oth = d[(d["domain"] == "Archaea") & (d["phylum"] != "Thermoplasmatota")]
bac = d[d["domain"] == "Bacteria"]
bac_p95 = bac["dark_fraction"].quantile(0.95)
pr(f"\nThermoplasmatota n={len(th)}, med={th['dark_fraction'].median():.4f}")
pr(f"Other Archaea n={len(ar_oth)}, med={ar_oth['dark_fraction'].median():.4f}")
pr(f"Bacteria n={len(bac)}, med={bac['dark_fraction'].median():.4f}, p95={bac_p95:.4f}")
pr(f"Thermo > bac p95: {(th['dark_fraction'] > bac_p95).sum()}/{len(th)}")

# First unnamed rank == 2 (class)
th_r2 = th[th["first_unnamed_rank"] == 2]
oth_r2 = d[(d["first_unnamed_rank"] == 2) & (d["phylum"] != "Thermoplasmatota")]
pr(f"Rank=2 (class): Thermo n={len(th_r2)} med={th_r2['dark_fraction'].median():.4f} vs Other n={len(oth_r2)} med={oth_r2['dark_fraction'].median():.4f}")

# 2. Thermoplasmatota ecology & metabolism
section("2. THERMOPLASMATOTA ECOLOGY & METABOLISM CHECK")
te = pd.read_csv(
    os.path.join(INP, "06_darkmatter/lineage_darkness_clean/results/v2_thermo_ecology_summary.tsv.gz"),
    sep="\t",
)
pr("v2_thermo_ecology_summary.tsv.gz (selected rows):")
pr(te[te["group"].isin(["Thermoplasmatota", "Halobacteriota", "all"])].to_string(index=False))

tm_rep_path = os.path.join(SNAP, "06_darkmatter/lineage_darkness_clean/results/v2_thermo_metabolism_report.txt")
if os.path.exists(tm_rep_path):
    pr("\n--- v2_thermo_metabolism_report.txt (full) ---")
    with open(tm_rep_path) as fh:
        pr(fh.read())

tm_sum_path = os.path.join(SNAP, "06_darkmatter/lineage_darkness_clean/results/v2_thermo_metabolism_summary.tsv")
if os.path.exists(tm_sum_path):
    tm_sum = pd.read_csv(tm_sum_path, sep="\t")
    # Check specific features cited in lines 286-298 of submitted_main.tex:
    # glycolysis (82%), TCA cycle through 2-oxoglutarate:ferredoxin oxidoreductase (96%),
    # glyoxylate cycle (79% vs 35% in matched Halobacteriota), complex I and II,
    # cytochrome c oxidase (59%), archaeal ATPase (91%),
    # branched-chain amino-acid transporters (79%), sugar transporters (2% vs 73% in matched Halobacteriota),
    # extracellular serine protease (96%), endoglucanase (79%), chitinase (44%),
    # ilv_carbon_skeleton (2% vs 91%), bcaa_aminotransferase (93%),
    # catalase (29% vs 85%), archaellum (100%), CRISPR-Cas (79%)
    focus_feats = [
        "M00001", "M00002", "M00009", "M00011", "M00012",
        "oxoglutarate_ferredoxin_OR", "complexI_nuo_half", "complexII_sdh_2sub",
        "complexIV_coxA", "complexIV_coxA_coxB", "ATPase_VA_half",
        "branched_chain_aa_ABC", "sugar_ABC_msm",
        "ilv_carbon_skeleton", "bcaa_aminotransferase",
        "catalase", "archaellum_fla", "crispr_any", "crispr_cas1_cas2",
        "COF_cobalamin", "COF_thiamine", "COF_F420"
    ]
    sub_tm = tm_sum[tm_sum["feature"].isin(focus_feats) & tm_sum["group"].isin(["thermo_237", "halo_matched", "archaea_other_matched", "bacteria_matched", "halo_all", "bacteria_all"])]
    pr("\nSelected rows from v2_thermo_metabolism_summary.tsv:")
    pr(sub_tm[["feature", "feature_name", "group", "n", "n_present", "share"]].to_string(index=False))

out_log.close()
