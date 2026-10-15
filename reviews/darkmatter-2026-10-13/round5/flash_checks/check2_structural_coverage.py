#!/usr/bin/env python3
"""Check 2: Exploratory analysis of structural query-score deficits vs alignment coverage.
Investigates whether query-score deficits (TM_query < 0.5) track sub-query alignment coverage
at fixed length/confidence rather than novel fold topology.
"""
import json
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, mannwhitneyu

BASE = Path(__file__).resolve().parents[4]
DM = BASE / "reproducibility" / "darkmatter"

df = pd.read_csv(DM / "structural_raw_audit" / "per_model.tsv", sep="\t")

ret = df[df.raw_rows > 0].copy()
ret["cov_ratio"] = ret["qtm"] / ret["alntm"]

dark_ret = ret[ret.arm == "dark"].copy()
ctrl_ret = ret[ret.arm == "control"].copy()

# 1. Decomposition of qtm < 0.5 among retrieved models
dark_qtm_below = dark_ret[dark_ret.qtm < 0.5]
dark_qtm_below_alntm_high = dark_qtm_below[dark_qtm_below.alntm >= 0.5]
dark_qtm_below_alntm_low = dark_qtm_below[dark_qtm_below.alntm < 0.5]

ctrl_qtm_below = ctrl_ret[ctrl_ret.qtm < 0.5]
ctrl_qtm_below_alntm_high = ctrl_qtm_below[ctrl_qtm_below.alntm >= 0.5]
ctrl_qtm_below_alntm_low = ctrl_qtm_below[ctrl_qtm_below.alntm < 0.5]

# 2. Length comparison between qtm < 0.5 and qtm >= 0.5 in dark arm
dark_qtm_above = dark_ret[dark_ret.qtm >= 0.5]
mwu_len = mannwhitneyu(dark_qtm_below.length, dark_qtm_above.length)

# 3. Correlation between length, pLDDT, and cov_ratio
corr_dark_len_cov, p_len_cov = spearmanr(dark_ret.length, dark_ret.cov_ratio)
corr_dark_plddt_cov, p_plddt_cov = spearmanr(dark_ret.plddt, dark_ret.cov_ratio)

# 4. Length-stratified breakdown of qtm < 0.5 vs alntm >= 0.5
bins = [100, 150, 200, 300, 500, 10000]
labels = ["100-149", "150-199", "200-299", "300-499", "500+"]
dark_ret["len_bin"] = pd.cut(dark_ret.length, bins=bins, labels=labels, right=False)

strata = []
for b, g in dark_ret.groupby("len_bin", observed=True):
    n = len(g)
    q_below = (g.qtm < 0.5).sum()
    aln_high_given_q_below = ((g.qtm < 0.5) & (g.alntm >= 0.5)).sum()
    med_cov = g.cov_ratio.median()
    med_qtm = g.qtm.median()
    med_alntm = g.alntm.median()
    strata.append({
        "len_bin": str(b),
        "n_models": int(n),
        "qtm_below_05": int(q_below),
        "pct_qtm_below": float(q_below / n * 100),
        "qtm_below_with_alntm_ge_05": int(aln_high_given_q_below),
        "pct_aln_high_of_q_below": float(aln_high_given_q_below / q_below * 100) if q_below else 0.0,
        "median_cov_ratio": float(med_cov),
        "median_qtm": float(med_qtm),
        "median_alntm": float(med_alntm)
    })

# 5. True structural outliers (alntm < 0.5 even with raw hits, plus no-hit retrieval failures)
dark_all = df[df.arm == "dark"]
no_rows = dark_all[dark_all.raw_rows == 0]
true_aln_novel = dark_all[(dark_all.raw_rows > 0) & (dark_all.alntm < 0.5)]

results = {
    "total_qualified_dark_models": len(dark_all),
    "dark_retrieved_models": len(dark_ret),
    "dark_retrieval_failures_no_rows": len(no_rows),
    "dark_qtm_below_05": {
        "count": len(dark_qtm_below),
        "alntm_ge_05_count": len(dark_qtm_below_alntm_high),
        "alntm_ge_05_pct": float(len(dark_qtm_below_alntm_high) / len(dark_qtm_below) * 100),
        "alntm_below_05_count": len(dark_qtm_below_alntm_low),
        "alntm_below_05_pct": float(len(dark_qtm_below_alntm_low) / len(dark_qtm_below) * 100),
        "alntm_median": float(dark_qtm_below.alntm.median()),
        "alntm_iqr": [float(dark_qtm_below.alntm.quantile(0.25)), float(dark_qtm_below.alntm.quantile(0.75))],
        "cov_ratio_median": float(dark_qtm_below.cov_ratio.median()),
        "cov_ratio_iqr": [float(dark_qtm_below.cov_ratio.quantile(0.25)), float(dark_qtm_below.cov_ratio.quantile(0.75))]
    },
    "dark_length_contrast": {
        "qtm_below_median_length": float(dark_qtm_below.length.median()),
        "qtm_below_iqr_length": [float(dark_qtm_below.length.quantile(0.25)), float(dark_qtm_below.length.quantile(0.75))],
        "qtm_above_median_length": float(dark_qtm_above.length.median()),
        "qtm_above_iqr_length": [float(dark_qtm_above.length.quantile(0.25)), float(dark_qtm_above.length.quantile(0.75))],
        "mann_whitney_u": float(mwu_len.statistic),
        "mann_whitney_pvalue": float(mwu_len.pvalue)
    },
    "correlation_spearman": {
        "length_vs_cov_ratio_rho": float(corr_dark_len_cov),
        "length_vs_cov_ratio_pvalue": float(p_len_cov),
        "plddt_vs_cov_ratio_rho": float(corr_dark_plddt_cov),
        "plddt_vs_cov_ratio_pvalue": float(p_plddt_cov)
    },
    "strata_by_query_length": strata,
    "true_novel_fold_candidates": {
        "retrieval_failure_count": len(no_rows),
        "alntm_below_05_count": len(true_aln_novel),
        "total_conservative_fold_novelty_upper_bound": len(no_rows) + len(true_aln_novel),
        "pct_of_all_qualified_dark_models": float((len(no_rows) + len(true_aln_novel)) / len(dark_all) * 100)
    }
}

out_dir = Path(__file__).resolve().parent
with open(out_dir / "check2_structural_coverage.json", "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
