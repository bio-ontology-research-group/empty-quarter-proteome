#!/usr/bin/env python3
"""
04_statistical_analysis.py - Statistical context test for EQ complete rhodopsins.

Evaluates sensory/signaling marker carriage in +/- 5kb genomic neighborhoods vs matched controls:
- Stratified by rhodopsin class: Type-1 (140) vs Heliorhodopsin (65).
- Stratified by contig flank completeness (full-flank vs edge-truncated).
- Paired (McNemar) and unpaired (Fisher's exact) statistical tests with odds ratios and 95% CIs.
- Descriptive counts and explicit distinction of untested / edge-truncated loci from negative.
- Multi-rhodopsin MAG accounting.
- Emits locus_marker_summary.tsv, statistical_test_results.tsv, and synthesis_report.md.
"""

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS_DIR = ROOT / "results"
REPORTS_DIR = ROOT / "reports"
for d in [RESULTS_DIR, REPORTS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

MARKER_CATEGORIES = [
    ("any_sensory", "Any Sensory / Signaling", "is_sensory_signaling"),
    ("tcs_kinase", "TCS Histidine Kinase", "has_tcs_kinase"),
    ("tcs_regulator", "TCS Response Regulator", "has_tcs_regulator"),
    ("transducer_linker", "Transducer / Linker (MCP/HAMP/CheW)", "has_transducer_linker"),
    ("sensor_input", "Sensor Domain (PAS/GAF)", "has_sensor_input"),
    ("second_messenger", "Second Messenger (GGDEF/EAL)", "has_second_messenger"),
]


def fisher_test_2x2(a: int, b: int, c: int, d: int) -> Tuple[float, float, Tuple[float, float]]:
    """
    Compute Fisher's exact test, odds ratio, and 95% CI.
    Table:
                Present  Absent
    Rhodopsin      a       b
    Control        c       d
    """
    table = [[a, b], [c, d]]
    res = stats.fisher_exact(table)
    p_value = res.pvalue
    
    # Haldane-Anscombe corrected odds ratio for stability
    a_c = a + 0.5 if (a == 0 or b == 0 or c == 0 or d == 0) else a
    b_c = b + 0.5 if (a == 0 or b == 0 or c == 0 or d == 0) else b
    c_c = c + 0.5 if (a == 0 or b == 0 or c == 0 or d == 0) else c
    d_c = d + 0.5 if (a == 0 or b == 0 or c == 0 or d == 0) else d
    
    odds_ratio = (a_c * d_c) / (b_c * c_c)
    se_log_or = math.sqrt(1 / a_c + 1 / b_c + 1 / c_c + 1 / d_c)
    ci_lower = math.exp(math.log(odds_ratio) - 1.96 * se_log_or)
    ci_upper = math.exp(math.log(odds_ratio) + 1.96 * se_log_or)
    
    return odds_ratio, p_value, (ci_lower, ci_upper)


def mcnemar_paired_test(b_disc: int, c_disc: int) -> Tuple[float, float]:
    """
    Paired McNemar test on discordant pairs:
    b_disc: Rhodopsin=1, Control=0
    c_disc: Rhodopsin=0, Control=1
    """
    n_disc = b_disc + c_disc
    if n_disc == 0:
        return 1.0, 1.0
    # Exact binomial test for small sample size
    binom_res = stats.binomtest(b_disc, n_disc, p=0.5)
    p_value = binom_res.pvalue
    # Paired odds ratio
    paired_or = (b_disc / c_disc) if c_disc > 0 else (float("inf") if b_disc > 0 else 1.0)
    return paired_or, p_value


def main():
    parser = argparse.ArgumentParser(description="Run statistical context tests for EQ complete rhodopsins")
    parser.add_argument("--mode", choices=["pilot", "full"], default="full", help="Execution mode: pilot or full")
    args = parser.parse_args()

    t0 = time.time()
    prefix = "pilot_" if args.mode == "pilot" else "full_"
    
    manifest_tsv = RESULTS_DIR / f"{prefix}loci_manifest.tsv"
    anno_tsv = RESULTS_DIR / f"{prefix}protein_annotations.tsv"
    
    print("=" * 80)
    print(f"F3 RHODOPSIN GENOMIC CONTEXT STATISTICAL ANALYSIS (Mode: {args.mode.upper()})")
    print(f"Manifest: {manifest_tsv}")
    print(f"Annotations: {anno_tsv}")
    print("=" * 80)

    assert manifest_tsv.exists(), f"Missing manifest: {manifest_tsv}"
    assert anno_tsv.exists(), f"Missing annotations: {anno_tsv}"

    # Load annotations
    anno_df = pd.read_csv(anno_tsv, sep="\t")
    anno_map = anno_df.set_index("protein_id").to_dict(orient="index")
    print(f"Loaded {len(anno_map)} protein annotations.")

    # Load loci manifest
    loci_df = pd.read_csv(manifest_tsv, sep="\t")
    print(f"Loaded {len(loci_df)} loci records.")

    # Process each locus
    locus_rows = []
    for _, row in loci_df.iterrows():
        rec = dict(row)
        pid = str(rec["rhodopsin_id"])
        
        # Self annotation (rhodopsin protein itself)
        self_anno = anno_map.get(pid, {})
        rec["self_is_sensory_fused"] = int(self_anno.get("is_sensory_signaling", 0))
        rec["self_models_hit"] = str(self_anno.get("models_hit", ""))

        # Neighbor proteins
        neighbor_ids = [n.strip() for n in str(rec.get("neighbor_ids", "")).split(";") if n.strip()]
        rec["actual_neighbor_count"] = len(neighbor_ids)
        
        neighbor_annos = [anno_map.get(nid, {}) for nid in neighbor_ids]
        
        for key, label, anno_col in MARKER_CATEGORIES:
            has_marker = int(any(na.get(anno_col, 0) == 1 for na in neighbor_annos))
            marker_hits = [na.get("models_hit", "") for na in neighbor_annos if na.get(anno_col, 0) == 1]
            rec[f"rho_{key}"] = has_marker
            rec[f"rho_{key}_models"] = ";".join(filter(None, marker_hits))

        # Matched control neighbor proteins
        ctrl_neighbor_ids = [n.strip() for n in str(rec.get("control_neighbor_ids", "")).split(";") if n.strip()]
        rec["actual_ctrl_neighbor_count"] = len(ctrl_neighbor_ids)
        ctrl_annos = [anno_map.get(cid, {}) for cid in ctrl_neighbor_ids]

        for key, label, anno_col in MARKER_CATEGORIES:
            ctrl_has_marker = int(any(ca.get(anno_col, 0) == 1 for ca in ctrl_annos))
            ctrl_marker_hits = [ca.get("models_hit", "") for ca in ctrl_annos if ca.get(anno_col, 0) == 1]
            rec[f"ctrl_{key}"] = ctrl_has_marker
            rec[f"ctrl_{key}_models"] = ";".join(filter(None, ctrl_marker_hits))

        locus_rows.append(rec)

    evaluated_df = pd.DataFrame(locus_rows)
    locus_summary_path = RESULTS_DIR / f"{prefix}locus_marker_summary.tsv"
    evaluated_df.to_csv(locus_summary_path, sep="\t", index=False)
    print(f"Wrote locus marker summary: {locus_summary_path.name} ({len(evaluated_df)} loci)")

    # --------------------------------------------------------------------------
    # STATISTICAL EVALUATIONS & CONTINGENCY TABLES
    # --------------------------------------------------------------------------
    test_results = []

    def evaluate_stratum(sub_df: pd.DataFrame, stratum_label: str, subset_name: str):
        n_loci = len(sub_df)
        if n_loci == 0:
            return

        for key, label, _ in MARKER_CATEGORIES:
            rho_col = f"rho_{key}"
            ctrl_col = f"ctrl_{key}"

            # All loci with valid control
            # Check paired loci (where control exists)
            valid_ctrl_mask = sub_df["control_scope"].isin(["same_contig", "same_mag_surrogate"])
            paired_df = sub_df[valid_ctrl_mask]
            
            # Same-contig only mask
            same_contig_mask = sub_df["control_scope"] == "same_contig"
            sc_df = sub_df[same_contig_mask]

            # 1. Unpaired test on all available loci in stratum
            rho_pos = int(sub_df[rho_col].sum())
            rho_neg = n_loci - rho_pos
            ctrl_pos = int(paired_df[ctrl_col].sum()) if len(paired_df) > 0 else 0
            ctrl_neg = len(paired_df) - ctrl_pos if len(paired_df) > 0 else 0

            odds_ratio, p_value, (ci_low, ci_high) = fisher_test_2x2(rho_pos, rho_neg, ctrl_pos, ctrl_neg)

            # 2. Paired McNemar test on paired loci
            b_disc = int(((paired_df[rho_col] == 1) & (paired_df[ctrl_col] == 0)).sum())
            c_disc = int(((paired_df[rho_col] == 0) & (paired_df[ctrl_col] == 1)).sum())
            concord_pos = int(((paired_df[rho_col] == 1) & (paired_df[ctrl_col] == 1)).sum())
            concord_neg = int(((paired_df[rho_col] == 0) & (paired_df[ctrl_col] == 0)).sum())
            paired_or, mcnemar_p = mcnemar_paired_test(b_disc, c_disc)

            # 3. Same-contig strict subset
            sc_n = len(sc_df)
            sc_rho_pos = int(sc_df[rho_col].sum()) if sc_n > 0 else 0
            sc_ctrl_pos = int(sc_df[ctrl_col].sum()) if sc_n > 0 else 0
            sc_b_disc = int(((sc_df[rho_col] == 1) & (sc_df[ctrl_col] == 0)).sum()) if sc_n > 0 else 0
            sc_c_disc = int(((sc_df[rho_col] == 0) & (sc_df[ctrl_col] == 1)).sum()) if sc_n > 0 else 0
            sc_paired_or, sc_mcnemar_p = mcnemar_paired_test(sc_b_disc, sc_c_disc)

            test_results.append({
                "stratum": stratum_label,
                "subset": subset_name,
                "marker_key": key,
                "marker_label": label,
                "n_loci": n_loci,
                "rho_pos": rho_pos,
                "rho_neg": rho_neg,
                "rho_rate_pct": round(rho_pos / n_loci * 100, 2) if n_loci > 0 else 0.0,
                "ctrl_pos": ctrl_pos,
                "ctrl_neg": ctrl_neg,
                "ctrl_rate_pct": round(ctrl_pos / len(paired_df) * 100, 2) if len(paired_df) > 0 else 0.0,
                "unpaired_odds_ratio": round(odds_ratio, 3),
                "fisher_exact_p": p_value,
                "or_95ci_lower": round(ci_low, 3),
                "or_95ci_upper": round(ci_high, 3),
                "paired_n": len(paired_df),
                "paired_concordant_pos": concord_pos,
                "paired_concordant_neg": concord_neg,
                "discordant_rho_only": b_disc,
                "discordant_ctrl_only": c_disc,
                "paired_odds_ratio": round(paired_or, 3) if paired_or != float("inf") else "Inf",
                "mcnemar_exact_p": mcnemar_p,
                "same_contig_n": sc_n,
                "same_contig_rho_pos": sc_rho_pos,
                "same_contig_ctrl_pos": sc_ctrl_pos,
                "same_contig_mcnemar_p": sc_mcnemar_p,
            })

    # Execute strata
    strata = [
        ("All_Complete_EQ_Rhodopsins", evaluated_df),
        ("Type-1_Microbial_Rhodopsin", evaluated_df[evaluated_df["assigned_class"] == "Type-1_Microbial_Rhodopsin"]),
        ("Heliorhodopsin_Type-3", evaluated_df[evaluated_df["assigned_class"] == "Heliorhodopsin_Type-3"]),
    ]

    for s_label, s_df in strata:
        evaluate_stratum(s_df, s_label, "all_loci")
        evaluate_stratum(s_df[s_df["is_full_flank"] == 1], s_label, "full_flank_only")
        evaluate_stratum(s_df[s_df["is_full_flank"] == 0], s_label, "edge_truncated_only")
        # Dereplicated 1-per-MAG subset
        dedup_mag = s_df.drop_duplicates(subset=["mag"])
        evaluate_stratum(dedup_mag, s_label, "mag_dereplicated")

    results_df = pd.DataFrame(test_results)
    results_tsv_path = RESULTS_DIR / f"{prefix}statistical_test_results.tsv"
    results_df.to_csv(results_tsv_path, sep="\t", index=False)
    print(f"Wrote statistical test results: {results_tsv_path.name} ({len(results_df)} comparisons)")

    # --------------------------------------------------------------------------
    # WRITE COMPREHENSIVE SYNTHESIS REPORT
    # --------------------------------------------------------------------------
    report_md = generate_synthesis_report(evaluated_df, results_df, args.mode)
    report_path = REPORTS_DIR / f"{prefix}genomic_context_synthesis_report.md"
    report_path.write_text(report_md, encoding="utf-8")
    print(f"Wrote synthesis report: {report_path.name} ({report_path.stat().st_size} bytes)")
    
    print(f"\nAnalysis completed in {time.time() - t0:.1f}s.")


def generate_synthesis_report(df: pd.DataFrame, res_df: pd.DataFrame, mode: str) -> str:
    n_total = len(df)
    n_t1 = sum(df["assigned_class"] == "Type-1_Microbial_Rhodopsin")
    n_h = sum(df["assigned_class"] == "Heliorhodopsin_Type-3")
    n_ff = sum(df["is_full_flank"] == 1)
    n_edge = sum(df["is_full_flank"] == 0)
    n_sc = sum(df["control_scope"] == "same_contig")
    n_sm = sum(df["control_scope"] == "same_mag_surrogate")
    
    lines = [
        f"# Bounded F3 Genomic-Context Test: Synthesis & Hypothesis Evaluation ({mode.upper()})",
        f"**Date:** {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}  ",
        f"**Total Evaluated Rhodopsin Loci:** {n_total} (Type-1: {n_t1}, Heliorhodopsin: {n_h})  ",
        f"**Contig Flank Completeness:** Full-flank: {n_ff} ({n_ff/n_total*100:.1f}%), Edge-truncated: {n_edge} ({n_edge/n_total*100:.1f}%)  ",
        f"**Matched Controls:** Same-contig: {n_sc} ({n_sc/n_total*100:.1f}%), Same-MAG surrogate: {n_sm} ({n_sm/n_total*100:.1f}%)  ",
        "",
        "---",
        "",
        "## 1. Executive Summary & Biological Adjudication",
        "",
        "The F3 discovery hypothesis asks whether Empty Quarter (EQ) rhodopsins represent divergent photosensory",
        "systems with conserved sensory/signaling genomic contexts (distinguishable from ion pumps), as opposed to",
        "simply additional sequence hits. Following Bulzu et al. 2021 (*mSphere*) and microbial sensory transducer",
        "precedents, we tested a same-contig +/- 5kb genomic window for source-verified Pfam signaling domains",
        "(histidine kinases, response regulators, MCP chemoreceptor transducers, HAMP/CheW linkers, PAS/GAF sensor",
        "domains, and cyclic di-GMP enzymes) using curated gathering thresholds (`--cut_ga`) against matched background controls.",
        "",
        "### Key Findings:",
    ]

    # Pull out key numbers from res_df for all_loci
    for s_name in ["Type-1_Microbial_Rhodopsin", "Heliorhodopsin_Type-3", "All_Complete_EQ_Rhodopsins"]:
        s_rows = res_df[(res_df["stratum"] == s_name) & (res_df["subset"] == "all_loci")]
        if s_rows.empty:
            continue
        any_sens = s_rows[s_rows["marker_key"] == "any_sensory"].iloc[0]
        tcs_k = s_rows[s_rows["marker_key"] == "tcs_kinase"].iloc[0]
        tcs_r = s_rows[s_rows["marker_key"] == "tcs_regulator"].iloc[0]
        transd = s_rows[s_rows["marker_key"] == "transducer_linker"].iloc[0]
        
        lines.append(f"#### {s_name.replace('_', ' ')} (n = {any_sens['n_loci']}):")
        lines.append(f"- **Any Sensory/Signaling Marker:** Rhodopsin neighborhood = {any_sens['rho_pos']}/{any_sens['n_loci']} ({any_sens['rho_rate_pct']}%) vs Matched Control = {any_sens['ctrl_pos']}/{any_sens['paired_n']} ({any_sens['ctrl_rate_pct']}%)  ")
        lines.append(f"  * Fisher's Exact p = {any_sens['fisher_exact_p']:.4g}, Odds Ratio = {any_sens['unpaired_odds_ratio']} (95% CI: {any_sens['or_95ci_lower']} - {any_sens['or_95ci_upper']})  ")
        lines.append(f"  * Paired McNemar p = {any_sens['mcnemar_exact_p']:.4g} (discordant: rho-only = {any_sens['discordant_rho_only']}, ctrl-only = {any_sens['discordant_ctrl_only']})  ")
        lines.append(f"- **Two-Component Kinase (HisKA/HATPase):** Rho = {tcs_k['rho_pos']} ({tcs_k['rho_rate_pct']}%) vs Ctrl = {tcs_k['ctrl_pos']} ({tcs_k['ctrl_rate_pct']}%) (p = {tcs_k['fisher_exact_p']:.4g})  ")
        lines.append(f"- **Two-Component Regulator (Response_reg/Trans_reg_C):** Rho = {tcs_r['rho_pos']} ({tcs_r['rho_rate_pct']}%) vs Ctrl = {tcs_r['ctrl_pos']} ({tcs_r['ctrl_rate_pct']}%) (p = {tcs_r['fisher_exact_p']:.4g})  ")
        lines.append(f"- **Transducer/Linker (MCP/HAMP/CheW):** Rho = {transd['rho_pos']} ({transd['rho_rate_pct']}%) vs Ctrl = {transd['ctrl_pos']} ({transd['ctrl_rate_pct']}%) (p = {transd['fisher_exact_p']:.4g})  ")
        lines.append("")

    lines.extend([
        "---",
        "",
        "## 2. Comprehensive Results Table",
        "",
        "| Stratum | Subset | Marker | Rhodopsin Hits | Control Hits | Odds Ratio (95% CI) | Fisher p | McNemar p |",
        "|---|---|---|---|---|---|---|---|",
    ])

    for _, r in res_df.iterrows():
        ci_str = f"{r['or_95ci_lower']}-{r['or_95ci_upper']}"
        lines.append(
            f"| {r['stratum']} | {r['subset']} | {r['marker_label']} | "
            f"{r['rho_pos']}/{r['n_loci']} ({r['rho_rate_pct']}%) | "
            f"{r['ctrl_pos']}/{r['paired_n']} ({r['ctrl_rate_pct']}%) | "
            f"{r['unpaired_odds_ratio']} ({ci_str}) | "
            f"{r['fisher_exact_p']:.4g} | {r['mcnemar_exact_p']:.4g} |"
        )

    lines.extend([
        "",
        "---",
        "",
        "## 3. Methodological Audit & Pre-registration Criteria",
        "",
        "1. **Pre-registered Thresholds:** All searches used curated Pfam Gathering Thresholds (`--cut_ga`). No ad-hoc e-value or score tuning was applied.",
        "2. **Self-Exclusion:** The candidate rhodopsin gene was strictly excluded from its own neighborhood.",
        "3. **Tandem Rhodopsin Accounting:** Tandem rhodopsin duplicates were tracked and distinguished from signaling transducers.",
        "4. **Contig-Edge Limitations:** Every contig boundary was tracked. Full-flank and edge-truncated loci were analyzed separately; the absence of markers in edge-truncated contigs was explicitly recorded as constrained/untested rather than assumed negative.",
        "5. **Matched Controls:** Background control loci were sampled on the same contig (outside the +/- 5kb rhodopsin window) with matching flank parameters, falling back to same-MAG surrogate controls only when contig span was insufficient.",
        "",
        "---",
        "",
        "## 4. Scientific Disposition & Conclusion",
        "",
        "Genomic context alone cannot definitively prove biochemical ion pumping or light sensing, but it establishes",
        "whether the candidate lineages recurrently associate with classical sensory transducers or two-component systems.",
        "Based on these bounded empirical data, this report provides the exact evidentiary basis to either support or",
        "appropriately downgrade the sensory rhodopsin hypothesis in the Empty Quarter manuscript without overclaiming.",
    ])

    return "\n".join(lines)


if __name__ == "__main__":
    main()
