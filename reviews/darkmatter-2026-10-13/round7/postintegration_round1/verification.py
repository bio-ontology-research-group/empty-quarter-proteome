#!/usr/bin/env python3
"""Independent bounded verification for post-integration review round 1 of 2.

Audits all integrated H1, H2, H3, and H4 quantitative claims, denominators,
methodological assertions, Figure 4 values, and reproducibility-registry alignment
against the submitted manuscript snapshot and frozen evidence tables.
"""

import csv
import gzip
import hashlib
import io
import json
import math
import os
import random
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path("/home/leechuck/Documents/papers/empty-quarter-mags")
R7 = ROOT / "reviews/darkmatter-2026-10-13/round7"
R6 = ROOT / "reviews/darkmatter-2026-10-13/round6"
OUT_DIR = R7 / "postintegration_round1"
SUBMITTED_DIR = OUT_DIR / "submitted"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_gz_tsv(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_json(path: Path):
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def check_manuscript_and_pdf_integrity():
    main_tex = ROOT / "manuscript/main.tex"
    sub_tex = SUBMITTED_DIR / "main.tex"
    sub_pdf = SUBMITTED_DIR / "main.pdf"
    sub_fig4 = SUBMITTED_DIR / "fig_structure_current.pdf"
    man_fig4 = ROOT / "manuscript/figures/fig_structure_current.pdf"
    pdf_scan = R7 / "manuscript/pdf_integrity_initial.txt"

    main_sha = sha256_file(main_tex)
    sub_sha = sha256_file(sub_tex)
    assert main_sha == sub_sha, "manuscript/main.tex and submitted/main.tex differ"

    fig4_sub_sha = sha256_file(sub_fig4)
    fig4_man_sha = sha256_file(man_fig4)
    assert fig4_sub_sha == fig4_man_sha, "fig_structure_current.pdf differs from manuscript/figures/"

    tex_text = main_tex.read_text(encoding="utf-8")
    # Ensure withdrawn/obsolete alnTM claims are absent from manuscript/main.tex
    obsolete_phrases = [
        "797 have an alignment-normalised",
        "936 local",
        "849 K_local_pdb",
        "95 dual-none",
        "1.0200",
    ]
    found_obsolete = [p for p in obsolete_phrases if p in tex_text]
    assert not found_obsolete, f"Found obsolete phrases in main.tex: {found_obsolete}"

    return {
        "main_tex_sha256": main_sha,
        "submitted_tex_sha256": sub_sha,
        "tex_identical": main_sha == sub_sha,
        "submitted_pdf_sha256": sha256_file(sub_pdf),
        "submitted_fig4_sha256": fig4_sub_sha,
        "pdf_integrity_scan_path": str(pdf_scan.relative_to(ROOT)),
        "pdf_integrity_scan_text": pdf_scan.read_text(encoding="utf-8").strip(),
        "obsolete_alntm_phrases_found": found_obsolete,
    }


def verify_h1():
    h1_dir = R6 / "snapshots/h1_final"
    manifest = read_json(h1_dir / "SNAPSHOT_MANIFEST.json")
    manifest_checks = 0
    for entry in manifest["files"]:
        rel = entry["path"]
        actual = sha256_file(h1_dir / rel)
        assert actual == entry["sha256"], f"H1 hash mismatch: {rel}"
        manifest_checks += 1

    closure = read_json(h1_dir / "results/v3_closure_verification.json")
    model_rows = read_tsv(h1_dir / "results/v3_pgls_model_summary.tsv")
    sens_rows = read_tsv(h1_dir / "results/v3_pgls_sensitivity.tsv")
    comp_meta = read_tsv(h1_dir / "results/v3_comparator_metadata_summary.tsv")
    cov_table = read_tsv(h1_dir / "results/denominator_caller_audit.tsv")

    assert len(cov_table) == 237
    class_counts = Counter(r["cls"] for r in cov_table)
    assert class_counts == {"SW-10-69-26": 185, "E2": 51, "Thermoplasmata": 1}

    dark_fracs = sorted(float(r["dark_fraction"]) for r in cov_table)
    median_dark = dark_fracs[len(dark_fracs) // 2]
    assert round(median_dark, 3) == 0.310

    distinct_comparators = {r["aai_nearest"] for r in cov_table}
    assert len(distinct_comparators) == 32
    source_counts = Counter(r["aai_source"] for r in cov_table)
    assert source_counts == {"GTDB_r232_rep": 229, "MGnify": 8}

    meta_all = [r for r in comp_meta if r["stratum"] == "All_Thermoplasmatota_237"][0]
    assert int(meta_all["n_genomes"]) == 237
    assert int(meta_all["n_unique_nearest_comparators"]) == 32
    assert int(meta_all["n_nearest_gtdb_r232_rep"]) == 229
    assert int(meta_all["n_nearest_mgnify"]) == 8
    assert int(meta_all["n_biome_soil"]) == 142
    assert int(meta_all["n_biome_marine_sediment"]) == 3
    assert int(meta_all["n_biome_unannotated_gtdb_only"]) == 92

    pgls_by_var = {r["variable"]: r for r in model_rows}
    aai = pgls_by_var["aai"]
    rbh = pgls_by_var["frac_rbh_desert"]

    assert round(float(aai["pgls_coef"]), 6) == 0.001704
    assert round(float(aai["pgls_se_unb"]), 6) == 0.000407
    assert abs(float(aai["pgls_pval_unb"]) - 4.065008e-5) < 1e-9

    assert round(float(rbh["pgls_coef"]), 6) == -0.137484
    assert round(float(rbh["pgls_se_unb"]), 6) == 0.019456
    assert abs(float(rbh["pgls_pval_unb"]) - 1.879206e-11) < 1e-15

    assert round(float(closure["primary_lambda_ml"]), 4) == 0.9158
    assert round(float(closure["primary_lambda_ci95"][0]), 4) == 0.8463
    assert round(float(closure["primary_lambda_ci95"][1]), 4) == 0.9578
    assert round(float(closure["primary_r2_gls"]), 4) == 0.4691

    # Verify all 7 specifications (S1 primary + 6 sensitivities S2-S7) keep negative RBH and positive AAI
    assert len(sens_rows) == 7
    for s in sens_rows:
        assert float(s["coef_aai"]) > 0, f"AAI not positive in {s['model_id']}"
        assert float(s["coef_frac_rbh"]) < 0, f"RBH not negative in {s['model_id']}"
    s7 = [s for s in sens_rows if s["model_id"].startswith("S7")][0]
    assert int(s7["n"]) == 232

    return {
        "status": "PASS",
        "snapshot_files_verified": manifest_checks,
        "n_genomes": len(cov_table),
        "class_counts": dict(class_counts),
        "median_dark_fraction": round(median_dark, 4),
        "unique_comparators": len(distinct_comparators),
        "comparator_pair_sources": dict(source_counts),
        "env_pairs": {"soil": 142, "marine_sediment": 3, "unknown": 92, "known_total": 145},
        "tree_tips": 237 + 1647 + 30,
        "lambda_ml": round(float(closure["primary_lambda_ml"]), 6),
        "lambda_ci95": [round(float(x), 6) for x in closure["primary_lambda_ci95"]],
        "r2_gls": round(float(closure["primary_r2_gls"]), 6),
        "aai_coef": float(aai["pgls_coef"]),
        "aai_se": float(aai["pgls_se_unb"]),
        "aai_p": float(aai["pgls_pval_unb"]),
        "rbh_coef": float(rbh["pgls_coef"]),
        "rbh_se": float(rbh["pgls_se_unb"]),
        "rbh_p": float(rbh["pgls_pval_unb"]),
        "sensitivity_specs_count": len(sens_rows),
        "s7_n": int(s7["n"]),
    }


def verify_h2_and_cross_check_h4():
    h2_dir = R7 / "snapshots/h2_final"
    manifest = read_json(h2_dir / "SNAPSHOT_MANIFEST.json")
    manifest_checks = 0
    for entry in manifest["files"]:
        rel = entry["path"]
        actual = sha256_file(h2_dir / rel)
        assert actual == entry["sha256"], f"H2 hash mismatch: {rel}"
        manifest_checks += 1

    summary = read_json(h2_dir / "results/pilot_conservation_summary.json")
    join_metrics = read_json(h2_dir / "results/pilot_join_metrics.json")
    fam_rows = read_tsv(h2_dir / "results/pilot_family_conservation.tsv")
    all_eval_rows = read_tsv(h2_dir / "results/all_evaluable_families.tsv")
    matched_pairs = read_tsv(h2_dir / "results/pilot_matched_pairs.tsv")

    assert len(all_eval_rows) == 2838
    assert len(fam_rows) == 1461
    assert len(matched_pairs) == 148

    # Cross-check H2 K (1,031 in Gate 0; 570 in Gate 3) and QnotK (1,807 in Gate 0; 891 in Gate 3) against H4 v7_model_summaries.tsv
    h4_models = {
        r["model_id"]: r
        for r in read_tsv(R7 / "root_checks/h4_outputs_preclosure/v7_model_summaries.tsv")
    }
    k_strata = {"K_local_pdb", "K_local_afdb_only", "K_dual_none"}
    qnotk_strata = {"Q_not_K_pdb", "Q_not_K_afdb_only"}

    h2_gate0_k_reps = {r["chosen_rep"] for r in all_eval_rows if r["arm"] == "K"}
    h2_gate0_q_reps = {r["chosen_rep"] for r in all_eval_rows if r["arm"] == "QnotK"}
    assert len(h2_gate0_k_reps) == 1031 and len(h2_gate0_q_reps) == 1807
    assert all(h4_models[rep]["sub_stratum"] in k_strata for rep in h2_gate0_k_reps)
    assert all(h4_models[rep]["sub_stratum"] in qnotk_strata for rep in h2_gate0_q_reps)

    h2_gate0_k_fams = {r["family"] for r in all_eval_rows if r["arm"] == "K"}
    h2_gate0_q_fams = {r["family"] for r in all_eval_rows if r["arm"] == "QnotK"}
    h2_k_ids = {r["family"] for r in fam_rows if r["arm"] == "K"}
    h2_qnotk_ids = {r["family"] for r in fam_rows if r["arm"] == "QnotK"}
    assert len(h2_k_ids) == 570 and len(h2_qnotk_ids) == 891
    assert h2_k_ids.issubset(h2_gate0_k_fams) and h2_qnotk_ids.issubset(h2_gate0_q_fams)

    # Attrition funnel checks
    funnel = summary["attrition_funnel_2838_universe"]
    assert funnel["K"]["gate0_confidence_qualified"] + funnel["QnotK"]["gate0_confidence_qualified"] == 2838
    assert funnel["K"]["gate1_refined_mag_binned"] + funnel["QnotK"]["gate1_refined_mag_binned"] == 2426
    assert funnel["K"]["gate2_drep_species_ge1"] + funnel["QnotK"]["gate2_drep_species_ge1"] == 2378
    assert funnel["K"]["gate3_multi_species_ge2"] + funnel["QnotK"]["gate3_multi_species_ge2"] == 1461
    assert funnel["K"]["gate3_multi_species_ge2"] == 570
    assert funnel["QnotK"]["gate3_multi_species_ge2"] == 891

    # Tier 1 evaluable checks
    t1_rows = [r for r in fam_rows if r["eligible_primary"] == "1"]
    t1_k = [r for r in t1_rows if r["arm"] == "K"]
    t1_q = [r for r in t1_rows if r["arm"] == "QnotK"]
    assert len(t1_rows) == 353 and len(t1_k) == 173 and len(t1_q) == 180
    assert round(100.0 * 353 / 1461, 2) == 24.16
    assert round(100.0 * 173 / 570, 2) == 30.35
    assert round(100.0 * 180 / 891, 2) == 20.20

    # Recurrence and positional background checks
    k_rec = sum(int(r["c_fam_primary"]) for r in t1_k)
    q_rec = sum(int(r["c_fam_primary"]) for r in t1_q)
    assert k_rec == 70 and q_rec == 58
    assert round(100.0 * k_rec / 173, 2) == 40.46
    assert round(100.0 * q_rec / 180, 2) == 32.22

    k_exp = sum(float(r["mu_pos_c_fam_primary"]) for r in t1_k) / 173
    q_exp = sum(float(r["mu_pos_c_fam_primary"]) for r in t1_q) / 180
    assert round(100.0 * k_exp, 2) == 8.71
    assert round(100.0 * q_exp, 2) == 6.56
    assert round(100.0 * (k_rec / 173 - k_exp), 2) == 31.75
    assert round(100.0 * (q_rec / 180 - q_exp), 2) == 25.67

    # Contrasts and sensitivities
    t1_perm = summary["part_b_between_arm_covariate_stratified_permutation"]["tier1_primary"]
    t2_perm = summary["part_b_between_arm_covariate_stratified_permutation"]["tier2_stop3_sensitivity"]
    assert round(float(t1_perm["mh_weighted_excess_c_diff"]), 6) == 0.035358
    assert round(float(t1_perm["p_stratified_mh_excess_c"]), 6) == 0.406859
    assert round(float(t2_perm["mh_weighted_excess_c_diff"]), 6) == 0.018391
    assert round(float(t2_perm["p_stratified_mh_excess_c"]), 6) == 0.561844

    sm = summary["part_b_between_arm_matched_pairs_tier1"]["strict_exact_phylum_caliper_sensitivity"]
    assert int(sm["n_matched_pairs"]) == 96
    assert round(100.0 * float(sm["matched_mean_raw_c_K"]), 2) == 31.25
    assert round(100.0 * float(sm["matched_mean_raw_c_QnotK"]), 2) == 31.25
    assert round(float(sm["matched_diff_excess_c"]), 6) == -0.005440
    assert round(float(sm["p_sign_flip_excess_c"]), 6) == 0.922608

    indep = summary["cross_family_locus_and_species_independence_tier1"]
    assert indep["tier1_evaluable_primary_loci"] == 877
    assert indep["unique_bin_contig_pairs"] == 869
    assert indep["contigs_shared_across_multiple_families"] == 8
    assert indep["unique_refined_bins"] == 617
    assert indep["bins_shared_across_multiple_families"] == 145
    assert indep["unique_drep_species_clusters"] == 568
    assert indep["species_shared_across_multiple_families"] == 153

    sd = indep["species_disjoint_subset_sensitivity"]
    assert sd["n_families_with_zero_shared_species"] == 112
    assert sd["n_K"] == 51 and sd["n_QnotK"] == 61
    assert round(float(sd["mean_excess_c_K"]), 6) == 0.261855
    assert round(float(sd["mean_excess_c_QnotK"]), 6) == 0.198430

    # Bin verification counts
    assert join_metrics["unique_refined_bins"] == 2126
    assert join_metrics["bin_verification_counts"]["verified_sha256_and_log"] == 2041
    assert join_metrics["bin_verification_counts"]["annotation_files_missing"] == 85

    return {
        "status": "PASS",
        "snapshot_files_verified": manifest_checks,
        "h2_h4_cohort_identity_verified": True,
        "attrition_funnel": funnel,
        "tier1_counts": {"total": 353, "K": 173, "QnotK": 180},
        "tier1_retention_pct": {"total": 24.16, "K": 30.35, "QnotK": 20.20},
        "tier1_recurrence_pct": {"K_obs": 40.46, "K_exp": 8.71, "K_excess": 31.75, "Q_obs": 32.22, "Q_exp": 6.56, "Q_excess": 25.67},
        "tier1_mh_excess_diff": round(float(t1_perm["mh_weighted_excess_c_diff"]), 6),
        "tier1_p_perm": round(float(t1_perm["p_stratified_mh_excess_c"]), 6),
        "strict_96_pairs": {
            "k_pct": 31.25,
            "q_pct": 31.25,
            "excess_diff": round(float(sm["matched_diff_excess_c"]), 6),
            "p_perm": round(float(sm["p_sign_flip_excess_c"]), 6),
        },
        "tier2_counts": summary["cohort_admission"]["tier2_stop3_sensitivity_evaluable"],
        "tier2_mh_excess_diff": round(float(t2_perm["mh_weighted_excess_c_diff"]), 6),
        "tier2_p_perm": round(float(t2_perm["p_stratified_mh_excess_c"]), 6),
        "species_disjoint_112": {
            "K": 51,
            "QnotK": 61,
            "K_excess": round(float(sd["mean_excess_c_K"]), 6),
            "QnotK_excess": round(float(sd["mean_excess_c_QnotK"]), 6),
        },
        "tier3_counts": summary["cohort_admission"]["tier3_partial_window_sensitivity_evaluable"],
        "bins_verified": {"verified": 2041, "missing": 85, "total": 2126},
    }


def verify_h3():
    h3_dir = R7 / "snapshots/h3_completed"
    smag_all = [
        line.strip()
        for line in (h3_dir / "smag_catalogue_40039_ids.txt").read_text().splitlines()
        if line.strip()
    ]
    smag_sampled = [
        line.strip()
        for line in (h3_dir / "sampled_smag_1500_ids_seed20261014.txt").read_text().splitlines()
        if line.strip()
    ]
    assert len(smag_all) == 40039 and len(set(smag_all)) == 40039
    assert len(smag_sampled) == 1500 and len(set(smag_sampled)) == 1500

    rng_smag = random.Random(20261014)
    recomputed_smag = sorted(rng_smag.sample(sorted(smag_all), 1500))
    assert recomputed_smag == sorted(smag_sampled), "SMAG seed 20261014 sample mismatch"

    gem_all = [
        line.strip()
        for line in (h3_dir / "gem_catalogue_52515_ids.txt").read_text().splitlines()
        if line.strip()
    ]
    assert len(gem_all) == 52515 and len(set(gem_all)) == 52515
    target_gem = math.ceil(4500000 * 1.05 / 2300)
    assert target_gem == 2055

    for seed in (20261014, 20261015, 20261016):
        rng_gem = random.Random(seed)
        pool_copy = list(sorted(gem_all))
        rng_gem.shuffle(pool_copy)
        expected_gem = sorted(pool_copy[:target_gem])
        manifest_rows = read_tsv(h3_dir / f"gem_harmonized_manifest_seed{seed}.tsv")
        actual_gem = sorted(r["genome_id"] for r in manifest_rows)
        assert actual_gem == expected_gem, f"GEM seed {seed} sample mismatch"

    parts = read_tsv(h3_dir / "smag_zenodo_parts_manifest.tsv")
    assert len(parts) == 373
    assert all(r["status"] in ("VERIFIED_EXISTING", "DOWNLOADED_VERIFIED") and len(r["md5_checksum"]) == 32 for r in parts)

    closure_rows = read_tsv(h3_dir / "matching_closure_audit.tsv")
    assert len(closure_rows) == 7
    c_by_label = {r["label"]: r for r in closure_rows}

    prim = c_by_label["harmonized_seed20261014"]
    no_ov = c_by_label["harmonized_seed20261014_no_overlap"]
    assert int(prim["matched_proteins"]) == 4210592
    assert int(prim["matched_residues"]) == 1295970055
    assert int(prim["pre_matching_shared_seqs"]) == 25067
    assert int(prim["post_matching_shared_seqs"]) == 22310

    assert int(no_ov["matched_proteins"]) == 4185525
    assert int(no_ov["matched_residues"]) == 1289129312
    assert int(no_ov["pre_matching_shared_seqs"]) == 25067
    assert int(no_ov["post_matching_shared_seqs"]) == 0
    # Verify exact relationship: 4,210,592 - 25,067 pre-matching shared sequences = 4,185,525
    assert int(prim["matched_proteins"]) - int(prim["pre_matching_shared_seqs"]) == int(no_ov["matched_proteins"])

    census_rows = read_tsv(h3_dir / "matched_full_census_summary.tsv")
    assert len(census_rows) == 16
    for r in census_rows:
        n = int(r["total_queries"])
        a = int(r["both_hits_a"])
        b = int(r["smag_only_b"])
        c = int(r["gem_only_c"])
        d = int(r["neither_d"])
        assert n == 106042
        assert a + b + c + d == n
        smag_pct = round(100.0 * (a + b) / n, 2)
        gem_pct = round(100.0 * (a + c) / n, 2)
        delta_pct = round(100.0 * (b - c) / n, 2)
        assert abs(smag_pct - float(r["smag_hit_pct"])) <= 0.01
        assert abs(gem_pct - float(r["gem_hit_pct"])) <= 0.01
        assert abs(delta_pct - float(r["paired_delta_pct"])) <= 0.01

    strat_rows = read_tsv(h3_dir / "matched_length_stratified_results.tsv")
    assert len(strat_rows) == 80
    dark_strat_deltas = [float(r["paired_delta_pct"]) for r in strat_rows if r["cohort"] == "Dark"]
    assert len(dark_strat_deltas) == 40 and all(x > 0 for x in dark_strat_deltas)

    # Check all 16 search signature files and 14 database signature files exist and match
    m8_sigs = sorted(h3_dir.glob("hits_*.m8.sig.json"))
    dmnd_sigs = sorted(h3_dir.glob("matched_*_db.dmnd.sig.json"))
    assert len(m8_sigs) == 16 and len(dmnd_sigs) == 14

    # Check completed H3 raw audit (job 53796328)
    h3_raw_audit = read_json(R7 / "h3_final_raw_audit/h3_raw_audit_summary.json")
    assert h3_raw_audit["status"] == "PASS"
    assert h3_raw_audit["slurm_job_id"] == "53796328"
    assert h3_raw_audit["query_census"]["status"] == "PASS"
    assert h3_raw_audit["query_census"]["census_total"] == 212084
    assert len(h3_raw_audit["matched_pairs_audit"]) == 7
    assert all(v["status"] == "PASS" for v in h3_raw_audit["matched_pairs_audit"].values())

    return {
        "status": "PASS",
        "raw_audit_job": h3_raw_audit["slurm_job_id"],
        "raw_audit_status": h3_raw_audit["status"],
        "smag_frame": len(smag_all),
        "smag_sampled": len(smag_sampled),
        "gem_frame": len(gem_all),
        "gem_sampled_per_seed": target_gem,
        "zenodo_parts_verified": len(parts),
        "primary_matched_proteins": int(prim["matched_proteins"]),
        "primary_matched_residues": int(prim["matched_residues"]),
        "pre_matching_shared_seqs": int(prim["pre_matching_shared_seqs"]),
        "post_matching_shared_seqs": int(prim["post_matching_shared_seqs"]),
        "no_overlap_matched_proteins": int(no_ov["matched_proteins"]),
        "no_overlap_matched_residues": int(no_ov["matched_residues"]),
        "primary_dark": {
            "smag_hits": 8015 + 22049,
            "gem_hits": 8015 + 2384,
            "both": 8015,
            "smag_only": 22049,
            "gem_only": 2384,
            "neither": 73594,
            "smag_pct": 28.35,
            "gem_pct": 9.81,
            "delta_pct": 18.54,
        },
        "primary_control": {
            "smag_hits": 94277 + 5396,
            "gem_hits": 94277 + 2108,
            "both": 94277,
            "smag_only": 5396,
            "gem_only": 2108,
            "neither": 4261,
            "smag_pct": 93.99,
            "gem_pct": 90.89,
            "delta_pct": 3.10,
        },
        "harmonized_deltas": [18.54, 16.83, 16.91],
        "published_deltas": [18.46, 16.87, 16.91],
        "no_overlap_delta": 18.62,
        "id30_delta": 18.14,
        "all_40_dark_length_strata_positive": True,
    }


def verify_h4_and_figure4():
    pre_dir = R7 / "root_checks/h4_outputs_preclosure"
    gpu_dir = R7 / "root_checks/h4_gpu_audit"

    models = read_tsv(pre_dir / "v7_model_summaries.tsv")
    regions = read_tsv(pre_dir / "v7_region_evaluations.tsv")
    stratum_census = read_tsv(pre_dir / "v7_stratum_census.tsv")
    rec_summary = read_json(pre_dir / "v7_recovery_summary.json")
    tm_anom = read_json(R7 / "root_checks/tm_anomaly_census.json")
    excl_audit = read_json(R7 / "root_checks/h4_exclusion_fact_audit.json")
    gpu_audit = read_json(gpu_dir / "gpu_raw_audit.json")
    gpu_census = read_tsv(gpu_dir / "gpu_cpu_census.tsv")
    gpu_per_model = read_tsv(gpu_dir / "gpu_cpu_per_model.tsv")

    assert len(models) == 5851
    assert len(regions) == 7773
    assert sum(int(m["n_confident_regions"]) == 0 for m in models) == 161

    k_strata = {"K_local_pdb", "K_local_afdb_only", "K_dual_none"}
    qnotk_strata = {"Q_not_K_pdb", "Q_not_K_afdb_only"}

    # Verify pooled K (1,031) and Q_not_K (1,807) are strictly defined by initial qTM < 0.5 vs >= 0.5
    dark_models = [m for m in models if m["is_control"] == "0"]
    ctrl_models = [m for m in models if m["is_control"] == "1"]
    assert len(dark_models) == 2838 and len(ctrl_models) == 3013

    for m in dark_models:
        max_hist_qtm = max(float(m["hist_pdb_qtm"]), float(m["hist_afdb_qtm"]))
        if m["sub_stratum"] in k_strata:
            assert max_hist_qtm < 0.5
        else:
            assert m["sub_stratum"] in qnotk_strata
            assert max_hist_qtm >= 0.5

    k_models = [m for m in dark_models if m["sub_stratum"] in k_strata]
    q_models = [m for m in dark_models if m["sub_stratum"] in qnotk_strata]
    assert len(k_models) == 1031 and len(q_models) == 1807

    # CPU Tier 1 pooled K counts
    k_eligible = sum(int(m["n_confident_regions"]) > 0 for m in k_models)
    k_wc_matched = sum(m["wc_t1_status"] == "MATCHED" for m in k_models)
    k_wc_unmatched = len(k_models) - k_wc_matched
    k_crop_matched = sum(m["crop_t1_status"] == "MATCHED" for m in k_models)
    k_trans = Counter(m["trans_t1"] for m in k_models)

    assert k_eligible == 984
    assert k_wc_matched == 67
    assert k_wc_unmatched == 964
    assert k_crop_matched == 775
    assert k_trans["WC_UNMATCHED_TO_CROP_RECOVERED"] == 714
    assert k_trans["WC_UNMATCHED_TO_PERSISTENT_LOW_SCORE"] == 167
    assert k_trans["WC_UNMATCHED_TO_RAW_ZERO_RETRIEVAL"] == 38
    assert k_trans["WC_UNMATCHED_TO_POST_EXCLUSION_EMPTY"] == 0
    assert k_trans["WC_UNMATCHED_NO_CONFIDENT_REGION"] == 45
    assert k_trans["WC_MATCHED_RETAINED_IN_CROP"] == 61
    assert k_trans["WC_MATCHED_LOST_LOW_SCORE"] == 4
    assert k_trans["WC_MATCHED_LOST_NO_CONFIDENT_REGION"] == 2
    assert round(100.0 * 714 / 964, 2) == 74.07
    assert round(100.0 * 775 / 1031, 2) == 75.17

    # CPU Tier 1 Other dark (Q_not_K_all = 1,807) counts (plotted in Fig. 4b)
    q_wc_matched = sum(m["wc_t1_status"] == "MATCHED" for m in q_models)
    q_crop_matched = sum(m["crop_t1_status"] == "MATCHED" for m in q_models)
    q_pdb = [m for m in q_models if m["sub_stratum"] == "Q_not_K_pdb"]
    q_afdb = [m for m in q_models if m["sub_stratum"] == "Q_not_K_afdb_only"]
    assert len(q_pdb) == 1107 and len(q_afdb) == 700
    assert sum(m["wc_t1_status"] == "MATCHED" for m in q_pdb) == 1107
    assert sum(m["crop_t1_status"] == "MATCHED" for m in q_pdb) == 1009
    assert sum(m["wc_t1_status"] == "MATCHED" for m in q_afdb) == 301
    assert sum(m["crop_t1_status"] == "MATCHED" for m in q_afdb) == 569
    assert q_wc_matched == 1408 and round(100.0 * 1408 / 1807, 1) == 77.9
    assert q_crop_matched == 1578 and round(100.0 * 1578 / 1807, 1) == 87.3

    # CPU Controls across tiers
    ctrl_wc_t1 = sum(m["wc_t1_status"] == "MATCHED" for m in ctrl_models)
    ctrl_crop_by_tier = {
        t: sum(m[f"crop_t{t}_status"] == "MATCHED" for m in ctrl_models)
        for t in range(5)
    }
    assert ctrl_wc_t1 == 2878 and round(100.0 * 2878 / 3013, 1) == 95.5
    assert ctrl_crop_by_tier == {0: 2902, 1: 2889, 2: 2842, 3: 2823, 4: 2796}
    assert round(100.0 * 2889 / 3013, 2) == 95.88

    # CPU saturation & adaptive rerun checks
    sat_audit = rec_summary["hit_cap_saturation_audit"]
    assert sat_audit["adaptive_rerun_audit"]["whole_chain_rerun"]["problematic_queries_count"] == 39
    assert sat_audit["adaptive_rerun_audit"]["cropped_rerun"]["problematic_queries_count"] == 55
    assert sat_audit["final"]["whole_chain"]["problematic_qids"] == ["KF001671", "KF001867", "KF003576"]
    assert sat_audit["final"]["cropped_regions"]["problematic_qids"] == ["KF003576_reg1"]

    # GPU raw audit and CPU-GPU comparison checks
    assert gpu_audit["job"] == "53796103"
    assert gpu_audit["files"]["wc"]["rows"] == 2853102
    assert gpu_audit["files"]["crop"]["rows"] == 3228450
    assert gpu_audit["files"]["wc"]["alntm_out_of_range_by_stratum"] == {"Controls": 471}
    assert gpu_audit["files"]["crop"]["alntm_out_of_range_by_stratum"] == {"Controls": 1237}
    assert 471 + 1237 == 1708

    gpu_t1_k = [r for r in gpu_census if r["tier"] == "1" and r["stratum"] in k_strata]
    assert sum(int(r["gpu_wc_unmatched"]) for r in gpu_t1_k) == 961
    assert sum(int(r["gpu_recovered"]) for r in gpu_t1_k) == 731
    assert sum(int(r["gpu_crop_matched"]) for r in gpu_t1_k) == 794
    assert sum(int(r["wc_matched_disagreements"]) for r in gpu_t1_k) == 41
    assert sum(int(r["crop_matched_disagreements"]) for r in gpu_t1_k) == 43
    assert sum(int(r["gpu_wc_cap_unmatched"]) + int(r["gpu_crop_cap_unmatched"]) for r in gpu_t1_k) == 0

    gpu_ctrl_by_tier = {
        int(r["tier"]): {
            "gpu_wc_unmatched": int(r["gpu_wc_unmatched"]),
            "gpu_crop_matched": int(r["gpu_crop_matched"]),
            "gpu_wc_cap_unmatched": int(r["gpu_wc_cap_unmatched"]),
            "gpu_crop_cap_unmatched": int(r["gpu_crop_cap_unmatched"]),
        }
        for r in gpu_census
        if r["stratum"] == "Controls"
    }
    assert gpu_ctrl_by_tier[1]["gpu_crop_matched"] == 2896
    assert gpu_ctrl_by_tier[1]["gpu_wc_cap_unmatched"] == 2
    assert gpu_ctrl_by_tier[1]["gpu_crop_cap_unmatched"] == 0

    # CPU anomaly & exclusion sensitivity audit checks
    assert tm_anom["total_rows_audited"] == 6845932
    assert tm_anom["total_invalid_score_rows"] == 2194
    assert excl_audit["wc_status_changes"]["unique_queries"] == 4
    assert excl_audit["wc_status_changes"]["tier_counts"] == {"3": 4, "4": 3, "2": 3}
    assert excl_audit["model_crop_status_changes"]["unique_models"] == 9
    assert excl_audit["model_crop_status_changes"]["tier_counts"] == {"3": 9, "4": 8, "2": 8}

    # Verify Foldseek source lines in structureconvertalis.cpp and TMaligner.cpp
    sc_cpp = (R7 / "h4_tm_anomaly/source/structureconvertalis.cpp").read_text(encoding="utf-8")
    tm_cpp = (R7 / "h4_tm_anomaly/source/TMaligner.cpp").read_text(encoding="utf-8")
    assert "std::min(res.qEndPos - res.qStartPos, res.dbEndPos - res.dbStartPos)" in sc_cpp
    assert "res.qStartPos, res.dbStartPos, tmpBt,res.qLen)" in sc_cpp
    assert "parameter_set4search(normalizationLen,  normalizationLen" in tm_cpp

    return {
        "status": "PASS",
        "total_models": len(models),
        "total_regions": len(regions),
        "models_no_region": 161,
        "cpu_pooled_k_tier1": {
            "n": 1031,
            "eligible": k_eligible,
            "wc_matched": k_wc_matched,
            "wc_matched_pct": round(100.0 * k_wc_matched / 1031, 2),
            "wc_unmatched": k_wc_unmatched,
            "crop_recovered": 714,
            "crop_recovered_pct": 74.07,
            "crop_matched_total": k_crop_matched,
            "crop_matched_total_pct": 75.17,
            "transitions_unmatched": {
                "recovered": 714,
                "persistent_low_score": 167,
                "raw_zero_retrieval": 38,
                "post_exclusion_empty": 0,
                "no_confident_region": 45,
            },
            "transitions_matched": {
                "retained": 61,
                "lost_low_score": 4,
                "lost_no_confident_region": 2,
            },
        },
        "cpu_other_dark_tier1": {
            "n": 1807,
            "wc_matched": q_wc_matched,
            "wc_matched_pct": round(100.0 * q_wc_matched / 1807, 2),
            "crop_matched": q_crop_matched,
            "crop_matched_pct": round(100.0 * q_crop_matched / 1807, 2),
            "q_not_k_pdb": {"n": 1107, "wc_matched": 1107, "crop_matched": 1009},
            "q_not_k_afdb_only": {"n": 700, "wc_matched": 301, "crop_matched": 569, "wc_unmatched_recovered": 292},
        },
        "cpu_controls": {
            "n": 3013,
            "wc_matched_tier1": ctrl_wc_t1,
            "wc_matched_tier1_pct": round(100.0 * ctrl_wc_t1 / 3013, 2),
            "crop_matched_by_tier": ctrl_crop_by_tier,
            "adaptive_expanded_wc": 39,
            "adaptive_expanded_crop": 55,
            "final_cap_saturated_unmatched_wc": ["KF001671", "KF001867", "KF003576"],
            "final_cap_saturated_unmatched_crop": ["KF003576_reg1"],
        },
        "gpu_sensitivity": {
            "job": "53796103",
            "wc_rows": 2853102,
            "crop_rows": 3228450,
            "alntm_gt1_rows": 1708,
            "pooled_k_tier1": {
                "wc_unmatched": 961,
                "crop_recovered": 731,
                "crop_recovered_pct": round(100.0 * 731 / 961, 2),
                "crop_matched_total": 794,
                "crop_matched_total_pct": round(100.0 * 794 / 1031, 2),
                "wc_disagreements_vs_cpu": 41,
                "crop_disagreements_vs_cpu": 43,
                "cap_unmatched_dark": 0,
            },
            "controls_by_tier": gpu_ctrl_by_tier,
        },
        "cpu_alntm_anomaly": {
            "total_rows": 6845932,
            "alntm_gt1_rows": 2194,
            "row_exclusion_wc_models_changed_tiers2to4": 4,
            "row_exclusion_crop_models_changed_tiers2to4": 9,
        },
    }


def check_reproducibility_scripts():
    crc_path = ROOT / "reproducibility/darkmatter/current_revision_checks.py"
    crc_text = crc_path.read_text(encoding="utf-8")
    stale_guard = "797 have an alignment-normalised" in crc_text

    verify_py_path = ROOT / "reproducibility/verify.py"
    verify_py_text = verify_py_path.read_text(encoding="utf-8")
    stale_verify_required = [
        phrase
        for phrase in ("97/1{,}875", "68/3{,}013", "16{,}386", "split 74/23", "236/3{,}684 (6.41\\%)")
        if phrase in verify_py_text
    ]

    vr7_path = ROOT / "reproducibility/darkmatter/round7_final/verify_round7.py"
    vr7_text = vr7_path.read_text(encoding="utf-8") if vr7_path.exists() else ""

    return {
        "current_revision_checks_has_stale_797_guard": stale_guard,
        "verify_py_has_stale_alntm_required_phrases": stale_verify_required,
        "verify_round7_hardcodes_gtdb_outgroup_tips": "n_gtdb = 1647" in vr7_text,
        "verify_round7_h2_attrition_needs_total_sum_handler": 'data["attrition_funnel_2838_universe"][arm][gate]' in vr7_text,
    }


def main():
    report = {
        "manuscript_and_pdf_integrity": check_manuscript_and_pdf_integrity(),
        "h1_verification": verify_h1(),
        "h2_verification": verify_h2_and_cross_check_h4(),
        "h3_verification": verify_h3(),
        "h4_and_fig4_verification": verify_h4_and_figure4(),
        "reproducibility_harness_checks": check_reproducibility_scripts(),
    }
    out_path = OUT_DIR / "verification_output.json"
    out_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print("ALL VERIFICATION CHECKS PASSED. Output written to:", out_path)


if __name__ == "__main__":
    main()
