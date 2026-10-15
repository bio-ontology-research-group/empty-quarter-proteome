"""verify_round7.py - Bounded reproducibility verifier helper for H1-H4.

This module is imported by reproducibility/verify.py. It evaluates active
quantitative claims for Hypotheses H1-H4 from pinned compact artifacts under
reproducibility/darkmatter/round7_final/{h1,h2,h3,h4}/.

Counts and rates are derived from original compact tables rather than static
assert comparisons where possible.
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
H1_DIR = os.path.join(HERE, "h1")
H2_DIR = os.path.join(HERE, "h2")
H3_DIR = os.path.join(HERE, "h3")
H4_DIR = os.path.join(HERE, "h4")


def _read_tsv(path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _read_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


# ===========================================================================
# H1: Within-Thermoplasmatota PGLS & Reference Proximity
# ===========================================================================
def h1_pgls_model(variable, field):
    """Extract coefficient or SE from v3_pgls_model_summary.tsv."""
    rows = _read_tsv(os.path.join(H1_DIR, "v3_pgls_model_summary.tsv"))
    for r in rows:
        if r["variable"] == variable:
            return float(r[field])
    raise KeyError((variable, field))


def h1_closure_val(key):
    """Extract scalar from v3_closure_verification.json."""
    data = _read_json(os.path.join(H1_DIR, "v3_closure_verification.json"))
    if key in ("reconciliation_pass", "status"):
        return 1.0 if data["status"] == "PASSED_ALL_ASSERTIONS" else 0.0
    return float(data[key])


def h1_closure_lambda_ci(index):
    """Extract lower (0) or upper (1) bound of 95% profile likelihood CI."""
    data = _read_json(os.path.join(H1_DIR, "v3_closure_verification.json"))
    return float(data["primary_lambda_ci95"][index])


def h1_comparator_meta(field, stratum="All_Thermoplasmatota_237"):
    """Extract counts from v3_comparator_metadata_summary.tsv."""
    rows = _read_tsv(os.path.join(H1_DIR, "v3_comparator_metadata_summary.tsv"))
    for r in rows:
        if r["stratum"] == stratum:
            return float(r[field])
    raise KeyError((stratum, field))


def h1_tree_taxa_count(subset="total"):
    """Count tips in tree dynamically from pinned thermo_tree_taxa.tsv."""
    taxa_path = os.path.join(H1_DIR, "thermo_tree_taxa.tsv")
    rows = _read_tsv(taxa_path)
    if subset == "total":
        return float(len(rows))
    elif subset in ("desert", "gtdb_thermoplasmatota", "outgroup"):
        return float(sum(1 for r in rows if r.get("set") == subset))
    raise KeyError(subset)


# ===========================================================================
# H2: Flanking-Orthologue Recurrence & Genomic Context
# ===========================================================================
def h2_attrition(arm, gate):
    """Extract family counts across attrition gates."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    return float(data["attrition_funnel_2838_universe"][arm][gate])


def h2_cohort_admission(tier, key):
    """Extract evaluable family count by tier and arm."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    return float(data["cohort_admission"][tier][key])


def h2_recurrence_stat(arm, stat):
    """Extract observed, expected, or excess recurrence."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    sec = data["tier1_primary_gate_unadjusted"]
    if stat == "observed_pct":
        k = f"{arm}_conserved"
        return round(100.0 * float(sec[k]["rate"]), 2)
    elif stat == "expected_pct":
        # expected = observed - excess
        obs = float(sec[f"{arm}_conserved"]["rate"])
        excess = float(sec["mean_excess_c_fam"][arm])
        return round(100.0 * (obs - excess), 2)
    elif stat == "excess_pp":
        return round(100.0 * float(sec["mean_excess_c_fam"][arm]), 2)
    raise KeyError((arm, stat))


def h2_contrast_stat(tier, stat):
    """Extract Mantel-Haenszel contrast or permutation p-value."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    sec = data["part_b_between_arm_covariate_stratified_permutation"][tier]
    if stat == "mh_excess_diff":
        return float(sec["mh_weighted_excess_c_diff"])
    elif stat == "p_perm":
        return float(sec["p_stratified_mh_excess_c"])
    raise KeyError((tier, stat))


def h2_strict_matching_stat(stat):
    """Extract strict matching statistics."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    sec = data["part_b_between_arm_matched_pairs_tier1"]["strict_exact_phylum_caliper_sensitivity"]
    if stat == "pairs":
        return float(sec["n_matched_pairs"])
    elif stat == "recurrence_pct":
        return round(100.0 * float(sec["matched_mean_raw_c_K"]), 2)
    elif stat == "excess_diff":
        return float(sec["matched_diff_excess_c"])
    elif stat == "p_perm":
        return float(sec["p_sign_flip_excess_c"])
    raise KeyError(stat)


def h2_species_disjoint_stat(stat):
    """Extract species-disjoint sensitivity statistics."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    sec = data["cross_family_locus_and_species_independence_tier1"]["species_disjoint_subset_sensitivity"]
    if stat in ("families", "n_families_with_zero_shared_species"):
        return float(sec["n_families_with_zero_shared_species"])
    elif stat in ("k_families", "n_K"):
        return float(sec["n_K"])
    elif stat in ("qnotk_families", "n_QnotK"):
        return float(sec["n_QnotK"])
    elif stat in ("k_excess", "mean_excess_c_K"):
        return float(sec["mean_excess_c_K"])
    elif stat in ("qnotk_excess", "mean_excess_c_QnotK"):
        return float(sec["mean_excess_c_QnotK"])
    elif stat in ("diff_excess", "diff_excess_c"):
        return float(sec["diff_excess_c"])
    raise KeyError(stat)


def h2_locus_independence_stat(stat):
    """Extract locus and species independence metrics."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    sec = data["cross_family_locus_and_species_independence_tier1"]
    alias_map = {
        "unique_bin_contigs": "unique_bin_contig_pairs",
        "unique_bins": "unique_refined_bins",
        "unique_drep_species": "unique_drep_species_clusters",
        "contigs_with_ge2_families": "contigs_shared_across_multiple_families",
        "bins_with_ge2_families": "bins_shared_across_multiple_families",
        "species_with_ge2_families": "species_shared_across_multiple_families",
    }
    key = alias_map.get(stat, stat)
    return float(sec[key])


def h2_bin_verification(stat):
    """Extract contig/bin join metrics."""
    data = _read_json(os.path.join(H2_DIR, "pilot_join_metrics.json"))
    if stat == "verified":
        return float(data["bin_verification_counts"]["verified_sha256_and_log"])
    elif stat == "missing":
        return float(data["bin_verification_counts"]["annotation_files_missing"])
    elif stat == "total":
        return float(data["unique_refined_bins"])
    raise KeyError(stat)


def h2_window_retention_pct(arm):
    """Derive complete-window retention percentage."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    evaluable = float(data["cohort_admission"]["tier1_primary_evaluable"][arm])
    total = float(data["cohort_admission"]["selected_by_arm"][arm])
    return round(100.0 * evaluable / total, 2)


def h2_multispecies_gate_retention_pct():
    """Derive percentage of multispecies families passing Tier 1."""
    data = _read_json(os.path.join(H2_DIR, "pilot_conservation_summary.json"))
    evaluable = float(data["cohort_admission"]["tier1_primary_evaluable"]["total"])
    total = float(data["cohort_admission"]["total_families_selected"])
    return round(100.0 * evaluable / total, 2)


def h2_family_table_stat(stat):
    """Derive evaluable family counts and loci sums directly from pilot_family_conservation.tsv."""
    rows = _read_tsv(os.path.join(H2_DIR, "pilot_family_conservation.tsv"))
    if stat == "tier1_evaluable_total":
        return float(len([r for r in rows if r["eligible_primary"] == "1"]))
    elif stat == "tier1_evaluable_k":
        return float(len([r for r in rows if r["eligible_primary"] == "1" and r["arm"] == "K"]))
    elif stat == "tier1_evaluable_qnotk":
        return float(len([r for r in rows if r["eligible_primary"] == "1" and r["arm"] == "QnotK"]))
    elif stat == "tier1_selected_loci":
        return float(sum(int(r["n_species_complete_window"]) for r in rows if r["eligible_primary"] == "1"))
    elif stat == "tier2_evaluable_total":
        return float(len([r for r in rows if r["eligible_stop3_sens"] == "1"]))
    elif stat == "tier2_evaluable_k":
        return float(len([r for r in rows if r["eligible_stop3_sens"] == "1" and r["arm"] == "K"]))
    elif stat == "tier2_evaluable_qnotk":
        return float(len([r for r in rows if r["eligible_stop3_sens"] == "1" and r["arm"] == "QnotK"]))
    elif stat == "tier3_evaluable_total":
        return float(len([r for r in rows if r["eligible_sensitivity"] == "1"]))
    raise KeyError(stat)


# ===========================================================================
# H3: Matched-Reference Comparison (SMAG vs GEM)
# ===========================================================================
def h3_census_val(run_label, cohort, col):
    """Extract or compute value from matched_full_census_summary.tsv."""
    rows = _read_tsv(os.path.join(H3_DIR, "matched_full_census_summary.tsv"))
    for r in rows:
        if r["run_label"] == run_label and r["cohort"] == cohort:
            if col == "smag_hits":
                return float(int(r["both_hits_a"]) + int(r["smag_only_b"]))
            elif col == "gem_hits":
                return float(int(r["both_hits_a"]) + int(r["gem_only_c"]))
            return float(r[col])
    raise KeyError((run_label, cohort, col))


def h3_closure_val(label, col):
    """Extract value from matching_closure_audit.tsv."""
    rows = _read_tsv(os.path.join(H3_DIR, "matching_closure_audit.tsv"))
    for r in rows:
        if r["label"] == label:
            return float(r[col])
    raise KeyError((label, col))


def h3_zenodo_parts_count():
    """Count archive parts in smag_zenodo_parts_manifest.tsv."""
    rows = _read_tsv(os.path.join(H3_DIR, "smag_zenodo_parts_manifest.tsv"))
    return float(len(rows))


# ===========================================================================
# H4: Paired PDB Structural Searches & Retrieval Sensitivity
# ===========================================================================
def _get_stratum_census():
    return _read_tsv(os.path.join(H4_DIR, "v7_stratum_census.tsv"))


def h4_census_pooled_k(stat):
    """Derive pooled K (K_local_pdb + K_local_afdb_only + K_dual_none) endpoints."""
    rows = _get_stratum_census()
    k_strata = {"K_local_pdb", "K_local_afdb_only", "K_dual_none"}
    k_rows = [r for r in rows if r["tier"] == "tier_1" and r["stratum"] in k_strata]
    if stat == "model_count":
        return float(sum(int(r["model_count"]) for r in k_rows))
    elif stat == "models_eligible":
        return float(sum(int(r["models_eligible"]) for r in k_rows))
    elif stat == "wc_matched":
        return float(sum(int(r["wc_matched"]) for r in k_rows))
    elif stat == "wc_unmatched":
        tot = sum(int(r["model_count"]) for r in k_rows)
        matched = sum(int(r["wc_matched"]) for r in k_rows)
        return float(tot - matched)
    elif stat == "crop_recovered":
        return float(sum(int(r["wc_unmatched_recovered_count"]) for r in k_rows))
    elif stat == "crop_recovery_pct":
        tot = sum(int(r["model_count"]) for r in k_rows)
        matched = sum(int(r["wc_matched"]) for r in k_rows)
        rec = sum(int(r["wc_unmatched_recovered_count"]) for r in k_rows)
        return round(100.0 * rec / (tot - matched), 2)
    elif stat == "crop_matched_total":
        return float(sum(int(r["crop_matched"]) for r in k_rows))
    raise KeyError(stat)


def h4_census_control_crops(tier):
    """Extract control crop matches for tier (tier_1 to tier_4)."""
    rows = _get_stratum_census()
    for r in rows:
        if r["stratum"] == "Controls" and r["tier"] == tier:
            return float(r["crop_matched"])
    raise KeyError(tier)


def h4_census_control_crops_pct(tier):
    """Extract control crop match percentage for tier."""
    rows = _get_stratum_census()
    for r in rows:
        if r["stratum"] == "Controls" and r["tier"] == tier:
            return round(100.0 * int(r["crop_matched"]) / int(r["model_count"]), 2)
    raise KeyError(tier)


def h4_total_models():
    """Derive total evaluated models across all strata."""
    rows = _get_stratum_census()
    tier1 = [r for r in rows if r["tier"] == "tier_1" and not r["stratum"].endswith("_all")]
    return float(sum(int(r["model_count"]) for r in tier1))


def h4_total_crops():
    """Derive total evaluated crops across all strata."""
    rows = _get_stratum_census()
    tier1 = [r for r in rows if r["tier"] == "tier_1" and not r["stratum"].endswith("_all")]
    return float(sum(int(r["total_regions"]) for r in tier1))


def h4_models_no_region():
    """Derive total models lacking an eligible region."""
    rows = _get_stratum_census()
    tier1 = [r for r in rows if r["tier"] == "tier_1" and not r["stratum"].endswith("_all")]
    return float(sum(int(r["models_no_region"]) for r in tier1))


def h4_gpu_census_pooled_k(stat):
    """Derive GPU pooled K endpoints from gpu_cpu_census.tsv."""
    rows = _read_tsv(os.path.join(H4_DIR, "gpu_cpu_census.tsv"))
    k_strata = {"K_local_pdb", "K_local_afdb_only", "K_dual_none"}
    t1_rows = [r for r in rows if r["tier"] == "1" and r["stratum"] in k_strata]
    if stat == "gpu_wc_unmatched":
        return float(sum(int(r["gpu_wc_unmatched"]) for r in t1_rows))
    elif stat == "gpu_recovered":
        return float(sum(int(r["gpu_recovered"]) for r in t1_rows))
    elif stat == "gpu_crop_matched":
        return float(sum(int(r["gpu_crop_matched"]) for r in t1_rows))
    elif stat == "wc_disagreements":
        return float(sum(int(r["wc_matched_disagreements"]) for r in t1_rows))
    elif stat == "crop_disagreements":
        return float(sum(int(r["crop_matched_disagreements"]) for r in t1_rows))
    elif stat == "gpu_cap_unmatched":
        return float(sum(int(r["gpu_wc_cap_unmatched"]) + int(r["gpu_crop_cap_unmatched"]) for r in t1_rows))
    raise KeyError(stat)


def h4_gpu_control_stat(stat):
    """Extract GPU control tier 1 statistics."""
    rows = _read_tsv(os.path.join(H4_DIR, "gpu_cpu_census.tsv"))
    for r in rows:
        if r["tier"] == "1" and r["stratum"] == "Controls":
            if stat == "crop_matched":
                return float(r["gpu_crop_matched"])
            elif stat == "wc_cap_unmatched":
                return float(r["gpu_wc_cap_unmatched"])
    raise KeyError(stat)


def h4_cap_audit(stat):
    """Extract CPU adaptive rerun hit cap expansion and saturation counts."""
    data = _read_json(os.path.join(H4_DIR, "v7_recovery_summary.json"))
    audit = data["hit_cap_saturation_audit"]
    if stat == "wc_expanded":
        return float(audit["adaptive_rerun_audit"]["whole_chain_rerun"]["problematic_queries_count"])
    elif stat == "crop_expanded":
        return float(audit["adaptive_rerun_audit"]["cropped_rerun"]["problematic_queries_count"])
    elif stat == "wc_saturated_unmatched":
        return float(len(audit["final"]["whole_chain"]["problematic_qids"]))
    elif stat == "crop_saturated_unmatched":
        return float(len(audit["final"]["cropped_regions"]["problematic_qids"]))
    raise KeyError(stat)


def h4_gpu_raw_stat(stat):
    """Extract alignment row counts and out-of-range rows from GPU raw audit."""
    data = _read_json(os.path.join(H4_DIR, "gpu_raw_audit.json"))
    if stat == "wc_rows":
        return float(data["files"]["wc"]["rows"])
    elif stat == "crop_rows":
        return float(data["files"]["crop"]["rows"])
    elif stat == "alntm_gt1":
        wc_anom = data["files"]["wc"]["alntm_out_of_range_by_stratum"].get("Controls", 0)
        crop_anom = data["files"]["crop"]["alntm_out_of_range_by_stratum"].get("Controls", 0)
        return float(wc_anom + crop_anom)
    raise KeyError(stat)


def h4_cpu_anomaly_stat(stat):
    """Extract CPU raw rows and alnTM > 1 anomalies."""
    data = _read_json(os.path.join(H4_DIR, "tm_anomaly_census.json"))
    if stat == "rows_audited":
        return float(data["total_rows_audited"])
    elif stat == "alntm_gt1":
        return float(data["total_invalid_score_rows"])
    raise KeyError(stat)


def h4_model_crop_status_changes():
    """Count unique control models whose crop status changed under alnTM>1 exclusion."""
    data = _read_json(os.path.join(H4_DIR, "h4_exclusion_fact_audit.json"))
    return float(data["model_crop_status_changes"]["unique_models"])



def h4_stratum_value(stratum, field, tier='tier_1'):
    return float(next(r[field] for r in _get_stratum_census() if r['stratum']==stratum and r['tier']==tier))


def h4_k_transition(field):
    return sum(h4_stratum_value(s,field) for s in ['K_local_pdb','K_local_afdb_only','K_dual_none'])


def h1_same_class_n():
    return float(next(r['n'] for r in _read_tsv(os.path.join(H1_DIR,'v3_pgls_sensitivity.tsv')) if r['model_id'].startswith('S6_')))


def h4_gpu_control_by_tier(tier):
    return float(next(r['gpu_crop_matched'] for r in _read_tsv(os.path.join(H4_DIR,'gpu_cpu_census.tsv')) if r['stratum']=='Controls' and int(r['tier'])==tier))
