#!/usr/bin/env python3
"""Independent verification of exported panel TSVs against raw source tables."""
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
PACKAGE_REPRO = (
    HERE.parent
    if (HERE.parent / "darkmatter").is_dir()
    else Path(os.environ.get("EQ_PROTEOME_REPRO_ROOT", "/home/leechuck/Public/software/empty-quarter-proteome/reproducibility"))
)
PANEL_DIR = Path(os.environ.get("EQ_PANEL_DIR", str(HERE / "panel_data")))
COMP_RES = Path(
    os.environ.get(
        "EQ_COMP_RES",
        str(
            (HERE / "inputs/compute_results")
            if (HERE / "inputs/compute_results").is_dir()
            else (PACKAGE_REPRO / "reframed_20261108/inputs/compute_results")
        ),
    )
)
REPRO_DM = Path(os.environ.get("EQ_REPRO_DM", str(PACKAGE_REPRO / "darkmatter")))
ARCH_RES = Path(os.environ.get("EQ_ARCH_RES", str(PACKAGE_REPRO / "integration_20261108/results")))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def test_panel_sources_sha256_manifest():
    manifest = json.loads((PANEL_DIR / "panel_sources_sha256.json").read_text())
    assert len(manifest["sources"]) == 18
    assert len(manifest["panels"]) == 14
    for key, rec in manifest["sources"].items():
        if rec.get("repro_rel_path") and (PACKAGE_REPRO / rec["repro_rel_path"]).exists():
            p = PACKAGE_REPRO / rec["repro_rel_path"]
        else:
            p = Path(rec["path"])
        assert p.exists(), f"Missing source {key}: {p}"
        assert sha256_file(p) == rec["sha256"], f"Hash mismatch for source {key}"
    for key, rec in manifest["panels"].items():
        p = PANEL_DIR / Path(rec["path"]).name
        assert p.exists(), f"Missing panel {key}: {p}"
        assert sha256_file(p) == rec["sha256"], f"Hash mismatch for panel {key}"
        df = pd.read_csv(p, sep="\t")
        assert len(df) == rec["rows"], f"Row count mismatch for panel {key}"


def test_fig1_panels_against_sources():
    p1a = pd.read_csv(PANEL_DIR / "panel_fig1a.tsv", sep="\t")
    src_map = pd.read_csv(ARCH_RES / "mapping_group_review.tsv", sep="\t")
    exp_1a = src_map[
        (src_map["mode"] == "all_primary")
        & (src_map["metric"] == "relative_abundance")
        & (src_map["group"] == "EQ_Thermoplasmatota")
        & (src_map["threshold"].isin([0.1, 0.2, 0.5]))
    ].sort_values("threshold")
    assert list(p1a["breadth_pct"]) == [10, 20, 50]
    assert list(p1a["sites"]) == [36, 36, 36]
    assert list(p1a["visits"]) == [68, 68, 68]
    assert list(p1a["positive_sites"]) == [35, 35, 28]
    assert np.allclose(p1a["mean"], exp_1a["mean"])
    assert np.allclose(p1a["median"], exp_1a["median"])

    p1b = pd.read_csv(PANEL_DIR / "panel_fig1b.tsv", sep="\t")
    src_site = pd.read_csv(REPRO_DM / "natcomms_20261020/results/ecology_site_contrasts.tsv", sep="\t")
    exp_1b = src_site[(src_site["metric"] == "thermo_read_percent") & (src_site["contrast"] == "Deep-Surface")]
    assert len(p1b) == 37
    assert int((p1b["difference"] > 0).sum()) == 36
    assert np.isclose(p1b["Surface"].median(), exp_1b["Surface"].median())
    assert np.isclose(p1b["Deep"].median(), exp_1b["Deep"].median())

    p1c = pd.read_csv(PANEL_DIR / "panel_fig1c.tsv", sep="\t").set_index("feature")
    src_met = pd.read_csv(REPRO_DM / "round4_corrections/inputs/v2_thermo_metabolism_summary.tsv.gz", sep="\t")
    src_met = src_met[src_met["group"] == "thermo_237"].set_index("feature")
    for feat, r in p1c.iterrows():
        assert r["n_genomes"] == 237
        assert r["n_present"] == int(src_met.loc[feat, "n_present"])
        assert np.isclose(r["pct_present"], 100.0 * src_met.loc[feat, "share"])


def test_fig2_panels_against_sources():
    p2a = pd.read_csv(PANEL_DIR / "panel_fig2a.tsv", sep="\t").set_index("taxon_group")
    src_rec = pd.read_csv(COMP_RES / "genecall_reconcile_20261102/per_genome_reconciliation.tsv", sep="\t")
    assert int(p2a["n_genomes"].sum()) == 5438
    assert int(p2a.loc["Thermoplasmatota", "n_genomes"]) == 237
    assert int(p2a.loc["Other archaea", "n_genomes"]) == 354
    assert int(p2a.loc["Bacteria", "n_genomes"]) == 4847
    thermo = src_rec[src_rec["gtdb_classification"].str.contains("p__Thermoplasmatota")]
    assert np.isclose(p2a.loc["Thermoplasmatota", "all_dark_median_pct"], 100.0 * thermo["dark_fraction"].median())
    assert np.isclose(p2a.loc["Thermoplasmatota", "comp100_dark_median_pct"], 100.0 * thermo["complete_ge100_dark_fraction"].median())

    # Verify per_genome_reconciliation.tsv uses n_dark / n_labelled, NOT n_dark / n_genes
    assert np.allclose(src_rec["dark_fraction"], src_rec["n_dark"] / src_rec["n_labelled"], atol=1e-6)
    assert np.allclose(
        src_rec["complete_ge100_dark_fraction"],
        src_rec["complete_ge100_n_dark"] / src_rec["complete_ge100_n_labelled"],
        atol=1e-6,
    )
    assert int(src_rec["n_genes"].sum()) == 20451986
    assert int(src_rec["n_labelled"].sum()) == 20216302
    assert not np.allclose(src_rec["dark_fraction"], src_rec["n_dark"] / src_rec["n_genes"], atol=1e-4)

    p2c = pd.read_csv(PANEL_DIR / "panel_fig2c.tsv", sep="\t").set_index("grp")
    dn = pd.read_csv(COMP_RES / "dnds_strata_20261103/dnds_with_strata.tsv.gz", sep="\t", low_memory=False)
    for grp in ["control", "dark_other_phylum", "dark_Thermoplasmatota"]:
        sub = dn[dn["grp"] == grp]["median_dNdS"].dropna()
        assert int(p2c.loc[grp, "n_families"]) == len(sub)
        assert np.isclose(p2c.loc[grp, "median_dNdS"], sub.median())

    p2d = pd.read_csv(PANEL_DIR / "panel_fig2d.tsv", sep="\t").set_index("grp")
    pn = pd.read_csv(COMP_RES / "pnps_20261104/family_pnps_strata.tsv.gz", sep="\t")
    for grp in ["control", "dark_other", "dark_Thermoplasmatota"]:
        sub = pn[pn["grp"] == grp]["median_pnps"].dropna()
        assert int(p2d.loc[grp, "n_families"]) == len(sub)
        assert np.isclose(p2d.loc[grp, "median_pnps"], sub.median())


def test_fig3_and_fig4_panels_against_sources():
    p3a = pd.read_csv(PANEL_DIR / "panel_fig3a.tsv", sep="\t").set_index("group")
    assert int(p3a.loc["dark", "n_genes"]) == 237251
    assert int(p3a.loc["nondark", "n_genes"]) == 529405
    assert int(p3a.loc["dark", "core_n"] + p3a.loc["dark", "shell_n"] + p3a.loc["dark", "cloud_n"]) == 237251
    assert int(p3a.loc["nondark", "core_n"] + p3a.loc["nondark", "shell_n"] + p3a.loc["nondark", "cloud_n"]) == 529405

    p3b = pd.read_csv(PANEL_DIR / "panel_fig3b.tsv", sep="\t").set_index("group")
    t2 = pd.read_csv(COMP_RES / "sp_tm_analysis_20261104/t2_thermo_dark_vs_nondark.tsv", sep="\t")
    t2_all = t2[t2["length_bin"] == "all"].set_index("group")
    assert int(p3b.loc["thermo_dark", "n"]) == 237251
    assert int(p3b.loc["thermo_nondark", "n"]) == 529405
    assert np.isclose(p3b.loc["thermo_dark", "pct_sp_any"], t2_all.loc["thermo_dark", "pct_sp_any"])
    assert np.isclose(p3b.loc["thermo_nondark", "pct_sp_any"], t2_all.loc["thermo_nondark", "pct_sp_any"])

    p3c = pd.read_csv(PANEL_DIR / "panel_fig3c.tsv", sep="\t").set_index(["cohort", "estimator"])
    exp_adj = pd.read_csv(
        COMP_RES / "science_20261108_completed/canonical_task_b_export_stage1b_within_genome_and_genus_sensitivity.tsv",
        sep="\t",
    )
    exp_sp = exp_adj[exp_adj["endpoint"] == "sp_any"].set_index(["cohort", "estimator"])
    for key_3c, r_3c in p3c.iterrows():
        src_r = exp_sp.loc[key_3c]
        assert int(r_3c["n_independent_lineage_units"]) == int(src_r["n_independent_lineage_units"])
        assert int(r_3c["n_informative_strata"]) == int(src_r["n_informative_strata"])
        assert np.isclose(r_3c["or_estimate"], src_r["or_estimate"])
        assert np.isclose(r_3c["or_ci95_low"], src_r["or_ci95_low"])
        assert np.isclose(r_3c["or_ci95_high"], src_r["or_ci95_high"])

    p4c = pd.read_csv(PANEL_DIR / "panel_fig4c.tsv", sep="\t").iloc[0]
    cand_eval = pd.read_csv(
        COMP_RES / "science_20261108_completed/canonical_task_c_targeted_family_functional_evaluation.tsv",
        sep="\t",
    )
    # Audited 47-family deep-dive disposition (Exp 4 targeted audit v2)
    assert int(p4c["n_evaluated_families"]) == 47
    assert int(p4c["n_total_member_loci_checked"]) == 428
    assert int(p4c["n_verified_complete_interior_loci"]) == 326
    assert int(p4c["disposition_uncertain_insufficient"]) == 38
    assert int(p4c["disposition_fold_or_motif_only"]) == 9
    assert int(p4c["disposition_biochemical_validated"]) == 0
    assert pd.isna(p4c["remaining_pending_blocker"]) or str(p4c["remaining_pending_blocker"]).lower() in ("none", "nan", "")
    # Historical 19-family screening baseline
    assert int(p4c["historical_n_evaluated_families"]) == len(cand_eval) == 19
    assert int(p4c["historical_n_total_member_loci_checked"]) == int(cand_eval["n_total_member_loci"].sum()) == 295
    assert int(p4c["historical_n_verified_complete_interior_loci"]) == int(cand_eval["n_verified_complete_interior_loci"].sum()) == 217

    p4a = pd.read_csv(PANEL_DIR / "panel_fig4a.tsv", sep="\t")
    c4 = pd.read_csv(COMP_RES / "c4_c5_c7_analysis_20261104/c4/c4_summary.tsv", sep="\t")
    assert int(p4a["control_1pct_n"].sum()) == 1061
    assert int(p4a["dark_n"].sum()) == 106042
    assert int(p4a["thermo_dark_n"].sum()) == 15289
    ctrl_all = c4[(c4["arm"] == "control_1pct") & (c4["stratum"] == "all")].iloc[0]
    assert np.isclose(100.0 * ctrl_all["either_p90"], 81.53, atol=0.05)

    p4b = pd.read_csv(PANEL_DIR / "panel_fig4b.tsv", sep="\t").set_index("group_label")
    c5 = pd.read_csv(COMP_RES / "c4_c5_c7_analysis_20261104/c5/c5_family_summary.tsv", sep="\t")
    assert int(p4b.loc["Control families", "n_families"]) == 106042
    assert int(p4b.loc["All dark families", "n_families"]) == 106042
    assert int(p4b.loc["Thermo. dark families", "n_families"]) == 15289
    c5_dark = c5[(c5["arm"] == "dark") & (c5["stratum"] == "all")].iloc[0]
    assert np.isclose(p4b.loc["All dark families", "rep_domain_signature_pct"], 100.0 * c5_dark["rep_domain_signature"])
    assert np.isclose(p4b.loc["All dark families", "rep_disorder_or_coil_only_pct"], 100.0 * c5_dark["rep_disorder_or_coil_only"])

    p_c17 = pd.read_csv(PANEL_DIR / "panel_fig4_c17_struct.tsv", sep="\t").set_index("cohort")
    assert int(p_c17.loc["ESMFold_dark", "models_total"]) == 105867
    assert int(p_c17.loc["ESMFold_dark", "qualified_models"]) == 19620
    assert int(p_c17.loc["ESMFold_dark", "qualified_pdb_match"]) == 18092
    assert int(p_c17.loc["ESMFold_dark", "qualified_afdb_match_overlapping"]) == 15629
    assert int(p_c17.loc["ESMFold_dark", "qualified_no_ref_match"]) == 542
    assert int(p_c17.loc["ESMFold_control", "models_total"]) == 105996
    assert int(p_c17.loc["ESMFold_control", "qualified_models"]) == 35600
    assert int(p_c17.loc["ColabFold_Thermo_dark_rescue", "models_total"]) == 9644
    assert int(p_c17.loc["ColabFold_Thermo_dark_rescue", "qualified_models"]) == 3453
    assert int(p_c17.loc["ColabFold_Thermo_dark_rescue", "qualified_pdb_match"]) == 2720
    assert int(p_c17.loc["ColabFold_Thermo_dark_rescue", "qualified_no_ref_match"]) == 552

