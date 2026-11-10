#!/usr/bin/env python3
"""Unit and regression tests for Experiment 4 (`exp4_surface`).

Exercises all scientific contract requirements and code-review fixes:
- `mechanism-20261109/EXP4-CONTRACT.txt`
- `mechanism-20261109/EXP4-CURATED-SOURCE-NOTE.txt`
- `mechanism-20261109/EXP4-EVIDENCE-REVIEW.txt`
- `mechanism-20261109/EXP4-CODE-REVIEW.txt` (Items 1-9)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

CODE_DIR = Path(__file__).resolve().parent.parent / "code"
sys.path.insert(0, str(CODE_DIR))

from build_curated_reference_panel import (
    ALL_INITIAL_PANEL_SPECS,
    EXCLUSION_MANIFEST,
    PANEL_SPECS,
    map_pdb_polymer_to_uniprot_coordinates,
)
from exp4_surface_helpers import (
    PANEL_SPECS_BY_ID,
    apply_bh_fdr,
    classify_export_and_surface_topology,
    classify_flanking_machinery_guilds,
    classify_flanking_machinery_guilds_detailed,
    compute_sequence_integrity_and_compositional_metrics,
    evaluate_competing_roles_for_family,
    evaluate_gspa_heldout_recovery_and_leakage_guard,
    extract_tmbed_segment_spans,
    resolve_pdb_chain_and_architecture_guards,
    verify_prediction_sequence_hash_match,
)


def test_curated_panel_roles_and_hypothesized_separation() -> None:
    """Verify EXP4-EVIDENCE-REVIEW.txt, EXP4-REFERENCE-IDENTITY-REVIEW.txt, and EXP4-AUDIT-FINAL.txt."""
    by_id = {s["ref_id"]: s for s in PANEL_SPECS}

    # Exactly 17 retained verified references and 7 quarantined entries (24 total audited)
    assert len(PANEL_SPECS) == 17
    assert len(EXCLUSION_MANIFEST) == 7
    assert len(ALL_INITIAL_PANEL_SPECS) == 24

    # 1. Gaines 2025 8REY CdiPB and 8RH5 ArlB are positive ONLY for surface_assembly_pilin_archaellin
    for rid in ("REF_PIL_8REY_CdiPB", "REF_PIL_8RH5_ArlB"):
        spec = by_id[rid]
        assert spec["experimentally_supported_roles"] == "surface_assembly_pilin_archaellin"
        assert "virus_interaction_or_entry" not in spec["experimentally_supported_roles"]
        assert "virus_interaction_or_entry" in spec["hypothesized_roles"]
        assert spec["hypothesized_role_evidence_tier"] == "untested"

    # 2. True OppA binders (SSO1288, APE_1583.1) are positive ONLY for substrate_binding_or_uptake (NOT peptidase)
    for rid in ("REF_SBP_SSO1288_TrueOppA", "REF_SBP_APE1583_TrueOppA"):
        spec = by_id[rid]
        assert spec["experimentally_supported_roles"] == "substrate_binding_or_uptake"
        assert "peptidase" not in spec["experimentally_supported_roles"]

    # 3. Viral-encoded virion proteins: 3J9X is retained; 4KFR and 2W77 are quarantined
    assert "REF_VIR_3J9X_SIRV2_gp26_Coat" in by_id
    spec_3j9x = by_id["REF_VIR_3J9X_SIRV2_gp26_Coat"]
    assert spec_3j9x["role_class"] == "viral_source_entry_or_capsid_protein"
    assert spec_3j9x["source_biological_compartment"] == "viral_virion_encoded_not_host_receptor"
    assert spec_3j9x["uniprot_acc"] == "Q8V9P2"
    assert "REF_VIR_4KFR_STIV_C381_PilusBinder" in EXCLUSION_MANIFEST
    assert "REF_VIR_2W77_AFV1_ORF132_Claw" in EXCLUSION_MANIFEST

    # 4. 9UES chain A propeptide is an explicit negative control with no supported functional role
    spec_9a = by_id["REF_CTRL_9UES_ChainA_Propeptide_NegativeCtrl"]
    assert spec_9a["experimentally_supported_roles"] == "none"
    assert spec_9a["control_type"] == "negative_propeptide_chain_control"

    # 5. 7PTR is Haloferax volcanii S-layer glycoprotein csg (P25062), NOT CsxA
    spec_7ptr = by_id["REF_COAT_7PTR_CsxA_SLayer"]
    assert "Haloferax volcanii" in spec_7ptr["organism"]
    assert "csg" in spec_7ptr["gene_or_locus"]
    assert spec_7ptr["uniprot_acc"] == "P25062"

    # 6. ArlI and ArlJ are assembly machinery references (context ONLY)
    for mach_id in ("REF_PIL_OXIME000406_ArlI", "REF_PIL_OXIME000405_ArlJ"):
        spec_mach = by_id[mach_id]
        assert spec_mach["control_type"] == "assembly_machinery_reference"


def test_pdb_chain_guards_no_tmscore_ratio_inference_and_strict_tcov() -> None:
    """Verify EXP4-CODE-REVIEW.txt items 1, 4, and 7."""
    q_seq = "M" + "A" * 158  # 159 aa

    # 1. Raw 8ce8-assembly1_E hit with qtm/ttm ratio matching 159 aa MUST still remain UNRESOLVED without exact target sequence hash
    res_8ce8_raw = resolve_pdb_chain_and_architecture_guards(
        "8ce8-assembly1_E",
        q_seq,
        qtmscore=0.769,
        ttmscore=0.769,
        qcov=0.97,
        tcov=0.97,
    )
    assert res_8ce8_raw["pdb_guard_status"].startswith("unresolved_8CE8_chain_E")
    assert res_8ce8_raw["structural_fold_compatible"] is False
    assert res_8ce8_raw["has_full_mechanism_support"] is False
    assert res_8ce8_raw["is_functional_role_positive"] is False

    # 2. Canonical REF_CTRL_8CE8_Entity5_AuthE_LabelI_CcmE resolves chain metadata directly,
    #    allows fold compatibility when qcov/tcov>=0.50, but BLOCKS automatic CcmE functional transfer
    res_8ce8_ref5 = resolve_pdb_chain_and_architecture_guards(
        "REF_CTRL_8CE8_Entity5_AuthE_LabelI_CcmE",
        q_seq,
        qtmscore=0.769,
        ttmscore=0.667,
        qcov=0.95,
        tcov=0.90,
        verified_alignment_active_site_conserved=False,
        has_verified_operon_context=False,
    )
    assert res_8ce8_ref5["resolved_entity_id"] == "5"
    assert res_8ce8_ref5["resolved_label_asym_id"] == "I"
    assert res_8ce8_ref5["resolved_auth_asym_id"] == "E"
    assert res_8ce8_ref5["resolved_uniprot_acc"] == "P69490"
    assert res_8ce8_ref5["structural_fold_compatible"] is True
    assert res_8ce8_ref5["has_full_mechanism_support"] is False
    assert res_8ce8_ref5["is_functional_role_positive"] is False
    assert "blocked_automatic_CcmE_transfer" in res_8ce8_ref5["pdb_guard_status"]

    # 3. 9UES chain A propeptide is always rejected
    res_9ues_a = resolve_pdb_chain_and_architecture_guards(
        "REF_CTRL_9UES_ChainA_Propeptide_NegativeCtrl",
        "M" * 350,
        qtmscore=0.88,
        ttmscore=0.88,
        qcov=0.95,
        tcov=0.95,
    )
    assert res_9ues_a["structural_fold_compatible"] is False
    assert res_9ues_a["is_functional_role_positive"] is False

    # 4. Missing tcov (NaN) MUST fail structural_fold_compatible on 9UES chain B and 8REY
    res_9ues_b_nan_tcov = resolve_pdb_chain_and_architecture_guards(
        "REF_PEP_9UES_ChainB_TK1689_Protease",
        "M" * 420,
        qtmscore=0.82,
        ttmscore=0.80,
        qcov=0.90,
        tcov=np.nan,
    )
    assert res_9ues_b_nan_tcov["structural_fold_compatible"] is False
    assert res_9ues_b_nan_tcov["is_functional_role_positive"] is False

    # 5. Coverage alone on 9UES chain B sets structural_fold_compatible=True but is_functional_role_positive=False
    res_9ues_b_cov_only = resolve_pdb_chain_and_architecture_guards(
        "REF_PEP_9UES_ChainB_TK1689_Protease",
        "M" * 420,
        qtmscore=0.82,
        ttmscore=0.80,
        qcov=0.90,
        tcov=0.85,
        verified_alignment_active_site_conserved=False,
    )
    assert res_9ues_b_cov_only["structural_fold_compatible"] is True
    assert res_9ues_b_cov_only["has_full_mechanism_support"] is False
    assert res_9ues_b_cov_only["is_functional_role_positive"] is False


def test_missing_tmbed_and_single_helix_position_handling() -> None:
    """Verify EXP4-CODE-REVIEW.txt items 2 and 3."""
    # 1. Missing TMbed (np.nan or None) with SignalP=SP MUST NOT produce freely_secreted or intracellular
    res_nan = classify_export_and_surface_topology("SP", np.nan, np.nan, np.nan)
    assert res_nan["is_exported_any"] is True
    assert res_nan["is_freely_secreted_candidate"] is None
    assert res_nan["is_surface_or_membrane_anchored"] is None
    assert res_nan["topology_status"] == "unknown_missing_tmbed"
    assert res_nan["surface_localization_class"] == "exported_signalp_known_topology_unknown_missing_tmbed"

    # 2. SignalP=OTHER with missing TMbed MUST NOT produce intracellular_soluble
    res_other_nan = classify_export_and_surface_topology("OTHER", None, None, None)
    assert res_other_nan["is_exported_any"] is None
    assert res_other_nan["is_freely_secreted_candidate"] is None
    assert res_other_nan["surface_localization_class"] == "signalp_other_topology_unknown_missing_tmbed"

    # 3. Unknown SignalP class MUST NOT be classified as Non_secretory
    res_unk_sp = classify_export_and_surface_topology("UNKNOWN_STATUS", False, 0, 0)
    assert res_unk_sp["is_exported_any"] is None
    assert res_unk_sp["secretion_pathway_class"] == "unknown_missing_or_unrecognized_signalp"

    # 4. H=1 without segment coordinates MUST be unresolved (not assumed N-terminal anchor)
    res_h1_no_coords = classify_export_and_surface_topology("SP", True, 1, 0)
    assert res_h1_no_coords["tm_anchor_position_status"] == "unresolved_without_segment_coordinates"
    assert res_h1_no_coords["is_freely_secreted_candidate"] is None
    assert res_h1_no_coords["surface_localization_class"] == "exported_single_tm_segment_position_unresolved_signal_vs_anchor"

    # 5. H=1 WITH verified segment coordinates distinguishes C-terminal anchor vs N-terminal anchor vs cleaved SP overlap
    res_cterm = classify_export_and_surface_topology(
        "SP", True, 1, 0, aa_length=250, signalp_cs_pos=22, tm_helix_spans_1based=[(220, 242)]
    )
    assert res_cterm["surface_localization_class"] == "c_terminal_anchored_surface_coat_candidate"
    assert res_cterm["is_surface_or_membrane_anchored"] is True
    assert res_cterm["is_freely_secreted_candidate"] is False

    res_nterm = classify_export_and_surface_topology(
        "SP", True, 1, 0, aa_length=250, signalp_cs_pos=22, tm_helix_spans_1based=[(26, 46)]
    )
    assert res_nterm["surface_localization_class"] == "n_terminal_anchored_single_pass_surface"

    res_sp_overlap = classify_export_and_surface_topology(
        "SP", True, 1, 0, aa_length=250, signalp_cs_pos=25, tm_helix_spans_1based=[(4, 22)]
    )
    assert res_sp_overlap["surface_localization_class"] == "freely_secreted_extracellular_candidate"
    assert res_sp_overlap["is_freely_secreted_candidate"] is True


def test_bounded_neighbourhood_keywords_and_overlapping_sequons() -> None:
    """Verify EXP4-CODE-REVIEW.txt items 5 and 6."""
    # Substring false positives ('copper' has 'opp', 'DNA modification' has 'modification', generic 'glycosyltransferase', 'integrase')
    g_false = classify_flanking_machinery_guilds(
        "copper resistance protein CopA; DNA modification methylase; glycosyltransferase family 2; tyrosine recombinase/integrase"
    )
    assert g_false == {"generic_any_known_annotation"}

    # True CRISPR defense vs viral capsid hallmark split
    g_def = classify_flanking_machinery_guilds("CRISPR-associated endonuclease Cas1")
    assert "host_defense_antiviral" in g_def
    assert "viral_or_mobile_element_origin" not in g_def

    g_vir = classify_flanking_machinery_guilds("archaeal virion major capsid protein")
    assert "viral_or_mobile_element_origin" in g_vir
    assert "host_defense_antiviral" not in g_vir

    # Overlapping N-X-[ST] sequons in 'MNNSTAA': 'NNS' at pos 2 and 'NST' at pos 3
    comp = compute_sequence_integrity_and_compositional_metrics("MNNSTAA")
    assert comp["predicted_nxst_sequon_count"] == 2
    assert comp["predicted_nxst_sequon_positions_1based"] == "2;3"
    assert comp["glycosylation_evidence_status"] == "predicted_sequence_motif_only_not_observed_modification"


def test_gspa_adversarial_spam_predictor_and_unverified_holdout_manifest() -> None:
    """Verify EXP4-CODE-REVIEW.txt item 9:
    - Adversarial spam predictor (1000 false GO terms + 1 true GO term) gets term_precision_macro_proxy = 1/1001 (~0.001), NOT 1.0!
    - Missing/unverified training split manifest marks calibration_status = BLOCKED_* and is_calibrated_heldout = False.
    - Synthetic gap-filling demands and untested hypothesized roles raise ValueError.
    """
    gt_df = pd.DataFrame(
        [
            {
                "protein_id": "P1",
                "family_id": "FAM_1",
                "tax_order": "o__SW-10-69-26",
                "true_go_terms": "GO:0004175",
                "is_synthetic_gapfill_demand": False,
                "is_hypothesized_untested_role": False,
            }
        ]
    )
    spam_rows = [
        {
            "protein_id": "P1",
            "arm": "context",
            "predicted_go_term": "GO:0004175",
            "score": 0.90,
            "is_leaked_self_or_same_family_hit": False,
        }
    ] + [
        {
            "protein_id": "P1",
            "arm": "context",
            "predicted_go_term": f"GO:{1000000 + i:07d}",
            "score": 0.85,
            "is_leaked_self_or_same_family_hit": False,
        }
        for i in range(1000)
    ]
    claims_df = pd.DataFrame(spam_rows)

    res = evaluate_gspa_heldout_recovery_and_leakage_guard(
        claims_df,
        gt_df,
        holdout_col="family_id",
        score_threshold=0.50,
        training_split_manifest_verified=False,
    )
    assert len(res) == 1
    row = res.iloc[0]
    assert bool(row["is_calibrated_heldout"]) is False
    assert str(row["calibration_status"]).startswith("BLOCKED_unverified_pretrained_train_split")
    assert float(row["protein_any_true_term_hit_rate"]) == pytest.approx(1.0)
    assert float(row["term_precision_macro_proxy"]) == pytest.approx(1.0 / 1001.0, rel=1e-4)
    assert float(row["term_precision_micro_proxy"]) == pytest.approx(1.0 / 1001.0, rel=1e-4)

    # Synthetic gapfill demand must raise ValueError
    gt_bad = gt_df.copy()
    gt_bad["is_synthetic_gapfill_demand"] = True
    with pytest.raises(ValueError, match="Synthetic gap-filling demands"):
        evaluate_gspa_heldout_recovery_and_leakage_guard(claims_df, gt_bad)


def test_tmbed_segment_extraction_and_pipeline_imports() -> None:
    """Verify 1-based TMbed helix segment extraction, header/hit gates (`EXP4-PIPELINE-REVIEW.txt`), and CLI `--help`."""
    import subprocess
    from exp4_surface_helpers import (
        PANEL_SPECS_BY_ID,
        evaluate_competing_roles_for_family,
        extract_tmbed_segments_1based,
        qualify_structural_hit_for_member,
        resolve_curated_reference_id_from_header,
    )
    import exp4_surface_pipeline as pipe

    lab = "iiHHHHHHHHHHHHHHHHHHHHoooohhhhhhhhhhiii"
    spans = extract_tmbed_segments_1based(lab, "Hh")
    assert spans == [(3, 22), (27, 36)]
    assert extract_tmbed_segments_1based("SSSSSooooo", "Hh") == []

    tax = pipe.parse_gtdb_lineage("d__Archaea;p__Thermoplasmatota;c__E2;o__DHVEG-1;f__F1;g__G1;s__")
    assert tax["c"] == "c__E2"
    assert tax["o"] == "o__DHVEG-1"

    # Item 1: Verify every actual header in curated_role_reference_panel.faa resolves to a known PANEL_SPECS_BY_ID entry
    faa_path = Path(__file__).resolve().parent.parent / "inputs" / "curated_role_reference_panel.faa"
    if not faa_path.exists():
        faa_path = Path(__file__).resolve().parent.parent / "compact_evidence" / "curated_role_reference_panel.faa"
    with open(faa_path, "rt") as fh:
        faa_headers = [ln[1:].strip() for ln in fh if ln.startswith(">")]
    assert len(faa_headers) == 17
    for hdr in faa_headers:
        rid = resolve_curated_reference_id_from_header(hdr, require_known=True)
        assert rid in PANEL_SPECS_BY_ID
        assert PANEL_SPECS_BY_ID[rid]["ref_id"] == rid

    # Verify Foldseek target variants (.pdb, .pdb.gz, chain suffixes, lookup map, polymer sequence)
    assert (
        resolve_curated_reference_id_from_header("REF_PIL_8REY_CdiPB.pdb_A", require_known=True)
        == "REF_PIL_8REY_CdiPB"
    )
    assert (
        resolve_curated_reference_id_from_header(
            "0",
            db_lookup_map={"0": "REF_PEP_9UES_ChainB_TK1689_Protease.pdb_B"},
            require_known=True,
        )
        == "REF_PEP_9UES_ChainB_TK1689_Protease"
    )
    with pytest.raises(ValueError, match="could not be mapped"):
        resolve_curated_reference_id_from_header("UNKNOWN_REF|foo", require_known=True)

    # Item 2: Verify per-hit, per-source exact sequence & pLDDT gate blocks cluster-propagated mismatches and low pLDDT
    seq_a = "M" + "A" * 150
    seq_b = "M" + "V" * 150
    raw_hit = {
        "ref_id": "REF_PIL_8REY_CdiPB",
        "model_source": "ColabFold_AF2",
        "qtmscore": 0.68,
        "ttmscore": 0.65,
        "qcov": 0.80,
        "tcov": 0.78,
        "experimentally_supported_roles": "surface_assembly_pilin_archaellin",
    }
    # Mismatched polymer sequence (cluster representative propagation) must be rejected
    q_mismatch = qualify_structural_hit_for_member(
        raw_hit,
        member_gene="G1",
        member_cat_id="C1",
        member_aa_seq=seq_a,
        model_info={"model_polymer_aa_seq": seq_b, "plddt_0to1": 0.88},
    )
    assert q_mismatch["qualified_for_role_evidence"] is False
    assert "mismatch" in q_mismatch["rejection_reason"]

    # Low pLDDT (< 0.70) must be rejected
    q_low_plddt = qualify_structural_hit_for_member(
        raw_hit,
        member_gene="G1",
        member_cat_id="C1",
        member_aa_seq=seq_a,
        model_info={"model_polymer_aa_seq": seq_a, "plddt_0to1": 0.62},
    )
    assert q_low_plddt["qualified_for_role_evidence"] is False
    assert "plddt_below_0.70" in q_low_plddt["rejection_reason"]

    # Exact polymer match + pLDDT >= 0.70 + length >= 100 + TM >= 0.50 + cov >= 0.50 passes
    q_valid = qualify_structural_hit_for_member(
        raw_hit,
        member_gene="G1",
        member_cat_id="C1",
        member_aa_seq=seq_a,
        model_info={"model_polymer_aa_seq": seq_a, "plddt_0to1": 0.85},
    )
    assert q_valid["qualified_for_role_evidence"] is True
    assert q_valid["rejection_reason"] == "none"

    # Verify evaluate_competing_roles_for_family ignores rejected structural hits
    comp_rej = evaluate_competing_roles_for_family(
        family_id="FAM_TEST",
        candidate_role="de_novo_exported_dark_candidate",
        is_dark=True,
        surface_localization_class="prepilin_or_archaellin_filament_subunit_candidate",
        frac_exported=1.0,
        curated_seq_hits=[],
        curated_struct_hits=[q_mismatch, q_low_plddt],
        own_domain_descriptions=[],
        own_c6_hits=[],
        pdb_guard_result={"pdb_guard_status": "no_pdb_hit", "has_full_mechanism_support": False},
        role_null_results={},
        tax_classes_present={"c__E2"},
    )
    assert comp_rej["n_struct_hits_qualified"] == 0
    assert comp_rej["n_struct_hits_rejected_gate"] == 2

    # Verify CLI --help runs cleanly
    pipe_script = Path(__file__).resolve().parent.parent / "code" / "exp4_surface_pipeline.py"
    proc = subprocess.run([sys.executable, str(pipe_script), "--help"], capture_output=True, text=True, check=True)
    assert "--mode" in proc.stdout


def test_driver_directory_validation_and_search_controls(tmp_path: Path) -> None:
    """Verify require_dir vs require_file semantics and driver search controls."""
    from common_science import require_dir, require_file
    import exp4_surface_pipeline as pipe
    from exp4_surface_helpers import PANEL_SPECS_BY_ID, resolve_pdb_chain_and_architecture_guards

    # 1. require_dir succeeds on real directories and rejects files/missing
    d = tmp_path / "test_dir"
    d.mkdir()
    f = tmp_path / "test_file.txt"
    f.write_text("hello", encoding="utf-8")

    assert require_dir(d, "valid directory") == d
    with pytest.raises(ValueError, match="Expected directory"):
        require_dir(f, "file as dir")
    with pytest.raises(FileNotFoundError, match="Missing required"):
        require_dir(tmp_path / "nonexistent", "missing dir")

    # 2. require_file rejects directories
    with pytest.raises(ValueError, match="Expected regular file"):
        require_file(d, "dir as file")

    # 3. Verify actual driver inputs/pdb_chains exists and is a directory
    pdb_chains_path = Path(__file__).resolve().parent.parent / "inputs" / "pdb_chains"
    if pdb_chains_path.exists():
        assert require_dir(pdb_chains_path, "inputs/pdb_chains").is_dir()
        n_pdbs = len(list(pdb_chains_path.glob("REF_*.pdb")))
        assert n_pdbs == 11
        n_quarantined = len(list((pdb_chains_path / "quarantined").glob("REF_*.pdb")))
        assert n_quarantined == 7

    # 4. Verify confusable negative control guards tested by driver
    for neg_rid in (
        "REF_CTRL_9UES_ChainA_Propeptide_NegativeCtrl",
        "REF_CTRL_8CE8_Entity5_AuthE_LabelI_CcmE",
        "REF_CTRL_8CE8_Entity1_LabelE_Autha_CcmA",
        "REF_SBP_SSO3043_SugarConfusableCtrl",
    ):
        spec = PANEL_SPECS_BY_ID[neg_rid]
        assert spec["ref_id"] == neg_rid
        g = resolve_pdb_chain_and_architecture_guards(
            neg_rid, "M" + "A" * 160, 0.72, 0.72, 0.85, 0.85, verified_n_terminal_helix=True
        )
        assert not g["has_full_mechanism_support"]


def test_assembly_machinery_context_not_direct_pilin() -> None:
    """Verify EXP4-AUDIT-FINAL.txt Item 4:
    ArlI/J are assembly machinery context ONLY (control_type assembly_machinery_reference).
    They must be excluded from positive seq_by_role / struct_by_role direct pilin assignment,
    and flagged separately as machinery context.
    """
    arli_hit = {
        "ref_id": "REF_PIL_OXIME000406_ArlI",
        "control_type": "assembly_machinery_reference",
        "experimentally_supported_roles": "surface_assembly_pilin_archaellin",
        "qcov": 0.95,
        "tcov": 0.95,
        "evalue": 1e-40,
        "bits": 350.0,
    }

    res = evaluate_competing_roles_for_family(
        family_id="FAM_ARLI_TEST",
        candidate_role="exported_machinery",
        is_dark=False,
        surface_localization_class="exported_signalp_known_topology_unknown_missing_tmbed",
        frac_exported=1.0,
        curated_seq_hits=[arli_hit],
        curated_struct_hits=[],
        own_domain_descriptions=["Type II/IV secretion system ATPase"],
        own_c6_hits=[],
        pdb_guard_result={"pdb_guard_status": "no_pdb_hit", "has_full_mechanism_support": False},
        role_null_results={},
        tax_classes_present={"c__Thermoplasmata"},
    )

    # Must NOT count as direct pilin/archaellin subunit positive
    assert res["tier_surface_assembly_pilin_archaellin"] == "none"
    assert res["score_surface_assembly_pilin_archaellin"] == 0.0

    # Must be recorded as assembly machinery context
    assert res["has_assembly_machinery_context"] is True
    assert "REF_PIL_OXIME000406_ArlI" in res["assembly_machinery_hits"]
    assert "REF_PIL_OXIME000406_ArlI" in res["assembly_machinery_context_status"]


def test_diamond_positive_fixture_execution(tmp_path: Path) -> None:
    """Verify EXP4-AUDIT-FINAL.txt Item 6:
    Actual DIAMOND execution on non-empty valid FASTA fixture succeeds cleanly without syntax errors.
    """
    import shutil
    import subprocess

    diamond_bin = shutil.which("diamond")
    if not diamond_bin:
        pytest.skip("diamond binary not installed in test environment")

    # Create small valid FASTA with 2 reference entries
    ref_faa = tmp_path / "ref.faa"
    ref_faa.write_text(
        ">REF_TEST_1|role:surface_assembly_pilin_archaellin\n"
        "MAKQTLKLVGLGAVGAAGLVGLAAVGAGVAAAATTAG\n"
        ">REF_TEST_2|role:substrate_binding_or_uptake\n"
        "MKKIILSLALVLSACSSNDVVEITGKAPAEVKKIKVE\n",
        encoding="ascii",
    )
    dmnd_db = tmp_path / "ref.dmnd"

    # Run makedb
    cmd_makedb = [diamond_bin, "makedb", "--in", str(ref_faa), "-d", str(dmnd_db), "--quiet"]
    proc_makedb = subprocess.run(cmd_makedb, capture_output=True, text=True)
    assert proc_makedb.returncode == 0, f"DIAMOND makedb failed: {proc_makedb.stderr}"
    assert dmnd_db.with_suffix(".dmnd").exists()

    # Query with exact match
    q_faa = tmp_path / "query.faa"
    q_faa.write_text(
        ">QUERY_1\n"
        "MAKQTLKLVGLGAVGAAGLVGLAAVGAGVAAAATTAG\n",
        encoding="ascii",
    )
    out_tsv = tmp_path / "diamond_out.tsv"
    cmd_blastp = [
        diamond_bin,
        "blastp",
        "-d",
        str(dmnd_db),
        "-q",
        str(q_faa),
        "-o",
        str(out_tsv),
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qlen",
        "slen",
        "--quiet",
    ]
    proc_blastp = subprocess.run(cmd_blastp, capture_output=True, text=True)
    assert proc_blastp.returncode == 0, f"DIAMOND blastp failed: {proc_blastp.stderr}"
    assert out_tsv.exists()
    hits = out_tsv.read_text().strip().splitlines()
    assert len(hits) >= 1
    assert "REF_TEST_1" in hits[0]


def test_exact_circular_null_convolution_vs_brute_force_and_no_p_floor() -> None:
    """Verify EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt Item A and EXP4-NULL-ACCEPTANCE.txt."""
    from common_science import within_contig_circular_null

    # 1. Two-contig fixture with known brute-force enumeration
    c1_len = 8
    c2_len = 10
    contigs = ["c1"] * c1_len + ["c2"] * c2_len
    target_mask = [False] * (c1_len + c2_len)
    target_mask[2] = True  # c1 target
    target_mask[c1_len + 3] = True  # c2 target

    flank_mask = [False] * (c1_len + c2_len)
    flank_mask[1] = True  # c1 flank near target
    flank_mask[c1_len + 4] = True  # c2 flank near target

    res = within_contig_circular_null(
        contigs,
        target_mask,
        flank_mask,
        window=2,
        use_exact_convolution=True,
    )
    assert res["obs_adjacent"] == 2
    assert res["null_calculation_method"] == "exact_discrete_convolution_no_monte_carlo_floor"
    assert 0.0 < res["p_exact_one_sided"] <= 1.0

    # Brute-force verification over all 8 x 10 = 80 shift combinations
    brute_counts = []
    for s1 in range(c1_len):
        cnt1 = 1 if abs(((2 + s1) % c1_len) - 1) <= 2 and ((2 + s1) % c1_len) != 1 else 0
        for s2 in range(c2_len):
            cnt2 = 1 if abs(((3 + s2) % c2_len) - 4) <= 2 and ((3 + s2) % c2_len) != 4 else 0
            brute_counts.append(cnt1 + cnt2)
    p_brute = float(np.mean(np.array(brute_counts) >= 2))
    assert np.isclose(res["p_exact_one_sided"], p_brute, atol=1e-12)

    # 2. Rare-effect multi-contig test: p-value strictly below the 1/1001 Monte Carlo floor
    c_rare = 5
    rare_contigs = []
    rare_targets = []
    rare_flanks = []
    for idx in range(c_rare):
        cid = f"contig_{idx}"
        rare_contigs.extend([cid] * 100)
        t_vec = [False] * 100
        t_vec[50] = True
        f_vec = [False] * 100
        f_vec[51] = True  # adjacent at obs shift 0
        rare_targets.extend(t_vec)
        rare_flanks.extend(f_vec)

    res_rare = within_contig_circular_null(
        rare_contigs,
        rare_targets,
        rare_flanks,
        window=5,
        use_exact_convolution=True,
    )
    assert res_rare["obs_adjacent"] == 5
    # Chance of hitting on 1 contig is 10/100 = 0.1; across 5 contigs (0.1)^5 = 1e-5
    # Monte Carlo with 1000 perms would floor at 0.000999; exact convolution calculates ~1e-5
    assert res_rare["p_exact_one_sided"] < 0.0005, f"Expected p < 0.0005, got {res_rare['p_exact_one_sided']}"
    assert res_rare["is_numerical_underflow"] is False


def test_genus_dereplicate_map_no_retained_loci_and_missing_map() -> None:
    """Verify EXP4-NULL-ACCEPTANCE.txt critical dereplication bug fixes."""
    from stage3_helpers import run_grouped_contig_circular_nulls

    # Contig record with 1 locus in genome G1, genus block B1
    recs = [
        {
            "genome": "G1",
            "contig": "c1",
            "genus_block": "B1",
            "cand_idx": 2,
            "contig_guilds": [{"generic_any_known_annotation"} for _ in range(10)],
        }
    ]
    recs[0]["contig_guilds"][1].add("nutrient_acquisition_peptidase")

    # 1. Genus dereplicate map selects G2 for B1 (eliminating G1)
    g_map_exclude = {"B1": "G2"}
    res_elim = run_grouped_contig_circular_nulls(
        recs,
        specific_guild="nutrient_acquisition_peptidase",
        genus_dereplicate_map=g_map_exclude,
    )
    gd = res_elim["genus_dereplicated_null"]
    assert gd["status"] == "nonestimable_no_retained_loci"
    assert gd["n_target_loci"] == 0
    assert np.isnan(gd["p_exact_one_sided"])
    assert gd["n_genomes_retained"] == 0

    # 2. Genus dereplicate map is None (sensitivity not run; must NOT copy full result)
    res_none = run_grouped_contig_circular_nulls(
        recs,
        specific_guild="nutrient_acquisition_peptidase",
        genus_dereplicate_map=None,
    )
    gd_none = res_none["genus_dereplicated_null"]
    assert gd_none["status"] == "sensitivity_not_run_missing_genus_map"
    assert gd_none["n_target_loci"] == 0
    assert np.isnan(gd_none["p_exact_one_sided"])

    # 3. Unmapped genus block in genus map: excluded, NOT silently defaulted to genome_id
    res_unmapped = run_grouped_contig_circular_nulls(
        recs,
        specific_guild="nutrient_acquisition_peptidase",
        genus_dereplicate_map={"B_OTHER": "G_OTHER"},  # B1 not in map
    )
    gd_unmapped = res_unmapped["genus_dereplicated_null"]
    assert gd_unmapped["status"] == "nonestimable_no_retained_loci"
    assert gd_unmapped["n_target_loci"] == 0


def test_negative_reference_and_unrelated_context_cannot_converge() -> None:
    """Verify EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt Item 4."""
    from exp4_surface_helpers import evaluate_same_role_evidence_and_convergence

    # Row with CcmE negative control fold and high peptidase context
    row_ccme = {
        "family_id": "test_fam_ccme",
        "is_dark": True,
        "cohort_group": "candidate_route_a2_exported_dark_eligible",
        "representative_gene": "gene_1",
        "surface_localization_class": "freely_secreted_extracellular_candidate",
        "topology_status": "verified_tmbed(S=1,H=0,B=0)",
        "is_exported_any": True,
        "model_source": "colabfold",
        "model_polymer_hash_status": "exact_match",
        "model_polymer_exact_match": True,
        "model_plddt_0to1": 0.85,
        "best_curated_seq_ref_id": "none",
        "best_curated_struct_ref_id": "REF_CTRL_8CE8_Entity5_AuthE_LabelI_CcmE",
        "best_curated_struct_qtmscore": 0.82,
        "best_curated_struct_qcov": 0.90,
        "best_curated_struct_tcov": 0.90,
        "pdb_guard_status": "blocked_automatic_CcmE_transfer__resolved_REF_CTRL_8CE8_Entity5_AuthE_LabelI_CcmE_fold_only",
        "rep_own_domains": "none",
        "rep_own_c6_hits": "none",
        # Highly significant peptidase context
        "null_obs_adj_nutrient_acquisition_peptidase": 4,
        "null_n_gblocks_nutrient_acquisition_peptidase": 3,
        "joint_matrix_q_bh_nutrient_acquisition_peptidase": 0.001,
        "global_q_bh_nutrient_acquisition_peptidase": 0.0005,
    }
    conv_ccme = evaluate_same_role_evidence_and_convergence(row_ccme)
    assert conv_ccme["supported_evidence_route"] != "convergent_sequence_structure_and_same_role_context"
    assert conv_ccme["supported_evidence_route"] != "qualified_structure_and_same_role_context_with_architecture"
    assert conv_ccme["is_computational_role_hypothesis"] is True
    assert conv_ccme["is_biochemical_assignment"] is False

    # Row with 9UES Chain A propeptide negative control
    row_9ues = dict(row_ccme)
    row_9ues["best_curated_struct_ref_id"] = "REF_CTRL_9UES_ChainA_Propeptide_NegativeCtrl"
    row_9ues["pdb_guard_status"] = "rejected_9UES_chainA_propeptide_negative_architecture_control"
    conv_9ues = evaluate_same_role_evidence_and_convergence(row_9ues)
    assert conv_9ues["supported_evidence_route"] != "convergent_sequence_structure_and_same_role_context"


def test_same_role_evidence_qualification_and_separate_routes() -> None:
    """Verify EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt Item 5."""
    from exp4_surface_helpers import evaluate_same_role_evidence_and_convergence

    # Row with positive peptidase sequence match + exported architecture + peptidase context, but MISSING model
    row_pep_seq = {
        "family_id": "test_pep_seq_only",
        "is_dark": True,
        "cohort_group": "candidate_route_a2_exported_dark_eligible",
        "representative_gene": "gene_2",
        "surface_localization_class": "freely_secreted_extracellular_candidate",
        "topology_status": "verified_tmbed(S=1,H=0,B=0)",
        "is_exported_any": True,
        "model_source": "none",
        "model_polymer_hash_status": "missing_model",
        "model_polymer_exact_match": False,
        "model_plddt_0to1": np.nan,
        "best_curated_seq_ref_id": "REF_PEP_1K32_TricornProtease",
        "best_curated_seq_evalue": 1e-12,
        "best_curated_seq_qcov": 0.85,
        "best_curated_seq_tcov": 0.80,
        "best_curated_struct_ref_id": "none",
        "pdb_guard_status": "none",
        "rep_own_domains": "none",
        "rep_own_c6_hits": "none",
        "null_obs_adj_nutrient_acquisition_peptidase": 3,
        "null_n_gblocks_nutrient_acquisition_peptidase": 2,
        "joint_matrix_q_bh_nutrient_acquisition_peptidase": 0.02,
        "global_q_bh_nutrient_acquisition_peptidase": 0.01,
    }
    conv_pep = evaluate_same_role_evidence_and_convergence(row_pep_seq)
    assert conv_pep["supported_evidence_route"] == "sequence_and_same_role_context_with_architecture"
    assert conv_pep["candidate_computational_role_hypothesis"] == "nutrient_acquisition_peptidase"
    assert conv_pep["is_computational_role_hypothesis"] is True
    assert conv_pep["is_biochemical_assignment"] is False


def test_coat_guild_key_unification_and_known_coat_neighbours() -> None:
    """Verify EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt Item 2 coat guild unification."""
    from exp4_surface_helpers import classify_flanking_machinery_guilds
    from common_science import within_contig_circular_null

    desc = "surface layer glycoprotein S-layer SlaA"
    guilds = classify_flanking_machinery_guilds(desc)
    assert "cell_coat_or_envelope" in guilds

    # Null testing cell_coat_or_envelope
    contigs = ["c1"] * 10
    tgt = [False] * 10
    tgt[2] = True
    flk = [False] * 10
    flk[3] = True  # adjacent coat gene

    res = within_contig_circular_null(
        contigs,
        tgt,
        flk,
        window=2,
        use_exact_convolution=True,
    )
    assert res["obs_adjacent"] == 1
    assert res["p_exact_one_sided"] < 1.0


def test_real_v1_row_fixture_and_schema_regression() -> None:
    """Verify EXP4-SCHEMA-REGRESSION-REVIEW.txt Items 1, 2, 4, and 5:
    - Real v1 row 60PRr1::60PRr1_285908_17 with synthetic same-role context passes structure route.
    - Actual emitted model_polymer_hash_status ('exact_sequence_hash_match') accepted.
    - Structured topology verified_tmbed validated; missing stays unknown and ineligible.
    - Domain profile match without direct sequence ref is domain_profile_and_same_role_context_with_architecture.
    - Role-specific hits from detail_{r} audited even when overall best hit is another role or negative control.
    """
    from exp4_surface_helpers import evaluate_same_role_evidence_and_convergence

    # 1. Real v1 row 60PRr1::60PRr1_285908_17 parameters from frozen baseline
    real_row_v1 = {
        "family_id": "60PRr1::60PRr1_285908_17",
        "is_dark": True,
        "cohort_group": "candidate_route_a2_exported_dark_eligible",
        "representative_gene": "60PRr1__single_SemiBin_79::60PRr1_285908_17",
        "surface_localization_class": "n_terminal_anchored_single_pass_surface",
        "topology_status": "verified_tmbed(S=1,H=1[22-37],B=0)",
        "tm_anchor_position_status": "verified_n_terminal_tm_anchor[22-37]",
        "is_exported_any": True,
        "model_source": "ESMFold",
        "model_plddt_0to1": 0.7463,
        "model_polymer_hash_status": "exact_sequence_hash_match",  # actual emitted string!
        "model_polymer_exact_match": True,
        "best_curated_seq_ref_id": "REF_PIL_8REY_CdiPB",
        "best_curated_seq_evalue": 1.05e-06,
        "best_curated_seq_qcov": 0.155,  # qcov < 0.50 -> does not qualify for direct seq ref
        "best_curated_seq_tcov": 0.219,
        "best_curated_struct_ref_id": "REF_PIL_5O4U_FlaB0",
        "best_curated_struct_qtmscore": 0.561,
        "best_curated_struct_ttmscore": 0.5588,
        "best_curated_struct_qcov": 0.922,
        "best_curated_struct_tcov": 0.923,
        "pdb_guard_status": "structural_and_n_term_helix_compatible_REF_PIL_5O4U_FlaB0_machinery_context_required_outside_Thermoplasmatales",
        "rep_own_domains": "Pfam:PF07790(Archaeal Type IV pilin, N-terminal)[18-67];HHsearch:Pilin_N(missing)",
        "rep_own_c6_hits": "none",
        "sbp_substrate_specificity_status": "no_sbp_match",
        "viral_origin_audit_status": "no_viral_virion_similarity_or_hallmark_flank",
        "max_role_evidence_score": 0.45,
        "top_competing_roles": "surface_assembly_pilin_archaellin",
    }

    # In actual v1 data, context p is 1.0 (obs_adj=0)
    actual_v1 = dict(real_row_v1)
    actual_v1["null_obs_adj_surface_assembly_pilin_archaellin"] = 0
    actual_v1["null_n_gblocks_surface_assembly_pilin_archaellin"] = 0
    actual_v1["joint_matrix_q_bh_surface_assembly_pilin_archaellin"] = 1.0
    actual_res = evaluate_same_role_evidence_and_convergence(actual_v1)
    assert actual_res["supported_evidence_route"] == "exploratory_fold_or_motif_candidate_no_operon_convergence"
    assert actual_res["candidate_computational_role_hypothesis"] == "surface_assembly_pilin_archaellin"
    assert actual_res["is_computational_role_hypothesis"] is True

    # With DELIBERATELY SYNTHETIC same-role context (q=.01, obs=3, genus=3) as specified by root review
    synth_v1 = dict(real_row_v1)
    synth_v1["null_obs_adj_surface_assembly_pilin_archaellin"] = 3
    synth_v1["null_n_gblocks_surface_assembly_pilin_archaellin"] = 3
    synth_v1["joint_matrix_q_bh_surface_assembly_pilin_archaellin"] = 0.01
    synth_res = evaluate_same_role_evidence_and_convergence(synth_v1)
    assert synth_res["supported_evidence_route"] == "qualified_structure_and_same_role_context_with_architecture"
    assert synth_res["candidate_computational_role_hypothesis"] == "surface_assembly_pilin_archaellin"
    assert synth_res["has_direct_positive_struct_ref_match"] is True
    assert synth_res["matched_positive_struct_ref_id"] == "REF_PIL_5O4U_FlaB0"

    # 2. Missing TMbed topology stays unknown and ineligible
    unver_topo = dict(synth_v1)
    unver_topo["topology_status"] = "unknown_missing_tmbed"
    unver_res = evaluate_same_role_evidence_and_convergence(unver_topo)
    assert unver_res["supported_evidence_route"] != "qualified_structure_and_same_role_context_with_architecture"

    # 3. Domain profile alone + context yields domain_profile route, NOT direct sequence ref route
    dom_only = {
        "family_id": "test_dom_only",
        "is_dark": True,
        "cohort_group": "candidate_route_a2_exported_dark_eligible",
        "surface_localization_class": "freely_secreted_extracellular_candidate",
        "topology_status": "verified_tmbed(S=1,H=0,B=0)",
        "is_exported_any": True,
        "model_source": "none",
        "model_polymer_hash_status": "missing_model",
        "model_polymer_exact_match": False,
        "model_plddt_0to1": np.nan,
        "best_curated_seq_ref_id": "none",
        "best_curated_struct_ref_id": "none",
        "pdb_guard_status": "none",
        "rep_own_domains": "Pfam:PF00082(Subtilisin-like protease)",
        "rep_own_c6_hits": "none",
        "null_obs_adj_nutrient_acquisition_peptidase": 3,
        "null_n_gblocks_nutrient_acquisition_peptidase": 2,
        "joint_matrix_q_bh_nutrient_acquisition_peptidase": 0.02,
        "global_q_bh_nutrient_acquisition_peptidase": 0.01,
    }
    dom_res = evaluate_same_role_evidence_and_convergence(dom_only)
    assert dom_res["supported_evidence_route"] == "domain_profile_and_same_role_context_with_architecture"
    assert dom_res["has_domain_profile_match"] is True
    assert dom_res["has_direct_positive_seq_ref_match"] is False

    # 4. Role-specific hit from detail_{r} inspected even if overall best hit is CcmE negative control
    masked_hit = dict(synth_v1)
    masked_hit["best_curated_struct_ref_id"] = "REF_CTRL_8CE8_Entity5_AuthE_LabelI_CcmE"
    masked_hit["detail_surface_assembly_pilin_archaellin"] = (
        "tier=tier2_qualified_architecture_and_machinery_context;"
        "seq_ref=none;struct_ref=REF_PIL_5O4U_FlaB0;ctx_adj=3(p=0.001,q=0.01)"
    )
    masked_res = evaluate_same_role_evidence_and_convergence(masked_hit)
    assert masked_res["supported_evidence_route"] == "qualified_structure_and_same_role_context_with_architecture"
    assert masked_res["candidate_computational_role_hypothesis"] == "surface_assembly_pilin_archaellin"


def _get_followup_review_dir() -> Path:
    res_dir = Path(__file__).resolve().parent.parent / "results" / "followup_v2_review"
    if not res_dir.exists():
        res_dir = Path(__file__).resolve().parent.parent / "compact_evidence"
    return res_dir


def test_followup_audited_model_availability_partition() -> None:
    """Verify EXP4-V2-RESULTS-AND-TARGETED-AUDIT.txt Item 4: Model availability partition counters.

    Asserts that exact, missing, and mismatch partition predicates sum exactly to cohort totals
    and globally to 1,460 (565 exact / 802 missing / 93 mismatch).
    """
    res_dir = _get_followup_review_dir()
    cohort_p = res_dir / "exp4_cohort_architecture_and_role_comparison_v2_audited.tsv"
    summary_p = res_dir / "exp4_surface_final_summary_v2_audited.json"

    if not cohort_p.exists() or not summary_p.exists():
        pytest.skip("Followup audited review files not present locally.")

    import json
    with open(summary_p) as fh:
        summary = json.load(fh)

    global_part = summary["model_availability_partition_global"]
    assert global_part["partition_valid"] is True
    assert global_part["exact_sequence_hash_match"] == 565
    assert global_part["missing_prediction_not_negative"] == 802
    assert global_part["invalid_sequence_mismatch"] == 93
    assert global_part["sum_partition"] == 1460

    cohort_df = pd.read_csv(cohort_p, sep="\t")
    assert cohort_df["n_families"].sum() == 1460
    assert cohort_df["n_models_exact_match"].sum() == 565
    assert cohort_df["n_models_missing"].sum() == 802
    assert cohort_df["n_models_mismatch"].sum() == 93

    for _, row in cohort_df.iterrows():
        total = row["n_families"]
        part_sum = row["n_models_exact_match"] + row["n_models_missing"] + row["n_models_mismatch"]
        assert part_sum == total, f"Cohort {row['cohort_group']} partition sum {part_sum} != {total}"


def test_followup_audited_genus_joint_fdr_sensitivity() -> None:
    """Verify EXP4-V2-RESULTS-AND-TARGETED-AUDIT.txt Items 2 & 5: Full-universe genus joint BH-FDR.

    Asserts that across all 8,760 tests (missing p treated as 1.0), min q == 1.0000.
    """
    res_dir = _get_followup_review_dir()
    eval_p = res_dir / "exp4_surface_family_evaluation_v2_audited.tsv"

    if not eval_p.exists():
        pytest.skip("Followup audited evaluation table not present locally.")

    df = pd.read_csv(eval_p, sep="\t")
    assert len(df) == 1460

    null_guilds = [
        "nutrient_acquisition_peptidase",
        "substrate_binding_or_uptake",
        "surface_assembly_pilin_archaellin",
        "cell_coat_or_envelope",
        "host_defense_antiviral",
        "viral_or_mobile_element_origin",
    ]

    for g in null_guilds:
        assert f"genus_derep_joint_matrix_q_bh_{g}" in df.columns

    min_q = df[[f"genus_derep_joint_matrix_q_bh_{g}" for g in null_guilds]].min().min()
    assert min_q >= 1.0 - 1e-9, f"Expected min genus joint q >= 1.0, got {min_q}"

    # Verify honest n_aligned_sequences: if deferred, must be 0!
    deferred = df[df["msa_status"] == "deferred_non_focal_no_role_hit"]
    assert (deferred["n_aligned_sequences"] == 0).all(), "Deferred MSAs must report 0 aligned sequences"


def test_followup_audited_lead_evidence_and_disposition() -> None:
    """Verify EXP4-TARGETED-AUDIT-CLAIM-BOUNDARIES.txt & EXP4-FINAL-PROSE-CORRECTION.txt: Bounded inference & exact checks for leads."""
    res_dir = _get_followup_review_dir()
    dossiers_p = res_dir / "exp4_47_candidate_dossiers.tsv"

    if not dossiers_p.exists():
        pytest.skip("Followup candidate dossiers not present locally.")

    df = pd.read_csv(dossiers_p, sep="\t")

    # Lead 1: 22PRr1::22PRr1_2531955_3
    r1 = df[df["family_id"] == "22PRr1::22PRr1_2531955_3"].iloc[0]
    assert r1["primary_lead_guild"] == "surface_assembly_pilin_archaellin"
    assert bool(r1["primary_significant_at_0_10"]) is True
    assert bool(r1["primary_significant_at_0_05"]) is True
    assert bool(r1["genus_significant_at_0_10"]) is False
    assert bool(r1["genus_significant_at_0_05"]) is False
    assert r1["genus_joint_q_bh"] >= 1.0 - 1e-9

    # Distinguish seed vs verified interior representative
    assert "22PRr1_2531955_3" in str(r1["seed_gene_id"])
    assert bool(r1["seed_is_verified_interior"]) is False  # contig-edge
    assert "1PRr2" in str(r1["representative_gene"])  # longest verified interior member
    assert bool(r1["rep_is_verified_interior"]) is True
    assert r1["rep_model_verification_status"] == "no_model_available"

    # Distinct counts
    assert int(r1["n_observed_flanking_loci"]) == 13
    assert int(r1["n_context_genus_blocks"]) == 6
    assert int(r1["n_genus_derep_retained_loci"]) == 5
    assert int(r1["n_genus_derep_retained_genus_blocks"]) == 5

    # Narrative guards: root-reviewed disposition wording from EXP4-FINAL-PROSE-CORRECTION.txt!
    disp1 = str(r1["scientific_disposition_summary"])
    assert "CRITICAL FALSIFICATION" not in disp1
    assert "absence of structural homology" not in disp1
    assert "completely lacks" not in disp1
    assert "experimental AlphaFold" not in disp1
    assert "genus-restricted sensitivity does not" in disp1
    assert "F37Dr2" in disp1
    assert "S20PRr1" in disp1
    assert "PilJ" in disp1
    assert "Co-transcription" in disp1
    assert "archaeal pilin identity" in disp1

    # Lead 2: V60PRr2::V60PRr2_1314979_4
    r2 = df[df["family_id"] == "V60PRr2::V60PRr2_1314979_4"].iloc[0]
    assert r2["primary_lead_guild"] == "cell_coat_or_envelope"
    assert bool(r2["primary_significant_at_0_10"]) is True
    assert bool(r2["primary_significant_at_0_05"]) is True
    assert bool(r2["genus_significant_at_0_10"]) is False
    assert bool(r2["genus_significant_at_0_05"]) is False
    assert r2["genus_joint_q_bh"] >= 1.0 - 1e-9

    # Distinguish seed vs verified interior representative
    assert "1314979_4" in str(r2["seed_gene_id"])
    assert bool(r2["seed_is_verified_interior"]) is False  # contig-edge
    assert bool(r2["rep_is_verified_interior"]) is True
    assert r2["rep_model_verification_status"] == "no_model_available"

    # Distinct counts
    assert int(r2["n_observed_flanking_loci"]) == 11
    assert int(r2["n_context_genus_blocks"]) == 5
    assert int(r2["n_genus_derep_retained_loci"]) == 4
    assert int(r2["n_genus_derep_retained_genus_blocks"]) == 4

    # Narrative guards: root-reviewed disposition wording from EXP4-FINAL-PROSE-CORRECTION.txt!
    disp2 = str(r2["scientific_disposition_summary"])
    assert "CRITICAL FALSIFICATION" not in disp2
    assert "absence of structural homology" not in disp2
    assert "completely lacks" not in disp2
    assert "experimental AlphaFold" not in disp2
    assert "genus-restricted sensitivity does not" in disp2
    assert "must not be interpreted as absence of homology" in disp2
    assert "Co-transcription" in disp2
    assert "S-layer membership" in disp2


def test_followup_audited_47_candidate_dossiers_consistency() -> None:
    """Verify EXP4-TARGETED-AUDIT-CLAIM-BOUNDARIES.txt Item 3: Honest 47 candidate dossiers."""
    res_dir = _get_followup_review_dir()
    dossiers_p = res_dir / "exp4_47_candidate_dossiers.tsv"

    if not dossiers_p.exists():
        pytest.skip("Followup candidate dossiers not present locally.")

    df = pd.read_csv(dossiers_p, sep="\t")
    assert len(df) == 47

    # Honest sequence count vs aligned count
    deferred = df[df["msa_status"] == "deferred_non_focal_no_role_hit"]
    assert (deferred["n_aligned_sequences"] == 0).all()
    assert (deferred["n_member_sequences"] > 0).all()

    # Never claim biochemical assignment
    assert (df["is_biochemical_assignment"] == False).all()


def test_followup_audited_member_and_flank_evidence_tables() -> None:
    """Verify EXP4-TARGETED-AUDIT-CLAIM-BOUNDARIES.txt Items 1, 3, & 5: Raw evidence tables."""
    res_dir = _get_followup_review_dir()
    mem_p = res_dir / "exp4_47_member_evidence.tsv"
    flank_p = res_dir / "exp4_47_flank_evidence.tsv"

    if not mem_p.exists() or not flank_p.exists():
        pytest.skip("Member or flank evidence table not present locally.")

    m_df = pd.read_csv(mem_p, sep="\t")
    f_df = pd.read_csv(flank_p, sep="\t")

    assert len(m_df) == 428
    assert len(f_df) == 2883

    # Reject identity-by-length: is_exact_model_match must require exact hash match status
    for _, row in m_df.iterrows():
        if row["is_exact_model_match"]:
            assert row["best_model_verification_status"] == "exact_sequence_hash_match"
        if row["is_verified_interior_with_exact_model"]:
            assert row["verified_complete_interior_member"] is True
            assert row["is_exact_model_match"] is True

    # Coordinate audit: coordinates incorporated where authoritative inputs exist; missing remain unknown
    assert "coordinate_provenance" in f_df.columns
    assert "annotated_guild_support" in f_df.columns
    assert (f_df["machinery_experimental_verification_status"] == "unverified_computational_annotation_confidence_unknown").all()

    # Lead 1 and Lead 2 context pairs: all 24 verified immediate adjacent, cooriented, tight intergenic gaps
    lead1_flanks = f_df[(f_df["family_id"] == "22PRr1::22PRr1_2531955_3") & (f_df["flanking_specific_guilds"].str.contains("surface_assembly_pilin_archaellin"))]
    lead2_flanks = f_df[(f_df["family_id"] == "V60PRr2::V60PRr2_1314979_4") & (f_df["flanking_specific_guilds"].str.contains("cell_coat_or_envelope"))]
    assert len(lead1_flanks) == 13
    assert len(lead2_flanks) == 11

    for l_df, gap_min, gap_max in [(lead1_flanks, -4, 37), (lead2_flanks, -4, 6)]:
        assert (l_df["same_strand"].astype(str) == "True").all()
        assert (l_df["relative_gene_offset"].abs() == 1).all()
        gaps = l_df["intergenic_distance_bp"].astype(int).tolist()
        assert min(gaps) >= gap_min and max(gaps) <= gap_max
        assert (l_df["coordinate_provenance"].str.startswith("verified_loci_input_join:sha256:")).all()

    # Unverified / missing flank coordinates must remain honestly unknown
    missing_flanks = f_df[f_df["coordinate_provenance"] == "unknown_not_in_verified_loci_input"]
    assert len(missing_flanks) == 2883 - 346
    assert (missing_flanks["intergenic_distance_bp"] == "unknown_not_in_verified_loci_input").all()


def test_independent_polymer_and_msa_replay_audit() -> None:
    """Verify EXP4-COMPLETE-REPRODUCIBILITY.txt Directive 1: Independent replay audit."""
    res_dir = _get_followup_review_dir()
    replay_json_p = res_dir / "exp4_independent_polymer_and_msa_replay_audit.json"
    replay_tsv_p = res_dir / "exp4_428_independent_polymer_replay.tsv"

    if not replay_json_p.exists() or not replay_tsv_p.exists():
        pytest.skip("Independent replay audit artifacts not present locally.")

    import json
    with open(replay_json_p) as fh:
        audit = json.load(fh)

    assert audit["n_shard_bundles_audited"] == 17
    assert audit["all_17_shards_present"] is True
    assert audit["n_deep_dive_families"] == 47
    assert audit["n_members_replayed"] == 428
    assert audit["n_members_exact_sequence_hash_match"] == 32
    assert audit["n_verified_interior_members_exact_match"] == 27
    assert audit["member_evidence_agreement_pass"] is True
    assert audit["all_msa_taxa_valid"] is True
    assert audit["all_ungapped_sequences_byte_match"] is True

    r_df = pd.read_csv(replay_tsv_p, sep="\t")
    assert len(r_df) == 428


def test_followup_audited_neighbourhood_panels_coverage() -> None:
    """Verify EXP4-TARGETED-AUDIT-CLAIM-BOUNDARIES.txt Item 2: Full neighbourhood extraction for leads."""
    res_dir = _get_followup_review_dir()
    neigh_p = res_dir / "exp4_surface_neighbourhood_panels_v2_audited.tsv"

    if not neigh_p.exists():
        pytest.skip("Followup audited neighbourhood panels not present locally.")

    df = pd.read_csv(neigh_p, sep="\t")
    n_22 = (df["family_id"] == "22PRr1::22PRr1_2531955_3").sum()
    n_v60 = (df["family_id"] == "V60PRr2::V60PRr2_1314979_4").sum()

    assert n_22 == 115, f"Expected 115 neighbourhood rows for 22PRr1, got {n_22}"
    assert n_v60 == 140, f"Expected 140 neighbourhood rows for V60PRr2, got {n_v60}"


def test_root_reviewed_lead_interpretations_and_independent_tree_audit() -> None:
    """Verify EXP4-FINAL-PROSE-CORRECTION.txt: Root reviewed interpretations and independent tree leaf audits."""
    res_dir = _get_followup_review_dir()
    interp_p = res_dir / "exp4_root_reviewed_lead_interpretations.json"
    tree_audit_p = res_dir / "exp4_independent_tree_leaf_sequence_audit.json"
    crosswalk_p = res_dir / "exp4_tree_leaf_name_crosswalk.tsv"
    join_p = res_dir / "exp4_replay_member_identity_join.json"

    if not interp_p.exists() or not tree_audit_p.exists():
        pytest.skip("Root reviewed artifacts not present locally.")

    import json
    with open(interp_p) as fh:
        interp = json.load(fh)
    assert "22PRr1::22PRr1_2531955_3" in interp
    assert "Ten partners belong to F37Dr2" in interp["22PRr1::22PRr1_2531955_3"]
    assert "three to S20PRr1" in interp["22PRr1::22PRr1_2531955_3"]

    with open(tree_audit_p) as fh:
        tree_audit = json.load(fh)
    assert len(tree_audit) == 3
    for entry in tree_audit:
        assert entry["checks"]["alignment_ids_match_members"] is True
        assert entry["checks"]["tree_leaves_unique"] is True
        assert entry["checks"]["sequence_hashes_match_raw_member_replay"] is True
        assert entry["checks"]["tree_leaf_ids_match_explicit_unique_colon_to_underscore_map"] is True

    cw_df = pd.read_csv(crosswalk_p, sep="\t")
    assert len(cw_df) == 45
    assert (cw_df["tree_leaf_id"] == cw_df["original_gene_id"].str.replace(":", "_")).all()

    with open(join_p) as fh:
        join_data = json.load(fh)
    assert join_data["n_members"] == 428
    assert join_data["checks"]["exact_member_identity_join"] is True
    assert join_data["checks"]["sequence_hash_agreement"] is True


def test_fresh_root_preflight_refuses_historical_and_existing_outputs(tmp_path: Path) -> None:
    """Verify EXP4-REANALYSIS-ENTRYPOINT-REVIEW.txt:
    Entrypoint wrapper preflight must refuse execution in historical root,
    refuse execution if pre-existing outputs or .done markers exist (no silent skip),
    and refuse execution if required input bundles are missing.
    """
    import os
    import subprocess

    array_script = CODE_DIR / "sbatch_exp4_array.sh"
    collect_script = CODE_DIR / "sbatch_exp4_collect.sh"
    audit_script = CODE_DIR / "sbatch_exp4_followup_audit.sh"

    assert array_script.exists()
    assert collect_script.exists()
    assert audit_script.exists()

    base_env = dict(os.environ)

    # 1. Array wrapper refuses execution in historical root
    env_hist = {**base_env, "EQ_EXP4_ROOT": "/scratch/project/k10018/eq-mechanism-20261109/exp4_surface"}
    proc_hist = subprocess.run(["bash", str(array_script)], env=env_hist, capture_output=True, text=True)
    assert proc_hist.returncode != 0
    assert "Refusing execution in historical root" in proc_hist.stderr

    # 2. Array wrapper refuses pre-existing .done marker (does NOT silently skip)
    run_root = tmp_path / "fresh_run_root"
    shards_dir = run_root / "results" / "shards"
    shards_dir.mkdir(parents=True, exist_ok=True)
    done_marker = shards_dir / "shard_0.v2.done"
    done_marker.touch()

    env_done = {**base_env, "EQ_EXP4_ROOT": str(run_root), "SLURM_ARRAY_TASK_ID": "0"}
    proc_done = subprocess.run(["bash", str(array_script)], env=env_done, capture_output=True, text=True)
    assert proc_done.returncode != 0
    assert "Pre-existing output" in proc_done.stderr or "will not silently skip" in proc_done.stderr

    # 3. Array wrapper refuses pre-existing evaluation output
    done_marker.unlink()
    eval_out = shards_dir / "shard_0_family_evaluation.tsv"
    eval_out.write_text("family_id\tguild\n")

    proc_eval = subprocess.run(["bash", str(array_script)], env=env_done, capture_output=True, text=True)
    assert proc_eval.returncode != 0
    assert "Pre-existing output" in proc_eval.stderr

    # 4. Array wrapper refuses missing input shard bundle
    eval_out.unlink()
    proc_no_in = subprocess.run(["bash", str(array_script)], env=env_done, capture_output=True, text=True)
    assert proc_no_in.returncode != 0
    assert "Required input bundle missing" in proc_no_in.stderr

    # 5. Collector wrapper refuses historical root
    proc_col_hist = subprocess.run(["bash", str(collect_script)], env=env_hist, capture_output=True, text=True)
    assert proc_col_hist.returncode != 0
    assert "Refusing execution in historical root" in proc_col_hist.stderr

    # 6. Collector wrapper refuses pre-existing collector output (actual v2 filename)
    col_out = run_root / "results" / "exp4_surface_family_evaluation_v2.tsv"
    col_out.parent.mkdir(parents=True, exist_ok=True)
    col_out.write_text("family_id\tguild\n")
    proc_col_pre = subprocess.run(["bash", str(collect_script)], env=env_done, capture_output=True, text=True)
    assert proc_col_pre.returncode != 0
    assert "Pre-existing output" in proc_col_pre.stderr

    # 7. Followup audit wrapper refuses historical root
    proc_audit_hist = subprocess.run(["bash", str(audit_script)], env=env_hist, capture_output=True, text=True)
    assert proc_audit_hist.returncode != 0
    assert "Refusing execution in historical root" in proc_audit_hist.stderr


def test_historical_root_guard_path_normalization(tmp_path: Path) -> None:
    """Verify EXP4-ENTRYPOINT-ACTUAL-OUTPUT-FIX.txt:
    Historical root guard must not be bypassed by trailing slash, dot segments, or symlinks.
    """
    import os
    import subprocess

    array_script = CODE_DIR / "sbatch_exp4_array.sh"
    collect_script = CODE_DIR / "sbatch_exp4_collect.sh"
    base_env = dict(os.environ)

    # 1. Trailing slash bypass attempt
    env_slash = {**base_env, "EQ_EXP4_ROOT": "/scratch/project/k10018/eq-mechanism-20261109/exp4_surface/"}
    p_slash = subprocess.run(["bash", str(array_script)], env=env_slash, capture_output=True, text=True)
    assert p_slash.returncode != 0
    assert "Refusing execution in historical root" in p_slash.stderr

    # 2. Dot-segment /.. bypass attempt
    env_dots = {**base_env, "EQ_EXP4_ROOT": "/scratch/project/k10018/foo/../eq-mechanism-20261109/exp4_surface"}
    p_dots = subprocess.run(["bash", str(collect_script)], env=env_dots, capture_output=True, text=True)
    assert p_dots.returncode != 0
    assert "Refusing execution in historical root" in p_dots.stderr

    # 3. Symlink bypass attempt
    symlink_root = tmp_path / "symlink_to_historical"
    symlink_root.symlink_to("/scratch/project/k10018/eq-mechanism-20261109/exp4_surface")
    env_sym = {**base_env, "EQ_EXP4_ROOT": str(symlink_root)}
    p_sym = subprocess.run(["bash", str(array_script)], env=env_sym, capture_output=True, text=True)
    assert p_sym.returncode != 0
    assert "Refusing execution in historical root" in p_sym.stderr


def test_collector_source_to_wrapper_contract_and_postcondition(tmp_path: Path) -> None:
    """Verify EXP4-ENTRYPOINT-ACTUAL-OUTPUT-FIX.txt:
    1. Collector wrapper must check the actual writer output filename (exp4_surface_family_evaluation_v2.tsv),
       not an invented filename.
    2. Collector wrapper postcondition succeeds when actual collector writes the v2 output.
    """
    import os
    import subprocess

    pipeline_code = (CODE_DIR / "exp4_surface_pipeline.py").read_text()
    collect_script = CODE_DIR / "sbatch_exp4_collect.sh"
    collect_code = collect_script.read_text()
    audit_code = (CODE_DIR / "sbatch_exp4_followup_audit.sh").read_text()

    # Source contract assertions:
    # 1. Pipeline writes exp4_surface_family_evaluation_v2.tsv
    assert 'all_eval.to_csv(P_RES / "exp4_surface_family_evaluation_v2.tsv"' in pipeline_code
    # 2. No invented exp4_surface_all_evaluated_families.tsv in wrappers or pipeline
    assert "exp4_surface_all_evaluated_families.tsv" not in pipeline_code
    assert "exp4_surface_all_evaluated_families.tsv" not in collect_code
    assert "exp4_surface_all_evaluated_families.tsv" not in audit_code
    # 3. Collector script checks actual v2 file
    assert 'COLLECT_OUT="$EQ_EXP4_ROOT/results/exp4_surface_family_evaluation_v2.tsv"' in collect_code
    # 4. Followup audit script checks actual v2 file
    assert 'COLLECT_EVAL="$EQ_EXP4_ROOT/results/exp4_surface_family_evaluation_v2.tsv"' in audit_code

    # Postcondition execution test with mocked collector:
    mock_root = tmp_path / "mock_collect_run"
    mock_code = mock_root / "code"
    mock_res = mock_root / "results"
    mock_shards = mock_res / "shards"
    mock_logs = mock_root / "logs"

    mock_code.mkdir(parents=True, exist_ok=True)
    mock_shards.mkdir(parents=True, exist_ok=True)
    mock_logs.mkdir(parents=True, exist_ok=True)

    # Stage 17 shard evaluation files
    for s_id in range(17):
        (mock_shards / f"shard_{s_id}_family_evaluation.tsv").write_text("family_id\tguild\n")

    # Mock exp4_surface_pipeline.py that writes actual v2 collector output
    mock_pipeline = mock_code / "exp4_surface_pipeline.py"
    mock_pipeline.write_text(
        "import sys, os\n"
        "from pathlib import Path\n"
        "root = Path(os.environ['EQ_EXP4_ROOT'])\n"
        "(root / 'results' / 'exp4_surface_family_evaluation_v2.tsv').write_text('family_id\\tguild\\n')\n"
        "print('[mock collect] wrote exp4_surface_family_evaluation_v2.tsv')\n"
    )

    # Copy sbatch_exp4_collect.sh into mock_code
    (mock_code / "sbatch_exp4_collect.sh").write_text(collect_code)

    env = {**os.environ, "EQ_EXP4_ROOT": str(mock_root), "SLURM_JOB_ID": "mock123"}
    proc = subprocess.run(["bash", str(mock_code / "sbatch_exp4_collect.sh")], env=env, capture_output=True, text=True)
    assert proc.returncode == 0, f"Collector failed unexpectedly: stderr={proc.stderr}"
    assert "Completed successfully" in proc.stdout
    assert (mock_res / "exp4_surface_family_evaluation_v2.tsv").exists()


def test_readme_and_fetch_instructions_staging_contract() -> None:
    """Verify EXP4-ENTRYPOINT-ACTUAL-OUTPUT-FIX.txt:
    README and fetch instructions must stage package code from LOCAL via rsync/scp (-J ws-ts),
    must NOT attempt nonexistent laptop paths on Shaheen login,
    must use staged relative manifest for pickle hash verification,
    and must stage/hash bulk data inside Slurm preflight rather than raw login node.
    """
    pkg_root = CODE_DIR.parent
    readme_text = (pkg_root / "README.md").read_text()
    fetch_text = (pkg_root / "manifests" / "institutional_fetch_instructions.txt").read_text()

    for doc in (readme_text, fetch_text):
        # Must document staging from local with jump host
        assert "rsync" in doc or "scp" in doc
        assert "ws-ts" in doc

        # Must not have ssh shaheen followed by cp /home/leechuck/Public
        # Check that after ssh shaheen, no command references /home/leechuck/Public
        shaheen_idx = doc.find("ssh shaheen")
        if shaheen_idx != -1:
            remote_section = doc[shaheen_idx:]
            # Ensure no local laptop paths in remote execution block
            assert "/home/leechuck/Public" not in remote_section, (
                "Found local laptop path /home/leechuck/Public inside remote Shaheen execution instructions"
            )

        # Must verify relative or staged manifest path on Shaheen
        assert "$EQ_EXP4_ROOT/manifests/upstream_shard_bundles.sha256" in doc or "manifests/upstream_shard_bundles.sha256" in doc

        # Must document Slurm preflight / compute allocation for bulk staging & hashing
        assert "preflight" in doc.lower() or "srun" in doc or "sbatch" in doc






