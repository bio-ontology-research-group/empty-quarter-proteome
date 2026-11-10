#!/usr/bin/env python3
"""Unit and scientific gate tests for Stage 3 candidate functional evaluation.

Covers all requirements in:
- `coordinated-revision-20261108/stage3-scientific-review.txt` (Item 6)
- `coordinated-revision-20261108/stage1b-final-review.txt` (Items 1-5)
"""
from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

CODE_DIR = Path(__file__).resolve().parent.parent / "code"
sys.path.insert(0, str(CODE_DIR))

from common_science import (  # noqa: E402
    PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES,
    scale_plddt_explicit_source,
)
from stage3_helpers import (  # noqa: E402
    STAGE1B_VERIFIED_PROPOSALS_REQUIRED_COLS,
    audit_gene_table_syntax_and_fasta,
    build_gspa_pilot_3genome_audit_table,
    build_manifest_genome_alias_resolver,
    calibrate_leave_order_out_threshold,
    classify_family_control_or_candidate_role,
    compute_alignment_summary_metrics,
    evaluate_reference_motif_conservation,
    evaluate_seven_evidence_axes,
    extract_verified_route_a1_reaction_pathway_context,
    is_exact_integer_value,
    qualify_structural_and_domain_evidence_strict,
    reconcile_and_summarize_order_metabolic_tables,
    run_grouped_contig_circular_nulls,
    select_curated_gapseq_reference_for_family,
    summarize_candidate_selection_outcome,
    verify_member_locus_strict,
)


def _revcomp(seq: str) -> str:
    tr = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(tr)[::-1]


def test_unknown_or_missing_sequence_and_noninteger_coords_fail_promotion() -> None:
    """Gate 1: Unknown/missing sequence flags, internal stops, or non-integer coords fail promotion."""
    assert is_exact_integer_value(10) is True
    assert is_exact_integer_value("10") is True
    assert is_exact_integer_value(10.0) is True
    assert is_exact_integer_value(10.5) is False
    assert is_exact_integer_value("10.5") is False
    assert is_exact_integer_value(None) is False
    assert is_exact_integer_value(np.nan) is False

    valid_aa = "M" + "A" * 109  # 110 aa
    valid_md5 = hashlib.md5(valid_aa.encode("ascii")).hexdigest()

    df = pd.DataFrame(
        [
            # Row 0: valid interior gene with matching FASTA AA sequence
            {
                "genome": "G1",
                "gene": "G1_c1_1",
                "contig": "c1",
                "start": 10,
                "end": 100,
                "strand": 1,
                "partial": "00",
                "length": 30,
                "internal_stop": False,
                "md5": "0" * 32,
                "contig_in_catalogue": True,
            },
            {
                "genome": "G1",
                "gene": "G1_c1_2",
                "contig": "c1",
                "start": 200,
                "end": 400,
                "strand": 1,
                "partial": "00",
                "length": 66,
                "internal_stop": False,
                "md5": "0" * 32,
                "contig_in_catalogue": True,
            },
            {
                "genome": "G1",
                "gene": "G1_c1_3_valid_interior",
                "contig": "c1",
                "start": 501,
                "end": 833,  # 333 nt = 110 aa + stop
                "strand": 1,
                "partial": "00",
                "length": 110,
                "internal_stop": False,
                "md5": valid_md5,
                "contig_in_catalogue": True,
            },
            # Row 3: interior gene missing from FASTA (must NOT be hardcoded has_full_seq=True)
            {
                "genome": "G1",
                "gene": "G1_c1_4_missing_seq",
                "contig": "c1",
                "start": 901,
                "end": 1233,
                "strand": 1,
                "partial": "00",
                "length": 110,
                "internal_stop": False,
                "md5": valid_md5,
                "contig_in_catalogue": True,
            },
            # Row 4: interior gene with internal '*' stop in FASTA
            {
                "genome": "G1",
                "gene": "G1_c1_5_internal_stop",
                "contig": "c1",
                "start": 1301,
                "end": 1633,
                "strand": 1,
                "partial": "00",
                "length": 110,
                "internal_stop": False,
                "md5": valid_md5,
                "contig_in_catalogue": True,
            },
            # Row 5: interior gene with non-integer start coordinate 1700.5
            {
                "genome": "G1",
                "gene": "G1_c1_6_nonint_coord",
                "contig": "c1",
                "start": 1700.5,
                "end": 2032,
                "strand": 1,
                "partial": "00",
                "length": 110,
                "internal_stop": False,
                "md5": valid_md5,
                "contig_in_catalogue": True,
            },
            {
                "genome": "G1",
                "gene": "G1_c1_7",
                "contig": "c1",
                "start": 2100,
                "end": 2250,
                "strand": 1,
                "partial": "00",
                "length": 50,
                "internal_stop": False,
                "md5": "0" * 32,
                "contig_in_catalogue": True,
            },
            {
                "genome": "G1",
                "gene": "G1_c1_8",
                "contig": "c1",
                "start": 2300,
                "end": 2450,
                "strand": 1,
                "partial": "00",
                "length": 50,
                "internal_stop": False,
                "md5": "0" * 32,
                "contig_in_catalogue": True,
            },
        ]
    )

    fasta_map = {
        "G1_c1_3_valid_interior": valid_aa + "*",
        "G1_c1_5_internal_stop": "MA*EK" + "A" * 105 + "*",
        "G1_c1_6_nonint_coord": valid_aa + "*",
    }
    audited = audit_gene_table_syntax_and_fasta(df, fasta_map).set_index("gene")

    assert bool(audited.loc["G1_c1_3_valid_interior", "has_full_seq"]) is True
    assert bool(audited.loc["G1_c1_3_valid_interior", "provisional_syntax_promotion_flag"]) is True

    # Missing sequence must fail promotion and record unknown_missing_aa_sequence
    assert bool(audited.loc["G1_c1_4_missing_seq", "has_full_seq"]) is False
    assert audited.loc["G1_c1_4_missing_seq", "fasta_sequence_audit_status"] == "unknown_missing_aa_sequence"
    assert bool(audited.loc["G1_c1_4_missing_seq", "provisional_syntax_promotion_flag"]) is False

    # Internal '*' in sequence must never be silently deleted and must fail promotion
    assert bool(audited.loc["G1_c1_5_internal_stop", "has_full_seq"]) is False
    assert audited.loc["G1_c1_5_internal_stop", "fasta_sequence_audit_status"] == "invalid_internal_stop_in_candidate_aa"
    assert bool(audited.loc["G1_c1_5_internal_stop", "provisional_syntax_promotion_flag"]) is False

    # Non-integer coordinate must fail coordinate_syntax_valid and promotion
    assert bool(audited.loc["G1_c1_6_nonint_coord", "coordinate_syntax_valid"]) is False
    assert bool(audited.loc["G1_c1_6_nonint_coord", "provisional_syntax_promotion_flag"]) is False


def test_verify_member_locus_strict_against_refined_contig_and_manifest() -> None:
    """Gate 1b: Member locus promotion requires exact refined-contig CDS translation and verified genome manifest."""
    cds_plus = "ATG" + ("GCT" * 109) + "TAA"  # M + 109 A + stop = 110 AA
    expected_aa = "M" + ("A" * 109)
    contig_seq = ("A" * 100) + cds_plus + ("C" * 100)
    contigs = {"c1": contig_seq}

    # 1. Verified exact complete interior member
    ok = verify_member_locus_strict(
        "verified_exact_match_genome_manifest",
        contigs,
        "c1",
        101,
        100 + len(cds_plus),
        1,
        "00",
        expected_aa + "*",
        "interior",
    )
    assert ok["member_final_promotion_gate_pass"] is True
    assert ok["verification_status"] == "verified_exact_complete_cds"

    # 2. Unverified genome manifest fails promotion even if CDS matches
    bad_man = verify_member_locus_strict(
        "sha256_or_count_mismatch_genome_manifest",
        contigs,
        "c1",
        101,
        100 + len(cds_plus),
        1,
        "00",
        expected_aa,
        "interior",
    )
    assert bad_man["member_final_promotion_gate_pass"] is False

    # 3. Missing candidate AA sequence fails promotion
    missing_aa = verify_member_locus_strict(
        "verified_exact_match_genome_manifest",
        contigs,
        "c1",
        101,
        100 + len(cds_plus),
        1,
        "00",
        None,
        "interior",
    )
    assert missing_aa["member_final_promotion_gate_pass"] is False
    assert missing_aa["verification_status"] == "unknown_missing_candidate_sequence"

    # 4. Inadmissible start codon (GGG with partial='00') fails promotion
    cds_bad_start = "GGG" + ("GCT" * 109) + "TAA"
    bad_start = verify_member_locus_strict(
        "verified_exact_match_genome_manifest",
        {"c1": ("A" * 100) + cds_bad_start + ("C" * 100)},
        "c1",
        101,
        100 + len(cds_bad_start),
        1,
        "00",
        "G" + ("A" * 109),
        "interior",
    )
    assert bad_start["has_admissible_start_codon"] is False
    assert bad_start["member_final_promotion_gate_pass"] is False
    assert bad_start["verification_status"] == "inadmissible_start_codon"


def test_nondark_family_not_accepted_as_positive_control_without_curated_evidence() -> None:
    """Gate 2: Non-dark family is NOT automatically an annotated positive control without curated reference evidence."""
    # 1. Prior heat-stress dark candidate examples cannot be positive controls
    for fid in PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES:
        res = classify_family_control_or_candidate_role(
            fid,
            is_dark=True,
            selection_route="prior_heat_stress_candidate_example",
        )
        assert res["candidate_role"] == "prior_heat_stress_candidate_example"
        with pytest.raises(ValueError, match="prior heat-stress dark candidate example"):
            classify_family_control_or_candidate_role(
                fid,
                is_dark=True,
                selection_route="prior_heat_stress_candidate_example",
                claim_annotated_positive_control=True,
            )

    # 2. Non-dark family without curated reference evidence is downgraded to sequence_annotated_comparator
    unverified_ctrl = classify_family_control_or_candidate_role(
        "F23BRr2::F23BRr2_291009_8",
        is_dark=False,
        selection_route="route_A1_non_dark_comparator",
        claim_annotated_positive_control=True,
        curated_reference_accession=None,
        curated_function_description="putative enzyme",
        curated_evidence_source=None,
        query_domain_signatures=["IPR000001"],
        reference_domain_signatures=None,
        independently_verified_curated_homology=False,
        query_to_curated_ref_cov=0.45,
        strict_reject_unverified_positive_control=False,
    )
    assert unverified_ctrl["candidate_role"] == "sequence_annotated_comparator"
    assert unverified_ctrl["positive_control_evidence_status"] == "positive_control_evidence_unavailable"
    assert unverified_ctrl["has_curated_positive_control_evidence"] is False

    # 2b. Negative fixture (stage3-final-integrity-review.txt Item 2):
    # Unrelated proteins forced-aligned at full coverage (1.0) plus an unrelated InterPro domain
    # on the query and a reaction SwissProt reference MUST NOT become an annotated positive control!
    forced_msa_unrelated = classify_family_control_or_candidate_role(
        "F23BRr2::F23BRr2_291009_8",
        is_dark=False,
        selection_route="route_A1_non_dark_comparator",
        claim_annotated_positive_control=True,
        curated_reference_accession="P58300",
        curated_function_description="Phosphoribosylformylglycinamidine synthase subunit PurQ",
        curated_evidence_source="SwissProt_Archaea_rev",
        query_domain_signatures=["Pfam:PF99999(Unrelated_domain)"],
        reference_domain_signatures=["Pfam:PF00117(GATase)"],
        independently_verified_curated_homology=False,
        query_to_curated_ref_cov=1.0,
        strict_reject_unverified_positive_control=False,
    )
    assert forced_msa_unrelated["candidate_role"] == "sequence_annotated_comparator"
    assert forced_msa_unrelated["positive_control_evidence_status"] == "positive_control_evidence_unavailable"
    assert forced_msa_unrelated["has_curated_positive_control_evidence"] is False

    # And raises ValueError when strict_reject_unverified_positive_control=True
    with pytest.raises(ValueError, match="Cannot promote non-dark family"):
        classify_family_control_or_candidate_role(
            "F23BRr2::F23BRr2_291009_8",
            is_dark=False,
            selection_route="route_A1_non_dark_comparator",
            claim_annotated_positive_control=True,
            curated_reference_accession="P58300",
            curated_function_description="PurQ",
            curated_evidence_source="SwissProt_Archaea_rev",
            query_domain_signatures=["Pfam:PF99999(Unrelated_domain)"],
            reference_domain_signatures=None,
            independently_verified_curated_homology=False,
            query_to_curated_ref_cov=1.0,
            strict_reject_unverified_positive_control=True,
        )

    # 3. Non-dark family WITH independently verified curated homology + matching reference domain architecture + coverage >= 0.70 passes
    verified_ctrl = classify_family_control_or_candidate_role(
        "F23BRr2::F23BRr2_291009_8",
        is_dark=False,
        selection_route="route_A1_curated_control",
        claim_annotated_positive_control=True,
        curated_reference_accession="P58300",
        curated_function_description="Phosphoribosylformylglycinamidine synthase subunit PurQ",
        curated_evidence_source="SwissProt_Archaea_rev",
        query_domain_signatures=["Pfam:PF00117(GATase)"],
        reference_domain_signatures=["Pfam:PF00117(GATase)"],
        independently_verified_curated_homology=True,
        query_to_curated_ref_cov=0.92,
        strict_reject_unverified_positive_control=True,
    )
    assert verified_ctrl["candidate_role"] == "annotated_positive_control"
    assert verified_ctrl["has_curated_positive_control_evidence"] is True


def test_mixed_locus_maxima_rejected_and_rep_own_separated() -> None:
    """Gate 3: Mixed-locus maxima are rejected; invalid/absent rep sequence never becomes synthetic valid or qualified."""
    rep_gene = "F24Dr2_413137_12"
    rep_seq = "M" + "A" * 150
    other_gene = "O24Dr3_99999_4"
    other_seq = "M" + "G" * 180

    # Negative fixture (stage3-final-integrity-review.txt Item 1):
    # Absent (None/"") or invalid ("MA*EK...") representative sequence must NEVER become a synthetic valid sequence
    # or produce a qualified structural/domain result even if pLDDT=0.95 and qTM=0.85 are supplied!
    for bad_seq in (None, "", "MA*EK" + "A" * 150):
        bad_res = qualify_structural_and_domain_evidence_strict(
            rep_gene,
            bad_seq,
            [
                {
                    "locus_gene": rep_gene,
                    "locus_aa_seq": bad_seq,
                    "claim_as_representative_own": True,
                    "plddt_raw": 0.95,
                    "plddt_source_scale": "0_to_1",
                    "pdb_target": "4l3p_A",
                    "pdb_qtmscore": 0.85,
                    "pdb_qcov": 0.92,
                    "hhsearch_pfam": "PF00001",
                    "hhsearch_prob": 99.0,
                    "hhsearch_qcov": 0.90,
                    "interpro_domains": ["Pfam:PF00001"],
                }
            ],
        )
        assert bad_res["rep_sequence_audit_status"] != "ok"
        assert bad_res["rep_aa_length"] == 0
        assert bad_res["rep_own_model_qualified"] is False
        assert bad_res["rep_own_qualified_pdb_fold_match"] is False
        assert bad_res["rep_own_has_qualified_whole_fold_support"] is False
        assert bad_res["rep_own_has_domain_or_profile_support"] is False

    # Attempting to claim another member's structural/domain record as representative-own raises ValueError
    with pytest.raises(ValueError, match="Mixed-locus evidence aggregation rejected"):
        qualify_structural_and_domain_evidence_strict(
            rep_gene,
            rep_seq,
            [
                {
                    "locus_gene": other_gene,
                    "locus_aa_seq": other_seq,
                    "claim_as_representative_own": True,
                    "plddt_raw": 0.88,
                    "plddt_source_scale": "0_to_1",
                    "pdb_target": "4l3p_A",
                    "pdb_qtmscore": 0.72,
                    "pdb_qcov": 0.85,
                }
            ],
        )

    # When properly passed without false claim, rep_own_* stays empty/false while other_member_* records the member
    res = qualify_structural_and_domain_evidence_strict(
        rep_gene,
        rep_seq,
        [
            {
                "locus_gene": rep_gene,
                "locus_aa_seq": rep_seq,
                "claim_as_representative_own": True,
                "model_source": "ColabFold_AF2",
                "model_polymer_aa_seq": rep_seq,
                "plddt_raw": 0.62,  # below 0.70
                "plddt_source_scale": "0_to_1",
                "pdb_hit_source": "ColabFold_AF2",
                "pdb_target": "1abc_A",
                "pdb_qtmscore": 0.55,
                "pdb_qcov": 0.80,
                "afdb_hit_source": "ColabFold_AF2",
                "afdb_target": "AF-P12345-F1",
                "afdb_evalue": 1e-12,
                "afdb_qcov": 0.85,
                "afdb_qtmscore": np.nan,  # AFDB has no qtmscore in c9 table
                "interpro_domains": [],
                "interpro_non_domain_flags": ["MobiDBLite:mobidb-lite"],
            },
            {
                "locus_gene": other_gene,
                "locus_aa_seq": other_seq,
                "claim_as_representative_own": False,
                "model_source": "ESMFold",
                "model_polymer_aa_seq": other_seq,
                "plddt_raw": 0.88,
                "plddt_source_scale": "0_to_1",
                "pdb_hit_source": "ESMFold",
                "pdb_target": "4l3p_A",
                "pdb_qtmscore": 0.72,
                "pdb_qcov": 0.85,
                "hhsearch_pfam": "PF00001",
                "hhsearch_prob": 95.0,
                "hhsearch_qcov": 0.75,
            },
        ],
    )
    assert res["rep_own_model_qualified"] is False
    assert res["rep_own_qualified_pdb_fold_match"] is False
    assert res["rep_own_qualified_afdb_fold_match"] is False
    assert res["rep_own_afdb_support_status"] == "afdb_3di_aa_alignment_only_qtm_unavailable"
    assert res["rep_own_has_qualified_whole_fold_support"] is False
    assert res["rep_own_has_domain_or_profile_support"] is False
    assert res["rep_own_interpro_domain_signatures"] == "none"
    assert "MobiDBLite" in res["rep_own_interpro_nondomain_disorder_coils"]
    assert other_gene in res["other_member_qualified_structural_evidence"]
    assert other_gene in res["other_member_qualified_profile_evidence"]


def test_adversarial_model_polymer_sequence_and_source_mismatch_rejected() -> None:
    """Gate 3b (stage3-final-integrity-review.txt Item 3): Same protein ID + same length but wrong model polymer sequence or mismatched model/hit source fails qualification."""
    rep_gene = "F24Dr2_413137_12"
    rep_seq = "M" + "A" * 150  # 151 AA
    wrong_poly_same_len = "M" + "G" * 150  # same length (151 AA), wrong sequence!

    # 1. Same protein ID and exact same length (151 AA), but model polymer sequence != audited locus sequence
    wrong_poly_res = qualify_structural_and_domain_evidence_strict(
        rep_gene,
        rep_seq,
        [
            {
                "locus_gene": rep_gene,
                "locus_aa_seq": rep_seq,
                "claim_as_representative_own": True,
                "model_source": "ColabFold_AF2",
                "model_polymer_aa_seq": wrong_poly_same_len,
                "plddt_raw": 0.91,
                "plddt_source_scale": "0_to_1",
                "pdb_hit_source": "ColabFold_AF2",
                "pdb_target": "4l3p_A",
                "pdb_qtmscore": 0.78,
                "pdb_qcov": 0.88,
            }
        ],
    )
    assert wrong_poly_res["rep_own_model_seq_provenance_status"] == "mismatch_model_polymer_sequence_same_length"
    assert wrong_poly_res["rep_own_model_qualified"] is False
    assert wrong_poly_res["rep_own_qualified_pdb_fold_match"] is False
    assert wrong_poly_res["rep_own_has_qualified_whole_fold_support"] is False
    # Raw hit remains visible as unverified candidate evidence
    assert wrong_poly_res["rep_own_pdb_target"] == "4l3p_A"

    # 2. Missing model polymer sequence (length equality alone cannot qualify!)
    missing_poly_res = qualify_structural_and_domain_evidence_strict(
        rep_gene,
        rep_seq,
        [
            {
                "locus_gene": rep_gene,
                "locus_aa_seq": rep_seq,
                "claim_as_representative_own": True,
                "model_source": "ColabFold_AF2",
                "model_polymer_aa_seq": None,
                "plddt_raw": 0.91,
                "plddt_source_scale": "0_to_1",
                "pdb_hit_source": "ColabFold_AF2",
                "pdb_target": "4l3p_A",
                "pdb_qtmscore": 0.78,
                "pdb_qcov": 0.88,
            }
        ],
    )
    assert missing_poly_res["rep_own_model_seq_provenance_status"] == "unverified_missing_model_polymer_sequence"
    assert missing_poly_res["rep_own_model_qualified"] is False
    assert missing_poly_res["rep_own_qualified_pdb_fold_match"] is False

    # 3. Matching model polymer sequence, but AF2 model pLDDT paired with an ESMFold Foldseek hit (cross-source pairing)
    cross_source_res = qualify_structural_and_domain_evidence_strict(
        rep_gene,
        rep_seq,
        [
            {
                "locus_gene": rep_gene,
                "locus_aa_seq": rep_seq,
                "claim_as_representative_own": True,
                "model_source": "ColabFold_AF2",
                "model_polymer_aa_seq": rep_seq,
                "plddt_raw": 0.91,
                "plddt_source_scale": "0_to_1",
                "pdb_hit_source": "ESMFold",
                "pdb_target": "4l3p_A",
                "pdb_qtmscore": 0.78,
                "pdb_qcov": 0.88,
            }
        ],
    )
    assert cross_source_res["rep_own_model_qualified"] is True
    assert cross_source_res["rep_own_qualified_pdb_fold_match"] is False
    assert cross_source_res["rep_own_has_qualified_whole_fold_support"] is False
    assert "mismatched_model_vs_hit_source" in cross_source_res["rep_own_model_seq_provenance_status"]


def test_plddt_explicit_scale_and_coverage_qualification() -> None:
    """Gate 4: Explicit per-source pLDDT scale (0, 1, 1.4, 70, 100) and qcov >= 0.50 qualification."""
    assert scale_plddt_explicit_source(0.0, "0_to_1") == pytest.approx(0.0)
    assert scale_plddt_explicit_source(1.0, "0_to_1") == pytest.approx(1.0)
    with pytest.raises(ValueError, match="out of bounds"):
        scale_plddt_explicit_source(1.4, "0_to_1")
    with pytest.raises(ValueError, match="out of bounds"):
        scale_plddt_explicit_source(70.0, "0_to_1")

    assert scale_plddt_explicit_source(70.0, "0_to_100") == pytest.approx(0.70)
    assert scale_plddt_explicit_source(100.0, "0_to_100") == pytest.approx(1.0)
    assert scale_plddt_explicit_source(1.4, "0_to_100") == pytest.approx(0.014)

    # High qTM-score (0.75) on a 350 AA protein with pLDDT 0.85 but low query coverage (qcov = 0.25 < 0.50)
    # must NOT qualify as whole-fold structural support
    rep_seq = "M" + "A" * 349
    low_cov = qualify_structural_and_domain_evidence_strict(
        "G_rep",
        rep_seq,
        [
            {
                "locus_gene": "G_rep",
                "locus_aa_seq": rep_seq,
                "claim_as_representative_own": True,
                "model_source": "ColabFold_AF2",
                "model_polymer_aa_seq": rep_seq,
                "plddt_raw": 85.0,
                "plddt_source_scale": "0_to_100",
                "pdb_hit_source": "ColabFold_AF2",
                "pdb_target": "2xyz_A",
                "pdb_qtmscore": 0.75,
                "pdb_qcov": 0.25,
            }
        ],
    )
    assert low_cov["rep_own_model_qualified"] is True
    assert low_cov["rep_own_qualified_pdb_fold_match"] is False


def test_zero_eligible_families_legitimate_vs_missing_input_failure() -> None:
    """Gate 4b (stage3-negative-outcome-review.txt): Zero eligible dark families after strict audit is valid; missing input fails."""
    with pytest.raises(ValueError, match="Missing or empty"):
        summarize_candidate_selection_outcome(None, None)

    r1_zero = pd.DataFrame(
        [
            {
                "family_id": "F_dark_1",
                "inherited_catalogue_dark": True,
                "final_candidate_promotion_eligible": False,
                "n_verified_complete_interior_member_loci": 1,
            }
        ]
    )
    r2_zero = pd.DataFrame(
        [
            {
                "family_id": "F_dark_2",
                "final_candidate_promotion_eligible": False,
                "n_verified_complete_interior_member_loci": 2,
            }
        ]
    )
    res_zero = summarize_candidate_selection_outcome(r1_zero, r2_zero)
    assert res_zero["n_total_selected_new_dark"] == 0
    assert res_zero["no_eligible_family"] is True
    assert res_zero["selection_outcome_status"] == "no_eligible_dark_family_after_strict_audit"


def test_leave_order_out_threshold_calibration_no_arbitrary_fallback() -> None:
    """Gate 5: Leave-order-out threshold calibration never falls back to 2.0 when target 0.50 is not reached."""
    rng = np.random.default_rng(42)
    train_scores = rng.uniform(0.5, 5.0, size=100)
    # Max recovery is ~0.20 everywhere, well below 0.50 target
    train_top1 = np.zeros(100, dtype="float64")
    train_top1[:15] = 1.0
    eval_scores = rng.uniform(0.5, 5.0, size=40)
    eval_top1 = np.zeros(40, dtype="float64")

    res = calibrate_leave_order_out_threshold(
        train_scores,
        train_top1,
        eval_scores,
        eval_top1,
        documented_target_top1_recovery=0.50,
        min_train_passing=20,
        is_held_out_split=True,
    )
    assert res["evaluation_mode"] == "leave_order_out_threshold_calibration"
    assert res["is_held_out_evaluation"] is True
    assert res["threshold_calibration_status"] == "no_qualifying_threshold_attained"
    assert pd.isna(res["calibrated_score_threshold"])
    assert pd.isna(res["eval_conditional_top1_recovery_at_threshold"])

    # Pooled in-sample descriptive mode must explicitly report is_held_out_evaluation=False
    pooled = calibrate_leave_order_out_threshold(
        train_scores,
        np.ones(100, dtype="float64"),
        train_scores,
        np.ones(100, dtype="float64"),
        documented_target_top1_recovery=0.50,
        min_train_passing=20,
        is_held_out_split=False,
    )
    assert pooled["evaluation_mode"] == "in_sample_pooled_descriptive_not_held_out"
    assert pooled["is_held_out_evaluation"] is False
    assert pooled["threshold_calibration_status"] == "qualifying_threshold_attained"


def test_motif_conservation_requires_stated_reference_and_seven_axes_do_not_overclaim() -> None:
    """Gate 6 & 7: Motif conservation requires stated reference; >=3 genera never auto-promotes to mechanism."""
    aln = {
        "m1": "MA-H-E-K",
        "m2": "MA-H-E-K",
    }
    no_ref = evaluate_reference_motif_conservation(aln, reference_id=None, documented_catalytic_sites=None)
    assert no_ref["status"] == "curated_reference_or_active_site_residues_absent_functional_specificity_pending"

    with_ref = evaluate_reference_motif_conservation(
        {"ref_sp": "MA-H-E-K", "m1": "MA-H-E-K", "m2": "MA-H-E-K"},
        reference_id="ref_sp",
        documented_catalytic_sites=[(3, "H", "catalytic_His"), (4, "DE", "catalytic_acid")],
    )
    assert with_ref["status"] == "conserved_at_all_2_documented_reference_sites"
    assert with_ref["n_sites_conserved_in_all_members"] == 2

    metrics = compute_alignment_summary_metrics({"m1": "MA-HEK", "m2": "MA-HDK"})
    assert metrics["n_aligned_sequences"] == 2
    assert metrics["alignment_cols"] == 6
    assert metrics["median_pairwise_aa_identity"] == pytest.approx(4.0 / 5.0)

    axes = evaluate_seven_evidence_axes(
        candidate_role="de_novo_dark_candidate_family",
        rep_final_promotion_gate_pass=True,
        n_verified_complete_loci=6,
        n_total_loci=6,
        rep_signalp_class="SP",
        rep_tmbed_s=True,
        rep_tm_helices=0,
        rep_tm_beta=0,
        frac_family_exported=1.0,
        rep_has_own_domain_or_profile=False,
        rep_own_domain_summary="none",
        rep_has_own_qualified_whole_fold=False,
        rep_own_fold_summary="none",
        other_member_evidence_summary="none",
        motif_conservation_status=no_ref["status"],
        median_alignment_coverage=0.98,
        n_verified_orders=2,
        n_verified_independent_genus_blocks=4,
        max_copies_per_genome=1,
        specific_hypothesis_name="specific_cofactor_salvage_flank",
        specific_context_obs_adj=4,
        specific_context_p_emp=0.012,
        n_verified_genus_blocks_with_specific_context=4,
        competing_assignments_summary="broad_repeat_fold_vs_unannotated_surface_protein",
        positive_control_evidence_status="not_applicable_dark_candidate",
    )
    assert "localisation_not_catalytic_or_S_layer_proof" in axes["axis2_localisation_status"]
    assert "posthoc_guild_selection_dependent_exploratory_never_establishes_mechanism" in axes["axis6_specific_context_support_status"]
    assert "no_automated_mechanistic_promotion" in axes["supervisor_interpretation_handoff"]


def test_route_a1_real_schema_handoff_to_curated_reference_and_reject_wrong_columns(tmp_path: Path) -> None:
    """Gate 6b (stage3-real-schema-followup.txt): Real Stage 1b proposal columns (`missing_reaction`, `pathway`)
    reach curated reference selection; non-existent `rxn_id`/`pathway_id` schema fails immediately.
    """
    rev_dir = tmp_path / "rev"
    rev_dir.mkdir(parents=True, exist_ok=True)
    (rev_dir / "rxn00123.fasta").write_text(
        ">sp|P58300|PURQ_METJA Phosphoribosylformylglycinamidine synthase subunit PurQ\n"
        "MKKIAVIGAGAVGSTLAFAAAQRGIAREIVLEDIAK\n"
    )

    real_prop_df = pd.DataFrame(
        [
            # Verified complete interior priority proposal -> must reach reaction/pathway context & reference selection
            {
                "genome": "G1",
                "candidate_gene": "G1_c1_3",
                "pathway": "|PWY-5100|",
                "missing_reaction": "rxn00123",
                "reaction_name": "Phosphoribosylformylglycinamidine synthase",
                "ec": "6.3.5.3",
                "review_priority": True,
                "inherited_catalogue_dark": False,
                "family_id": "FAM_PURQ",
                "contig_edge_status": "interior",
                "final_promotion_gate_pass": True,
            },
            # Unverified or edge row -> must NOT pollute verified context
            {
                "genome": "G2",
                "candidate_gene": "G2_c9_1",
                "pathway": "|PWY-UNVERIFIED|",
                "missing_reaction": "rxn99999",
                "reaction_name": "Unverified edge reaction",
                "ec": "1.1.1.1",
                "review_priority": True,
                "inherited_catalogue_dark": True,
                "family_id": "FAM_PURQ",
                "contig_edge_status": "edge",
                "final_promotion_gate_pass": False,
            },
        ]
    )

    assert set(STAGE1B_VERIFIED_PROPOSALS_REQUIRED_COLS).issubset(real_prop_df.columns)
    ctx_map = extract_verified_route_a1_reaction_pathway_context(real_prop_df)
    assert "FAM_PURQ" in ctx_map
    assert ctx_map["FAM_PURQ"]["missing_reactions"] == ["rxn00123"]
    assert ctx_map["FAM_PURQ"]["pathways"] == ["|PWY-5100|"]
    assert ctx_map["FAM_PURQ"]["reaction_names"] == ["Phosphoribosylformylglycinamidine synthase"]
    assert ctx_map["FAM_PURQ"]["ecs"] == ["6.3.5.3"]

    cref = select_curated_gapseq_reference_for_family("FAM_PURQ", ctx_map, tmp_path)
    assert cref["has_reviewed_ref"] is True
    assert cref["rxn_id"] == "rxn00123"
    assert cref["accession"] == "P58300"
    assert cref["sequence"].startswith("MKKIAVIGAG")

    # Negative fixture: wrong schema (`rxn_id`, `pathway_id` instead of `missing_reaction`, `pathway`) must raise ValueError
    wrong_schema_df = real_prop_df.rename(columns={"missing_reaction": "rxn_id", "pathway": "pathway_id"})
    with pytest.raises(ValueError, match="missing required columns.*pathway.*missing_reaction"):
        extract_verified_route_a1_reaction_pathway_context(wrong_schema_df)


def test_grouped_contig_circular_nulls_real_api_and_merge_same_contig_targets() -> None:
    """Gate 6c: Exercise real `within_contig_circular_null(contig_ids, target_mask, flank_target_mask)` API,
    merge multiple target members on the same `(genome, contig)`, and separate specific vs generic nulls.
    """
    spec_guild = "extracellular_peptidase_or_proteolysis"
    gen_guild = "generic_any_known_annotation"

    # Contig G1::c1 has 12 genes and TWO target family members at indices 2 and 8 (same contig!)
    c1_guilds: list[set[str]] = [set() for _ in range(12)]
    c1_guilds[1] = {spec_guild, gen_guild}  # adjacent to target at 2
    c1_guilds[3] = {gen_guild}              # generic-only adjacent to target at 2
    c1_guilds[7] = {spec_guild, gen_guild}  # adjacent to target at 8

    # Contig G2::c5 has 10 genes and ONE target family member at index 4
    c2_guilds: list[set[str]] = [set() for _ in range(10)]
    c2_guilds[5] = {gen_guild}              # generic-only adjacent to target at 4

    records = [
        {"genome": "G1", "contig": "c1", "cand_idx": 2, "contig_guilds": c1_guilds},
        {"genome": "G1", "contig": "c1", "cand_idx": 8, "contig_guilds": c1_guilds},
        {"genome": "G2", "contig": "c5", "cand_idx": 4, "contig_guilds": c2_guilds},
    ]

    bundle = run_grouped_contig_circular_nulls(
        records,
        specific_guild=spec_guild,
        specific_guild_selection_mode="exploratory_posthoc_observed_neighbourhood_winner_selection_dependent",
        window=2,
        n_perm=200,
        seed=20261108,
    )
    assert bundle["n_member_records"] == 3
    assert bundle["n_unique_contigs"] == 2
    assert bundle["n_merged_same_contig_target_loci"] == 1
    assert bundle["specific_null_inference_scope"] == "exploratory_posthoc_selection_dependent_not_independent_confirmation"

    spec_res = bundle["specific_null"]
    gen_res = bundle["generic_null"]
    assert spec_res["n_target_loci"] == 3
    assert gen_res["n_target_loci"] == 3
    # Target at G1::c1[2] and G1::c1[8] have spec_guild within window=2; G2::c5[4] has only generic_guild
    assert spec_res["obs_adjacent"] == 2
    assert gen_res["obs_adjacent"] == 3
    assert 0.0 < float(spec_res["p_empirical_one_sided"]) <= 1.0
    assert 0.0 < float(gen_res["p_empirical_one_sided"]) <= 1.0


def test_gspa_dotted_columns_preserved_and_missing_schema_rejected() -> None:
    """Gate 7c (stage3-gspa-summary-schema-review.txt & stage3-final-summary-review.txt):
    Verify dotted GSPA columns (`nonroot_0.5_proteins`, `nonroot_0.5_dark`, etc.) survive without
    falling back to 0 via `itertuples`/`getattr`, preserve `context` vs `shuffled` separation,
    and use the non-inferential `scope_and_limitation_note`.
    """
    g1 = "F23PRr1__single_SemiBin_76"
    thermo_meta = pd.DataFrame(
        [
            {
                "genome": g1,
                "tax_order": "o__JACQPN01",
                "genus_singleton_block": "g__JAJPHT01",
                "checkm_completeness": 92.5,
            }
        ]
    )
    alias_map = build_manifest_genome_alias_resolver([g1])

    gspa_sum_df = pd.DataFrame(
        [
            {
                "genome": g1,
                "mode": "fast",
                "proteins": 3427,
                "dark": 1062,
                "strict_dark": 920,
                "claims": 904316,
                "nonroot_0.5_proteins": 1350,
                "nonroot_0.5_dark": 0,
                "nonroot_0.5_strict_dark": 0,
                "mf_bp_depth4_0.5_proteins": 1315,
                "mf_bp_depth4_0.5_dark": 0,
                "mf_bp_depth4_0.5_strict_dark": 0,
                "nonroot_0.8_proteins": 1348,
                "nonroot_0.8_dark": 0,
                "nonroot_0.8_strict_dark": 0,
                "mf_bp_depth4_0.8_proteins": 1146,
                "mf_bp_depth4_0.8_dark": 0,
                "mf_bp_depth4_0.8_strict_dark": 0,
            },
            {
                "genome": g1,
                "mode": "cnn",
                "proteins": 3427,
                "dark": 1062,
                "strict_dark": 920,
                "claims": 7614290,
                "nonroot_0.5_proteins": 3425,
                "nonroot_0.5_dark": 1060,
                "nonroot_0.5_strict_dark": 918,
                "mf_bp_depth4_0.5_proteins": 3347,
                "mf_bp_depth4_0.5_dark": 1014,
                "mf_bp_depth4_0.5_strict_dark": 878,
                "nonroot_0.8_proteins": 1776,
                "nonroot_0.8_dark": 172,
                "nonroot_0.8_strict_dark": 154,
                "mf_bp_depth4_0.8_proteins": 1643,
                "mf_bp_depth4_0.8_dark": 107,
                "mf_bp_depth4_0.8_strict_dark": 98,
            },
            {
                "genome": g1,
                "mode": "cnn_baseline",
                "proteins": 3427,
                "dark": 1062,
                "strict_dark": 920,
                "claims": 102489,
                "nonroot_0.5_proteins": 1768,
                "nonroot_0.5_dark": 169,
                "nonroot_0.5_strict_dark": 151,
                "mf_bp_depth4_0.5_proteins": 1602,
                "mf_bp_depth4_0.5_dark": 100,
                "mf_bp_depth4_0.5_strict_dark": 92,
                "nonroot_0.8_proteins": 0,
                "nonroot_0.8_dark": 0,
                "nonroot_0.8_strict_dark": 0,
                "mf_bp_depth4_0.8_proteins": 0,
                "mf_bp_depth4_0.8_dark": 0,
                "mf_bp_depth4_0.8_strict_dark": 0,
            },
            {
                "genome": g1,
                "mode": "context",
                "proteins": 3427,
                "dark": 1062,
                "strict_dark": 920,
                "claims": 38976,
                "nonroot_0.5_proteins": 1348,
                "nonroot_0.5_dark": 0,
                "nonroot_0.5_strict_dark": 0,
                "mf_bp_depth4_0.5_proteins": 1123,
                "mf_bp_depth4_0.5_dark": 0,
                "mf_bp_depth4_0.5_strict_dark": 0,
                "nonroot_0.8_proteins": 0,
                "nonroot_0.8_dark": 0,
                "nonroot_0.8_strict_dark": 0,
                "mf_bp_depth4_0.8_proteins": 0,
                "mf_bp_depth4_0.8_dark": 0,
                "mf_bp_depth4_0.8_strict_dark": 0,
            },
            {
                "genome": g1,
                "mode": "shuffled",
                "proteins": 3427,
                "dark": 1062,
                "strict_dark": 920,
                "claims": 38976,
                "nonroot_0.5_proteins": 1348,
                "nonroot_0.5_dark": 0,
                "nonroot_0.5_strict_dark": 0,
                "mf_bp_depth4_0.5_proteins": 1123,
                "mf_bp_depth4_0.5_dark": 0,
                "mf_bp_depth4_0.5_strict_dark": 0,
                "nonroot_0.8_proteins": 0,
                "nonroot_0.8_dark": 0,
                "nonroot_0.8_strict_dark": 0,
                "mf_bp_depth4_0.8_proteins": 0,
                "mf_bp_depth4_0.8_dark": 0,
                "mf_bp_depth4_0.8_strict_dark": 0,
            },
        ]
    )
    gspa_ctx_df = pd.DataFrame(
        [
            {"genome": g1, "variant": "fast", "mode": "context", "increased_claims": 5, "increased_proteins": 5, "crossed_05": 5, "dark_crossed_05": 0, "max_increase": 0.0928},
            {"genome": g1, "variant": "fast", "mode": "shuffled", "increased_claims": 2, "increased_proteins": 2, "crossed_05": 2, "dark_crossed_05": 0, "max_increase": 0.1134},
            {"genome": g1, "variant": "cnn", "mode": "context", "increased_claims": 1, "increased_proteins": 1, "crossed_05": 1, "dark_crossed_05": 0, "max_increase": 0.0821},
            {"genome": g1, "variant": "cnn", "mode": "shuffled", "increased_claims": 0, "increased_proteins": 0, "crossed_05": 0, "dark_crossed_05": 0, "max_increase": 0.0},
        ]
    )

    out = build_gspa_pilot_3genome_audit_table(gspa_sum_df, gspa_ctx_df, thermo_meta, alias_map)
    by_mode = out.set_index("gspa_mode")
    assert int(by_mode.loc["fast", "gspa_nonroot_05_proteins"]) == 1350
    assert int(by_mode.loc["fast", "gspa_nonroot_05_dark_proteins"]) == 0
    assert int(by_mode.loc["cnn", "gspa_nonroot_05_proteins"]) == 3425
    assert int(by_mode.loc["cnn", "gspa_nonroot_05_dark_proteins"]) == 1060
    assert int(by_mode.loc["cnn_baseline", "gspa_nonroot_05_proteins"]) == 1768
    assert int(by_mode.loc["cnn_baseline", "gspa_nonroot_05_dark_proteins"]) == 169
    assert pd.isna(by_mode.loc["fast", "mode_crossed_05"])
    assert int(by_mode.loc["context", "mode_crossed_05"]) == 5
    assert int(by_mode.loc["shuffled", "mode_crossed_05"]) == 2
    assert (
        by_mode.loc["cnn", "scope_and_limitation_note"]
        == "three-genome exploratory pilot; observed and shuffled counts reported separately; prediction counts do not establish functional precision or statistical separation"
    )

    # Missing dotted column or NaN must raise ValueError (never fall back to 0)
    bad_missing = gspa_sum_df.drop(columns=["nonroot_0.5_dark"])
    with pytest.raises(ValueError, match="nonroot_0.5_dark"):
        build_gspa_pilot_3genome_audit_table(bad_missing, gspa_ctx_df, thermo_meta, alias_map)

    bad_nan = gspa_sum_df.copy()
    bad_nan.loc[1, "nonroot_0.5_dark"] = np.nan
    with pytest.raises(ValueError, match="nonroot_0.5_dark"):
        build_gspa_pilot_3genome_audit_table(bad_nan, gspa_ctx_df, thermo_meta, alias_map)


def test_c12_c18_manifest_strict_id_and_reject_unregistered_unsuffixed_parent() -> None:
    """Gate 7d (stage3-final-summary-review.txt Item 1 & stage3-gspa-summary-schema-review.txt):
    - Do not create unsuffixed parent alias merely because canonical mag ends `_sub` or `refined_sub == 1`.
    - Reject unregistered unsuffixed parent even when `refined_sub == 1` exists.
    - Enforce 1:1 medium selection (do not count multiple media as distinct genomes) and keep missing models explicit.
    """
    g_sub = "15Dr1__self_SemiBin_98_sub"
    g_parent_unregistered = "15Dr1__self_SemiBin_98"
    g_missing_model = "S22PRr1__self_SemiBin_29"

    thermo_meta = pd.DataFrame(
        [
            {"genome": g_sub, "tax_order": "o__JACQPN01", "genus_singleton_block": "g__JAIVGT01", "checkm_completeness": 60.25},
            {"genome": g_missing_model, "tax_order": "o__JACPAO01", "genus_singleton_block": "g__JAHFTW01", "checkm_completeness": 72.40},
        ]
    )
    manifest_df = pd.DataFrame(
        [
            {"mag": g_sub, "path": f"genomes/{g_sub}.fa", "bytes": 100, "sha256": "a" * 64, "contigs": 10, "bases": 1000},
            {"mag": g_missing_model, "path": f"genomes/{g_missing_model}.fa", "bytes": 100, "sha256": "b" * 64, "contigs": 10, "bases": 1000},
        ]
    )
    canonical_rep_df = pd.DataFrame(
        [
            {"mag": g_sub, "drep_fasta_basename": f"{g_sub}.fa", "refined_sub": 1},
            {"mag": g_missing_model, "drep_fasta_basename": f"{g_missing_model}.fa", "refined_sub": 0},
        ]
    )

    alias_map = build_manifest_genome_alias_resolver(
        [g_sub, g_missing_model],
        manifest_df=manifest_df,
        canonical_rep_df=canonical_rep_df,
    )
    assert alias_map[g_sub] == g_sub
    assert g_parent_unregistered not in alias_map, "Must not infer unsuffixed parent alias from _sub or refined_sub=1"

    mod_df = pd.DataFrame(
        [
            {"genome": g_sub, "eq_thermoplasmatota": True, "pathways": 415, "predicted_pathways": 31, "good_blast_reactions": 183, "transporter_hits": 4729},
            {"genome": g_missing_model, "eq_thermoplasmatota": True, "pathways": 415, "predicted_pathways": 30, "good_blast_reactions": 168, "transporter_hits": 4300},
        ]
    )
    c12_df = pd.DataFrame(
        [
            {"genome": g_sub, "medium": "native", "filled": "TRUE", "n_rxn": 870, "n_draft": 565, "origin_counts": "0:435", "n_genes": 109, "growth": 0.015, "fba_status": 0},
            {"genome": g_sub, "medium": "mmglu", "filled": "FALSE", "n_rxn": np.nan, "n_draft": np.nan, "origin_counts": "", "n_genes": np.nan, "growth": np.nan, "fba_status": np.nan},
            {"genome": g_sub, "medium": "desert", "filled": "TRUE", "n_rxn": 862, "n_draft": 565, "origin_counts": "0:435", "n_genes": 109, "growth": 0.043, "fba_status": 0},
            {"genome": "GCA_021163455.1", "medium": "native", "filled": "TRUE", "n_rxn": 900, "n_draft": 600, "origin_counts": "0:500", "n_genes": 120, "growth": 0.02, "fba_status": 0},
        ]
    )
    c18_df = pd.DataFrame(
        [
            {"genome": g_sub, "pathway": "ARGASEDEG-PWY", "has_hole": 0, "n_window_genes": 1, "n_dark": 1, "dark_families": "F1"},
            {"genome": g_sub, "pathway": "TCA", "has_hole": 1, "n_window_genes": 4, "n_dark": 2, "dark_families": "F2"},
        ]
    )
    prop_df = pd.DataFrame(
        [
            {
                "genome": g_sub,
                "candidate_gene": "G1_1",
                "pathway": "TCA",
                "missing_reaction": "RXN1",
                "reaction_name": "rxn1",
                "ec": "1.1.1.1",
                "review_priority": True,
                "inherited_catalogue_dark": True,
                "family_id": "F2",
                "contig_edge_status": "interior",
                "final_promotion_gate_pass": True,
            }
        ]
    )

    order_df = reconcile_and_summarize_order_metabolic_tables(
        thermo_meta, mod_df, c12_df, c18_df, prop_df, alias_map, primary_medium="native"
    )
    assert int(order_df["n_genomes"].sum()) == 2
    by_ord = order_df.set_index("tax_order")
    assert int(by_ord.loc["o__JACQPN01", "n_genomes"]) == 1
    assert int(by_ord.loc["o__JACQPN01", "n_genomes_with_gapseq_model"]) == 1
    assert float(by_ord.loc["o__JACQPN01", "median_gapfilled_rxn_added"]) == pytest.approx(870.0 - 565.0)
    assert float(by_ord.loc["o__JACQPN01", "median_desert_gapfilled_rxn_added"]) == pytest.approx(862.0 - 565.0)
    assert int(by_ord.loc["o__JACPAO01", "n_genomes_missing_gapseq_model"]) == 1
    assert by_ord.loc["o__JACPAO01", "missing_gapseq_model_genomes"] == g_missing_model
    assert pd.isna(by_ord.loc["o__JACPAO01", "median_total_model_rxn"])

    # Unregistered unsuffixed parent ID must fail even when refined_sub=1 exists in canonical_rep_df
    bad_c12 = c12_df.copy()
    bad_c12.loc[0, "genome"] = g_parent_unregistered
    with pytest.raises(ValueError, match=g_parent_unregistered):
        reconcile_and_summarize_order_metabolic_tables(
            thermo_meta, mod_df, bad_c12, c18_df, prop_df, alias_map, primary_medium="native"
        )


def test_stage3_real_outputs_if_present() -> None:
    """Gate 8: Pending real-output integration check (runs only when Stage 3 output TSVs exist)."""
    import json

    res_dir = Path(os.environ.get("EQ_SCIENCE_OUT", str(CODE_DIR.parent))) / "results"
    eval_path = res_dir / "task_c_targeted_family_functional_evaluation.tsv"
    loo_path = res_dir / "task_c_masked_benchmark_leave_order_out_threshold_calibration.tsv"
    mem_path = res_dir / "task_c_stage3_member_locus_manifest_verified.tsv"
    sum_path = res_dir / "stage3_summary.json"
    if not eval_path.exists() or not loo_path.exists():
        pytest.skip("PENDING_REAL_OUTPUT_INTEGRATION: Stage 3 real output TSVs not yet generated in results/")

    assert mem_path.exists(), "Missing task_c_stage3_member_locus_manifest_verified.tsv"
    assert sum_path.exists(), "Missing stage3_summary.json"

    ev_df = pd.read_csv(eval_path, sep="\t")
    assert ev_df["family_id"].is_unique
    n_new_dark = int((ev_df["candidate_role"] == "de_novo_dark_candidate_family").sum())
    assert 0 <= n_new_dark <= 10, f"Expected 0..10 new dark candidate families, got {n_new_dark}"

    summary = json.loads(sum_path.read_text())
    if n_new_dark == 0:
        assert summary.get("selection_outcome_status") == "no_eligible_dark_family_after_strict_audit"
        assert summary.get("no_eligible_family") is True

    # Prior heat-stress candidates must never be labelled annotated_positive_control
    prior_rows = ev_df[ev_df["family_id"].isin(PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES)]
    assert (prior_rows["candidate_role"] == "prior_heat_stress_candidate_example").all()

    # Any row labelled annotated_positive_control must have has_curated_positive_control_evidence == True
    pos_rows = ev_df[ev_df["candidate_role"] == "annotated_positive_control"]
    if not pos_rows.empty:
        assert pos_rows["has_curated_positive_control_evidence"].all()

    # Verify all selected supporting member loci were independently validated
    mem_df = pd.read_csv(mem_path, sep="\t")
    assert not mem_df.empty
    assert "member_final_promotion_gate_pass" in mem_df.columns
    assert "genome_manifest_provenance_status" in mem_df.columns

    # Check leave-order-out threshold calibration table
    loo_df = pd.read_csv(loo_path, sep="\t")
    pooled_rows = loo_df[loo_df["held_out_order"] == "ALL_POOLED_IN_SAMPLE"]
    assert not pooled_rows.empty
    assert (pooled_rows["is_held_out_evaluation"] == False).all()  # noqa: E712
    held_rows = loo_df[loo_df["held_out_order"] != "ALL_POOLED_IN_SAMPLE"]
    assert (held_rows["is_held_out_evaluation"] == True).all()  # noqa: E712
    no_qual = loo_df[loo_df["threshold_calibration_status"] == "no_qualifying_threshold_attained"]
    if not no_qual.empty:
        assert no_qual["calibrated_score_threshold"].isna().all()

    # Check order-stratified metabolic summary and 3-genome GSPA pilot audit
    ord_path = res_dir / "task_c_order_stratified_metabolic_and_gspa_descriptive_summary.tsv"
    gspa_path = res_dir / "task_c_gspa_pilot_3genome_descriptive_audit.tsv"
    assert ord_path.exists(), "Missing task_c_order_stratified_metabolic_and_gspa_descriptive_summary.tsv"
    assert gspa_path.exists(), "Missing task_c_gspa_pilot_3genome_descriptive_audit.tsv"

    ord_df = pd.read_csv(ord_path, sep="\t")
    assert int(ord_df["n_genomes"].sum()) == 237, f"Expected 237 total genomes across orders, got {int(ord_df['n_genomes'].sum())}"
    assert int(ord_df["n_genomes_with_gapseq_model"].sum()) == 236
    assert int(ord_df["n_genomes_missing_gapseq_model"].sum()) == 1

    gspa_df = pd.read_csv(gspa_path, sep="\t")
    f23_cnn = gspa_df[(gspa_df["pilot_genome"] == "F23PRr1__single_SemiBin_76") & (gspa_df["gspa_mode"] == "cnn")]
    assert len(f23_cnn) == 1
    assert int(f23_cnn["gspa_nonroot_05_proteins"].iloc[0]) == 3425
    assert int(f23_cnn["gspa_nonroot_05_dark_proteins"].iloc[0]) == 1060
    f23_fast = gspa_df[(gspa_df["pilot_genome"] == "F23PRr1__single_SemiBin_76") & (gspa_df["gspa_mode"] == "fast")]
    assert int(f23_fast["gspa_nonroot_05_proteins"].iloc[0]) == 1350
    assert int(f23_fast["gspa_nonroot_05_dark_proteins"].iloc[0]) == 0
    f23_base = gspa_df[(gspa_df["pilot_genome"] == "F23PRr1__single_SemiBin_76") & (gspa_df["gspa_mode"] == "cnn_baseline")]
    assert int(f23_base["gspa_nonroot_05_proteins"].iloc[0]) == 1768
    assert int(f23_base["gspa_nonroot_05_dark_proteins"].iloc[0]) == 169

