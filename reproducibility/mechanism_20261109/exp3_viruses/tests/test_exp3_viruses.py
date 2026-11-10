#!/usr/bin/env python3
"""Unit and adversarial scientific tests for exp3_viruses.

Covers all contract-mandated positive and negative fixtures from:
- `EXP3-CONTRACT.txt`
- `EXP3-STAGING-NOTE.txt`
- `EXP3-EVIDENCE-REVIEW.txt`
- `EXP3-CODE-REVIEW.txt`
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

CODE_DIR = Path(__file__).resolve().parent.parent / "code"
INPUTS_DIR = Path(__file__).resolve().parent.parent / "inputs"
sys.path.insert(0, str(CODE_DIR))

from exp3_common import (
    are_same_namespaced_contig,
    canonical_spacer_seq_rc,
    classify_host_link_evidence_tier,
    evaluate_viral_identity_and_decoys,
    is_valid_crispr_repeat_and_spacer,
    qualify_namespaced_contig_id,
    reconcile_host_alias_strict,
    revcomp_dna,
    score_spacer_alignment,
    sha256_bytes,
    summarize_mag_host_disposition,
    verify_and_reconcile_fasta,
)
from stage1_crispr_recovery import (
    deduplicate_mag_spacers,
    parse_minced_txt,
)
from stage2_viral_identity_and_spacers import (
    parse_checkv_and_genomad_dirs,
)


def test_toy_spacer_mismatch_end_coverage_and_gap_event_vs_residue():
    """Verify exact, 1-mismatch, 2-mismatch, end-coverage rejection, and gapopen event vs residue handling."""
    s0 = score_spacer_alignment(
        qlen=32, qstart=1, qend=32, sstart=100, send=131,
        length=32, mismatch=0, gapopen=0, pident=100.0,
    )
    assert s0["match_tier"] == "exact_0mm_full"
    assert s0["passes_strict_spacer_gate"] is True
    assert s0["effective_mismatches"] == 0
    assert s0["full_spacer_coverage"] == 1.0

    s1 = score_spacer_alignment(
        qlen=32, qstart=1, qend=32, sstart=500, send=469,
        length=32, mismatch=1, gapopen=0, pident=96.875,
    )
    assert s1["match_tier"] == "strict_1mm_end2"
    assert s1["passes_strict_spacer_gate"] is True
    assert s1["effective_mismatches"] == 1

    # Partial end-coverage failure: 7 bases unaligned at 3' end with 100% BLAST internal identity
    s_end_fail = score_spacer_alignment(
        qlen=32, qstart=1, qend=25, sstart=200, send=224,
        length=25, mismatch=0, gapopen=0, pident=100.0,
    )
    assert s_end_fail["uncovered_3p"] == 7
    assert s_end_fail["effective_mismatches"] == 7
    assert s_end_fail["match_tier"] == "rejected_mismatch_or_end_coverage"
    assert s_end_fail["passes_strict_spacer_gate"] is False
    assert s_end_fail["passes_exploratory_spacer_gate"] is False

    # 3 internal mismatches -> rejected
    s_mm_fail = score_spacer_alignment(
        qlen=32, qstart=1, qend=32, sstart=200, send=231,
        length=32, mismatch=3, gapopen=0, pident=90.625,
    )
    assert s_mm_fail["match_tier"] == "rejected_mismatch_or_end_coverage"
    assert s_mm_fail["passes_strict_spacer_gate"] is False

    # Exploratory 2-mismatch hit
    s_exp = score_spacer_alignment(
        qlen=32, qstart=1, qend=32, sstart=200, send=231,
        length=32, mismatch=2, gapopen=0, pident=93.75,
    )
    assert s_exp["match_tier"] == "near_2mm_exploratory"
    assert s_exp["passes_strict_spacer_gate"] is False
    assert s_exp["passes_exploratory_spacer_gate"] is True

    # Gapped alignment with ONLY gapopen (gap event count, no gap residue count `gaps`):
    # Must NOT claim numerical full_spacer_identity from gapopen!
    s_gap_no_res = score_spacer_alignment(
        qlen=32, qstart=1, qend=32, sstart=200, send=227,
        length=32, mismatch=0, gapopen=1, pident=87.5, gaps=None,
    )
    assert s_gap_no_res["match_tier"] == "rejected_gapped_alignment_unsupported"
    assert s_gap_no_res["effective_mismatches"] is None
    assert s_gap_no_res["full_spacer_identity"] is None
    assert s_gap_no_res["passes_strict_spacer_gate"] is False
    assert s_gap_no_res["passes_exploratory_spacer_gate"] is False

    # Gapped alignment WITH total gap residues (`gaps=4`, `gapopen=1`):
    # Must use total gap residues (4), never gapopen event count (1)!
    s_gap_with_res = score_spacer_alignment(
        qlen=32, qstart=1, qend=32, sstart=200, send=227,
        length=32, mismatch=0, gapopen=1, pident=87.5, gaps=4,
    )
    assert s_gap_with_res["match_tier"] == "rejected_gapped_alignment"
    assert s_gap_with_res["gap_residues"] == 4
    assert s_gap_with_res["effective_mismatches"] == 4
    assert s_gap_with_res["full_spacer_identity"] == round((32 - 4) / 32.0, 6)
    assert s_gap_with_res["passes_strict_spacer_gate"] is False


def test_crispr_repeat_and_spacer_low_complexity_rejection():
    """Verify low-complexity repeats/spacers and homopolymers are rejected."""
    ok, _ = is_valid_crispr_repeat_and_spacer("G", "ACGT" * 8)
    assert not ok

    ok2, reason2 = is_valid_crispr_repeat_and_spacer(
        "GTTTGAGAGTTGTGATTTGCAAGTACGT",
        "ACGTAAAAAAAACGTACGTACGTACGTACGT",
    )
    assert not ok2
    assert reason2 == "rejected_homopolymer_run"

    ok3, reason3 = is_valid_crispr_repeat_and_spacer(
        "GTTTGAGAGTTGTGATTTGCAAGTACGT",
        "ATATATATATATATATATATATATATATATAT",
    )
    assert not ok3
    assert reason3 in ("rejected_low_complexity_entropy", "rejected_dinucleotide_tandem_repeat")

    ok4, reason4 = is_valid_crispr_repeat_and_spacer(
        "CTTTCAATTCTATAGTAGATTAGC",
        "GCTAGCTTGACCTGATCGATCGTAGCTAGCTA",
    )
    assert ok4
    assert reason4 == "valid_crispr_spacer"


def test_wrong_host_alias_no_heuristic_suffix_stripping():
    """Verify heuristic suffix stripping (_sub, self_/single_) is strictly rejected."""
    canonical_591 = {"15Dr1__self_SemiBin_98_sub", "41PRr2__single_SemiBin_236", "44Dr2__single_SemiBin_185"}
    crosswalk = {"43PRr1__single_SemiBin_144": "44Dr2__single_SemiBin_185"}

    st1, m1 = reconcile_host_alias_strict("15Dr1__self_SemiBin_98_sub", canonical_591, crosswalk)
    assert st1 == "canonical_591_direct"
    assert m1 == "15Dr1__self_SemiBin_98_sub"

    st2, m2 = reconcile_host_alias_strict("15Dr1__self_SemiBin_98", canonical_591, crosswalk)
    assert st2 == "rejected_unreconciled_host_alias"
    assert m2 is None

    st3, m3 = reconcile_host_alias_strict("41PRr2__SemiBin_236", canonical_591, crosswalk)
    assert st3 == "rejected_unreconciled_host_alias"
    assert m3 is None

    st4, m4 = reconcile_host_alias_strict("43PRr1__single_SemiBin_144", canonical_591, crosswalk)
    assert st4 == "nonrep_refined_bin_mapped_to_species_rep"
    assert m4 == "44Dr2__single_SemiBin_185"


def test_unknown_viral_identity_vs_external_viral_controls_and_curated_nonviral_decoys():
    """Verify EXP3-EVIDENCE-REVIEW and EXP3-CODE-REVIEW rules:
    1. CheckV Not-determined + no hallmark (0 host genes, not curated control) -> Unknown_unverified_viral_identity (NEVER Decoy/Verified_nonviral).
    2. External habitat reference viruses (`decoys_derep.fna`) are STILL viral sequences (`is_verified_viral=True`, `is_decoy=False`),
       recorded separately as `External_habitat_viral_specificity_control` (`External_habitat_viral_control_match_retained`).
    3. Only deliberately selected verified non-viral plasmid/cellular controls get `Decoy_curated_*` labels.
    4. Missing annotations (`None`/`NA`) are not conflated with scored zero (`0`), and overlapping tool counts are not summed.
    """
    # 1. CheckV Not-determined + 0 viral hallmarks + 0 host genes -> Unknown_unverified_viral_identity
    unk_eval = evaluate_viral_identity_and_decoys(
        target_id="59PRr1::59PRr1_225981",
        target_hit_lo=1200,
        target_hit_hi=1232,
        host_mag="S54PRr1__self_SemiBin_338",
        host_contig_id="S54PRr1_88",
        host_array_lo=200,
        host_array_hi=800,
        target_meta={
            "in_plasmid_pool": False,
            "is_curated_plasmid_control": False,
            "is_curated_cellular_control": False,
            "in_cellular_mag_contig_index": False,
            "plasmid_hallmarks": 0,
            "virus_hallmarks": 0,
            "archaeal_viral_hallmarks": 0,
            "checkv_viral_genes": 0,
            "checkv_host_genes": 0,
            "checkv_quality": "Not-determined",
        },
    )
    assert unk_eval["viral_identity_status"] == "Unknown_unverified_viral_identity"
    assert unk_eval["is_unknown_viral_identity"] is True
    assert unk_eval["is_decoy"] is False
    assert unk_eval["is_verified_viral"] is False

    tier_unk = classify_host_link_evidence_tier(
        host_reconciliation_status="canonical_591_direct",
        viral_eval=unk_eval,
        n_strict_spacers=1,
        n_exact_0mm_spacers=1,
        n_exploratory_spacers=1,
        has_cas_operon_in_mag=True,
        has_verified_provirus_junction_in_mag=False,
        has_classifier_corroboration=False,
        n_distinct_host_genera_for_virus=1,
        n_distinct_host_phyla_for_virus=1,
    )
    assert tier_unk["evidence_tier"] == "Tier3_unverified_viral_identity_retained"
    assert tier_unk["eligible_stage4_primary"] is False
    assert tier_unk["eligible_stage4_sensitivity"] is False

    # 2. External habitat reference virus (e.g. RefSeq phage from non-desert habitat):
    # STILL a viral sequence (is_verified_viral=True, is_decoy=False), recorded separately for mapping specificity!
    ext_viral_eval = evaluate_viral_identity_and_decoys(
        target_id="NC_001416.1_Enterobacteria_phage_lambda",
        target_hit_lo=100,
        target_hit_hi=132,
        host_mag="S54PRr1__self_SemiBin_338",
        host_contig_id="S54PRr1_88",
        host_array_lo=200,
        host_array_hi=800,
        target_meta={
            "is_external_habitat_viral_control": True,
            "virus_hallmarks": 5,
            "checkv_viral_genes": 20,
            "checkv_quality": "Complete",
        },
    )
    assert ext_viral_eval["viral_identity_status"] == "Verified_viral_external_habitat_reference"
    assert ext_viral_eval["habitat_specificity_role"] == "External_habitat_viral_specificity_control"
    assert ext_viral_eval["is_verified_viral"] is True
    assert ext_viral_eval["is_decoy"] is False
    assert ext_viral_eval["is_external_viral_specificity_control"] is True

    tier_ext_viral = classify_host_link_evidence_tier(
        host_reconciliation_status="canonical_591_direct",
        viral_eval=ext_viral_eval,
        n_strict_spacers=1,
        n_exact_0mm_spacers=1,
        n_exploratory_spacers=1,
        has_cas_operon_in_mag=True,
        has_verified_provirus_junction_in_mag=False,
        has_classifier_corroboration=False,
        n_distinct_host_genera_for_virus=1,
        n_distinct_host_phyla_for_virus=1,
    )
    assert tier_ext_viral["evidence_tier"] == "External_habitat_viral_control_match_retained"
    assert tier_ext_viral["eligible_stage4_primary"] is False

    # 3. Deliberately selected curated plasmid control -> Decoy_curated_plasmid_control
    plasmid_eval = evaluate_viral_identity_and_decoys(
        target_id="12Dr1::12Dr1_99999",
        target_hit_lo=500,
        target_hit_hi=532,
        host_mag="15Dr1__self_SemiBin_98_sub",
        host_contig_id="15Dr1_100",
        host_array_lo=1000,
        host_array_hi=1500,
        target_meta={
            "is_curated_plasmid_control": True,
            "plasmid_hallmarks": 2,
            "virus_hallmarks": 0,
            "checkv_viral_genes": 0,
            "checkv_host_genes": 5,
            "checkv_quality": "Not-determined",
        },
    )
    assert plasmid_eval["is_decoy"] is True
    assert plasmid_eval["viral_identity_status"] == "Decoy_curated_plasmid_control"

    # 4. Overlapping Foldseek + geNomad + CheckV counts on a single ORF (each = 1, distinct_viral_orf_count = 1):
    # Must NOT sum to 3 independent viral genes!
    overlap_eval = evaluate_viral_identity_and_decoys(
        target_id="44Dr2::44Dr2_fragment_99",
        target_hit_lo=100,
        target_hit_hi=132,
        host_mag="41PRr2__single_SemiBin_236",
        host_contig_id="41PRr2_10",
        host_array_lo=500,
        host_array_hi=800,
        target_meta={
            "virus_hallmarks": 1,
            "archaeal_viral_hallmarks": 1,
            "foldseek_viral_hallmarks": 1,
            "checkv_viral_genes": 1,
            "checkv_host_genes": 0,
            "checkv_quality": "Low-quality",
            "genomad_virus_score": 0.85,
        },
    )
    assert overlap_eval["min_independent_viral_orfs"] == 1
    assert overlap_eval["is_verified_viral"] is False
    assert overlap_eval["is_provisional_viral"] is True
    assert overlap_eval["viral_identity_status"] == "Provisional_viral_low_completeness_fragment"

    # 5. Missing/unscored annotations (None / NA) vs scored zero:
    missing_eval = evaluate_viral_identity_and_decoys(
        target_id="99Dr1::99Dr1_unbinned_contig_7",
        target_hit_lo=100,
        target_hit_hi=132,
        host_mag="41PRr2__single_SemiBin_236",
        host_contig_id="41PRr2_10",
        host_array_lo=500,
        host_array_hi=800,
        target_meta={
            "in_cellular_mag_contig_index": True,
            "virus_hallmarks": None,
            "checkv_viral_genes": "NA",
            "checkv_host_genes": None,
            "checkv_quality": "NA",
        },
    )
    assert missing_eval["is_unknown_viral_identity"] is True
    assert missing_eval["is_cellular_like"] is False
    assert missing_eval["status_detail"] == "missing_unscored_viral_and_checkv_annotations_not_zero"


def test_namespaced_sample_contig_collision():
    """Verify same bare contig suffix in different samples NEVER collides (`sampleA::contig_1` != `sampleB::contig_1`)."""
    s1, b1, q1, _, _ = qualify_namespaced_contig_id("sampleA::contig_1|provirus_100_900")
    s2, b2, q2, _, _ = qualify_namespaced_contig_id("contig_1", fallback_sample="sampleB")
    assert (s1, b1, q1) == ("sampleA", "contig_1", "sampleA::contig_1")
    assert (s2, b2, q2) == ("sampleB", "contig_1", "sampleB::contig_1")

    same_ab, qa, qb = are_same_namespaced_contig(
        target_id="sampleA::contig_1|provirus_34656_72876",
        host_contig_id="contig_1",
        host_mag="sampleB__single_SemiBin_136",
    )
    assert same_ab is False
    assert qa == "sampleA::contig_1"
    assert qb == "sampleB::contig_1"

    # Even if coordinates overlap `host_array_lo..host_array_hi`, a hit on `sampleA::contig_1`
    # from a spacer on `sampleB::contig_1` is in a DIFFERENT sample and must NOT be masked as
    # a same-contig technical self-CRISPR-array match!
    cross_sample_eval = evaluate_viral_identity_and_decoys(
        target_id="sampleA::contig_1|provirus_34656_72876",
        target_hit_lo=1001,  # parent coord = 35656
        target_hit_hi=1032,
        host_mag="sampleB__single_SemiBin_136",
        host_contig_id="contig_1",
        host_array_lo=35500,
        host_array_hi=36200,
        target_meta={
            "contig_in_same_host_mag": False,
            "virus_hallmarks": 2,
            "checkv_viral_genes": 8,
            "checkv_host_genes": 3,
            "checkv_quality": "High-quality",
            "has_verified_host_flanks": True,
            "verified_viral_region_lo": 34656,
            "verified_viral_region_hi": 72876,
        },
    )
    assert cross_sample_eval["is_same_namespaced_contig"] is False
    assert cross_sample_eval["is_technical_crispr_array_self_match"] is False
    assert cross_sample_eval["is_verified_viral"] is True


def test_same_chromosomal_contig_provirus_and_technical_self_array_mask():
    """Verify technical self-CRISPR mask, chromosomal non-provirus, unexcised provirus candidate, and valid provirus."""
    # Case 1: Technical self-CRISPR array match on same namespaced contig -> masked as Technical_self_crispr_array_masked
    self_hit = evaluate_viral_identity_and_decoys(
        target_id="V45Dr1::V45Dr1_606357|provirus_34656_72876",
        target_hit_lo=1001,  # parent coord = 35656, overlaps array 35500..36200
        target_hit_hi=1032,
        host_mag="V45Dr1__single_SemiBin_136",
        host_contig_id="V45Dr1_606357",
        host_array_lo=35500,
        host_array_hi=36200,
        target_meta={
            "contig_in_same_host_mag": True,
            "virus_hallmarks": 1,
            "checkv_viral_genes": 3,
            "checkv_host_genes": 5,
            "checkv_quality": "Medium-quality",
            "has_verified_host_flanks": True,
            "verified_viral_region_lo": 34656,
            "verified_viral_region_hi": 72876,
        },
    )
    assert self_hit["is_same_namespaced_contig"] is True
    assert self_hit["is_technical_crispr_array_self_match"] is True
    assert self_hit["is_decoy"] is False
    assert self_hit["viral_identity_status"] == "Technical_self_crispr_array_masked"

    # Case 2: Same chromosomal contig in host MAG with ZERO viral genes/hallmarks -> Cellular_like_chromosomal_contig_not_provirus
    chrom_not_prov = evaluate_viral_identity_and_decoys(
        target_id="V45Dr1::V45Dr1_606357",
        target_hit_lo=10000,
        target_hit_hi=10032,
        host_mag="V45Dr1__single_SemiBin_136",
        host_contig_id="V45Dr1_999999",
        host_array_lo=500,
        host_array_hi=900,
        target_meta={
            "contig_in_same_host_mag": True,
            "virus_hallmarks": 0,
            "checkv_viral_genes": 0,
            "checkv_host_genes": 40,
            "checkv_quality": "Not-determined",
            "has_verified_host_flanks": False,
            "verified_viral_region_lo": None,
            "verified_viral_region_hi": None,
        },
    )
    assert chrom_not_prov["is_cellular_like"] is True
    assert chrom_not_prov["viral_identity_status"] == "Cellular_like_chromosomal_contig_not_provirus"

    # Case 3: Unexcised contig in host MAG carrying viral hallmarks AND host genes -> Unresolved_unexcised_provirus_candidate
    unexcised_prov = evaluate_viral_identity_and_decoys(
        target_id="V45Dr1::V45Dr1_606357",
        target_hit_lo=10000,
        target_hit_hi=10032,
        host_mag="V45Dr1__single_SemiBin_136",
        host_contig_id="V45Dr1_999999",
        host_array_lo=500,
        host_array_hi=900,
        target_meta={
            "contig_in_same_host_mag": True,
            "virus_hallmarks": 2,
            "checkv_viral_genes": 5,
            "checkv_host_genes": 40,
            "checkv_quality": "Low-quality",
            "has_verified_host_flanks": False,
            "verified_viral_region_lo": None,
            "verified_viral_region_hi": None,
        },
    )
    assert unexcised_prov["is_provisional_viral"] is True
    assert unexcised_prov["is_cellular_like"] is False
    assert unexcised_prov["viral_identity_status"] == "Unresolved_unexcised_provirus_candidate"

    # Case 4: Spacer hit lands in flanking host region outside excised proviral interval -> Cellular_like_flanking_host_region_outside_provirus
    flank_hit = evaluate_viral_identity_and_decoys(
        target_id="V45Dr1::V45Dr1_606357",
        target_hit_lo=10000,
        target_hit_hi=10032,
        host_mag="V45Dr1__single_SemiBin_136",
        host_contig_id="V45Dr1_999999",
        host_array_lo=500,
        host_array_hi=900,
        target_meta={
            "contig_in_same_host_mag": True,
            "virus_hallmarks": 3,
            "checkv_viral_genes": 10,
            "checkv_host_genes": 20,
            "checkv_quality": "High-quality",
            "has_verified_host_flanks": True,
            "verified_viral_region_lo": 34656,
            "verified_viral_region_hi": 72876,
        },
    )
    assert flank_hit["is_cellular_like"] is True
    assert flank_hit["viral_identity_status"] == "Cellular_like_flanking_host_region_outside_provirus"

    # Case 5: Valid excised provirus with host flanks, viral hallmarks, and interior hit far from CRISPR array
    valid_prov = evaluate_viral_identity_and_decoys(
        target_id="V45Dr1::V45Dr1_606357|provirus_34656_72876",
        target_hit_lo=15000,
        target_hit_hi=15032,
        host_mag="V45Dr1__single_SemiBin_136",
        host_contig_id="V45Dr1_606357",
        host_array_lo=1000,
        host_array_hi=1800,
        target_meta={
            "contig_in_same_host_mag": True,
            "virus_hallmarks": 3,
            "checkv_viral_genes": 12,
            "checkv_host_genes": 8,
            "checkv_quality": "High-quality",
            "has_verified_host_flanks": True,
            "verified_viral_region_lo": 34656,
            "verified_viral_region_hi": 72876,
        },
    )
    assert valid_prov["is_decoy"] is False
    assert valid_prov["is_verified_viral"] is True
    assert valid_prov["provirus_boundary_status"] == "verified_provirus_with_host_flanks"


def test_multiple_targets_and_multihost_ambiguity_tiering():
    """Verify ambiguous multi-host links are retained in Tier 3 and excluded from Stage 4."""
    verified_viral = {
        "viral_identity_status": "Verified_viral_independent_hallmark_and_checkv",
        "is_verified_viral": True,
        "is_provisional_viral": False,
        "is_unknown_viral_identity": False,
        "is_cellular_like": False,
        "is_technical_crispr_array_self_match": False,
        "is_decoy": False,
        "is_external_viral_specificity_control": False,
    }
    t1 = classify_host_link_evidence_tier(
        host_reconciliation_status="canonical_591_direct",
        viral_eval=verified_viral,
        n_strict_spacers=2,
        n_exact_0mm_spacers=1,
        n_exploratory_spacers=2,
        has_cas_operon_in_mag=True,
        has_verified_provirus_junction_in_mag=False,
        has_classifier_corroboration=False,
        n_distinct_host_genera_for_virus=1,
        n_distinct_host_phyla_for_virus=1,
    )
    assert t1["evidence_tier"] == "Tier1_verified_direct"
    assert t1["eligible_stage4_primary"] is True
    assert t1["eligible_stage4_sensitivity"] is True

    t3_multi = classify_host_link_evidence_tier(
        host_reconciliation_status="canonical_591_direct",
        viral_eval=verified_viral,
        n_strict_spacers=3,
        n_exact_0mm_spacers=2,
        n_exploratory_spacers=3,
        has_cas_operon_in_mag=True,
        has_verified_provirus_junction_in_mag=False,
        has_classifier_corroboration=True,
        n_distinct_host_genera_for_virus=2,
        n_distinct_host_phyla_for_virus=1,
    )
    assert t3_multi["evidence_tier"] == "Tier3_ambiguous_multihost_retained"
    assert t3_multi["eligible_stage4_primary"] is False
    assert t3_multi["eligible_stage4_sensitivity"] is False


def test_adversarial_fasta_reconciliation_non_circular_and_same_length_mutation():
    """Verify non-circular canonical FASTA verification and adversarial same-length sequence mutation rejection."""
    ref_fasta = b">contig_1\nACGTACGT\n>contig_2\nTTTTGGGG\n"
    ref_sha = sha256_bytes(ref_fasta)

    # 1. Primary canonical check (verified_canonical_reference_bytes=None) passes on exact match
    v_exact = verify_and_reconcile_fasta("sampleA__self_SemiBin_1", ref_fasta, ref_sha, 2, 16, None)
    assert v_exact["is_verified_canonical"] is True
    assert v_exact["status"] == "verified_exact_match_genome_manifest"

    # 2. Adversarial same-length sequence mutation (T -> A at end of contig_1; same headers, 2 contigs, 16 bases):
    # Must FAIL primary canonical check!
    mutated_same_len_fasta = b">contig_1\nACGTACGA\n>contig_2\nTTTTGGGG\n"
    v_mut_primary = verify_and_reconcile_fasta(
        "sampleA__self_SemiBin_1", mutated_same_len_fasta, ref_sha, 2, 16, None
    )
    assert v_mut_primary["is_verified_canonical"] is False
    assert v_mut_primary["status"] == "rejected_canonical_sha256_sequence_mismatch"

    # 3. Adversarial circularity check: attempting to pass `mutated_same_len_fasta` as its own
    # `verified_canonical_reference_bytes` against `ref_sha` MUST raise ValueError!
    with pytest.raises(ValueError, match="failed independent manifest SHA-256"):
        verify_and_reconcile_fasta(
            "sampleA__self_SemiBin_1",
            mutated_same_len_fasta,
            ref_sha,
            2,
            16,
            verified_canonical_reference_bytes=mutated_same_len_fasta,
        )

    # 4. Even when true `ref_fasta` is passed as `verified_canonical_reference_bytes`,
    # `mutated_same_len_fasta` MUST still fail because per-contig sequence SHA-256 differs!
    v_mut_with_ref = verify_and_reconcile_fasta(
        "sampleA__self_SemiBin_1",
        mutated_same_len_fasta,
        ref_sha,
        2,
        16,
        verified_canonical_reference_bytes=ref_fasta,
    )
    assert v_mut_with_ref["is_verified_canonical"] is False
    assert v_mut_with_ref["status"] == "rejected_canonical_sha256_sequence_mismatch"

    # 5. Registered header transform (`sampleA::contig_1`) with identical sequences passes ONLY when
    # independently hash-verified `ref_fasta` is provided:
    prefixed_fasta = b">sampleA::contig_1 description\nACGT\nACGT\n>sampleA::contig_2\nTTTTGGGG\n"
    v_trans = verify_and_reconcile_fasta(
        "sampleA__self_SemiBin_1",
        prefixed_fasta,
        ref_sha,
        2,
        16,
        verified_canonical_reference_bytes=ref_fasta,
    )
    assert v_trans["is_verified_canonical"] is True
    assert v_trans["status"] == "verified_identical_contigs_registered_header_transform"
    assert v_trans["header_transform_applied"] == "strip_sample_double_colon"

    # 6. Wrong sample prefix (`sampleB::contig_1` on `sampleA__self_SemiBin_1`) is rejected!
    wrong_sample_prefix_fasta = b">sampleB::contig_1\nACGTACGT\n>sampleB::contig_2\nTTTTGGGG\n"
    v_wrong_sample = verify_and_reconcile_fasta(
        "sampleA__self_SemiBin_1",
        wrong_sample_prefix_fasta,
        ref_sha,
        2,
        16,
        verified_canonical_reference_bytes=ref_fasta,
    )
    assert v_wrong_sample["is_verified_canonical"] is False

    # 7. Altered contig membership (extra contig) is rejected!
    unrefined_fasta = b">contig_1\nACGTACGT\n>contig_2\nTTTTGGGG\n>contig_3_extra\nCCCCAAAA\n"
    v_wrong = verify_and_reconcile_fasta(
        "sampleA__self_SemiBin_1",
        unrefined_fasta,
        ref_sha,
        2,
        16,
        verified_canonical_reference_bytes=ref_fasta,
    )
    assert v_wrong["is_verified_canonical"] is False
    assert v_wrong["status"] == "rejected_wrong_product_contig_membership_differs"


def test_canonical_591_archaeal_mags_and_dispositions():
    """Verify staged canonical taxonomy has exactly 591 archaeal MAGs (237 Thermoplasmatota) and dispositions."""
    tax_path = INPUTS_DIR / "taxonomy_reconciliation_all5438.tsv.gz"
    assert tax_path.exists(), f"Missing staged taxonomy file: {tax_path}"
    df = pd.read_csv(tax_path, sep="\t")
    assert len(df) == 5438
    arch = df[df["canonical_domain"] == "d__Archaea"]
    assert len(arch) == 591
    thermo = arch[arch["canonical_phylum"] == "p__Thermoplasmatota"]
    assert len(thermo) == 237

    assert summarize_mag_host_disposition("m1", False, False, 0, 0, []) == "Unknown_unscreened"
    assert summarize_mag_host_disposition("m2", True, True, 0, 0, []) == "Zero_spacers_detected"
    assert summarize_mag_host_disposition("m3", True, True, 2, 14, []) == "Spacers_detected_no_viral_hit"
    assert (
        summarize_mag_host_disposition("m4", True, True, 1, 5, ["Tier3_unverified_viral_identity_retained"])
        == "Candidate_link_unknown_viral_identity_retained"
    )
    assert (
        summarize_mag_host_disposition("m5", True, True, 1, 5, ["Technical_self_crispr_array_masked"])
        == "Spacers_detected_only_technical_self_array_masked"
    )
    assert (
        summarize_mag_host_disposition("m6", True, True, 2, 8, ["Tier1_verified_direct", "Decoy_control_rejected"])
        == "Linked_Tier1_verified_direct"
    )
    assert (
        summarize_mag_host_disposition("m7", True, True, 1, 4, ["External_habitat_viral_control_match_retained"])
        == "External_habitat_viral_reference_match_only"
    )
    assert (
        summarize_mag_host_disposition("m8", True, False, 0, 0, [], caller_failed=True)
        == "Failed_caller_execution"
    )


def test_stage1_spacer_dedup_rc_and_multi_array_loci_and_header_namespace():
    """Verify EXP3-STAGE1-REVIEW rules:
    1. Embedded contig header sample prefix (`41PRr2_100`) takes precedence over MAG fallback sample (`99Dr1`).
    2. Repeated spacer in 2 distant arrays (100..300 and 50000..50200) + reverse complement in a 3rd array
       deduplicates to ONE independent spacer (`n_distinct_arrays_in_mag=3`) and preserves individual array
       intervals (`all_array_loci`) rather than collapsing to `100..50200`.
    3. Filename `|provirus_start_end` without `verified_viral_region_lo`/`has_verified_host_flanks` is
       `predicted_excised_provirus_interval_unverified_flanks`, not verified host flanks.
    """
    # 1. Embedded sample prefix priority over fallback sample
    samp, base, qual, _, _ = qualify_namespaced_contig_id("41PRr2_100", fallback_sample="99Dr1")
    assert (samp, base, qual) == ("41PRr2", "41PRr2_100", "41PRr2::41PRr2_100")

    # 2. MinCED parser + deduplication across multiple arrays & reverse complement
    spc1 = "GCTAGCTTGACCTGATCGATCGTAGCTAGCTA"
    spc1_rc = revcomp_dna(spc1)
    can1, _ = canonical_spacer_seq_rc(spc1)
    can1_rc, _ = canonical_spacer_seq_rc(spc1_rc)
    assert can1 == can1_rc

    minced_sample_txt = f"""Sequence '41PRr2_100' (80000 bp)

CRISPR 1   Range: 100 - 300
POSITION\tREPEAT\t\t\t\tSPACER
--------\t----------------------------\t--------------------------------
100\t\tGTTTGAGAGTTGTGATTTGCAAGTACGT\t{spc1}\t[ 28, 32 ]
160\t\tGTTTGAGAGTTGTGATTTGCAAGTACGT
--------\t----------------------------\t--------------------------------

CRISPR 2   Range: 50000 - 50200
POSITION\tREPEAT\t\t\t\tSPACER
--------\t----------------------------\t--------------------------------
50000\t\tGTTTGAGAGTTGTGATTTGCAAGTACGT\t{spc1_rc}\t[ 28, 32 ]
50060\t\tGTTTGAGAGTTGTGATTTGCAAGTACGT
--------\t----------------------------\t--------------------------------
"""
    arrays, spacers = parse_minced_txt(minced_sample_txt, "99Dr1__self_SemiBin_1", "99Dr1", True)
    assert len(arrays) == 2
    assert arrays[0]["namespaced_contig_id"] == "41PRr2::41PRr2_100"
    assert len(spacers) == 2

    dedup = deduplicate_mag_spacers(spacers)
    assert len(dedup) == 1
    assert dedup[0]["n_distinct_arrays_in_mag"] == 2
    assert dedup[0]["observed_orientations"] == "fwd+rc"
    assert dedup[0]["all_array_loci"] == "41PRr2::41PRr2_100:100-300;41PRr2::41PRr2_100:50000-50200"
    # Primary array interval must be the first individual array (100..300), NEVER synthetic 100..50200!
    assert (dedup[0]["primary_array_start"], dedup[0]["primary_array_end"]) == (100, 300)

    # 3. Filename `|provirus_start_end` without verified host flanks is NOT promoted to verified_provirus_with_host_flanks
    unver_flank = evaluate_viral_identity_and_decoys(
        target_id="V45Dr1::V45Dr1_606357|provirus_34656_72876",
        target_hit_lo=15000,
        target_hit_hi=15032,
        host_mag="S54PRr1__self_SemiBin_338",
        host_contig_id="S54PRr1_10",
        host_array_lo=100,
        host_array_hi=500,
        target_meta={
            "virus_hallmarks": 2,
            "checkv_viral_genes": 5,
            "checkv_host_genes": 0,
            "checkv_quality": "Medium-quality",
            "has_verified_host_flanks": False,
        },
    )
    assert unver_flank["provirus_boundary_status"] == "predicted_excised_provirus_interval_unverified_flanks"
    assert unver_flank["read_junction_support_status"] == "untested_missing"


def test_stage2_sha256_bound_annotation_and_edges_schema(tmp_path: Path):
    """Verify EXP3-STAGE2-REVIEW rules:
    1. CheckV + geNomad annotations are bound strictly by `seq_sha256` (`bound_seq_sha256`),
       and `has_verified_host_flanks` requires BOTH host + viral regions in `contamination.tsv`
       AND `checkv_host_genes > 0` AND `checkv_viral_genes > 0` (never filename coordinates alone).
    2. `edges.tsv` uses `votu_rep`, `host_mag`, `host_phylum`, `n_spacers` and has 8 Thermoplasmatota rows.
    """
    edges_p = INPUTS_DIR / "edges.tsv"
    assert edges_p.exists()
    edf = pd.read_csv(edges_p, sep="\t")
    assert {"host_mag", "votu_rep", "host_phylum", "n_spacers"} <= set(edf.columns)
    thermo_e = edf[edf["host_phylum"].isin(["Thermoplasmatota", "p__Thermoplasmatota"])]
    assert thermo_e["host_mag"].nunique() == 8
    assert len(thermo_e) == 14

    # Test parse_checkv_and_genomad_dirs sha256 binding and flank verification
    cv_dir = tmp_path / "checkv"
    cv_dir.mkdir()
    (cv_dir / "quality_summary.tsv").write_text(
        "contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n"
        "S1::S1_100|provirus_1000_9000\t8001\tNo\tNA\t10\t4\t0\tMedium-quality\tGenome-fragment\t55.0\tAAI\t0.0\t1.0\t\n"
        "S2::S2_200|provirus_2000_12000\t10001\tYes\t8000\t15\t6\t3\tHigh-quality\tGenome-fragment\t92.0\tAAI\t10.0\t1.0\t\n"
    )
    (cv_dir / "contamination.tsv").write_text(
        "contig_id\tcontig_length\ttotal_genes\tviral_genes\thost_genes\tprovirus\tproviral_length\thost_length\tregion_types\tregion_lengths\tregion_coords_bp\tregion_coords_genes\tregion_viral_genes\tregion_host_genes\n"
        "S1::S1_100|provirus_1000_9000\t8001\t10\t4\t0\tNo\tNA\tNA\tviral\t8001\t1-8001\t1-10\t4\t0\n"
        "S2::S2_200|provirus_2000_12000\t10001\t15\t6\t3\tYes\t8000\t2001\thost,viral\t2000,8001\t1-2000,2001-10001\t1-3,4-15\t0,6\t3,0\n"
    )
    sha1 = sha256_bytes(b"ACGT" * 200)
    sha2 = sha256_bytes(b"TGCA" * 250)
    seq_sha_map = {
        "S1::S1_100|provirus_1000_9000": sha1,
        "S2::S2_200|provirus_2000_12000": sha2,
    }
    parsed = parse_checkv_and_genomad_dirs(cv_dir, None, None, seq_sha_map, "unit_test")
    assert sha1 in parsed and sha2 in parsed
    # S1 has |provirus_1000_9000 in filename, but NO host region / host genes -> has_verified_host_flanks must NOT be True!
    assert parsed[sha1].get("has_verified_host_flanks", False) is False
    # S2 has both host and viral regions in contamination.tsv + checkv_host_genes=3 + checkv_viral_genes=6 -> True!
    assert parsed[sha2].get("has_verified_host_flanks", False) is True


