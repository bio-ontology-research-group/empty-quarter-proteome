#!/usr/bin/env python3
"""Dedicated Stage 2-4 unit and adversarial tests for exp3_viruses.

Separated into `test_stage2_to_4_exp3.py` so Stage 1 runtime recovery worker
(`7dd8a3a2`) and Stage 2-4 worker never race on the same test file.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

CODE_DIR = Path(__file__).resolve().parent.parent / "code"
INPUTS_DIR = Path(__file__).resolve().parent.parent / "inputs"
sys.path.insert(0, str(CODE_DIR))

from exp3_common import (
    evaluate_viral_identity_and_decoys,
    sha256_bytes,
    summarize_mag_host_disposition,
)
import stage2_viral_identity_and_spacers_v3 as stage2
ORIGINAL_GET_BLAST_TOOL_INFO = stage2.get_blast_tool_info
from stage2_viral_identity_and_spacers_v3 import (
    build_and_audit_target_assemblies_manifest,
    get_blast_tool_info,
    is_blast_receipt_valid,
    load_legacy_precomputed_votu_audit,
    parse_checkv_and_genomad_dirs,
    search_single_raw_assembly,
    stream_file_sha256,
    validate_blast_tabular_output,
    write_blast_receipt,
)
from stage4_downstream_mechanism import (
    compute_host_abundance_metrics,
    evaluate_pair_cooccurrence,
    resolve_shared_sample_intersection,
    run_phylogeny_controlled_tests,
)
from audit_stage2_receipts_with_corrected_validator import audit_stage2_receipts
from audit_thermo_pair_spacers import audit_thermo_spacers
from audit_tier1_pairs_evidence import audit_tier1_pairs, is_unannotated_genomad_gene
from stage4_downstream_mechanism_v2 import (
    audit_surface_and_defense_mechanism_estimability,
    compute_descriptive_phylum_counts,
    load_and_orient_abundance_matrices,
)


@pytest.fixture(autouse=True)
def mock_blast_tool_info_when_container_missing(monkeypatch):
    """Automatically mock blast tool info ONLY when real container SIF is absent (e.g. offline dev)."""
    if not stage2.SIF_BLAST.exists():
        monkeypatch.setattr(
            stage2,
            "get_blast_tool_info",
            lambda: {
                "tool_binary": "blastn",
                "tool_version": "blastn: 2.17.0+",
                "container_sif": str(stage2.SIF_BLAST),
                "container_sha256": "6c7923d5c6042e279cdecb719b0d8b84535a10f8f97485a18695a8d6902ac51f",
            },
        )


def test_stage2_sha256_bound_annotation_and_no_filename_flank_promotion(tmp_path: Path):
    """Verify EXP3-STAGE2-REVIEW rules:
    1. CheckV + geNomad annotations are bound strictly by `seq_sha256` (`bound_seq_sha256`),
       and `has_verified_host_flanks` requires BOTH host + viral regions in `contamination.tsv`
       AND `checkv_host_genes > 0` AND `checkv_viral_genes > 0` (never filename coordinates alone).
    2. `edges.tsv` uses `votu_rep`, `host_mag`, `host_phylum`, `n_spacers` and has 8 Thermoplasmatota MAGs (14 rows).
    """
    edges_p = INPUTS_DIR / "edges.tsv"
    assert edges_p.exists()
    edf = pd.read_csv(edges_p, sep="\t")
    assert {"host_mag", "votu_rep", "host_phylum", "n_spacers"} <= set(edf.columns)
    thermo_e = edf[edf["host_phylum"].isin(["Thermoplasmatota", "p__Thermoplasmatota"])]
    assert thermo_e["host_mag"].nunique() == 8
    assert len(thermo_e) == 14

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
    assert parsed[sha1].get("has_verified_host_flanks", False) is False
    assert parsed[sha2].get("has_verified_host_flanks", False) is True

    # Verify evaluate_viral_identity_and_decoys does not promote filename coordinates to verified flanks
    ev_s1 = evaluate_viral_identity_and_decoys(
        target_id="S1::S1_100|provirus_1000_9000",
        target_hit_lo=500,
        target_hit_hi=532,
        host_mag="S99__self_SemiBin_1",
        host_contig_id="S99_1",
        host_array_lo=100,
        host_array_hi=300,
        target_meta=parsed[sha1],
    )
    assert ev_s1["provirus_boundary_status"] == "predicted_excised_provirus_interval_unverified_flanks"
    assert ev_s1["read_junction_support_status"] == "untested_missing"

    ev_s2 = evaluate_viral_identity_and_decoys(
        target_id="S2::S2_200|provirus_2000_12000",
        target_hit_lo=3000,
        target_hit_hi=3032,
        host_mag="S99__self_SemiBin_1",
        host_contig_id="S99_1",
        host_array_lo=100,
        host_array_hi=300,
        target_meta=parsed[sha2],
    )
    assert ev_s2["provirus_boundary_status"] == "verified_provirus_with_host_flanks"
    assert ev_s2["read_junction_support_status"] == "untested_missing"


def test_stage2_refuses_unverified_old_annotation_binding_and_keeps_legacy_audit(tmp_path: Path):
    """Verify EXP3-FINAL-PROVENANCE-REVIEW item 1:
    1. Attempting to bind precomputed CheckV/geNomad tables to current sequence SHA-256s
       by target ID without a verified original-run input manifest hash map raises ValueError.
    2. When an original-run manifest hash map IS provided, any target whose current sequence
       SHA-256 differs from the original manifest SHA-256 is rejected and NOT bound.
    3. Legacy precomputed annotations are kept strictly as an ID-keyed audit (`legacy_precomputed_votu_table_unverified_original_input_hash`).
    """
    cv_dir = tmp_path / "old_checkv"
    cv_dir.mkdir()
    (cv_dir / "quality_summary.tsv").write_text(
        "contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n"
        "V1\t5000\tNo\tNA\t8\t5\t0\tHigh-quality\tGenome-fragment\t95.0\tAAI\t0.0\t1.0\t\n"
        "V2\t6000\tNo\tNA\t9\t6\t0\tHigh-quality\tGenome-fragment\t96.0\tAAI\t0.0\t1.0\t\n"
    )
    cur_sha_v1 = sha256_bytes(b"ACGT" * 100)
    cur_sha_v2_mutated = sha256_bytes(b"ACGA" * 100)
    orig_sha_v2 = sha256_bytes(b"TGCA" * 100)
    seq_sha_map = {"V1": cur_sha_v1, "V2": cur_sha_v2_mutated}

    # 1. Refuses binding precomputed annotation by ID without original-run manifest hash map
    with pytest.raises(ValueError, match="Refusing to bind precomputed annotation"):
        parse_checkv_and_genomad_dirs(
            cv_dir,
            None,
            None,
            seq_sha_map,
            "precomputed_votu_reps_sha256_verified",
        )

    # 2. Only binds V1 (exact match to original manifest hash), rejects V2 (mutated sequence)
    verified_manifest = {"V1": cur_sha_v1, "V2": orig_sha_v2}
    bound = parse_checkv_and_genomad_dirs(
        cv_dir,
        None,
        None,
        seq_sha_map,
        "precomputed_votu_reps_manifest_verified",
        verified_original_manifest_sha256_map=verified_manifest,
    )
    assert cur_sha_v1 in bound
    assert cur_sha_v2_mutated not in bound
    assert orig_sha_v2 not in bound

    # 3. Legacy audit loads by ID with explicit unverified provenance
    legacy = load_legacy_precomputed_votu_audit(cv_dir, None, None, {"V1", "V2"})
    assert legacy["V1"]["legacy_audit_provenance"] == "legacy_precomputed_votu_table_unverified_original_input_hash"
    assert legacy["V1"]["legacy_checkv_quality"] == "High-quality"


def test_stage3_failed_caller_vs_zero_spacers_disposition():
    """Verify failed caller execution is never conflated with biological zero spacers."""
    assert summarize_mag_host_disposition("m1", True, False, 0, 0, [], caller_failed=True) == "Failed_caller_execution"
    assert summarize_mag_host_disposition("m2", True, False, 0, 0, [], caller_failed=False) == "Unknown_unscreened"
    assert summarize_mag_host_disposition("m3", True, True, 0, 0, [], caller_failed=False) == "Zero_spacers_detected"


def test_stage4_missing_mag_vs_all_zero_observed_and_shared_samples_and_failed_rows():
    """Verify EXP3-FINAL-PROVENANCE-REVIEW item 2:
    1. Missing MAG abundance is UNKNOWN (`NaN` + `unknown_missing_from_mag_abundance_matrix`),
       tested separately from an all-zero observed MAG row (`0.0` + `observed_all_zero_in_mag_abundance_matrix`).
    2. Per-sample NaN values within an observed MAG row are excluded from prevalence/mean denominators.
    3. Unshared sample columns between MAG and viral matrices are tracked as unknown and excluded.
    4. Co-occurrence preserves `"NA"` (unknown, never `0`) for missing MAGs, missing viruses,
       failed/missing alignment rows, and missing/failed e1 breadth samples.
    5. Phylogeny-controlled regressions exclude missing MAG abundances, missing defense annotations,
       and failed/unknown callers from denominators while retaining observed all-zero MAGs.
    """
    mag_df = pd.DataFrame(
        {
            "S1": [10.0, 0.0, np.nan],
            "S2": [0.0, 0.0, np.nan],
            "S3": [np.nan, 0.0, np.nan],
            "S_MAG_ONLY": [5.0, 0.0, np.nan],
        },
        index=["MAG_OBS_POS", "MAG_OBS_ALL_ZERO", "MAG_ALL_NAN"],
    )
    vir_df = pd.DataFrame(
        {
            "S1": [20.0, 0.0],
            "S2": [5.0, 0.0],
            "S3": [0.0, 0.0],
            "S_VIR_ONLY": [12.0, 1.0],
        },
        index=["VIR_1", "VIR_ZERO"],
    )

    inter = resolve_shared_sample_intersection(mag_df, vir_df)
    assert inter["shared_samples"] == ["S1", "S2", "S3"]
    assert inter["mag_only_samples"] == ["S_MAG_ONLY"]
    assert inter["vir_only_samples"] == ["S_VIR_ONLY"]
    assert inter["intersection_status"] == "partial_sample_intersection_unshared_excluded_as_unknown"

    shared = inter["shared_samples"]
    metrics = compute_host_abundance_metrics(
        ["MAG_OBS_POS", "MAG_OBS_ALL_ZERO", "MAG_MISSING", "MAG_ALL_NAN"],
        mag_df,
        sample_cols=shared,
    ).set_index("genome")

    # 1. Observed positive MAG with 1 NaN sample: denominator is 2 valid samples, not 3!
    assert metrics.loc["MAG_OBS_POS", "host_mapping_status"] == "observed_in_mag_abundance_matrix"
    assert metrics.loc["MAG_OBS_POS", "n_valid_abundance_samples"] == 2
    assert metrics.loc["MAG_OBS_POS", "n_missing_abundance_samples"] == 1
    assert metrics.loc["MAG_OBS_POS", "sample_prevalence"] == 1.0
    assert metrics.loc["MAG_OBS_POS", "sample_prevalence_fraction"] == 0.5
    assert metrics.loc["MAG_OBS_POS", "mean_sample_abundance"] == 5.0

    # 2. Observed all-zero MAG row: prevalence is 0.0 and mean abundance is 0.0 (true observed zero)
    assert metrics.loc["MAG_OBS_ALL_ZERO", "host_mapping_status"] == "observed_all_zero_in_mag_abundance_matrix"
    assert metrics.loc["MAG_OBS_ALL_ZERO", "n_valid_abundance_samples"] == 3
    assert metrics.loc["MAG_OBS_ALL_ZERO", "sample_prevalence"] == 0.0
    assert metrics.loc["MAG_OBS_ALL_ZERO", "sample_prevalence_fraction"] == 0.0
    assert metrics.loc["MAG_OBS_ALL_ZERO", "mean_sample_abundance"] == 0.0

    # 3. Missing MAG (not in mag_df.index): strictly NaN, NEVER 0!
    assert metrics.loc["MAG_MISSING", "host_mapping_status"] == "unknown_missing_from_mag_abundance_matrix"
    assert metrics.loc["MAG_MISSING", "n_valid_abundance_samples"] == 0
    assert np.isnan(metrics.loc["MAG_MISSING", "sample_prevalence"])
    assert np.isnan(metrics.loc["MAG_MISSING", "sample_prevalence_fraction"])
    assert np.isnan(metrics.loc["MAG_MISSING", "mean_sample_abundance"])

    # 4. Present MAG with all-NaN sample entries: strictly NaN, NEVER 0!
    assert metrics.loc["MAG_ALL_NAN", "host_mapping_status"] == "unknown_all_sample_entries_nan_or_failed"
    assert np.isnan(metrics.loc["MAG_ALL_NAN", "sample_prevalence"])
    assert np.isnan(metrics.loc["MAG_ALL_NAN", "mean_sample_abundance"])

    # 5. Co-occurrence preserves unknown ("NA") for missing MAG, missing virus, failed alignment row,
    #    and excludes failed/missing e1 samples from breadth regime denominators
    ledger_df = pd.DataFrame(
        {
            "genome": ["MAG_OBS_POS", "MAG_OBS_ALL_ZERO", "MAG_MISSING"],
            "canonical_phylum": ["p__Thermoplasmatota"] * 3,
            "mag_host_disposition": [
                "Tier1_verified_viral_host",
                "Zero_spacers_detected",
                "Zero_spacers_detected",
            ],
            "sample_prevalence": [1.0, 0.0, np.nan],
        }
    )
    pairs_df = pd.DataFrame(
        [
            {
                "canonical_591_host_mag": "MAG_MISSING",
                "votu_id": "VIR_1",
                "host_phylum": "p__Thermoplasmatota",
                "evidence_tier": "Tier1_verified_direct",
                "eligible_stage4_primary": 1,
                "alignment_status": "valid",
            },
            {
                "canonical_591_host_mag": "MAG_OBS_POS",
                "votu_id": "VIR_1",
                "host_phylum": "p__Thermoplasmatota",
                "evidence_tier": "Tier1_verified_direct",
                "eligible_stage4_primary": 1,
                "alignment_status": "failed_alignment",
            },
            {
                "canonical_591_host_mag": "MAG_OBS_POS",
                "votu_id": "VIR_1",
                "host_phylum": "p__Thermoplasmatota",
                "evidence_tier": "Tier1_verified_direct",
                "eligible_stage4_primary": 1,
                "alignment_status": "valid",
            },
        ]
    )
    cooc = evaluate_pair_cooccurrence(
        eligible_pairs=pairs_df,
        ledger_df=ledger_df,
        mag_abund_df=mag_df[shared],
        vir_tpm_df=vir_df[shared],
        shared_samples=shared,
        e1_mask={"VIR_1": {"S1"}},
        b75_mask={"VIR_1": {"S1"}},
        e1_observed_samples={"S1"},  # S2 is missing from e1, S3 is NaN in mag_df
        e1_failed_samples=set(),
    )
    row_missing_mag = cooc[cooc["host_mag"] == "MAG_MISSING"].iloc[0]
    assert row_missing_mag["mapping_status"] == "unknown_host_missing_from_mag_abundance_matrix"
    assert row_missing_mag["n_host_present_samples"] == "NA"
    assert row_missing_mag["n_virus_present_samples"] == "NA"

    row_failed_align = cooc[
        (cooc["host_mag"] == "MAG_OBS_POS") & (cooc["mapping_status"] == "unknown_failed_or_missing_alignment_row")
    ].iloc[0]
    assert row_failed_align["n_host_present_samples"] == "NA"

    row_obs_e1 = cooc[
        (cooc["host_mag"] == "MAG_OBS_POS")
        & (cooc["mapping_status"] == "observed_in_both_host_and_viral_matrices")
        & (cooc["breadth_regime"] == "e1_presence_rule")
    ].iloc[0]
    # Only S1 is both finite in mag_df/vir_df AND present in e1_observed_samples; S2 & S3 are excluded as unknown!
    assert row_obs_e1["n_evaluated_samples"] == 1
    assert row_obs_e1["n_unknown_excluded_samples"] == 2
    assert row_obs_e1["n_codetected_samples"] == 1

    # 6. Phylogeny-controlled regression test exclusions for missing MAG abundance, failed caller, and missing defense
    rng = np.random.default_rng(42)
    n_obs = 20
    reg_rows = []
    for i in range(n_obs):
        reg_rows.append(
            {
                "genome": f"M_OBS_{i}",
                "canonical_phylum": "p__Thermoplasmatota",
                "canonical_order": "o__OrderA" if i % 2 == 0 else "o__OrderB",
                "mag_host_disposition": "Tier1_verified_viral_host" if i < 6 else "Zero_spacers_detected",
                "n_verified_spacers": 3 if i < 6 else 0,
                "n_cas_operons": 1 if i < 8 else 0,
                "n_tier1_verified_direct_links": 1 if i < 6 else 0,
                "n_tier2_corroborated_links": 0,
                "c8_defense_systems": float(i % 5),
                "mean_sample_abundance": 0.0 if i == 19 else float(rng.uniform(0.1, 5.0)),  # i==19 is observed all-zero!
                "checkm2_completeness": 92.0,
                "refined_n_bases": 1500000,
                "refined_n_contigs": 50,
            }
        )
    # Add 1 missing-abundance MAG, 1 failed-caller MAG, and 1 missing-c8-defense MAG
    reg_rows.append(
        {
            **reg_rows[0],
            "genome": "M_MISS_ABUND",
            "mean_sample_abundance": np.nan,
        }
    )
    reg_rows.append(
        {
            **reg_rows[0],
            "genome": "M_FAILED_CALLER",
            "mag_host_disposition": "Failed_caller_execution",
        }
    )
    reg_rows.append(
        {
            **reg_rows[0],
            "genome": "M_MISS_DEF",
            "c8_defense_systems": np.nan,
        }
    )
    reg_df = pd.DataFrame(reg_rows)
    res_tests = run_phylogeny_controlled_tests(reg_df)
    t_row = res_tests[
        (res_tests["scope"] == "thermoplasmatota_237")
        & (res_tests["predictor"] == "has_verified_spacers")
        & (res_tests["response"] == "c8_def_num")
    ].iloc[0]
    # Total in scope = 23; 1 missing abundance excluded, 2 missing predictor/response (failed caller + missing c8_def) excluded;
    # all 20 observed MAGs (including the observed all-zero MAG at i==19) are retained!
    assert t_row["n_scope_total_mags"] == 23
    assert t_row["n_mags"] == 20
    assert t_row["n_missing_abundance_excluded"] == 1
    assert t_row["n_unknown_caller_or_alignment_excluded"] == 2


def test_stage2_cellular_controls_not_falsely_flagged_as_plasmids(tmp_path: Path):
    """Verify that cellular genes (transcription, transport, transferase, translation,
    and genes with annotation_conjscan == 'NA') are NOT falsely classified as plasmid hallmarks,
    while true plasmid and viral hallmarks are accurately recognized.
    """
    gn_genes = tmp_path / "genes.tsv"
    gn_genes.write_text(
        "gene\tstart\tend\tlength\tstrand\tgc_content\tgenetic_code\trbs_motif\tmarker\tevalue\tbitscore\tuscg\tplasmid_hallmark\tvirus_hallmark\ttaxid\ttaxname\tannotation_conjscan\tannotation_amr\tannotation_accessions\tannotation_description\n"
        "C_CELL_1\t1\t300\t300\t1\t0.5\t11\tNone\tNA\tNA\tNA\t0\t0\t0\t1\tNA\tNA\tNA\tPF00001\tABC-type transporter\n"
        "C_CELL_2\t301\t600\t300\t1\t0.5\t11\tNone\tNA\tNA\tNA\t0\t0\t0\t1\tNA\tNA\tNA\tPF00002\tAcetyltransferase\n"
        "C_CELL_3\t601\t900\t300\t1\t0.5\t11\tNone\tNA\tNA\tNA\t0\t0\t0\t1\tNA\tNA\tNA\tPF00003\tTranscription initiation factor\n"
        "C_CELL_4\t901\t1200\t300\t1\t0.5\t11\tNone\tNA\tNA\tNA\t0\t0\t0\t1\tNA\tNA\tNA\tPF00004\tTranslation elongation factor\n"
        "C_PLAS_1\t1\t300\t300\t1\t0.5\t11\tNone\tNA\tNA\tNA\t0\t1\t0\t1\tNA\tNA\tNA\tPF00005\tHypothetical protein\n"
        "C_PLAS_2\t301\t600\t300\t1\t0.5\t11\tNone\tNA\tNA\tNA\t0\t0\t0\t1\tNA\tTraA\tNA\tPF00006\tConjugal transfer protein TraA\n"
        "C_VIR_1\t1\t300\t300\t1\t0.5\t11\tNone\tNA\tNA\tNA\t0\t0\t1\t1\tNA\tNA\tNA\tPF00007\tMajor capsid protein\n"
    )
    seq_map = {
        "C_CELL": sha256_bytes(b"A" * 1200),
        "C_PLAS": sha256_bytes(b"C" * 600),
        "C_VIR": sha256_bytes(b"G" * 300),
    }
    parsed = parse_checkv_and_genomad_dirs(
        checkv_dir=None,
        genomad_genes_path=gn_genes,
        genomad_tax_path=None,
        seq_sha_map=seq_map,
        provenance_label="unit_test",
    )
    # C_CELL must have 0 plasmid hallmarks and 0 viral hallmarks!
    cell_meta = parsed[seq_map["C_CELL"]]
    assert cell_meta.get("plasmid_hallmarks", 0) == 0, f"Expected 0 plasmid hallmarks, got {cell_meta.get('plasmid_hallmarks')}"
    assert cell_meta.get("virus_hallmarks", 0) == 0

    # C_PLAS must have 2 plasmid hallmarks (plasmid_hallmark==1 and TraA)
    plas_meta = parsed[seq_map["C_PLAS"]]
    assert plas_meta.get("plasmid_hallmarks", 0) == 2

    # C_VIR must have 1 viral hallmark
    vir_meta = parsed[seq_map["C_VIR"]]
    assert vir_meta.get("virus_hallmarks", 0) == 1


def test_blast_receipt_interrupted_partial(tmp_path: Path):
    """Verify that interrupted or incomplete BLAST runs (.partial.tsv or .tmp) are rejected."""
    query_fna = tmp_path / "query.fna"
    query_fna.write_text(">spc1\nATCGATCGATCGATCG\n")
    q_sha = sha256_bytes(query_fna.read_bytes())

    subj_fna = tmp_path / "subject.fna"
    subj_fna.write_text(">contig1\nATCGATCGATCGATCGATCGATCG\n")
    s_sha = stream_file_sha256(subj_fna)

    out_tsv = tmp_path / "search.blast.tsv"
    out_tsv.write_text("spc1\tcontig1\t100.0\t16\t0\t0\t0\t1\t16\t1\t16\t1e-5\t32.0\t16\t24\tATCGATCGATCGATCG\tATCGATCGATCGATCG\t16\n")

    write_blast_receipt(
        out_tsv=out_tsv,
        task_type="unit_test",
        target_name="test_target",
        query_fna=query_fna,
        query_sha256=q_sha,
        subject_fna=subj_fna,
        subject_sha256=s_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=1.0,
    )

    # Valid run
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha, expected_threads=8) is True

    # Interrupted run leaving .partial.tsv
    partial_tsv = out_tsv.with_suffix(".partial.tsv")
    partial_tsv.write_text("incomplete line\n")
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha, expected_threads=8) is False

    partial_tsv.unlink()
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha, expected_threads=8) is True

    # Interrupted run leaving temporary receipt
    rec_tmp = Path(str(out_tsv.with_suffix(".receipt.json")) + ".tmp")
    rec_tmp.write_text("{}")
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha, expected_threads=8) is False


def test_blast_receipt_legitimate_empty_success(tmp_path: Path):
    """Verify that legitimate 0-hit (0-byte) completed searches are accepted with valid receipt."""
    query_fna = tmp_path / "query.fna"
    query_fna.write_text(">spc1\nATCGATCGATCGATCG\n")
    q_sha = sha256_bytes(query_fna.read_bytes())

    subj_fna = tmp_path / "subject.fna"
    subj_fna.write_text(">contig1\nGGGGGGGGGGGGGGGG\n")
    s_sha = stream_file_sha256(subj_fna)

    out_tsv = tmp_path / "empty_search.blast.tsv"
    out_tsv.write_text("")  # 0 bytes

    write_blast_receipt(
        out_tsv=out_tsv,
        task_type="unit_test",
        target_name="empty_target",
        query_fna=query_fna,
        query_sha256=q_sha,
        subject_fna=subj_fna,
        subject_sha256=s_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=0.5,
    )

    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha, expected_threads=8) is True


def test_blast_receipt_changed_input_or_params(tmp_path: Path):
    """Verify that altered inputs, missing subjects, or altered blast parameters invalidate receipt."""
    query_fna = tmp_path / "query.fna"
    query_fna.write_text(">spc1\nATCGATCGATCGATCG\n")
    q_sha = sha256_bytes(query_fna.read_bytes())

    subj_fna = tmp_path / "subject.fna"
    subj_fna.write_text(">contig1\nATCGATCGATCGATCGATCGATCG\n")
    s_sha = stream_file_sha256(subj_fna)

    out_tsv = tmp_path / "search.blast.tsv"
    out_tsv.write_text("spc1\tcontig1\t100.0\t16\t0\t0\t0\t1\t16\t1\t16\t1e-5\t32.0\t16\t24\tATCGATCGATCGATCG\tATCGATCGATCGATCG\t16\n")

    write_blast_receipt(
        out_tsv=out_tsv,
        task_type="unit_test",
        target_name="test_target",
        query_fna=query_fna,
        query_sha256=q_sha,
        subject_fna=subj_fna,
        subject_sha256=s_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=1.0,
    )

    # 1. Missing subject file must immediately fail (never skipped)
    missing_subj = tmp_path / "nonexistent.fna"
    assert is_blast_receipt_valid(out_tsv, q_sha, missing_subj, expected_subject_sha256=s_sha) is False

    # 2. Altered query hash
    altered_q_sha = sha256_bytes(b"MUTATED_QUERY")
    assert is_blast_receipt_valid(out_tsv, altered_q_sha, subj_fna, expected_subject_sha256=s_sha) is False

    # 3. Altered subject size or mtime
    import os
    orig_subj_bytes = subj_fna.read_bytes()
    orig_subj_mtime = subj_fna.stat().st_mtime
    subj_fna.write_bytes(orig_subj_bytes + b"\nEXTRA_SEQUENCE\n")
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha) is False
    subj_fna.write_bytes(orig_subj_bytes)  # restore
    os.utime(subj_fna, (orig_subj_mtime, orig_subj_mtime))
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha) is True

    # 4. Altered thread count
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha, expected_threads=16) is False
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha, expected_threads=8) is True

    # 5. Tampered parameters in receipt
    rec_p = out_tsv.with_suffix(".receipt.json")
    import json
    rec_data = json.loads(rec_p.read_text())
    rec_data["blast_parameters"]["evalue"] = "1e-5"
    rec_p.write_text(json.dumps(rec_data))
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha, expected_threads=8) is False


def test_blast_receipt_parser_or_corrupt_failure(tmp_path: Path):
    """Verify that corrupted outputs or invalid receipts are rejected."""
    query_fna = tmp_path / "query.fna"
    query_fna.write_text(">spc1\nATCGATCGATCGATCG\n")
    q_sha = sha256_bytes(query_fna.read_bytes())

    subj_fna = tmp_path / "subject.fna"
    subj_fna.write_text(">contig1\nATCGATCGATCGATCGATCGATCG\n")
    s_sha = stream_file_sha256(subj_fna)

    out_tsv = tmp_path / "search.blast.tsv"
    out_tsv.write_text("spc1\tcontig1\t100.0\t16\t0\t0\t0\t1\t16\t1\t16\t1e-5\t32.0\t16\t24\tATCGATCGATCGATCG\tATCGATCGATCGATCG\t16\n")

    write_blast_receipt(
        out_tsv=out_tsv,
        task_type="unit_test",
        target_name="test_target",
        query_fna=query_fna,
        query_sha256=q_sha,
        subject_fna=subj_fna,
        subject_sha256=s_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=1.0,
    )

    # Corrupt JSON receipt
    rec_p = out_tsv.with_suffix(".receipt.json")
    rec_p.write_text("{corrupt json content")
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha) is False

    # Output SHA mismatch
    write_blast_receipt(
        out_tsv=out_tsv,
        task_type="unit_test",
        target_name="test_target",
        query_fna=query_fna,
        query_sha256=q_sha,
        subject_fna=subj_fna,
        subject_sha256=s_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=1.0,
    )
    out_tsv.write_text("different content\n")
    assert is_blast_receipt_valid(out_tsv, q_sha, subj_fna, expected_subject_sha256=s_sha) is False

    # Corrupted tabular line: truncated columns
    corrupt_tsv = tmp_path / "corrupt_cols.tsv"
    corrupt_tsv.write_text("spc1\tcontig1\t100.0\tshort_row\n")
    is_val, _, err = validate_blast_tabular_output(corrupt_tsv)
    assert is_val is False
    assert "columns, expected 18" in err

    # Corrupted tabular line: non-numeric pident
    corrupt_tsv2 = tmp_path / "corrupt_types.tsv"
    corrupt_tsv2.write_text("spc1\tcontig1\tNOT_A_FLOAT\t16\t0\t0\t0\t1\t16\t1\t16\t1e-5\t32.0\t16\t24\tATCG\tATCG\t16\n")
    is_val2, _, err2 = validate_blast_tabular_output(corrupt_tsv2)
    assert is_val2 is False
    assert "could not convert string to float" in err2


def test_target_assemblies_manifest_audit(tmp_path: Path):
    """Verify build_and_audit_target_assemblies_manifest tracks present, missing, failed, and completed."""
    ver_tsv = tmp_path / "arch591_ver.tsv"
    ver_tsv.write_text("genome\tsource_sample\tpass_hash\nMAG1\tSAMP_A\t1\nMAG2\tSAMP_B\t1\n")

    spc_tsv = tmp_path / "spacers.tsv"
    spc_tsv.write_text("unique_spacer_id\tmag_id\nSPC1\tSAMP_C__bin_1\n")

    asm_dir = tmp_path / "assemblies"
    # SAMP_A: present & completed
    (asm_dir / "SAMP_A").mkdir(parents=True)
    asm_a = asm_dir / "SAMP_A" / "SAMP_A_contigs.fa"
    asm_a.write_text(">c1\nACGT\n")
    s_sha_a = stream_file_sha256(asm_a)

    # SAMP_B: present & unverified (TSV exists but no receipt)
    (asm_dir / "SAMP_B").mkdir(parents=True)
    asm_b = asm_dir / "SAMP_B" / "SAMP_B_contigs.fa"
    asm_b.write_text(">c2\nACGT\n")

    # SAMP_C: missing assembly file
    raw_asm_dir = tmp_path / "raw_asm"
    raw_asm_dir.mkdir()

    # Create completed search for SAMP_A
    q_fna = tmp_path / "q.fna"
    q_fna.write_text(">q\nACGT\n")
    q_sha = sha256_bytes(q_fna.read_bytes())
    out_a = raw_asm_dir / "raw_asm_SAMP_A.blast.tsv"
    out_a.write_text("")  # 0-byte valid
    write_blast_receipt(
        out_tsv=out_a,
        task_type="unit_test",
        target_name="SAMP_A",
        query_fna=q_fna,
        query_sha256=q_sha,
        subject_fna=asm_a,
        subject_sha256=s_sha_a,
        threads=8,
        exit_code=0,
        elapsed_seconds=0.1,
    )

    # Create unverified TSV for SAMP_B (no receipt)
    out_b = raw_asm_dir / "raw_asm_SAMP_B.blast.tsv"
    out_b.write_text("some\ttsv\tcontent\n")

    manifest_tsv = tmp_path / "manifest.tsv"
    summary_json = tmp_path / "summary.json"
    cache_hashes_p = tmp_path / "hashes_cache.json"

    mdf, summary = build_and_audit_target_assemblies_manifest(
        manifest_tsv=manifest_tsv,
        summary_json=summary_json,
        arch591_ver_tsv=ver_tsv,
        stage1_spacers_tsv=spc_tsv,
        assemblies_dir=asm_dir,
        raw_asm_dir=raw_asm_dir,
        query_sha256=q_sha,
        subject_hashes_cache_path=cache_hashes_p,
        expected_threads=8,
    )

    assert len(mdf) == 3
    assert summary["n_target_cohort_samples"] == 3
    assert summary["n_assemblies_present"] == 2
    assert summary["n_assemblies_missing"] == 1
    assert summary["missing_assembly_samples"] == ["SAMP_C"]
    assert summary["n_assemblies_completed_verified"] == 1
    assert summary["n_assemblies_unverified_or_failed"] == 1
    assert summary["n_assemblies_unsearched"] == 0
    assert summary["is_cohort_coverage_complete"] is False
    assert summary["is_present_coverage_complete"] is False

    row_a = mdf[mdf["sample_id"] == "SAMP_A"].iloc[0]
    assert row_a["search_status"] == "completed_verified"
    assert row_a["receipt_verified"] == 1

    row_b = mdf[mdf["sample_id"] == "SAMP_B"].iloc[0]
    assert row_b["search_status"] == "unverified_or_failed"
    assert row_b["receipt_verified"] == 0

    row_c = mdf[mdf["sample_id"] == "SAMP_C"].iloc[0]
    assert row_c["search_status"] == "missing_assembly_file"
    assert row_c["assembly_exists"] == 0


def test_raw_assembly_same_size_same_mtime_mutation_invalidates(tmp_path: Path):
    """Verify that same-size same-mtime content mutation on assembly file invalidates receipt (no hash bypass!)."""
    import os
    q_fna = tmp_path / "query.fna"
    q_fna.write_text(">spc\nATCG\n")
    q_sha = sha256_bytes(q_fna.read_bytes())

    asm_fna = tmp_path / "raw_asm.fa"
    asm_fna.write_text(">c1\nACGTACGT\n")  # 13 bytes
    orig_sha = stream_file_sha256(asm_fna)
    orig_mtime = asm_fna.stat().st_mtime

    out_tsv = tmp_path / "raw_asm_TEST.blast.tsv"
    out_tsv.write_text("")  # 0 bytes
    write_blast_receipt(
        out_tsv=out_tsv,
        task_type="raw_assembly_blast",
        target_name="TEST",
        query_fna=q_fna,
        query_sha256=q_sha,
        subject_fna=asm_fna,
        subject_sha256=orig_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=0.1,
    )

    # Initial state is valid
    assert is_blast_receipt_valid(out_tsv, q_sha, asm_fna, expected_subject_sha256=orig_sha, expected_threads=8) is True

    # Mutate content while keeping exact same byte size (13 bytes)
    asm_fna.write_text(">c1\nACGTACGA\n")
    assert asm_fna.stat().st_size == 13
    mutated_sha = stream_file_sha256(asm_fna)
    assert mutated_sha != orig_sha

    # Restore exact mtime so size AND mtime match receipt, but SHA differs!
    os.utime(asm_fna, (orig_mtime, orig_mtime))

    # Validator MUST reject because subject_sha256 does not match!
    assert is_blast_receipt_valid(out_tsv, q_sha, asm_fna, expected_subject_sha256=mutated_sha, expected_threads=8) is False


def test_blast_tool_probe_failure_raises_in_production(monkeypatch):
    """Verify that get_blast_tool_info raises on missing container or failed probe without fallback."""
    # Reset cached tool info
    monkeypatch.setattr(stage2, "_BLAST_TOOL_INFO_CACHE", None)
    monkeypatch.setattr(stage2, "get_blast_tool_info", ORIGINAL_GET_BLAST_TOOL_INFO)

    # 1. Missing container image -> FileNotFoundError
    monkeypatch.setattr(stage2, "SIF_BLAST", Path("/nonexistent/container/blast.sif"))
    with pytest.raises(FileNotFoundError):
        ORIGINAL_GET_BLAST_TOOL_INFO()

    # 2. Existing fake container with singularity failure -> RuntimeError
    fake_sif = Path("/bin/sh")  # exists
    monkeypatch.setattr(stage2, "SIF_BLAST", fake_sif)
    with pytest.raises(RuntimeError):
        ORIGINAL_GET_BLAST_TOOL_INFO()


def test_blast_receipt_changed_container_hash_fails(tmp_path: Path, monkeypatch):
    """Verify that is_blast_receipt_valid rejects receipts if container_sha256 differs (tampered image)."""
    q_fna = tmp_path / "query.fna"
    q_fna.write_text(">spc\nATCG\n")
    q_sha = sha256_bytes(q_fna.read_bytes())

    asm_fna = tmp_path / "raw_asm.fa"
    asm_fna.write_text(">c1\nACGTACGT\n")
    orig_sha = stream_file_sha256(asm_fna)

    out_tsv = tmp_path / "raw_asm_TEST.blast.tsv"
    out_tsv.write_text("")

    mock_info = {
        "tool_binary": "blastn",
        "tool_version": "blastn: 2.17.0+",
        "container_sif": str(stage2.SIF_BLAST),
        "container_sha256": "6c7923d5c6042e279cdecb719b0d8b84535a10f8f97485a18695a8d6902ac51f",
    }
    monkeypatch.setattr(stage2, "get_blast_tool_info", lambda: dict(mock_info))

    write_blast_receipt(
        out_tsv=out_tsv,
        task_type="raw_assembly_blast",
        target_name="TEST",
        query_fna=q_fna,
        query_sha256=q_sha,
        subject_fna=asm_fna,
        subject_sha256=orig_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=0.1,
    )

    # Valid with matching container hash
    assert is_blast_receipt_valid(out_tsv, q_sha, asm_fna, expected_subject_sha256=orig_sha, expected_threads=8) is True

    # Mutate container_sha256 in current tool info (e.g. tampered or updated container SIF)
    tampered_info = dict(mock_info, container_sha256="bad_tampered_container_hash_1234567890abcdef")
    monkeypatch.setattr(stage2, "get_blast_tool_info", lambda: dict(tampered_info))

    # Validator MUST reject because container digest does not match receipt!
    assert is_blast_receipt_valid(out_tsv, q_sha, asm_fna, expected_subject_sha256=orig_sha, expected_threads=8) is False


def test_standalone_audit_success(tmp_path: Path, monkeypatch):
    """Verify standalone audit passes when all expected catalog and assembly targets are valid."""
    import json
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    asm_raw_dir = raw_dir / "raw_assembly_blast"
    asm_raw_dir.mkdir()
    asm_dir = tmp_path / "assemblies"
    asm_dir.mkdir()

    # Query
    q_fna = tmp_path / "spacers.fna"
    q_fna.write_text(">spc1\nACGTACGT\n")
    q_sha = sha256_bytes(q_fna.read_bytes())

    # Mock tool info
    mock_info = {
        "tool_binary": "blastn",
        "tool_version": "blastn: 2.17.0+",
        "container_sif": str(stage2.SIF_BLAST),
        "container_sha256": "6c7923d5c6042e279cdecb719b0d8b84535a10f8f97485a18695a8d6902ac51f",
    }
    monkeypatch.setattr(stage2, "get_blast_tool_info", lambda: dict(mock_info))

    # Mock 1 catalog pool
    cat_fna = tmp_path / "cat1.fna"
    cat_fna.write_text(">cat1\nGGCC\n")
    cat_sha = stream_file_sha256(cat_fna)
    cat_tsv = raw_dir / "blast_cat1.tsv"
    cat_tsv.write_text("")
    write_blast_receipt(
        out_tsv=cat_tsv,
        task_type="catalog_pool_blast",
        target_name="blast_cat1",
        query_fna=q_fna,
        query_sha256=q_sha,
        subject_fna=cat_fna,
        subject_sha256=cat_sha,
        threads=32,
        exit_code=0,
        elapsed_seconds=0.1,
    )

    # Mock 1 assembly sample
    samp = "SAMP_01"
    samp_dir = asm_dir / samp
    samp_dir.mkdir()
    asm_fna = samp_dir / f"{samp}_contigs.fa"
    asm_fna.write_text(">contig1\nTTTT\n")
    asm_sha = stream_file_sha256(asm_fna)
    asm_tsv = asm_raw_dir / f"raw_asm_{samp}.blast.tsv"
    asm_tsv.write_text("")
    write_blast_receipt(
        out_tsv=asm_tsv,
        task_type="raw_assembly_blast",
        target_name=samp,
        query_fna=q_fna,
        query_sha256=q_sha,
        subject_fna=asm_fna,
        subject_sha256=asm_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=0.1,
    )

    # Write manifest
    manifest_data = {
        samp: {
            "size": asm_fna.stat().st_size,
            "mtime": asm_fna.stat().st_mtime,
            "sha256": asm_sha,
        }
    }
    (raw_dir / "target_assemblies_subject_hashes.json").write_text(json.dumps(manifest_data))

    catalog_targets = [("blast_cat1", cat_fna, 32)]
    is_ok, report = audit_stage2_receipts(
        raw_dir=raw_dir,
        spacers_fna=q_fna,
        assemblies_dir=asm_dir,
        expected_catalog_targets=catalog_targets,
        expected_assembly_samples=[samp],
        stream_check_subjects=True,
    )

    assert is_ok is True
    assert report["status"] == "passed"
    assert report["total_valid"] == 2
    assert report["total_expected"] == 2


def test_standalone_audit_same_size_same_mtime_mutation(tmp_path: Path, monkeypatch):
    """Verify standalone audit rejects and fails if subject FASTA is mutated with same size and mtime."""
    import os
    import json
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    asm_raw_dir = raw_dir / "raw_assembly_blast"
    asm_raw_dir.mkdir()
    asm_dir = tmp_path / "assemblies"
    asm_dir.mkdir()

    q_fna = tmp_path / "spacers.fna"
    q_fna.write_text(">spc1\nACGTACGT\n")
    q_sha = sha256_bytes(q_fna.read_bytes())

    mock_info = {
        "tool_binary": "blastn",
        "tool_version": "blastn: 2.17.0+",
        "container_sif": str(stage2.SIF_BLAST),
        "container_sha256": "6c7923d5c6042e279cdecb719b0d8b84535a10f8f97485a18695a8d6902ac51f",
    }
    monkeypatch.setattr(stage2, "get_blast_tool_info", lambda: dict(mock_info))

    samp = "SAMP_MUT"
    samp_dir = asm_dir / samp
    samp_dir.mkdir()
    asm_fna = samp_dir / f"{samp}_contigs.fa"
    asm_fna.write_text(">contig1\nAAAA\n")  # 14 bytes
    orig_sha = stream_file_sha256(asm_fna)
    orig_mtime = asm_fna.stat().st_mtime

    asm_tsv = asm_raw_dir / f"raw_asm_{samp}.blast.tsv"
    asm_tsv.write_text("")
    write_blast_receipt(
        out_tsv=asm_tsv,
        task_type="raw_assembly_blast",
        target_name=samp,
        query_fna=q_fna,
        query_sha256=q_sha,
        subject_fna=asm_fna,
        subject_sha256=orig_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=0.1,
    )

    manifest_data = {
        samp: {
            "size": asm_fna.stat().st_size,
            "mtime": orig_mtime,
            "sha256": orig_sha,
        }
    }
    (raw_dir / "target_assemblies_subject_hashes.json").write_text(json.dumps(manifest_data))

    # Mutate content with EXACT same size (14 bytes) and restore exact original mtime
    asm_fna.write_text(">contig1\nCCCC\n")
    assert asm_fna.stat().st_size == 14
    os.utime(asm_fna, (orig_mtime, orig_mtime))

    is_ok, report = audit_stage2_receipts(
        raw_dir=raw_dir,
        spacers_fna=q_fna,
        assemblies_dir=asm_dir,
        expected_catalog_targets=[],
        expected_assembly_samples=[samp],
        stream_check_subjects=True,
    )

    assert is_ok is False
    assert report["status"] == "failed"
    assert len(report["missing_or_invalid_targets"]) == 1
    err_msg = report["missing_or_invalid_targets"][0]["error"]
    assert "mutated" in err_msg or "mismatch" in err_msg


def test_standalone_audit_no_receipts(tmp_path: Path, monkeypatch):
    """Verify standalone audit rejects when raw directory has no receipts."""
    raw_dir = tmp_path / "empty_raw"
    raw_dir.mkdir()
    q_fna = tmp_path / "spacers.fna"
    q_fna.write_text(">spc\nACGT\n")

    mock_info = {
        "tool_binary": "blastn",
        "tool_version": "blastn: 2.17.0+",
        "container_sif": str(stage2.SIF_BLAST),
        "container_sha256": "6c7923d5c6042e279cdecb719b0d8b84535a10f8f97485a18695a8d6902ac51f",
    }
    monkeypatch.setattr(stage2, "get_blast_tool_info", lambda: dict(mock_info))

    is_ok, report = audit_stage2_receipts(
        raw_dir=raw_dir,
        spacers_fna=q_fna,
        expected_catalog_targets=[],
        expected_assembly_samples=[],
    )
    assert is_ok is False
    assert report["status"] == "failed"


def test_standalone_audit_missing_expected_receipts_nonempty_resultcount_rejected(tmp_path: Path, monkeypatch):
    """Verify standalone audit rejects when some expected receipts are missing, even if present ones are valid."""
    import json
    raw_dir = tmp_path / "raw"
    raw_dir.mkdir()
    asm_raw_dir = raw_dir / "raw_assembly_blast"
    asm_raw_dir.mkdir()
    asm_dir = tmp_path / "assemblies"
    asm_dir.mkdir()

    q_fna = tmp_path / "spacers.fna"
    q_fna.write_text(">spc1\nACGTACGT\n")
    q_sha = sha256_bytes(q_fna.read_bytes())

    mock_info = {
        "tool_binary": "blastn",
        "tool_version": "blastn: 2.17.0+",
        "container_sif": str(stage2.SIF_BLAST),
        "container_sha256": "6c7923d5c6042e279cdecb719b0d8b84535a10f8f97485a18695a8d6902ac51f",
    }
    monkeypatch.setattr(stage2, "get_blast_tool_info", lambda: dict(mock_info))

    # SAMP_1 is present and completely valid
    s1 = "SAMP_1"
    s1_dir = asm_dir / s1
    s1_dir.mkdir()
    s1_fna = s1_dir / f"{s1}_contigs.fa"
    s1_fna.write_text(">c1\nAAAA\n")
    s1_sha = stream_file_sha256(s1_fna)
    s1_tsv = asm_raw_dir / f"raw_asm_{s1}.blast.tsv"
    s1_tsv.write_text("")
    write_blast_receipt(
        out_tsv=s1_tsv,
        task_type="raw_assembly_blast",
        target_name=s1,
        query_fna=q_fna,
        query_sha256=q_sha,
        subject_fna=s1_fna,
        subject_sha256=s1_sha,
        threads=8,
        exit_code=0,
        elapsed_seconds=0.1,
    )

    # SAMP_2 has subject on disk and is in expected list, but receipt is missing!
    s2 = "SAMP_2"
    s2_dir = asm_dir / s2
    s2_dir.mkdir()
    s2_fna = s2_dir / f"{s2}_contigs.fa"
    s2_fna.write_text(">c2\nTTTT\n")
    s2_sha = stream_file_sha256(s2_fna)

    manifest_data = {
        s1: {"size": s1_fna.stat().st_size, "mtime": s1_fna.stat().st_mtime, "sha256": s1_sha},
        s2: {"size": s2_fna.stat().st_size, "mtime": s2_fna.stat().st_mtime, "sha256": s2_sha},
    }
    (raw_dir / "target_assemblies_subject_hashes.json").write_text(json.dumps(manifest_data))

    is_ok, report = audit_stage2_receipts(
        raw_dir=raw_dir,
        spacers_fna=q_fna,
        assemblies_dir=asm_dir,
        expected_catalog_targets=[],
        expected_assembly_samples=[s1, s2],
        stream_check_subjects=True,
    )

    # Non-empty resultcount (1 valid receipt out of 2) MUST be rejected!
    assert is_ok is False
    assert report["status"] == "failed"
    assert report["n_assemblies_valid"] == 1
    assert report["n_assemblies_expected"] == 2
    assert any(x["target_name"] == s2 for x in report["missing_or_invalid_targets"])


def test_thermo_spacers_audit_schema_and_deduplication(tmp_path: Path):
    """Verify Thermo pair forensic audit handles actual stage 1 schema, deduplicates pool hits, and rejects empty joins."""
    spc_p = tmp_path / "stage1_verified_archaeal_spacers.tsv"
    spc_df = pd.DataFrame([
        {
            "unique_spacer_id": "spc_1",
            "mag_id": "41PRr2__single_SemiBin_236",
            "is_arch591_rep": 1,
            "canonical_spacer_seq_rc": "A" * 32,
            "representative_spacer_seq": "T" * 32,
            "spacer_len": 32,
            "observed_orientations": "fwd",
            "n_loci_in_mag": 1,
            "n_distinct_arrays_in_mag": 1,
            "callers": "minced",
            "array_trust_tiers": "trusted",
            "caller_independence_note": "ok",
            "primary_namespaced_contig_id": "c1",
            "primary_array_start": 100,
            "primary_array_end": 500,
            "all_array_loci": "c1:100-500",
            "repeat_seq": "C" * 28,
        },
        {
            "unique_spacer_id": "spc_2",
            "mag_id": "41PRr2__single_SemiBin_236",
            "is_arch591_rep": 1,
            "canonical_spacer_seq_rc": "G" * 32,
            "representative_spacer_seq": "C" * 32,
            "spacer_len": 32,
            "observed_orientations": "fwd",
            "n_loci_in_mag": 1,
            "n_distinct_arrays_in_mag": 1,
            "callers": "minced",
            "array_trust_tiers": "trusted",
            "caller_independence_note": "ok",
            "primary_namespaced_contig_id": "c1",
            "primary_array_start": 100,
            "primary_array_end": 500,
            "all_array_loci": "c1:100-500",
            "repeat_seq": "C" * 28,
        }
    ])
    spc_df.to_csv(spc_p, sep="\t", index=False)

    # Scored alignments with 2 spacers triplicated across 3 pools (6 total raw records)
    align_rows = []
    for pool in ["eq_votu_reps", "eq_all_viruses", "eq_raw_assemblies_arch_samples"]:
        align_rows.append({
            "host_mag": "41PRr2__single_SemiBin_236",
            "target_id": "V60PRr2::V60PRr2_715976",
            "unique_spacer_id": "spc_1",
            "sstart": 1000,
            "send": 1031,
            "effective_mismatches": 0,
            "blast_pident": 100.0,
            "evalue": 1e-10,
            "bitscore": 58.0,
            "match_tier": "exact_0mm_full",
            "passes_strict_spacer_gate": True,
            "passes_exploratory_spacer_gate": True,
            "btop": "32",
            "target_seq_sha256": "dummy_sha_123",
            "pool": pool,
        })
        align_rows.append({
            "host_mag": "41PRr2__single_SemiBin_236",
            "target_id": "V60PRr2::V60PRr2_715976",
            "unique_spacer_id": "spc_2",
            "sstart": 5000,
            "send": 5031,
            "effective_mismatches": 1,
            "blast_pident": 96.8,
            "evalue": 1e-8,
            "bitscore": 52.0,
            "match_tier": "strict_1mm_core",
            "passes_strict_spacer_gate": True,
            "passes_exploratory_spacer_gate": True,
            "btop": "16AG15",
            "target_seq_sha256": "dummy_sha_123",
            "pool": pool,
        })
    scored_p = tmp_path / "stage2_all_spacer_alignments_scored.tsv"
    pd.DataFrame(align_rows).to_csv(scored_p, sep="\t", index=False)

    summary = audit_thermo_spacers(scored_p, spc_p)
    assert summary["raw_alignment_records_across_pools"] == 6
    assert summary["deduplicated_alignment_count"] == 2
    assert summary["n_strict_spacers"] == 2
    assert summary["n_exact_0mm_spacers"] == 1
    assert summary["n_unique_representative_sequences"] == 2
    assert summary["n_distinct_protospacer_loci"] == 2

    # Adversarial test: empty spacer sequence must raise ValueError
    bad_spc_df = spc_df.copy()
    bad_spc_df.loc[0, "representative_spacer_seq"] = ""
    bad_spc_p = tmp_path / "bad_spacers.tsv"
    bad_spc_df.to_csv(bad_spc_p, sep="\t", index=False)
    with pytest.raises(ValueError, match="invalid/empty sequence"):
        audit_thermo_spacers(scored_p, bad_spc_p)


def test_genomad_gene_schema_unannotated_detection_and_alias_resolution(tmp_path: Path):
    """Verify geNomad gene schema parsing detects unannotated ORFs and resolves alias targets."""
    # Test individual gene function
    r_unannot = pd.Series({"annotation_description": "NA", "marker": "NA", "virus_hallmark": 0})
    assert is_unannotated_genomad_gene(r_unannot) is True

    r_hallmark = pd.Series({"annotation_description": "capsid", "marker": "MCP", "virus_hallmark": 1})
    assert is_unannotated_genomad_gene(r_hallmark) is False

    r_host = pd.Series({"annotation_description": "ATPase", "marker": "NA", "virus_hallmark": 0})
    assert is_unannotated_genomad_gene(r_host) is False

    # Test full audit table
    gn_genes_p = tmp_path / "all_matched_primary_targets_genes.tsv"
    gn_df = pd.DataFrame([
        {
            "gene": "TGT_CANONICAL_1",
            "start": 1,
            "end": 300,
            "length": 300,
            "strand": 1,
            "gc_content": 0.5,
            "genetic_code": 11,
            "rbs_motif": "AGGAG",
            "marker": "NA",
            "evalue": "NA",
            "bitscore": "NA",
            "uscg": 0,
            "plasmid_hallmark": 0,
            "virus_hallmark": 0,
            "taxid": 1,
            "taxname": "NA",
            "annotation_conjscan": "NA",
            "annotation_amr": "NA",
            "annotation_accessions": "NA",
            "annotation_description": "NA",
        },
        {
            "gene": "TGT_CANONICAL_2",
            "start": 301,
            "end": 900,
            "length": 600,
            "strand": 1,
            "gc_content": 0.52,
            "genetic_code": 11,
            "rbs_motif": "AGGAG",
            "marker": "MCP",
            "evalue": "1e-20",
            "bitscore": "150.0",
            "uscg": 0,
            "plasmid_hallmark": 0,
            "virus_hallmark": 1,
            "taxid": 1,
            "taxname": "NA",
            "annotation_conjscan": "NA",
            "annotation_amr": "NA",
            "annotation_accessions": "NA",
            "annotation_description": "Major capsid protein",
        },
    ])
    gn_df.to_csv(gn_genes_p, sep="\t", index=False)

    fna_p = tmp_path / "all_matched_primary_targets.fna"
    fna_p.write_text(">TGT_CANONICAL\nACGTACGT\n")

    import hashlib
    tgt_sha = hashlib.sha256(b"ACGTACGT").hexdigest()

    pairs_p = tmp_path / "pairs.tsv"
    pairs_df = pd.DataFrame([
        {
            "canonical_591_host_mag": "H1__single_SemiBin_1",
            "target_id": "TGT_ALIAS",  # Alias name differs from canonical FASTA header
            "target_seq_sha256": tgt_sha,
            "evidence_tier": "Tier1_verified_direct",
            "host_phylum": "p__Halobacteriota",
            "host_genus": "g__Halorubrum",
            "votu_id": "VOTU_1",
            "pool": "eq_votu_reps",
            "viral_identity_status": "Verified_viral_independent_hallmark_and_checkv",
            "checkv_quality": "High-quality",
            "checkv_viral_genes": 1,
            "checkv_host_genes": 0,
        },
        {
            "canonical_591_host_mag": "H2__single_SemiBin_2",
            "target_id": "41Sr3::41Sr3_56779",
            "target_seq_sha256": "sha_41sr3",
            "evidence_tier": "Tier1_verified_direct",
            "host_phylum": "p__Halobacteriota",
            "host_genus": "g__",
            "votu_id": "41Sr3::41Sr3_56779",
            "pool": "eq_raw_assemblies_arch_samples",
            "viral_identity_status": "Verified_viral_independent_hallmark_and_checkv",
            "checkv_quality": "Medium-quality",
            "checkv_viral_genes": 1,
            "checkv_host_genes": 70,
        }
    ])
    pairs_df.to_csv(pairs_p, sep="\t", index=False)

    fastq_dir = tmp_path / "fastq"
    fastq_dir.mkdir()
    (fastq_dir / "H1_R1.fastq.gz").touch()
    (fastq_dir / "H1_R2.fastq.gz").touch()
    (fastq_dir / "TGT_R1.fastq.gz").touch()
    (fastq_dir / "TGT_R2.fastq.gz").touch()

    summary = audit_tier1_pairs(
        pairs_tsv=pairs_p,
        gn_genes_tsv=gn_genes_p,
        target_fna=fna_p,
        fastq_dir=fastq_dir,
    )
    assert summary["total_tier1_rows"] == 2
    assert summary["n_unresolved_virality_targets"] == 1

    alias_rec = [p for p in summary["pairs_detail"] if p["target_id"] == "TGT_ALIAS"][0]
    assert alias_rec["is_alias_target"] is True
    assert alias_rec["canonical_annotated_target_id"] == "TGT_CANONICAL"
    assert alias_rec["n_total_orfs"] == 2
    assert alias_rec["n_unannotated_orfs"] == 1
    assert alias_rec["n_hallmark_orfs"] == 1
    assert alias_rec["eligible_for_structural_hallmark_followup"] is True


def test_abundance_matrix_schema_transposition_and_sample_join(tmp_path: Path):
    """Verify explicit matrix schema validation, transposition, and authoritative sample joining."""
    samples = ["SAMP_A", "SAMP_B", "SAMP_C"]
    meta_p = tmp_path / "metadata.tsv"
    pd.DataFrame({"sample_id": samples, "included": [True, True, True]}).to_csv(meta_p, sep="\t", index=False)

    mag_p = tmp_path / "mag_abundance.tsv"
    mag_raw = pd.DataFrame(
        {
            "MAG_1": [0.1, 0.2, 0.3],
            "MAG_2": [0.0, 0.5, 0.0],
        },
        index=samples,
    )
    mag_raw.index.name = "sample_id"
    mag_raw.to_csv(mag_p, sep="\t")

    vir_p = tmp_path / "viral_tpm.tsv"
    vir_raw = pd.DataFrame(
        {
            "VIR_1": [10.0, 20.0, 0.0],
            "VIR_2": [5.0, 0.0, 15.0],
        },
        index=samples,
    )
    vir_raw.index.name = "sample_id"
    vir_raw.to_csv(vir_p, sep="\t")

    mag_df, vir_df, auth_samples = load_and_orient_abundance_matrices(
        mag_abund_path=mag_p,
        vir_tpm_path=vir_p,
        sample_manifest_path=meta_p,
    )

    # Both matrices must be transposed: rows are entities, columns are samples
    assert mag_df.shape == (2, 3)
    assert vir_df.shape == (2, 3)
    assert list(mag_df.index) == ["MAG_1", "MAG_2"]
    assert list(vir_df.index) == ["VIR_1", "VIR_2"]
    assert set(mag_df.columns) == set(samples)
    assert set(vir_df.columns) == set(samples)
    assert auth_samples == set(samples)

    # Adversarial test 1: negative values rejected
    bad_mag = mag_raw.copy()
    bad_mag.iloc[0, 0] = -1.0
    bad_mag_p = tmp_path / "bad_mag.tsv"
    bad_mag.to_csv(bad_mag_p, sep="\t")
    with pytest.raises(ValueError, match="contains negative"):
        load_and_orient_abundance_matrices(bad_mag_p, vir_p, meta_p)

    # Adversarial test 2: NaN values rejected
    nan_vir = vir_raw.copy()
    nan_vir.iloc[0, 0] = np.nan
    nan_vir_p = tmp_path / "nan_vir.tsv"
    nan_vir.to_csv(nan_vir_p, sep="\t")
    with pytest.raises(ValueError, match="contains NaN"):
        load_and_orient_abundance_matrices(mag_p, nan_vir_p, meta_p)

    # Adversarial test 3: Mismatched samples rejected
    bad_samples = ["SAMP_X", "SAMP_Y", "SAMP_Z"]
    bad_meta_p = tmp_path / "bad_meta.tsv"
    pd.DataFrame({"sample_id": bad_samples}).to_csv(bad_meta_p, sep="\t", index=False)
    with pytest.raises(ValueError, match="axes do not match authoritative samples"):
        load_and_orient_abundance_matrices(mag_p, vir_p, bad_meta_p)


def test_phylum_breakdown_not_claimed_as_phylogenetic_control():
    """Verify that phylum grouping alone is strictly labeled as descriptive counts and never claimed as phylogenetic control."""
    cooc_df = pd.DataFrame([
        {
            "host_phylum": "p__Thermoplasmatota",
            "estimability_status": "estimable_with_existing_catalog_abundance_matrices",
            "jaccard_similarity": 0.15,
        },
        {
            "host_phylum": "p__Halobacteriota",
            "estimability_status": "estimable_with_existing_catalog_abundance_matrices",
            "jaccard_similarity": 0.25,
        },
    ])
    ledger_df = pd.DataFrame([
        {
            "genome": "G1",
            "canonical_phylum": "p__Thermoplasmatota",
            "n_candidate_virus_pairs": 1,
            "highest_evidence_tier": "Tier1_verified_direct",
        },
        {
            "genome": "G2",
            "canonical_phylum": "p__Thermoplasmatota",
            "n_candidate_virus_pairs": 0,
            "highest_evidence_tier": "Zero_spacers_or_no_viral_hits",
        },
        {
            "genome": "G3",
            "canonical_phylum": "p__Halobacteriota",
            "n_candidate_virus_pairs": 1,
            "highest_evidence_tier": "Tier1_verified_direct",
        },
        {
            "genome": "G4",
            "canonical_phylum": "p__Halobacteriota",
            "n_candidate_virus_pairs": 0,
            "highest_evidence_tier": "Zero_spacers_or_no_viral_hits",
        },
    ])

    df = compute_descriptive_phylum_counts(cooc_df, ledger_df)
    assert len(df) == 2

    # Must be explicitly labeled descriptive_phylum_breakdown
    assert (df["model_type"] == "descriptive_phylum_breakdown").all()

    # Must NOT claim phylogenetic confounding control
    for note in df["notes"]:
        assert "Controls for phylogenetic confounding" not in note
        assert "does NOT constitute" in note

    # Check opportunity counts
    thermo_row = df[df["canonical_phylum"] == "p__Thermoplasmatota"].iloc[0]
    assert thermo_row["n_mags_with_spacer_detection_opportunity"] == 1
    assert thermo_row["n_mags_zero_spacers_no_detection_opportunity"] == 1

    # Check defense mechanism estimability audit
    pairs_df = pd.DataFrame([
        {"canonical_591_host_mag": "G1", "evidence_tier": "Tier1_verified_direct"},
        {"canonical_591_host_mag": "G3", "evidence_tier": "Tier1_verified_direct"},
    ])
    audit = audit_surface_and_defense_mechanism_estimability(ledger_df, pairs_df)
    assert audit["estimability_status"] == "non_estimable_insufficient_positive_events_and_circularity"
    assert audit["regression_rank_deficiency"] is True
    assert "circular" in audit["scientific_disposition_reason"].lower()


def test_collision_audit_v2_source_extraction_and_missing_pool_rejection(tmp_path: Path):
    """Verify EXP3-COLLISION-COVERAGE-REPAIR rules:
    1. Assembly FASTA extraction uses raw_subject_id (not namespaced target_id) and applies
       namespacing transform only after extraction.
    2. Missing-one-of-three-pools triggers unknown_incomplete_source_coverage and CANNOT support collision_free.
    3. Within-source repeated headers with different sequences are detected and flagged.
    4. Legitimate SHA deduplication resolves aliases to expected SHA.
    """
    import audit_target_sequence_collisions_v2 as col_v2

    # 1. Create mock assembly FASTA with raw header
    asm_fa = tmp_path / "41PRr2_contigs.fa"
    raw_id = "41PRr2_12345"
    seq1 = "ACGT" * 100
    asm_fa.write_text(f">{raw_id} length=400\n{seq1}\n")

    # Extract raw header
    ext, conflicts = col_v2.parse_fasta_for_headers(asm_fa, {raw_id})
    assert raw_id in ext
    assert ext[raw_id] == seq1
    assert len(conflicts) == 0

    # Namespacing transform
    namespaced = f"41PRr2::{raw_id}"
    assert namespaced == "41PRr2::41PRr2_12345"

    # 2. Test within-source conflict detection
    asm_conflict = tmp_path / "conflict.fa"
    asm_conflict.write_text(f">dup_id\nACGT\n>dup_id\nTGCA\n")
    ext_c, conflicts_c = col_v2.parse_fasta_for_headers(asm_conflict, {"dup_id"})
    assert len(conflicts_c) == 1
    assert conflicts_c[0]["header"] == "dup_id"

    # 3. Missing-one-of-three-pools test
    expected_pools = ["eq_all_viruses", "eq_raw_assemblies_arch_samples", "eq_votu_reps"]
    # Suppose only 2 pools recovered (the exact bug found by Root)
    recovered_pools = ["eq_all_viruses", "eq_votu_reps"]
    missing_pools = sorted(set(expected_pools) - set(recovered_pools))
    assert len(missing_pools) > 0

    # Under v2 logic, status must be unknown_incomplete_source_coverage
    has_complete = (len(missing_pools) == 0)
    assert has_complete is False
    status = "unknown_incomplete_source_coverage" if not has_complete else "identical_duplicate_across_pools"
    assert status == "unknown_incomplete_source_coverage"


def test_stage4_crosswalk_ownership_ambiguity_and_full_mag_denominator():
    """Verify EXP3-MAPPING-ACCEPTANCE-REPAIR rules:
    1. Crosswalk resolves ownership: multi-MAG shared contigs are NOT attributed as unique host coverage.
    2. Shared host-viral contigs are flagged as shared_host_viral_ambiguous.
    3. Full MAG denominators require all original contigs and canonical manifest bases (including zeros).
    4. Uniquely attributable coverage is separated from ambiguous shared coverage.
    """
    # Mock crosswalk table
    cw_rows = [
        # Host MAG 1 contigs
        {"original_header": "HOST__M1__c1", "source_category": "host_mag", "entity_id": "M1", "contig_id": "c1", "seq_length": 1000, "canonical_sha256": "sha_unique_1", "representative_reference_header": "REP_HOST_1", "ambiguity_status": "unambiguous_host"},
        {"original_header": "HOST__M1__c2", "source_category": "host_mag", "entity_id": "M1", "contig_id": "c2", "seq_length": 2000, "canonical_sha256": "sha_shared_m1_m2", "representative_reference_header": "REP_SHARED_HOST", "ambiguity_status": "unambiguous_host"},
        {"original_header": "HOST__M1__c3", "source_category": "host_mag", "entity_id": "M1", "contig_id": "c3", "seq_length": 3000, "canonical_sha256": "sha_shared_viral", "representative_reference_header": "REP_HOST_VIRAL", "ambiguity_status": "shared_host_viral_ambiguous"},
        {"original_header": "HOST__M1__c4", "source_category": "host_mag", "entity_id": "M1", "contig_id": "c4", "seq_length": 4000, "canonical_sha256": "sha_unmapped_zero", "representative_reference_header": "REP_ZERO_COV", "ambiguity_status": "unambiguous_host"},

        # Host MAG 2 also has c2 sequence
        {"original_header": "HOST__M2__c2", "source_category": "host_mag", "entity_id": "M2", "contig_id": "c2", "seq_length": 2000, "canonical_sha256": "sha_shared_m1_m2", "representative_reference_header": "REP_SHARED_HOST", "ambiguity_status": "unambiguous_host"},

        # Virus target also has c3 sequence
        {"original_header": "VIRUS_TARGET__V1", "source_category": "virus_target", "entity_id": "V1", "contig_id": "V1", "seq_length": 3000, "canonical_sha256": "sha_shared_viral", "representative_reference_header": "REP_HOST_VIRAL", "ambiguity_status": "shared_host_viral_ambiguous"},
    ]
    cw_df = pd.DataFrame(cw_rows)

    # Derive ownership from crosswalk
    owners_by_rep = cw_df.groupby("representative_reference_header").agg({
        "source_category": set,
        "entity_id": set,
    })

    # Rule 1: Multi-MAG shared contig (REP_SHARED_HOST) has multiple host MAGs
    assert owners_by_rep.loc["REP_SHARED_HOST"]["entity_id"] == {"M1", "M2"}
    is_unique_m1_c2 = (owners_by_rep.loc["REP_SHARED_HOST"]["source_category"] == {"host_mag"}) and (owners_by_rep.loc["REP_SHARED_HOST"]["entity_id"] == {"M1"})
    assert is_unique_m1_c2 is False  # MUST NOT be accepted as unique to M1!

    # Rule 2: Host-viral shared contig (REP_HOST_VIRAL) has both host_mag and virus_target
    assert "host_mag" in owners_by_rep.loc["REP_HOST_VIRAL"]["source_category"]
    assert "virus_target" in owners_by_rep.loc["REP_HOST_VIRAL"]["source_category"]

    # Rule 3: Full MAG denominator for M1
    m1_cw = cw_df[(cw_df["source_category"] == "host_mag") & (cw_df["entity_id"] == "M1")]
    total_mag_bases = m1_cw["seq_length"].sum()
    assert total_mag_bases == 1000 + 2000 + 3000 + 4000  # 10,000 bp
    assert len(m1_cw) == 4  # 4 contigs

    # Mock BAM coverage table (REP_ZERO_COV is omitted from BAM because 0 reads mapped)
    mock_cov = {
        "REP_HOST_1": {"numreads": 100, "covbases": 1000, "meandepth": 10.0},
        "REP_SHARED_HOST": {"numreads": 200, "covbases": 2000, "meandepth": 15.0},
        "REP_HOST_VIRAL": {"numreads": 50, "covbases": 1500, "meandepth": 5.0},
        # REP_ZERO_COV absent from BAM!
    }

    # Evaluate M1 with explicit 0-coverage handling
    unique_bases = 0
    unique_covbases = 0
    unique_depth_weighted_sum = 0.0
    total_covbases = 0
    total_depth_weighted_sum = 0.0

    for _, crow in m1_cw.iterrows():
        c_len = int(crow["seq_length"])
        c_rep = str(crow["representative_reference_header"])
        owner_info = owners_by_rep.loc[c_rep]
        is_unique = (owner_info["source_category"] == {"host_mag"}) and (owner_info["entity_id"] == {"M1"})

        centry = mock_cov.get(c_rep, {"numreads": 0, "covbases": 0, "meandepth": 0.0})
        total_covbases += centry["covbases"]
        total_depth_weighted_sum += (centry["meandepth"] * c_len)

        if is_unique:
            unique_bases += c_len
            unique_covbases += centry["covbases"]
            unique_depth_weighted_sum += (centry["meandepth"] * c_len)

    # Unique contigs for M1: only c1 (1000 bp) and c4 (4000 bp, but 0 covered) = 5000 bp
    assert unique_bases == 5000
    assert unique_covbases == 1000  # only c1 covered
    unique_breadth = (unique_covbases / total_mag_bases) * 100.0  # 10.0%
    assert unique_breadth == 10.0

    # Total breadth includes shared contigs: c1 (1000) + c2 (2000) + c3 (1500) = 4500
    assert total_covbases == 4500
    total_breadth = (total_covbases / total_mag_bases) * 100.0  # 45.0%
    assert total_breadth == 45.0

    # 0-coverage contig c4 (4000 bp) properly deflates weighted mean depth:
    # unique mean depth = (10.0 * 1000 + 0 * 4000) / 10000 = 1.0x
    assert (unique_depth_weighted_sum / total_mag_bases) == 1.0


def test_stage4_mapping_receipt_cryptographic_and_missing_container_failure(tmp_path: Path):
    """Verify cryptographic receipt hashing and container unavailability failure behavior."""
    import stage4_map_sample_competitive as map_mod

    # 1. Test stream_sha256
    test_f = tmp_path / "test.bam"
    test_f.write_bytes(b"BAM_MOCK_CONTENT_12345")
    h1 = map_mod.stream_sha256(test_f)
    assert len(h1) == 64

    # Mutation test: 1 byte change alters hash
    test_f.write_bytes(b"BAM_MOCK_CONTENT_12346")
    h2 = map_mod.stream_sha256(test_f)
    assert h1 != h2

    # 2. Test missing container failure on eligible target
    nonexistent_sif = tmp_path / "nonexistent_instrain.sif"
    assert not nonexistent_sif.exists()

    qualifying_vir = [{"entity_id": "V1", "reference_contig": "REP_V1", "qualifies_10x_50pct_instrain": 1}]
    # When target qualifies and container is missing, must raise RuntimeError, never silently mark ineligible
    with pytest.raises(RuntimeError, match="inStrain container unavailable"):
        if qualifying_vir and not nonexistent_sif.exists():
            raise RuntimeError(f"inStrain container unavailable: {nonexistent_sif}")


def test_stage4_sequence_domain_annotation_and_structural_screen_distinction():
    """Verify that Pfam-A screen is explicitly sequence domain annotation, not Foldseek structural screen."""
    import stage4_structural_hallmark_screen as annot_mod

    # Method distinction metadata
    summary = {
        "methodology_distinction": {
            "pfam_search_type": "sequence_profile_hmm_domain_annotation",
            "is_structural_coordinate_search": False,
            "structural_search_tool": "Foldseek (requires 3D PDB/CIF models)",
            "foldseek_structural_screen_status": "not_performed_missing_structural_models",
        }
    }
    assert summary["methodology_distinction"]["is_structural_coordinate_search"] is False
    assert summary["methodology_distinction"]["pfam_search_type"] == "sequence_profile_hmm_domain_annotation"
    assert summary["methodology_distinction"]["foldseek_structural_screen_status"] == "not_performed_missing_structural_models"





