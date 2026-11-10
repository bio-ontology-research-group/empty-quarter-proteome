#!/usr/bin/env python3
"""7-Gate Standalone Reproducibility & Replay Verifier for Exp3 Raw-vs-MAPQ20 Coverage QC (`exp3_viruses/coverage_qc`)."""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))

import stage4_raw_vs_mapq20_focal_coverage as raw_qc  # noqa: E402
from stage4_competitive_common import load_analysis_pairs, stream_sha256  # noqa: E402
from stage4_map_sample_competitive import COVERAGE_COLUMNS, compute_coverage_metrics  # noqa: E402


def assert_true(cond: bool, msg: str) -> None:
    if not cond:
        raise AssertionError(msg)


def verify_all(root: Path = ROOT, verify_9sr2_benchmark: bool = True) -> dict[str, Any]:
    in_dir = root / "inputs"
    res_dir = root / "results"
    man_dir = root / "manifests"

    # [1/7] SHA256 Package Manifest & Executed Code Hash Integrity
    sha_manifest_path = man_dir / "sha256_manifest.tsv"
    assert_true(sha_manifest_path.is_file(), f"Missing {sha_manifest_path}")
    sha_manifest = pd.read_csv(sha_manifest_path, sep="\t")
    assert_true(
        list(sha_manifest.columns) == ["relative_path", "bytes", "sha256"],
        f"Unexpected sha256_manifest columns: {list(sha_manifest.columns)}",
    )
    assert_true(
        len(sha_manifest) >= 60,
        f"Expected >=60 packaged files in sha256_manifest.tsv, got {len(sha_manifest)}",
    )
    for row in sha_manifest.itertuples(index=False):
        target_path = root / str(row.relative_path)
        assert_true(target_path.is_file(), f"Missing packaged file: {row.relative_path}")
        actual_bytes = int(target_path.stat().st_size)
        assert_true(
            actual_bytes == int(row.bytes),
            f"Byte size mismatch for {row.relative_path}: {actual_bytes} != {row.bytes}",
        )
        actual_sha = stream_sha256(target_path)
        assert_true(
            actual_sha == str(row.sha256),
            f"SHA256 mismatch for {row.relative_path}: {actual_sha} != {row.sha256}",
        )

    code_hashes_path = man_dir / "stage4_raw_vs_mapq20_code_hashes.tsv"
    assert_true(code_hashes_path.is_file(), f"Missing {code_hashes_path}")
    expected_code_hashes: dict[str, str] = {}
    for raw_line in code_hashes_path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        sha_val, rel_str = line.split(None, 1)
        expected_code_hashes[rel_str.strip()] = sha_val.strip()

    for rel_key in [
        "code/stage4_raw_vs_mapq20_focal_coverage.py",
        "code/sbatch_stage4_raw_vs_mapq20.sh",
        "code/stage4_competitive_common.py",
        "code/stage4_map_sample_competitive.py",
        "tests/test_stage4_raw_vs_mapq20.py",
    ]:
        assert_true(rel_key in expected_code_hashes, f"Missing {rel_key} in stage4_raw_vs_mapq20_code_hashes.tsv")
        assert_true(
            stream_sha256(root / rel_key) == expected_code_hashes[rel_key],
            f"Code hash mismatch for {rel_key}",
        )
    print(f"[1/7] PASS: SHA256 package manifest ({len(sha_manifest)} files) and code hashes verified.")

    # [2/7] Exact 14-Sample x 7-Target Scope & Stage 2 vOTU Context
    sub_receipt = in_dir / "stage4_competitive_pipeline_submission.json"
    pair_ledger = in_dir / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    boundary_cfg = in_dir / "stage4_boundary_sensitive_intervals.json"
    stage2_table = in_dir / "stage2_target_viral_identity_table.tsv"

    samples, sub_payload = raw_qc.load_expected_samples(sub_receipt)
    assert_true(len(samples) == 14, f"Expected 14 samples, got {len(samples)}")
    assert_true(sub_payload.get("job_validate_ref") == "15905724", "Expected job_validate_ref == 15905724")
    assert_true(sub_payload.get("job_map_array") == "15905725", "Expected job_map_array == 15905725")

    pairs, reviewed, boundary = load_analysis_pairs(pair_ledger, boundary_cfg)
    assert_true(len(reviewed) == 6 and len(boundary) == 1 and len(pairs) == 7, "Expected 6+1=7 pairs")
    pair_meta_by_sha = raw_qc.build_pair_metadata_by_sha(pairs)
    assert_true(len(pair_meta_by_sha) == 7, f"Expected 7 unique target_seq_sha256, got {len(pair_meta_by_sha)}")

    role_counts_targets = pd.Series([m["analysis_role"] for m in pair_meta_by_sha.values()]).value_counts().to_dict()
    assert_true(role_counts_targets.get("primary") == 5, f"Expected 5 primary targets, got {role_counts_targets}")
    assert_true(
        role_counts_targets.get("sensitivity_unresolved_virality") == 1,
        f"Expected 1 sensitivity_unresolved_virality target, got {role_counts_targets}",
    )
    assert_true(
        role_counts_targets.get("boundary_sensitive_provisional") == 1,
        f"Expected 1 boundary_sensitive_provisional target, got {role_counts_targets}",
    )

    votu_ctx = raw_qc.load_stage2_votu_context(stage2_table, pairs)
    assert_true(len(votu_ctx) == 7, f"Expected 7 vOTU context entries, got {len(votu_ctx)}")
    sha_9sr2 = "96ff8fd33cdd8a3b3dd1db91ba5645a14797b63723f26e1f5ff2d1be838ef32d"
    assert_true(
        votu_ctx[sha_9sr2]["stage2_is_votu_cluster_rep"] == 0
        and votu_ctx[sha_9sr2]["stage2_cluster_votu_rep_id"] == "55Sr3::55Sr3_129891",
        "Unexpected Stage 2 vOTU metadata for 9Sr2::9Sr2_284237",
    )
    print("[2/7] PASS: Exact 14-sample x 7-target scope (98 rows) and Stage 2 vOTU context metadata verified.")

    # [3/7] Reference Validation Sentinel & 14 Sample Mapping Receipt Bindings
    pair_ledger_sha = stream_sha256(pair_ledger)
    boundary_cfg_sha = stream_sha256(boundary_cfg)
    ref_dir = in_dir / "stage4_competitive_reference"
    expected_ref_validation = raw_qc.verify_reference_validation_contract(
        ref_dir=ref_dir,
        expected_pair_ledger_sha256=pair_ledger_sha,
        expected_boundary_config_sha256=boundary_cfg_sha,
    )

    map_dir = in_dir / "stage4_competitive_mapping"
    existing_sample_dirs = sorted(p.parent.name for p in map_dir.glob("*/*_mapping_receipt.json"))
    assert_true(
        existing_sample_dirs == sorted(samples),
        f"Sample receipt directories {existing_sample_dirs} != expected {sorted(samples)}",
    )

    eval_by_sample_sha: dict[tuple[str, str], dict[str, Any]] = {}
    sample_receipt_meta: dict[str, dict[str, Any]] = {}

    for sample in samples:
        s_dir = map_dir / sample
        r_path = s_dir / f"{sample}_mapping_receipt.json"
        d_path = s_dir / "MAPPING_DONE"
        e_path = s_dir / f"{sample}_focal_targets_evaluation.tsv"
        assert_true(r_path.is_file() and d_path.is_file() and e_path.is_file(), f"Missing mapping inputs for {sample}")

        r_sha = stream_sha256(r_path)
        d_sha = stream_sha256(d_path)
        d_payload = json.loads(d_path.read_text())
        assert_true(d_payload.get("receipt_sha256") == r_sha, f"MAPPING_DONE receipt_sha256 mismatch for {sample}")

        r_payload = json.loads(r_path.read_text())
        assert_true(r_payload.get("status") == "completed_verified", f"Receipt status not completed_verified for {sample}")
        assert_true(int(r_payload.get("schema_version", 0)) == 3, f"Receipt schema_version != 3 for {sample}")
        assert_true(r_payload.get("stage") == "stage4_map_sample_competitive", f"Unexpected stage for {sample}")
        assert_true(r_payload.get("sample") == sample, f"Receipt sample mismatch for {sample}: {r_payload.get('sample')}")

        r_inputs = r_payload.get("inputs", {})
        assert_true(r_inputs.get("pair_ledger_sha256") == pair_ledger_sha, f"pair_ledger_sha256 mismatch for {sample}")
        assert_true(r_inputs.get("boundary_config_sha256") == boundary_cfg_sha, f"boundary_config_sha256 mismatch for {sample}")

        ref_binding = r_payload.get("reference_binding", {})
        for bind_key in ["validation_receipt_sha256", "reference_fna_sha256", "crosswalk_sha256", "build_receipt_sha256"]:
            assert_true(
                ref_binding.get(bind_key) == expected_ref_validation[bind_key],
                f"reference_binding.{bind_key} mismatch for {sample}",
            )
        assert_true(
            ref_binding.get("bowtie2_index_sha256") == expected_ref_validation["bowtie2_index_sha256"],
            f"reference_binding.bowtie2_index_sha256 mismatch for {sample}",
        )

        validated_metrics = raw_qc.validate_sample_mapping_metrics(r_payload.get("mapping_metrics"), sample)

        outputs = r_payload.get("outputs", {})
        for sha_key in ["raw_bam_sha256", "filtered_bam_sha256", "coverage_tsv_sha256", "focal_evaluation_tsv_sha256"]:
            raw_qc.require_valid_sha256_hex(outputs.get(sha_key), f"{sample}.outputs.{sha_key}")
        assert_true(
            stream_sha256(e_path) == outputs["focal_evaluation_tsv_sha256"],
            f"focal_evaluation_tsv_sha256 mismatch for {sample}",
        )

        # Verify both receipt focal_evaluations and packaged focal_targets_evaluation.tsv have exact 7 virus_target rows
        receipt_virus_evals = [row for row in r_payload.get("focal_evaluations", []) if row.get("entity_type") == "virus_target"]
        assert_true(len(receipt_virus_evals) == 7, f"Expected 7 virus_target rows in {sample} receipt, got {len(receipt_virus_evals)}")
        receipt_by_sha = {str(row["entity_sha256"]): row for row in receipt_virus_evals}
        assert_true(set(receipt_by_sha) == set(pair_meta_by_sha), f"Receipt virus_target SHA mismatch in {sample}")

        ev_df = pd.read_csv(e_path, sep="\t")
        ev_virus_df = ev_df[ev_df["entity_type"] == "virus_target"].copy()
        assert_true(len(ev_virus_df) == 7, f"Expected 7 virus_target rows in {e_path}, got {len(ev_virus_df)}")
        assert_true(
            set(ev_virus_df["entity_sha256"].astype(str)) == set(pair_meta_by_sha),
            f"Evaluation TSV virus_target SHA set mismatch in {sample}",
        )
        for _, ev_row in ev_virus_df.iterrows():
            t_sha = str(ev_row["entity_sha256"])
            rec_row = receipt_by_sha[t_sha]
            assert_true(
                int(ev_row["filtered_reads_mapped"]) == int(rec_row["filtered_reads_mapped"]),
                f"filtered_reads_mapped mismatch between TSV and receipt for {sample} {t_sha}",
            )
            assert_true(
                math.isclose(float(ev_row["horizontal_breadth_pct"]), float(rec_row["horizontal_breadth_pct"]), rel_tol=1e-9, abs_tol=1e-9),
                f"horizontal_breadth_pct mismatch between TSV and receipt for {sample} {t_sha}",
            )
            assert_true(
                math.isclose(float(ev_row["mean_depth"]), float(rec_row["mean_depth"]), rel_tol=1e-9, abs_tol=1e-9),
                f"mean_depth mismatch between TSV and receipt for {sample} {t_sha}",
            )
            eval_by_sample_sha[(sample, t_sha)] = rec_row

        sample_receipt_meta[sample] = {
            "mapping_receipt_sha256": r_sha,
            "done_sha256": d_sha,
            "raw_bam_sha256": outputs["raw_bam_sha256"],
            "filtered_bam_sha256": outputs["filtered_bam_sha256"],
            "coverage_tsv_sha256": outputs["coverage_tsv_sha256"],
            "focal_evaluation_tsv_sha256": outputs["focal_evaluation_tsv_sha256"],
            "validated_metrics": validated_metrics,
        }

    print("[3/7] PASS: Reference validation sentinel and all 14 sample mapping receipts & bindings verified.")

    # [4/7] Independent Replay of 98-Row Table from Packaged samtools Tabular Outputs & Producer Evaluations
    focal_raw_df = pd.read_csv(res_dir / "stage4_raw_vs_mapq20_focal_coverage_98rows.tsv", sep="\t", keep_default_na=False)
    focal_df = pd.read_csv(res_dir / "stage4_raw_vs_mapq20_focal_coverage_98rows.tsv", sep="\t")
    raw_tab_df = pd.read_csv(res_dir / "stage4_raw_vs_mapq20_samtools_raw_tabular_98rows.tsv", sep="\t")
    filt_tab_df = pd.read_csv(res_dir / "stage4_raw_vs_mapq20_samtools_filtered_tabular_98rows.tsv", sep="\t")
    target_sum_df = pd.read_csv(res_dir / "stage4_raw_vs_mapq20_target_summary.tsv", sep="\t")
    sample_sum_df = pd.read_csv(res_dir / "stage4_raw_vs_mapq20_sample_summary.tsv", sep="\t")

    expected_tab_cols = ["sample", "#rname", "startpos", "endpos", "numreads", "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"]
    assert_true(list(raw_tab_df.columns) == expected_tab_cols, f"Unexpected raw_tab_df columns: {list(raw_tab_df.columns)}")
    assert_true(list(filt_tab_df.columns) == expected_tab_cols, f"Unexpected filt_tab_df columns: {list(filt_tab_df.columns)}")

    assert_true(len(focal_df) == 98, f"Expected 98 focal rows, got {len(focal_df)}")
    assert_true(len(raw_tab_df) == 98, f"Expected 98 raw tabular rows, got {len(raw_tab_df)}")
    assert_true(len(filt_tab_df) == 98, f"Expected 98 filtered tabular rows, got {len(filt_tab_df)}")
    assert_true(len(target_sum_df) == 7, f"Expected 7 target summary rows, got {len(target_sum_df)}")
    assert_true(len(sample_sum_df) == 14, f"Expected 14 sample summary rows, got {len(sample_sum_df)}")

    # Enforce exact Cartesian product of 14 samples x 7 target SHAs and unique keys
    expected_cartesian = {(s, t) for s in samples for t in pair_meta_by_sha}
    observed_pairs = list(zip(focal_df["sample"].astype(str), focal_df["raw_oriented_target_sha256"].astype(str)))
    assert_true(len(observed_pairs) == 98 and len(set(observed_pairs)) == 98, "Duplicate (sample, raw_oriented_target_sha256) keys")
    assert_true(set(observed_pairs) == expected_cartesian, "Observed 98 rows do not match exact 14 x 7 Cartesian product")

    observed_contig_keys = list(zip(focal_df["sample"].astype(str), focal_df["reference_contig"].astype(str)))
    assert_true(len(set(observed_contig_keys)) == 98, "Duplicate (sample, reference_contig) keys in focal_df")
    assert_true(
        list(zip(raw_tab_df["sample"].astype(str), raw_tab_df["#rname"].astype(str))) == observed_contig_keys,
        "raw_tab_df (sample, #rname) order/keys do not match focal_df",
    )
    assert_true(
        list(zip(filt_tab_df["sample"].astype(str), filt_tab_df["#rname"].astype(str))) == observed_contig_keys,
        "filt_tab_df (sample, #rname) order/keys do not match focal_df",
    )

    # Role counts across 98 rows
    role_counts_98 = focal_df["analysis_role"].value_counts().to_dict()
    assert_true(role_counts_98 == {"primary": 70, "sensitivity_unresolved_virality": 14, "boundary_sensitive_provisional": 14}, f"Unexpected 98-row role counts: {role_counts_98}")

    raw_by_key = {(str(r["sample"]), str(r["#rname"])): r for _, r in raw_tab_df.iterrows()}
    filt_by_key = {(str(r["sample"]), str(r["#rname"])): r for _, r in filt_tab_df.iterrows()}

    for idx in range(98):
        f_row = focal_df.iloc[idx]
        f_raw_str = focal_raw_df.iloc[idx]
        s_id = str(f_row["sample"])
        t_sha = str(f_row["raw_oriented_target_sha256"])
        ref_contig = str(f_row["reference_contig"])
        r_row = raw_by_key[(s_id, ref_contig)]
        m_row = filt_by_key[(s_id, ref_contig)]
        feval = eval_by_sample_sha[(s_id, t_sha)]
        pmeta = pair_meta_by_sha[t_sha]
        s_meta = sample_receipt_meta[s_id]

        length_bp = int(f_row["length_bp"])
        assert_true(length_bp == int(feval["length_bp"]) == int(r_row["endpos"]) == int(m_row["endpos"]), f"Row {idx} length_bp mismatch")
        assert_true(int(f_row["startpos"]) == 1 == int(r_row["startpos"]) == int(m_row["startpos"]), f"Row {idx} startpos != 1")

        # Raw coverage fields vs samtools raw tabular output
        raw_reads = int(r_row["numreads"])
        raw_covbases = int(r_row["covbases"])
        raw_samtools_pct = float(r_row["coverage"])
        recomputed_raw_breadth = 100.0 * raw_covbases / length_bp
        raw_depth = float(r_row["meandepth"])
        raw_baseq = float(r_row["meanbaseq"])
        raw_mapq = float(r_row["meanmapq"])

        assert_true(int(f_row["raw_numreads"]) == raw_reads, f"Row {idx} raw_numreads mismatch")
        assert_true(int(f_row["raw_covbases"]) == raw_covbases, f"Row {idx} raw_covbases mismatch")
        assert_true(math.isclose(float(f_row["raw_coverage_samtools_pct"]), raw_samtools_pct, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} raw_coverage_samtools_pct mismatch")
        assert_true(math.isclose(float(f_row["raw_horizontal_breadth_pct"]), recomputed_raw_breadth, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} raw_horizontal_breadth_pct mismatch")
        assert_true(math.isclose(float(f_row["raw_meandepth"]), raw_depth, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} raw_meandepth mismatch")
        assert_true(math.isclose(float(f_row["raw_meanbaseq"]), raw_baseq, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} raw_meanbaseq mismatch")
        assert_true(math.isclose(float(f_row["raw_meanmapq"]), raw_mapq, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} raw_meanmapq mismatch")
        assert_true(int(f_row["raw_reaches_10x_50pct_qc_only"]) == int(raw_depth >= 10.0 and recomputed_raw_breadth >= 50.0), f"Row {idx} raw_reaches_10x_50pct_qc_only mismatch")

        # Filtered coverage fields vs samtools filtered tabular output and producer focal_evaluation
        filt_reads = int(m_row["numreads"])
        filt_covbases = int(m_row["covbases"])
        filt_samtools_pct = float(m_row["coverage"])
        recomputed_filt_breadth = 100.0 * filt_covbases / length_bp
        filt_depth = float(m_row["meandepth"])
        filt_baseq = float(m_row["meanbaseq"])
        filt_mapq = float(m_row["meanmapq"])

        assert_true(int(f_row["filtered_numreads"]) == filt_reads == int(feval["filtered_reads_mapped"]), f"Row {idx} filtered_numreads mismatch")
        assert_true(int(f_row["unique_filtered_reads_mapped"]) == int(feval["unique_filtered_reads_mapped"]), f"Row {idx} unique_filtered_reads_mapped mismatch")
        assert_true(int(f_row["ambiguous_filtered_reads_mapped"]) == int(feval["ambiguous_filtered_reads_mapped"]), f"Row {idx} ambiguous_filtered_reads_mapped mismatch")
        assert_true(
            int(f_row["unique_filtered_reads_mapped"]) + int(f_row["ambiguous_filtered_reads_mapped"]) == filt_reads,
            f"Row {idx} unique + ambiguous != filtered_numreads",
        )
        assert_true(int(f_row["filtered_covbases"]) == filt_covbases, f"Row {idx} filtered_covbases mismatch")
        assert_true(math.isclose(float(f_row["filtered_coverage_samtools_pct"]), filt_samtools_pct, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} filtered_coverage_samtools_pct mismatch")
        assert_true(
            math.isclose(float(f_row["filtered_horizontal_breadth_pct"]), recomputed_filt_breadth, rel_tol=1e-9, abs_tol=1e-9)
            and math.isclose(float(f_row["filtered_horizontal_breadth_pct"]), float(feval["horizontal_breadth_pct"]), rel_tol=1e-9, abs_tol=1e-9),
            f"Row {idx} filtered_horizontal_breadth_pct mismatch",
        )
        assert_true(
            math.isclose(float(f_row["filtered_meandepth"]), filt_depth, rel_tol=1e-7, abs_tol=1e-7)
            and math.isclose(float(f_row["filtered_meandepth"]), float(feval["mean_depth"]), rel_tol=1e-7, abs_tol=1e-7),
            f"Row {idx} filtered_meandepth mismatch",
        )
        assert_true(math.isclose(float(f_row["filtered_meanbaseq"]), filt_baseq, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} filtered_meanbaseq mismatch")
        assert_true(math.isclose(float(f_row["filtered_meanmapq"]), filt_mapq, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} filtered_meanmapq mismatch")

        # Receipt hash bindings per row
        assert_true(str(f_row["mapping_receipt_sha256"]) == s_meta["mapping_receipt_sha256"], f"Row {idx} mapping_receipt_sha256 mismatch")
        assert_true(str(f_row["raw_bam_sha256"]) == s_meta["raw_bam_sha256"], f"Row {idx} raw_bam_sha256 mismatch")
        assert_true(str(f_row["filtered_bam_sha256"]) == s_meta["filtered_bam_sha256"], f"Row {idx} filtered_bam_sha256 mismatch")

        # Replay classification
        qualifies_primary = int(feval["qualifies_10x_50pct_instrain"])
        source_cat_unique = bool(feval["uniquely_attributable"])
        assert_true(int(f_row["qualifies_primary_10x_50pct_instrain"]) == qualifies_primary, f"Row {idx} qualifies_primary mismatch")
        assert_true(bool(f_row["source_category_uniquely_viral"]) == source_cat_unique, f"Row {idx} source_category_uniquely_viral mismatch")

        re_class, re_reason, re_note = raw_qc.classify_recruitment_qc(
            raw_numreads=raw_reads,
            raw_meandepth=raw_depth,
            raw_breadth_pct=recomputed_raw_breadth,
            filtered_numreads=filt_reads,
            filtered_meandepth=float(feval["mean_depth"]),
            filtered_breadth_pct=float(feval["horizontal_breadth_pct"]),
            qualifies_primary=qualifies_primary,
            analysis_role=str(pmeta["analysis_role"]),
            source_category_uniquely_viral=source_cat_unique,
        )
        assert_true(str(f_row["recruitment_qc_class"]) == re_class, f"Row {idx} recruitment_qc_class replay mismatch")
        assert_true(str(f_row["non_estimability_reason"]) == re_reason, f"Row {idx} non_estimability_reason replay mismatch")
        assert_true(str(f_row["scientific_qualification"]) == re_note, f"Row {idx} scientific_qualification replay mismatch")

        # Verify zero-denominator NA semantics (empty string in TSV, NaN in pandas, never 0)
        if raw_reads == 0:
            assert_true(str(f_raw_str["retained_read_fraction"]) == "", f"Row {idx} expected empty retained_read_fraction when raw_numreads==0")
            assert_true(str(f_raw_str["lost_to_mapq20_read_fraction"]) == "", f"Row {idx} expected empty lost_to_mapq20_read_fraction when raw_numreads==0")
            assert_true(pd.isna(f_row["retained_read_fraction"]) and pd.isna(f_row["lost_to_mapq20_read_fraction"]), f"Row {idx} expected NaN read fractions")
        else:
            assert_true(math.isclose(float(f_row["retained_read_fraction"]), filt_reads / raw_reads, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} retained_read_fraction mismatch")
            assert_true(math.isclose(float(f_row["lost_to_mapq20_read_fraction"]), (raw_reads - filt_reads) / raw_reads, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} lost_to_mapq20_read_fraction mismatch")

        if raw_covbases == 0:
            assert_true(str(f_raw_str["retained_breadth_fraction"]) == "" and pd.isna(f_row["retained_breadth_fraction"]), f"Row {idx} expected empty retained_breadth_fraction when raw_covbases==0")
        else:
            assert_true(math.isclose(float(f_row["retained_breadth_fraction"]), filt_covbases / raw_covbases, rel_tol=1e-9, abs_tol=1e-9), f"Row {idx} retained_breadth_fraction mismatch")

        if raw_depth == 0.0:
            assert_true(str(f_raw_str["retained_meandepth_fraction"]) == "" and pd.isna(f_row["retained_meandepth_fraction"]), f"Row {idx} expected empty retained_meandepth_fraction when raw_meandepth==0")
        else:
            assert_true(math.isclose(float(f_row["retained_meandepth_fraction"]), float(feval["mean_depth"]) / raw_depth, rel_tol=1e-7, abs_tol=1e-7), f"Row {idx} retained_meandepth_fraction mismatch")

    # Replay target_summary_df and sample_summary_df arithmetic
    replayed_target_df = raw_qc.summarize_by_target(focal_df)
    pd.testing.assert_frame_equal(replayed_target_df, target_sum_df, check_dtype=False, rtol=1e-9, atol=1e-9)

    for _, s_row in sample_sum_df.iterrows():
        s_id = str(s_row["sample"])
        sub_f = focal_df[focal_df["sample"] == s_id]
        assert_true(len(sub_f) == 7, f"Expected 7 rows for sample {s_id}")
        assert_true(int(s_row["focal7_raw_reads_sum"]) == int(sub_f["raw_numreads"].sum()), f"focal7_raw_reads_sum mismatch for {s_id}")
        assert_true(int(s_row["focal7_filtered_reads_sum"]) == int(sub_f["filtered_numreads"].sum()), f"focal7_filtered_reads_sum mismatch for {s_id}")
        assert_true(int(s_row["focal7_lost_to_mapq20_reads_sum"]) == int(sub_f["lost_to_mapq20_reads"].sum()), f"focal7_lost_to_mapq20_reads_sum mismatch for {s_id}")
        assert_true(int(s_row["focal7_n_targets_raw_reads_gt0"]) == int((sub_f["raw_numreads"] > 0).sum()), f"focal7_n_targets_raw_reads_gt0 mismatch for {s_id}")
        assert_true(int(s_row["focal7_n_targets_raw_reads_eq0"]) == int((sub_f["raw_numreads"] == 0).sum()), f"focal7_n_targets_raw_reads_eq0 mismatch for {s_id}")
        assert_true(
            int(s_row["focal7_n_targets_raw_reads_gt0"]) + int(s_row["focal7_n_targets_raw_reads_eq0"]) == 7,
            f"Raw gt0 + eq0 != 7 for {s_id}",
        )

    print("[4/7] PASS: Independent 98-row, 7-target summary, and 14-sample summary replay verified.")

    # [5/7] Arithmetic Invariants & Non-Causal Scientific Classification Audit
    for idx, f_row in focal_df.iterrows():
        raw_n = int(f_row["raw_numreads"])
        filt_n = int(f_row["filtered_numreads"])
        raw_cb = int(f_row["raw_covbases"])
        filt_cb = int(f_row["filtered_covbases"])
        raw_d = float(f_row["raw_meandepth"])
        filt_d = float(f_row["filtered_meandepth"])
        raw_b = float(f_row["raw_horizontal_breadth_pct"])
        filt_b = float(f_row["filtered_horizontal_breadth_pct"])

        assert_true(raw_n >= filt_n >= 0, f"Row {idx}: raw_numreads < filtered_numreads")
        assert_true(raw_cb >= filt_cb >= 0, f"Row {idx}: raw_covbases < filtered_covbases")
        assert_true(raw_d + 1e-6 >= filt_d >= 0.0, f"Row {idx}: raw_meandepth < filtered_meandepth")
        assert_true(0.0 <= raw_b <= 100.0 and 0.0 <= filt_b <= 100.0, f"Row {idx}: breadth out of [0,100]")
        assert_true(int(f_row["lost_to_mapq20_reads"]) == raw_n - filt_n, f"Row {idx}: lost_to_mapq20_reads mismatch")
        assert_true(int(f_row["lost_to_mapq20_covbases"]) == raw_cb - filt_cb, f"Row {idx}: lost_to_mapq20_covbases mismatch")

        if str(f_row["recruitment_qc_class"]) == "high_raw_support_10x_50pct_lost_at_mapq20":
            assert_true(
                str(f_row["non_estimability_reason"])
                == "high_raw_support_lost_at_mapq20_exact_reference_attribution_unresolved",
                f"Row {idx}: unexpected non_estimability_reason for high_raw_support_10x_50pct_lost_at_mapq20",
            )
            assert_true(
                "exact-reference attribution unresolved; compatible with competitive ambiguity or weak alignment discrimination"
                in str(f_row["scientific_qualification"]),
                f"Row {idx}: unexpected scientific_qualification for high_raw_support_10x_50pct_lost_at_mapq20",
            )

    print("[5/7] PASS: All 98 rows satisfy monotonic raw>=MAPQ20 arithmetic invariants and non-causal QC wording.")

    # [6/7] Authoritative 9Sr2 Source-Target Benchmark Verification
    auth_log = in_dir / "exp3_9Sr2_source_target_unfiltered_coverage.log"
    assert_true(auth_log.is_file(), f"Missing {auth_log}")
    auth_frame = compute_coverage_metrics(auth_log.read_text(), {"VIRUS_TARGET__9Sr2::9Sr2_284237": 56873})
    assert_true(len(auth_frame) == 1, "Expected 1 row in authoritative 9Sr2 log")
    auth_check = raw_qc.verify_authoritative_9sr2_row(focal_df, verify_benchmark=verify_9sr2_benchmark)
    if verify_9sr2_benchmark:
        assert_true(bool(auth_check.get("verified")), "Authoritative 9Sr2 verification failed")
    print("[6/7] PASS: Authoritative 9Sr2 source-target benchmark verified.")

    # [7/7] Sentinel, Summary JSON & Slurm Execution Receipt Verification
    summary_path = res_dir / "stage4_raw_vs_mapq20_summary.json"
    done_path = res_dir / "STAGE4_RAW_VS_MAPQ20_DONE"
    slurm_path = man_dir / "slurm_execution_receipts.tsv"
    assert_true(summary_path.is_file() and done_path.is_file() and slurm_path.is_file(), "Missing summary/done/slurm receipt")

    summary_payload = json.loads(summary_path.read_text())
    done_payload = json.loads(done_path.read_text())
    assert_true(done_payload.get("summary_sha256") == stream_sha256(summary_path), "DONE summary_sha256 mismatch")
    assert_true(summary_payload.get("status") == "completed_verified", "Summary status != completed_verified")
    assert_true(int(summary_payload.get("schema_version", 0)) == 1, "Summary schema_version != 1")
    assert_true(summary_payload.get("stage") == "stage4_raw_vs_mapq20_focal_coverage", "Summary stage mismatch")
    if verify_9sr2_benchmark:
        assert_true(bool(summary_payload.get("authoritative_9Sr2_verification", {}).get("verified")), "Summary 9Sr2 verification false")

    scope_counts = summary_payload.get("scope_counts", {})
    assert_true(int(scope_counts.get("n_samples_expected_and_verified", 0)) == 14, "n_samples_expected_and_verified != 14")
    assert_true(int(scope_counts.get("n_unique_viral_targets", 0)) == 7, "n_unique_viral_targets != 7")
    assert_true(int(scope_counts.get("n_sample_target_rows", 0)) == 98, "n_sample_target_rows != 98")
    assert_true(
        int(scope_counts.get("n_rows_measured_zero_raw_reads", -1)) + int(scope_counts.get("n_rows_raw_reads_gt0", -1)) == 98,
        "scope_counts raw zero + raw gt0 != 98",
    )
    assert_true(
        sum(int(v) for v in summary_payload.get("recruitment_qc_class_counts_all98", {}).values()) == 98,
        "recruitment_qc_class_counts_all98 sum != 98",
    )

    out_meta = summary_payload.get("outputs", {})
    for base_key, filename in [
        ("focal_coverage_98rows_tsv", "stage4_raw_vs_mapq20_focal_coverage_98rows.tsv"),
        ("samtools_raw_tabular_98rows_tsv", "stage4_raw_vs_mapq20_samtools_raw_tabular_98rows.tsv"),
        ("samtools_filtered_tabular_98rows_tsv", "stage4_raw_vs_mapq20_samtools_filtered_tabular_98rows.tsv"),
        ("target_summary_tsv", "stage4_raw_vs_mapq20_target_summary.tsv"),
        ("sample_summary_tsv", "stage4_raw_vs_mapq20_sample_summary.tsv"),
    ]:
        p = res_dir / filename
        assert_true(p.is_file(), f"Missing output file {p}")
        assert_true(stream_sha256(p) == str(out_meta.get(f"{base_key}_sha256")), f"Output SHA256 mismatch for {base_key}")
        assert_true(int(p.stat().st_size) == int(out_meta.get(f"{base_key}_bytes", -1)), f"Output bytes mismatch for {base_key}")

    slurm_df = pd.read_csv(slurm_path, sep="\t")
    assert_true(len(slurm_df) >= 1, "Expected >=1 row in slurm_execution_receipts.tsv")
    raw_qc_rows = slurm_df[slurm_df["pipeline_stage"] == "stage4_raw_vs_mapq20_focal_coverage"]
    assert_true(len(raw_qc_rows) == 1, "Expected 1 stage4_raw_vs_mapq20_focal_coverage row in slurm_execution_receipts.tsv")
    assert_true(str(raw_qc_rows.iloc[0]["state"]) == "COMPLETED", f"Slurm job state != COMPLETED: {raw_qc_rows.iloc[0]['state']}")
    assert_true(str(raw_qc_rows.iloc[0]["exit_code"]) == "0:0", f"Slurm job exit_code != 0:0: {raw_qc_rows.iloc[0]['exit_code']}")

    print("[7/7] PASS: STAGE4_RAW_VS_MAPQ20_DONE sentinel, summary JSON, and Slurm execution receipts verified.")
    print("ALL 7 GATES PASSED.")
    return {
        "status": "ALL_7_GATES_PASSED",
        "manifest_files": len(sha_manifest),
        "focal_rows": len(focal_df),
        "qc_class_counts": summary_payload.get("recruitment_qc_class_counts_all98", {}),
    }


if __name__ == "__main__":
    verify_all()
