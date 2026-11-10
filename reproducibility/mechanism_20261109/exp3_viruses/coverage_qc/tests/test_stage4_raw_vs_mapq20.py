from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

EXP3_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = EXP3_ROOT / "code"
sys.path.insert(0, str(CODE_DIR))

import stage4_raw_vs_mapq20_focal_coverage as raw_qc  # noqa: E402
from stage4_competitive_common import load_analysis_pairs, stream_sha256  # noqa: E402
from stage4_map_sample_competitive import compute_coverage_metrics  # noqa: E402


def run_cmd(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def resolve_test_input(*candidates: Path) -> Path:
    for cand in candidates:
        if cand.is_file():
            return cand
    raise FileNotFoundError(f"None of candidate input paths exist: {[str(c) for c in candidates]}")



def make_valid_reference_validation_fixture(
    ref_dir: Path,
    pair_ledger_sha: str,
    boundary_config_sha: str,
) -> dict[str, object]:
    ref_dir.mkdir(parents=True, exist_ok=True)
    index_files = [
        {"name": f"unified_competitive_reference{suffix}.bt2l", "sha256": f"{idx + 1}" * 64}
        for idx, suffix in enumerate([".1", ".2", ".3", ".4", ".rev.1", ".rev.2"])
    ]
    val_receipt = {
        "status": "completed_verified",
        "schema_version": 3,
        "stage": "stage4_validate_competitive_reference",
        "pair_ledger": {"sha256": pair_ledger_sha},
        "boundary_config": {"sha256": boundary_config_sha},
        "build_receipt": {"sha256": "b" * 64},
        "reference_files": {
            "reference_fna": {"sha256": "c" * 64},
            "reference_fai": {"sha256": "d" * 64},
            "crosswalk_tsv": {"sha256": "e" * 64},
        },
        "bowtie2_index": {"files": index_files},
        "checks": {
            "all_required_sources_present_hashed_and_count_verified": True,
            "exact_5438_mag_ids_contigs_bases_and_fasta_byte_hashes": True,
            "unique_reference_bam_headers": True,
            "reference_fasta_fai_crosswalk_sha_header_consistency": True,
            "six_actual_index_files_hashed_and_builder_receipt_matched": True,
            "full_focal_target_alias_coverage": True,
        },
    }
    val_path = ref_dir / "stage4_reference_validation_receipt.json"
    val_path.write_text(json.dumps(val_receipt, indent=2, sort_keys=True) + "\n")
    val_sha = stream_sha256(val_path)
    (ref_dir / "REFERENCE_VALIDATION_DONE").write_text(json.dumps({"receipt_sha256": val_sha}) + "\n")
    return {
        "validation_receipt_sha256": val_sha,
        "reference_fna_sha256": "c" * 64,
        "crosswalk_sha256": "e" * 64,
        "build_receipt_sha256": "b" * 64,
        "bowtie2_index_sha256": {row["name"]: row["sha256"] for row in index_files},
    }


def test_actual_stage2_votu_context_and_pair_ledger_integration():
    """Verify real workspace pair ledger, boundary config, submission scope, and Stage 2 vOTU metadata."""
    samples, payload = raw_qc.load_expected_samples(
        resolve_test_input(
            EXP3_ROOT / "provenance/stage4_competitive_pipeline_submission.json",
            EXP3_ROOT / "inputs/stage4_competitive_pipeline_submission.json",
        )
    )
    assert len(samples) == 14
    assert samples == [
        "41PRr2", "V60PRr2", "55Dr2", "9Sr2", "55Sr3",
        "S49Dr1", "48Sr2", "F48Sr2", "V48Dr1", "S48PRr1",
        "41Sr3", "10Dr2", "43PRr3", "44Dr2",
    ]
    assert payload["job_validate_ref"] == "15905724"
    assert payload["job_map_array"] == "15905725"

    pairs, reviewed, boundary = load_analysis_pairs(
        resolve_test_input(
            EXP3_ROOT / "results/stage3_tier1_pairs_individual_evidence_audit.tsv",
            EXP3_ROOT / "inputs/stage3_tier1_pairs_individual_evidence_audit.tsv",
        ),
        resolve_test_input(
            EXP3_ROOT / "inputs/stage4_boundary_sensitive_intervals.json",
        ),
    )
    assert len(reviewed) == 6
    assert len(boundary) == 1
    assert len(pairs) == 7

    pair_meta = raw_qc.build_pair_metadata_by_sha(pairs)
    assert len(pair_meta) == 7

    votu_ctx = raw_qc.load_stage2_votu_context(
        resolve_test_input(
            EXP3_ROOT / "results/stage2_target_viral_identity_table.tsv",
            EXP3_ROOT / "inputs/stage2_target_viral_identity_table.tsv",
        ),
        pairs,
    )
    assert len(votu_ctx) == 7

    # 9Sr2::9Sr2_284237 is a non-representative member pointing to 55Sr3::55Sr3_129891
    sha_9sr2 = "96ff8fd33cdd8a3b3dd1db91ba5645a14797b63723f26e1f5ff2d1be838ef32d"
    assert votu_ctx[sha_9sr2]["canonical_annotated_target_id"] == "9Sr2::9Sr2_284237"
    assert votu_ctx[sha_9sr2]["stage2_is_votu_cluster_rep"] == 0
    assert votu_ctx[sha_9sr2]["stage2_cluster_votu_rep_id"] == "55Sr3::55Sr3_129891"
    assert "context_only_not_demonstrated_cause_of_ambiguity" in votu_ctx[sha_9sr2]["votu_metadata_role"]

    # Exact-sequence aliases 48Sr2::48Sr2_190170 (is_rep=0) and coasm_48_S::k127_1210458 (is_rep=1)
    sha_48sr2 = "ca78a3b6ce5eca9cb6608e7bd98477ff1c6b486ee21ba2ace1d45437ece707f1"
    assert votu_ctx[sha_48sr2]["stage2_is_votu_cluster_rep"] == 0
    assert votu_ctx[sha_48sr2]["stage2_cluster_votu_rep_id"] == "48PRr2::48PRr2_1185466"
    assert (
        votu_ctx[sha_48sr2]["stage2_alias_votu_metadata"]
        == "48Sr2::48Sr2_190170:rep=48PRr2::48PRr2_1185466:is_rep=0;"
        "coasm_48_S::k127_1210458:rep=coasm_48_S::k127_1210458:is_rep=1"
    )


def test_authoritative_9sr2_log_matches_benchmark():
    """Verify the root authoritative 9Sr2 samtools coverage log parses with actual samtools schema."""
    log_path = resolve_test_input(
        EXP3_ROOT / "inputs/exp3_9Sr2_source_target_unfiltered_coverage.log",
        EXP3_ROOT.parent / "review-20261109/exp3_9Sr2_source_target_unfiltered_coverage.log",
    )
    assert log_path.is_file()
    frame = compute_coverage_metrics(log_path.read_text(), {"VIRUS_TARGET__9Sr2::9Sr2_284237": 56873})
    assert len(frame) == 1
    row = frame.iloc[0]
    assert row["rname"] == raw_qc.AUTHORITATIVE_9SR2_BENCHMARK["reference_contig"]
    assert int(row["startpos"]) == 1
    assert int(row["endpos"]) == 56873
    assert int(row["numreads"]) == 12065
    assert int(row["covbases"]) == 56860
    assert float(row["coverage"]) == pytest.approx(99.9771)
    assert float(row["meandepth"]) == pytest.approx(30.919)
    assert float(row["meanbaseq"]) == pytest.approx(35.6)
    assert float(row["meanmapq"]) == pytest.approx(1.51)


def test_reference_validation_contract_and_mapping_metrics_fail_closed(tmp_path: Path):
    """Verify missing/mismatched REFERENCE_VALIDATION_DONE and missing/nonfinite mapping_metrics fail closed."""
    ref_dir = tmp_path / "ref"
    p_sha = "1" * 64
    b_sha = "2" * 64

    with pytest.raises(FileNotFoundError, match="reference validation receipt missing"):
        raw_qc.verify_reference_validation_contract(ref_dir, p_sha, b_sha)

    make_valid_reference_validation_fixture(ref_dir, p_sha, b_sha)
    verified = raw_qc.verify_reference_validation_contract(ref_dir, p_sha, b_sha)
    assert verified["reference_fna_sha256"] == "c" * 64

    # Corrupt REFERENCE_VALIDATION_DONE sentinel hash
    done_path = ref_dir / "REFERENCE_VALIDATION_DONE"
    good_done = done_path.read_text()
    done_path.write_text(json.dumps({"receipt_sha256": "0" * 64}) + "\n")
    with pytest.raises(ValueError, match="REFERENCE_VALIDATION_DONE receipt_sha256 mismatch"):
        raw_qc.verify_reference_validation_contract(ref_dir, p_sha, b_sha)

    # Missing REFERENCE_VALIDATION_DONE sentinel
    done_path.unlink()
    with pytest.raises(FileNotFoundError, match="REFERENCE_VALIDATION_DONE sentinel missing"):
        raw_qc.verify_reference_validation_contract(ref_dir, p_sha, b_sha)
    done_path.write_text(good_done)

    # mapping_metrics validation: missing key or nonfinite value must fail closed (never default to 0)
    valid_metrics = {
        "raw_total_records": 100,
        "raw_primary_mapped_records": 80,
        "filtered_high_quality_mapq20_reads": 40,
        "filtered_fraction_of_primary": 0.5,
        "coverage_references_verified": 3173693,
        "filter": "MAPQ>=20; exclude UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
    }
    parsed = raw_qc.validate_sample_mapping_metrics(valid_metrics, "S1")
    assert parsed["raw_total_records"] == 100

    for key in valid_metrics:
        missing = dict(valid_metrics)
        missing.pop(key)
        with pytest.raises(ValueError, match="missing required keys"):
            raw_qc.validate_sample_mapping_metrics(missing, "S1")

    bad_frac = dict(valid_metrics, filtered_fraction_of_primary=math.nan)
    with pytest.raises(ValueError, match="non-finite filtered_fraction_of_primary"):
        raw_qc.validate_sample_mapping_metrics(bad_frac, "S1")

    bad_bool = dict(valid_metrics, raw_total_records=True)
    with pytest.raises(ValueError, match="expected integer"):
        raw_qc.validate_sample_mapping_metrics(bad_bool, "S1")


@pytest.mark.skipif(
    subprocess.run(["bash", "-c", "command -v samtools"], capture_output=True).returncode != 0,
    reason="samtools required",
)
def test_actual_samtools_raw_vs_mapq20_fixture_and_flags(tmp_path: Path):
    """Exercise actual samtools coverage with --ff 3844 on MAPQ mix, zero, duplicate, supplementary, and missing."""
    sam = tmp_path / "mix.sam"
    seq60 = "A" * 60
    qual60 = "I" * 60
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:VIRUS_TARGET__9Sr2::9Sr2_284237\tLN:100",
        "@SQ\tSN:VIRUS_TARGET__V60PRr2::V60PRr2_715976\tLN:100",
        "@SQ\tSN:VIRUS_TARGET__55Dr2::55Dr2_75333\tLN:100",
    ]
    # 20 low-MAPQ primary reads (MAPQ=1) + 1 high-MAPQ read (MAPQ=40) on 9Sr2 target -> raw >10x/60%, filtered 0.6x/60%
    for idx in range(20):
        lines.append(f"low_{idx}\t0\tVIRUS_TARGET__9Sr2::9Sr2_284237\t1\t1\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    lines.append(f"high_1\t0\tVIRUS_TARGET__9Sr2::9Sr2_284237\t1\t40\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    # Flagged records that --ff 3844 must exclude even in raw (-q 0) mode:
    # 256=SECONDARY, 512=QCFAIL, 1024=DUP, 2048=SUPPLEMENTARY, 4=UNMAP
    for flag in [256, 512, 1024, 2048]:
        lines.append(f"flagged_{flag}\t{flag}\tVIRUS_TARGET__9Sr2::9Sr2_284237\t41\t60\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    lines.append(f"unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t{seq60}\t{qual60}")

    # 20 high-MAPQ primary reads on V60 target -> both raw and filtered reach >=10x and >=50%
    for idx in range(20):
        lines.append(f"v60_{idx}\t0\tVIRUS_TARGET__V60PRr2::V60PRr2_715976\t1\t42\t60M\t*\t0\t0\t{seq60}\t{qual60}")

    # 55Dr2 has 0 reads -> measured_zero_assigned_reads
    sam.write_text("\n".join(lines) + "\n")

    raw_bam = tmp_path / "raw.bam"
    filt_bam = tmp_path / "filt.bam"
    run_cmd(["samtools", "sort", "-o", str(raw_bam), str(sam)])
    run_cmd(["samtools", "index", str(raw_bam)])
    run_cmd(["samtools", "view", "-b", "-q", "20", "-F", "3844", "-o", str(filt_bam), str(raw_bam)])
    run_cmd(["samtools", "index", str(filt_bam)])

    # Check 9Sr2 target: flagged records at 41..100 must NOT inflate covbases beyond 60 or numreads beyond 21
    raw_9s, _, _ = raw_qc.run_samtools_region_coverage(raw_bam, "VIRUS_TARGET__9Sr2::9Sr2_284237", 0, 100)
    filt_9s, _, _ = raw_qc.run_samtools_region_coverage(filt_bam, "VIRUS_TARGET__9Sr2::9Sr2_284237", 20, 100)
    assert raw_9s["numreads"] == 21
    assert raw_9s["covbases"] == 60
    assert raw_9s["meandepth"] == pytest.approx(12.6)
    assert filt_9s["numreads"] == 1
    assert filt_9s["covbases"] == 60
    assert filt_9s["meandepth"] == pytest.approx(0.6)

    qc_class, reason, qual = raw_qc.classify_recruitment_qc(
        raw_numreads=raw_9s["numreads"],
        raw_meandepth=raw_9s["meandepth"],
        raw_breadth_pct=60.0,
        filtered_numreads=filt_9s["numreads"],
        filtered_meandepth=filt_9s["meandepth"],
        filtered_breadth_pct=60.0,
        qualifies_primary=0,
        analysis_role="primary",
        source_category_uniquely_viral=True,
    )
    assert qc_class == "high_raw_support_10x_50pct_lost_at_mapq20"
    assert reason == "high_raw_support_lost_at_mapq20_exact_reference_attribution_unresolved"
    assert "exact-reference attribution unresolved; compatible with competitive ambiguity or weak alignment discrimination" in qual
    assert "due to" not in qual

    # Check V60 target: primary eligible
    raw_v60, _, _ = raw_qc.run_samtools_region_coverage(raw_bam, "VIRUS_TARGET__V60PRr2::V60PRr2_715976", 0, 100)
    filt_v60, _, _ = raw_qc.run_samtools_region_coverage(filt_bam, "VIRUS_TARGET__V60PRr2::V60PRr2_715976", 20, 100)
    assert raw_v60["numreads"] == 20 and filt_v60["numreads"] == 20
    qc_v60, _, _ = raw_qc.classify_recruitment_qc(
        raw_numreads=20, raw_meandepth=12.0, raw_breadth_pct=60.0,
        filtered_numreads=20, filtered_meandepth=12.0, filtered_breadth_pct=60.0,
        qualifies_primary=1, analysis_role="primary", source_category_uniquely_viral=True,
    )
    assert qc_v60 == "primary_mapq20_10x_50pct_eligible"

    # Check 55Dr2 target: measured zero assigned reads
    raw_zero, _, _ = raw_qc.run_samtools_region_coverage(raw_bam, "VIRUS_TARGET__55Dr2::55Dr2_75333", 0, 100)
    assert raw_zero["numreads"] == 0 and raw_zero["covbases"] == 0 and raw_zero["meandepth"] == 0.0
    qc_zero, reason_zero, _ = raw_qc.classify_recruitment_qc(
        raw_numreads=0, raw_meandepth=0.0, raw_breadth_pct=0.0,
        filtered_numreads=0, filtered_meandepth=0.0, filtered_breadth_pct=0.0,
        qualifies_primary=0, analysis_role="primary", source_category_uniquely_viral=True,
    )
    assert qc_zero == "measured_zero_assigned_reads"
    assert reason_zero == "measured_zero_assigned_reads_to_reference_not_proved_biological_absence"

    # Wrong expected length or missing reference fails closed
    with pytest.raises(ValueError, match="length mismatch"):
        raw_qc.run_samtools_region_coverage(raw_bam, "VIRUS_TARGET__55Dr2::55Dr2_75333", 0, 999)
    with pytest.raises(RuntimeError, match="samtools coverage failed"):
        raw_qc.run_samtools_region_coverage(raw_bam, "VIRUS_TARGET__NONEXISTENT", 0, 100)


@pytest.mark.skipif(
    subprocess.run(["bash", "-c", "command -v samtools"], capture_output=True).returncode != 0,
    reason="samtools required",
)
def test_full_14x7_pipeline_and_fail_closed_corruptions(tmp_path: Path):
    """Build a complete 14-sample x 7-target fixture with actual BAMs and test end-to-end + corruptions."""
    pair_ledger = resolve_test_input(
        EXP3_ROOT / "results/stage3_tier1_pairs_individual_evidence_audit.tsv",
        EXP3_ROOT / "inputs/stage3_tier1_pairs_individual_evidence_audit.tsv",
    )
    boundary_config = resolve_test_input(
        EXP3_ROOT / "inputs/stage4_boundary_sensitive_intervals.json",
    )
    stage2_table = resolve_test_input(
        EXP3_ROOT / "results/stage2_target_viral_identity_table.tsv",
        EXP3_ROOT / "inputs/stage2_target_viral_identity_table.tsv",
    )
    submission_receipt = resolve_test_input(
        EXP3_ROOT / "provenance/stage4_competitive_pipeline_submission.json",
        EXP3_ROOT / "inputs/stage4_competitive_pipeline_submission.json",
    )

    samples, _ = raw_qc.load_expected_samples(submission_receipt)
    pairs, _, _ = load_analysis_pairs(pair_ledger, boundary_config)
    pair_meta = raw_qc.build_pair_metadata_by_sha(pairs)

    pair_ledger_sha = stream_sha256(pair_ledger)
    boundary_config_sha = stream_sha256(boundary_config)

    ref_dir = tmp_path / "ref"
    ref_binding = make_valid_reference_validation_fixture(ref_dir, pair_ledger_sha, boundary_config_sha)

    map_dir = tmp_path / "map"
    map_dir.mkdir()

    targets_sorted = [pair_meta[sha] for sha in sorted(pair_meta)]
    sq_lines = ["@HD\tVN:1.6\tSO:coordinate"]
    for t in targets_sorted:
        sq_lines.append(f"@SQ\tSN:VIRUS_TARGET__{t['entity_id']}\tLN:100")

    seq60 = "G" * 60
    qual60 = "I" * 60
    sam_records = list(sq_lines)
    # Give 9Sr2::9Sr2_284237 20 low-MAPQ reads + 1 high-MAPQ read
    for i in range(20):
        sam_records.append(f"r_{i}\t0\tVIRUS_TARGET__9Sr2::9Sr2_284237\t1\t1\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    sam_records.append(f"r_hi\t0\tVIRUS_TARGET__9Sr2::9Sr2_284237\t1\t35\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    # Give 55Dr2::55Dr2_75333 1 low-MAPQ read + 1 high-MAPQ read (low total support with attrition)
    sam_records.append(f"d_lo1\t0\tVIRUS_TARGET__55Dr2::55Dr2_75333\t1\t2\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    sam_records.append(f"d_hi1\t0\tVIRUS_TARGET__55Dr2::55Dr2_75333\t1\t30\t60M\t*\t0\t0\t{seq60}\t{qual60}")

    template_sam = tmp_path / "template.sam"
    template_sam.write_text("\n".join(sam_records) + "\n")
    template_raw = tmp_path / "template_raw.bam"
    template_filt = tmp_path / "template_filt.bam"
    run_cmd(["samtools", "sort", "-o", str(template_raw), str(template_sam)])
    run_cmd(["samtools", "index", str(template_raw)])
    run_cmd(["samtools", "view", "-b", "-q", "20", "-F", "3844", "-o", str(template_filt), str(template_raw)])
    run_cmd(["samtools", "index", str(template_filt)])
    cov_stdout = run_cmd(["samtools", "coverage", str(template_filt)]).stdout
    cov_by_rname = compute_coverage_metrics(cov_stdout).set_index("rname").to_dict(orient="index")

    for idx, sample in enumerate(samples):
        sdir = map_dir / sample
        sdir.mkdir()
        raw_bam = sdir / f"{sample}_raw.bam"
        filt_bam = sdir / f"{sample}_mapq20.bam"
        cov_tsv = sdir / f"{sample}_filtered_coverage.tsv"
        eval_tsv = sdir / f"{sample}_focal_targets_evaluation.tsv"

        raw_bam.write_bytes(template_raw.read_bytes())
        Path(str(raw_bam) + ".bai").write_bytes(Path(str(template_raw) + ".bai").read_bytes())
        filt_bam.write_bytes(template_filt.read_bytes())
        Path(str(filt_bam) + ".bai").write_bytes(Path(str(template_filt) + ".bai").read_bytes())
        cov_tsv.write_text(cov_stdout)

        focal_evals = []
        for sha in sorted(pair_meta):
            pm = pair_meta[sha]
            rname = f"VIRUS_TARGET__{pm['entity_id']}"
            c_entry = cov_by_rname[rname]
            reads = int(c_entry["numreads"])
            breadth = 100.0 * int(c_entry["covbases"]) / 100
            depth = float(c_entry["meandepth"])
            is_boundary = pm["analysis_role"] == "boundary_sensitive_provisional"
            unique = not is_boundary
            focal_evals.append({
                "sample": sample,
                "entity_id": pm["entity_id"],
                "entity_type": "virus_target",
                "entity_sha256": sha,
                "raw_oriented_target_sha256": sha,
                "canonical_sha256": sha,
                "target_aliases": ";".join(pm["target_aliases"]),
                "host_mags": ";".join(pm["host_mags"]),
                "reference_contig": rname,
                "reference_orientation_to_ledger": "same_as_reference",
                "length_bp": 100,
                "filtered_reads_mapped": reads,
                "unique_filtered_reads_mapped": reads if unique else 0,
                "ambiguous_filtered_reads_mapped": 0 if unique else reads,
                "horizontal_breadth_pct": breadth,
                "mean_depth": depth,
                "ambiguity_status": (
                    "boundary_sensitive_provisional; whole contig is not claimed viral; see interval rows"
                    if is_boundary else "unique_viral_sequence_aliases_collapsed"
                ),
                "uniquely_attributable": unique,
                "analysis_role": pm["analysis_role"],
                "qualifies_10x_50pct_instrain": 0,
                "instrain_action": "not_eligible_coverage_ambiguity_or_sensitivity",
                "individual_status": pm["individual_virality_status"],
            })
        pd.DataFrame(focal_evals).to_csv(eval_tsv, sep="\t", index=False)

        receipt = {
            "status": "completed_verified",
            "schema_version": 3,
            "stage": "stage4_map_sample_competitive",
            "sample": sample,
            "execution_environment": {"slurm_job_id": "15905725", "slurm_array_task_id": str(idx)},
            "inputs": {
                "pair_ledger_sha256": pair_ledger_sha,
                "boundary_config_sha256": boundary_config_sha,
            },
            "reference_binding": ref_binding,
            "mapping_metrics": {
                "raw_total_records": 23,
                "raw_primary_mapped_records": 23,
                "filtered_high_quality_mapq20_reads": 2,
                "filtered_fraction_of_primary": 2 / 23,
                "coverage_references_verified": 7,
                "filter": "MAPQ>=20; exclude UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
            },
            "focal_evaluations": focal_evals,
            "outputs": {
                "raw_bam": str(raw_bam),
                "raw_bam_bytes": raw_bam.stat().st_size,
                "raw_bam_sha256": stream_sha256(raw_bam),
                "filtered_bam": str(filt_bam),
                "filtered_bam_sha256": stream_sha256(filt_bam),
                "coverage_tsv": str(cov_tsv),
                "coverage_tsv_sha256": stream_sha256(cov_tsv),
                "focal_evaluation_tsv": str(eval_tsv),
                "focal_evaluation_tsv_sha256": stream_sha256(eval_tsv),
            },
        }
        rpath = sdir / f"{sample}_mapping_receipt.json"
        rpath.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
        (sdir / "MAPPING_DONE").write_text(json.dumps({"receipt_sha256": stream_sha256(rpath)}) + "\n")

    out_dir = tmp_path / "out"
    summary = raw_qc.run_raw_vs_mapq20_qc(
        map_dir=map_dir,
        ref_dir=ref_dir,
        pair_ledger=pair_ledger,
        boundary_config=boundary_config,
        stage2_identity_table=stage2_table,
        submission_receipt=submission_receipt,
        out_dir=out_dir,
        threads=4,
        verify_9sr2_benchmark=False,
    )
    assert summary["status"] == "completed_verified"
    assert summary["scope_counts"]["n_samples_expected_and_verified"] == 14
    assert summary["scope_counts"]["n_unique_viral_targets"] == 7
    assert summary["scope_counts"]["n_sample_target_rows"] == 98

    df98 = pd.read_csv(out_dir / "stage4_raw_vs_mapq20_focal_coverage_98rows.tsv", sep="\t")
    assert len(df98) == 98
    assert set(df98["recruitment_qc_class"]) == {
        "high_raw_support_10x_50pct_lost_at_mapq20",
        "low_total_read_support_with_mapq20_attrition",
        "measured_zero_assigned_reads",
    }

    # Fail-closed test 1: corrupt MAPPING_DONE hash binding
    done_9s = map_dir / "9Sr2" / "MAPPING_DONE"
    orig_done = done_9s.read_text()
    done_9s.write_text(json.dumps({"receipt_sha256": "0" * 64}) + "\n")
    with pytest.raises(ValueError, match="MAPPING_DONE receipt_sha256 mismatch"):
        raw_qc.run_raw_vs_mapq20_qc(
            map_dir=map_dir, ref_dir=ref_dir, pair_ledger=pair_ledger,
            boundary_config=boundary_config, stage2_identity_table=stage2_table,
            submission_receipt=submission_receipt, out_dir=out_dir, threads=2,
            verify_9sr2_benchmark=False,
        )
    done_9s.write_text(orig_done)

    # Fail-closed test 2: missing MAPPING_DONE sentinel
    done_9s.unlink()
    with pytest.raises(FileNotFoundError, match="missing MAPPING_DONE"):
        raw_qc.run_raw_vs_mapq20_qc(
            map_dir=map_dir, ref_dir=ref_dir, pair_ledger=pair_ledger,
            boundary_config=boundary_config, stage2_identity_table=stage2_table,
            submission_receipt=submission_receipt, out_dir=out_dir, threads=2,
            verify_9sr2_benchmark=False,
        )
    done_9s.write_text(orig_done)

    # Fail-closed test 3: corrupt raw BAM file bytes without updating receipt
    raw_9s_path = map_dir / "9Sr2" / "9Sr2_raw.bam"
    orig_raw_bytes = raw_9s_path.read_bytes()
    raw_9s_path.write_bytes(orig_raw_bytes + b"corrupt")
    with pytest.raises(ValueError, match="receipted output hash mismatch|raw_bam_bytes mismatch"):
        raw_qc.run_raw_vs_mapq20_qc(
            map_dir=map_dir, ref_dir=ref_dir, pair_ledger=pair_ledger,
            boundary_config=boundary_config, stage2_identity_table=stage2_table,
            submission_receipt=submission_receipt, out_dir=out_dir, threads=2,
            verify_9sr2_benchmark=False,
        )
    raw_9s_path.write_bytes(orig_raw_bytes)
