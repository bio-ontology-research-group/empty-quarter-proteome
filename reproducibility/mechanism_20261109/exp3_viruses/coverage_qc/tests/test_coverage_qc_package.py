from __future__ import annotations

import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

PKG_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PKG_ROOT))
sys.path.insert(0, str(PKG_ROOT / "code"))
sys.path.insert(0, str(PKG_ROOT / "tests"))

import stage4_raw_vs_mapq20_focal_coverage as raw_qc  # noqa: E402
from stage4_competitive_common import load_analysis_pairs, stream_sha256  # noqa: E402
from stage4_map_sample_competitive import compute_coverage_metrics  # noqa: E402
from test_stage4_raw_vs_mapq20 import make_valid_reference_validation_fixture, run_cmd  # noqa: E402
import verify_coverage_qc  # noqa: E402


def rebuild_sha256_manifest(pkg_dir: Path) -> None:
    rows = []
    for p in sorted(pkg_dir.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(pkg_dir).as_posix()
        if rel == "manifests/sha256_manifest.tsv" or ".pytest" in rel or "__pycache__" in rel:
            continue
        rows.append({"relative_path": rel, "bytes": p.stat().st_size, "sha256": stream_sha256(p)})
    pd.DataFrame(rows).to_csv(pkg_dir / "manifests/sha256_manifest.tsv", sep="\t", index=False)


def test_real_shaheen_reference_validation_and_available_sample_receipts():
    """Validate the real Shaheen reference validation receipt/sentinel and all available real sample receipts/evaluations."""
    in_dir = PKG_ROOT / "inputs"
    pair_ledger = in_dir / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    boundary_cfg = in_dir / "stage4_boundary_sensitive_intervals.json"
    sub_receipt = in_dir / "stage4_competitive_pipeline_submission.json"

    samples, _ = raw_qc.load_expected_samples(sub_receipt)
    pairs, _, _ = load_analysis_pairs(pair_ledger, boundary_cfg)
    pair_meta_by_sha = raw_qc.build_pair_metadata_by_sha(pairs)

    pair_ledger_sha = stream_sha256(pair_ledger)
    boundary_cfg_sha = stream_sha256(boundary_cfg)
    expected_ref_validation = raw_qc.verify_reference_validation_contract(
        in_dir / "stage4_competitive_reference",
        pair_ledger_sha,
        boundary_cfg_sha,
    )
    assert expected_ref_validation["validation_receipt_sha256"] == "589d5a28598ea36f5237f198f687dc6bc9cbf02434c0717956419625eece758f"

    map_dir = in_dir / "stage4_competitive_mapping"
    verified_receipts = 0
    verified_eval_tsvs = 0
    for sample in samples:
        s_dir = map_dir / sample
        e_path = s_dir / f"{sample}_focal_targets_evaluation.tsv"
        if e_path.is_file():
            ev_df = pd.read_csv(e_path, sep="\t")
            ev_virus = ev_df[ev_df["entity_type"] == "virus_target"]
            assert len(ev_virus) == 7
            assert set(ev_virus["entity_sha256"].astype(str)) == set(pair_meta_by_sha)
            verified_eval_tsvs += 1

        r_path = s_dir / f"{sample}_mapping_receipt.json"
        d_path = s_dir / "MAPPING_DONE"
        if r_path.is_file() and d_path.is_file():
            r_sha = stream_sha256(r_path)
            assert json.loads(d_path.read_text()).get("receipt_sha256") == r_sha
            r_payload = json.loads(r_path.read_text())
            assert r_payload["sample"] == sample
            assert r_payload["status"] == "completed_verified"
            assert r_payload["schema_version"] == 3
            raw_qc.validate_sample_mapping_metrics(r_payload["mapping_metrics"], sample)
            for k in ["validation_receipt_sha256", "reference_fna_sha256", "crosswalk_sha256", "build_receipt_sha256"]:
                assert r_payload["reference_binding"][k] == expected_ref_validation[k]
            assert r_payload["reference_binding"]["bowtie2_index_sha256"] == expected_ref_validation["bowtie2_index_sha256"]
            assert stream_sha256(e_path) == r_payload["outputs"]["focal_evaluation_tsv_sha256"]
            verified_receipts += 1

    assert verified_eval_tsvs == 14
    assert verified_receipts >= 13


@pytest.mark.skipif(
    subprocess.run(["bash", "-c", "command -v samtools"], capture_output=True).returncode != 0,
    reason="samtools required",
)
def test_synthetic_package_fixture_passes_all_7_gates_and_fails_closed(tmp_path: Path):
    """Build a complete synthetic 14x7 coverage_qc package fixture using actual producer output and verify all 7 gates."""
    pkg = tmp_path / "coverage_qc_synth"
    for d in ["code", "inputs/stage4_competitive_reference", "inputs/stage4_competitive_mapping", "results", "manifests", "tests"]:
        (pkg / d).mkdir(parents=True, exist_ok=True)

    for fname in [
        "stage4_raw_vs_mapq20_focal_coverage.py",
        "sbatch_stage4_raw_vs_mapq20.sh",
        "stage4_competitive_common.py",
        "stage4_map_sample_competitive.py",
    ]:
        shutil.copy2(PKG_ROOT / "code" / fname, pkg / "code" / fname)
    shutil.copy2(PKG_ROOT / "tests/test_stage4_raw_vs_mapq20.py", pkg / "tests/test_stage4_raw_vs_mapq20.py")

    for in_name in [
        "stage3_tier1_pairs_individual_evidence_audit.tsv",
        "stage4_boundary_sensitive_intervals.json",
        "stage2_target_viral_identity_table.tsv",
        "stage4_competitive_pipeline_submission.json",
        "exp3_9Sr2_source_target_unfiltered_coverage.log",
    ]:
        shutil.copy2(PKG_ROOT / "inputs" / in_name, pkg / "inputs" / in_name)

    pair_ledger = pkg / "inputs/stage3_tier1_pairs_individual_evidence_audit.tsv"
    boundary_cfg = pkg / "inputs/stage4_boundary_sensitive_intervals.json"
    stage2_table = pkg / "inputs/stage2_target_viral_identity_table.tsv"
    sub_receipt = pkg / "inputs/stage4_competitive_pipeline_submission.json"

    samples, _ = raw_qc.load_expected_samples(sub_receipt)
    pairs, _, _ = load_analysis_pairs(pair_ledger, boundary_cfg)
    pair_meta = raw_qc.build_pair_metadata_by_sha(pairs)

    pair_ledger_sha = stream_sha256(pair_ledger)
    boundary_cfg_sha = stream_sha256(boundary_cfg)

    ref_dir = pkg / "inputs/stage4_competitive_reference"
    ref_binding = make_valid_reference_validation_fixture(ref_dir, pair_ledger_sha, boundary_cfg_sha)

    # Use length 103 (prime) so 100*covbases/length_bp differs from 4-decimal rounded samtools coverage column!
    bam_work = tmp_path / "bam_work"
    bam_work.mkdir()
    targets_sorted = [pair_meta[sha] for sha in sorted(pair_meta)]
    sq_lines = ["@HD\tVN:1.6\tSO:coordinate"]
    for t in targets_sorted:
        sq_lines.append(f"@SQ\tSN:VIRUS_TARGET__{t['entity_id']}\tLN:103")

    seq60 = "G" * 60
    qual60 = "I" * 60
    sam_records = list(sq_lines)
    for i in range(20):
        sam_records.append(f"r_{i}\t0\tVIRUS_TARGET__9Sr2::9Sr2_284237\t1\t1\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    sam_records.append(f"r_hi\t0\tVIRUS_TARGET__9Sr2::9Sr2_284237\t1\t35\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    sam_records.append(f"d_lo1\t0\tVIRUS_TARGET__55Dr2::55Dr2_75333\t1\t2\t60M\t*\t0\t0\t{seq60}\t{qual60}")
    sam_records.append(f"d_hi1\t0\tVIRUS_TARGET__55Dr2::55Dr2_75333\t1\t30\t60M\t*\t0\t0\t{seq60}\t{qual60}")

    template_sam = bam_work / "template.sam"
    template_sam.write_text("\n".join(sam_records) + "\n")
    template_raw = bam_work / "template_raw.bam"
    template_filt = bam_work / "template_filt.bam"
    run_cmd(["samtools", "sort", "-o", str(template_raw), str(template_sam)])
    run_cmd(["samtools", "index", str(template_raw)])
    run_cmd(["samtools", "view", "-b", "-q", "20", "-F", "3844", "-o", str(template_filt), str(template_raw)])
    run_cmd(["samtools", "index", str(template_filt)])
    cov_stdout = run_cmd(["samtools", "coverage", str(template_filt)]).stdout
    cov_by_rname = compute_coverage_metrics(cov_stdout).set_index("rname").to_dict(orient="index")

    map_dir = pkg / "inputs/stage4_competitive_mapping"
    for idx, sample in enumerate(samples):
        sdir = map_dir / sample
        sdir.mkdir()
        raw_bam = bam_work / f"{sample}_raw.bam"
        filt_bam = bam_work / f"{sample}_mapq20.bam"
        cov_tsv = bam_work / f"{sample}_filtered_coverage.tsv"
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
            breadth = 100.0 * int(c_entry["covbases"]) / 103
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
                "length_bp": 103,
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
        focal_evals.append({
            "sample": sample,
            "entity_id": "44Dr2::44Dr2_4002037|catalog_viral_core",
            "entity_type": "boundary_interval",
            "entity_sha256": "6293eefbf49a4cf8cb58cb6aa2c6dd877e3e81450917bd018483580bf3028bb1",
            "raw_oriented_target_sha256": "",
            "canonical_sha256": "",
            "target_aliases": "",
            "host_mags": "43PRr3__self_SemiBin_469",
            "reference_contig": "VIRUS_TARGET__44Dr2::44Dr2_4002037",
            "reference_orientation_to_ledger": "",
            "length_bp": 50,
            "filtered_reads_mapped": 0,
            "unique_filtered_reads_mapped": 0,
            "ambiguous_filtered_reads_mapped": 0,
            "horizontal_breadth_pct": 0.0,
            "mean_depth": 0.0,
            "ambiguity_status": "nested_catalog_core_competitor; low MAPQ-filtered coverage is not absence",
            "uniquely_attributable": False,
            "analysis_role": "boundary_sensitive_provisional",
            "qualifies_10x_50pct_instrain": 0,
            "instrain_action": "not_eligible_boundary_sensitive_provisional",
            "individual_status": "boundary_sensitive_provisional_not_primary_host_link",
        })
        focal_evals.append({
            "sample": sample,
            "entity_id": "41PRr2__single_SemiBin_236",
            "entity_type": "host_mag",
            "entity_sha256": "",
            "raw_oriented_target_sha256": "",
            "canonical_sha256": "",
            "target_aliases": "",
            "host_mags": "",
            "reference_contig": "multiple_original_contigs_see_crosswalk",
            "reference_orientation_to_ledger": "",
            "length_bp": 3372843,
            "filtered_reads_mapped": 10,
            "unique_filtered_reads_mapped": 10,
            "ambiguous_filtered_reads_mapped": 0,
            "horizontal_breadth_pct": 0.1,
            "mean_depth": 0.01,
            "ambiguity_status": "ownership_derived_from_full_crosswalk",
            "uniquely_attributable": False,
            "analysis_role": "host_context",
            "qualifies_10x_50pct_instrain": 0,
            "instrain_action": "host_context_only_not_profiled",
            "individual_status": "",
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
                "boundary_config_sha256": boundary_cfg_sha,
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

    raw_qc.run_raw_vs_mapq20_qc(
        map_dir=map_dir,
        ref_dir=ref_dir,
        pair_ledger=pair_ledger,
        boundary_config=boundary_cfg,
        stage2_identity_table=stage2_table,
        submission_receipt=sub_receipt,
        out_dir=pkg / "results",
        threads=2,
        verify_9sr2_benchmark=False,
    )

    code_hash_lines = [
        f"{stream_sha256(pkg / rel)}  {rel}"
        for rel in [
            "code/stage4_raw_vs_mapq20_focal_coverage.py",
            "code/sbatch_stage4_raw_vs_mapq20.sh",
            "code/stage4_competitive_common.py",
            "code/stage4_map_sample_competitive.py",
            "tests/test_stage4_raw_vs_mapq20.py",
        ]
    ]
    (pkg / "manifests/stage4_raw_vs_mapq20_code_hashes.tsv").write_text("\n".join(code_hash_lines) + "\n")
    pd.DataFrame([
        {
            "job_id": "15910231",
            "pipeline_stage": "stage4_raw_vs_mapq20_focal_coverage",
            "job_name": "exp3_raw_cov_qc",
            "user": "hohndor",
            "partition": "workq",
            "account": "k10018",
            "nodes": "nid00001",
            "n_tasks": 1,
            "state": "COMPLETED",
            "exit_code": "0:0",
            "elapsed_hms": "00:01:00",
            "start_local": "2026-11-09T22:05:00",
            "end_local": "2026-11-09T22:06:00",
            "submission_context": "Synthetic test execution receipt",
        }
    ]).to_csv(pkg / "manifests/slurm_execution_receipts.tsv", sep="\t", index=False)
    rebuild_sha256_manifest(pkg)

    res = verify_coverage_qc.verify_all(pkg, verify_9sr2_benchmark=False)
    assert res["status"] == "ALL_7_GATES_PASSED"
    assert res["focal_rows"] == 98

    # Fail-closed check 1: Tamper with STAGE4_RAW_VS_MAPQ20_DONE summary_sha256
    done_path = pkg / "results/STAGE4_RAW_VS_MAPQ20_DONE"
    orig_done = done_path.read_text()
    done_path.write_text(json.dumps({"summary_sha256": "0" * 64}) + "\n")
    rebuild_sha256_manifest(pkg)
    with pytest.raises(AssertionError, match="DONE summary_sha256 mismatch"):
        verify_coverage_qc.verify_all(pkg, verify_9sr2_benchmark=False)
    done_path.write_text(orig_done)

    # Fail-closed check 2: Replace empty zero-denominator retained_read_fraction with 0.0
    rows_path = pkg / "results/stage4_raw_vs_mapq20_focal_coverage_98rows.tsv"
    orig_rows = rows_path.read_text()
    df_bad = pd.read_csv(rows_path, sep="\t", keep_default_na=False)
    zero_idx = df_bad.index[df_bad["raw_numreads"].astype(int) == 0][0]
    df_bad.loc[zero_idx, "retained_read_fraction"] = "0.0"
    df_bad.to_csv(rows_path, sep="\t", index=False)
    rebuild_sha256_manifest(pkg)
    with pytest.raises(AssertionError, match="expected empty retained_read_fraction when raw_numreads==0"):
        verify_coverage_qc.verify_all(pkg, verify_9sr2_benchmark=False)
    rows_path.write_text(orig_rows)

    # Fail-closed check 3: Replace precise raw_horizontal_breadth_pct with rounded samtools coverage pct (for length 103, 60/103*100 = 58.25242718... != 58.2524)
    df_round = pd.read_csv(rows_path, sep="\t", keep_default_na=False)
    nonzero_idx = df_round.index[df_round["raw_covbases"].astype(int) > 0][0]
    df_round.loc[nonzero_idx, "raw_horizontal_breadth_pct"] = df_round.loc[nonzero_idx, "raw_coverage_samtools_pct"]
    df_round.to_csv(rows_path, sep="\t", index=False)
    rebuild_sha256_manifest(pkg)
    with pytest.raises(AssertionError, match="raw_horizontal_breadth_pct mismatch"):
        verify_coverage_qc.verify_all(pkg, verify_9sr2_benchmark=False)
    rows_path.write_text(orig_rows)

    # Fail-closed check 4: Duplicate a (sample, raw_oriented_target_sha256) row (keeping 98 rows total)
    df_dup = pd.read_csv(rows_path, sep="\t", keep_default_na=False)
    df_dup.iloc[1] = df_dup.iloc[0]
    df_dup.to_csv(rows_path, sep="\t", index=False)
    rebuild_sha256_manifest(pkg)
    with pytest.raises(AssertionError, match="Duplicate \\(sample, raw_oriented_target_sha256\\) keys"):
        verify_coverage_qc.verify_all(pkg, verify_9sr2_benchmark=False)
    rows_path.write_text(orig_rows)
    rebuild_sha256_manifest(pkg)


@pytest.mark.skipif(
    not (PKG_ROOT / "results/STAGE4_RAW_VS_MAPQ20_DONE").is_file(),
    reason="Real Shaheen results not yet synced into PKG_ROOT/results",
)
def test_real_packaged_coverage_qc_passes_all_7_gates():
    """Verify the real packaged coverage_qc directory passes all 7 verification gates including 9Sr2 benchmark."""
    res = verify_coverage_qc.verify_all(PKG_ROOT, verify_9sr2_benchmark=True)
    assert res["status"] == "ALL_7_GATES_PASSED"
    assert res["focal_rows"] == 98
    assert res["manifest_files"] >= 60
