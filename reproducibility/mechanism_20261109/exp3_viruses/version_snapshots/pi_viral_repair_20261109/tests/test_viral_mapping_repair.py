from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

CODE = Path(__file__).resolve().parents[1] / "code"
sys.path.insert(0, str(CODE))

import audit_target_sequence_collisions_v2 as collision
import stage4_aggregate_competitive_mapping as aggregate
import stage4_map_sample_competitive as mapping
import stage4_structural_hallmark_screen as sequence_followup
from stage4_competitive_common import (
    CROSSWALK_COLUMNS,
    PAIR_COLUMNS,
    VERIFIED,
    UNRESOLVED,
    canonical_sha256,
    load_boundary_candidates,
    load_reviewed_pairs,
    reference_owners,
    resolve_focal_targets,
    sha256_sequence,
    stream_sha256,
    validate_crosswalk_schema,
)


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, check=True, capture_output=True, text=True)


def cw_row(header: str, category: str, entity: str, contig: str, length: int, canonical: str, rep: str) -> dict:
    return {
        "original_header": header,
        "source_category": category,
        "entity_id": entity,
        "contig_id": contig,
        "seq_length": length,
        "canonical_sha256": canonical,
        "is_reference_representative": int(header == rep),
        "representative_header": rep,
        "dedup_reason": "primary_representative" if header == rep else "duplicate_sequence_of_earlier_contig",
        "ambiguity_status": "not_used_for_ownership",
    }


def pair_row(host: str, target: str, canonical_target: str, sha: str, status: str = VERIFIED) -> dict:
    values = {
        "host_mag": host,
        "host_phylum": "p__test",
        "host_sample": host.split("__")[0],
        "has_host_fastq": "True",
        "target_id": target,
        "canonical_annotated_target_id": canonical_target,
        "votu_id": canonical_target,
        "virus_sample": target.split("::")[0],
        "has_vir_fastq": "True",
        "target_seq_sha256": sha,
        "individual_virality_status": status,
        "eligible_for_joint_competitive_mapping": "True",
        "eligible_for_structural_hallmark_followup": "True",
    }
    assert set(values) == set(PAIR_COLUMNS)
    return values


@pytest.mark.skipif(subprocess.run(["bash", "-c", "command -v samtools"], capture_output=True).returncode != 0, reason="samtools required")
def test_actual_samtools_coverage_filter_zero_duplicate_lowmapq_and_histogram(tmp_path: Path):
    """Exercise actual BAM and samtools, not a hand-written coverage TSV mock."""
    sam = tmp_path / "fixture.sam"
    seq50 = "A" * 50
    qual50 = "I" * 50
    sam.write_text(
        "@HD\tVN:1.6\tSO:coordinate\n"
        "@SQ\tSN:HOST__M1__c1\tLN:100\n"
        "@SQ\tSN:VIRUS_TARGET__V1\tLN:80\n"
        "@SQ\tSN:HOST__M1__zero\tLN:60\n"
        f"mapped\t0\tHOST__M1__c1\t1\t60\t50M\t*\t0\t0\t{seq50}\t{qual50}\n"
        f"lowmapq\t0\tHOST__M1__c1\t51\t10\t50M\t*\t0\t0\t{seq50}\t{qual50}\n"
        f"duplicate\t1024\tHOST__M1__c1\t1\t60\t50M\t*\t0\t0\t{seq50}\t{qual50}\n"
        f"viral\t0\tVIRUS_TARGET__V1\t1\t60\t40M\t*\t0\t0\t{'C' * 40}\t{'I' * 40}\n"
        f"unmapped\t4\t*\t0\t0\t*\t*\t0\t0\t{'N' * 20}\t{'I' * 20}\n"
    )
    raw = tmp_path / "raw.bam"
    filtered = tmp_path / "filtered.bam"
    run(["samtools", "view", "-b", "-o", str(raw), str(sam)])
    run(["samtools", "view", "-b", "-q", "20", "-F", "3844", "-o", str(filtered), str(raw)])
    run(["samtools", "quickcheck", str(filtered)])
    run(["samtools", "index", str(filtered)])
    lengths = mapping.bam_reference_lengths(run(["samtools", "view", "-H", str(filtered)]).stdout)
    actual = run(["samtools", "coverage", str(filtered)]).stdout
    frame = mapping.compute_coverage_metrics(actual, lengths)
    assert set(frame["rname"]) == set(lengths)
    by_name = frame.set_index("rname")
    assert by_name.loc["HOST__M1__c1", "numreads"] == 1
    assert by_name.loc["VIRUS_TARGET__V1", "numreads"] == 1
    assert by_name.loc["HOST__M1__zero", "numreads"] == 0
    assert by_name.loc["HOST__M1__zero", "covbases"] == 0
    assert by_name.loc["HOST__M1__zero", "meandepth"] == 0
    histogram = run(["samtools", "coverage", "-m", str(filtered)]).stdout
    with pytest.raises(ValueError, match="unexpected samtools coverage header"):
        mapping.compute_coverage_metrics(histogram, lengths)
    with pytest.raises(ValueError, match="empty output"):
        mapping.compute_coverage_metrics("")
    with pytest.raises(ValueError, match="absent from BAM header"):
        mapping.compute_coverage_metrics(actual, {name: length for name, length in lengths.items() if name != "HOST__M1__zero"})
    depth = run(["samtools", "depth", "-a", "-d", "0", "-Q", "20", "-r", "VIRUS_TARGET__V1:41-80", str(filtered)]).stdout
    interval = mapping.compute_interval_coverage(depth, "VIRUS_TARGET__V1", 41, 80)
    assert interval == {"length_bp": 40, "callable_positions_depth_ge1": 0, "horizontal_breadth_pct": 0.0, "mean_depth": 0.0}
    with pytest.raises(ValueError, match="missing|noncontiguous"):
        mapping.compute_interval_coverage("VIRUS_TARGET__V1\t41\t0\n", "VIRUS_TARGET__V1", 41, 80)


def test_actual_builder_crosswalk_schema_and_alias_hash_domains():
    sequence = "AACCGTAACG"
    from stage4_competitive_common import reverse_complement
    oriented = reverse_complement(sequence)
    raw_sha = sha256_sequence(oriented)
    canonical = canonical_sha256(sequence)
    assert raw_sha != canonical  # distinct domains in this non-palindromic fixture
    crosswalk = pd.DataFrame([
        cw_row("VIRUS_TARGET__V1", "virus_target", "V1", "V1", len(sequence), canonical, "VIRUS_TARGET__V1"),
        cw_row("VIRUS_COASM__V1_alias", "coasm_viruses", "V1_alias", "V1_alias", len(sequence), canonical, "VIRUS_TARGET__V1"),
    ])
    assert list(crosswalk.columns) == CROSSWALK_COLUMNS
    validate_crosswalk_schema(crosswalk)
    pairs = pd.DataFrame([{
        "target_seq_sha256": raw_sha,
        "target_aliases": [["V1", "V1_alias"]],
        "canonical_annotated_target_id": "V1",
    }])
    resolved = resolve_focal_targets(pairs, crosswalk, lambda rep: sequence)
    assert resolved[raw_sha]["reference_contig"] == "VIRUS_TARGET__V1"
    assert resolved[raw_sha]["reference_orientation_to_ledger"] == "reverse_complement_of_reference"
    with pytest.raises(ValueError, match="absent from original-source crosswalk"):
        resolve_focal_targets(pairs, crosswalk.iloc[[0]], lambda rep: sequence)
    with pytest.raises(ValueError, match="raw oriented target SHA mismatch"):
        bad_pairs = pairs.copy()
        bad_pairs.loc[:, "target_seq_sha256"] = "0" * 64
        resolve_focal_targets(bad_pairs, crosswalk, lambda rep: sequence)
    with pytest.raises(ValueError, match="canonical hash mismatch"):
        resolve_focal_targets(pairs, crosswalk, lambda rep: "TTTTTTTTTT")


def test_reviewed_pair_ledger_alias_collapse_and_sensitivity(tmp_path: Path):
    sha1 = "1" * 64
    sha2 = "2" * 64
    rows = [
        pair_row("H1__bin", "V1", "V1", sha1),
        pair_row("H1__bin", "V1_alias", "V1", sha1),
        pair_row("H2__bin", "V2", "V2", sha2, UNRESOLVED),
    ]
    path = tmp_path / "pairs.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    pairs = load_reviewed_pairs(path)
    assert len(pairs) == 2
    primary = pairs[pairs["analysis_role"] == "primary"].iloc[0]
    assert primary["target_aliases"] == ["V1", "V1_alias"]
    assert len(primary["target_aliases"]) == 2
    assert (pairs["analysis_role"] == "sensitivity_unresolved_virality").sum() == 1
    bad = pd.DataFrame(rows).drop(columns=["target_seq_sha256"])
    bad.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="missing required columns"):
        load_reviewed_pairs(path)


def test_host_full_denominator_shared_owners_zero_and_duplicate_reads():
    shas = [f"{i}" * 64 for i in range(1, 5)]
    rows = [
        cw_row("HOST__M1__c1", "host_mag", "M1", "c1", 100, shas[0], "HOST__M1__c1"),
        cw_row("HOST__M1__c1copy", "host_mag", "M1", "c1copy", 100, shas[0], "HOST__M1__c1"),
        cw_row("HOST__M1__shared", "host_mag", "M1", "shared", 200, shas[1], "HOST__M1__shared"),
        cw_row("HOST__M2__shared", "host_mag", "M2", "shared", 200, shas[1], "HOST__M1__shared"),
        cw_row("HOST__M1__viral", "host_mag", "M1", "viral", 80, shas[2], "HOST__M1__viral"),
        cw_row("VIRUS_TARGET__V1", "virus_target", "V1", "V1", 80, shas[2], "HOST__M1__viral"),
        cw_row("HOST__M1__zero", "host_mag", "M1", "zero", 60, shas[3], "HOST__M1__zero"),
    ]
    crosswalk = pd.DataFrame(rows)
    owners = reference_owners(crosswalk)
    coverage = {
        "HOST__M1__c1": {"numreads": 2, "covbases": 100, "meandepth": 10.0},
        "HOST__M1__shared": {"numreads": 4, "covbases": 100, "meandepth": 5.0},
        "HOST__M1__viral": {"numreads": 3, "covbases": 40, "meandepth": 2.0},
        "HOST__M1__zero": {"numreads": 0, "covbases": 0, "meandepth": 0.0},
    }
    host_rows = crosswalk[(crosswalk["source_category"] == "host_mag") & (crosswalk["entity_id"] == "M1")]
    result = mapping.evaluate_host_coverage("M1", host_rows, owners, coverage, 5, 540, "f" * 64)
    assert result["filtered_reads_mapped"] == 9
    assert result["unique_filtered_reads_mapped"] == 2  # c1copy cannot inflate reads
    assert result["ambiguous_filtered_reads_mapped"] == 7
    assert result["length_bp"] == 540
    assert result["unique_mag_bases"] == 260
    assert result["ambiguous_mag_bases"] == 280
    assert result["unique_mag_covbases"] == 200
    with pytest.raises(ValueError, match="unknown is not zero"):
        incomplete = dict(coverage)
        incomplete.pop("HOST__M1__zero")
        mapping.evaluate_host_coverage("M1", host_rows, owners, incomplete, 5, 540, "f" * 64)


def test_collision_truthful_gate_identity_within_source_and_final_mismatch():
    sha = hashlib.sha256(b"ACGT").hexdigest()
    pools = {"p1": sha, "p2": sha}
    assert collision.determine_target_status([], pools, sha, sha, [], []) == "identical_duplicate_across_pools"
    assert collision.determine_target_status(["p3"], pools, sha, sha, [], []) == "unknown_incomplete_source_coverage"
    assert collision.determine_target_status([], pools, "0" * 64, sha, [], []) == "conflict_with_identity_table"
    assert collision.determine_target_status([], pools, sha, "0" * 64, [], []) == "conflict_with_final_merged_fasta"
    assert collision.determine_target_status([], pools, sha, sha, [{"header": "V1"}], []) == "conflict_within_source_repeated_header"
    good_row = {"status": "identical_duplicate_across_pools", "has_complete_source_coverage": True, "matches_identity_table": True, "final_fna_sha": sha, "identity_table_sha": sha}
    assert collision.collision_free_verdict([good_row], [], [])
    assert not collision.collision_free_verdict([good_row], [{"header": "other"}], [])
    assert not collision.collision_free_verdict([good_row], [], [{"header": "V1"}])


def test_missing_instrain_container_and_callable_positions(tmp_path: Path):
    with pytest.raises(RuntimeError, match="container/tool unavailable"):
        mapping.require_instrain_available(tmp_path / "missing.sif", None)
    output = tmp_path / "profile" / "output"
    output.mkdir(parents=True)
    pd.DataFrame([{"scaffold": "VIRUS_TARGET__V1", "length": 100, "breadth": 0.5}]).to_csv(output / "profile_scaffold_info.tsv", sep="\t", index=False)
    (output / "profile_SNVs.tsv").write_text("scaffold\tposition\nVIRUS_TARGET__V1\t10\n")
    report = mapping.collect_instrain_callable_positions(output.parent, "VIRUS_TARGET__V1")
    assert report["callable_positions"] == 50
    assert report["variant_table_rows"] == 1
    assert "not evidence of recombination" in report["variant_interpretation"]


def test_aggregate_receipt_mutation_fails(tmp_path: Path):
    sample = "S1"
    out = tmp_path / sample
    out.mkdir()
    outputs = {}
    for key in ["raw_bam", "filtered_bam", "coverage_tsv", "focal_evaluation_tsv"]:
        path = out / key
        if key == "focal_evaluation_tsv":
            path.write_text("sample\tentity_id\nS1\tV1\n")
        else:
            path.write_bytes(key.encode())
        outputs[key] = str(path)
        outputs[key + "_sha256"] = stream_sha256(path)
    receipt = {
        "status": "completed_verified", "schema_version": 3, "sample": sample,
        "outputs": outputs,
    }
    receipt_path = out / f"{sample}_mapping_receipt.json"
    receipt_path.write_text(json.dumps(receipt))
    (out / "MAPPING_DONE").write_text(json.dumps({"receipt_sha256": stream_sha256(receipt_path)}))
    assert aggregate.verify_sample_receipt(receipt_path)["sample"] == sample
    Path(outputs["focal_evaluation_tsv"]).write_text("sample\tentity_id\nS1\tMUTATED\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        aggregate.verify_sample_receipt(receipt_path)


def test_strict_gene_sequence_id_and_exact_extraction(tmp_path: Path):
    assert sequence_followup.gene_sequence_id("V1|provirus_1_20_3") == "V1|provirus_1_20"
    assert sequence_followup.gene_sequence_id("V10_2") == "V10"
    with pytest.raises(ValueError):
        sequence_followup.gene_sequence_id("V1_without_suffix")
    proteins = tmp_path / "proteins.faa"
    proteins.write_text(">V1_1\nAAAA\n>V10_1\nCCCC\n")
    output = tmp_path / "candidate.faa"
    count, ids = sequence_followup.extract_candidate_proteins({"V1_1"}, proteins, output)
    assert (count, ids) == (1, ["V1_1"])
    with pytest.raises(ValueError, match="exact protein extraction failed"):
        sequence_followup.extract_candidate_proteins({"V1_1", "V1_2"}, proteins, output)


def test_root_boundary_candidate_stays_provisional_with_three_intervals():
    config = Path(__file__).resolve().parents[1] / "inputs/stage4_boundary_sensitive_intervals.json"
    boundary = load_boundary_candidates(config)
    assert len(boundary) == 1
    candidate = boundary.iloc[0]
    assert candidate["host_sample"] == "43PRr3"
    assert candidate["virus_samples"] == ["44Dr2"]
    assert candidate["analysis_role"] == "boundary_sensitive_provisional"
    assert [row["name"] for row in candidate["boundary_intervals"]] == [
        "catalog_viral_core", "disputed_checkv_extension", "checkv_host_flank",
    ]
    assert sum(row["end"] - row["start"] + 1 for row in candidate["boundary_intervals"]) == 67924
    assert candidate["nested_competitor_ids"] == ["44Dr2::44Dr2_4002037|provirus_1_34828"]


def test_sbatch_and_launcher_shaheen_memory_array_and_dependency_contract():
    val_text = (CODE / "sbatch_stage4_validate_ref.sh").read_text()
    map_text = (CODE / "sbatch_stage4_map_array.sh").read_text()
    agg_text = (CODE / "sbatch_stage4_aggregate.sh").read_text()
    seq_text = (CODE / "sbatch_stage4_sequence_followup.sh").read_text()
    launch_text = (CODE / "submit_stage4_competitive_mapping_pipeline.sh").read_text()

    assert "#SBATCH --mem=350G" in val_text
    assert "#SBATCH --mem=350G" in map_text
    assert "370G" not in val_text
    assert "370G" not in map_text
    assert "#SBATCH --array=0-13\n" in map_text
    assert "%4" not in map_text
    for sample in [
        "41PRr2", "V60PRr2", "55Dr2", "9Sr2", "55Sr3",
        "S49Dr1", "48Sr2", "F48Sr2", "V48Dr1", "S48PRr1",
        "41Sr3", "10Dr2", "43PRr3", "44Dr2",
    ]:
        assert f'"{sample}"' in map_text
    for script_text in [val_text, map_text, agg_text, seq_text]:
        assert 'export TMPDIR="${EXP3_ROOT}/tmp/${SLURM_JOB_ID}"' in script_text
        assert "/tmp" not in script_text.replace("${EXP3_ROOT}/tmp", "")
    assert "sbatch --dependency=afterok:${JOB_REF}" in launch_text
    assert "sbatch --dependency=afterok:${JOB_VAL}" in launch_text
    assert "sbatch --dependency=afterok:${JOB_MAP}" in launch_text
    assert "sbatch --dependency=afterok:${JOB_AGGR}" in launch_text
