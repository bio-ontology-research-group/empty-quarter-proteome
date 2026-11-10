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
    assert 'VAL_DEP="afterok:${JOB_REF}"' in launch_text
    assert 'VAL_DEP="satisfied_completed_${JOB_REF}"' in launch_text
    assert "sbatch --dependency=afterok:${JOB_VAL}" in launch_text
    assert "sbatch --dependency=afterok:${JOB_MAP}" in launch_text
    assert "completed_independent_of_mapping_15894662" in launch_text


def test_preserved_executed_code_hashes_and_orf10_official_sources():
    import stage4_posthoc_sequence_domain_qualification as posthoc_qual

    assert stream_sha256(CODE / "stage4_structural_hallmark_screen.py") == "8fb33e18af4dc6b165acc7a1b932c1ff5ec3cf1f09e1b4d2bdcac1c25a3abbd9"
    assert stream_sha256(CODE / "stage2_viral_identity_and_spacers_v3.py") == "e6c50836d24b926bd7d0bf089f71ab623109fed221691020b344b5b5de68562e"
    assert stream_sha256(CODE / "submit_stage4_competitive_mapping_pipeline.sh") == "3485410bd6a34e3191611f8191560af4eeb732e4ae12afd6b2c5ca51df7467ee"
    for hist_snap in [
        CODE.parent / "version_snapshots/pi_viral_repair_20261109/code/submit_stage4_competitive_mapping_pipeline.sh",
        CODE.parent / "provenance/version_snapshots/pi_viral_repair_20261109/code/submit_stage4_competitive_mapping_pipeline.sh",
    ]:
        if hist_snap.exists():
            assert stream_sha256(hist_snap) == "9bd3887b4dc91163d9c73650e6e573299a757822f06854f40fa8d88412c86ce3"
    for sent_snap in [
        CODE.parent / "version_snapshots/pi_viral_sentinel_repair_20261109/code/submit_stage4_competitive_mapping_pipeline.sh",
        CODE.parent / "provenance/version_snapshots/pi_viral_sentinel_repair_20261109/code/submit_stage4_competitive_mapping_pipeline.sh",
    ]:
        if sent_snap.exists():
            assert stream_sha256(sent_snap) == "67e5ee0e881c98f4e8b1af6554fb7e5435a58eb95cd2a3b77d23c4266e936fa0"

    orf10_refs = posthoc_qual.parse_local_orf10_reference_sources()
    assert orf10_refs["PF09979"]["name"] == "Pa193 gp24 scaffolding protein"
    assert orf10_refs["PF09979"]["short_name"] == "Pa193_gp24"
    assert orf10_refs["PF09979"]["interpro_integrated"] == "IPR016913"
    assert orf10_refs["K09960"]["name"] == "uncharacterized protein"
    assert "COG3566" in orf10_refs["K09960"]["dblinks"]
    assert "Uncharacterized conserved protein, DUF2213 domain" in orf10_refs["COG3566"]["heading"]


def _resolve_evidence_dir() -> Path:
    root = Path(__file__).resolve().parents[1]
    if (root / "results/stage4_unannotated_proteins_pfam_domtbl.txt").exists():
        return root / "results"
    return root / "compact_evidence"


def test_domtblout_coordinate_coverage_ievalue_and_overlap_qualification():
    import stage4_posthoc_sequence_domain_qualification as posthoc_qual

    res_dir = _resolve_evidence_dir()
    domtbl_path = res_dir / "stage4_unannotated_proteins_pfam_domtbl.txt"
    extract_path = res_dir / "stage4_target_genomad_checkv_pfam_ga_extract.json"
    ext = json.loads(extract_path.read_text())

    assert ext["schema_version"] == 2
    prov = ext["extraction_provenance"]
    assert prov["upstream_slurm_jobs"]["sequence_domain_annotation_job_id"] == "15894662"
    assert prov["source_files"]["pfam_db_sha256_from_job_15894662_summary"] == "a78a42d6faf265b6bfca59e8f062d06fae6083ce2c6e335d7b381f20b82b7903"
    assert prov["source_files"]["raw_domtbl_sha256"] == "e5ee9a96f18395a85518f04acbde489a81866454ff29f58002f9485c9f8348f5"

    raw_df = posthoc_qual.parse_domtblout_full_coordinates(domtbl_path, domain_evalue_cutoff=1e-3)
    assert len(raw_df) == 325
    assert raw_df["query_gene"].nunique() == 80
    assert int(raw_df["passes_full_seq_e_1e3"].sum()) == 325

    ie_df = raw_df[raw_df["passes_independent_domain_ie_1e3"]]
    assert len(ie_df) == 153
    assert ie_df["query_gene"].nunique() == 75

    qualified = posthoc_qual.annotate_overlapping_pfam_intervals(raw_df, pfam_ga_map=ext["pfam_ga_map"])
    ga_df = qualified[qualified["passes_pfam_ga"] == True]  # noqa: E712
    assert len(ga_df) == 82
    assert ga_df["query_gene"].nunique() == 50
    assert bool(ga_df["passes_independent_domain_ie_1e3"].all())

    orf12 = qualified[qualified["query_gene"] == "V60PRr2::V60PRr2_715976_12"].set_index("pfam_target_name")
    assert set(orf12.index) == {"Linocin_M18", "DUF6260", "Phage_capsid"}
    assert orf12["overlap_cluster_id"].nunique() == 1
    assert bool(orf12["is_competing_overlapping_hit"].all())
    assert orf12.loc["Linocin_M18", "best_ievalue_model_in_cluster"] == "Linocin_M18(PF04454.19)"
    assert bool(orf12.loc["Linocin_M18", "passes_independent_domain_ie_1e3"]) is True
    assert bool(orf12.loc["Linocin_M18", "passes_pfam_ga"]) is True
    assert (int(orf12.loc["Linocin_M18", "hmm_from"]), int(orf12.loc["Linocin_M18", "hmm_to"]), int(orf12.loc["Linocin_M18", "tlen"])) == (7, 234, 253)
    assert (int(orf12.loc["Linocin_M18", "ali_from"]), int(orf12.loc["Linocin_M18", "ali_to"]), int(orf12.loc["Linocin_M18", "qlen"])) == (22, 268, 305)
    assert bool(orf12.loc["DUF6260", "passes_independent_domain_ie_1e3"]) is True
    assert bool(orf12.loc["DUF6260", "passes_pfam_ga"]) is False
    assert bool(orf12.loc["Phage_capsid", "passes_independent_domain_ie_1e3"]) is False
    note12 = str(orf12.loc["Linocin_M18", "functional_qualification_note"])
    assert "terminase_and_capsid_assembly_related_profile_neighborhood" in note12
    assert "terminase_protease_neighborhood" not in note12
    assert "not_peroxidase_or_cargo_enzyme" in note12

    orf41 = qualified[qualified["query_gene"] == "V60PRr2::V60PRr2_715976_41"]
    assert len(orf41) == 8
    assert int(orf41["passes_independent_domain_ie_1e3"].sum()) == 5
    assert int((orf41["passes_pfam_ga"] == True).sum()) == 2  # noqa: E712


def test_legacy_archaeal_viral_hallmarks_keyword_semantics_and_posthoc_tables():
    import stage4_posthoc_sequence_domain_qualification as posthoc_qual

    res_dir = _resolve_evidence_dir()
    summary = posthoc_qual.run_posthoc_qualification(res_dir)
    assert summary["schema_version"] == 2
    assert summary["executed_sequence_script_sha256"] == "8fb33e18af4dc6b165acc7a1b932c1ff5ec3cf1f09e1b4d2bdcac1c25a3abbd9"
    assert summary["threshold_hierarchy"]["raw_full_sequence_evalue_lte_1e3"]["domain_rows"] == 325
    assert summary["threshold_hierarchy"]["raw_full_sequence_evalue_lte_1e3"]["unique_proteins"] == 80
    assert summary["threshold_hierarchy"]["posthoc_independent_domain_ievalue_lte_1e3"]["domain_rows"] == 153
    assert summary["threshold_hierarchy"]["posthoc_independent_domain_ievalue_lte_1e3"]["unique_proteins"] == 75
    assert summary["threshold_hierarchy"]["posthoc_pfam_gathering_threshold_ga"]["domain_rows"] == 82
    assert summary["threshold_hierarchy"]["posthoc_pfam_gathering_threshold_ga"]["unique_proteins"] == 50
    assert summary["threshold_hierarchy"]["posthoc_pfam_gathering_threshold_ga"]["domain_rows_also_passing_ie_1e3"] == 82
    assert summary["structural_followup_qualification"]["status"] == "not_performed_no_exact_sequence_matched_existing_models"
    assert summary["structural_followup_qualification"]["n_exact_sequence_matched_models"] == 0

    orf10_audit = summary["focal_target_V60PRr2_715976_audit"]["orf8_to_orf12_genomic_neighborhood"]["intervening_orf10"]
    assert orf10_audit["marker"] == "GENOMAD.101057.VV"
    assert orf10_audit["annotation_description"] == "NA"
    assert orf10_audit["official_reference_sources_2026_11_09"]["PF09979"]["name"] == "Pa193 gp24 scaffolding protein"
    assert "NOT a confirmed biochemical prohead protease" in orf10_audit["scientific_qualification"]

    summary_text = json.dumps(summary)
    assert "terminase_protease_neighborhood" not in summary_text
    assert "immediately downstream of ORF8" not in summary_text

    focal_table = pd.read_csv(res_dir / "stage4_focal_V60PRr2_715976_orf_evidence_table.tsv", sep="\t")
    assert len(focal_table) == 51
    assert int(focal_table["was_screened_as_unannotated_candidate"].sum()) == 45
    assert int((focal_table["n_posthoc_ie_1e3_pfam_rows"] > 0).sum()) == 6
    orf10_row = focal_table[focal_table["gene"] == "V60PRr2::V60PRr2_715976_10"].iloc[0]
    assert "Pa193_gp24_capsid_assembly_scaffolding_related_profile" in str(orf10_row["qualified_scientific_interpretation"])
    assert "not_confirmed_protease" in str(orf10_row["qualified_scientific_interpretation"])

    orf8_sem = posthoc_qual.audit_genomad_gene_hallmark_and_keyword_semantics({
        "gene": "V60PRr2::V60PRr2_715976_8",
        "virus_hallmark": "1",
        "taxid": "1246",
        "taxname": "Caudoviricetes",
        "annotation_description": "phage terminase, large subunit, PBSX family",
    })
    assert orf8_sem["is_genomad_virus_hallmark"] is True
    assert orf8_sem["matches_legacy_arch_viral_kw_re"] is True
    assert orf8_sem["matches_generic_morphogenetic_kw"] is True
    assert orf8_sem["matches_archaea_specific_viral_kw"] is False
    assert orf8_sem["taxname_is_broad_caudoviricetes_not_host_domain_specific"] is True
    assert "not_archaea_specific_proof" in orf8_sem["legacy_column_semantics_note"]


def test_host_array_contig_provenance_and_genomad_marker_semantics():
    import stage4_host_array_provenance_audit as host_prov

    summary = host_prov.run_host_array_provenance_audit(CODE.parent)
    assert summary["schema_version"] == 2
    assert summary["status"] == "host_array_contig_provenance_completed_verified"
    assert summary["source_extract"]["extract_sha256"] == (
        "b486422974e76f504f5a98b6ac2a859377b96615362d745c57f3207dd324328d"
    )

    # Regression guard: official geNomad marker classes are CC, CP, CV, PC, PP, PV, VC, VP, VV
    # (never invented .A / .B domain suffixes) and batch_taxonomy.tsv is viral ICTV only.
    with pytest.raises(ValueError, match="chromosome/plasmid/virus specificity classes"):
        host_prov.summarize_genomad_marker_classes([{"marker": "GENOMAD.000001.A"}])
    with pytest.raises(ValueError, match="chromosome/plasmid/virus specificity classes"):
        host_prov.summarize_genomad_marker_classes([{"marker": "GENOMAD.000002.B"}])
    with pytest.raises(ValueError, match="viral ICTV lineages only"):
        host_prov.validate_genomad_viral_ictv_taxonomy({"lineage": "d__Archaea;p__Thermoplasmatota", "n_genes_with_taxonomy": 5})
    with pytest.raises(ValueError, match="viral ICTV lineages only"):
        host_prov.validate_genomad_viral_ictv_taxonomy({"lineage": "Bacteria;Proteobacteria", "n_genes_with_taxonomy": 3})

    summary_json_text = json.dumps(summary)
    assert ".A=" not in summary_json_text
    assert ".B=" not in summary_json_text

    # Focal Thermoplasmatota host 41PRr2__single_SemiBin_236 reconciliation
    rec = summary["focal_thermoplasmatota_host_41PRr2_SemiBin_236_reconciliation"]
    assert rec["focal_host_mag"] == "41PRr2__single_SemiBin_236"
    assert rec["focal_array_contig"] == "41PRr2::41PRr2_1588631"
    assert rec["focal_contig_length_bp"] == 651471
    assert rec["focal_contig_length_bp"] >= 438682
    assert rec["focal_contig_length_rank_in_mag"] == 2
    assert rec["focal_contig_fraction_of_mag_bases"] == pytest.approx(0.193152, abs=1e-6)
    assert rec["focal_contig_raw_seq_sha256"] == "0f98513b8967dffe42f9fd040ab970a13e7e33bd7050dd597994d4dab4c2206f"
    assert rec["focal_contig_canonical_seq_sha256"] == "0b93476eaf5b87f468e7ac8680429f8642d395f256b02e606bd17dfc1763bc1f"
    assert rec["n_physical_crispr_array_loci"] == 1
    assert rec["n_crispr_caller_records"] == 2
    assert rec["array_interval_1based"] == [433346, 438682]
    assert rec["adjacent_type_ie_cas_operon"] == "41PRr2_1588631@8:I-E:424731..432985"
    assert rec["n_raw_extracted_matching_spacers"] == 7
    assert rec["n_strict_unique_spacers_verified"] == 6
    assert rec["n_exact_0mm_unique_spacers_verified"] == 4
    assert rec["n_exploratory_only_2mm_spacers"] == 1
    assert rec["exploratory_only_spacer_ids"] == ["archspc_000625"]
    assert rec["all_spacers_verified_literal_in_contig_fasta"] is True
    assert rec["root_independent_join_verified"] is True

    focal_ctg = summary["hosts"]["41PRr2__single_SemiBin_236"]["array_bearing_contigs"]["41PRr2_1588631"]
    assert focal_ctg["genomad_contig_chromosome_score"] == "0.9521"
    assert focal_ctg["genomad_contig_virus_score"] == "0.0148"
    assert focal_ctg["genomad_contig_total_orfs"] == 600
    assert focal_ctg["genomad_contig_marker_orfs_count"] == 179
    assert focal_ctg["genomad_contig_uscg_count"] == 23
    assert focal_ctg["genomad_marker_class_counts_summary"] == "CP:98;CC:35;CV:20;VC:13;PC:6;VV:6;VP:1"
    assert focal_ctg["crosswalk_row_count"] == 2
    assert focal_ctg["crosswalk_distinct_host_mag_owners"] == "41PRr2__single_SemiBin_236"
    assert focal_ctg["crosswalk_other_source_entries"] == "virus_target:41PRr2::41PRr2_1588631"
    assert focal_ctg["other_entry_in_reviewed_viral_pairs"] is False
    assert focal_ctg["qualified_source_alias_interpretation"] == (
        "source_category_duplication_same_contig_in_raw_search_pool_unique_mag_owner"
    )
    corrob_summary = focal_ctg["corroborative_coding_annotations_summary"]
    assert "41PRr2_1588631_359(397553..398164" in corrob_summary and "K03049" in corrob_summary
    assert "41PRr2_1588631_361(400408..401073" in corrob_summary and "K02967" in corrob_summary
    assert "41PRr2_1588631_408(449433..451175" in corrob_summary and "K07557" in corrob_summary
    assert "41PRr2_1588631_137(132048..132620" in corrob_summary and "K19664" in corrob_summary
    assert "Ribosomal protein S24E" in corrob_summary
    assert "Ribosomal protein L37E" in corrob_summary
    assert "conditional on MAG bin membership" in focal_ctg["scientific_qualification"]

    res_dir = _resolve_evidence_dir()
    ctg_table = pd.read_csv(res_dir / "stage4_host_array_contig_provenance_table.tsv", sep="\t")
    assert len(ctg_table) == 8
    assert set(ctg_table["host_mag"]) == {
        "41PRr2__single_SemiBin_236",
        "41PRr2__single_SemiBin_97",
        "S49Dr1__self_SemiBin_1106",
        "V48Dr1__self_SemiBin_3215",
        "43PRr3__self_SemiBin_469",
        "41Sr3__self_SemiBin_108",
    }
    assert (ctg_table["crosswalk_row_count"] == 2).all()
    assert (ctg_table["other_entry_in_reviewed_viral_pairs"] == False).all()  # noqa: E712
    assert (
        ctg_table["qualified_source_alias_interpretation"]
        == "source_category_duplication_same_contig_in_raw_search_pool_unique_mag_owner"
    ).all()
    for qual_text in ctg_table["scientific_qualification"]:
        assert "conditional on MAG bin membership" in str(qual_text)


def test_host_contig_mapper_coverage_join_by_representative_header_and_fail_closed(tmp_path):
    import stage4_host_array_provenance_audit as host_prov

    cov_tsv = tmp_path / "41PRr2_filtered_coverage.tsv"
    cov_tsv.write_text(
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
        "HOST__41PRr2__single_SemiBin_236__41PRr2_1588631\t1\t651471\t12000\t600000\t92.099265\t2.75\t35.1\t41.2\n"
        "HOST__41PRr2__single_SemiBin_236__41PRr2_other1\t1\t348529\t6000\t300000\t86.076051\t2.50\t35.0\t40.8\n",
        encoding="utf-8",
    )
    parsed = host_prov.parse_samtools_filtered_coverage_tsv(cov_tsv)
    assert set(parsed.keys()) == {
        "HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
        "HOST__41PRr2__single_SemiBin_236__41PRr2_other1",
    }

    host_cw = [
        {
            "source_category": "host_mag",
            "entity_id": "41PRr2__single_SemiBin_236",
            "contig_id": "41PRr2_1588631",
            "representative_header": "HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
            "seq_length": "651471",
        },
        {
            "source_category": "host_mag",
            "entity_id": "41PRr2__single_SemiBin_236",
            "contig_id": "41PRr2_other1",
            "representative_header": "HOST__41PRr2__single_SemiBin_236__41PRr2_other1",
            "seq_length": "348529",
        },
    ]
    eval_row = {
        "unique_mag_contigs": "1.0",
        "unique_mag_bases": "348529.0",
        "unique_filtered_reads_mapped": "6000.0",
        "unique_horizontal_breadth_pct": "30.0",
        "unique_mean_depth": "0.871323",
    }
    res = host_prov.evaluate_host_array_contig_coverage_from_mapper(
        host_mag="41PRr2__single_SemiBin_236",
        array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
        array_contig_length_bp=651471,
        mag_total_contigs=2,
        mag_total_bases=1000000,
        host_crosswalk_rows=host_cw,
        coverage_by_rname=parsed,
        focal_eval_host_row=eval_row,
    )
    assert res["sample_mapping_coverage_status"] == "completed_sample_mapping_integrated_via_representative_header"
    assert res["sample_mapping_contig_rname"] == "HOST__41PRr2__single_SemiBin_236__41PRr2_1588631"
    assert res["sample_mapping_contig_mean_depth"] == pytest.approx(2.75)
    assert res["sample_mapping_remaining_mag_mean_depth"] == pytest.approx(2.50)
    assert res["sample_mapping_total_mag_numreads"] == 18000
    assert res["mapper_raw_conservative_unique_mag_contigs"] == 1
    assert res["mapper_raw_conservative_unique_mag_bases"] == 348529
    assert res["mapper_raw_conservative_unique_mag_reads"] == 6000

    with pytest.raises(ValueError, match="non-integer float"):
        host_prov.evaluate_host_array_contig_coverage_from_mapper(
            host_mag="41PRr2__single_SemiBin_236",
            array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
            array_contig_length_bp=651471,
            mag_total_contigs=2,
            mag_total_bases=1000000,
            host_crosswalk_rows=host_cw,
            coverage_by_rname=parsed,
            focal_eval_host_row={**eval_row, "unique_mag_contigs": "1.5"},
        )

    # Fail-closed: missing != zero when a host representative_header is absent from coverage
    incomplete_cov = {"HOST__41PRr2__single_SemiBin_236__41PRr2_1588631": parsed["HOST__41PRr2__single_SemiBin_236__41PRr2_1588631"]}
    with pytest.raises(ValueError, match="missing != zero"):
        host_prov.evaluate_host_array_contig_coverage_from_mapper(
            host_mag="41PRr2__single_SemiBin_236",
            array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
            array_contig_length_bp=651471,
            mag_total_contigs=2,
            mag_total_bases=1000000,
            host_crosswalk_rows=host_cw,
            coverage_by_rname=incomplete_cov,
        )


def test_validator_accepts_zero_byte_reference_build_done_sentinel_and_rejects_missing_or_empty_data(tmp_path):
    import stage4_validate_competitive_reference as validator

    build_done = tmp_path / "REFERENCE_BUILD_DONE"
    ref_fna = tmp_path / "unified_competitive_reference.fna"
    crosswalk_path = tmp_path / "stage4_reference_sequence_crosswalk.tsv"
    build_receipt_path = tmp_path / "unified_competitive_reference.receipt.json"
    refined_tar = tmp_path / "refined_mags.tar"

    ref_fna.write_text(">HOST__m1__c1\nACGT\n", encoding="utf-8")
    crosswalk_path.write_text("representative_header\tseq_length\nHOST__m1__c1\t4\n", encoding="utf-8")
    build_receipt_path.write_text(json.dumps({"status": "completed"}), encoding="utf-8")
    refined_tar.write_bytes(b"nonempty_tar_placeholder")

    # 1. Missing REFERENCE_BUILD_DONE sentinel is rejected
    with pytest.raises(FileNotFoundError, match="required build sentinel missing"):
        validator.verify_build_inputs_and_sentinel(
            build_done=build_done,
            nonempty_paths=[ref_fna, crosswalk_path, build_receipt_path, refined_tar],
            build_receipt_path=build_receipt_path,
        )

    # 2. Zero-byte touch-created REFERENCE_BUILD_DONE sentinel is accepted when all data/receipt files are nonempty
    build_done.touch()
    assert build_done.stat().st_size == 0
    receipt = validator.verify_build_inputs_and_sentinel(
        build_done=build_done,
        nonempty_paths=[ref_fna, crosswalk_path, build_receipt_path, refined_tar],
        build_receipt_path=build_receipt_path,
    )
    assert receipt["status"] == "completed"

    # 3. Empty data artifact (e.g., 0-byte ref_fna) is still rejected
    empty_fna = tmp_path / "empty.fna"
    empty_fna.touch()
    with pytest.raises(FileNotFoundError, match="required build/source input missing or empty"):
        validator.verify_build_inputs_and_sentinel(
            build_done=build_done,
            nonempty_paths=[empty_fna, crosswalk_path, build_receipt_path, refined_tar],
            build_receipt_path=build_receipt_path,
        )

    # 4. Non-completed build receipt status is still rejected
    bad_receipt = tmp_path / "bad_receipt.json"
    bad_receipt.write_text(json.dumps({"status": "running"}), encoding="utf-8")
    with pytest.raises(ValueError, match="builder receipt is not completed"):
        validator.verify_build_inputs_and_sentinel(
            build_done=build_done,
            nonempty_paths=[ref_fna, crosswalk_path, bad_receipt, refined_tar],
            build_receipt_path=bad_receipt,
        )


def test_boundary_json_load_analysis_pairs_and_validator_consumer_contract():
    from stage4_competitive_common import load_analysis_pairs
    import stage4_validate_competitive_reference as validator

    res_dir = _resolve_evidence_dir()
    pair_path = res_dir / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    boundary_path = CODE.parent / "inputs/stage4_boundary_sensitive_intervals.json"

    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_path, boundary_path)
    assert len(boundary_pairs) == 1
    assert "target_length" in boundary_pairs.columns
    candidate = boundary_pairs.to_dict(orient="records")[0]
    assert int(candidate["target_length"]) == 67924
    assert candidate["target_id"] == "44Dr2::44Dr2_4002037"
    assert candidate["catalog_viral_interval"] == [1, 34828]
    assert candidate["nested_competitor_ids"] == ["44Dr2::44Dr2_4002037|provirus_1_34828"]
    assert candidate["individual_virality_status"] == "boundary_sensitive_provisional_not_primary_host_link"
    assert candidate["analysis_role"] == "boundary_sensitive_provisional"

    target_sha = candidate["target_seq_sha256"]
    full_rep = "VIRUS_ALL__44Dr2::44Dr2_4002037"
    nested_rep = "VIRUS_VOTU__44Dr2::44Dr2_4002037|provirus_1_34828"
    full_seq = ("ACGT" * 16981)[:67924]
    nested_seq = full_seq[:34828]
    focal_sequences = {full_rep: full_seq, nested_rep: nested_seq}
    resolved = {
        target_sha: {
            "reference_contig": full_rep,
            "length_bp": 67924,
            "reference_orientation_to_ledger": "same_as_reference",
        }
    }
    crosswalk = pd.DataFrame([
        {
            "source_category": "all_viruses",
            "entity_id": candidate["target_id"],
            "contig_id": candidate["target_id"],
            "original_header": full_rep,
            "representative_header": full_rep,
            "canonical_sha256": canonical_sha256(full_seq),
            "seq_length": 67924,
            "is_reference_representative": 1,
        },
        {
            "source_category": "votu_reps",
            "entity_id": "44Dr2::44Dr2_4002037|provirus_1_34828",
            "contig_id": "44Dr2::44Dr2_4002037|provirus_1_34828",
            "original_header": nested_rep,
            "representative_header": nested_rep,
            "canonical_sha256": canonical_sha256(nested_seq),
            "seq_length": 34828,
            "is_reference_representative": 1,
        },
    ])

    nested_val = validator.validate_boundary_candidates_and_nested_competitors(
        boundary_pairs=boundary_pairs,
        resolved=resolved,
        focal_sequences=focal_sequences,
        crosswalk=crosswalk,
    )
    assert nested_val == {
        "44Dr2::44Dr2_4002037|provirus_1_34828": {
            "reference_contig": nested_rep,
            "length_bp": 34828,
            "exact_catalog_core_match_forward_or_reverse_complement": True,
        }
    }

    # Verify fail-closed if reference length disagrees with candidate["target_length"]
    bad_resolved = {
        target_sha: {
            "reference_contig": full_rep,
            "length_bp": 67923,
            "reference_orientation_to_ledger": "same_as_reference",
        }
    }
    with pytest.raises(ValueError, match="boundary target length differs from root disposition"):
        validator.validate_boundary_candidates_and_nested_competitors(
            boundary_pairs=boundary_pairs,
            resolved=bad_resolved,
            focal_sequences=focal_sequences,
            crosswalk=crosswalk,
        )


def test_focal_host_taxonomy_and_quality_reconciliation_across_extract_table_summary_and_report():
    res_dir = _resolve_evidence_dir()
    extract = json.loads((res_dir / "stage4_host_array_contig_provenance_extract.json").read_text())
    summary = json.loads((res_dir / "stage4_host_array_provenance_summary.json").read_text())
    table = pd.read_csv(res_dir / "stage4_host_array_contig_provenance_table.tsv", sep="\t")

    focal_mag = "41PRr2__single_SemiBin_236"
    expected_tax = "d__Archaea;p__Thermoplasmatota;c__E2;o__JACPAO01;f__JAHWKR01;g__;s__"

    raw_ref = extract["hosts"][focal_mag]["mag_taxonomy_and_quality"]["refined_representatives_record"]
    raw_gunc = extract["hosts"][focal_mag]["mag_taxonomy_and_quality"]["gunc_per_genome_record"]
    assert raw_ref["gtdb_r232_classification"] == expected_tax
    assert float(raw_ref["Completeness"]) == pytest.approx(93.92)
    assert float(raw_ref["Contamination"]) == pytest.approx(0.54)
    assert str(raw_gunc["pass.GUNC"]) == "True"
    assert float(raw_gunc["clade_separation_score"]) == pytest.approx(0.02)

    sum_host = summary["hosts"][focal_mag]
    assert sum_host["mag_gtdb_r232_taxonomy"] == expected_tax
    assert float(sum_host["mag_checkm2_completeness"]) == pytest.approx(93.92)
    assert float(sum_host["mag_checkm2_contamination"]) == pytest.approx(0.54)
    gunc_sum = sum_host["mag_gunc_summary"]
    assert str(gunc_sum["pass.GUNC"]) == "True"
    assert gunc_sum["taxonomic_level"] == "kingdom"
    assert float(gunc_sum["clade_separation_score"]) == pytest.approx(0.02)
    assert float(gunc_sum["contamination_portion"]) == pytest.approx(0.45)
    assert float(gunc_sum["reference_representation_score"]) == pytest.approx(0.32)

    row = table[table["host_mag"] == focal_mag].iloc[0]
    assert row["mag_gtdb_r232_taxonomy"] == expected_tax
    assert row["mag_order"] == "o__JACPAO01"
    assert row["mag_family"] == "f__JAHWKR01"
    assert str(row["mag_genus"]) in ("g__", "nan", "")
    assert float(row["mag_checkm2_completeness"]) == pytest.approx(93.92)
    assert float(row["mag_checkm2_contamination"]) == pytest.approx(0.54)
    assert float(row["mag_gunc_css"]) == pytest.approx(0.02)

    report_path = CODE.parent.parent / "review-20261109/pi-viral-repair-report.txt"
    if report_path.exists():
        rep_text = report_path.read_text(encoding="utf-8")
        assert expected_tax in rep_text
        assert "93.92%" in rep_text
        assert "0.54%" in rep_text
        assert "CSS=0.02" in rep_text
        for stale_str in ("JACQPQ01", "92.71%", "2.70%", "CSS=0.05"):
            assert stale_str not in rep_text


def test_posthoc_instrain_callable_semantics_breadth_vs_breadth_mincov_and_missing_snv_table(tmp_path):
    from stage4_posthoc_instrain_and_mapping_followup import (
        INSTRAIN_PROFILE_FINISHED_BANNER,
        qualify_instrain_profile_output,
    )

    rep = "VIRUS_ALL__V60PRr2::V60PRr2_715976"
    other_rep = "VIRUS_ALL__OTHER_CONTIG"

    # Case 1: breadth (0.80 -> 8000) != breadth_minCov (0.55 -> 5500), with SNVs.tsv present.
    # Also verifies target-only row counting (ignores other_rep row) and SNV_count vs SNS_count distinction.
    prof1 = tmp_path / "prof1"
    (prof1 / "output").mkdir(parents=True)
    (prof1 / "log").mkdir(parents=True)
    (prof1 / "log/log.log").write_text(f"startup\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")
    (prof1 / "output/prof1_scaffold_info.tsv").write_text(
        f"scaffold\tlength\tbreadth\tbreadth_minCov\tcoverage\tnucl_diversity\tSNS_count\tSNV_count\tdivergent_site_count\n"
        f"{rep}\t10000\t0.80\t0.55\t14.2\t0.0012\t1\t1\t2\n",
        encoding="utf-8",
    )
    (prof1 / "output/prof1_SNVs.tsv").write_text(
        "scaffold\tposition\tref_base\tcon_base\tvar_base\tallele_count\tclass\tcryptic\n"
        f"{rep}\t100\tA\tG\tG\t1\tSNS\tFalse\n"
        f"{rep}\t250\tC\tC\tT\t2\tSNV\tFalse\n"
        f"{other_rep}\t50\tT\tC\tC\t1\tSNS\tFalse\n",
        encoding="utf-8",
    )
    q1 = qualify_instrain_profile_output(prof1, rep, min_cov=5)
    assert q1["read_covered_positions_ge1x"] == 8000
    assert q1["inferred_depth_qualified_positions_ge_mincov"] == 5500
    assert q1["depth_qualified_positions_ge_mincov"] == 5500
    assert q1["read_covered_positions_ge1x"] != q1["inferred_depth_qualified_positions_ge_mincov"]
    assert q1["sns_count_scaffold_info"] == 1
    assert q1["snv_count_scaffold_info"] == 1
    assert q1["divergent_site_count_scaffold_info"] == 2
    assert q1["snv_table_status"] == "present_verified_target_rows"
    assert q1["snv_target_row_count"] == 2
    assert q1["snv_target_non_cryptic_row_count"] == 2
    assert q1["exported_non_cryptic_sns_count"] == 1
    assert q1["exported_non_cryptic_snv_count"] == 1

    # Negative case: exported SNVs.tsv containing cryptic=True row is invalid vendor output and must fail closed
    prof1_bad_cryptic = tmp_path / "prof1_bad_cryptic"
    (prof1_bad_cryptic / "output").mkdir(parents=True)
    (prof1_bad_cryptic / "log").mkdir(parents=True)
    (prof1_bad_cryptic / "log/log.log").write_text(f"startup\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")
    (prof1_bad_cryptic / "output/prof1_scaffold_info.tsv").write_text(
        f"scaffold\tlength\tbreadth\tbreadth_minCov\tcoverage\tnucl_diversity\tSNS_count\tSNV_count\tdivergent_site_count\n"
        f"{rep}\t10000\t0.80\t0.55\t14.2\t0.0012\t1\t1\t2\n",
        encoding="utf-8",
    )
    (prof1_bad_cryptic / "output/prof1_SNVs.tsv").write_text(
        "scaffold\tposition\tref_base\tcon_base\tcryptic\n"
        f"{rep}\t100\tA\tG\tFalse\n"
        f"{rep}\t250\tC\tT\tFalse\n"
        f"{rep}\t300\tG\tA\tTrue\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="contains 1 cryptic=True row"):
        qualify_instrain_profile_output(prof1_bad_cryptic, rep, min_cov=5)

    # Case 2: breadth_minCov absent from scaffold_info.tsv AND SNVs.tsv missing on disk -> explicit unknown (never breadth fallback or silent 0)
    prof2 = tmp_path / "prof2"
    (prof2 / "output").mkdir(parents=True)
    (prof2 / "log").mkdir(parents=True)
    (prof2 / "log/log.log").write_text(f"startup\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")
    (prof2 / "output/prof2_scaffold_info.tsv").write_text(
        f"scaffold\tlength\tbreadth\tSNV_count\tSNS_count\tdivergent_site_count\n{rep}\t10000\t0.80\t0\t0\t0\n",
        encoding="utf-8",
    )
    q2 = qualify_instrain_profile_output(prof2, rep, min_cov=5)
    assert q2["read_covered_positions_ge1x"] == 8000
    assert q2["breadth_mincov_fraction"] == "unknown_breadth_minCov_column_absent"
    assert q2["inferred_depth_qualified_positions_ge_mincov"] == "unknown_breadth_minCov_column_absent"
    assert q2["snv_table_status"] == "unknown_snv_table_missing"
    assert q2["snv_target_row_count"] == "unknown_snv_table_missing"

    # Case 3: inStrain 1.9.1 writes 1-byte newline (or header-only) SNVs.tsv when 0 variants + scaffold_info confirms 0 + terminal log banner
    prof3 = tmp_path / "prof3"
    (prof3 / "output").mkdir(parents=True)
    (prof3 / "log").mkdir(parents=True)
    (prof3 / "log/log.log").write_text(f"startup\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")
    (prof3 / "output/prof3_scaffold_info.tsv").write_text(
        f"scaffold\tlength\tbreadth\tbreadth_minCov\tSNS_count\tSNV_count\tdivergent_site_count\n"
        f"{rep}\t10000\t0.75\t0.60\t0\t0\t0\n",
        encoding="utf-8",
    )
    (prof3 / "output/prof3_SNVs.tsv").write_text("\n", encoding="utf-8")
    q3 = qualify_instrain_profile_output(prof3, rep, min_cov=5)
    assert q3["read_covered_positions_ge1x"] == 7500
    assert q3["inferred_depth_qualified_positions_ge_mincov"] == 6000
    assert q3["snv_table_status"] == "present_empty_dataframe_newline_verified_zero_variants_with_scaffold_info_and_log"
    assert q3["snv_target_row_count"] == 0

    # Case 4: Generic startup log strings ("inStrain profile", "Finished", "scaffold_info") without terminal banner MUST fail closed
    prof4 = tmp_path / "prof4"
    (prof4 / "output").mkdir(parents=True)
    (prof4 / "log").mkdir(parents=True)
    (prof4 / "log/log.log").write_text("inStrain profile ... Finished scaffold_info step\n", encoding="utf-8")
    (prof4 / "output/prof4_scaffold_info.tsv").write_text(
        f"scaffold\tlength\tbreadth\tbreadth_minCov\tSNV_count\n{rep}\t10000\t0.75\t0.60\t0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="lacks terminal completion banner"):
        qualify_instrain_profile_output(prof4, rep, min_cov=5)


def test_load_planned_samples_on_actual_real_inputs_reconciles_coassembly_and_copositive_samples():
    from stage4_competitive_common import load_analysis_pairs
    from stage4_posthoc_instrain_and_mapping_followup import PLANNED_14_SAMPLES, load_planned_samples

    res_dir = _resolve_evidence_dir()
    pair_ledger = res_dir / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    boundary_cfg = CODE.parent / "inputs/stage4_boundary_sensitive_intervals.json"
    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_ledger, boundary_cfg)
    planned = load_planned_samples(CODE.parent, pairs)
    assert planned == list(PLANNED_14_SAMPLES)
    assert len(planned) == 14


def test_posthoc_completion_contract_14_samples_eligibility_stale_dirs_and_compare_integration(tmp_path):
    from stage4_competitive_common import atomic_json, atomic_text, load_analysis_pairs, sha256_sequence, stream_sha256
    from stage4_posthoc_instrain_and_mapping_followup import (
        INSTRAIN_COMPARE_FINISHED_BANNER,
        INSTRAIN_PROFILE_FINISHED_BANNER,
        PLANNED_14_SAMPLES,
        run_crosssample_sensitivity_and_comparisons_if_eligible,
    )

    # Copy actual pair ledger, boundary config, and submission manifest into tmp_path
    res_src = _resolve_evidence_dir()
    (tmp_path / "results/stage4_competitive_reference").mkdir(parents=True)
    (tmp_path / "results/stage4_competitive_mapping").mkdir(parents=True)
    (tmp_path / "inputs").mkdir(parents=True)
    (tmp_path / "provenance").mkdir(parents=True)

    pair_ledger = tmp_path / "results/stage3_tier1_pairs_individual_evidence_audit.tsv"
    boundary_cfg = tmp_path / "inputs/stage4_boundary_sensitive_intervals.json"
    pair_ledger.write_text((res_src / "stage3_tier1_pairs_individual_evidence_audit.tsv").read_text(), encoding="utf-8")
    boundary_cfg.write_text((CODE.parent / "inputs/stage4_boundary_sensitive_intervals.json").read_text(), encoding="utf-8")
    atomic_json(
        tmp_path / "provenance/stage4_competitive_pipeline_submission.json",
        {"samples": list(PLANNED_14_SAMPLES)},
    )

    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_ledger, boundary_cfg)
    orig_focal_sha = str(
        reviewed_pairs.loc[reviewed_pairs["target_id"] == "V60PRr2::V60PRr2_715976", "target_seq_sha256"].iloc[0]
    )
    test_seq = "ACGT" * 2500
    focal_sha = sha256_sequence(test_seq)
    pair_ledger.write_text(
        pair_ledger.read_text(encoding="utf-8").replace(orig_focal_sha, focal_sha),
        encoding="utf-8",
    )
    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_ledger, boundary_cfg)
    focal_rep = "VIRUS_ALL__V60PRr2::V60PRr2_715976"

    # Reference validation receipt fixture
    ref_dir = tmp_path / "results/stage4_competitive_reference"
    val_receipt_path = ref_dir / "stage4_reference_validation_receipt.json"
    val_payload = {
        "status": "completed_verified",
        "schema_version": 3,
        "reference_files": {
            "reference_fna": {"path": str(ref_dir / "unified_competitive_reference.fna"), "sha256": "a" * 64},
            "crosswalk_tsv": {"path": str(ref_dir / "stage4_reference_sequence_crosswalk.tsv"), "sha256": "b" * 64},
        },
        "build_receipt": {"path": str(ref_dir / "unified_competitive_reference.receipt.json"), "sha256": "c" * 64},
    }
    atomic_json(val_receipt_path, val_payload)
    atomic_text(ref_dir / "REFERENCE_VALIDATION_DONE", "ok\n")
    val_receipt_sha = stream_sha256(val_receipt_path)

    def _write_sample_bundle(
        sample: str,
        qualifies_focal: bool = False,
        omit_eligible_profile: bool = False,
        write_valid_sensitivity_receipt: bool = False,
    ) -> Path:
        sdir = tmp_path / "results/stage4_competitive_mapping" / sample
        sdir.mkdir(parents=True, exist_ok=True)
        cov_tsv = sdir / f"{sample}_filtered_coverage.tsv"
        eval_tsv = sdir / f"{sample}_focal_targets_evaluation.tsv"
        raw_bam = sdir / f"{sample}_raw.bam"
        filt_bam = sdir / f"{sample}_mapq20.bam"

        cov_tsv.write_text("#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n", encoding="utf-8")
        raw_bam.write_bytes(f"raw_bam_{sample}".encode("utf-8"))
        filt_bam.write_bytes(f"filt_bam_{sample}".encode("utf-8"))

        eval_rows = []
        seen_shas = set()
        for p in pairs.to_dict(orient="records"):
            tsha = p["target_seq_sha256"]
            if tsha in seen_shas:
                continue
            seen_shas.add(tsha)
            is_qual = bool(qualifies_focal and tsha == focal_sha and p["analysis_role"] == "primary")
            rep_name = focal_rep if tsha == focal_sha else f"VIRUS_ALL__{p['target_id']}"
            eval_rows.append({
                "sample": sample,
                "entity_id": p["target_id"],
                "entity_type": "virus_target",
                "entity_sha256": tsha,
                "reference_contig": rep_name,
                "analysis_role": p["analysis_role"],
                "filtered_reads_mapped": 5000 if is_qual else 0,
                "horizontal_breadth_pct": 85.0 if is_qual else 0.0,
                "mean_depth": 22.5 if is_qual else 0.0,
                "uniquely_attributable": True if p["analysis_role"] == "primary" else False,
                "ambiguity_status": "unique_viral_sequence_aliases_collapsed",
                "qualifies_10x_50pct_instrain": 1 if is_qual else 0,
            })
        for host_mag in sorted(set(pairs["host_mag"])):
            eval_rows.append({
                "sample": sample,
                "entity_id": host_mag,
                "entity_type": "host_mag",
                "entity_sha256": "",
                "reference_contig": "multiple",
                "analysis_role": "host_context",
                "filtered_reads_mapped": 100,
                "horizontal_breadth_pct": 10.0,
                "mean_depth": 1.0,
                "uniquely_attributable": False,
                "ambiguity_status": "ownership_derived_from_full_crosswalk",
                "qualifies_10x_50pct_instrain": 0,
            })
        pd.DataFrame(eval_rows).to_csv(eval_tsv, sep="\t", index=False)

        instrain_results = []
        commands = []
        if qualifies_focal and not omit_eligible_profile:
            pdir = sdir / f"inStrain_source_primary_{focal_sha[:16]}"
            (pdir / "output").mkdir(parents=True, exist_ok=True)
            (pdir / "log").mkdir(parents=True, exist_ok=True)
            (pdir / "log/log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")
            (pdir / "output/prof_scaffold_info.tsv").write_text(
                f"scaffold\tlength\tbreadth\tbreadth_minCov\tcoverage\tnucl_diversity\tSNS_count\tSNV_count\tdivergent_site_count\n"
                f"{focal_rep}\t10000\t0.85\t0.70\t22.5\t0.0005\t0\t1\t1\n",
                encoding="utf-8",
            )
            (pdir / "output/prof_SNVs.tsv").write_text(
                f"scaffold\tposition\tref_base\tcon_base\tcryptic\n{focal_rep}\t500\tA\tG\tFalse\n",
                encoding="utf-8",
            )
            out_hashes = {
                str(p.relative_to(pdir)): stream_sha256(p)
                for p in sorted(pdir.rglob("*")) if p.is_file()
            }
            instrain_results.append({
                "target_raw_sha256": focal_sha,
                "reference_contig": focal_rep,
                "profile_dir": str(pdir),
                "status": "profiled_successfully",
                "tool_version": "inStrain version 1.9.1",
                "min_mapq_strict_greater_than": 19,
                "effective_minimum_integer_mapq": 20,
                "min_read_ani_strict_greater_than": 0.95,
                "analysis_mode": "source_primary",
                "callable_positions": int(round(10000 * 0.85)),
                "output_sha256": out_hashes,
            })
            commands.append({
                "description": f"inStrain source_primary profile {focal_sha}",
                "command": ["inStrain", "profile", str(filt_bam), f"focal_{focal_sha[:16]}.fna", "-o", str(pdir), "--min_mapq", "19", "--min_read_ani", "0.95"],
                "returncode": 0,
            })

            # Write focal FASTA whose sequence SHA matches focal_sha in our test fixture
            #Wait: _verify_focal_fasta_against_target_sha checks sha256_sequence(seq) == target_sha!
            #So for our test fixture, let's use a real sequence whose sha256 we know, or patch _verify_focal_fasta_against_target_sha / use focal_fasta with matching sha!

        receipt_path = sdir / f"{sample}_mapping_receipt.json"
        receipt_payload = {
            "status": "completed_verified",
            "schema_version": 3,
            "sample": sample,
            "reference_binding": {
                "validation_receipt_sha256": val_receipt_sha,
                "reference_fna_sha256": "a" * 64,
                "crosswalk_sha256": "b" * 64,
                "build_receipt_sha256": "c" * 64,
            },
            "commands": commands,
            "instrain": {
                "performed": bool(instrain_results),
                "results": instrain_results,
            },
            "outputs": {
                "coverage_tsv": str(cov_tsv),
                "coverage_tsv_sha256": stream_sha256(cov_tsv),
                "focal_evaluation_tsv": str(eval_tsv),
                "focal_evaluation_tsv_sha256": stream_sha256(eval_tsv),
                "raw_bam": str(raw_bam),
                "raw_bam_sha256": stream_sha256(raw_bam),
                "filtered_bam": str(filt_bam),
                "filtered_bam_sha256": stream_sha256(filt_bam),
            },
        }
        atomic_json(receipt_path, receipt_payload)
        atomic_text(sdir / "MAPPING_DONE", json.dumps({"receipt_sha256": stream_sha256(receipt_path)}) + "\n")
        return receipt_path

    # 1. Incomplete 14-sample receipt set (only 13 of 14 samples exist) -> fails closed when require_complete=True,
    #    and returns pending_incomplete_sample_mapping (never completed_verified or non_estimable) when require_complete=False
    for s in PLANNED_14_SAMPLES[:-1]:
        _write_sample_bundle(s, qualifies_focal=False)

    with pytest.raises(RuntimeError, match="incomplete 14-sample competitive mapping receipt set"):
        run_crosssample_sensitivity_and_comparisons_if_eligible(tmp_path, require_complete=True)

    pending_res = run_crosssample_sensitivity_and_comparisons_if_eligible(tmp_path, require_complete=False)
    assert pending_res["status"] == "pending_incomplete_sample_mapping"
    assert pending_res["missing_planned_samples"] == [PLANNED_14_SAMPLES[-1]]
    assert pending_res["population_comparisons"] == []

    # 2. Add the 14th sample, but with an eligible target whose source_primary profile is MISSING -> fails closed!
    _write_sample_bundle(PLANNED_14_SAMPLES[-1], qualifies_focal=True, omit_eligible_profile=True)
    with pytest.raises(ValueError, match="lacks an inStrain profile in receipt"):
        run_crosssample_sensitivity_and_comparisons_if_eligible(tmp_path, require_complete=True)

    # 3. Corrupt MAPPING_DONE hash on a sample -> fails closed!
    _write_sample_bundle(PLANNED_14_SAMPLES[-1], qualifies_focal=False)
    atomic_text(
        tmp_path / "results/stage4_competitive_mapping" / PLANNED_14_SAMPLES[-1] / "MAPPING_DONE",
        json.dumps({"receipt_sha256": "0" * 64}) + "\n",
    )
    with pytest.raises(ValueError, match="MAPPING_DONE does not match current receipt"):
        run_crosssample_sensitivity_and_comparisons_if_eligible(tmp_path, require_complete=True)

    # 4. All 14 samples valid, 2 samples qualify for source_primary, but runtime tool is missing
    #    and a stale partial sensitivity directory exists without receipt -> preserves stale dir and fails closed!
    s0_dir = tmp_path / "results/stage4_competitive_mapping/41PRr2"
    for s in PLANNED_14_SAMPLES:
        _write_sample_bundle(s, qualifies_focal=(s in {"V60PRr2", "41PRr2"}))
    for s in ("V60PRr2", "41PRr2"):
        sdir = tmp_path / "results/stage4_competitive_mapping" / s
        atomic_text(sdir / f"focal_{focal_sha[:16]}.fna", f">{focal_rep}\n{test_seq}\n")

    stale_sens = s0_dir / f"inStrain_crosssample_sensitivity_{focal_sha[:16]}"
    (stale_sens / "output").mkdir(parents=True)
    (stale_sens / "output/partial.txt").write_text("stale partial", encoding="utf-8")

    with pytest.raises(RuntimeError, match="eligible crosssample_sensitivity profile required"):
        run_crosssample_sensitivity_and_comparisons_if_eligible(
            tmp_path,
            require_complete=True,
            instrain_sif=tmp_path / "nonexistent_instrain.sif",
            singularity_bin=None,
        )
    preserved = list(s0_dir.glob(f"inStrain_crosssample_sensitivity_{focal_sha[:16]}.stale_preserved_*"))
    assert len(preserved) == 1
    assert (preserved[0] / "output/partial.txt").exists()

    # 5. Now write valid receipted crosssample_sensitivity profiles for V60PRr2 and 41PRr2,
    #    and a stale partial compare directory without receipt -> preserves stale compare dir and fails closed when tool missing!
    for s in ("V60PRr2", "41PRr2"):
        sdir = tmp_path / "results/stage4_competitive_mapping" / s
        sens_dir = sdir / f"inStrain_crosssample_sensitivity_{focal_sha[:16]}"
        (sens_dir / "output").mkdir(parents=True, exist_ok=True)
        (sens_dir / "log").mkdir(parents=True, exist_ok=True)
        (sens_dir / "log/log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")
        (sens_dir / "output/sens_scaffold_info.tsv").write_text(
            f"scaffold\tlength\tbreadth\tbreadth_minCov\tcoverage\tnucl_diversity\tSNS_count\tSNV_count\tdivergent_site_count\n"
            f"{focal_rep}\t10000\t0.88\t0.74\t24.0\t0.0008\t0\t2\t2\n",
            encoding="utf-8",
        )
        (sens_dir / "output/sens_SNVs.tsv").write_text(
            f"scaffold\tposition\tref_base\tcon_base\tcryptic\n"
            f"{focal_rep}\t500\tA\tG\tFalse\n"
            f"{focal_rep}\t600\tC\tT\tFalse\n",
            encoding="utf-8",
        )
        out_hashes = {
            str(p.relative_to(sens_dir)): stream_sha256(p)
            for p in sorted(sens_dir.rglob("*")) if p.is_file()
        }
        rec_path = sens_dir / "posthoc_instrain_profile_receipt.json"
        atomic_json(
            rec_path,
            {
                "status": "completed_verified",
                "schema_version": 1,
                "sample": s,
                "target_raw_sha256": focal_sha,
                "reference_contig": focal_rep,
                "analysis_mode": "crosssample_sensitivity",
                "tool_version": "inStrain version 1.9.1",
                "parameters": {"min_mapq": 19, "min_read_ani": 0.92, "min_cov": 5},
                "inputs": {
                    "filtered_bam_sha256": stream_sha256(sdir / f"{s}_mapq20.bam"),
                    "focal_fasta_sha256": stream_sha256(sdir / f"focal_{focal_sha[:16]}.fna"),
                },
                "output_sha256": out_hashes,
            },
        )
        atomic_text(sens_dir / "POSTHOC_INSTRAIN_PROFILE_DONE", json.dumps({"receipt_sha256": stream_sha256(rec_path)}) + "\n")

    stale_comp = tmp_path / "results/stage4_competitive_mapping" / f"inStrain_compare_source_primary_{focal_sha[:16]}"
    (stale_comp / "output").mkdir(parents=True)
    (stale_comp / "output/incomplete.tsv").write_text("partial", encoding="utf-8")

    with pytest.raises(RuntimeError, match="eligible inStrain compare required"):
        run_crosssample_sensitivity_and_comparisons_if_eligible(
            tmp_path,
            require_complete=True,
            instrain_sif=tmp_path / "nonexistent_instrain.sif",
            singularity_bin=None,
        )
    preserved_comp = list((tmp_path / "results/stage4_competitive_mapping").glob(f"inStrain_compare_source_primary_{focal_sha[:16]}.stale_preserved_*"))
    assert len(preserved_comp) == 1

    # 6. Now write valid receipted inStrain compare outputs for both source_primary and crosssample_sensitivity
    for mode, scaff_rel in [
        ("source_primary", "output/prof_scaffold_info.tsv"),
        ("crosssample_sensitivity", "output/sens_scaffold_info.tsv"),
    ]:
        comp_dir = tmp_path / "results/stage4_competitive_mapping" / f"inStrain_compare_{mode}_{focal_sha[:16]}"
        (comp_dir / "output").mkdir(parents=True, exist_ok=True)
        (comp_dir / "log").mkdir(parents=True, exist_ok=True)
        (comp_dir / "log/log.log").write_text(f"compare\n{INSTRAIN_COMPARE_FINISHED_BANNER}\n", encoding="utf-8")
        (comp_dir / "output/comp_comparisonsTable.tsv").write_text(
            "scaffold\tname1\tname2\tcoverage_overlap\tcompared_bases_count\tpercent_genome_compared\tlength\tconsensus_SNPs\tpopulation_SNPs\tpopANI\tconANI\n"
            f"{focal_rep}\t41PRr2_mapq20.bam\tV60PRr2_mapq20.bam\t0.82\t6800\t0.68\t10000\t1\t0\t1.0\t0.9998529411764706\n",
            encoding="utf-8",
        )
        out_hashes = {
            str(p.relative_to(comp_dir)): stream_sha256(p)
            for p in sorted(comp_dir.rglob("*")) if p.is_file()
        }
        scaff_hashes = {
            s: stream_sha256(
                tmp_path / "results/stage4_competitive_mapping" / s / f"inStrain_{mode}_{focal_sha[:16]}" / scaff_rel
            )
            for s in ("41PRr2", "V60PRr2")
        }
        crec_path = comp_dir / "posthoc_instrain_compare_receipt.json"
        atomic_json(
            crec_path,
            {
                "status": "completed_verified",
                "schema_version": 1,
                "target_seq_sha256": focal_sha,
                "reference_contig": focal_rep,
                "analysis_mode": mode,
                "qualifying_samples": ["41PRr2", "V60PRr2"],
                "input_profile_scaffold_info_sha256": scaff_hashes,
                "output_sha256": out_hashes,
            },
        )
        atomic_text(comp_dir / "POSTHOC_INSTRAIN_COMPARE_DONE", json.dumps({"receipt_sha256": stream_sha256(crec_path)}) + "\n")

    full_res = run_crosssample_sensitivity_and_comparisons_if_eligible(
        tmp_path,
        require_complete=True,
        instrain_sif=tmp_path / "nonexistent_instrain.sif",
        singularity_bin=None,
    )
    assert full_res["status"] == "completed_verified"
    assert full_res["n_sample_receipts_verified"] == 14
    assert full_res["n_source_primary_profiles_qualified"] == 2
    assert full_res["n_crosssample_sensitivity_profiles_qualified"] == 2
    comp_by_key = {(c["target_seq_sha256"], c["analysis_mode"]): c for c in full_res["population_comparisons"]}
    focal_comp_pri = comp_by_key[(focal_sha, "source_primary")]
    assert focal_comp_pri["comparison_status"] == "completed_verified_instrain_compare"
    assert focal_comp_pri["min_compared_bases_count"] == 6800
    assert focal_comp_pri["comparison_rows"][0]["popANI"] == pytest.approx(1.0)
    assert len(focal_comp_pri["measured_14_sample_evaluation_rows"]) == 14


def test_instrain_191_runtime_contract_probe_and_comparison_parser_strict_validation(tmp_path):
    from stage4_competitive_common import stream_sha256
    from stage4_posthoc_instrain_and_mapping_followup import (
        INSTRAIN_COMPARE_FINISHED_BANNER,
        INSTRAIN_PROFILE_FINISHED_BANNER,
        parse_instrain_comparisons_table,
        qualify_instrain_profile_output,
    )

    res_dir = _resolve_evidence_dir()
    probe = json.loads((res_dir / "stage4_instrain_191_runtime_contract_probe.json").read_text(encoding="utf-8"))
    assert probe["status"] == "completed_verified_instrain_1.9.1_runtime_contract_probe"

    # Reconstruct real Shaheen inStrain 1.9.1 outputs from the probe JSON and verify exact SHA-256 hashes
    prof_a = tmp_path / "profA"
    prof_b = tmp_path / "profB"
    comp_ab = tmp_path / "compAB"
    for d in (prof_a, prof_b, comp_ab):
        (d / "output").mkdir(parents=True, exist_ok=True)
        (d / "raw_data").mkdir(parents=True, exist_ok=True)
        (d / "log").mkdir(parents=True, exist_ok=True)

    (prof_a / "output/profA_scaffold_info.tsv").write_text(probe["profA_zero_variants"]["scaffold_info_tsv"], encoding="utf-8")
    (prof_a / "output/profA_SNVs.tsv").write_text("\n", encoding="utf-8")
    (prof_a / "raw_data/attributes.tsv").write_text(
        f"Attribute\tvalue\ttype\tdescription\n{probe['profA_zero_variants']['bam_loc_attribute']}\n",
        encoding="utf-8",
    )
    (prof_a / "log/log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")

    (prof_b / "output/profB_scaffold_info.tsv").write_text(probe["profB_one_sns_one_snv"]["scaffold_info_tsv"], encoding="utf-8")
    (prof_b / "output/profB_SNVs.tsv").write_text(probe["profB_one_sns_one_snv"]["snvs_tsv"], encoding="utf-8")
    (prof_b / "raw_data/attributes.tsv").write_text(
        f"Attribute\tvalue\ttype\tdescription\n{probe['profB_one_sns_one_snv']['bam_loc_attribute']}\n",
        encoding="utf-8",
    )
    (prof_b / "log/log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")

    (comp_ab / "output/compAB_comparisonsTable.tsv").write_text(
        probe["compAB_positive_and_zero_overlap"]["comparisons_table_tsv"],
        encoding="utf-8",
    )
    (comp_ab / "log/log.log").write_text(f"start\n{INSTRAIN_COMPARE_FINISHED_BANNER}\n", encoding="utf-8")

    assert stream_sha256(prof_a / "output/profA_scaffold_info.tsv") == probe["profA_zero_variants"]["scaffold_info_sha256"]
    assert stream_sha256(prof_a / "output/profA_SNVs.tsv") == probe["profA_zero_variants"]["snvs_tsv_sha256"]
    assert stream_sha256(prof_b / "output/profB_scaffold_info.tsv") == probe["profB_one_sns_one_snv"]["scaffold_info_sha256"]
    assert stream_sha256(prof_b / "output/profB_SNVs.tsv") == probe["profB_one_sns_one_snv"]["snvs_tsv_sha256"]
    assert stream_sha256(comp_ab / "output/compAB_comparisonsTable.tsv") == probe["compAB_positive_and_zero_overlap"]["comparisons_table_sha256"]

    # Verify profA (0 variants, 1-byte newline SNVs.tsv) and profB (1 SNS + 1 SNV)
    qa = qualify_instrain_profile_output(prof_a, "TINY_VIRUS_REF", min_cov=5)
    assert qa["snv_count_scaffold_info"] == 0
    assert qa["sns_count_scaffold_info"] == 0
    assert qa["divergent_site_count_scaffold_info"] == 0
    assert qa["snv_table_status"] == "present_empty_dataframe_newline_verified_zero_variants_with_scaffold_info_and_log"
    assert qa["snv_target_row_count"] == 0

    qb = qualify_instrain_profile_output(prof_b, "TINY_VIRUS_REF", min_cov=5)
    assert qb["snv_count_scaffold_info"] == 1
    assert qb["sns_count_scaffold_info"] == 1
    assert qb["divergent_site_count_scaffold_info"] == 2
    assert qb["snv_table_status"] == "present_verified_target_rows"
    assert qb["snv_target_row_count"] == 2
    assert qb["snv_target_non_cryptic_row_count"] == 2

    qual_profiles = [
        {"sample": "sampleA", "profile_dir": str(prof_a), "instrain_length_bp": 1000},
        {"sample": "sampleB", "profile_dir": str(prof_b), "instrain_length_bp": 1000},
    ]

    # Positive overlap (TINY_VIRUS_REF: 790 compared bases, 2 consensus SNPs, 1 population SNP)
    comp_pos = parse_instrain_comparisons_table(comp_ab, "TINY_VIRUS_REF", qual_profiles)
    assert comp_pos["n_pairwise_comparisons"] == 1
    assert comp_pos["n_positive_overlap_pairs"] == 1
    assert comp_pos["n_zero_overlap_non_estimable_pairs"] == 0
    row_pos = comp_pos["comparison_rows"][0]
    assert row_pos["compared_bases_count"] == 790
    assert row_pos["consensus_SNPs"] == 2
    assert row_pos["population_SNPs"] == 1
    assert row_pos["conANI"] == pytest.approx((790 - 2) / 790)
    assert row_pos["popANI"] == pytest.approx((790 - 1) / 790)
    assert row_pos["pairwise_estimability_status"] == "estimable_positive_overlap_ge_50pct_genome_compared"

    # Zero overlap (TINY_DISJOINT_REF: 0 compared bases -> NA_zero_coverage_overlap, never numeric 0)
    comp_zero = parse_instrain_comparisons_table(comp_ab, "TINY_DISJOINT_REF", qual_profiles)
    assert comp_zero["n_pairwise_comparisons"] == 1
    assert comp_zero["n_positive_overlap_pairs"] == 0
    assert comp_zero["n_zero_overlap_non_estimable_pairs"] == 1
    row_zero = comp_zero["comparison_rows"][0]
    assert row_zero["compared_bases_count"] == 0
    assert row_zero["consensus_SNPs"] == "NA_zero_coverage_overlap"
    assert row_zero["population_SNPs"] == "NA_zero_coverage_overlap"
    assert row_zero["conANI"] == "NA_zero_coverage_overlap"
    assert row_zero["popANI"] == "NA_zero_coverage_overlap"
    assert row_zero["pairwise_estimability_status"] == "non_estimable_zero_coverage_overlap"

    # Empty SNVs.tsv with missing SNS_count or divergent_site_count must report explicit unknown (never 0)
    prof_missing_sns = tmp_path / "prof_missing_sns"
    (prof_missing_sns / "output").mkdir(parents=True)
    (prof_missing_sns / "log").mkdir(parents=True)
    (prof_missing_sns / "log/log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")
    (prof_missing_sns / "output/prof_scaffold_info.tsv").write_text(
        "scaffold\tlength\tbreadth\tbreadth_minCov\tSNV_count\nTINY_VIRUS_REF\t1000\t0.79\t0.79\t0\n",
        encoding="utf-8",
    )
    (prof_missing_sns / "output/prof_SNVs.tsv").write_text("\n", encoding="utf-8")
    q_missing_sns = qualify_instrain_profile_output(prof_missing_sns, "TINY_VIRUS_REF", min_cov=5)
    assert q_missing_sns["snv_table_status"] == "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
    assert q_missing_sns["snv_target_row_count"] == "unknown_scaffold_info_missing_sns_or_divergent_site_columns"

    # Empty SNVs.tsv with SNV_count=0 but SNS_count=1 / divergent_site_count=1 must fail closed
    prof_contradict_sns = tmp_path / "prof_contradict_sns"
    (prof_contradict_sns / "output").mkdir(parents=True)
    (prof_contradict_sns / "log").mkdir(parents=True)
    (prof_contradict_sns / "log/log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")
    (prof_contradict_sns / "output/prof_scaffold_info.tsv").write_text(
        "scaffold\tlength\tbreadth\tbreadth_minCov\tSNS_count\tSNV_count\tdivergent_site_count\n"
        "TINY_VIRUS_REF\t1000\t0.79\t0.79\t1\t0\t1\n",
        encoding="utf-8",
    )
    (prof_contradict_sns / "output/prof_SNVs.tsv").write_text("\n", encoding="utf-8")
    with pytest.raises(ValueError, match="empty SNVs.tsv .* contradicts scaffold_info counts"):
        qualify_instrain_profile_output(prof_contradict_sns, "TINY_VIRUS_REF", min_cov=5)

    # Self-pair or wrong pair in comparisonsTable.tsv must fail closed
    bad_comp = tmp_path / "bad_comp"
    (bad_comp / "output").mkdir(parents=True)
    (bad_comp / "log").mkdir(parents=True)
    (bad_comp / "log/log.log").write_text(f"start\n{INSTRAIN_COMPARE_FINISHED_BANNER}\n", encoding="utf-8")
    (bad_comp / "output/bad_comparisonsTable.tsv").write_text(
        "scaffold\tname1\tname2\tcoverage_overlap\tcompared_bases_count\tpercent_genome_compared\tlength\tconsensus_SNPs\tpopulation_SNPs\tpopANI\tconANI\n"
        "TINY_VIRUS_REF\tsampleA_mapq20.bam\tsampleA_mapq20.bam\t1.0\t790\t0.79\t1000\t0\t0\t1.0\t1.0\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="self-comparison row"):
        parse_instrain_comparisons_table(bad_comp, "TINY_VIRUS_REF", qual_profiles)


def test_instrain_191_namespaced_header_overlay_and_v60_recovery_contract():
    from stage4_posthoc_instrain_and_mapping_followup import (
        INSTRAIN_191_ORIGINAL_SNVPROFILE_SHA256,
        INSTRAIN_191_PATCH_DIFF_SHA256,
        INSTRAIN_191_PATCHED_SNVPROFILE_SHA256,
        verify_and_resolve_instrain_namespaced_header_overlay,
    )
    import stage4_recover_v60prr2_competitive_mapping as v60_rec

    overlay = verify_and_resolve_instrain_namespaced_header_overlay(root=CODE.parent)
    assert overlay["tool_version"] == "1.9.1"
    assert overlay["original_snvprofile_py_sha256"] == INSTRAIN_191_ORIGINAL_SNVPROFILE_SHA256
    assert overlay["patched_snvprofile_py_sha256"] == INSTRAIN_191_PATCHED_SNVPROFILE_SHA256
    assert overlay["patch_diff_sha256"] == INSTRAIN_191_PATCH_DIFF_SHA256
    assert overlay["singularity_bind_arg"].endswith(":/usr/local/lib/python3.8/site-packages/inStrain/SNVprofile.py:ro")

    # Verify the exact 3-line patch against namespaced :: headers
    ns_scaff = "VIRUS_TARGET__V60PRr2::V60PRr2_715976"
    covt_key = f"{ns_scaff}::0"
    with pytest.raises(ValueError, match="too many values to unpack"):
        _s, _m = covt_key.split("::")
    scaff_out, mm_out = covt_key.rsplit("::", 1)
    assert scaff_out == ns_scaff
    assert int(mm_out) == 0

    pair_key = f"{ns_scaff}::sampleA_mapq20.bam-sampleB_mapq20.bam::2"
    with pytest.raises(ValueError, match="too many values to unpack"):
        _s2, _p2, _m2 = pair_key.split("::")
    scaff2, pair2, mm2 = pair_key.rsplit("::", 2)
    assert scaff2 == ns_scaff
    assert pair2 == "sampleA_mapq20.bam-sampleB_mapq20.bam"
    assert int(mm2) == 2

    # Verify stage4_recover_v60prr2_competitive_mapping source contracts from EXP3-V60-RECOVERY-ROOT-REVIEW.txt
    rec_src = Path(v60_rec.__file__).read_text(encoding="utf-8")
    assert "present_empty_dataframe_newline_verified_zero_variants_with_scaffold_info_and_log" in rec_src
    assert "verified_empty_no_header_written_by_instrain_zero_variants" not in rec_src
    assert "2,866,184" not in rec_src
    assert "1_recovered_from_15905725_1" not in rec_src
    assert "shutil.rmtree(ns_probe_dir)" not in rec_src


def test_completed_14_sample_competitive_mapping_posthoc_instrain_and_host_coverage_outputs():
    res_dir = _resolve_evidence_dir()

    # 1. Completed 14-sample competitive mapping aggregation (15905726)
    aggr_sum_path = res_dir / "stage4_competitive_mapping_all_pairs_summary.json"
    aggr_done_path = res_dir / "STAGE4_COMPETITIVE_MAPPING_DONE"
    aggr_sum = json.loads(aggr_sum_path.read_text(encoding="utf-8"))
    aggr_done = json.loads(aggr_done_path.read_text(encoding="utf-8"))
    assert aggr_done["summary_sha256"] == stream_sha256(aggr_sum_path)
    assert stream_sha256(aggr_sum_path) == "558eab32f5dbc43279f77eae2eee868acc3f904a27e05c7aa2c799cd714222a6"
    assert aggr_sum["status"] == "completed_verified"
    assert aggr_sum["n_samples_mapped"] == 14
    assert aggr_sum["n_pair_sample_combinations"] == 98
    assert aggr_sum["n_boundary_interval_sample_rows"] == 42
    assert aggr_sum["n_qualifying_unique_viral_sha_sample_combinations"] == 1
    assert aggr_sum["n_instrain_profiles_completed"] == 1
    assert aggr_sum["inputs"]["mapping_receipts"]["V60PRr2"] == (
        "2c00a7c260d31fd03682e04216dfbedc328d862ca4fc574c3c65d27ab62ef498"
    )

    pairs_df = pd.read_csv(res_dir / "stage4_competitive_mapping_all_pairs_summary.tsv", sep="\t")
    eval_df = pd.read_csv(res_dir / "stage4_competitive_mapping_all_contigs_evaluation.tsv", sep="\t")
    bound_df = pd.read_csv(res_dir / "stage4_competitive_mapping_all_boundary_intervals.tsv", sep="\t")
    assert (len(pairs_df), len(eval_df), len(bound_df)) == (98, 224, 42)

    # 2. Completed V60PRr2 recovery provenance (15910797 preserving 15905725_1 and 15910426)
    v60_rec_path = res_dir / "stage4_competitive_mapping/V60PRr2/V60PRr2_mapping_receipt.json"
    v60_prov_path = res_dir / "stage4_competitive_mapping/V60PRr2/V60PRr2_recovery_provenance.json"
    assert stream_sha256(v60_rec_path) == "2c00a7c260d31fd03682e04216dfbedc328d862ca4fc574c3c65d27ab62ef498"
    v60_prov = json.loads(v60_prov_path.read_text(encoding="utf-8"))
    assert v60_prov["status"] == "completed_verified_v60prr2_instrain_recovery"
    assert v60_prov["recovery_execution"]["slurm_job_id"] == "15910797"
    assert v60_prov["recovery_execution"]["remapped_fastqs"] is False
    assert v60_prov["recovery_execution"]["rebuilt_reference"] is False
    prior_attempts = v60_prov["recovery_execution"]["preserved_failed_attempts"]
    assert [a["source_slurm_job_id"] for a in prior_attempts] == ["15905725_1", "15910426"]

    # 3. Completed posthoc inStrain primary (>0.95) and relaxed sensitivity (>0.92) profiles (15910957)
    post_sum_path = res_dir / "stage4_instrain_posthoc_callable_and_sensitivity_summary.json"
    post_tbl_path = res_dir / "stage4_instrain_posthoc_profiles_table.tsv"
    post_sum = json.loads(post_sum_path.read_text(encoding="utf-8"))
    post_tbl = pd.read_csv(post_tbl_path, sep="\t")
    assert post_sum["status"] == "completed_verified"
    assert post_sum["n_sample_receipts_verified"] == 14
    assert post_sum["n_source_primary_profiles_qualified"] == 1
    assert post_sum["n_crosssample_sensitivity_profiles_qualified"] == 1
    assert len(post_tbl) == 2

    src_prof = post_sum["source_primary_profiles"][0]
    sens_prof = post_sum["crosssample_sensitivity_profiles"][0]
    assert src_prof["inferred_depth_qualified_positions_ge_mincov"] == 38585
    assert src_prof["snv_count_scaffold_info_unfiltered_cumulative"] == 24
    assert src_prof["sns_count_scaffold_info_unfiltered_cumulative"] == 0
    assert src_prof["exported_non_cryptic_snv_count"] == 22
    assert src_prof["exported_non_cryptic_sns_count"] == 0
    assert src_prof["cryptic_variant_site_count_excluded_from_exported_snv_table"] == 2

    assert sens_prof["inferred_depth_qualified_positions_ge_mincov"] == 38589
    assert sens_prof["snv_count_scaffold_info_unfiltered_cumulative"] == 24
    assert sens_prof["sns_count_scaffold_info_unfiltered_cumulative"] == 1
    assert sens_prof["divergent_site_count_scaffold_info_unfiltered_cumulative"] == 25
    assert sens_prof["exported_non_cryptic_snv_count"] == 21
    assert sens_prof["exported_non_cryptic_sns_count"] == 1
    assert sens_prof["cryptic_variant_site_count_excluded_from_exported_snv_table"] == 3

    src_pos = set(src_prof["raw_cumulative_variant_reconciliation"]["exported_noncryptic_positions_zero_based"])
    sens_pos = set(sens_prof["raw_cumulative_variant_reconciliation"]["exported_noncryptic_positions_zero_based"])
    assert sorted(sens_pos - src_pos) == [3]
    assert sorted(src_pos - sens_pos) == [22974]

    # 4. Completed host-array sample coverage integration across all 6 hosts / 8 contigs (no pending placeholders)
    host_tbl_path = res_dir / "stage4_host_array_contig_provenance_table.tsv"
    host_tbl = pd.read_csv(host_tbl_path, sep="\t")
    assert len(host_tbl) == 8
    assert not host_tbl.astype(str).apply(lambda s: s.str.contains("NA_pending_")).any().any()
    assert (
        host_tbl["sample_mapping_coverage_status"]
        == "completed_sample_mapping_integrated_via_representative_header"
    ).all()
    focal_host_row = host_tbl[host_tbl["raw_contig_id"] == "41PRr2_1588631"].iloc[0]
    assert float(focal_host_row["sample_mapping_contig_mean_depth"]) == pytest.approx(0.871882, abs=1e-6)
    assert float(focal_host_row["sample_mapping_contig_breadth_pct"]) == pytest.approx(8.149096, abs=1e-6)
    assert float(focal_host_row["sample_mapping_remaining_mag_mean_depth"]) == pytest.approx(13.369586, abs=1e-6)
    assert float(focal_host_row["sample_mapping_remaining_mag_breadth_pct"]) == pytest.approx(71.472404, abs=1e-6)





