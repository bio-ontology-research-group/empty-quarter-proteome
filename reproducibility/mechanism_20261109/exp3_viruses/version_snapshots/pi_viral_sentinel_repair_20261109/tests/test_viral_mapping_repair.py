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
    assert stream_sha256(CODE / "submit_stage4_competitive_mapping_pipeline.sh") == "67e5ee0e881c98f4e8b1af6554fb7e5435a58eb95cd2a3b77d23c4266e936fa0"
    for hist_snap in [
        CODE.parent / "version_snapshots/pi_viral_repair_20261109/code/submit_stage4_competitive_mapping_pipeline.sh",
        CODE.parent / "provenance/version_snapshots/pi_viral_repair_20261109/code/submit_stage4_competitive_mapping_pipeline.sh",
    ]:
        if hist_snap.exists():
            assert stream_sha256(hist_snap) == "9bd3887b4dc91163d9c73650e6e573299a757822f06854f40fa8d88412c86ce3"

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
        "unique_mag_contigs": "1",
        "unique_mag_bases": "348529",
        "unique_filtered_reads_mapped": "6000",
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
    assert res["mapper_raw_conservative_unique_mag_reads"] == 6000

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


