#!/usr/bin/env python3
"""Recover V60PRr2 competitive mapping completion and execute the namespaced-header
inStrain 1.9.1 runtime regression probe without remapping FASTQs or rebuilding the
unified competitive reference.

Root-directed design (EXP3-INSTRAIN-NAMESPACED-HEADER-REPAIR.txt):
1. Verify the version/hash-bound local `SNVprofile.py` overlay patch (`rsplit('::', 1)` /
   `rsplit('::', 2)`) against the immutable `/project/k10018/sw/images/eq-annot/instrain.sif`.
2. Execute a real tiny inStrain 1.9.1 profile + compare regression probe with namespaced
   `::` scaffold headers (`VIRUS_TARGET__TINY_VIRUS::REF_1` and `VIRUS_TARGET__TINY_DISJOINT::REF_2`),
   demonstrating that unpatched inStrain 1.9.1 fails at `SNVprofile.py:766` (`thing.split('::')`)
   while the bind-mounted overlay succeeds with exact terminal completion banners and identical
   SNV/SNS/overlap outcomes as the non-namespaced fixture.
3. Recover sample `V60PRr2` (failed in Slurm array task `15905725_1` at inStrain Step 4 after
   completing Bowtie2 alignment, BAM sorting/indexing, MAPQ>=20 filtering/indexing, full
   3,173,693-contig `samtools coverage`, and focal target evaluation):
   - Independently verify reference validation receipt, 6 `.bt2l` index files, crosswalk, `.fai`,
     and build receipt hashes;
   - Independently verify `V60PRr2` FASTQ inputs and compute their SHA-256 hashes;
   - Independently verify existing `V60PRr2_raw.bam` + `.bai` and `V60PRr2_mapq20.bam` + `.bai`
     via `samtools quickcheck`, `@SQ` header reference lengths, and read counts;
   - Independently verify existing `V60PRr2_filtered_coverage.tsv` and recompute all 7 viral
     targets + 3 boundary intervals + 6 host MAGs to confirm 100% numerical/categorical identity
     with `V60PRr2_focal_targets_evaluation.tsv`;
   - Preserve the failed profile directory `inStrain_source_primary_6b19afac4f46beca` from
     `15905725_1` under `inStrain_source_primary_6b19afac4f46beca.failed_15905725_1_unpatched_split_bug`;
   - Execute a fresh `inStrain profile` for `VIRUS_TARGET__V60PRr2::V60PRr2_715976` with the
     bind-mounted `SNVprofile.py` overlay, verify terminal completion banner and output hashes,
     and write schema-v3 `V60PRr2_mapping_receipt.json` + `MAPPING_DONE` with explicit recovery
     provenance bound to `15905725_1`.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd

from stage4_competitive_common import (
    VIRAL_CATEGORIES,
    atomic_json,
    atomic_text,
    load_analysis_pairs,
    reference_owners,
    require_columns,
    resolve_focal_targets,
    stream_sha256,
    validate_crosswalk_schema,
)
from stage4_map_sample_competitive import (
    FASTQ_DIR,
    SIF_INSTRAIN,
    CommandRunner,
    bam_reference_lengths,
    collect_instrain_callable_positions,
    compute_coverage_metrics,
    compute_interval_coverage,
    evaluate_host_coverage,
    fetch_faidx_sequence,
    require_instrain_available,
    tool_version,
    verify_validation_receipt,
)
from stage4_posthoc_instrain_and_mapping_followup import (
    INSTRAIN_COMPARE_FINISHED_BANNER,
    INSTRAIN_PROFILE_FINISHED_BANNER,
    _verify_focal_fasta_against_target_sha,
    parse_instrain_comparisons_table,
    qualify_instrain_profile_output,
    verify_and_resolve_instrain_namespaced_header_overlay,
)

BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))


def run_namespaced_header_instrain_191_runtime_probe(
    root: Path,
    threads: int = 8,
    instrain_sif: Path = SIF_INSTRAIN,
) -> dict[str, Any]:
    """Run the real inStrain 1.9.1 runtime regression probe with namespaced `::` headers
    (`VIRUS_TARGET__TINY_VIRUS::REF_1` and `VIRUS_TARGET__TINY_DISJOINT::REF_2`).
    """
    singularity_exec = shutil.which("singularity")
    require_instrain_available(instrain_sif, singularity_exec)
    assert singularity_exec is not None

    overlay_info = verify_and_resolve_instrain_namespaced_header_overlay(
        root=root,
        instrain_sif=instrain_sif,
        singularity_exec=singularity_exec,
        verify_container_runtime=True,
    )

    base_probe_dir = root / "scratch_tmp/instrain_191_contract_probe"
    ns_probe_dir = root / "scratch_tmp/instrain_191_namespaced_header_probe"
    preserved_prior_probe_dirs: list[str] = []
    if ns_probe_dir.exists():
        stamp = f"{int(time.time())}_{os.getpid()}"
        preserved_probe = ns_probe_dir.with_name(f"{ns_probe_dir.name}.preserved_{stamp}")
        counter = 1
        while preserved_probe.exists():
            preserved_probe = ns_probe_dir.with_name(f"{ns_probe_dir.name}.preserved_{stamp}_{counter}")
            counter += 1
        os.replace(ns_probe_dir, preserved_probe)
        preserved_prior_probe_dirs.append(str(preserved_probe))
    ns_probe_dir.mkdir(parents=True, exist_ok=True)

    ns_virus_ref = "VIRUS_TARGET__TINY_VIRUS::REF_1"
    ns_disjoint_ref = "VIRUS_TARGET__TINY_DISJOINT::REF_2"

    for fname in ("ref.fna", "sampleA.sam", "sampleB.sam"):
        src = base_probe_dir / fname
        if not src.exists():
            raise FileNotFoundError(f"base contract probe fixture file missing: {src}")
        txt = (
            src.read_text(encoding="utf-8")
            .replace("TINY_VIRUS_REF", ns_virus_ref)
            .replace("TINY_DISJOINT_REF", ns_disjoint_ref)
        )
        (ns_probe_dir / fname).write_text(txt, encoding="utf-8")

    ref_fna = ns_probe_dir / "ref.fna"
    bam_a = ns_probe_dir / "sampleA_mapq20.bam"
    bam_b = ns_probe_dir / "sampleB_mapq20.bam"
    for sam_name, bam_path in [("sampleA.sam", bam_a), ("sampleB.sam", bam_b)]:
        sam_path = ns_probe_dir / sam_name
        subprocess.run(["samtools", "sort", "-o", str(bam_path), str(sam_path)], check=True)
        subprocess.run(["samtools", "index", str(bam_path)], check=True)

    # 1. Demonstrate that unpatched inStrain 1.9.1 fails on namespaced :: headers at Step 4
    unpatched_dir = ns_probe_dir / "profA_unpatched_expected_failure"
    unpatched_cmd = [
        singularity_exec,
        "exec",
        "-B",
        "/project",
        "-B",
        "/scratch",
        str(instrain_sif),
        "inStrain",
        "profile",
        str(bam_a),
        str(ref_fna),
        "-o",
        str(unpatched_dir),
        "-p",
        str(threads),
        "--min_mapq",
        "19",
        "--min_read_ani",
        "0.95",
        "--min_cov",
        "5",
    ]
    unpatched_proc = subprocess.run(unpatched_cmd, capture_output=True, text=True)
    if unpatched_proc.returncode == 0:
        raise RuntimeError("expected unpatched inStrain 1.9.1 to fail on namespaced :: header, but it returned 0")
    if "scaff, mm = thing.split('::')" not in unpatched_proc.stderr or "too many values to unpack (expected 2)" not in unpatched_proc.stderr:
        raise RuntimeError(f"unexpected failure stderr from unpatched inStrain 1.9.1: {unpatched_proc.stderr}")

    # 2. Run patched inStrain profile on sampleA (0 variants) and sampleB (1 SNS + 1 SNV)
    prof_a_dir = ns_probe_dir / "profA_patched"
    prof_b_dir = ns_probe_dir / "profB_patched"
    for bam_path, out_dir in [(bam_a, prof_a_dir), (bam_b, prof_b_dir)]:
        cmd = [
            singularity_exec,
            "exec",
            "-B",
            "/project",
            "-B",
            "/scratch",
            "-B",
            overlay_info["singularity_bind_arg"],
            str(instrain_sif),
            "inStrain",
            "profile",
            str(bam_path),
            str(ref_fna),
            "-o",
            str(out_dir),
            "-p",
            str(threads),
            "--min_mapq",
            "19",
            "--min_read_ani",
            "0.95",
            "--min_cov",
            "5",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"patched inStrain profile failed for {bam_path} (rc={proc.returncode}): {proc.stderr}")

    qual_a_virus = qualify_instrain_profile_output(prof_a_dir, ns_virus_ref, min_cov=5)
    qual_b_virus = qualify_instrain_profile_output(prof_b_dir, ns_virus_ref, min_cov=5)

    if (
        qual_a_virus["snv_count_scaffold_info"] != 0
        or qual_a_virus["sns_count_scaffold_info"] != 0
        or qual_a_virus["snv_table_status"] != "present_empty_dataframe_newline_verified_zero_variants_with_scaffold_info_and_log"
        or qual_a_virus["inferred_depth_qualified_positions_ge_mincov"] != 790
    ):
        raise ValueError(f"unexpected patched profA qualification metrics: {qual_a_virus}")

    if (
        qual_b_virus["snv_count_scaffold_info"] != 1
        or qual_b_virus["sns_count_scaffold_info"] != 1
        or qual_b_virus["divergent_site_count_scaffold_info"] != 2
        or qual_b_virus["snv_target_row_count"] != 2
        or qual_b_virus["inferred_depth_qualified_positions_ge_mincov"] != 790
    ):
        raise ValueError(f"unexpected patched profB qualification metrics: {qual_b_virus}")

    # 3. Run patched inStrain compare on profA_patched and profB_patched
    comp_dir = ns_probe_dir / "compAB_patched"
    comp_cmd = [
        singularity_exec,
        "exec",
        "-B",
        "/project",
        "-B",
        "/scratch",
        "-B",
        overlay_info["singularity_bind_arg"],
        str(instrain_sif),
        "inStrain",
        "compare",
        "-i",
        str(prof_a_dir),
        str(prof_b_dir),
        "-o",
        str(comp_dir),
        "-p",
        str(threads),
    ]
    comp_proc = subprocess.run(comp_cmd, capture_output=True, text=True)
    if comp_proc.returncode != 0:
        raise RuntimeError(f"patched inStrain compare failed (rc={comp_proc.returncode}): {comp_proc.stderr}")

    qual_profiles_meta = [
        {"sample": "sampleA", "profile_dir": str(prof_a_dir), "instrain_length_bp": 1000},
        {"sample": "sampleB", "profile_dir": str(prof_b_dir), "instrain_length_bp": 1000},
    ]
    parsed_comp_virus = parse_instrain_comparisons_table(comp_dir, ns_virus_ref, qual_profiles_meta)
    parsed_comp_disjoint = parse_instrain_comparisons_table(comp_dir, ns_disjoint_ref, qual_profiles_meta)

    v_row = parsed_comp_virus["comparison_rows"][0]
    d_row = parsed_comp_disjoint["comparison_rows"][0]
    if (
        v_row["compared_bases_count"] != 790
        or v_row["consensus_SNPs"] != 2
        or v_row["population_SNPs"] != 1
        or abs(float(v_row["popANI"]) - (789 / 790)) > 1e-6
        or abs(float(v_row["conANI"]) - (788 / 790)) > 1e-6
    ):
        raise ValueError(f"unexpected patched compare positive-overlap metrics: {v_row}")
    if (
        d_row["compared_bases_count"] != 0
        or d_row["conANI"] != "NA_zero_coverage_overlap"
        or d_row["popANI"] != "NA_zero_coverage_overlap"
        or d_row["pairwise_estimability_status"] != "non_estimable_zero_coverage_overlap"
    ):
        raise ValueError(f"unexpected patched compare zero-overlap metrics: {d_row}")

    prof_a_scaff = prof_a_dir / "output/profA_patched_scaffold_info.tsv"
    prof_a_snvs = prof_a_dir / "output/profA_patched_SNVs.tsv"
    prof_a_log = prof_a_dir / "log/log.log"
    prof_b_scaff = prof_b_dir / "output/profB_patched_scaffold_info.tsv"
    prof_b_snvs = prof_b_dir / "output/profB_patched_SNVs.tsv"
    prof_b_log = prof_b_dir / "log/log.log"
    comp_table = comp_dir / "output/compAB_patched_comparisonsTable.tsv"
    comp_log = comp_dir / "log/log.log"

    probe_receipt = {
        "schema_version": 1,
        "status": "completed_verified_instrain_1.9.1_namespaced_header_runtime_probe",
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "instrain_version": tool_version([singularity_exec, "exec", str(instrain_sif), "inStrain", "profile", "--version"]),
        "container_path": str(instrain_sif),
        "instrain_namespaced_header_overlay": overlay_info,
        "probe_dir_on_shaheen": str(ns_probe_dir),
        "preserved_prior_probe_dirs": preserved_prior_probe_dirs,
        "namespaced_scaffolds": {
            "positive_overlap_scaffold": ns_virus_ref,
            "zero_overlap_scaffold": ns_disjoint_ref,
        },
        "inputs": {
            "ref_fna_sha256": stream_sha256(ref_fna),
            "sampleA_mapq20_bam_sha256": stream_sha256(bam_a),
            "sampleB_mapq20_bam_sha256": stream_sha256(bam_b),
        },
        "unpatched_failure_regression": {
            "returncode": unpatched_proc.returncode,
            "failed_at_snvprofile_line_766_split_bug": True,
            "stderr_excerpt": "\n".join(unpatched_proc.stderr.strip().splitlines()[-12:]),
        },
        "profA_patched_zero_variants": {
            "log_has_profile_finished_banner": INSTRAIN_PROFILE_FINISHED_BANNER in prof_a_log.read_text(encoding="utf-8"),
            "log_sha256": stream_sha256(prof_a_log),
            "scaffold_info_sha256": stream_sha256(prof_a_scaff),
            "scaffold_info_tsv": prof_a_scaff.read_text(encoding="utf-8"),
            "snvs_tsv_bytes": prof_a_snvs.stat().st_size,
            "snvs_tsv_raw_repr": repr(prof_a_snvs.read_text(encoding="utf-8")),
            "snvs_tsv_sha256": stream_sha256(prof_a_snvs),
            "qualified_summary": qual_a_virus,
        },
        "profB_patched_one_sns_one_snv": {
            "log_has_profile_finished_banner": INSTRAIN_PROFILE_FINISHED_BANNER in prof_b_log.read_text(encoding="utf-8"),
            "log_sha256": stream_sha256(prof_b_log),
            "scaffold_info_sha256": stream_sha256(prof_b_scaff),
            "scaffold_info_tsv": prof_b_scaff.read_text(encoding="utf-8"),
            "snvs_tsv": prof_b_snvs.read_text(encoding="utf-8"),
            "snvs_tsv_sha256": stream_sha256(prof_b_snvs),
            "qualified_summary": qual_b_virus,
        },
        "compAB_patched_positive_and_zero_overlap": {
            "log_has_compare_finished_banner": INSTRAIN_COMPARE_FINISHED_BANNER in comp_log.read_text(encoding="utf-8"),
            "log_sha256": stream_sha256(comp_log),
            "comparisons_table_sha256": stream_sha256(comp_table),
            "comparisons_table_tsv": comp_table.read_text(encoding="utf-8"),
            "parsed_positive_overlap": parsed_comp_virus,
            "parsed_zero_overlap": parsed_comp_disjoint,
        },
    }
    out_probe_json = root / "results/stage4_instrain_191_namespaced_header_runtime_probe.json"
    atomic_json(out_probe_json, probe_receipt)
    print(f"[Namespaced-header inStrain 1.9.1 probe verified] {out_probe_json}", flush=True)
    return probe_receipt


def _evaluations_match_table(recomputed: list[dict[str, Any]], existing_df: pd.DataFrame) -> None:
    """Verify that recomputed focal evaluation rows match existing `<sample>_focal_targets_evaluation.tsv`."""
    recomp_df = pd.DataFrame(recomputed)
    if list(recomp_df.columns) != list(existing_df.columns):
        raise ValueError(f"focal evaluation column mismatch: {list(recomp_df.columns)} != {list(existing_df.columns)}")
    if len(recomp_df) != len(existing_df):
        raise ValueError(f"focal evaluation row count mismatch: {len(recomp_df)} != {len(existing_df)}")
    for idx in range(len(recomp_df)):
        r_new = recomp_df.iloc[idx]
        r_old = existing_df.iloc[idx]
        for col in recomp_df.columns:
            v_new = r_new[col]
            v_old = r_old[col]
            if pd.isna(v_new) and pd.isna(v_old):
                continue
            if isinstance(v_new, (float, int)) and not isinstance(v_new, bool) and pd.notna(v_old):
                try:
                    f_new = float(v_new)
                    f_old = float(v_old)
                    if math.isclose(f_new, f_old, rel_tol=1e-5, abs_tol=1e-7):
                        continue
                except (ValueError, TypeError):
                    pass
            if str(v_new).strip() != str(v_old).strip():
                raise ValueError(
                    f"focal evaluation mismatch at row {idx} col {col!r}: recomputed={v_new!r} vs existing={v_old!r}"
                )


def recover_v60prr2_sample_mapping(
    root: Path,
    threads: int = 32,
    fastq_dir: Path = FASTQ_DIR,
    instrain_sif: Path = SIF_INSTRAIN,
) -> dict[str, Any]:
    """Recover V60PRr2 sample mapping completion from verified existing BAMs and coverage TSV."""
    start = time.time()
    sample = "V60PRr2"
    res_dir = root / "results"
    out_dir = res_dir / "stage4_competitive_mapping" / sample
    ref_dir = res_dir / "stage4_competitive_reference"
    ref_fna = ref_dir / "unified_competitive_reference.fna"
    crosswalk_path = ref_dir / "stage4_reference_sequence_crosswalk.tsv"
    pair_path = res_dir / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    boundary_config = root / "inputs/stage4_boundary_sensitive_intervals.json"
    manifest_path = root / "inputs/genome_manifest.tsv"
    original_mapper_script = root / "code/stage4_map_sample_competitive.py"
    common_script = root / "code/stage4_competitive_common.py"

    done_path = out_dir / "MAPPING_DONE"
    done_path.unlink(missing_ok=True)
    runner = CommandRunner()

    for tool in ["bowtie2", "samtools"]:
        if shutil.which(tool) is None:
            raise RuntimeError(f"required mapping tool unavailable: {tool}")
    singularity_exec = shutil.which("singularity")
    require_instrain_available(instrain_sif, singularity_exec)
    assert singularity_exec is not None

    overlay_info = verify_and_resolve_instrain_namespaced_header_overlay(
        root=root,
        instrain_sif=instrain_sif,
        singularity_exec=singularity_exec,
        verify_container_runtime=True,
    )
    instrain_version = tool_version([singularity_exec, "exec", str(instrain_sif), "inStrain", "profile", "--version"])
    versions = {
        "bowtie2": tool_version(["bowtie2", "--version"]),
        "samtools": tool_version(["samtools", "--version"]),
        "instrain": instrain_version,
    }

    validation_receipt, validation_binding = verify_validation_receipt(ref_dir)
    print("[Recovery Preflight] Reference, six index files, crosswalk, FAI and build receipt hashes verified", flush=True)

    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_path, boundary_config)
    crosswalk = pd.read_csv(crosswalk_path, sep="\t")
    validate_crosswalk_schema(crosswalk)
    owners = reference_owners(crosswalk)
    resolved = resolve_focal_targets(
        pairs,
        crosswalk,
        lambda reference: fetch_faidx_sequence(runner, ref_fna, reference),
    )
    manifest = pd.read_csv(manifest_path, sep="\t")
    require_columns(manifest, ["mag", "sha256", "contigs", "bases"], "genome manifest")
    manifest = manifest.set_index("mag", drop=False)
    for host in sorted(set(pairs["host_mag"])):
        if host not in manifest.index:
            raise ValueError(f"focal host absent from canonical manifest: {host}")

    # Verify original FASTQ inputs
    r1 = fastq_dir / f"{sample}_R1.fastq.gz"
    r2 = fastq_dir / f"{sample}_R2.fastq.gz"
    se = fastq_dir / f"{sample}_SE.fastq.gz"
    if not r1.exists() or not r2.exists():
        raise FileNotFoundError(f"paired FASTQ missing for {sample}: {r1} / {r2}")
    read_paths = {"r1": r1, "r2": r2}
    if se.exists():
        read_paths["se"] = se
    read_hashes = {f"{key}_sha256": stream_sha256(path) for key, path in read_paths.items()}
    read_hashes["se_sha256"] = read_hashes.get("se_sha256", "not_present")

    # Verify existing BAMs, BAI indexes, and original alignment stderr logs from 15905725_1
    raw_bam = out_dir / f"{sample}_raw.bam"
    raw_bai = out_dir / f"{sample}_raw.bam.bai"
    filtered_bam = out_dir / f"{sample}_mapq20.bam"
    filtered_bai = out_dir / f"{sample}_mapq20.bam.bai"
    coverage_path = out_dir / f"{sample}_filtered_coverage.tsv"
    evaluation_path = out_dir / f"{sample}_focal_targets_evaluation.tsv"
    for req_file in [raw_bam, raw_bai, filtered_bam, filtered_bai, coverage_path, evaluation_path]:
        if not req_file.exists() or req_file.stat().st_size == 0:
            raise FileNotFoundError(f"required existing V60PRr2 artifact missing or empty: {req_file}")

    runner.run(["samtools", "quickcheck", str(raw_bam)], "quickcheck existing raw BAM")
    runner.run(["samtools", "quickcheck", str(filtered_bam)], "quickcheck existing filtered BAM")

    raw_total = int(runner.run(["samtools", "view", "-c", str(raw_bam)], "count all raw records").stdout.strip())
    raw_primary = int(runner.run(["samtools", "view", "-c", "-F", "2820", str(raw_bam)], "count primary mapped raw records").stdout.strip())
    filtered_count = int(runner.run(["samtools", "view", "-c", str(filtered_bam)], "count filtered records").stdout.strip())
    header = runner.run(["samtools", "view", "-H", str(filtered_bam)], "read filtered BAM header").stdout
    bam_lengths = bam_reference_lengths(header)
    if set(bam_lengths) != set(owners):
        raise ValueError("BAM references differ from crosswalk representatives")
    for name, length in bam_lengths.items():
        if owners[name]["length_bp"] != length:
            raise ValueError(f"BAM/crosswalk length mismatch for {name}")

    # Verify existing coverage TSV against BAM header and recompute all focal evaluations
    coverage_frame = compute_coverage_metrics(coverage_path.read_text(encoding="utf-8"), bam_lengths)
    coverage = {row["rname"]: row for row in coverage_frame.to_dict(orient="records")}

    evaluations: list[dict[str, Any]] = []
    for target_sha, group in pairs.groupby("target_seq_sha256", sort=True):
        info = resolved[target_sha]
        representative = info["reference_contig"]
        entry = coverage[representative]
        categories = owners[representative]["categories"]
        unique = bool(categories) and categories <= VIRAL_CATEGORIES
        if "host_mag" in categories:
            ambiguity = "shared_host_viral_ambiguous"
        elif "all_plasmids" in categories:
            ambiguity = "shared_viral_plasmid_ambiguous"
        elif "cellular_controls" in categories:
            ambiguity = "shared_viral_decoy_ambiguous"
        elif unique:
            ambiguity = "unique_viral_sequence_aliases_collapsed"
        else:
            ambiguity = "unresolved_nonviral_ownership"
        roles = set(group["analysis_role"])
        if len(roles) != 1:
            raise ValueError(f"conflicting analysis roles for {target_sha}: {roles}")
        role = next(iter(roles))
        depth = float(entry["meandepth"])
        breadth = 100.0 * int(entry["covbases"]) / info["length_bp"]
        qualifies = bool(depth >= 10 and breadth >= 50 and unique and role == "primary")
        if role == "boundary_sensitive_provisional":
            unique = False
            ambiguity = "boundary_sensitive_provisional; whole contig is not claimed viral; see interval rows"
            qualifies = False
        reads = int(entry["numreads"])
        evaluations.append({
            "sample": sample,
            "entity_id": group.iloc[0]["target_id"],
            "entity_type": "virus_target",
            "entity_sha256": target_sha,
            "raw_oriented_target_sha256": target_sha,
            "canonical_sha256": info["canonical_sha256"],
            "target_aliases": ";".join(info["target_aliases"]),
            "host_mags": ";".join(sorted(set(group["host_mag"]))),
            "reference_contig": representative,
            "reference_orientation_to_ledger": info["reference_orientation_to_ledger"],
            "length_bp": info["length_bp"],
            "filtered_reads_mapped": reads,
            "unique_filtered_reads_mapped": reads if unique else 0,
            "ambiguous_filtered_reads_mapped": 0 if unique else reads,
            "horizontal_breadth_pct": breadth,
            "mean_depth": depth,
            "ambiguity_status": ambiguity,
            "uniquely_attributable": unique,
            "analysis_role": role,
            "qualifies_10x_50pct_instrain": int(qualifies),
            "instrain_action": "eligible" if qualifies else "not_eligible_coverage_ambiguity_or_sensitivity",
            "individual_status": group.iloc[0]["individual_virality_status"],
        })

    for candidate in boundary_pairs.to_dict(orient="records"):
        target_sha = candidate["target_seq_sha256"]
        representative = resolved[target_sha]["reference_contig"]
        catalog_start, catalog_end = candidate["catalog_viral_interval"]
        for nested_id in candidate["nested_competitor_ids"]:
            nested_rows = crosswalk[crosswalk["entity_id"] == nested_id]
            if nested_rows.empty:
                raise ValueError(f"nested boundary competitor absent from crosswalk: {nested_id}")
            if not set(nested_rows["representative_header"]) <= set(bam_lengths):
                raise ValueError(f"nested boundary competitor absent from BAM: {nested_id}")
        for interval in candidate["boundary_intervals"]:
            start_position, end_position = int(interval["start"]), int(interval["end"])
            if resolved[target_sha]["reference_orientation_to_ledger"] == "reverse_complement_of_reference":
                reference_start = resolved[target_sha]["length_bp"] - end_position + 1
                reference_end = resolved[target_sha]["length_bp"] - start_position + 1
            else:
                reference_start, reference_end = start_position, end_position
            region = f"{representative}:{reference_start}-{reference_end}"
            depth_result = runner.run([
                "samtools", "depth", "-a", "-d", "0", "-Q", "20", "-r", region, str(filtered_bam),
            ], f"boundary interval depth {region}")
            metrics = compute_interval_coverage(depth_result.stdout, representative, reference_start, reference_end)
            interval_reads = int(runner.run([
                "samtools", "view", "-c", "-q", "20", "-F", "3844", str(filtered_bam), region,
            ], f"boundary interval reads {region}").stdout.strip())
            nested_ambiguous = start_position <= catalog_end and end_position >= catalog_start
            evaluations.append({
                "sample": sample,
                "entity_id": f"{candidate['target_id']}|{interval['name']}",
                "entity_type": "boundary_interval",
                "entity_sha256": target_sha,
                "host_mag": candidate["host_mag"],
                "host_mags": candidate["host_mag"],
                "reference_contig": representative,
                "interval_name": interval["name"],
                "interval_start_1based": start_position,
                "interval_end_1based_inclusive": end_position,
                "reference_interval_start_1based": reference_start,
                "reference_interval_end_1based_inclusive": reference_end,
                "length_bp": metrics["length_bp"],
                "filtered_reads_mapped": interval_reads,
                "unique_filtered_reads_mapped": 0,
                "ambiguous_filtered_reads_mapped": interval_reads,
                "horizontal_breadth_pct": metrics["horizontal_breadth_pct"],
                "mean_depth": metrics["mean_depth"],
                "callable_positions_depth_ge1": metrics["callable_positions_depth_ge1"],
                "nested_competitor_ambiguous": nested_ambiguous,
                "nested_competitor_ids": ";".join(candidate["nested_competitor_ids"]),
                "ambiguity_status": (
                    "nested_catalog_core_competitor; low MAPQ-filtered coverage is not absence"
                    if nested_ambiguous else "boundary-sensitive interval; not uniquely viral-attributable"
                ),
                "uniquely_attributable": False,
                "analysis_role": "boundary_sensitive_provisional",
                "qualifies_10x_50pct_instrain": 0,
                "instrain_action": "not_eligible_boundary_sensitive_provisional",
                "individual_status": candidate["individual_virality_status"],
                "interpretation": interval["interpretation"],
            })

    for host in sorted(set(pairs["host_mag"])):
        host_rows = crosswalk[(crosswalk["source_category"] == "host_mag") & (crosswalk["entity_id"] == host)]
        row = evaluate_host_coverage(
            host, host_rows, owners, coverage,
            int(manifest.loc[host, "contigs"]), int(manifest.loc[host, "bases"]),
            str(manifest.loc[host, "sha256"]),
        )
        row["sample"] = sample
        row["analysis_role"] = "host_context"
        evaluations.append(row)

    existing_eval_df = pd.read_csv(evaluation_path, sep="\t")
    _evaluations_match_table(evaluations, existing_eval_df)
    print(f"[Recovery Verified] Existing {evaluation_path.name} matches recomputed 7 viral + 3 boundary + 6 host rows", flush=True)

    # Preserve failed profile directory from 15905725_1 before running patched profile
    min_read_ani = 0.95
    instrain_results: list[dict[str, Any]] = []
    preserved_failed_attempts: list[dict[str, Any]] = []
    qualifying = [row for row in evaluations if row["entity_type"] == "virus_target" and row["qualifies_10x_50pct_instrain"]]
    if len(qualifying) != 1:
        raise RuntimeError(f"expected exactly 1 qualifying primary viral target in V60PRr2, found={len(qualifying)}")

    for target in qualifying:
        target_sha = str(target["entity_sha256"])
        representative = str(target["reference_contig"])
        focal_fasta = out_dir / f"focal_{target_sha[:16]}.fna"
        if not focal_fasta.exists():
            sequence = fetch_faidx_sequence(runner, ref_fna, representative)
            atomic_text(focal_fasta, f">{representative}\n{sequence}\n")
        focal_fasta_sha = _verify_focal_fasta_against_target_sha(focal_fasta, representative, target_sha)

        profile_dir = out_dir / f"inStrain_source_primary_{target_sha[:16]}"
        preserved_dir = out_dir / f"inStrain_source_primary_{target_sha[:16]}.failed_15905725_1_unpatched_split_bug"
        if preserved_dir.exists():
            failed_hashes = {
                str(p.relative_to(preserved_dir)): stream_sha256(p)
                for p in sorted(preserved_dir.rglob("*")) if p.is_file()
            }
            preserved_failed_attempts.append({
                "original_profile_dir": str(profile_dir),
                "preserved_failed_profile_dir": str(preserved_dir),
                "source_slurm_job_id": "15905725_1",
                "failure_reason": (
                    "inStrain 1.9.1 SNVprofile.py _load_special line 766 `scaff, mm = thing.split('::')` "
                    "raised `ValueError: too many values to unpack (expected 2)` during Step 4 "
                    "`profile_genome_wide` on namespaced reference header `VIRUS_TARGET__V60PRr2::V60PRr2_715976`."
                ),
                "preserved_file_sha256": failed_hashes,
            })
        if profile_dir.exists():
            log_file = profile_dir / "log/log.log"
            has_banner = log_file.exists() and (INSTRAIN_PROFILE_FINISHED_BANNER in log_file.read_text(encoding="utf-8", errors="replace"))
            if not has_banner:
                if preserved_dir.exists():
                    stamp = f"{int(time.time())}_{os.getpid()}"
                    alt_preserved = out_dir / f"inStrain_source_primary_{target_sha[:16]}.failed_preserved_{stamp}"
                    os.replace(profile_dir, alt_preserved)
                    preserved_target = alt_preserved
                else:
                    os.replace(profile_dir, preserved_dir)
                    preserved_target = preserved_dir
                failed_hashes = {
                    str(p.relative_to(preserved_target)): stream_sha256(p)
                    for p in sorted(preserved_target.rglob("*")) if p.is_file()
                }
                preserved_failed_attempts.append({
                    "original_profile_dir": str(profile_dir),
                    "preserved_failed_profile_dir": str(preserved_target),
                    "source_slurm_job_id": "15905725_1",
                    "failure_reason": (
                        "inStrain 1.9.1 SNVprofile.py _load_special line 766 `scaff, mm = thing.split('::')` "
                        "raised `ValueError: too many values to unpack (expected 2)` during Step 4 "
                        "`profile_genome_wide` on namespaced reference header `VIRUS_TARGET__V60PRr2::V60PRr2_715976`."
                    ),
                    "preserved_file_sha256": failed_hashes,
                })
            else:
                # Archive prior completed patched profile (e.g. from job 15910426) before fresh run to ensure clean single provenance
                stamp = f"{int(time.time())}_{os.getpid()}"
                prior_completed = out_dir / f"inStrain_source_primary_{target_sha[:16]}.prior_completed_{stamp}"
                os.replace(profile_dir, prior_completed)
                prior_hashes = {
                    str(p.relative_to(prior_completed)): stream_sha256(p)
                    for p in sorted(prior_completed.rglob("*")) if p.is_file()
                }
                log_15910426 = root / "logs/stage4_recover_v60prr2_15910426.log"
                err_15910426 = root / "logs/stage4_recover_v60prr2_15910426.err"
                preserved_failed_attempts.append({
                    "original_profile_dir": str(profile_dir),
                    "preserved_failed_profile_dir": str(prior_completed),
                    "source_slurm_job_id": "15910426",
                    "failure_reason": (
                        "Patched inStrain 1.9.1 profile completed (`..:: inStrain profile finished ::..`), "
                        "but pre-reconciliation posthoc qualification failed with `ValueError: SNVs.tsv non-cryptic row count (22) "
                        "for VIRUS_TARGET__V60PRr2::V60PRr2_715976 does not match scaffold_info divergent_site_count (24)` "
                        "before raw cumulative variant reconciliation (22 non-cryptic exported + 2 cryptic excluded = 24 summary sites)."
                    ),
                    "slurm_stdout_log": str(log_15910426) if log_15910426.exists() else "missing",
                    "slurm_stdout_log_sha256": stream_sha256(log_15910426) if log_15910426.exists() else None,
                    "slurm_stderr_log": str(err_15910426) if err_15910426.exists() else "missing",
                    "slurm_stderr_log_sha256": stream_sha256(err_15910426) if err_15910426.exists() else None,
                    "preserved_file_sha256": prior_hashes,
                })

        runner.run([
            singularity_exec,
            "exec",
            "-B",
            "/project",
            "-B",
            "/scratch",
            "-B",
            overlay_info["singularity_bind_arg"],
            str(instrain_sif),
            "inStrain",
            "profile",
            str(filtered_bam),
            str(focal_fasta),
            "-o",
            str(profile_dir),
            "-p",
            str(threads),
            "--min_mapq",
            "19",
            "--min_read_ani",
            str(min_read_ani),
        ], f"inStrain source_primary profile {target_sha} (with namespaced-header rsplit overlay)")

        # Strictly verify completion banner and qualify output
        qual_posthoc = qualify_instrain_profile_output(
            profile_dir,
            representative,
            min_cov=5,
            min_cov_provenance="instrain_1.9.1_verified_default_min_cov_5_from_argumentParser",
        )
        callable_report = collect_instrain_callable_positions(profile_dir, representative)
        output_hashes = {
            str(path.relative_to(profile_dir)): stream_sha256(path)
            for path in sorted(profile_dir.rglob("*")) if path.is_file()
        }
        instrain_results.append({
            "target_raw_sha256": target_sha,
            "reference_contig": representative,
            "profile_dir": str(profile_dir),
            "status": "profiled_successfully",
            "tool_version": instrain_version,
            "min_mapq_strict_greater_than": 19,
            "effective_minimum_integer_mapq": 20,
            "min_read_ani_strict_greater_than": min_read_ani,
            "analysis_mode": "source_primary",
            "focal_fasta_sha256": focal_fasta_sha,
            "instrain_namespaced_header_overlay": overlay_info,
            "posthoc_qualified_metrics": qual_posthoc,
            **callable_report,
            "output_sha256": output_hashes,
        })

    # Compute hashes of existing BAMs, coverage TSV, evaluation TSV, and original logs from 15905725_1
    raw_bam_sha = stream_sha256(raw_bam)
    raw_bai_sha = stream_sha256(raw_bai)
    filtered_bam_sha = stream_sha256(filtered_bam)
    filtered_bai_sha = stream_sha256(filtered_bai)
    coverage_tsv_sha = stream_sha256(coverage_path)
    focal_eval_tsv_sha = stream_sha256(evaluation_path)

    orig_log_path = root / "logs/stage4_map_15905725_1.log"
    orig_err_path = root / "logs/stage4_map_15905725_1.err"
    bt2_stderr_path = out_dir / "alignment_bowtie2.stderr.log"
    sort_stderr_path = out_dir / "alignment_samtools_sort.stderr.log"
    view_stderr_path = out_dir / "alignment_samtools_view.stderr.log"

    recovery_provenance = {
        "schema_version": 1,
        "status": "completed_verified_v60prr2_instrain_recovery",
        "sample": sample,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "source_mapping_job": {
            "slurm_job_id": "15905725_1",
            "slurm_array_job_id": "15905725",
            "slurm_array_task_id": "1",
            "state": "FAILED",
            "exit_code": "1:0",
            "elapsed": "01:11:23",
            "completed_stages_in_source_job": [
                "Preflight reference/index/crosswalk/build-receipt hash verification",
                "Bowtie2 --very-sensitive alignment streamed to sorted V60PRr2_raw.bam + .bai",
                "samtools quickcheck on V60PRr2_raw.bam",
                "MAPQ>=20 and flag 3844 filtering to V60PRr2_mapq20.bam + .bai",
                "samtools quickcheck on V60PRr2_mapq20.bam",
                f"Full {len(coverage):,}-contig samtools coverage to {sample}_filtered_coverage.tsv",
                "Boundary interval depth/read evaluations and V60PRr2_focal_targets_evaluation.tsv generation",
                "samtools faidx extraction of focal_6b19afac4f46beca.fna",
                "inStrain profile Step 1 (Filter reads) and Step 2 (Profile scaffolds)",
            ],
            "failed_stage_in_source_job": (
                "inStrain profile Step 4 (profile_genome_wide): SNVprofile.py line 766 "
                "`scaff, mm = thing.split('::')` raised `ValueError: too many values to unpack (expected 2)` "
                "on namespaced reference contig `VIRUS_TARGET__V60PRr2::V60PRr2_715976`."
            ),
            "original_mapping_script_path": str(original_mapper_script),
            "original_mapping_script_sha256": stream_sha256(original_mapper_script),
            "original_common_script_sha256": stream_sha256(common_script),
            "original_slurm_stdout_log": str(orig_log_path) if orig_log_path.exists() else "missing",
            "original_slurm_stdout_log_sha256": stream_sha256(orig_log_path) if orig_log_path.exists() else None,
            "original_slurm_stderr_log": str(orig_err_path) if orig_err_path.exists() else "missing",
            "original_slurm_stderr_log_sha256": stream_sha256(orig_err_path) if orig_err_path.exists() else None,
            "alignment_bowtie2_stderr_sha256": stream_sha256(bt2_stderr_path) if bt2_stderr_path.exists() else None,
            "alignment_samtools_sort_stderr_sha256": stream_sha256(sort_stderr_path) if sort_stderr_path.exists() else None,
            "alignment_samtools_view_stderr_sha256": stream_sha256(view_stderr_path) if view_stderr_path.exists() else None,
        },
        "recovery_execution": {
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", "local"),
            "hostname": platform.node(),
            "recovery_script_path": str(Path(__file__).resolve()),
            "recovery_script_sha256": stream_sha256(Path(__file__).resolve()),
            "remapped_fastqs": False,
            "rebuilt_reference": False,
            "verified_existing_bams_and_coverage": True,
            "instrain_namespaced_header_overlay": overlay_info,
            "preserved_failed_attempts": preserved_failed_attempts,
            "commands_executed_in_recovery": runner.commands,
        },
        "verified_artifact_sha256": {
            "raw_bam_sha256": raw_bam_sha,
            "raw_bai_sha256": raw_bai_sha,
            "filtered_bam_sha256": filtered_bam_sha,
            "filtered_bai_sha256": filtered_bai_sha,
            "coverage_tsv_sha256": coverage_tsv_sha,
            "focal_evaluation_tsv_sha256": focal_eval_tsv_sha,
        },
    }
    recovery_prov_path = out_dir / f"{sample}_recovery_provenance.json"
    atomic_json(recovery_prov_path, recovery_provenance)

    receipt = {
        "status": "completed_verified",
        "schema_version": 3,
        "stage": "stage4_map_sample_competitive",
        "sample": sample,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": time.time() - start,
        "execution_environment": {
            "hostname": platform.node(),
            "kernel": platform.release(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", "local"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", "not_array"),
            "source_mapping_slurm_job_id": "15905725_1",
            "source_mapping_slurm_array_job_id": "15905725",
            "source_mapping_slurm_array_task_id": "1",
            "source_mapping_slurm_state": "FAILED",
            "source_mapping_slurm_exit_code": "1:0",
            "python_version": sys.version.split()[0],
            "tool_versions": versions,
        },
        "inputs": {
            "r1": str(r1),
            "r1_bytes": r1.stat().st_size,
            "r1_sha256": read_hashes["r1_sha256"],
            "r2": str(r2),
            "r2_bytes": r2.stat().st_size,
            "r2_sha256": read_hashes["r2_sha256"],
            "se": str(se) if se.exists() else "not_present",
            "se_sha256": read_hashes["se_sha256"],
            "pair_ledger": str(pair_path),
            "pair_ledger_sha256": stream_sha256(pair_path),
            "boundary_config": str(boundary_config),
            "boundary_config_sha256": stream_sha256(boundary_config),
            "mapping_script": str(original_mapper_script),
            "mapping_script_sha256": stream_sha256(original_mapper_script),
            "recovery_script": str(Path(__file__).resolve()),
            "recovery_script_sha256": stream_sha256(Path(__file__).resolve()),
            "common_script_sha256": stream_sha256(common_script),
            "manifest_sha256": stream_sha256(manifest_path),
        },
        "reference_binding": {
            "validation_receipt_sha256": validation_binding["receipt_sha256"],
            "validation_receipt_path": validation_binding["receipt_path"],
            "reference_fna_sha256": validation_receipt["reference_files"]["reference_fna"]["sha256"],
            "crosswalk_sha256": validation_receipt["reference_files"]["crosswalk_tsv"]["sha256"],
            "build_receipt_sha256": validation_receipt["build_receipt"]["sha256"],
            "bowtie2_index_sha256": {row["name"]: row["sha256"] for row in validation_receipt["bowtie2_index"]["files"]},
        },
        "commands": runner.commands,
        "mapping_metrics": {
            "raw_total_records": raw_total,
            "raw_primary_mapped_records": raw_primary,
            "filtered_high_quality_mapq20_reads": filtered_count,
            "filtered_fraction_of_primary": filtered_count / raw_primary if raw_primary else None,
            "filter": "MAPQ>=20; exclude UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
            "coverage_references_verified": len(coverage),
        },
        "focal_pair_scope": {
            "n_primary_pairs": int((reviewed_pairs["analysis_role"] == "primary").sum()),
            "n_sensitivity_pairs": int((reviewed_pairs["analysis_role"] != "primary").sum()),
            "n_boundary_sensitive_provisional_pairs": len(boundary_pairs),
            "n_unique_viral_sha": len(resolved),
        },
        "instrain": {
            "performed": bool(instrain_results),
            "status": "profiles_completed" if instrain_results else "not_run_no_eligible_targets",
            "analysis_mode": "source_primary",
            "min_read_ani_strict_greater_than": min_read_ani,
            "relaxed_crosssample_sensitivity_is_separate": True,
            "instrain_namespaced_header_overlay": overlay_info,
            "results": instrain_results,
            "recombination_analysis": "not_performed; SNVs/sequence polymorphism are not recombination evidence",
        },
        "recovery_provenance": recovery_provenance,
        "recovery_provenance_path": str(recovery_prov_path),
        "recovery_provenance_sha256": stream_sha256(recovery_prov_path),
        "focal_evaluations": evaluations,
        "outputs": {
            "raw_bam": str(raw_bam),
            "raw_bam_bytes": raw_bam.stat().st_size,
            "raw_bam_sha256": raw_bam_sha,
            "filtered_bam": str(filtered_bam),
            "filtered_bam_sha256": filtered_bam_sha,
            "coverage_tsv": str(coverage_path),
            "coverage_tsv_sha256": coverage_tsv_sha,
            "focal_evaluation_tsv": str(evaluation_path),
            "focal_evaluation_tsv_sha256": focal_eval_tsv_sha,
        },
    }
    mapping_receipt_path = out_dir / f"{sample}_mapping_receipt.json"
    atomic_json(mapping_receipt_path, receipt)
    atomic_text(done_path, json.dumps({"receipt_sha256": stream_sha256(mapping_receipt_path)}) + "\n")
    print(f"[Recovery Completed Verified] {sample}: {mapping_receipt_path}", flush=True)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=BRANCH_ROOT)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--skip-probe", action="store_true")
    parser.add_argument("--probe-only", action="store_true")
    args = parser.parse_args()

    if not args.skip_probe:
        run_namespaced_header_instrain_191_runtime_probe(args.root, threads=min(args.threads, 8))
    if not args.probe_only:
        recover_v60prr2_sample_mapping(args.root, threads=args.threads)
    return 0


if __name__ == "__main__":
    sys.exit(main())
