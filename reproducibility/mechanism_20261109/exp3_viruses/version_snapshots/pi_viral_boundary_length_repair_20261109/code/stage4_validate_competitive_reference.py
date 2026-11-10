#!/usr/bin/env python3
"""Independent post-build validation of the immutable competitive reference.

This consumer validates the running builder's actual outputs; it never rebuilds
or edits them.  All bulk hashing, archive streaming and FASTA scanning are meant
for the wrapping Slurm job.  The validator fails closed on missing sources,
source/crosswalk disagreement, non-unique reference headers, stale indexes, or
any unresolved reviewed focal alias.
"""
from __future__ import annotations

import collections
import hashlib
import io
import json
import os
import platform
import subprocess
import sys
import tarfile
import time
from pathlib import Path
from typing import Any, Iterator

import pandas as pd

from stage4_competitive_common import (
    VIRAL_CATEGORIES,
    atomic_json,
    atomic_text,
    canonical_sha256,
    load_analysis_pairs,
    reference_owners,
    require_columns,
    resolve_focal_targets,
    reverse_complement,
    stream_sha256,
    validate_crosswalk_schema,
)

BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
RES_DIR = BRANCH_ROOT / "results"
REF_DIR = RES_DIR / "stage4_competitive_reference"
CODE_DIR = BRANCH_ROOT / "code"
P_REFINED_TAR = Path("/project/k10018/empty-quarter/data/mobilome/tars/refined_bins_all.tar")
SOURCE_SPECS = [
    ("virus_target", "VIRUS_TARGET", RES_DIR / "stage2_raw_outputs_v2/all_matched_primary_targets.fna"),
    ("votu_reps", "VIRUS_VOTU", Path("/project/k10018/empty-quarter/data/mobilome/catalogues/votu_vclust_1.3.1_ani95_af85/votu_reps.fna")),
    ("all_viruses", "VIRUS_ALL", Path("/project/k10018/empty-quarter/data/mobilome/mobilome_inputs/all_viruses.fna")),
    ("coasm_viruses", "VIRUS_COASM", Path("/scratch/project/k10018/work/mobilome_hyp_20261104/e4/derep40/coasm_viruses.fna")),
    ("all_plasmids", "PLASMID", Path("/project/k10018/empty-quarter/data/mobilome/mobilome_inputs/all_plasmids.fna")),
    ("cellular_controls", "DECOY", RES_DIR / "stage2_raw_outputs_v2/curated_cellular_controls_verified.fna"),
]
PREFIX_BY_CATEGORY = {category: prefix for category, prefix, _ in SOURCE_SPECS}
PREFIX_BY_CATEGORY["host_mag"] = "HOST"


def run_command(cmd: list[str]) -> str:
    print(f"[CMD] {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {cmd}: {result.stderr}")
    return result.stdout


def parse_fasta(handle: io.TextIOBase) -> Iterator[tuple[str, str]]:
    header: str | None = None
    chunks: list[str] = []
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks).upper()
            header = line[1:].strip().split()[0]
            if not header:
                raise ValueError("empty FASTA header")
            chunks = []
        else:
            if header is None:
                if line.strip():
                    raise ValueError("sequence data before first FASTA header")
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks).upper()


def validate_crosswalk_internal(crosswalk: pd.DataFrame, manifest: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Validate SHA/header/ownership contracts using representative_header."""
    validate_crosswalk_schema(crosswalk)
    require_columns(manifest, ["mag", "sha256", "contigs", "bases"], "canonical genome manifest")
    if manifest["mag"].duplicated().any() or len(manifest) != 5438:
        raise ValueError(f"canonical manifest must contain exactly 5438 unique MAGs, found={len(manifest)}")
    flags = pd.to_numeric(crosswalk["is_reference_representative"], errors="raise")
    if not flags.isin([0, 1]).all():
        raise ValueError("is_reference_representative must be 0/1")
    representative_rows = crosswalk[flags == 1]
    if (representative_rows["original_header"] != representative_rows["representative_header"]).any():
        raise ValueError("representative row original_header differs from representative_header")
    rep_flag_counts = representative_rows.groupby("representative_header").size()
    if set(rep_flag_counts.index) != set(crosswalk["representative_header"]):
        raise ValueError("one or more crosswalk representatives have no representative row")
    if (rep_flag_counts != 1).any():
        raise ValueError("a representative header has zero/multiple representative rows")
    nonrepresentative = crosswalk[flags == 0]
    if (nonrepresentative["original_header"] == nonrepresentative["representative_header"]).any():
        raise ValueError("non-representative alias has original_header equal to representative_header")

    representatives: dict[str, dict[str, Any]] = {}
    for canonical, group in crosswalk.groupby("canonical_sha256", sort=False):
        reps = set(group["representative_header"])
        lengths = set(group["seq_length"].astype(int))
        if len(reps) != 1 or len(lengths) != 1:
            raise ValueError(f"canonical SHA has inconsistent representative/length: {canonical}")
        representative = next(iter(reps))
        previous = representatives.get(representative)
        record = {"canonical_sha256": canonical, "length_bp": next(iter(lengths))}
        if previous is not None and previous != record:
            raise ValueError(f"representative header maps to multiple canonical SHAs: {representative}")
        representatives[representative] = record

    for category, prefix in PREFIX_BY_CATEGORY.items():
        subset = crosswalk[crosswalk["source_category"] == category]
        if subset.empty:
            raise ValueError(f"required source category absent from crosswalk: {category}")
        if category == "host_mag":
            expected = "HOST__" + subset["entity_id"].astype(str) + "__" + subset["contig_id"].astype(str)
        else:
            expected = prefix + "__" + subset["entity_id"].astype(str)
            if (subset["entity_id"].astype(str) != subset["contig_id"].astype(str)).any():
                raise ValueError(f"non-host entity_id/contig_id mismatch in {category}")
        if (subset["original_header"].astype(str).reset_index(drop=True) != expected.reset_index(drop=True)).any():
            raise ValueError(f"original header/source ownership mismatch in {category}")

    hosts = crosswalk[crosswalk["source_category"] == "host_mag"]
    host_counts = hosts.groupby("entity_id").agg(contigs=("contig_id", "size"), bases=("seq_length", "sum"))
    expected = manifest.set_index("mag")
    if set(host_counts.index) != set(expected.index):
        raise ValueError(
            f"crosswalk canonical MAG IDs differ from manifest; missing={sorted(set(expected.index)-set(host_counts.index))[:10]}, "
            f"extra={sorted(set(host_counts.index)-set(expected.index))[:10]}"
        )
    joined = expected[["contigs", "bases"]].join(host_counts, lsuffix="_manifest", rsuffix="_crosswalk")
    bad = joined[
        (joined["contigs_manifest"].astype(int) != joined["contigs_crosswalk"].astype(int))
        | (joined["bases_manifest"].astype(int) != joined["bases_crosswalk"].astype(int))
    ]
    if not bad.empty:
        raise ValueError(f"canonical MAG contig/base totals differ from manifest: {bad.head().to_dict(orient='index')}")
    return representatives


def expected_source_counter(crosswalk: pd.DataFrame, category: str) -> collections.Counter[tuple[str, str]]:
    subset = crosswalk[crosswalk["source_category"] == category]
    return collections.Counter(zip(subset["original_header"].astype(str), subset["canonical_sha256"].astype(str)))


def audit_source_record(
    observed: collections.Counter[tuple[str, str]],
    expected: collections.Counter[tuple[str, str]],
    original_header: str,
    sequence: str,
    source: Path,
) -> None:
    if not sequence:
        raise ValueError(f"empty sequence {original_header} in {source}")
    key = (original_header, canonical_sha256(sequence))
    if key not in expected:
        raise ValueError(f"source record absent/mismatched in crosswalk: {source}: {original_header}")
    observed[key] += 1
    if observed[key] > expected[key]:
        raise ValueError(f"source contains more copies than crosswalk: {source}: {original_header}")


def parse_fai(path: Path) -> dict[str, int]:
    references: dict[str, int] = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"malformed FAI row: {line!r}")
            name, length = fields[0], int(fields[1])
            if not name or name in references or length <= 0:
                raise ValueError(f"duplicate/invalid FAI reference: {line!r}")
            references[name] = length
    if not references:
        raise ValueError("empty FAI")
    return references


def validate_reference_fasta(
    fasta: Path,
    representatives: dict[str, dict[str, Any]],
    focal_representatives: set[str],
) -> tuple[int, dict[str, str]]:
    observed: set[str] = set()
    focal_sequences: dict[str, str] = {}
    with open(fasta, errors="strict") as handle:
        for header, sequence in parse_fasta(handle):
            if any(character.isspace() for character in header):
                raise ValueError(f"BAM-incompatible whitespace in reference token: {header!r}")
            if header in observed:
                raise ValueError(f"duplicate reference FASTA/BAM header: {header}")
            observed.add(header)
            if header not in representatives:
                raise ValueError(f"reference header absent from crosswalk: {header}")
            expected = representatives[header]
            if len(sequence) != expected["length_bp"] or canonical_sha256(sequence) != expected["canonical_sha256"]:
                raise ValueError(f"reference sequence/SHA/length mismatch for {header}")
            if header in focal_representatives:
                focal_sequences[header] = sequence
    if observed != set(representatives):
        raise ValueError(
            f"reference FASTA headers differ from crosswalk: missing={sorted(set(representatives)-observed)[:10]}, "
            f"extra={sorted(observed-set(representatives))[:10]}"
        )
    if focal_representatives - set(focal_sequences):
        raise ValueError("one or more focal representative sequences were not recovered")
    return len(observed), focal_sequences


def verify_build_inputs_and_sentinel(
    build_done: Path,
    nonempty_paths: list[Path],
    build_receipt_path: Path,
) -> dict[str, Any]:
    """Require touch-created REFERENCE_BUILD_DONE sentinel to exist as a file (0 bytes allowed),
    while requiring all data/receipt/archive inputs to be nonempty files and build_receipt status=='completed'.
    """
    if not build_done.exists() or not build_done.is_file():
        raise FileNotFoundError(f"required build sentinel missing or not a file: {build_done}")
    for path in nonempty_paths:
        if not path.exists() or not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"required build/source input missing or empty: {path}")
    build_receipt = json.loads(build_receipt_path.read_text())
    if build_receipt.get("status") != "completed":
        raise ValueError("builder receipt is not completed")
    return build_receipt


def validate_boundary_candidates_and_nested_competitors(
    boundary_pairs: pd.DataFrame,
    resolved: dict[str, dict[str, Any]],
    focal_sequences: dict[str, str],
    crosswalk: pd.DataFrame,
) -> dict[str, Any]:
    nested_validation: dict[str, Any] = {}
    for candidate in boundary_pairs.to_dict(orient="records"):
        full_info = resolved[candidate["target_seq_sha256"]]
        if full_info["length_bp"] != int(candidate["target_length"]):
            raise ValueError(f"boundary target length differs from root disposition: {candidate['target_id']}")
        full_sequence = focal_sequences[full_info["reference_contig"]]
        if full_info["reference_orientation_to_ledger"] == "reverse_complement_of_reference":
            full_sequence = reverse_complement(full_sequence)
        for nested_id in candidate["nested_competitor_ids"]:
            alias_rows = crosswalk[crosswalk["entity_id"] == nested_id]
            if alias_rows.empty:
                raise ValueError(f"nested competitor missing from original-source crosswalk: {nested_id}")
            nested_reps = set(alias_rows["representative_header"])
            if len(nested_reps) != 1:
                raise ValueError(f"nested competitor resolves to multiple references: {nested_id}")
            nested_rep = next(iter(nested_reps))
            nested_sequence = focal_sequences[nested_rep]
            core_start, core_end = candidate["catalog_viral_interval"]
            expected_core = full_sequence[core_start - 1:core_end]
            if expected_core != nested_sequence and expected_core != reverse_complement(nested_sequence):
                raise ValueError(f"nested competitor is not the exact catalog core of the full target: {nested_id}")
            nested_validation[nested_id] = {
                "reference_contig": nested_rep,
                "length_bp": len(nested_sequence),
                "exact_catalog_core_match_forward_or_reverse_complement": True,
            }
    return nested_validation


def main() -> int:
    start = time.time()
    done = REF_DIR / "REFERENCE_VALIDATION_DONE"
    done.unlink(missing_ok=True)
    build_done = REF_DIR / "REFERENCE_BUILD_DONE"
    ref_fna = REF_DIR / "unified_competitive_reference.fna"
    crosswalk_path = REF_DIR / "stage4_reference_sequence_crosswalk.tsv"
    build_receipt_path = REF_DIR / "unified_competitive_reference.receipt.json"
    build_receipt = verify_build_inputs_and_sentinel(
        build_done=build_done,
        nonempty_paths=[ref_fna, crosswalk_path, build_receipt_path, P_REFINED_TAR],
        build_receipt_path=build_receipt_path,
    )

    manifest_path = BRANCH_ROOT / "inputs" / "genome_manifest.tsv"
    pair_path = RES_DIR / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    manifest = pd.read_csv(manifest_path, sep="\t")
    boundary_path = BRANCH_ROOT / "inputs/stage4_boundary_sensitive_intervals.json"
    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_path, boundary_path)
    print("[Crosswalk] Loading actual builder schema", flush=True)
    crosswalk = pd.read_csv(crosswalk_path, sep="\t")
    representatives = validate_crosswalk_internal(crosswalk, manifest)
    owners = reference_owners(crosswalk)
    if set(owners) != set(representatives):
        raise ValueError("ownership and representative indexes differ")

    focal_aliases = {alias for aliases in pairs["target_aliases"] for alias in aliases}
    nested_competitor_ids = {nested for values in boundary_pairs["nested_competitor_ids"] for nested in values}
    focal_aliases |= nested_competitor_ids
    focal_rows = crosswalk[crosswalk["source_category"].isin(VIRAL_CATEGORIES) & crosswalk["entity_id"].isin(focal_aliases)]
    focal_representatives = set(focal_rows["representative_header"])

    ref_fai = Path(str(ref_fna) + ".fai")
    if not ref_fai.exists() or ref_fai.stat().st_size == 0:
        run_command(["samtools", "faidx", str(ref_fna)])
    fai = parse_fai(ref_fai)
    if fai != {name: record["length_bp"] for name, record in representatives.items()}:
        raise ValueError("FAI reference names/lengths do not exactly match crosswalk")
    print(f"[Reference] Hashing and scanning {ref_fna}", flush=True)
    ref_sha = stream_sha256(ref_fna)
    reference_count, focal_sequences = validate_reference_fasta(ref_fna, representatives, focal_representatives)
    resolved = resolve_focal_targets(pairs, crosswalk, lambda name: focal_sequences[name])
    nested_validation = validate_boundary_candidates_and_nested_competitors(
        boundary_pairs=boundary_pairs,
        resolved=resolved,
        focal_sequences=focal_sequences,
        crosswalk=crosswalk,
    )

    index_files: list[dict[str, Any]] = []
    for suffix in [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"]:
        candidates = [REF_DIR / f"unified_competitive_reference{suffix}.bt2l", REF_DIR / f"unified_competitive_reference{suffix}.bt2"]
        existing = [path for path in candidates if path.exists() and path.stat().st_size > 0]
        if len(existing) != 1:
            raise ValueError(f"expected exactly one nonempty Bowtie2 index for {suffix}: {candidates}")
        path = existing[0]
        index_files.append({"name": path.name, "path": str(path), "bytes": path.stat().st_size, "sha256": stream_sha256(path)})
    builder_indices = build_receipt.get("index_files", {})
    for index in index_files:
        if index["name"] not in builder_indices or builder_indices[index["name"]].get("sha256") != index["sha256"]:
            raise ValueError(f"index differs from immutable builder receipt: {index['name']}")

    source_receipts: dict[str, Any] = {}
    expected_mag_bytes = manifest.set_index("mag")["sha256"].astype(str).to_dict()
    expected_mag_counts = manifest.set_index("mag")[["contigs", "bases"]].astype(int).to_dict(orient="index")
    host_expected = expected_source_counter(crosswalk, "host_mag")
    host_observed: collections.Counter[tuple[str, str]] = collections.Counter()
    found_mags: set[str] = set()
    mag_hashes: dict[str, str] = {}
    print(f"[Sources] Hashing and streaming canonical archive {P_REFINED_TAR}", flush=True)
    tar_sha = stream_sha256(P_REFINED_TAR)
    with tarfile.open(P_REFINED_TAR, "r:") as archive:
        for member in archive:
            if not member.isfile() or not member.name.endswith(".fa"):
                continue
            mag = Path(member.name).name[:-3]
            if mag not in expected_mag_bytes:
                continue
            if mag in found_mags:
                raise ValueError(f"duplicate canonical MAG member in archive: {mag}")
            extracted = archive.extractfile(member)
            if extracted is None:
                raise ValueError(f"cannot extract canonical MAG: {mag}")
            raw = extracted.read()
            raw_sha = hashlib.sha256(raw).hexdigest()
            if raw_sha != expected_mag_bytes[mag]:
                raise ValueError(f"canonical FASTA byte SHA mismatch for {mag}")
            records = list(parse_fasta(io.StringIO(raw.decode("utf-8"))))
            if len(records) != expected_mag_counts[mag]["contigs"] or sum(len(seq) for _, seq in records) != expected_mag_counts[mag]["bases"]:
                raise ValueError(f"source FASTA contig/base mismatch for {mag}")
            for contig, sequence in records:
                audit_source_record(host_observed, host_expected, f"HOST__{mag}__{contig}", sequence, P_REFINED_TAR)
            found_mags.add(mag)
            mag_hashes[mag] = raw_sha
            if len(found_mags) % 500 == 0:
                print(f"  verified {len(found_mags)}/5438 canonical MAGs", flush=True)
    if found_mags != set(expected_mag_bytes) or host_observed != host_expected:
        raise ValueError(f"canonical archive/crosswalk source coverage mismatch; found {len(found_mags)}/5438 MAGs")
    source_receipts["canonical_mag_archive"] = {
        "path": str(P_REFINED_TAR), "sha256": tar_sha, "bytes": P_REFINED_TAR.stat().st_size,
        "canonical_mag_count": len(found_mags), "canonical_contig_count": int(manifest["contigs"].sum()),
        "canonical_base_count": int(manifest["bases"].sum()), "per_mag_fasta_byte_sha256": mag_hashes,
    }

    for category, prefix, path in SOURCE_SPECS:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"required builder source missing/empty: {category}: {path}")
        print(f"[Sources] Validating {category}: {path}", flush=True)
        expected = expected_source_counter(crosswalk, category)
        observed: collections.Counter[tuple[str, str]] = collections.Counter()
        count = 0
        bases = 0
        with open(path, errors="strict") as handle:
            for raw_header, sequence in parse_fasta(handle):
                audit_source_record(observed, expected, f"{prefix}__{raw_header}", sequence, path)
                count += 1
                bases += len(sequence)
        if observed != expected or count != sum(expected.values()):
            missing = expected - observed
            extra = observed - expected
            raise ValueError(f"source/crosswalk count mismatch for {category}: missing={missing.most_common(3)}, extra={extra.most_common(3)}")
        source_receipts[category] = {
            "path": str(path), "sha256": stream_sha256(path), "bytes": path.stat().st_size,
            "record_count": count, "base_count": bases, "crosswalk_records_verified": sum(expected.values()),
        }

    ownership_counts: dict[str, int] = collections.Counter()
    for owner in owners.values():
        categories = owner["categories"]
        if "host_mag" in categories and len(categories) > 1:
            ownership_counts["host_shared_with_nonhost"] += 1
        elif categories == {"host_mag"} and len(owner["entities"]) > 1:
            ownership_counts["host_host_shared"] += 1
        elif categories == {"host_mag"}:
            ownership_counts["single_host_only"] += 1
        if categories <= VIRAL_CATEGORIES:
            ownership_counts["viral_only"] += 1
        if "all_plasmids" in categories and categories & VIRAL_CATEGORIES:
            ownership_counts["viral_plasmid_shared"] += 1
        if "cellular_controls" in categories and categories & VIRAL_CATEGORIES:
            ownership_counts["viral_decoy_shared"] += 1

    builder_script = CODE_DIR / "stage4_build_combined_competitive_reference.py"
    receipt = {
        "status": "completed_verified",
        "schema_version": 3,
        "stage": "stage4_validate_competitive_reference",
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": time.time() - start,
        "execution": {
            "hostname": platform.node(), "slurm_job_id": os.environ.get("SLURM_JOB_ID", "local"),
            "samtools_version": run_command(["samtools", "--version"]).splitlines()[0],
            "validator_script_sha256": stream_sha256(Path(__file__)),
            "common_script_sha256": stream_sha256(Path(__file__).with_name("stage4_competitive_common.py")),
        },
        "builder_script_observed": {
            "path": str(builder_script), "sha256": stream_sha256(builder_script),
            "timing": "observed_at_postbuild_validation_not_claimed_as_launch_hash",
        },
        "build_receipt": {"path": str(build_receipt_path), "sha256": stream_sha256(build_receipt_path)},
        "reference_files": {
            "reference_fna": {"path": str(ref_fna), "bytes": ref_fna.stat().st_size, "sha256": ref_sha},
            "reference_fai": {"path": str(ref_fai), "bytes": ref_fai.stat().st_size, "sha256": stream_sha256(ref_fai)},
            "crosswalk_tsv": {"path": str(crosswalk_path), "bytes": crosswalk_path.stat().st_size, "sha256": stream_sha256(crosswalk_path)},
            "contig_count": reference_count,
            "crosswalk_record_count": len(crosswalk),
        },
        "bowtie2_index": {"prefix": str(REF_DIR / "unified_competitive_reference"), "files": index_files},
        "input_sources": source_receipts,
        "pair_ledger": {"path": str(pair_path), "sha256": stream_sha256(pair_path)},
        "boundary_config": {"path": str(boundary_path), "sha256": stream_sha256(boundary_path)},
        "canonical_manifest": {"path": str(manifest_path), "sha256": stream_sha256(manifest_path)},
        "focal_alias_validation": {
            "n_logical_pairs": len(pairs),
            "n_primary_pairs": int((reviewed_pairs["analysis_role"] == "primary").sum()),
            "n_sensitivity_pairs": int((reviewed_pairs["analysis_role"] != "primary").sum()),
            "n_boundary_sensitive_provisional_pairs": len(boundary_pairs),
            "nested_boundary_competitors": nested_validation,
            "n_unique_raw_target_sha": len(resolved),
            "all_original_aliases_resolved": True,
            "raw_oriented_and_canonical_hash_domains_verified_separately": True,
            "targets": resolved,
        },
        "ownership_audit": dict(ownership_counts),
        "checks": {
            "all_required_sources_present_hashed_and_count_verified": True,
            "exact_5438_mag_ids_contigs_bases_and_fasta_byte_hashes": True,
            "unique_reference_bam_headers": True,
            "reference_fasta_fai_crosswalk_sha_header_consistency": True,
            "six_actual_index_files_hashed_and_builder_receipt_matched": True,
            "full_focal_target_alias_coverage": True,
        },
    }
    receipt_path = REF_DIR / "stage4_reference_validation_receipt.json"
    atomic_json(receipt_path, receipt)
    atomic_text(done, json.dumps({"receipt_sha256": stream_sha256(receipt_path)}) + "\n")
    print(f"[REFERENCE VALIDATION COMPLETED VERIFIED] {receipt_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
