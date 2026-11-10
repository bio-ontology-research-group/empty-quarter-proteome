#!/usr/bin/env python3
"""Shared, fail-closed contracts for competitive mapping and reference validation.

The pair ledger stores SHA256(raw oriented sequence).  The reference crosswalk
stores SHA256(min(sequence, reverse_complement(sequence))).  These hash domains
are deliberately never joined directly.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any, Callable

import pandas as pd

CROSSWALK_COLUMNS = [
    "original_header",
    "source_category",
    "entity_id",
    "contig_id",
    "seq_length",
    "canonical_sha256",
    "is_reference_representative",
    "representative_header",
    "dedup_reason",
    "ambiguity_status",
]
PAIR_COLUMNS = [
    "host_mag",
    "host_phylum",
    "host_sample",
    "has_host_fastq",
    "target_id",
    "canonical_annotated_target_id",
    "votu_id",
    "virus_sample",
    "has_vir_fastq",
    "target_seq_sha256",
    "individual_virality_status",
    "eligible_for_joint_competitive_mapping",
    "eligible_for_structural_hallmark_followup",
]
VIRAL_CATEGORIES = {"virus_target", "all_viruses", "votu_reps", "coasm_viruses"}
VERIFIED = "verified_direct_viral_target"
UNRESOLVED = "unresolved_virality_pending_boundary_review"
REVCOMP_TABLE = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SHA256_HEX = frozenset("0123456789abcdef")


def stream_sha256(path: Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


def sha256_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode("ascii")).hexdigest()


def reverse_complement(sequence: str) -> str:
    return sequence.translate(REVCOMP_TABLE)[::-1].upper()


def canonical_sha256(sequence: str) -> str:
    sequence = sequence.upper()
    rc = reverse_complement(sequence)
    return sha256_sequence(min(sequence, rc))


def atomic_json(path: Path, payload: Any) -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def atomic_text(path: Path, text: str = "") -> None:
    temporary = path.with_name(path.name + f".tmp.{os.getpid()}")
    temporary.write_text(text)
    os.replace(temporary, path)


def require_columns(frame: pd.DataFrame, required: list[str], label: str) -> None:
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} missing required columns: {missing}; actual={list(frame.columns)}")
    if frame.empty:
        raise ValueError(f"{label} is empty")


def validate_crosswalk_schema(frame: pd.DataFrame) -> None:
    """Require the actual builder schema, including representative_header."""
    require_columns(frame, CROSSWALK_COLUMNS, "reference crosswalk")
    for column in ["original_header", "source_category", "entity_id", "canonical_sha256", "representative_header"]:
        if frame[column].isna().any() or (frame[column].astype(str).str.len() == 0).any():
            raise ValueError(f"reference crosswalk has missing values in {column}")
    lengths = pd.to_numeric(frame["seq_length"], errors="raise")
    if (lengths <= 0).any():
        raise ValueError("reference crosswalk contains non-positive sequence lengths")
    shas = frame["canonical_sha256"].astype(str)
    invalid = shas[(shas.str.len() != 64) | ~shas.apply(lambda value: set(value) <= SHA256_HEX)]
    if not invalid.empty:
        raise ValueError(f"invalid canonical SHA256 values, examples={invalid.head().tolist()}")


def parse_bool(value: Any, label: str) -> bool:
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid Boolean for {label}: {value!r}")


def load_reviewed_pairs(path: Path) -> pd.DataFrame:
    """Load all reviewed pairs and collapse aliases by logical host + raw SHA.

    Verified targets are primary irrespective of legacy sample-eligibility flags;
    sample availability is handled separately.  The unresolved 41Sr3 candidate is
    retained only as sensitivity analysis.  No target count is hard-coded.
    """
    raw = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    require_columns(raw, PAIR_COLUMNS, "reviewed pair ledger")
    raw = raw.apply(lambda column: column.str.strip())
    unknown = sorted(set(raw["individual_virality_status"]) - {VERIFIED, UNRESOLVED})
    if unknown:
        raise ValueError(f"unreviewed individual_virality_status values: {unknown}")
    for index, row in raw.iterrows():
        sha = row["target_seq_sha256"]
        if len(sha) != 64 or not set(sha) <= SHA256_HEX:
            raise ValueError(f"invalid raw oriented target SHA at ledger row {index}: {sha!r}")
        for column in ["host_mag", "target_id", "canonical_annotated_target_id", "host_sample"]:
            if not row[column]:
                raise ValueError(f"missing {column} at ledger row {index}")
        for column in ["has_host_fastq", "has_vir_fastq", "eligible_for_joint_competitive_mapping", "eligible_for_structural_hallmark_followup"]:
            parse_bool(row[column], f"ledger row {index} {column}")

    records: list[dict[str, Any]] = []
    for (host, target_sha), group in raw.groupby(["host_mag", "target_seq_sha256"], sort=True):
        statuses = set(group["individual_virality_status"])
        if len(statuses) != 1:
            raise ValueError(f"conflicting review statuses for {host} / {target_sha}")
        for column in ["host_phylum", "host_sample", "canonical_annotated_target_id"]:
            if set(group[column]) != {group.iloc[0][column]}:
                raise ValueError(f"conflicting {column} within host+SHA aliases for {host} / {target_sha}")
        aliases = sorted(set(group["target_id"]))
        canonical_target = group.iloc[0]["canonical_annotated_target_id"]
        representative = canonical_target if canonical_target in aliases else aliases[0]
        status = group.iloc[0]["individual_virality_status"]
        records.append({
            "host_mag": host,
            "host_phylum": group.iloc[0]["host_phylum"],
            "host_sample": group.iloc[0]["host_sample"],
            "target_id": representative,
            "canonical_annotated_target_id": canonical_target,
            "target_aliases": aliases,
            "target_seq_sha256": target_sha,
            "individual_virality_status": status,
            "analysis_role": "primary" if status == VERIFIED else "sensitivity_unresolved_virality",
            "votu_ids": sorted(set(group["votu_id"])),
            "virus_samples": sorted(set(group["virus_sample"])),
            "has_host_fastq": any(parse_bool(v, "has_host_fastq") for v in group["has_host_fastq"]),
            "has_vir_fastq": any(parse_bool(v, "has_vir_fastq") for v in group["has_vir_fastq"]),
            "eligible_for_structural_hallmark_followup": any(
                parse_bool(v, "eligible_for_structural_hallmark_followup")
                for v in group["eligible_for_structural_hallmark_followup"]
            ),
            "source_pair_rows": len(group),
        })
    result = pd.DataFrame(records)
    if result.empty:
        raise ValueError("reviewed pair ledger produced no logical host+SHA pairs")
    return result


def load_boundary_candidates(path: Path) -> pd.DataFrame:
    """Load root-reviewed boundary candidates without promoting them to primary."""
    payload = json.loads(path.read_text())
    if payload.get("schema_version") != 1 or not payload.get("candidates"):
        raise ValueError("boundary candidate configuration is empty or has an unsupported schema")
    records: list[dict[str, Any]] = []
    for candidate in payload["candidates"]:
        required = {
            "host_mag", "host_sample", "target_id", "target_seq_sha256", "target_length",
            "virus_samples", "disposition", "catalog_viral_interval",
            "fresh_checkv_viral_interval", "fresh_checkv_host_interval",
            "nested_competitor_ids", "intervals",
        }
        if not required <= set(candidate):
            raise ValueError(f"boundary candidate missing fields: {sorted(required-set(candidate))}")
        if candidate["disposition"] != "boundary_sensitive_provisional_not_primary_host_link":
            raise ValueError(f"unapproved boundary disposition: {candidate['disposition']}")
        target_sha = candidate["target_seq_sha256"]
        if len(target_sha) != 64 or not set(target_sha) <= SHA256_HEX:
            raise ValueError("invalid boundary raw target SHA")
        if len(candidate["intervals"]) != 3:
            raise ValueError("boundary follow-up requires core, disputed extension and host flank")
        intervals = sorted(candidate["intervals"], key=lambda row: row["start"])
        position = 1
        for interval in intervals:
            if interval["start"] != position or interval["end"] < interval["start"]:
                raise ValueError(f"boundary intervals are not contiguous 1-based partitions: {intervals}")
            position = interval["end"] + 1
        if position - 1 != candidate["target_length"]:
            raise ValueError("boundary intervals do not sum to full target length")
        records.append({
            "host_mag": candidate["host_mag"],
            "host_phylum": candidate.get("host_phylum", "p__Thermoplasmatota"),
            "host_sample": candidate["host_sample"],
            "target_id": candidate["target_id"],
            "canonical_annotated_target_id": candidate["target_id"],
            "target_aliases": [candidate["target_id"]],
            "target_seq_sha256": target_sha,
            "target_length": int(candidate["target_length"]),
            "individual_virality_status": candidate["disposition"],
            "analysis_role": "boundary_sensitive_provisional",
            "votu_ids": [candidate["target_id"]],
            "virus_samples": list(candidate["virus_samples"]),
            "has_host_fastq": True,
            "has_vir_fastq": True,
            "eligible_for_structural_hallmark_followup": False,
            "source_pair_rows": 0,
            "boundary_intervals": intervals,
            "nested_competitor_ids": list(candidate["nested_competitor_ids"]),
            "catalog_viral_interval": list(candidate["catalog_viral_interval"]),
            "fresh_checkv_viral_interval": list(candidate["fresh_checkv_viral_interval"]),
            "fresh_checkv_host_interval": list(candidate["fresh_checkv_host_interval"]),
            "boundary_config_sha256": stream_sha256(path),
        })
    result = pd.DataFrame(records)
    if result["target_seq_sha256"].duplicated().any():
        raise ValueError("duplicate boundary target SHA")
    return result


def load_analysis_pairs(pair_ledger: Path, boundary_config: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reviewed = load_reviewed_pairs(pair_ledger)
    boundary = load_boundary_candidates(boundary_config)
    overlap = set(reviewed["target_seq_sha256"]) & set(boundary["target_seq_sha256"])
    if overlap:
        raise ValueError(f"boundary target duplicated in reviewed primary/sensitivity ledger: {sorted(overlap)}")
    combined = pd.concat([reviewed, boundary], ignore_index=True, sort=False)
    return combined, reviewed, boundary


def resolve_focal_targets(
    pairs: pd.DataFrame,
    crosswalk: pd.DataFrame,
    fetch_reference_sequence: Callable[[str], str],
) -> dict[str, dict[str, Any]]:
    """Resolve target IDs through original aliases, then verify oriented SHA.

    Alias lookup uses source IDs/headers.  Only after retrieving the representative
    sequence do we compare SHA256(reference) and SHA256(reverse_complement) to the
    ledger's raw-oriented SHA.  The canonical hash is independently recomputed.
    """
    validate_crosswalk_schema(crosswalk)
    resolved: dict[str, dict[str, Any]] = {}
    for target_sha, pair_group in pairs.groupby("target_seq_sha256", sort=True):
        aliases: set[str] = set()
        def collect_aliases(value: Any) -> None:
            if isinstance(value, str):
                aliases.add(value)
            else:
                for nested in value:
                    collect_aliases(nested)
        for values in pair_group["target_aliases"]:
            collect_aliases(values)
        aliases = sorted(aliases)
        alias_rows = crosswalk[
            crosswalk["source_category"].isin(VIRAL_CATEGORIES)
            & crosswalk["entity_id"].isin(aliases)
        ]
        missing_aliases = sorted(set(aliases) - set(alias_rows["entity_id"]))
        if missing_aliases:
            raise ValueError(f"focal aliases absent from original-source crosswalk: {missing_aliases}")
        canonical_values = set(alias_rows["canonical_sha256"])
        representatives = set(alias_rows["representative_header"])
        if len(canonical_values) != 1 or len(representatives) != 1:
            raise ValueError(
                f"inconsistent reference alias resolution for {aliases}: "
                f"canonical={sorted(canonical_values)}, representatives={sorted(representatives)}"
            )
        canonical = next(iter(canonical_values))
        representative = next(iter(representatives))
        lengths = set(pd.to_numeric(alias_rows["seq_length"], errors="raise").astype(int))
        if len(lengths) != 1:
            raise ValueError(f"inconsistent focal sequence lengths for {aliases}: {lengths}")
        sequence = fetch_reference_sequence(representative).upper()
        if not sequence:
            raise ValueError(f"empty reference sequence for focal alias {representative}")
        if len(sequence) != next(iter(lengths)):
            raise ValueError(f"reference length mismatch for {representative}")
        observed_canonical = canonical_sha256(sequence)
        if observed_canonical != canonical:
            raise ValueError(f"canonical hash mismatch for {representative}")
        direct_sha = sha256_sequence(sequence)
        rc_sha = sha256_sequence(reverse_complement(sequence))
        if target_sha == direct_sha == rc_sha:
            orientation = "palindromic"
        elif target_sha == direct_sha:
            orientation = "same_as_reference"
        elif target_sha == rc_sha:
            orientation = "reverse_complement_of_reference"
        else:
            raise ValueError(
                f"raw oriented target SHA mismatch for {aliases}: ledger={target_sha}, "
                f"reference={direct_sha}, reverse_complement={rc_sha}"
            )
        resolved[target_sha] = {
            "target_aliases": aliases,
            "canonical_annotated_target_ids": sorted(set(pair_group["canonical_annotated_target_id"])),
            "canonical_sha256": canonical,
            "reference_contig": representative,
            "reference_orientation_to_ledger": orientation,
            "length_bp": len(sequence),
            "reference_raw_oriented_sha256": direct_sha,
            "ledger_raw_oriented_sha256": target_sha,
        }
    return resolved


def reference_owners(crosswalk: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """Derive ownership from all rows; never trust the builder's category label."""
    owners: dict[str, dict[str, Any]] = {}
    grouped = crosswalk.groupby("representative_header", sort=False)
    for representative, group in grouped:
        canonical_values = set(group["canonical_sha256"])
        lengths = set(pd.to_numeric(group["seq_length"], errors="raise").astype(int))
        if len(canonical_values) != 1 or len(lengths) != 1:
            raise ValueError(f"inconsistent crosswalk representative {representative}")
        owners[representative] = {
            "categories": set(group["source_category"]),
            "entities": set(group["entity_id"]),
            "canonical_sha256": next(iter(canonical_values)),
            "length_bp": next(iter(lengths)),
        }
    return owners


def fasta_sequence_by_id(path: Path, wanted: set[str]) -> dict[str, str]:
    """Small-fixture/source helper; rejects conflicting or missing exact IDs."""
    found: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if current is None or current not in wanted:
            return
        sequence = "".join(chunks).upper()
        if current in found and found[current] != sequence:
            raise ValueError(f"conflicting duplicate FASTA ID {current} in {path}")
        found[current] = sequence

    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                finish()
                current = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    finish()
    missing = wanted - set(found)
    if missing:
        raise ValueError(f"FASTA IDs missing from {path}: {sorted(missing)}")
    return found
