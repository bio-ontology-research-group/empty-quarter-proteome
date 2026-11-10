#!/usr/bin/env python3
"""Map one sample to the validated unified competitive reference.

Fail-closed implementation contract:
* consume the builder's actual ``representative_header`` crosswalk column;
* resolve ledger target IDs through original aliases and verify raw oriented SHA
  separately from canonical (strand-insensitive) SHA;
* use actual samtools coverage tabular output (never ``-m`` histogram output),
  and require one valid row for every BAM reference, including measured zeros;
* retain full-MAG denominators while deriving sequence ownership from all
  crosswalk rows and counting each representative's reads at most once per MAG;
* conditionally run inStrain once per unique viral SHA using MAPQ >19 (>=20)
  and source-read ANI >0.95; relaxed >0.92 is a separately named sensitivity.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import shlex
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

BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
RES_DIR = BRANCH_ROOT / "results"
FASTQ_DIR = Path("/project/k10018/empty-quarter/data/reads/fastq")
SIF_INSTRAIN = Path("/project/k10018/sw/images/eq-annot/instrain.sif")
COVERAGE_COLUMNS = [
    "rname", "startpos", "endpos", "numreads", "covbases", "coverage",
    "meandepth", "meanbaseq", "meanmapq",
]


class CommandRunner:
    def __init__(self) -> None:
        self.commands: list[dict[str, Any]] = []

    def run(self, cmd: list[str], description: str) -> subprocess.CompletedProcess[str]:
        printable = shlex.join(cmd)
        print(f"[CMD] {description}: {printable}", flush=True)
        result = subprocess.run(cmd, capture_output=True, text=True)
        self.commands.append({"description": description, "command": cmd, "returncode": result.returncode})
        if result.returncode != 0:
            raise RuntimeError(f"{description} failed ({result.returncode}): {result.stderr}")
        return result

    def alignment_pipeline(
        self,
        bowtie_cmd: list[str],
        raw_bam: Path,
        threads: int,
        stderr_dir: Path,
    ) -> None:
        """Stream SAM through checked Popen processes; never create a giant SAM."""
        view_cmd = ["samtools", "view", "-u", "-"]
        temporary = raw_bam.with_name(raw_bam.name + f".tmp.{os.getpid()}")
        sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(temporary), "-"]
        commands = [bowtie_cmd, view_cmd, sort_cmd]
        print("[PIPELINE] " + " | ".join(shlex.join(command) for command in commands), flush=True)
        logs = [open(stderr_dir / f"alignment_{name}.stderr.log", "w") for name in ["bowtie2", "samtools_view", "samtools_sort"]]
        processes: list[subprocess.Popen[bytes]] = []
        try:
            bowtie = subprocess.Popen(bowtie_cmd, stdout=subprocess.PIPE, stderr=logs[0])
            processes.append(bowtie)
            view = subprocess.Popen(view_cmd, stdin=bowtie.stdout, stdout=subprocess.PIPE, stderr=logs[1])
            processes.append(view)
            assert bowtie.stdout is not None
            bowtie.stdout.close()
            sort = subprocess.Popen(sort_cmd, stdin=view.stdout, stderr=logs[2])
            processes.append(sort)
            assert view.stdout is not None
            view.stdout.close()
            returncodes = [process.wait() for process in processes]
            self.commands.append({
                "description": "stream Bowtie2 SAM to sorted BAM",
                "pipeline": commands,
                "returncodes": returncodes,
            })
            if any(returncodes):
                raise RuntimeError(f"alignment pipeline return codes {returncodes}; inspect {stderr_dir}")
            os.replace(temporary, raw_bam)
        except Exception:
            for process in processes:
                if process.poll() is None:
                    process.kill()
            for process in processes:
                process.wait()
            temporary.unlink(missing_ok=True)
            raise
        finally:
            for log in logs:
                log.close()


def tool_version(cmd: list[str]) -> str:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"tool version probe failed: {shlex.join(cmd)}: {result.stderr}")
    lines = (result.stdout or result.stderr).strip().splitlines()
    if not lines:
        raise RuntimeError(f"empty tool version output: {shlex.join(cmd)}")
    return lines[0].strip()


def compute_coverage_metrics(cov_tsv_text: str, expected_lengths: dict[str, int] | None = None) -> pd.DataFrame:
    """Strictly parse actual samtools coverage table, including explicit zero rows."""
    lines = cov_tsv_text.splitlines()
    if not lines or not lines[0].strip():
        raise ValueError("samtools coverage returned empty output")
    header = lines[0].split("\t")
    if header and header[0] == "#rname":
        header[0] = "rname"
    if header != COVERAGE_COLUMNS:
        raise ValueError(f"unexpected samtools coverage header: {header}; expected={COVERAGE_COLUMNS}")
    if len(lines) == 1:
        raise ValueError("samtools coverage returned a header but no reference rows")
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    for line_number, line in enumerate(lines[1:], start=2):
        parts = line.split("\t")
        if len(parts) != 9:
            raise ValueError(f"malformed samtools coverage row {line_number}: {len(parts)} fields")
        try:
            row: dict[str, Any] = {
                "rname": parts[0],
                "startpos": int(parts[1]),
                "endpos": int(parts[2]),
                "numreads": int(parts[3]),
                "covbases": int(parts[4]),
                "coverage": float(parts[5]),
                "meandepth": float(parts[6]),
                "meanbaseq": float(parts[7]),
                "meanmapq": float(parts[8]),
            }
        except ValueError as error:
            raise ValueError(f"non-numeric samtools coverage row {line_number}: {line}") from error
        name = row["rname"]
        if not name or name in seen:
            raise ValueError(f"empty or duplicate coverage reference at row {line_number}: {name!r}")
        seen.add(name)
        numeric = [value for key, value in row.items() if key != "rname"]
        if not all(math.isfinite(float(value)) for value in numeric):
            raise ValueError(f"non-finite coverage value at row {line_number}")
        if row["startpos"] != 1 or row["endpos"] <= 0 or row["numreads"] < 0:
            raise ValueError(f"invalid coordinates/read count at coverage row {line_number}: {row}")
        if not (0 <= row["covbases"] <= row["endpos"] and 0 <= row["coverage"] <= 100 and row["meandepth"] >= 0):
            raise ValueError(f"invalid coverage values at row {line_number}: {row}")
        if expected_lengths is not None:
            if name not in expected_lengths:
                raise ValueError(f"coverage reference absent from BAM header: {name}")
            if row["endpos"] != expected_lengths[name]:
                raise ValueError(f"coverage/BAM length mismatch for {name}")
        rows.append(row)
    if expected_lengths is not None and seen != set(expected_lengths):
        missing = sorted(set(expected_lengths) - seen)
        extra = sorted(seen - set(expected_lengths))
        raise ValueError(
            f"coverage rows do not exactly match BAM references; missing={missing[:10]} "
            f"(n={len(missing)}), extra={extra[:10]} (n={len(extra)}). Missing is not zero."
        )
    return pd.DataFrame(rows, columns=COVERAGE_COLUMNS)


def compute_interval_coverage(depth_text: str, reference: str, start: int, end: int) -> dict[str, Any]:
    """Parse samtools depth -a output and require every interval position."""
    if start < 1 or end < start:
        raise ValueError(f"invalid interval: {reference}:{start}-{end}")
    depths: list[int] = []
    for expected_position, line in enumerate(depth_text.splitlines(), start=start):
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError(f"malformed samtools depth row: {line!r}")
        if fields[0] != reference or int(fields[1]) != expected_position:
            raise ValueError(f"missing/noncontiguous depth at {reference}:{expected_position}")
        depth = int(fields[2])
        if depth < 0:
            raise ValueError(f"negative depth at {reference}:{expected_position}")
        depths.append(depth)
    expected_length = end - start + 1
    if len(depths) != expected_length:
        raise ValueError(f"depth rows missing for {reference}:{start}-{end}; missing is not zero")
    covered = sum(depth > 0 for depth in depths)
    return {
        "length_bp": expected_length,
        "callable_positions_depth_ge1": covered,
        "horizontal_breadth_pct": 100.0 * covered / expected_length,
        "mean_depth": sum(depths) / expected_length,
    }


def bam_reference_lengths(header_text: str) -> dict[str, int]:
    references: dict[str, int] = {}
    for line in header_text.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        if "SN" not in fields or "LN" not in fields:
            raise ValueError(f"malformed @SQ line: {line}")
        name = fields["SN"]
        length = int(fields["LN"])
        if not name or name in references or length <= 0:
            raise ValueError(f"duplicate/invalid BAM reference: {line}")
        references[name] = length
    if not references:
        raise ValueError("BAM header has no @SQ references")
    return references


def evaluate_host_coverage(
    host: str,
    host_rows: pd.DataFrame,
    owners: dict[str, dict[str, Any]],
    coverage: dict[str, dict[str, Any]],
    manifest_contigs: int,
    manifest_bases: int,
    manifest_fasta_sha256: str,
) -> dict[str, Any]:
    """Aggregate all original contigs without double-counting representative reads."""
    if host_rows.empty:
        raise ValueError(f"host {host} absent from crosswalk")
    total_contigs = len(host_rows)
    total_bases = int(host_rows["seq_length"].astype(int).sum())
    if (total_contigs, total_bases) != (manifest_contigs, manifest_bases):
        raise ValueError(
            f"host denominator mismatch for {host}: crosswalk={(total_contigs, total_bases)}, "
            f"manifest={(manifest_contigs, manifest_bases)}"
        )
    total_covbases = 0
    total_depth_weighted = 0.0
    unique_contigs = unique_bases = unique_covbases = 0
    unique_depth_weighted = 0.0
    unique_representatives: set[str] = set()
    ambiguous_representatives: set[str] = set()
    all_representatives: set[str] = set()
    for row in host_rows.to_dict(orient="records"):
        representative = str(row["representative_header"])
        if representative not in coverage:
            raise ValueError(f"missing coverage row for {representative}; unknown is not zero")
        length = int(row["seq_length"])
        entry = coverage[representative]
        depth = float(entry["meandepth"])
        covered = int(entry["covbases"])
        owner = owners[representative]
        is_unique = owner["categories"] == {"host_mag"} and owner["entities"] == {host}
        total_covbases += covered
        total_depth_weighted += depth * length
        all_representatives.add(representative)
        if is_unique:
            unique_representatives.add(representative)
            unique_contigs += 1
            unique_bases += length
            unique_covbases += covered
            unique_depth_weighted += depth * length
        else:
            ambiguous_representatives.add(representative)
    ambiguous_bases = total_bases - unique_bases
    ambiguous_covbases = total_covbases - unique_covbases
    ambiguous_depth_weighted = total_depth_weighted - unique_depth_weighted
    read_sum = lambda names: int(sum(coverage[name]["numreads"] for name in names))
    return {
        "entity_id": host,
        "entity_type": "host_mag",
        "host_mag": host,
        "reference_contig": "multiple_original_contigs_see_crosswalk",
        "length_bp": total_bases,
        "manifest_contigs": manifest_contigs,
        "manifest_bases": manifest_bases,
        "manifest_fasta_sha256": manifest_fasta_sha256,
        "denominator_verified_manifest": True,
        "filtered_reads_mapped": read_sum(all_representatives),
        "horizontal_breadth_pct": 100.0 * total_covbases / total_bases,
        "mean_depth": total_depth_weighted / total_bases,
        "unique_mag_contigs": unique_contigs,
        "unique_mag_bases": unique_bases,
        "unique_mag_covbases": unique_covbases,
        "unique_filtered_reads_mapped": read_sum(unique_representatives),
        "unique_horizontal_breadth_pct": 100.0 * unique_covbases / total_bases,
        "unique_mean_depth": unique_depth_weighted / total_bases,
        "ambiguous_mag_contigs": total_contigs - unique_contigs,
        "ambiguous_mag_bases": ambiguous_bases,
        "ambiguous_mag_covbases": ambiguous_covbases,
        "ambiguous_filtered_reads_mapped": read_sum(ambiguous_representatives),
        "ambiguous_horizontal_breadth_pct": 100.0 * ambiguous_covbases / total_bases,
        "ambiguous_mean_depth": ambiguous_depth_weighted / total_bases,
        "ambiguity_status": "ownership_derived_from_full_crosswalk",
        "uniquely_attributable": ambiguous_bases == 0,
        "qualifies_10x_50pct_instrain": 0,
        "instrain_action": "host_context_only_not_profiled",
    }


def fetch_faidx_sequence(runner: CommandRunner, fasta: Path, reference: str) -> str:
    result = runner.run(["samtools", "faidx", str(fasta), reference], f"extract exact reference {reference}")
    lines = result.stdout.splitlines()
    if not lines or not lines[0].startswith(">") or lines[0][1:].split()[0] != reference:
        raise ValueError(f"samtools faidx returned wrong reference for {reference}")
    return "".join(line.strip() for line in lines[1:] if not line.startswith(">"))


def verify_validation_receipt(ref_dir: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    receipt_path = ref_dir / "stage4_reference_validation_receipt.json"
    if not (ref_dir / "REFERENCE_VALIDATION_DONE").exists() or not receipt_path.exists():
        raise FileNotFoundError("validated reference receipt/DONE missing; mapping is gated on post-build validation")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("status") != "completed_verified" or receipt.get("schema_version") != 3:
        raise ValueError("reference validation receipt is not a completed schema-v3 receipt")
    files = receipt["reference_files"]
    checks = [
        (Path(files["reference_fna"]["path"]), files["reference_fna"]["sha256"]),
        (Path(files["reference_fai"]["path"]), files["reference_fai"]["sha256"]),
        (Path(files["crosswalk_tsv"]["path"]), files["crosswalk_tsv"]["sha256"]),
        (Path(receipt["build_receipt"]["path"]), receipt["build_receipt"]["sha256"]),
        (receipt_path, None),
    ]
    for index in receipt["bowtie2_index"]["files"]:
        checks.append((Path(index["path"]), index["sha256"]))
    for path, expected in checks:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"validated reference input missing/empty: {path}")
        if expected is not None:
            observed = stream_sha256(path)
            if observed != expected:
                raise ValueError(f"validated reference hash changed for {path}: {observed} != {expected}")
    return receipt, {"receipt_path": str(receipt_path), "receipt_sha256": stream_sha256(receipt_path)}


def verify_raw_bam_reuse(raw_bam: Path, receipt_path: Path, read_hashes: dict[str, str], reference_binding: str) -> None:
    if not raw_bam.exists():
        return
    if not receipt_path.exists():
        raise ValueError(f"existing raw BAM {raw_bam} has no cryptographic receipt; refusing blind reuse or overwrite")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("status") != "completed_verified":
        raise ValueError("existing raw BAM receipt is not completed_verified")
    inputs = receipt.get("inputs", {})
    for key, value in read_hashes.items():
        if inputs.get(key) != value:
            raise ValueError(f"existing raw BAM receipt does not match current {key}")
    if receipt.get("reference_binding", {}).get("validation_receipt_sha256") != reference_binding:
        raise ValueError("existing raw BAM was made against a different reference validation receipt")
    if receipt.get("outputs", {}).get("raw_bam_sha256") != stream_sha256(raw_bam):
        raise ValueError("existing raw BAM hash differs from its receipt")


def require_instrain_available(container: Path, singularity_executable: str | None) -> None:
    if not container.exists() or singularity_executable is None:
        raise RuntimeError(f"eligible inStrain targets exist but container/tool unavailable: {container}")


def collect_instrain_callable_positions(output_dir: Path, representative: str) -> dict[str, Any]:
    candidates = sorted(output_dir.glob("**/*scaffold_info.tsv"))
    if len(candidates) != 1:
        raise ValueError(f"expected one inStrain scaffold_info table in {output_dir}, found={candidates}")
    frame = pd.read_csv(candidates[0], sep="\t")
    require_columns(frame, ["scaffold", "length", "breadth"], "inStrain scaffold_info")
    selected = frame[frame["scaffold"].astype(str) == representative]
    if len(selected) != 1:
        raise ValueError(f"inStrain scaffold_info lacks exactly one row for {representative}")
    row = selected.iloc[0]
    length = int(row["length"])
    breadth = float(row["breadth"])
    if length <= 0 or not 0 <= breadth <= 1:
        raise ValueError(f"invalid inStrain length/breadth for {representative}")
    snv_tables = sorted(output_dir.glob("**/*SNVs.tsv"))
    variant_rows = 0
    if snv_tables:
        variant_rows = sum(max(0, sum(1 for _ in open(path)) - 1) for path in snv_tables)
    return {
        "scaffold_info_path": str(candidates[0]),
        "callable_positions": int(round(length * breadth)),
        "callable_position_method": "round(inStrain scaffold_info length * breadth)",
        "instrain_length_bp": length,
        "instrain_breadth_fraction": breadth,
        "variant_table_rows": variant_rows,
        "variant_interpretation": "sequence polymorphism only; variants are not evidence of recombination",
    }


def main() -> int:
    start = time.time()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--fastq-dir", type=Path, default=FASTQ_DIR)
    parser.add_argument("--boundary-config", type=Path, default=BRANCH_ROOT / "inputs/stage4_boundary_sensitive_intervals.json")
    parser.add_argument("--instrain-mode", choices=["source_primary", "crosssample_sensitivity"], default="source_primary")
    args = parser.parse_args()
    if args.threads < 1:
        raise ValueError("threads must be positive")
    sample = args.sample
    out_dir = args.out_dir or RES_DIR / "stage4_competitive_mapping" / sample
    ref_dir = RES_DIR / "stage4_competitive_reference"
    ref_fna = ref_dir / "unified_competitive_reference.fna"
    crosswalk_path = ref_dir / "stage4_reference_sequence_crosswalk.tsv"
    bt2_prefix = ref_dir / "unified_competitive_reference"
    pair_path = RES_DIR / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    manifest_path = BRANCH_ROOT / "inputs" / "genome_manifest.tsv"
    out_dir.mkdir(parents=True, exist_ok=True)
    done_path = out_dir / "MAPPING_DONE"
    done_path.unlink(missing_ok=True)
    runner = CommandRunner()

    for tool in ["bowtie2", "samtools"]:
        if shutil.which(tool) is None:
            raise RuntimeError(f"required mapping tool unavailable: {tool}")
    versions = {
        "bowtie2": tool_version(["bowtie2", "--version"]),
        "samtools": tool_version(["samtools", "--version"]),
    }
    validation_receipt, validation_binding = verify_validation_receipt(ref_dir)
    print("[Preflight] Reference, six index files, crosswalk, FAI and build receipt hashes verified", flush=True)

    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_path, args.boundary_config)
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

    r1 = args.fastq_dir / f"{sample}_R1.fastq.gz"
    r2 = args.fastq_dir / f"{sample}_R2.fastq.gz"
    se = args.fastq_dir / f"{sample}_SE.fastq.gz"
    if not r1.exists() or not r2.exists():
        raise FileNotFoundError(f"paired FASTQ missing for {sample}: {r1} / {r2}")
    read_paths = {"r1": r1, "r2": r2}
    if se.exists():
        read_paths["se"] = se
    read_hashes = {f"{key}_sha256": stream_sha256(path) for key, path in read_paths.items()}
    read_hashes["se_sha256"] = read_hashes.get("se_sha256", "not_present")
    raw_bam = out_dir / f"{sample}_raw.bam"
    filtered_bam = out_dir / f"{sample}_mapq20.bam"
    mapping_receipt_path = out_dir / f"{sample}_mapping_receipt.json"
    verify_raw_bam_reuse(raw_bam, mapping_receipt_path, read_hashes, validation_binding["receipt_sha256"])

    if not raw_bam.exists():
        bowtie_cmd = [
            "bowtie2", "-x", str(bt2_prefix), "-1", str(r1), "-2", str(r2),
            "-p", str(args.threads), "--very-sensitive",
        ]
        if se.exists():
            bowtie_cmd.extend(["-U", str(se)])
        runner.alignment_pipeline(bowtie_cmd, raw_bam, args.threads, out_dir)
        runner.run(["samtools", "index", "-@", str(args.threads), str(raw_bam)], "index raw BAM")
    runner.run(["samtools", "quickcheck", str(raw_bam)], "quickcheck raw BAM")

    temporary_filtered = filtered_bam.with_name(filtered_bam.name + f".tmp.{os.getpid()}")
    runner.run([
        "samtools", "view", "-b", "-q", "20", "-F", "3844", "-@", str(args.threads),
        "-o", str(temporary_filtered), str(raw_bam),
    ], "filter MAPQ>=20 and UNMAP/SECONDARY/QCFAIL/DUP/SUPPLEMENTARY")
    os.replace(temporary_filtered, filtered_bam)
    runner.run(["samtools", "index", "-@", str(args.threads), str(filtered_bam)], "index filtered BAM")
    runner.run(["samtools", "quickcheck", str(filtered_bam)], "quickcheck filtered BAM")

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

    coverage_path = out_dir / f"{sample}_filtered_coverage.tsv"
    coverage_result = runner.run(["samtools", "coverage", str(filtered_bam)], "compute actual tabular coverage")
    coverage_frame = compute_coverage_metrics(coverage_result.stdout, bam_lengths)
    atomic_text(coverage_path, coverage_result.stdout)
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
    # Root-reviewed boundary candidates get independent contiguous interval evidence.
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
    evaluation_frame = pd.DataFrame(evaluations)
    evaluation_path = out_dir / f"{sample}_focal_targets_evaluation.tsv"
    evaluation_frame.to_csv(evaluation_path, sep="\t", index=False)

    min_read_ani = 0.95 if args.instrain_mode == "source_primary" else 0.92
    instrain_results: list[dict[str, Any]] = []
    qualifying = [row for row in evaluations if row["entity_type"] == "virus_target" and row["qualifies_10x_50pct_instrain"]]
    if qualifying:
        require_instrain_available(SIF_INSTRAIN, shutil.which("singularity"))
        instrain_version = tool_version([
            "singularity", "exec", str(SIF_INSTRAIN), "inStrain", "profile", "--version",
        ])
        for target in qualifying:
            target_sha = target["entity_sha256"]
            representative = target["reference_contig"]
            focal_fasta = out_dir / f"focal_{target_sha[:16]}.fna"
            sequence = fetch_faidx_sequence(runner, ref_fna, representative)
            atomic_text(focal_fasta, f">{representative}\n{sequence}\n")
            profile_dir = out_dir / f"inStrain_{args.instrain_mode}_{target_sha[:16]}"
            if profile_dir.exists():
                raise FileExistsError(f"inStrain output already exists; refusing unverified overwrite: {profile_dir}")
            runner.run([
                "singularity", "exec", "-B", "/project", "-B", "/scratch", str(SIF_INSTRAIN),
                "inStrain", "profile", str(filtered_bam), str(focal_fasta), "-o", str(profile_dir),
                "-p", str(args.threads), "--min_mapq", "19", "--min_read_ani", str(min_read_ani),
            ], f"inStrain {args.instrain_mode} profile {target_sha}")
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
                "analysis_mode": args.instrain_mode,
                **callable_report,
                "output_sha256": output_hashes,
            })

    receipt = {
        "status": "completed_verified",
        "schema_version": 3,
        "stage": "stage4_map_sample_competitive",
        "sample": sample,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": time.time() - start,
        "execution_environment": {
            "hostname": platform.node(), "kernel": platform.release(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", "local"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID", "local"),
            "python_version": sys.version.split()[0], "tool_versions": versions,
        },
        "inputs": {
            "r1": str(r1), "r1_bytes": r1.stat().st_size, "r1_sha256": read_hashes["r1_sha256"],
            "r2": str(r2), "r2_bytes": r2.stat().st_size, "r2_sha256": read_hashes["r2_sha256"],
            "se": str(se) if se.exists() else "not_present", "se_sha256": read_hashes["se_sha256"],
            "pair_ledger": str(pair_path), "pair_ledger_sha256": stream_sha256(pair_path),
            "boundary_config": str(args.boundary_config), "boundary_config_sha256": stream_sha256(args.boundary_config),
            "mapping_script": str(Path(__file__)), "mapping_script_sha256": stream_sha256(Path(__file__)),
            "common_script_sha256": stream_sha256(Path(__file__).with_name("stage4_competitive_common.py")),
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
            "analysis_mode": args.instrain_mode,
            "min_read_ani_strict_greater_than": min_read_ani,
            "relaxed_crosssample_sensitivity_is_separate": True,
            "results": instrain_results,
            "recombination_analysis": "not_performed; SNVs/sequence polymorphism are not recombination evidence",
        },
        "focal_evaluations": evaluations,
        "outputs": {
            "raw_bam": str(raw_bam), "raw_bam_bytes": raw_bam.stat().st_size, "raw_bam_sha256": stream_sha256(raw_bam),
            "filtered_bam": str(filtered_bam), "filtered_bam_sha256": stream_sha256(filtered_bam),
            "coverage_tsv": str(coverage_path), "coverage_tsv_sha256": stream_sha256(coverage_path),
            "focal_evaluation_tsv": str(evaluation_path), "focal_evaluation_tsv_sha256": stream_sha256(evaluation_path),
        },
    }
    atomic_json(mapping_receipt_path, receipt)
    atomic_text(done_path, json.dumps({"receipt_sha256": stream_sha256(mapping_receipt_path)}) + "\n")
    print(f"[Completed verified] {sample}: {mapping_receipt_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
