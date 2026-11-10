#!/usr/bin/env python3
"""Extract reviewed target geNomad/CheckV rows and Pfam-A.hmm GA thresholds with full provenance.

Reads the completed Stage 2 fresh target annotations (`stage2_raw_outputs_v2/fresh_target_annotations`)
and the completed Slurm job 15894662 Pfam domtblout (`stage4_unannotated_proteins_pfam_domtbl.txt`),
verifies the production Pfam-A.hmm SHA-256 recorded by job 15894662
(`a78a42d6faf265b6bfca59e8f062d06fae6083ce2c6e335d7b381f20b82b7903`), and writes
`results/stage4_target_genomad_checkv_pfam_ga_extract.json` with full file hashes and
HMM profile header provenance (`EXP3-ORF10-REFERENCE-CORRECTION.txt`).
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
from typing import Any

EXPECTED_PRODUCTION_PFAM_DB_SHA256 = "a78a42d6faf265b6bfca59e8f062d06fae6083ce2c6e335d7b381f20b82b7903"

REVIEWED_TARGETS = [
    "41Sr3::41Sr3_56779",
    "44Dr2::44Dr2_4002037",
    "44Dr2::44Dr2_4002037|provirus_1_34828",
    "48Sr2::48Sr2_190170",
    "55Dr2::55Dr2_75333",
    "9Sr2::9Sr2_284237",
    "V60PRr2::V60PRr2_715976",
    "coasm_48_PR::k127_4172870",
    "coasm_48_S::k127_1210458",
]


def stream_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_pfam_ga_thresholds(pfam_path: Path, target_accessions: set[str]) -> tuple[str, dict[str, dict[str, Any]]]:
    """Extract HMM header version and GA/TC/NC cutoffs for all matched Pfam accessions."""
    ga_map: dict[str, dict[str, Any]] = {}
    first_header: str | None = None
    cur_header = ""
    cur_name = ""
    cur_acc = ""
    cur_desc = ""
    cur_leng: int | None = None
    cur_date = ""
    cur_ga: tuple[float, float] | None = None
    cur_tc: tuple[float, float] | None = None
    cur_nc: tuple[float, float] | None = None

    def flush_record() -> None:
        if cur_acc and cur_acc in target_accessions and cur_ga is not None:
            ga_map[cur_acc] = {
                "name": cur_name,
                "acc": cur_acc,
                "desc": cur_desc,
                "model_length": cur_leng,
                "hmm_format_header": cur_header,
                "model_date": cur_date,
                "ga_seq_score": cur_ga[0],
                "ga_dom_score": cur_ga[1],
                "tc_seq_score": cur_tc[0] if cur_tc else None,
                "tc_dom_score": cur_tc[1] if cur_tc else None,
                "nc_seq_score": cur_nc[0] if cur_nc else None,
                "nc_dom_score": cur_nc[1] if cur_nc else None,
            }

    with open(pfam_path, errors="replace") as handle:
        for line in handle:
            if line.startswith("HMMER"):
                cur_header = line.strip()
                if first_header is None:
                    first_header = cur_header
                cur_name = cur_acc = cur_desc = cur_date = ""
                cur_leng = None
                cur_ga = cur_tc = cur_nc = None
            elif line.startswith("NAME "):
                cur_name = line.split(None, 1)[1].strip()
            elif line.startswith("ACC "):
                cur_acc = line.split(None, 1)[1].strip()
            elif line.startswith("DESC "):
                cur_desc = line.split(None, 1)[1].strip()
            elif line.startswith("LENG "):
                cur_leng = int(line.split(None, 1)[1].strip())
            elif line.startswith("DATE "):
                cur_date = line.split(None, 1)[1].strip()
            elif line.startswith("GA ") and cur_acc in target_accessions:
                parts = line.split()[1:]
                cur_ga = (float(parts[0].rstrip(";")), float(parts[1].rstrip(";")))
            elif line.startswith("TC ") and cur_acc in target_accessions:
                parts = line.split()[1:]
                cur_tc = (float(parts[0].rstrip(";")), float(parts[1].rstrip(";")))
            elif line.startswith("NC ") and cur_acc in target_accessions:
                parts = line.split()[1:]
                cur_nc = (float(parts[0].rstrip(";")), float(parts[1].rstrip(";")))
            elif line.startswith("//"):
                flush_record()

    if not first_header:
        raise ValueError(f"Could not read HMMER header from {pfam_path}")
    missing = sorted(target_accessions - set(ga_map))
    if missing:
        raise ValueError(f"Missing Pfam GA entries for accessions: {missing}")
    return first_header, ga_map


def main() -> int:
    exp3_root = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
    res_dir = exp3_root / "results"
    raw_dir = res_dir / "stage2_raw_outputs_v2/fresh_target_annotations"
    genes_path = raw_dir / "genomad_out/all_matched_primary_targets_annotate/all_matched_primary_targets_genes.tsv"
    tax_path = raw_dir / "genomad_out/all_matched_primary_targets_annotate/all_matched_primary_targets_taxonomy.tsv"
    qual_path = raw_dir / "checkv_out/quality_summary.tsv"
    contam_path = raw_dir / "checkv_out/contamination.tsv"
    domtbl_path = res_dir / "stage4_unannotated_proteins_pfam_domtbl.txt"
    raw_summary_path = res_dir / "stage4_sequence_domain_annotation_summary.json"
    pfam_path = Path(os.environ.get("PFAM_DB", "/project/k10018/sw/db/pfam-ibex/pfam/Pfam-A.hmm"))

    raw_summary = json.loads(raw_summary_path.read_text())
    recorded_pfam_sha = raw_summary["tool_provenance"]["pfam_db_sha256"]
    if recorded_pfam_sha != EXPECTED_PRODUCTION_PFAM_DB_SHA256:
        raise ValueError(
            f"Recorded Pfam DB SHA-256 in job 15894662 summary ({recorded_pfam_sha}) "
            f"does not match expected production SHA-256 ({EXPECTED_PRODUCTION_PFAM_DB_SHA256})"
        )

    target_set = set(REVIEWED_TARGETS)
    gene_rows: list[dict[str, str]] = []
    with open(genes_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cid = row["gene"].rsplit("_", 1)[0]
            if cid in target_set:
                gene_rows.append(row)

    tax_rows: list[dict[str, str]] = []
    with open(tax_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["seq_name"] in target_set:
                tax_rows.append(row)

    qual_rows: list[dict[str, str]] = []
    with open(qual_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["contig_id"] in target_set:
                qual_rows.append(row)

    contam_rows: list[dict[str, str]] = []
    with open(contam_path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["contig_id"] in target_set:
                contam_rows.append(row)

    matched_accs: set[str] = set()
    with open(domtbl_path) as handle:
        for line in handle:
            if not line.startswith("#") and line.strip():
                matched_accs.add(line.split()[1])

    hmm_header, ga_map = extract_pfam_ga_thresholds(pfam_path, matched_accs)

    out = {
        "schema_version": 2,
        "extraction_provenance": {
            "extraction_mode": "posthoc_metadata_extract_from_completed_stage2_and_slurm_15894662_artifacts",
            "no_new_slurm_or_hmmscan_job_submitted": True,
            "upstream_slurm_jobs": {
                "sequence_domain_annotation_job_id": "15894662",
                "sequence_domain_annotation_job_name": "exp3_seq_domains",
                "sequence_domain_annotation_node": "nid00193",
                "sequence_domain_annotation_elapsed": "00:02:07",
                "sequence_domain_annotation_exit_code": "0:0",
                "executed_sequence_script_relpath": "code/stage4_structural_hallmark_screen.py",
                "executed_sequence_script_sha256": raw_summary["script_sha256"],
            },
            "extract_script_relpath": "code/extract_stage4_target_genomad_checkv_pfam_ga.py",
            "extract_script_sha256": stream_sha256(Path(__file__)),
            "source_files": {
                "pfam_db_production_path": str(pfam_path),
                "pfam_db_sha256_from_job_15894662_summary": recorded_pfam_sha,
                "pfam_db_sha256_verified_expected": EXPECTED_PRODUCTION_PFAM_DB_SHA256,
                "pfam_hmm_format_header": hmm_header,
                "pfam_declared_release": raw_summary["tool_provenance"]["pfam"]["declared_release"],
                "genomad_genes_relpath": "results/stage2_raw_outputs_v2/fresh_target_annotations/genomad_out/all_matched_primary_targets_annotate/all_matched_primary_targets_genes.tsv",
                "genomad_genes_sha256": stream_sha256(genes_path),
                "genomad_taxonomy_relpath": "results/stage2_raw_outputs_v2/fresh_target_annotations/genomad_out/all_matched_primary_targets_annotate/all_matched_primary_targets_taxonomy.tsv",
                "genomad_taxonomy_sha256": stream_sha256(tax_path),
                "checkv_quality_summary_relpath": "results/stage2_raw_outputs_v2/fresh_target_annotations/checkv_out/quality_summary.tsv",
                "checkv_quality_summary_sha256": stream_sha256(qual_path),
                "checkv_contamination_relpath": "results/stage2_raw_outputs_v2/fresh_target_annotations/checkv_out/contamination.tsv",
                "checkv_contamination_sha256": stream_sha256(contam_path),
                "raw_domtbl_relpath": "results/stage4_unannotated_proteins_pfam_domtbl.txt",
                "raw_domtbl_sha256": stream_sha256(domtbl_path),
                "raw_sequence_summary_relpath": "results/stage4_sequence_domain_annotation_summary.json",
                "raw_sequence_summary_sha256": stream_sha256(raw_summary_path),
            },
        },
        "targets": sorted(target_set),
        "n_gene_rows": len(gene_rows),
        "n_matched_pfam_accessions": len(ga_map),
        "gene_rows": gene_rows,
        "tax_rows": tax_rows,
        "qual_rows": qual_rows,
        "contam_rows": contam_rows,
        "pfam_ga_map": dict(sorted(ga_map.items())),
    }
    out_path = res_dir / "stage4_target_genomad_checkv_pfam_ga_extract.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {out_path} ({len(gene_rows)} gene rows, {len(ga_map)} Pfam GA entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
