#!/usr/bin/env python3
"""Pfam sequence-domain annotation plus audited existing-model Foldseek follow-up.

HMMER/Pfam is sequence profile annotation, never a structural search.  All
unannotated ORFs on credible unique target sequences are extracted by exact
sequence ID.  Structural follow-up uses only exact-sequence matches to existing
ColabFold/ESMFold model DBs; no new structure predictions are created.  If no
exact model exists, the audited status is explicitly not performed.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import pandas as pd

from stage4_competitive_common import atomic_json, atomic_text, load_reviewed_pairs, stream_sha256

BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
RES_DIR = BRANCH_ROOT / "results"
RAW_DIR = RES_DIR / "stage2_raw_outputs_v2"
SIF_HMMER = Path("/project/k10018/sw/images/eq-annot/hmmer.sif")
PFAM_DB = Path("/project/k10018/sw/db/pfam-ibex/pfam/Pfam-A.hmm")
SIF_FOLDSEEK = Path("/project/k10018/sw/images/eq-annot/foldseek.sif")
MODEL_ROOT = Path("/project/k10018/work/eq-annot-20261030")
CURATED_STRUCTURE_DB = Path("/scratch/project/k10018/eq-mechanism-20261109/exp4_surface/tmp/foldseek_curated/curated_pdb_db")
MODEL_SOURCES = [
    ("ColabFold_AF2", MODEL_ROOT / "results/c10_structclust/colabfold/members.tsv.gz", MODEL_ROOT / "results/c10_structclust/colabfold/db/all"),
    ("ESMFold", MODEL_ROOT / "results/c10_structclust/esmfold/members.tsv.gz", MODEL_ROOT / "results/c10_structclust/esmfold/db/all"),
]
GENE_SUFFIX = re.compile(r"^(.*)_([0-9]+)$")


def run_command(cmd: list[str], description: str) -> subprocess.CompletedProcess[str]:
    print(f"[CMD] {description}: {shlex.join(cmd)}", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"{description} failed ({result.returncode}): {result.stderr}")
    return result


def probe_hmmer_version(sif_path: Path) -> str:
    result = run_command(["singularity", "exec", str(sif_path), "hmmscan", "-h"], "probe HMMER")
    lines = (result.stdout or result.stderr).splitlines()
    version = next((line.strip() for line in lines if line.startswith("# HMMER")), None)
    if not version:
        raise ValueError("could not determine actual HMMER version")
    return version


def probe_pfam_version(pfam_path: Path) -> dict[str, str]:
    with open(pfam_path, errors="strict") as handle:
        header = handle.readline().strip()
    if not header.startswith("HMMER"):
        raise ValueError(f"invalid Pfam HMM header: {header}")
    # This database has no RELEASE line; its cryptographic hash is the version ID.
    return {"declared_release": "not_declared_in_database", "hmm_format_header": header}


def is_unannotated_genomad_gene(row: pd.Series) -> bool:
    try:
        hallmark = int(float(row.get("virus_hallmark", 0)))
    except (ValueError, TypeError):
        hallmark = 0
    if hallmark == 1:
        return False
    for column in ["annotation_description", "marker", "annotation_conjscan", "annotation_amr", "annotation_accessions"]:
        value = str(row.get(column, "")).strip()
        if value not in {"", "nan", "NA", "None", "-"}:
            return False
    return True


def gene_sequence_id(gene_id: str) -> str:
    """Strip only geNomad's terminal numeric ORF suffix; never substring-match."""
    match = GENE_SUFFIX.fullmatch(gene_id)
    if not match:
        raise ValueError(f"geNomad gene lacks terminal numeric ORF suffix: {gene_id}")
    return match.group(1)


def extract_candidate_proteins(candidate_gene_ids: set[str], all_proteins_faa: Path, out_faa: Path) -> tuple[int, list[str]]:
    if not all_proteins_faa.exists():
        raise FileNotFoundError(f"all proteins FASTA missing: {all_proteins_faa}")
    extracted: list[str] = []
    write_current = False
    with open(all_proteins_faa) as source, open(out_faa, "w") as output:
        for line in source:
            if line.startswith(">"):
                gene_id = line[1:].strip().split()[0]
                write_current = gene_id in candidate_gene_ids
                if write_current:
                    if gene_id in extracted:
                        raise ValueError(f"duplicate candidate protein FASTA ID: {gene_id}")
                    extracted.append(gene_id)
                    output.write(line)
            elif write_current:
                output.write(line)
    if len(extracted) != len(candidate_gene_ids) or set(extracted) != candidate_gene_ids:
        missing = sorted(candidate_gene_ids - set(extracted))
        extra = sorted(set(extracted) - candidate_gene_ids)
        raise ValueError(f"exact protein extraction failed: missing={missing[:20]}, extra={extra[:20]}")
    return len(extracted), extracted


def fasta_sequences(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        if current is None:
            return
        sequence = "".join(chunks).upper().rstrip("*")
        if not sequence or current in sequences:
            raise ValueError(f"empty/duplicate sequence ID in {path}: {current}")
        sequences[current] = sequence

    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                finish()
                current = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    finish()
    return sequences


def parse_domtblout(domtbl_path: Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    with open(domtbl_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 23:
                raise ValueError(f"malformed HMMER domtblout row: {line}")
            rows.append({
                "query_gene": parts[3], "pfam_target_name": parts[0], "pfam_accession": parts[1],
                "full_evalue": float(parts[6]), "full_score": float(parts[7]),
                "domain_evalue": float(parts[12]), "domain_score": float(parts[13]),
                "pfam_description": " ".join(parts[22:]),
            })
    return pd.DataFrame(rows)


def load_model_polymer_sequences(db_prefix: Path, target_snames: set[str]) -> tuple[dict[str, str], dict[str, str]]:
    """Read only required model records from an existing Foldseek DB."""
    lookup_path = Path(str(db_prefix) + ".lookup")
    index_path = Path(str(db_prefix) + ".index")
    for path in [db_prefix, lookup_path, index_path]:
        if not path.exists():
            raise FileNotFoundError(f"existing model DB component missing: {path}")
    key_to_sname: dict[str, str] = {}
    with open(lookup_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[1] in target_snames:
                key_to_sname[fields[0]] = fields[1]
    spans: list[tuple[str, int, int]] = []
    with open(index_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] in key_to_sname:
                spans.append((key_to_sname[fields[0]], int(fields[1]), int(fields[2])))
    sequences: dict[str, str] = {}
    with open(db_prefix, "rb") as handle:
        for sname, offset, length in spans:
            handle.seek(offset)
            sequence = handle.read(length).rstrip(b"\x00\r\n ").decode("ascii").strip().upper().rstrip("*")
            if sname in sequences and sequences[sname] != sequence:
                raise ValueError(f"conflicting model polymer sequence for {sname}")
            sequences[sname] = sequence
    return sequences, key_to_sname


def audit_existing_models(candidate_sequences: dict[str, str]) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, Any]]:
    records: list[dict[str, Any]] = []
    exact: list[dict[str, Any]] = []
    provenance: dict[str, Any] = {}
    for source, members_path, db_prefix in MODEL_SOURCES:
        if not members_path.exists():
            raise FileNotFoundError(f"model inventory unavailable: {members_path}")
        memberships: list[tuple[str, str]] = []
        with gzip.open(members_path, "rt") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 6:
                    raise ValueError(f"malformed model inventory row in {members_path}: {line}")
                gene_id, sname = fields[1], fields[5]
                if gene_id in candidate_sequences:
                    memberships.append((gene_id, sname))
        polymers, key_map = load_model_polymer_sequences(db_prefix, {sname for _, sname in memberships})
        for gene_id in sorted(candidate_sequences):
            gene_memberships = [(g, s) for g, s in memberships if g == gene_id]
            if not gene_memberships:
                records.append({"gene_id": gene_id, "model_source": source, "model_id": pd.NA, "status": "no_existing_model_inventory_entry"})
                continue
            for _, sname in gene_memberships:
                polymer = polymers.get(sname, "")
                candidate = candidate_sequences[gene_id]
                is_exact = bool(polymer) and polymer == candidate
                records.append({
                    "gene_id": gene_id, "model_source": source, "model_id": sname,
                    "status": "exact_sequence_match" if is_exact else "model_sequence_mismatch_or_unavailable",
                    "candidate_sha256": hashlib.sha256(candidate.encode()).hexdigest(),
                    "model_sha256": hashlib.sha256(polymer.encode()).hexdigest() if polymer else pd.NA,
                })
                if is_exact:
                    keys = [key for key, value in key_map.items() if value == sname]
                    if len(keys) != 1:
                        raise ValueError(f"model DB key resolution failed for {source}/{sname}")
                    exact.append({"gene_id": gene_id, "model_source": source, "model_id": sname, "db_key": keys[0], "db_prefix": str(db_prefix)})
        provenance[source] = {
            "members_path": str(members_path), "members_sha256": stream_sha256(members_path),
            "db_prefix": str(db_prefix), "db_lookup_sha256": stream_sha256(Path(str(db_prefix) + ".lookup")),
            "candidate_inventory_entries": len(memberships),
        }
    return pd.DataFrame(records), exact, provenance


def run_foldseek_exact_models(exact_models: list[dict[str, Any]], threads: int) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    if not exact_models:
        return pd.DataFrame(), []
    if not SIF_FOLDSEEK.exists() or shutil.which("singularity") is None:
        raise RuntimeError("exact structural models exist but Foldseek/container is unavailable")
    for suffix in ["", ".lookup", ".index"]:
        if not Path(str(CURATED_STRUCTURE_DB) + suffix).exists():
            raise FileNotFoundError(f"curated Foldseek reference DB missing: {CURATED_STRUCTURE_DB}{suffix}")
    version = run_command(["singularity", "exec", str(SIF_FOLDSEEK), "foldseek", "version"], "probe Foldseek").stdout.strip()
    all_hits: list[pd.DataFrame] = []
    run_receipts: list[dict[str, Any]] = []
    columns = ["query", "target", "fident", "alnlen", "qlen", "tlen", "qcov", "tcov", "qtmscore", "ttmscore", "alntmscore", "evalue", "bits"]
    for source in sorted({row["model_source"] for row in exact_models}):
        models = [row for row in exact_models if row["model_source"] == source]
        db_prefix = Path(models[0]["db_prefix"])
        subset_list = RES_DIR / f"stage4_foldseek_{source}_exact_model_keys.txt"
        subset_list.write_text("".join(f"{row['db_key']}\n" for row in models))
        query_db = RES_DIR / f"stage4_foldseek_{source}_exact_models_db"
        alignment_db = RES_DIR / f"stage4_foldseek_{source}_alignment_db"
        scratch = RES_DIR / f"stage4_foldseek_{source}_work"
        scratch.mkdir(exist_ok=True)
        base_cmd = ["singularity", "exec", "-B", "/project", "-B", "/scratch", str(SIF_FOLDSEEK), "foldseek"]
        run_command(base_cmd + ["createsubdb", str(subset_list), str(db_prefix), str(query_db)], f"subset existing {source} models")
        run_command(base_cmd + ["search", str(query_db), str(CURATED_STRUCTURE_DB), str(alignment_db), str(scratch), "-a", "--threads", str(threads)], f"Foldseek {source}")
        hit_path = RES_DIR / f"stage4_foldseek_{source}_hits.tsv"
        run_command(base_cmd + [
            "convertalis", str(query_db), str(CURATED_STRUCTURE_DB), str(alignment_db), str(hit_path),
            "--format-output", ",".join(columns),
        ], f"format Foldseek {source} hits")
        hits = pd.read_csv(hit_path, sep="\t", names=columns) if hit_path.stat().st_size > 0 else pd.DataFrame(columns=columns)
        hits["model_source"] = source
        all_hits.append(hits)
        run_receipts.append({
            "model_source": source, "foldseek_version": version, "n_query_models": len(models),
            "query_db": str(query_db), "reference_db": str(CURATED_STRUCTURE_DB),
            "reference_lookup_sha256": stream_sha256(Path(str(CURATED_STRUCTURE_DB) + ".lookup")),
            "hits_path": str(hit_path), "hits_sha256": stream_sha256(hit_path), "n_hits": len(hits),
        })
    return (pd.concat(all_hits, ignore_index=True) if all_hits else pd.DataFrame()), run_receipts


def main() -> int:
    start = time.time()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--evalue-cutoff", type=float, default=1e-3)
    args = parser.parse_args()
    done = RES_DIR / "STAGE4_SEQUENCE_DOMAIN_ANNOTATION_DONE"
    done.unlink(missing_ok=True)
    pairs = load_reviewed_pairs(RES_DIR / "stage3_tier1_pairs_individual_evidence_audit.tsv")
    eligible = pairs[(pairs["analysis_role"] == "primary") & pairs["eligible_for_structural_hallmark_followup"]]
    target_to_sha = dict(zip(eligible["canonical_annotated_target_id"], eligible["target_seq_sha256"]))
    if not target_to_sha:
        raise ValueError("no credible exact targets approved for sequence follow-up")
    print(f"[Candidates] {len(target_to_sha)} unique credible target sequences", flush=True)

    annotation_dir = RAW_DIR / "fresh_target_annotations/genomad_out/all_matched_primary_targets_annotate"
    genes_path = annotation_dir / "all_matched_primary_targets_genes.tsv"
    proteins_path = annotation_dir / "all_matched_primary_targets_proteins.faa"
    genes = pd.read_csv(genes_path, sep="\t")
    if "gene" not in genes or genes["gene"].duplicated().any():
        raise ValueError("geNomad genes require unique gene IDs")
    sequence_ids = genes["gene"].astype(str).map(gene_sequence_id)
    target_genes = genes[sequence_ids.isin(target_to_sha)].copy()
    target_genes["exact_target_sequence_id"] = sequence_ids[sequence_ids.isin(target_to_sha)]
    unannotated = target_genes[target_genes.apply(is_unannotated_genomad_gene, axis=1)]
    candidate_ids = set(unannotated["gene"].astype(str))
    genes_per_target = {target: int((unannotated["exact_target_sequence_id"] == target).sum()) for target in target_to_sha}
    candidate_faa = RES_DIR / "stage4_unannotated_viral_candidates.faa"
    n_extracted, extracted_ids = extract_candidate_proteins(candidate_ids, proteins_path, candidate_faa)
    candidate_sequences = fasta_sequences(candidate_faa)
    if set(candidate_sequences) != candidate_ids or n_extracted != len(candidate_ids):
        raise ValueError("candidate sequence count/ID contract failed")

    model_audit, exact_models, model_provenance = audit_existing_models(candidate_sequences)
    model_audit_path = RES_DIR / "stage4_existing_structural_model_exact_sequence_audit.tsv"
    model_audit.to_csv(model_audit_path, sep="\t", index=False)
    foldseek_hits, foldseek_runs = run_foldseek_exact_models(exact_models, args.threads)
    foldseek_path = RES_DIR / "stage4_existing_model_foldseek_results.tsv"
    foldseek_hits.to_csv(foldseek_path, sep="\t", index=False)
    structural_status = "completed_existing_exact_models" if exact_models else "not_performed_no_exact_sequence_matched_existing_models"

    for path in [SIF_HMMER, PFAM_DB] + [Path(str(PFAM_DB) + suffix) for suffix in [".h3f", ".h3i", ".h3m", ".h3p"]]:
        if not path.exists() or path.stat().st_size == 0:
            raise FileNotFoundError(f"HMMER/Pfam required file missing: {path}")
    if shutil.which("singularity") is None:
        raise RuntimeError("singularity unavailable for HMMER")
    hmmer_version = probe_hmmer_version(SIF_HMMER)
    pfam_version = probe_pfam_version(PFAM_DB)
    domtbl = RES_DIR / "stage4_unannotated_proteins_pfam_domtbl.txt"
    tbl = RES_DIR / "stage4_unannotated_proteins_pfam_tbl.txt"
    command = [
        "singularity", "exec", "-B", "/project", "-B", "/scratch", str(SIF_HMMER), "hmmscan",
        "--cpu", str(args.threads), "-E", str(args.evalue_cutoff), "--domtblout", str(domtbl),
        "--tblout", str(tbl), str(PFAM_DB), str(candidate_faa),
    ]
    run_command(command, "Pfam sequence-domain annotation")
    hits = parse_domtblout(domtbl)
    if not hits.empty and not set(hits["query_gene"]) <= candidate_ids:
        raise ValueError("HMMER output contains non-candidate query IDs")
    hits_path = RES_DIR / "stage4_sequence_domain_annotation_results.tsv"
    hits.to_csv(hits_path, sep="\t", index=False)

    summary = {
        "status": "sequence_domain_annotation_completed_verified",
        "schema_version": 3,
        "stage": "stage4_sequence_domain_annotation_and_existing_model_structural_followup",
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "elapsed_seconds": time.time() - start,
        "methodology_distinction": {
            "pfam_search_type": "sequence_profile_hmm_domain_annotation",
            "is_structural_coordinate_search": False,
            "foldseek_structural_screen_status": structural_status,
            "no_new_structure_predictions": True,
        },
        "candidate_selection": {
            "credible_targets": target_to_sha, "n_candidate_orfs": len(candidate_ids),
            "n_candidate_proteins_extracted": n_extracted, "extracted_ids_sha256": hashlib.sha256("\n".join(sorted(extracted_ids)).encode()).hexdigest(),
            "genes_per_target": genes_per_target, "all_unannotated_orfs_no_truncation": True,
        },
        "tool_provenance": {
            "hmmer_version": hmmer_version, "hmmer_container": str(SIF_HMMER), "hmmer_container_sha256": stream_sha256(SIF_HMMER),
            "pfam": pfam_version, "pfam_db": str(PFAM_DB), "pfam_db_sha256": stream_sha256(PFAM_DB),
            "evalue_cutoff": args.evalue_cutoff, "command": command,
        },
        "pfam_results": {
            "n_hits": len(hits), "n_unique_proteins_with_hits": int(hits["query_gene"].nunique()) if not hits.empty else 0,
            "results_path": str(hits_path), "results_sha256": stream_sha256(hits_path),
        },
        "structural_followup": {
            "status": structural_status, "candidate_model_audit_path": str(model_audit_path),
            "candidate_model_audit_sha256": stream_sha256(model_audit_path),
            "n_exact_sequence_matched_models": len(exact_models), "exact_models": exact_models,
            "model_inventory_provenance": model_provenance, "foldseek_runs": foldseek_runs,
            "foldseek_results_path": str(foldseek_path), "foldseek_results_sha256": stream_sha256(foldseek_path),
        },
        "script_sha256": stream_sha256(Path(__file__)),
    }
    summary_path = RES_DIR / "stage4_sequence_domain_annotation_summary.json"
    atomic_json(summary_path, summary)
    atomic_text(done, json.dumps({"summary_sha256": stream_sha256(summary_path)}) + "\n")
    print(f"[Sequence annotation completed; structural follow-up={structural_status}] {summary_path}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
