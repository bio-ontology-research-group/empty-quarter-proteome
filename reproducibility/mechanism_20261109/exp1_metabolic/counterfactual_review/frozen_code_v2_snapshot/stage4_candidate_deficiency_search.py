#!/usr/bin/env python3
"""Stage 4 of `exp1_metabolic`: Candidate-Conditioned Missing-Enzyme & Transporter Verification
(`EXP1-CONTRACT.txt` Item 4 & `EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt` Item B).

Audits candidate metabolic deficiencies across 5 transparent evidence axes:
1. Reference Sequence Availability: explicit audit of reaction reference protein sets (no silent fallback).
2. Refined MAG Homology: DIAMOND blastx against refined MAG contigs, tracking contig-edge distance,
   alignment length, and subject coverage to distinguish true truncated/split genes from complete near-edge genes.
3. Raw Metagenome Assembly Context: DIAMOND blastx against unbinned sample assemblies, classifying hits
   without confirmed host linkage as 'unresolved_sample_community_homolog_host_unassigned' rather than
   unwarranted MAG dropout.
4. Synteny & Pathway-Hole Genomic Context: pre-computed pathway-hole dark gene candidates.
5. Physiological & Epistemic Qualification: tested no-hits are recorded as conditional negative evidence,
   preserving candidate-conditioned search scope without asserting biological auxotrophy.

Authoritative raw search TSVs and provenance receipts are preserved under results/provenance/.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import tarfile
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def compute_file_sha256(path: Path) -> str:
    """Compute sha256 checksum of a file."""
    if not path.exists():
        return ""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(65536):
            h.update(chunk)
    return h.hexdigest()


def load_pathway_hole_synteny_candidates(
    gap_context_tsv_gz: Path,
    c18_candidates_tsv: Path,
) -> tuple[dict[tuple[str, str], list[str]], set[tuple[str, str]]]:
    """Load pre-computed pathway-hole neighbourhood dark-gene candidates keyed by `(genome, pathway_or_rxn)`."""
    synteny_by_genome_pwy: dict[tuple[str, str], list[str]] = {}
    synteny_genome_any: set[str] = set()

    if gap_context_tsv_gz.exists():
        with gzip.open(gap_context_tsv_gz, "rt", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                g = (row.get("genome") or row.get("mag") or "").strip()
                pwy = (row.get("pathway") or row.get("pathway_id") or "").strip()
                gene = (row.get("gene_id") or row.get("candidate_gene") or row.get("orf_id") or "").strip()
                if g and pwy:
                    synteny_by_genome_pwy.setdefault((g, pwy), []).append(gene)
                    synteny_genome_any.add(g)

    if c18_candidates_tsv.exists():
        with open(c18_candidates_tsv, "r", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                g = (row.get("genome") or row.get("mag") or "").strip()
                pwy = (row.get("pathway") or row.get("pathway_id") or row.get("rxn") or "").strip()
                gene = (row.get("gene_id") or row.get("candidate_gene") or row.get("protein_id") or "").strip()
                if g and pwy:
                    synteny_by_genome_pwy.setdefault((g, pwy), []).append(gene)
                    synteny_genome_any.add(g)

    return synteny_by_genome_pwy, {(g, p) for (g, p) in synteny_by_genome_pwy}


def build_candidate_rxn_fasta_and_inventory(
    rxn_base_ids: list[str],
    gapseq_rxn_fasta_dir: Path,
    out_fasta: Path,
) -> tuple[int, dict[str, dict[str, Any]]]:
    """Combine reference protein sequences from `/project/k10018/sw/db/gapseq_seq/Archaea/rxn/<rxn>.fasta`
    for all target reactions, tracking sequence availability and checksum per reaction.
    """
    n_seqs = 0
    ref_inventory: dict[str, dict[str, Any]] = {}
    seen_rxns: set[str] = set()

    with open(out_fasta, "w", encoding="utf-8") as out_fh:
        for bid in sorted(set(rxn_base_ids)):
            clean_rxn = re.sub(r"_[ce]\d*$", "", str(bid))
            if clean_rxn in seen_rxns:
                continue
            seen_rxns.add(clean_rxn)

            fa_path = gapseq_rxn_fasta_dir / f"{clean_rxn}.fasta"
            if not fa_path.exists():
                ref_inventory[clean_rxn] = {
                    "reference_available": False,
                    "n_reference_sequences": 0,
                    "fasta_path": str(fa_path),
                    "fasta_sha256": "",
                }
                continue

            rxn_seq_count = 0
            for line in fa_path.read_text(encoding="utf-8", errors="replace").splitlines():
                if line.startswith(">"):
                    hdr = line[1:].strip().split()[0]
                    out_fh.write(f">{clean_rxn}|{hdr}\n")
                    n_seqs += 1
                    rxn_seq_count += 1
                else:
                    seq_clean = re.sub(r"[^A-Za-z*]", "", line.strip())
                    if seq_clean:
                        out_fh.write(seq_clean + "\n")

            ref_inventory[clean_rxn] = {
                "reference_available": rxn_seq_count > 0,
                "n_reference_sequences": rxn_seq_count,
                "fasta_path": str(fa_path),
                "fasta_sha256": compute_file_sha256(fa_path),
            }

    return n_seqs, ref_inventory


def extract_refined_mags_and_run_diamond_blastx(
    refined_mags_tar: Path,
    db_dmnd: Path,
    raw_searches_dir: Path,
    tmp_dir: Path,
    threads: int,
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[tuple[str, str], int], str]:
    """Extract refined MAG contigs, record contig lengths, and run `diamond blastx` against reference database.
    Retains authoritative raw TSV in `raw_searches_dir` and captures alignment details (coverage, coordinates).
    """
    mags_fasta = tmp_dir / "all_refined_mags.fna"
    contig_lengths: dict[tuple[str, str], int] = {}

    with tarfile.open(refined_mags_tar, "r:gz") as tf, open(mags_fasta, "w", encoding="utf-8") as out_fh:
        for member in tf.getmembers():
            if not member.isfile() or not (member.name.endswith(".fa") or member.name.endswith(".fna") or member.name.endswith(".fasta")):
                continue
            g_name = Path(member.name).stem
            fobj = tf.extractfile(member)
            if fobj is None:
                continue
            text = fobj.read().decode("utf-8", errors="replace")
            curr_cid = ""
            curr_len = 0
            for ln in text.splitlines():
                if ln.startswith(">"):
                    if curr_cid:
                        contig_lengths[(g_name, curr_cid)] = curr_len
                    curr_cid = ln[1:].strip().split()[0]
                    curr_len = 0
                    out_fh.write(f">{g_name}||{curr_cid}\n")
                else:
                    s = ln.strip()
                    curr_len += len(s)
                    out_fh.write(s + "\n")
            if curr_cid:
                contig_lengths[(g_name, curr_cid)] = curr_len

    out_tsv = raw_searches_dir / "mag_contigs_vs_target_rxns.tsv"
    cmd = [
        "diamond",
        "blastx",
        "--query",
        str(mags_fasta),
        "--db",
        str(db_dmnd),
        "--out",
        str(out_tsv),
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "qstart",
        "qend",
        "qlen",
        "sstart",
        "send",
        "slen",
        "evalue",
        "bitscore",
        "--evalue",
        "1e-3",
        "--max-target-seqs",
        "100",
        "--threads",
        str(threads),
        "--tmpdir",
        str(tmp_dir),
    ]
    subprocess.run(cmd, check=True)

    best_hits: dict[tuple[str, str], dict[str, Any]] = {}
    if out_tsv.exists():
        with open(out_tsv, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                qseqid, sseqid, pident, length, qstart, qend, qlen, sstart, send, slen, evalue, bitscore = parts
                if "||" not in qseqid or "|" not in sseqid:
                    continue
                g_name, cid = qseqid.split("||", 1)
                rxn_id = sseqid.split("|", 1)[0]
                bs = float(bitscore)
                qs, qe, ql = int(qstart), int(qend), int(qlen)
                ss, se, sl = int(sstart), int(send), int(slen)
                aln_len = int(length)

                lo, hi = min(qs, qe), max(qs, qe)
                dist_edge = min(lo - 1, max(0, ql - hi))
                is_near_edge = dist_edge <= 300
                subj_cov = float(aln_len) / float(max(1, sl))
                is_truncated_near_edge = is_near_edge and (subj_cov < 0.70)
                is_full_length_near_edge = is_near_edge and (subj_cov >= 0.70)

                key = (g_name, rxn_id)
                prev = best_hits.get(key)
                if prev is None or bs > prev["bitscore"]:
                    best_hits[key] = {
                        "contig_id": cid,
                        "bitscore": bs,
                        "pident": float(pident),
                        "evalue": float(evalue),
                        "alignment_length": aln_len,
                        "subject_length": sl,
                        "subject_coverage": subj_cov,
                        "qstart": qs,
                        "qend": qe,
                        "qlen": ql,
                        "sstart": ss,
                        "send": se,
                        "dist_to_contig_edge_bp": dist_edge,
                        "is_near_contig_edge": is_near_edge,
                        "is_truncated_near_edge": is_truncated_near_edge,
                        "is_full_length_near_edge": is_full_length_near_edge,
                    }

    tsv_hash = compute_file_sha256(out_tsv)
    return best_hits, contig_lengths, tsv_hash


def search_sample_assemblies_for_community_homologs(
    samples: list[str],
    assemblies_dir: Path,
    db_dmnd: Path,
    raw_searches_dir: Path,
    tmp_dir: Path,
    threads: int,
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, bool], dict[str, str]]:
    """Search unbinned sample assemblies against target reactions database.
    Retains raw TSVs and audits assembly existence explicitly per sample.
    """
    sample_hits: dict[tuple[str, str], dict[str, Any]] = {}
    sample_availability: dict[str, bool] = {}
    sample_tsv_hashes: dict[str, str] = {}

    for sample in sorted(set(samples)):
        contigs_fa = assemblies_dir / sample / f"{sample}_contigs.fa"
        if not contigs_fa.exists():
            sample_availability[sample] = False
            continue
        sample_availability[sample] = True

        out_tsv = raw_searches_dir / f"sample_{sample}_vs_target_rxns.tsv"
        cmd = [
            "diamond",
            "blastx",
            "--query",
            str(contigs_fa),
            "--db",
            str(db_dmnd),
            "--out",
            str(out_tsv),
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "evalue",
            "bitscore",
            "--evalue",
            "1e-5",
            "--max-target-seqs",
            "100",
            "--threads",
            str(threads),
            "--tmpdir",
            str(tmp_dir),
        ]
        subprocess.run(cmd, check=True)
        sample_tsv_hashes[sample] = compute_file_sha256(out_tsv)

        if out_tsv.exists():
            with open(out_tsv, "r", encoding="utf-8") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 6:
                        continue
                    qseqid, sseqid, pident, aln_len, evalue, bitscore = parts
                    if "|" not in sseqid:
                        continue
                    rxn_id = sseqid.split("|", 1)[0]
                    bs = float(bitscore)
                    key = (sample, rxn_id)
                    prev = sample_hits.get(key)
                    if prev is None or bs > prev["sample_assembly_bitscore"]:
                        sample_hits[key] = {
                            "sample_contig_id": qseqid,
                            "sample_assembly_bitscore": bs,
                            "sample_assembly_pident": float(pident),
                            "sample_assembly_evalue": float(evalue),
                            "sample_alignment_length": int(aln_len),
                        }

    return sample_hits, sample_availability, sample_tsv_hashes


def main() -> None:
    ap = argparse.ArgumentParser(description="Stage 4: Candidate-Conditioned Deficiency Verification.")
    ap.add_argument("--instances-tsv", type=Path, required=True)
    ap.add_argument("--summary-tsv", type=Path, required=True)
    ap.add_argument("--gapseq-rxn-fasta-dir", type=Path, default=Path("/project/k10018/sw/db/gapseq_seq/Archaea/rxn"))
    ap.add_argument("--refined-mags-tar", type=Path, default=Path("/scratch/project/k10018/eq-science-20261108/inputs/thermo237_refined_mags.tar.gz"))
    ap.add_argument("--assemblies-dir", type=Path, default=Path("/project/k10018/empty-quarter/data/assemblies"))
    ap.add_argument("--gap-context-tsv-gz", type=Path, default=Path("/scratch/project/k10018/eq-archaeal-function-20261106/results/gap_context_candidates.tsv.gz"))
    ap.add_argument("--c18-candidates-tsv", type=Path, default=Path("/project/k10018/work/eq-annot-20261030/results/c18_pathway_holes/candidates.tsv"))
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--tmp-dir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=64)
    args = ap.parse_args()

    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.tmp_dir.mkdir(parents=True, exist_ok=True)
    prov_dir = args.out_dir / "provenance"
    prov_dir.mkdir(parents=True, exist_ok=True)
    raw_searches_dir = prov_dir / "raw_searches"
    raw_searches_dir.mkdir(parents=True, exist_ok=True)

    inst_df = pd.read_csv(args.instances_tsv, sep="\t") if args.instances_tsv.exists() and args.instances_tsv.stat().st_size > 10 else pd.DataFrame()

    # Empty candidate handling per Root Item B: emit schema-valid completed nonestimable outputs
    if inst_df.empty:
        empty_inst_cols = [
            "scenario_id", "rxn_id", "base_id", "rxn_name", "is_transport", "biomass_sensitivity_classification",
            "followup_selection_tier", "epistemic_status", "genome", "cohort", "canonical_order",
            "canonical_family", "canonical_genus", "completeness", "contamination", "sample_id",
            "clean_rxn_id", "reference_sequence_status", "mag_search_status", "mag_tblastn_bitscore",
            "mag_tblastn_contig_id", "mag_dist_to_contig_edge_bp", "mag_hit_near_contig_edge",
            "mag_hit_truncated_near_edge", "sample_assembly_status", "sample_assembly_search_status",
            "sample_assembly_bitscore", "sample_assembly_contig_id", "verified_host_linkage",
            "n_synteny_dark_candidates", "synteny_dark_candidates_csv", "deficiency_verification_disposition"
        ]
        pd.DataFrame(columns=empty_inst_cols).to_csv(args.out_dir / "exp1_stage4_deficiency_verification_instances.tsv", sep="\t", index=False)
        pd.DataFrame(columns=["clean_rxn_id", "rxn_name", "deficiency_verification_disposition", "n_instances"]).to_csv(args.out_dir / "exp1_stage4_deficiency_verification_summary.tsv", sep="\t", index=False)
        with open(prov_dir / "stage4_execution_receipts.json", "w", encoding="utf-8") as fh:
            json.dump({"status": "completed_zero_candidates", "n_candidates": 0, "elapsed_seconds": round(time.time() - t0, 2)}, fh, indent=2)
        print("[Stage 4 Complete] Zero candidates to verify; emitted schema-valid empty outputs.")
        return

    synteny_map, _ = load_pathway_hole_synteny_candidates(args.gap_context_tsv_gz, args.c18_candidates_tsv)
    target_base_ids = inst_df["base_id"].dropna().astype(str).unique().tolist()
    ref_fasta = args.tmp_dir / "stage4_target_rxns.faa"
    n_seqs, ref_inventory = build_candidate_rxn_fasta_and_inventory(target_base_ids, args.gapseq_rxn_fasta_dir, ref_fasta)

    mag_hits: dict[tuple[str, str], dict[str, Any]] = {}
    sample_hits: dict[tuple[str, str], dict[str, Any]] = {}
    sample_availability: dict[str, bool] = {}
    mag_tsv_hash = ""
    sample_tsv_hashes: dict[str, str] = {}

    if n_seqs > 0 and args.refined_mags_tar.exists():
        db_dmnd = args.tmp_dir / "stage4_target_rxns.dmnd"
        subprocess.run(
            ["diamond", "makedb", "--in", str(ref_fasta), "--db", str(db_dmnd), "--threads", str(args.threads)],
            check=True,
        )
        mag_hits, _, mag_tsv_hash = extract_refined_mags_and_run_diamond_blastx(
            refined_mags_tar=args.refined_mags_tar,
            db_dmnd=db_dmnd,
            raw_searches_dir=raw_searches_dir,
            tmp_dir=args.tmp_dir,
            threads=args.threads,
        )
        samples = sorted({
            str(row.get("sample_id") or (str(row["genome"]).split("__")[0] if "__" in str(row["genome"]) else str(row["genome"])))
            for _, row in inst_df.iterrows()
        })
        sample_hits, sample_availability, sample_tsv_hashes = search_sample_assemblies_for_community_homologs(
            samples=samples,
            assemblies_dir=args.assemblies_dir,
            db_dmnd=db_dmnd,
            raw_searches_dir=raw_searches_dir,
            tmp_dir=args.tmp_dir,
            threads=args.threads,
        )

    verified_rows: list[dict[str, Any]] = []
    for _, row in inst_df.iterrows():
        g = str(row["genome"])
        sample = str(row.get("sample_id") or (g.split("__")[0] if "__" in g else g))
        bid = str(row["base_id"])
        clean_rxn = re.sub(r"_[ce]\d*$", "", str(bid))
        pwy = str(row.get("pathway_id", "") or "").strip("|")

        # Axis 1: Reference Availability
        ref_meta = ref_inventory.get(clean_rxn, {})
        has_ref = ref_meta.get("reference_available", False)

        # Axis 2: Refined MAG Contigs Search
        mh = mag_hits.get((g, clean_rxn), {})
        mag_evaluated = bool(has_ref and args.refined_mags_tar.exists())
        mag_bs = float(mh.get("bitscore", 0.0)) if mag_evaluated else math.nan
        near_edge = bool(mh.get("is_near_contig_edge", False)) if mag_evaluated else False
        trunc_near_edge = bool(mh.get("is_truncated_near_edge", False)) if mag_evaluated else False
        full_near_edge = bool(mh.get("is_full_length_near_edge", False)) if mag_evaluated else False
        dist_edge = mh.get("dist_to_contig_edge_bp", -1) if mag_evaluated else -1
        subj_cov = mh.get("subject_coverage", math.nan) if mag_evaluated else math.nan

        # Axis 3: Unbinned Sample Metagenome Assembly Search
        sample_fa_exists = sample_availability.get(sample, False)
        sample_evaluated = bool(has_ref and sample_fa_exists)
        sh = sample_hits.get((sample, clean_rxn), {}) if sample_evaluated else {}
        samp_bs = float(sh.get("sample_assembly_bitscore", 0.0)) if sample_evaluated else math.nan
        verified_host_linkage = False  # Contig provenance unverified from raw assembly

        # Axis 4: Synteny & Genomic Context Dark Genes
        synt_genes = synteny_map.get((g, pwy), []) or synteny_map.get((g, f"|{pwy}|"), []) or synteny_map.get((g, clean_rxn), [])

        # Axis 5: Qualified Competing-Hypothesis Disposition (`EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt` Item B)
        if not has_ref:
            disposition = "not_evaluable_missing_reference_sequences"
        elif not args.refined_mags_tar.exists():
            disposition = "not_evaluable_missing_mag_assembly"
        elif not math.isnan(mag_bs) and mag_bs >= 35.0:
            if trunc_near_edge:
                disposition = "candidate_truncated_or_split_gene_at_contig_edge"
            elif full_near_edge:
                disposition = "complete_homolog_near_contig_edge"
            elif mag_bs >= 50.0 and mh.get("pident", 0.0) >= 40.0:
                disposition = "uncalled_or_divergent_homolog_detected_in_mag"
            else:
                disposition = "weak_homology_or_conserved_domain_in_mag"
        elif not math.isnan(samp_bs) and samp_bs >= 65.0:
            if verified_host_linkage:
                disposition = "demonstrated_assembly_or_binning_dropout_with_host_linkage"
            else:
                disposition = "unresolved_sample_community_homolog_host_unassigned"
        elif len(synt_genes) > 0:
            disposition = "candidate_novel_or_non_orthologous_replacement_in_synteny"
        elif not sample_fa_exists:
            disposition = "tested_no_homolog_detected_in_mag_sample_assembly_unevaluated"
        else:
            disposition = "tested_no_homolog_detected_conditional_negative_evidence"

        out_r = dict(row)
        out_r["sample_id"] = sample
        out_r["clean_rxn_id"] = clean_rxn
        out_r["reference_sequence_status"] = "reference_available" if has_ref else "reference_unavailable"
        out_r["n_reference_sequences"] = ref_meta.get("n_reference_sequences", 0)
        out_r["mag_search_status"] = "evaluated" if mag_evaluated else "not_evaluable"
        out_r["mag_tblastn_bitscore"] = mag_bs
        out_r["mag_tblastn_contig_id"] = mh.get("contig_id", "")
        out_r["mag_tblastn_subject_coverage"] = subj_cov
        out_r["mag_dist_to_contig_edge_bp"] = dist_edge
        out_r["mag_hit_near_contig_edge"] = near_edge
        out_r["mag_hit_truncated_near_edge"] = trunc_near_edge
        out_r["sample_assembly_status"] = "sample_assembly_available" if sample_fa_exists else "sample_assembly_unavailable"
        out_r["sample_assembly_search_status"] = "evaluated" if sample_evaluated else "not_evaluable"
        out_r["sample_assembly_bitscore"] = samp_bs
        out_r["sample_assembly_contig_id"] = sh.get("sample_contig_id", "")
        out_r["verified_host_linkage"] = verified_host_linkage
        out_r["n_synteny_dark_candidates"] = len(synt_genes)
        out_r["synteny_dark_candidates_csv"] = ",".join(synt_genes[:5])
        out_r["deficiency_verification_disposition"] = disposition
        verified_rows.append(out_r)

    ver_df = pd.DataFrame(verified_rows)
    ver_df.to_csv(args.out_dir / "exp1_stage4_deficiency_verification_instances.tsv", sep="\t", index=False)

    if not ver_df.empty:
        ver_sum = (
            ver_df.groupby(
                ["clean_rxn_id", "rxn_name", "biomass_sensitivity_classification", "deficiency_verification_disposition"],
                dropna=False,
            )
            .agg(
                n_instances=("genome", "count"),
                n_unique_genomes=("genome", "nunique"),
                n_unique_families=("canonical_family", "nunique") if "canonical_family" in ver_df.columns else ("genome", "count"),
                median_mag_bitscore=("mag_tblastn_bitscore", "median"),
                median_sample_assembly_bitscore=("sample_assembly_bitscore", "median"),
            )
            .reset_index()
            .sort_values(["n_unique_genomes", "n_instances"], ascending=[False, False])
        )
    else:
        ver_sum = pd.DataFrame()
    ver_sum.to_csv(args.out_dir / "exp1_stage4_deficiency_verification_summary.tsv", sep="\t", index=False)

    receipts = {
        "status": "completed",
        "n_candidate_instances": len(ver_df),
        "n_reactions_evaluated": int(ver_df["clean_rxn_id"].nunique()) if not ver_df.empty else 0,
        "n_reference_sequences_total": n_seqs,
        "mag_search_raw_tsv_sha256": mag_tsv_hash,
        "sample_search_tsv_hashes": sample_tsv_hashes,
        "ref_inventory": ref_inventory,
        "elapsed_seconds": round(time.time() - t0, 2),
    }
    with open(prov_dir / "stage4_execution_receipts.json", "w", encoding="utf-8") as fh:
        json.dump(receipts, fh, indent=2)

    dt = time.time() - t0
    print(f"[Stage 4 Complete in {dt:.1f}s] Verified {len(ver_df)} candidate deficiency instances across {int(ver_df['clean_rxn_id'].nunique()) if not ver_df.empty else 0} reactions.")


if __name__ == "__main__":
    main()
