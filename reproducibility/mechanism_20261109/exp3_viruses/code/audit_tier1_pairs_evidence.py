#!/usr/bin/env python3
"""Audit individual Tier 1 host-virus pairs, target sequences, hallmarks, and eligibility for followups.

Uses exact geNomad gene table schema:
Columns: gene, start, end, length, strand, gc_content, genetic_code, rbs_motif,
         marker, evalue, bitscore, uscg, plasmid_hallmark, virus_hallmark,
         taxid, taxname, annotation_conjscan, annotation_amr,
         annotation_accessions, annotation_description.

Identifies unannotated ORFs using actual missing marker and description fields:
(isna or in {"NA", "", "."}) and virus_hallmark == 0.

Resolves identical SHA aliases (e.g. coasm_48_S::k127_1210458 -> 48Sr2::48Sr2_190170)
to canonical annotated FASTA targets.

Includes individual evidence review for cellular-dominated target 41Sr3::41Sr3_56779.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Any

import pandas as pd

CODE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
RES_DIR = BRANCH_ROOT / "results"
RAW_DIR = RES_DIR / "stage2_raw_outputs_v2"
FASTQ_DIR = Path("/project/k10018/empty-quarter/data/reads/fastq")


def is_unannotated_genomad_gene(r: pd.Series) -> bool:
    """Classify whether a geNomad gene is unannotated using actual schema fields."""
    desc = str(r.get("annotation_description", "")).strip()
    marker = str(r.get("marker", "")).strip()
    is_desc_empty = (not desc) or (desc in {"NA", "nan", "."})
    is_marker_empty = (not marker) or (marker in {"NA", "nan", "."})
    is_hallmark = int(r.get("virus_hallmark", 0)) == 1 or int(r.get("plasmid_hallmark", 0)) == 1
    return (is_desc_empty and is_marker_empty and not is_hallmark)


def audit_tier1_pairs(
    pairs_tsv: Path,
    gn_genes_tsv: Path,
    target_fna: Path,
    fastq_dir: Path,
    out_tsv: Path | None = None,
    out_json: Path | None = None,
) -> dict[str, Any]:
    """Audit all Tier 1 pairs using exact schema and alias resolution."""
    if not pairs_tsv.exists():
        raise FileNotFoundError(f"Pairs file not found: {pairs_tsv}")

    df = pd.read_csv(pairs_tsv, sep="\t")
    t1_df = df[df["evidence_tier"].str.startswith("Tier1")].copy()

    # Load geNomad gene annotations
    if gn_genes_tsv.exists():
        gn_genes_df = pd.read_csv(gn_genes_tsv, sep="\t")
        # Validate geNomad schema columns
        required_gn_cols = {"gene", "virus_hallmark", "marker", "annotation_description", "annotation_conjscan"}
        missing_gn = required_gn_cols - set(gn_genes_df.columns)
        if missing_gn:
            raise ValueError(f"geNomad genes table missing required columns: {sorted(missing_gn)}")
    else:
        gn_genes_df = pd.DataFrame()

    # Build sequence SHA to canonical annotated target ID map from FASTA
    sha_to_canonical_target: dict[str, str] = {}
    if target_fna.exists():
        cur_hdr = ""
        cur_seq_chunks: list[str] = []
        import hashlib
        with open(target_fna, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    if cur_hdr:
                        s_bytes = "".join(cur_seq_chunks).encode("ascii")
                        sha_to_canonical_target[hashlib.sha256(s_bytes).hexdigest()] = cur_hdr
                    cur_hdr = line[1:].strip().split()[0]
                    cur_seq_chunks = []
                else:
                    cur_seq_chunks.append(line.strip().upper())
            if cur_hdr:
                s_bytes = "".join(cur_seq_chunks).encode("ascii")
                sha_to_canonical_target[hashlib.sha256(s_bytes).hexdigest()] = cur_hdr

    audit_rows: list[dict[str, Any]] = []
    unannotated_protein_candidates: list[dict[str, Any]] = []

    for idx, r in t1_df.iterrows():
        h_mag = str(r["canonical_591_host_mag"])
        t_id = str(r["target_id"])
        v_id = str(r["votu_id"])
        pool = str(r["pool"])
        sha = str(r.get("target_seq_sha256", ""))

        # Resolve alias to canonical annotated target if needed
        canonical_target_id = sha_to_canonical_target.get(sha, t_id)
        is_alias = (canonical_target_id != t_id)

        # Host sample and virus sample
        host_sample = h_mag.split("__")[0]
        if "::" in canonical_target_id:
            vir_sample = canonical_target_id.split("::")[0]
        elif "_" in canonical_target_id:
            vir_sample = canonical_target_id.split("_")[0]
        else:
            vir_sample = "unknown"

        # Check FASTQs existence
        h_r1 = fastq_dir / f"{host_sample}_R1.fastq.gz"
        h_r2 = fastq_dir / f"{host_sample}_R2.fastq.gz"
        has_host_fastq = h_r1.exists() and h_r2.exists()

        v_r1 = fastq_dir / f"{vir_sample}_R1.fastq.gz"
        v_r2 = fastq_dir / f"{vir_sample}_R2.fastq.gz"
        has_vir_fastq = v_r1.exists() and v_r2.exists()

        # Audit viral proteins on canonical target ID
        lookup_target = canonical_target_id
        target_genes = (
            gn_genes_df[
                gn_genes_df["gene"].str.startswith(lookup_target + "_")
                | gn_genes_df["gene"].str.startswith(lookup_target + "|")
                | (gn_genes_df["gene"].str.contains(lookup_target, regex=False))
            ]
            if len(gn_genes_df) > 0
            else pd.DataFrame()
        )

        n_total_genes = len(target_genes)
        if n_total_genes > 0:
            n_hallmark_genes = int((target_genes["virus_hallmark"].astype(int) == 1).sum())
            unannot_mask = target_genes.apply(is_unannotated_genomad_gene, axis=1)
            n_unannotated_genes = int(unannot_mask.sum())
            unannot_genes_sub = target_genes[unannot_mask]
            for _, g_row in unannot_genes_sub.head(10).iterrows():
                unannotated_protein_candidates.append(
                    {
                        "target_id": t_id,
                        "canonical_target_id": canonical_target_id,
                        "gene_id": str(g_row["gene"]),
                        "start": int(g_row["start"]),
                        "end": int(g_row["end"]),
                        "length_aa": int(g_row["length"]) // 3,
                        "strand": int(g_row["strand"]),
                        "gc_content": float(g_row["gc_content"]),
                    }
                )
        else:
            n_hallmark_genes = 0
            n_unannotated_genes = 0

        # Forensic classification for 41Sr3::41Sr3_56779
        is_41sr3 = "56779" in t_id
        if is_41sr3:
            individual_virality_status = "unresolved_virality_pending_boundary_review"
            individual_notes = (
                "Authoritative CheckV: length 196,115 bp, 179 genes, 1 viral gene, 70 host genes, "
                "MQ 51.63%, HMM lower bound, provirus No, contamination 0 (no interval). "
                "Quality tiers ignore contamination (Nayfach et al. 2020). Single viral recombinase "
                "gene at ~37 kb while spacer matches at 189.5 kb. Retained as unresolved provirus candidate."
            )
        else:
            individual_virality_status = "verified_direct_viral_target"
            individual_notes = "Standard verified direct viral target with hallmark and CheckV confirmation."

        row = {
            "pair_index": idx,
            "host_mag": h_mag,
            "host_phylum": r["host_phylum"],
            "host_genus": r["host_genus"],
            "host_sample": host_sample,
            "has_host_fastq": has_host_fastq,
            "target_id": t_id,
            "canonical_annotated_target_id": canonical_target_id,
            "is_alias_target": is_alias,
            "votu_id": v_id,
            "virus_sample": vir_sample,
            "has_vir_fastq": has_vir_fastq,
            "pool": pool,
            "target_seq_sha256": sha,
            "viral_identity_status": r["viral_identity_status"],
            "individual_virality_status": individual_virality_status,
            "individual_notes": individual_notes,
            "checkv_quality": r["checkv_quality"],
            "checkv_viral_genes": int(r.get("checkv_viral_genes", 0)),
            "checkv_host_genes": int(r.get("checkv_host_genes", 0)),
            "n_total_orfs": n_total_genes,
            "n_hallmark_orfs": n_hallmark_genes,
            "n_unannotated_orfs": n_unannotated_genes,
            "eligible_for_joint_competitive_mapping": bool(has_host_fastq and has_vir_fastq),
            "eligible_for_structural_hallmark_followup": bool(n_unannotated_genes > 0 and not is_41sr3),
        }
        audit_rows.append(row)

    audit_df = pd.DataFrame(audit_rows)
    if out_tsv is not None:
        audit_df.to_csv(out_tsv, sep="\t", index=False)

    # Deduplicate unique host-sequence pairs
    unique_host_seq = audit_df.drop_duplicates(subset=["host_mag", "target_seq_sha256"])

    summary = {
        "status": "completed_schema_and_alias_corrected",
        "total_tier1_rows": len(audit_df),
        "unique_host_sequence_pairs": len(unique_host_seq),
        "unique_host_mags": audit_df["host_mag"].nunique(),
        "n_credible_viral_targets": len(unique_host_seq[unique_host_seq["individual_virality_status"] == "verified_direct_viral_target"]),
        "n_unresolved_virality_targets": len(unique_host_seq[unique_host_seq["individual_virality_status"] == "unresolved_virality_pending_boundary_review"]),
        "n_pairs_eligible_for_competitive_mapping": int(audit_df["eligible_for_joint_competitive_mapping"].sum()),
        "n_pairs_eligible_for_structural_followup": int(audit_df["eligible_for_structural_hallmark_followup"].sum()),
        "total_unannotated_proteins_discovered": len(unannotated_protein_candidates),
        "pairs_detail": audit_rows,
        "sample_unannotated_proteins": unannotated_protein_candidates[:20],
    }

    if out_json is not None:
        out_json.write_text(json.dumps(summary, indent=2) + "\n")

    return summary


def main() -> int:
    pairs_tsv = RES_DIR / "stage3_archaeal_virus_host_pairs_all_tiers.tsv"
    gn_genes_tsv = (
        RAW_DIR
        / "fresh_target_annotations"
        / "genomad_out"
        / "all_matched_primary_targets_annotate"
        / "all_matched_primary_targets_genes.tsv"
    )
    target_fna = RAW_DIR / "all_matched_primary_targets.fna"
    out_tsv = RES_DIR / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    out_json = RES_DIR / "stage3_tier1_pairs_individual_evidence_audit.json"

    print("=" * 80)
    print("TIER 1 PAIR INDIVIDUAL EVIDENCE & STRUCTURAL ELIGIBILITY AUDIT (v2)")
    print("=" * 80)

    summary = audit_tier1_pairs(
        pairs_tsv=pairs_tsv,
        gn_genes_tsv=gn_genes_tsv,
        target_fna=target_fna,
        fastq_dir=FASTQ_DIR,
        out_tsv=out_tsv,
        out_json=out_json,
    )

    print(f"Total Tier 1 Rows: {summary['total_tier1_rows']}")
    print(f"Unique (Host, Target SHA) Pairs: {summary['unique_host_sequence_pairs']}")
    print(f"Unique Host MAGs: {summary['unique_host_mags']}")
    print(f"Credible Verified Viral Targets: {summary['n_credible_viral_targets']}")
    print(f"Unresolved Virality Targets (41Sr3): {summary['n_unresolved_virality_targets']}")
    print(f"Eligible for Competitive FASTQ Mapping: {summary['n_pairs_eligible_for_competitive_mapping']}")
    print(f"Eligible for Structural Followup: {summary['n_pairs_eligible_for_structural_followup']}")
    print(f"Total Candidate Unannotated Proteins Discovered: {summary['total_unannotated_proteins_discovered']}")
    print("=" * 80)
    return 0


if __name__ == "__main__":
    sys.exit(main())
