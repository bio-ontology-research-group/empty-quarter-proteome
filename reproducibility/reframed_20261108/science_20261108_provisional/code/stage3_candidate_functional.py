#!/usr/bin/env python3
"""Stage 3: Targeted Candidate Family Functional Evaluation, Specific-Context Nulls, and Honest Benchmark/GSPA Reporting.

Implements all requirements from:
- `coordinated-revision-20261108/stage3-scientific-review.txt` (Items 1-6)
- `coordinated-revision-20261108/stage3-negative-outcome-review.txt`
- `coordinated-revision-20261108/stage1b-final-review.txt` (Items 1-5)
- `coordinated-revision-20261108/scientific-design.txt`
- `coordinated-revision-20261108/compute-followup-design.txt`
- `science-20261108/STAGE1B-CONTRACT.txt`
"""
from __future__ import annotations

import collections
import glob
import gzip
import json
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_science import (
    PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES,
    clean_candidate_aa_sequence,
    normalize_genus_groups,
    parse_and_verify_refined_fasta_bytes,
    require_file,
    scale_plddt_explicit_source,
    sha256_file,
    validate_primary_key,
    within_contig_circular_null,
)
from stage1_audit_and_export import (
    SPANY_CLASSES,
    load_all_sp_tm,
    load_c6_thermo_gene_hits,
    parse_gtdb_lineage,
)
from stage3_helpers import (
    BROAD_REPEAT_KEYWORDS,
    CAZY_KEYWORDS,
    CHAPERONE_REPAIR_KEYWORDS,
    NON_DOMAIN_INTERPRO_ANALYSES,
    PEPTIDASE_KEYWORDS,
    SLAYER_KEYWORDS,
    STAGE1B_VERIFIED_PROPOSALS_REQUIRED_COLS,
    TRANSPORTER_KEYWORDS,
    audit_gene_table_syntax_and_fasta,
    build_gspa_pilot_3genome_audit_table,
    build_manifest_genome_alias_resolver,
    calibrate_leave_order_out_threshold,
    classify_family_control_or_candidate_role,
    compute_alignment_summary_metrics,
    evaluate_reference_motif_conservation,
    evaluate_seven_evidence_axes,
    extract_verified_route_a1_reaction_pathway_context,
    load_curated_swissprot_gapseq_reference_from_root,
    qualify_structural_and_domain_evidence_strict,
    reconcile_and_summarize_order_metabolic_tables,
    run_grouped_contig_circular_nulls,
    select_curated_gapseq_reference_for_family,
    summarize_candidate_selection_outcome,
    verify_member_locus_strict,
)

P_ANNOT = Path(os.environ.get("EQ_ANNOT_ROOT", "/project/k10018/work/eq-annot-20261030"))
P_ARCH = Path(os.environ.get("EQ_ARCH_ROOT", "/scratch/project/k10018/eq-archaeal-function-20261106"))
P_PROT_INT = Path(os.environ.get("EQ_PROT_INT_ROOT", "/scratch/project/k10018/eq-proteome-integration-20261106"))
P_GSPA_PILOT = Path(os.environ.get("EQ_GSPA_PILOT_ROOT", "/scratch/project/k10018/eq-proteome-gspa-pilot-20261105"))
P_GAPSEQ_SEQ = Path(os.environ.get("EQ_GAPSEQ_SEQ_ROOT", "/project/k10018/sw/envs/gapseq/share/gapseq/dat/seq/Archaea"))
P_SIF_MAFFT = Path(os.environ.get("EQ_SIF_MAFFT", "/project/k10018/sw/images/eq-annot/mafft.sif"))
P_SIF_IQTREE = Path(os.environ.get("EQ_SIF_IQTREE", "/project/k10018/sw/images/eq-annot/iqtree.sif"))

P_OUT_ROOT = Path(os.environ.get("EQ_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261108"))
P_IN = Path(os.environ.get("EQ_SCIENCE_IN_ROOT", str(P_OUT_ROOT / "inputs")))
P_RES = P_OUT_ROOT / "results"
P_DOSSIERS = P_RES / "stage3_family_dossiers"


def load_verified_refined_genomes_subset(
    target_genomes: set[str],
) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    """Verify and load refined genome contigs for `target_genomes` against `inputs/genome_manifest.tsv`."""
    manifest_path = require_file(P_IN / "genome_manifest.tsv", "inputs/genome_manifest.tsv")
    man_df = pd.read_csv(manifest_path, sep="\t")
    validate_primary_key(man_df, ["mag"], "genome_manifest.tsv")
    man_lookup = man_df.set_index("mag")

    tar_path = require_file(P_IN / "thermo237_refined_mags.tar.gz", "inputs/thermo237_refined_mags.tar.gz")
    prov_map: dict[str, str] = {}
    contigs_map: dict[str, dict[str, str]] = {}

    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf:
            if not member.isfile():
                continue
            bname = Path(member.name).name
            if not bname.endswith(".fa"):
                continue
            gid = bname[:-3]
            if gid not in target_genomes:
                continue
            fobj = tf.extractfile(member)
            if fobj is None:
                raise ValueError(f"Cannot extract {member.name} from {tar_path}")
            raw_bytes = fobj.read()
            if gid not in man_lookup.index:
                prov_map[gid] = "missing_from_genome_manifest"
                continue
            m_row = man_lookup.loc[gid]
            ver = parse_and_verify_refined_fasta_bytes(
                gid,
                raw_bytes,
                expected_sha256=str(m_row["sha256"]),
                expected_contigs=int(m_row["contigs"]),
                expected_bases=int(m_row["bases"]),
            )
            prov_map[gid] = str(ver["provenance_status"])
            contigs_map[gid] = ver["contigs_dict"]  # type: ignore[assignment]

    return prov_map, contigs_map


def load_raw_archaeal_aa_map(target_genes: set[str] | None = None) -> dict[str, str]:
    """Load raw AA sequences from `all_archaea.faa` without stripping internal '*' stops."""
    faa_path = require_file(P_ARCH / "proteins/all_archaea.faa", "all_archaea.faa")
    seqs: dict[str, str] = {}
    cur_id: str | None = None
    chunks: list[str] = []

    def _commit(gid: str | None, parts: list[str]) -> None:
        if gid is None:
            return
        if target_genes is not None and gid not in target_genes:
            return
        seqs[gid] = "".join(parts).strip()

    with open(faa_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                _commit(cur_id, chunks)
                hdr = line[1:].strip().split()[0]
                cur_id = hdr.split("|", 1)[1] if "|" in hdr else hdr
                chunks = []
            else:
                if target_genes is None or cur_id in target_genes:
                    chunks.append(line.strip())
        _commit(cur_id, chunks)
    return seqs


def load_curated_swissprot_gapseq_reference(rxn_id: str) -> dict[str, Any]:
    """Inspect `gapseq` `Archaea/rev/{rxn_id}.fasta` for a reviewed SwissProt reference sequence."""
    return load_curated_swissprot_gapseq_reference_from_root(rxn_id, P_GAPSEQ_SEQ)


def _load_foldseek_db_polymer_sequences(db_prefix: Path, target_snames: set[str]) -> dict[str, str]:
    """Read exact polymer AA sequences from a Foldseek DB (`{db_prefix}.lookup`, `{db_prefix}.index`, `{db_prefix}`)."""
    lookup_p = Path(f"{db_prefix}.lookup")
    index_p = Path(f"{db_prefix}.index")
    data_p = Path(f"{db_prefix}")
    if not target_snames or not (lookup_p.exists() and index_p.exists() and data_p.exists()):
        return {}

    key_to_sname: dict[str, str] = {}
    with open(lookup_p, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1] in target_snames:
                key_to_sname[parts[0]] = parts[1]
    if not key_to_sname:
        return {}

    spans: list[tuple[str, int, int]] = []
    with open(index_p, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[0] in key_to_sname:
                spans.append((key_to_sname[parts[0]], int(parts[1]), int(parts[2])))

    out: dict[str, str] = {}
    with open(data_p, "rb") as fh:
        for sname, offset, length in spans:
            fh.seek(offset)
            raw = fh.read(length).rstrip(b"\x00\r\n ")
            out[sname] = raw.decode("ascii", errors="replace").strip()
    return out


def load_per_protein_raw_evidence(target_cat_ids: set[str], target_genes: set[str]) -> dict[str, Any]:
    """Load raw per-protein structural, profile, InterProScan, C6, and DeepFRI evidence for `target_cat_ids`.

    Preserves separate model records and Foldseek hits per source (`ColabFold_AF2` vs `ESMFold`)
    and extracts the exact model polymer AA sequence from `c10_structclust/{colabfold,esmfold}/db/all`.
    """
    # 1. ColabFold & ESMFold model metadata (`plddt` on 0_to_1 scale) + polymer AA sequences from Foldseek DB
    struct_models: dict[tuple[str, str], dict[str, Any]] = {}
    for src_name, rel_p, db_rel in (
        ("ColabFold_AF2", "results/c10_structclust/colabfold/members.tsv.gz", "results/c10_structclust/colabfold/db/all"),
        ("ESMFold", "results/c10_structclust/esmfold/members.tsv.gz", "results/c10_structclust/esmfold/db/all"),
    ):
        p = require_file(P_ANNOT / rel_p, rel_p)
        snames_needed: set[str] = set()
        temp_rows: list[tuple[str, int, float, float, str]] = []
        with gzip.open(p, "rt") as fh:
            for line in fh:
                arm, pid, len_s, plddt_s, ptm_s, sname = line.rstrip("\n").split("\t")
                norm_pid = pid.replace("__", "::", 1) if "::" not in pid and "__" in pid else pid
                if norm_pid not in target_cat_ids:
                    continue
                plddt_01 = scale_plddt_explicit_source(float(plddt_s), "0_to_1")
                ptm_01 = scale_plddt_explicit_source(float(ptm_s), "0_to_1")
                m_len = int(len_s) if len_s.isdigit() else -1
                temp_rows.append((norm_pid, m_len, plddt_01, ptm_01, sname))
                snames_needed.add(sname)

        poly_map = _load_foldseek_db_polymer_sequences(P_ANNOT / db_rel, snames_needed)
        for norm_pid, m_len, plddt_01, ptm_01, sname in temp_rows:
            poly_seq = poly_map.get(sname, "")
            key = (norm_pid, src_name)
            prev = struct_models.get(key)
            if prev is None or plddt_01 > prev["plddt_0to1"]:
                struct_models[key] = {
                    "model_source": src_name,
                    "model_length": len(poly_seq) if poly_seq else m_len,
                    "model_polymer_aa_seq": poly_seq,
                    "plddt_0to1": plddt_01,
                    "ptm_0to1": ptm_01,
                    "plddt_source_scale": "0_to_1",
                    "struct_cluster": sname,
                }

    # 2. Foldseek PDB hits (`qtmscore`, `qcov`) keyed by `(pid, hit_source)` so AF2 and ESMFold are never cross-paired
    pdb_hits: dict[tuple[str, str], list[dict[str, Any]]] = collections.defaultdict(list)
    pdb_files = [
        ("ColabFold_AF2", str(P_ANNOT / "results/c9cf_colabfold/pdb_whole.tsv.gz")),
        ("ColabFold_AF2", str(P_ANNOT / "results/c9cf_colabfold/pdb_regions.tsv.gz")),
    ] + [("ESMFold", s) for s in sorted(glob.glob(str(P_ANNOT / "results/c9_foldseek/c9pdb_*/chunk_*.pdb_*.tsv.gz")))]
    for hit_src, p_str in pdb_files:
        p = Path(p_str)
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                vals = line.rstrip("\n").split("\t")
                if not vals:
                    continue
                raw_pid = vals[0]
                pid = raw_pid.replace("__", "::", 1) if "::" not in raw_pid and "__" in raw_pid else raw_pid
                if pid not in target_cat_ids:
                    continue
                row = dict(zip(hdr, vals))
                try:
                    pdb_hits[(pid, hit_src)].append(
                        {
                            "hit_source": hit_src,
                            "target": row.get("target", "missing"),
                            "region": row.get("region", "whole"),
                            "qtmscore": float(row["qtmscore"]),
                            "ttmscore": float(row.get("ttmscore", "nan")),
                            "qcov": float(row.get("qcov", "nan")),
                            "tcov": float(row.get("tcov", "nan")),
                            "evalue": float(row.get("evalue", "nan")),
                        }
                    )
                except ValueError:
                    continue

    # 3. Foldseek AFDB hits keyed by `(pid, hit_source)` (has evalue, qcov; lacks qtmscore)
    afdb_hits: dict[tuple[str, str], list[dict[str, Any]]] = collections.defaultdict(list)
    afdb_files = [
        ("ColabFold_AF2", str(P_ANNOT / "results/c9cf_colabfold/afdb_whole.tsv.gz")),
    ] + [("ESMFold", s) for s in sorted(glob.glob(str(P_ANNOT / "results/c9_foldseek/c9afdb_*/chunk_*.afdb_whole.tsv.gz")))]
    for hit_src, p_str in afdb_files:
        p = Path(p_str)
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                vals = line.rstrip("\n").split("\t")
                if not vals:
                    continue
                raw_pid = vals[0]
                pid = raw_pid.replace("__", "::", 1) if "::" not in raw_pid and "__" in raw_pid else raw_pid
                if pid not in target_cat_ids:
                    continue
                row = dict(zip(hdr, vals))
                try:
                    ev = float(row["evalue"])
                    q_s = float(row.get("qstart", "nan"))
                    q_e = float(row.get("qend", "nan"))
                    q_l = float(row.get("qlen", "nan"))
                    qcov = (abs(q_e - q_s) + 1.0) / q_l if (np.isfinite(q_s) and np.isfinite(q_e) and q_l > 0) else np.nan
                    afdb_hits[(pid, hit_src)].append(
                        {
                            "hit_source": hit_src,
                            "target": row.get("target", "missing"),
                            "evalue": ev,
                            "qcov": qcov,
                            "qtmscore": np.nan,
                        }
                    )
                except ValueError:
                    continue

    # 4. HHsearch per-query hits (`c4_per_query.tsv.gz` + `c4_hhsuite/dark/task_*.tophits.tsv.gz` for exact qcov)
    hh_hits: dict[str, dict[str, Any]] = {}
    c4_pq = require_file(P_ANNOT / "results/_agent_analysis/c4_c5_c7/c4/c4_per_query.tsv.gz", "c4_per_query.tsv.gz")
    with gzip.open(c4_pq, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            q = row.get("query", "")
            if q not in target_cat_ids:
                continue
            hh_hits[q] = {
                "pfam_id": row.get("pfam_id") or "missing",
                "pfam_prob": float(row["pfam_prob"]) if row.get("pfam_prob") not in ("", "NA", None) else np.nan,
                "pfam_desc": row.get("pfam_desc") or "missing",
                "pdb_id": row.get("pdb_id") or "missing",
                "pdb_prob": float(row["pdb_prob"]) if row.get("pdb_prob") not in ("", "NA", None) else np.nan,
                "pdb_desc": row.get("pdb_desc") or "missing",
                "qcov": np.nan,
            }

    for p_str in sorted(glob.glob(str(P_ANNOT / "results/c4_hhsuite/*/task_*.tophits.tsv.gz"))):
        with gzip.open(p_str, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                db_name, q_raw, qlen_s, _, t_acc, prob_s, _, _, q_range, _, _, desc = parts[:12]
                q = q_raw.replace("__", "::", 1) if "::" not in q_raw and "__" in q_raw else q_raw
                if q not in target_cat_ids or db_name != "pfam":
                    continue
                try:
                    qlen = float(qlen_s)
                    q_s, q_e = [float(x) for x in q_range.split("-")]
                    qcov = (abs(q_e - q_s) + 1.0) / qlen if qlen > 0 else np.nan
                    prob = float(prob_s)
                except ValueError:
                    continue
                if q not in hh_hits or pd.isna(hh_hits[q].get("qcov")) or prob >= float(hh_hits[q].get("pfam_prob", 0.0)):
                    hh_hits[q] = {
                        "pfam_id": t_acc,
                        "pfam_prob": prob,
                        "pfam_desc": desc,
                        "pdb_id": hh_hits.get(q, {}).get("pdb_id", "missing"),
                        "pdb_prob": hh_hits.get(q, {}).get("pdb_prob", np.nan),
                        "pdb_desc": hh_hits.get(q, {}).get("pdb_desc", "missing"),
                        "qcov": qcov,
                    }

    # 5. InterProScan per-protein intervals (`c5_interpro/{families,proteomes}/chunk_*.tsv.gz`)
    ipr_domains: dict[str, list[str]] = collections.defaultdict(list)
    ipr_nondomain: dict[str, list[str]] = collections.defaultdict(list)
    ipr_intervals: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    for p_str in sorted(glob.glob(str(P_ANNOT / "results/c5_interpro/*/chunk_*.tsv.gz"))):
        with gzip.open(p_str, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 11:
                    continue
                pid = parts[0]
                if pid not in target_cat_ids:
                    continue
                analysis = parts[3]
                sig_acc = parts[4]
                sig_desc = parts[5]
                start_s, end_s = parts[6], parts[7]
                ipr_acc = parts[11] if len(parts) > 11 else "-"
                ipr_desc = parts[12] if len(parts) > 12 else "-"
                tag = f"{analysis}:{sig_acc}({sig_desc if sig_desc != '-' else ipr_desc})[{start_s}-{end_s}]"
                if analysis in NON_DOMAIN_INTERPRO_ANALYSES:
                    if tag not in ipr_nondomain[pid]:
                        ipr_nondomain[pid].append(tag)
                else:
                    if tag not in ipr_domains[pid]:
                        ipr_domains[pid].append(tag)
                    try:
                        ipr_intervals[pid].append(
                            {
                                "analysis": analysis,
                                "sig_acc": sig_acc,
                                "sig_desc": sig_desc,
                                "ipr_acc": ipr_acc,
                                "start": int(start_s),
                                "end": int(end_s),
                                "seq_len": int(parts[2]) if parts[2].isdigit() else -1,
                            }
                        )
                    except ValueError:
                        pass

    # 6. C6 CAZy / MEROPS / TCDB per-gene hits
    c6_hits: dict[str, list[str]] = collections.defaultdict(list)
    for tag, fname in (
        ("TCDB", "tcdb_thermo.tsv.gz"),
        ("MEROPS", "merops_thermo.tsv.gz"),
        ("CAZy", "dbcan_thermo_overview.tsv.gz"),
    ):
        p = P_ANNOT / "results/c6_cazy_merops_tcdb" / fname
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if not parts:
                    continue
                q_raw = parts[0]
                q_gene = q_raw.split("::", 1)[1] if "::" in q_raw else (q_raw.split("|", 1)[1] if "|" in q_raw else q_raw)
                if q_raw in target_cat_ids or q_gene in target_genes:
                    hit_str = f"{tag}:{parts[1]}" if len(parts) > 1 else tag
                    c6_hits[q_gene].append(hit_str)

    # 7. DeepFRI per-protein GO predictions (prob >= 0.5)
    deepfri_hits: dict[str, list[str]] = collections.defaultdict(list)
    for p_str in sorted(glob.glob(str(P_ANNOT / "results/c9_deepfri/*/chunk_*.preds.tsv.gz"))):
        with gzip.open(p_str, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5 or parts[0] == "Protein":
                    continue
                raw_pid, ont, go_id, score_s, go_name = parts[:5]
                pid = raw_pid.replace("__", "::", 1) if "::" not in raw_pid and "__" in raw_pid else raw_pid
                if pid not in target_cat_ids:
                    continue
                try:
                    sc = float(score_s)
                except ValueError:
                    continue
                if sc >= 0.50:
                    deepfri_hits[pid].append(f"{ont}:{go_id}({go_name}:{sc:.2f})")

    return {
        "struct_models": struct_models,
        "pdb_hits": pdb_hits,
        "afdb_hits": afdb_hits,
        "hh_hits": hh_hits,
        "ipr_domains": ipr_domains,
        "ipr_nondomain": ipr_nondomain,
        "ipr_intervals": ipr_intervals,
        "c6_hits": c6_hits,
        "deepfri_hits": deepfri_hits,
    }


def run_mafft_and_iqtree_for_family(
    family_id: str,
    member_seqs: dict[str, str],
    curated_ref: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run MAFFT (`mafft.sif`) and IQ-TREE 3 (`iqtree.sif`, `LG+G4`) on non-redundant verified sequences."""
    P_DOSSIERS.mkdir(parents=True, exist_ok=True)
    safe_fid = family_id.replace("::", "__").replace("/", "_")
    faa_path = P_DOSSIERS / f"{safe_fid}.faa"
    aln_path = P_DOSSIERS / f"{safe_fid}.aln.faa"
    tree_prefix = P_DOSSIERS / f"{safe_fid}.iqtree"

    # Deduplicate exact AA sequences while keeping representative and reference distinct for alignment metrics
    all_seqs: dict[str, str] = {}
    ref_key: str | None = None
    if curated_ref and curated_ref.get("has_reviewed_ref") and curated_ref.get("sequence"):
        ref_key = f"REF__{curated_ref['accession']}"
        all_seqs[ref_key] = str(curated_ref["sequence"]).strip()
    for mid, mseq in member_seqs.items():
        if mseq:
            all_seqs[mid] = mseq

    if len(all_seqs) == 0:
        return {
            "msa_status": "no_verified_member_sequences",
            "n_aligned_sequences": 0,
            "alignment_cols": 0,
            "median_ungapped_coverage": np.nan,
            "min_pairwise_aa_identity": np.nan,
            "median_pairwise_aa_identity": np.nan,
            "max_pairwise_aa_identity": np.nan,
            "tree_status": "no_sequences",
            "tree_file": "none",
            "aligned_sequences": {},
            "reference_id": None,
        }

    with open(faa_path, "wt") as fh:
        for sid, sseq in all_seqs.items():
            fh.write(f">{sid}\n{sseq}\n")

    aligned_map: dict[str, str] = {}
    if len(all_seqs) == 1:
        aligned_map = dict(all_seqs)
        with open(aln_path, "wt") as fh:
            for sid, sseq in aligned_map.items():
                fh.write(f">{sid}\n{sseq}\n")
        msa_status = "single_sequence_no_msa_needed"
    elif P_SIF_MAFFT.exists():
        cmd = [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch,/dev/shm",
            str(P_SIF_MAFFT),
            "mafft",
            "--quiet",
            "--auto",
            "--thread",
            "8",
            str(faa_path),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if proc.returncode == 0 and proc.stdout.strip():
            with open(aln_path, "wt") as fh:
                fh.write(proc.stdout)
            cur_k: str | None = None
            chunks: list[str] = []
            for ln in proc.stdout.splitlines():
                if ln.startswith(">"):
                    if cur_k is not None:
                        aligned_map[cur_k] = "".join(chunks)
                    cur_k = ln[1:].strip().split()[0]
                    chunks = []
                else:
                    chunks.append(ln.strip())
            if cur_k is not None:
                aligned_map[cur_k] = "".join(chunks)
            msa_status = "mafft_completed"
        else:
            msa_status = f"mafft_failed_rc{proc.returncode}"
    else:
        msa_status = "mafft_sif_unavailable"

    # Compute alignment summary metrics on member sequences (excluding external reference from intra-family PID)
    member_aln = {k: v for k, v in aligned_map.items() if k != ref_key}
    metrics = compute_alignment_summary_metrics(member_aln if member_aln else aligned_map)

    # Deduplicate identical aligned sequences before running IQ-TREE (requires >= 4 distinct sequences)
    unique_by_seq: dict[str, str] = {}
    for sid, aseq in aligned_map.items():
        if aseq not in unique_by_seq.values():
            unique_by_seq[sid] = aseq

    tree_status = "fewer_than_4_nonredundant_sequences_tree_skipped"
    tree_file_str = "none"
    if len(unique_by_seq) >= 4 and P_SIF_IQTREE.exists() and aln_path.exists():
        nr_aln_path = P_DOSSIERS / f"{safe_fid}.nr.aln.faa"
        with open(nr_aln_path, "wt") as fh:
            for sid, aseq in unique_by_seq.items():
                fh.write(f">{sid}\n{aseq}\n")
        cmd_iq = [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch,/dev/shm",
            str(P_SIF_IQTREE),
            "iqtree3",
            "-s",
            str(nr_aln_path),
            "-m",
            "LG+G4",
            "-fast",
            "-nt",
            "4",
            "-redo",
            "-pre",
            str(tree_prefix),
        ]
        proc_iq = subprocess.run(cmd_iq, capture_output=True, text=True, check=False)
        tf = Path(f"{tree_prefix}.treefile")
        if proc_iq.returncode == 0 and tf.exists():
            tree_status = f"iqtree3_LG_G4_completed(n_nr={len(unique_by_seq)})"
            tree_file_str = str(tf.relative_to(P_OUT_ROOT)) if tf.is_relative_to(P_OUT_ROOT) else str(tf)
        else:
            tree_status = f"iqtree3_failed_or_constant_alignment(n_nr={len(unique_by_seq)})"

    return {
        "msa_status": msa_status,
        **metrics,
        "tree_status": tree_status,
        "tree_file": tree_file_str,
        "aligned_sequences": aligned_map,
        "reference_id": ref_key,
    }


def classify_flanking_gene_guild(desc: str, profiles: list[str]) -> set[str]:
    """Classify a flanking gene's functional guild(s) for specific-context vs generic known-gene testing."""
    guilds: set[str] = set()
    d_low = str(desc or "").lower()
    p_low = " ".join(profiles).lower()
    combined = f"{d_low} {p_low}".strip()
    if not combined or combined in ("missing", "nan", "hypothetical protein", "uncharacterized protein"):
        return guilds

    guilds.add("generic_any_known_annotation")
    if any(k in combined for k in PEPTIDASE_KEYWORDS):
        guilds.add("extracellular_peptidase_or_proteolysis")
    if any(k in combined for k in CAZY_KEYWORDS):
        guilds.add("carbohydrate_active_or_glycosidase")
    if any(k in combined for k in TRANSPORTER_KEYWORDS):
        guilds.add("solute_or_abc_transporter")
    if any(k in combined for k in CHAPERONE_REPAIR_KEYWORDS):
        guilds.add("heat_shock_chaperone_or_dna_repair")
    if any(k in combined for k in SLAYER_KEYWORDS):
        guilds.add("cell_envelope_or_s_layer_biogenesis")
    return guilds


def run_leave_order_out_benchmark_evaluation(tax_map: dict[str, str], genus_map: dict[str, str]) -> pd.DataFrame:
    """Evaluate `masked_context_benchmark.tsv` with honest leave-order-out threshold calibration."""
    bench_path = require_file(
        P_ARCH / "results/masked_context_benchmark.tsv",
        "eq-archaeal-function-20261106/results/masked_context_benchmark.tsv",
    )
    bm = pd.read_csv(bench_path, sep="\t")
    req_cols = {"genome", "pathway", "target_gene", "anchors", "candidates", "context_score", "top1", "top5", "within_contig_null_top5"}
    missing = sorted(req_cols - set(bm.columns))
    if missing:
        raise ValueError(f"masked_context_benchmark.tsv missing columns: {missing}")

    bm["tax_order"] = bm["genome"].map(lambda g: tax_map.get(str(g), "unclassified_order"))
    bm["genus_block"] = bm["genome"].map(lambda g: genus_map.get(str(g), f"UNCLASSIFIED_GENUS_SINGLETON::{g}"))
    bm["context_score_num"] = pd.to_numeric(bm["context_score"], errors="coerce").fillna(0.0)
    bm["top1_num"] = bm["top1"].astype(bool).astype("float64")
    bm["top5_num"] = bm["top5"].astype(bool).astype("float64")
    bm["null_top5_num"] = pd.to_numeric(bm["within_contig_null_top5"], errors="coerce")
    bm["candidates_num"] = pd.to_numeric(bm["candidates"], errors="coerce")

    orders = sorted(bm["tax_order"].unique())
    rows: list[dict[str, Any]] = []

    for ord_name in orders:
        eval_mask = bm["tax_order"] == ord_name
        train_df = bm.loc[~eval_mask]
        eval_df = bm.loc[eval_mask]
        if eval_df.empty or train_df.empty:
            continue
        cal = calibrate_leave_order_out_threshold(
            train_df["context_score_num"].to_numpy(),
            train_df["top1_num"].to_numpy(),
            eval_df["context_score_num"].to_numpy(),
            eval_df["top1_num"].to_numpy(),
            documented_target_top1_recovery=0.50,
            min_train_passing=20,
            is_held_out_split=True,
        )
        rows.append(
            {
                "held_out_order": ord_name,
                "n_eval_masked_instances": int(len(eval_df)),
                "n_eval_genomes": int(eval_df["genome"].nunique()),
                "n_eval_genus_blocks": int(eval_df["genus_block"].nunique()),
                "unfiltered_top1_recovery": float(eval_df["top1_num"].mean()),
                "unfiltered_top5_recovery": float(eval_df["top5_num"].mean()),
                "within_contig_null_top5_mean": float(eval_df["null_top5_num"].mean()),
                "top5_enrichment_over_within_contig_null": (
                    float(eval_df["top5_num"].mean() / eval_df["null_top5_num"].mean())
                    if float(eval_df["null_top5_num"].mean()) > 0
                    else np.nan
                ),
                "median_window_candidates": float(eval_df["candidates_num"].median()),
                **cal,
            }
        )

    # Pooled in-sample descriptive row (explicitly marked is_held_out_evaluation=False)
    cal_pool = calibrate_leave_order_out_threshold(
        bm["context_score_num"].to_numpy(),
        bm["top1_num"].to_numpy(),
        bm["context_score_num"].to_numpy(),
        bm["top1_num"].to_numpy(),
        documented_target_top1_recovery=0.50,
        min_train_passing=20,
        is_held_out_split=False,
    )
    rows.append(
        {
            "held_out_order": "ALL_POOLED_IN_SAMPLE",
            "n_eval_masked_instances": int(len(bm)),
            "n_eval_genomes": int(bm["genome"].nunique()),
            "n_eval_genus_blocks": int(bm["genus_block"].nunique()),
            "unfiltered_top1_recovery": float(bm["top1_num"].mean()),
            "unfiltered_top5_recovery": float(bm["top5_num"].mean()),
            "within_contig_null_top5_mean": float(bm["null_top5_num"].mean()),
            "top5_enrichment_over_within_contig_null": (
                float(bm["top5_num"].mean() / bm["null_top5_num"].mean())
                if float(bm["null_top5_num"].mean()) > 0
                else np.nan
            ),
            "median_window_candidates": float(bm["candidates_num"].median()),
            **cal_pool,
        }
    )

    return pd.DataFrame(rows)


def build_order_stratified_metabolic_and_gspa_tables(
    meta_df: pd.DataFrame,
    proposals_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build order-stratified metabolic/gap-filling summary and 3-genome GSPA pilot descriptive audit."""
    mod_sum_path = require_file(P_ARCH / "results/metabolic_model_summary.tsv", "metabolic_model_summary.tsv")
    c12_sum_path = require_file(P_ANNOT / "results/c11_gapseq/tables/c12_summary.tsv.gz", "c12_summary.tsv.gz")
    c18_win_path = require_file(P_ANNOT / "results/c18_pathway_holes/windows.tsv", "c18_pathway_holes/windows.tsv")
    manifest_path = require_file(P_IN / "genome_manifest.tsv", "inputs/genome_manifest.tsv")

    mod_df = pd.read_csv(mod_sum_path, sep="\t")
    c12_df = pd.read_csv(c12_sum_path, sep="\t")
    c18_df = pd.read_csv(c18_win_path, sep="\t")
    man_df = pd.read_csv(manifest_path, sep="\t")
    rep_path = P_IN / "refined_representatives.tsv"
    rep_df = pd.read_csv(rep_path, sep="\t") if rep_path.exists() else None

    thermo_meta = meta_df[meta_df["is_thermoplasmatota"].astype(bool)].copy()
    alias_map = build_manifest_genome_alias_resolver(
        thermo_meta["genome"].astype(str).tolist(),
        manifest_df=man_df,
        canonical_rep_df=rep_df,
    )

    order_df = reconcile_and_summarize_order_metabolic_tables(
        thermo_meta,
        mod_df,
        c12_df,
        c18_df,
        proposals_df,
        alias_map,
        primary_medium="native",
    )

    # 3-genome GSPA pilot descriptive audit (never extrapolated to all 237 genomes)
    gspa_sum_path = require_file(P_GSPA_PILOT / "results/summary.tsv", "eq-proteome-gspa-pilot-20261105/results/summary.tsv")
    gspa_ctx_path = require_file(P_GSPA_PILOT / "results/context-changes.tsv", "eq-proteome-gspa-pilot-20261105/results/context-changes.tsv")
    gspa_sum_df = pd.read_csv(gspa_sum_path, sep="\t")
    gspa_ctx_df = pd.read_csv(gspa_ctx_path, sep="\t")

    pilot_df = build_gspa_pilot_3genome_audit_table(
        gspa_sum_df,
        gspa_ctx_df,
        thermo_meta,
        alias_map,
        gspa_gap_pilot_root=P_ARCH / "results/gspa_gap_pilot",
    )
    return order_df, pilot_df


def refresh_stage3_metabolic_and_gspa_summaries_only() -> None:
    """Safely recompute and write the order-stratified metabolic summary and 3-genome GSPA pilot audit
    from authoritative Stage 1b/3 inputs without overwriting completed family dossiers or null evaluations.
    """
    P_RES.mkdir(parents=True, exist_ok=True)
    prop_ver_path = require_file(
        P_RES / "task_a_stage1b_route1_gapfill_proposals_verified.tsv.gz",
        "Stage 1b verified Route A1 proposals",
    )
    meta_path = require_file(P_ARCH / "inputs/thermoplasmatota.tsv", "thermoplasmatota.tsv")
    prop_df = pd.read_csv(prop_ver_path, sep="\t")
    meta_df = pd.read_csv(meta_path, sep="\t")
    tax_parsed = [parse_gtdb_lineage(x) for x in meta_df["gtdb_classification"]]
    for rk, col_name in [("c", "tax_class"), ("o", "tax_order"), ("f", "tax_family"), ("g", "tax_genus")]:
        meta_df[col_name] = [d[rk] for d in tax_parsed]
    g_blocks, is_unclass = normalize_genus_groups(meta_df["genome"].tolist(), meta_df["tax_genus"].tolist())
    meta_df["genus_singleton_block"] = g_blocks
    meta_df["is_unclassified_genus"] = is_unclass
    meta_df["is_thermoplasmatota"] = True

    order_met_df, pilot_gspa_df = build_order_stratified_metabolic_and_gspa_tables(meta_df, prop_df)
    order_met_df.to_csv(P_RES / "task_c_order_stratified_metabolic_and_gspa_descriptive_summary.tsv", sep="\t", index=False)
    pilot_gspa_df.to_csv(P_RES / "task_c_gspa_pilot_3genome_descriptive_audit.tsv", sep="\t", index=False)
    print(
        f"[Stage 3 Summary Refresh] Wrote {len(order_met_df)} order rows and {len(pilot_gspa_df)} GSPA pilot rows to {P_RES}"
    )


def main() -> None:
    P_RES.mkdir(parents=True, exist_ok=True)
    P_DOSSIERS.mkdir(parents=True, exist_ok=True)

    # Pre-flight check for singularity binary and container images before long loads
    if shutil.which("singularity") is None:
        raise RuntimeError("singularity binary not found in PATH; load singularity module before running Stage 3")
    require_file(P_SIF_MAFFT, "mafft.sif container")
    require_file(P_SIF_IQTREE, "iqtree.sif container")

    # 1. Require Stage 1b verified outputs (never fall back to Stage 1 provisional promotion Booleans)
    r1_ver_path = require_file(
        P_RES / "task_a_stage1b_route1_gapfill_family_ranking_verified.tsv",
        "Stage 1b verified Route A1 family ranking",
    )
    r2_ver_path = require_file(
        P_RES / "task_a_stage1b_route2_exported_dark_family_ranking_verified.tsv",
        "Stage 1b verified Route A2 family ranking",
    )
    mem_ver_path = require_file(
        P_RES / "task_a_stage1b_all_candidate_member_loci_verified.tsv.gz",
        "Stage 1b all-member locus verification table",
    )
    prop_ver_path = require_file(
        P_RES / "task_a_stage1b_route1_gapfill_proposals_verified.tsv.gz",
        "Stage 1b verified Route A1 proposals",
    )
    neg_ctrl_path = require_file(
        P_RES / "task_a_route2_annotation_negative_exported_controls.tsv",
        "Route A2 matched annotation-negative exported controls",
    )

    r1_df = pd.read_csv(r1_ver_path, sep="\t")
    r2_df = pd.read_csv(r2_ver_path, sep="\t")
    prop_df = pd.read_csv(prop_ver_path, sep="\t")
    neg_ctrl_df = pd.read_csv(neg_ctrl_path, sep="\t")

    # 2. Load Thermoplasmatota genome metadata, archaeal_genes.tsv, family memberships, and SignalP/TMbed predictions
    meta_path = require_file(P_ARCH / "inputs/thermoplasmatota.tsv", "thermoplasmatota.tsv")
    genes_path = require_file(P_ARCH / "inputs/archaeal_genes.tsv", "archaeal_genes.tsv")
    prof_hits_path = require_file(P_PROT_INT / "results/profile_gene_hits.tsv.gz", "profile_gene_hits.tsv.gz")
    rnacode_path = require_file(P_ANNOT / "results/c14_rnacode/rnacode_long_all.tsv.gz", "rnacode_long_all.tsv.gz")
    pnps_path = require_file(P_ANNOT / "results/c16_pnps/agg/family_pnps.tsv", "family_pnps.tsv")

    meta_df = pd.read_csv(meta_path, sep="\t")
    tax_parsed = [parse_gtdb_lineage(x) for x in meta_df["gtdb_classification"]]
    for rk, col_name in [("c", "tax_class"), ("o", "tax_order"), ("f", "tax_family"), ("g", "tax_genus")]:
        meta_df[col_name] = [d[rk] for d in tax_parsed]
    g_blocks, is_unclass = normalize_genus_groups(meta_df["genome"].tolist(), meta_df["tax_genus"].tolist())
    meta_df["genus_singleton_block"] = g_blocks
    meta_df["is_unclassified_genus"] = is_unclass
    meta_df["is_thermoplasmatota"] = True

    tax_order_map = dict(zip(meta_df["genome"].astype(str), meta_df["tax_order"].astype(str)))
    tax_genus_map = dict(zip(meta_df["genome"].astype(str), meta_df["tax_genus"].astype(str)))
    genus_block_map = dict(zip(meta_df["genome"].astype(str), meta_df["genus_singleton_block"].astype(str)))
    thermo_set = set(meta_df["genome"].astype(str))

    # 3. Select 0 <= n <= 10 new dark candidate families + controls + prior heat-stress examples
    sel_outcome = summarize_candidate_selection_outcome(r1_df, r2_df, max_route1_dark=5, max_route2_dark=5)
    sel_r1_ids = sel_outcome["selected_route1_family_ids"]
    sel_r2_ids = sel_outcome["selected_route2_family_ids"]

    # Record complete disposition of all screened dark families across Route A1 and Route A2
    disp_rows: list[dict[str, Any]] = []
    for r in r1_df.itertuples(index=False):
        fid = str(r.family_id)
        is_dk = bool(getattr(r, "inherited_catalogue_dark", True))
        elig = bool(r.final_candidate_promotion_eligible)
        disp_rows.append(
            {
                "selection_route": "Route_A1_pathway_gapfill",
                "family_id": fid,
                "inherited_catalogue_dark": is_dk,
                "rep_final_promotion_gate_pass": bool(r.rep_final_promotion_gate_pass),
                "n_total_family_member_loci": int(r.n_total_family_member_loci),
                "n_verified_complete_interior_member_loci": int(r.n_verified_complete_interior_member_loci),
                "n_verified_independent_genus_blocks": int(r.n_verified_independent_genus_blocks),
                "final_candidate_promotion_eligible": elig,
                "stage3_disposition": (
                    "selected_new_dark_candidate_route_A1"
                    if fid in sel_r1_ids
                    else ("eligible_dark_not_in_top5_cap" if (is_dk and elig) else ("non_dark_screened" if not is_dk else "rejected_failed_strict_gate_or_recurrence"))
                ),
            }
        )
    for r in r2_df.itertuples(index=False):
        fid = str(r.family_id)
        elig = bool(r.final_candidate_promotion_eligible)
        disp_rows.append(
            {
                "selection_route": "Route_A2_exported_dark",
                "family_id": fid,
                "inherited_catalogue_dark": True,
                "rep_final_promotion_gate_pass": bool(r.rep_final_promotion_gate_pass),
                "n_total_family_member_loci": int(r.n_total_family_member_loci),
                "n_verified_complete_interior_member_loci": int(r.n_verified_complete_interior_member_loci),
                "n_verified_independent_genus_blocks": int(r.n_verified_independent_genus_blocks),
                "final_candidate_promotion_eligible": elig,
                "stage3_disposition": (
                    "selected_new_dark_candidate_route_A2"
                    if fid in sel_r2_ids
                    else ("eligible_dark_not_in_top5_cap" if elig else "rejected_failed_strict_gate_or_recurrence")
                ),
            }
        )
    pd.DataFrame(disp_rows).to_csv(P_RES / "task_c_stage3_screened_family_disposition.tsv", sep="\t", index=False)

    # Select 2 non-dark Route A1 comparator/control families (passing verified complete interior gate)
    r1_nondark = r1_df[
        (~r1_df["inherited_catalogue_dark"].astype(bool))
        & r1_df["rep_final_promotion_gate_pass"].astype(bool)
        & (r1_df["n_verified_independent_genus_blocks"] >= 2)
    ].head(2)
    nondark_ctrl_ids = [str(x) for x in r1_nondark["family_id"].tolist()]

    # Prior heat-stress candidate examples (3 families)
    prior_ids = [fid for fid in sorted(PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES)]

    # Matched annotation-negative exported controls from Route A2 (up to 5 distinct families)
    neg_ctrl_col = "family_id" if "family_id" in neg_ctrl_df.columns else "control_family_id"
    neg_ctrl_ids: list[str] = []
    for fid in neg_ctrl_df[neg_ctrl_col].astype(str).unique():
        if fid not in set(sel_r1_ids + sel_r2_ids + prior_ids + nondark_ctrl_ids):
            neg_ctrl_ids.append(fid)
        if len(neg_ctrl_ids) >= 5:
            break

    target_families_ordered: list[tuple[str, str, bool, bool]] = []
    for fid in sel_r1_ids:
        target_families_ordered.append((fid, "Route_A1_pathway_gapfill_dark", True, False))
    for fid in sel_r2_ids:
        target_families_ordered.append((fid, "Route_A2_exported_dark", True, False))
    for fid in prior_ids:
        if fid not in {x[0] for x in target_families_ordered}:
            target_families_ordered.append((fid, "prior_heat_stress_candidate_example", True, False))
    for fid in nondark_ctrl_ids:
        if fid not in {x[0] for x in target_families_ordered}:
            target_families_ordered.append((fid, "Route_A1_non_dark_comparator", False, False))
    for fid in neg_ctrl_ids:
        if fid not in {x[0] for x in target_families_ordered}:
            target_families_ordered.append((fid, "Route_A2_matched_negative_exported_control", False, True))

    target_fam_set = {x[0] for x in target_families_ordered}

    # 4. Load gene table & family memberships for all Thermoplasmatota genomes to audit all member loci and +/-5 gene neighbourhoods
    genes_df = pd.read_csv(genes_path, sep="\t", low_memory=False)
    genes_thermo = genes_df[genes_df["genome"].astype(str).isin(thermo_set)].copy().reset_index(drop=True)

    mem2fam: dict[str, str] = {}
    for fn in ("member_fam_final.tsv", "ctrl_member_fam.tsv"):
        with open(require_file(P_ANNOT / "inputs/fam" / fn, fn)) as fh:
            for line in fh:
                m, r = line.rstrip("\n").split("\t")[:2]
                mem2fam[m] = r

    def _resolve_fam(c_id: str, s_rep: str, c_rep: str) -> str:
        for k in (c_id, s_rep, c_rep):
            if isinstance(k, str) and k in mem2fam:
                return mem2fam[k]
        if isinstance(s_rep, str) and s_rep and s_rep != "nan":
            return s_rep
        if isinstance(c_rep, str) and c_rep and c_rep != "nan":
            return c_rep
        return "unassigned"

    genes_thermo["family_id"] = [
        _resolve_fam(str(c), str(s), str(cr))
        for c, s, cr in zip(genes_thermo["cat_id"], genes_thermo["seq_cluster_rep"], genes_thermo["cat_cluster_rep"])
    ]
    genes_thermo["is_dark"] = pd.to_numeric(genes_thermo["dark"], errors="coerce").fillna(1.0) == 1.0

    # Load SignalP 6.0h & TMbed predictions and C6 hits for Thermoplasmatota genes
    tm_all, sp_all = load_all_sp_tm(n_workers=32)
    cazy_gene, mer_gene, tc_gene = load_c6_thermo_gene_hits(thermo_set, n_workers=16)

    # Load all_archaea.faa for all Thermoplasmatota genes in target families and their contigs
    target_member_rows = genes_thermo[genes_thermo["family_id"].astype(str).isin(target_fam_set)].copy()
    target_contigs_by_genome = set(
        zip(
            target_member_rows["genome"].astype(str),
            target_member_rows["contig"].astype(str),
        )
    )
    contig_genes_df = genes_thermo[
        list(zip(genes_thermo["genome"].astype(str), genes_thermo["contig"].astype(str)))
    ].copy() if False else genes_thermo[
        [
            (str(g), str(c)) in target_contigs_by_genome
            for g, c in zip(genes_thermo["genome"], genes_thermo["contig"])
        ]
    ].copy()

    # Attach SignalP & TMbed columns onto contig_genes_df and target_member_rows
    def _attach_sp_tm(df_in: pd.DataFrame) -> pd.DataFrame:
        sp_cls: list[str] = []
        sp_any_l: list[bool] = []
        tm_s_l: list[bool] = []
        tm_h_l: list[int] = []
        tm_b_l: list[int] = []
        for g in df_in["gene"].astype(str):
            sp_e = sp_all.get(("T", g)) or sp_all.get(("D", g))
            tm_e = tm_all.get(("T", g)) or tm_all.get(("D", g))
            sc = sp_e[0] if sp_e is not None else "missing"
            sp_cls.append(sc)
            sp_any_l.append(sc in SPANY_CLASSES)
            tm_s_l.append(bool(tm_e[1]) if tm_e is not None else False)
            tm_h_l.append(int(tm_e[2]) if tm_e is not None else 0)
            tm_b_l.append(int(tm_e[3]) if tm_e is not None else 0)
        df_out = df_in.copy()
        df_out["signalp_class"] = sp_cls
        df_out["sp_any"] = sp_any_l
        df_out["tmbed_S"] = tm_s_l
        df_out["tmbed_H"] = tm_h_l
        df_out["tmbed_B"] = tm_b_l
        return df_out

    target_member_rows = _attach_sp_tm(target_member_rows)
    contig_genes_df = _attach_sp_tm(contig_genes_df)

    raw_aa_map = load_raw_archaeal_aa_map(set(contig_genes_df["gene"].astype(str)))
    audited_contig_genes = audit_gene_table_syntax_and_fasta(contig_genes_df, raw_aa_map)
    audited_gene_lookup = audited_contig_genes.set_index("gene")

    # Verify authoritative refined genome FASTAs against genome_manifest.tsv for all genomes containing target family members
    member_genomes_set = set(target_member_rows["genome"].astype(str))
    prov_status_map, refined_contigs_map = load_verified_refined_genomes_subset(member_genomes_set)

    # Independently verify EVERY member locus across all evaluated families against the refined genome contigs
    member_manifest_rows: list[dict[str, Any]] = []
    target_cat_ids: set[str] = set()
    target_gene_ids: set[str] = set()

    for fam_id, sel_route, is_dk, is_neg_ctrl in target_families_ordered:
        fam_annot = target_member_rows[target_member_rows["family_id"].astype(str) == fam_id]
        for a_row in fam_annot.itertuples(index=False):
            gid = str(a_row.gene)
            g_mag = str(a_row.genome)
            cat_id = str(a_row.cat_id) if pd.notna(a_row.cat_id) and str(a_row.cat_id).strip() else f"{g_mag.split('__', 1)[0]}::{gid}"
            target_cat_ids.add(cat_id)
            target_gene_ids.add(gid)

            if gid not in audited_gene_lookup.index:
                continue
            g_info = audited_gene_lookup.loc[gid]
            prov_st = prov_status_map.get(g_mag, "missing_from_refined_archive")
            c_dict = refined_contigs_map.get(g_mag, {})
            raw_seq = raw_aa_map.get(gid)
            edge_st = str(g_info["contig_edge_status"])

            vres = verify_member_locus_strict(
                prov_st,
                c_dict,
                str(g_info["contig"]),
                g_info["start"],
                g_info["end"],
                g_info["strand"],
                g_info["partial"],
                raw_seq,
                edge_st,
            )
            verified_interior = bool(vres["member_final_promotion_gate_pass"] and edge_st == "interior")
            member_manifest_rows.append(
                {
                    "family_id": fam_id,
                    "selection_route": sel_route,
                    "gene": gid,
                    "cat_id": cat_id,
                    "genome": g_mag,
                    "tax_order": tax_order_map.get(g_mag, "unknown"),
                    "tax_genus": tax_genus_map.get(g_mag, "unknown"),
                    "genus_singleton_block": genus_block_map.get(g_mag, f"UNCLASSIFIED_GENUS_SINGLETON::{g_mag}"),
                    "contig": str(g_info["contig"]),
                    "start": g_info["start"],
                    "end": g_info["end"],
                    "strand": g_info["strand"],
                    "partial": str(g_info["partial"]),
                    "length": int(g_info["length"]),
                    "contig_rank": int(g_info["contig_rank"]),
                    "contig_n_genes": int(g_info["contig_n_genes"]),
                    "contig_edge_status": edge_st,
                    "fasta_sequence_audit_status": str(g_info["fasta_sequence_audit_status"]),
                    "genome_manifest_provenance_status": vres["genome_manifest_provenance_status"],
                    "refined_contig_present": bool(vres["refined_contig_present"]),
                    "locus_within_contig_bounds": bool(vres["locus_within_contig_bounds"]),
                    "has_admissible_start_codon": bool(vres["has_admissible_start_codon"]),
                    "has_terminal_stop_codon": bool(vres["has_terminal_stop_codon"]),
                    "cds_translation_exact_match": bool(vres["cds_translation_exact_match"]),
                    "n_internal_stops_translated": int(vres["n_internal_stops_translated"]),
                    "n_internal_stops_in_candidate": int(vres["n_internal_stops_in_candidate"]),
                    "verification_status": str(vres["verification_status"]),
                    "member_final_promotion_gate_pass": bool(vres["member_final_promotion_gate_pass"]),
                    "verified_complete_interior_member": verified_interior,
                    "signalp_class": str(a_row.signalp_class) if pd.notna(a_row.signalp_class) else "missing",
                    "sp_any": bool(a_row.sp_any) if pd.notna(a_row.sp_any) else False,
                    "tmbed_S": bool(a_row.tmbed_S) if pd.notna(a_row.tmbed_S) else False,
                    "tmbed_H": int(a_row.tmbed_H) if pd.notna(a_row.tmbed_H) else 0,
                    "tmbed_B": int(a_row.tmbed_B) if pd.notna(a_row.tmbed_B) else 0,
                    "audited_aa_sequence": str(g_info["audited_aa_sequence"]) if verified_interior else "",
                }
            )

    member_manifest_df = pd.DataFrame(member_manifest_rows)
    member_manifest_out = member_manifest_df.drop(columns=["audited_aa_sequence"], errors="ignore")
    member_manifest_out.to_csv(P_RES / "task_c_stage3_member_locus_manifest_verified.tsv", sep="\t", index=False)

    # 5. Load raw per-protein structural, profile, InterProScan, C6, DeepFRI, RNAcode, and pN/pS evidence
    raw_ev = load_per_protein_raw_evidence(target_cat_ids, target_gene_ids)

    rnacode_df = pd.read_csv(rnacode_path, sep="\t")
    rnacode_map = rnacode_df.set_index("family").to_dict(orient="index")
    pnps_df = pd.read_csv(pnps_path, sep="\t")
    pnps_map = pnps_df.set_index("family").to_dict(orient="index")

    # Build lookup for flanking gene annotations (C6, HHsearch c4_per_query, and microTrait/stress profile hits) on target contigs
    prof_df = pd.read_csv(prof_hits_path, sep="\t")
    prof_sub = prof_df[prof_df["gene"].isin(set(audited_contig_genes["gene"]))].copy()
    gene_profiles_map: dict[str, list[str]] = collections.defaultdict(list)
    for pr in prof_sub.itertuples(index=False):
        gene_profiles_map[str(pr.gene)].append(str(pr.profile))

    flank_fams_needed = set(audited_contig_genes["family_id"].astype(str)) | set(audited_contig_genes["cat_id"].astype(str))
    flank_hh_desc: dict[str, str] = {}
    with gzip.open(require_file(P_ANNOT / "results/_agent_analysis/c4_c5_c7/c4/c4_per_query.tsv.gz", "c4_per_query.tsv.gz"), "rt") as fh:
        hdr_c4 = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r_c4 = dict(zip(hdr_c4, line.rstrip("\n").split("\t")))
            q_c4 = r_c4.get("query", "").replace("__", "::", 1)
            if q_c4 in flank_fams_needed:
                bp = r_c4.get("best_pfam") or ""
                pd_d = r_c4.get("pdb_desc") or ""
                if bp or pd_d:
                    flank_hh_desc[q_c4] = f"Pfam:{bp} {pd_d}".strip()

    flank_annot_map: dict[str, dict[str, Any]] = {}
    for cg_row in audited_contig_genes.itertuples(index=False):
        fg = str(cg_row.gene)
        ffid = str(cg_row.family_id)
        fcid = str(cg_row.cat_id)
        is_dk_fg = bool(cg_row.is_dark)
        desc_parts: list[str] = []
        if fg in tc_gene:
            desc_parts.append(f"TC:{tc_gene[fg]}")
        if fg in mer_gene:
            desc_parts.append(f"MER:{mer_gene[fg]}")
        if fg in cazy_gene:
            desc_parts.append(f"CAZy:{cazy_gene[fg]}")
        hh_d = flank_hh_desc.get(ffid) or flank_hh_desc.get(fcid)
        if hh_d:
            desc_parts.append(hh_d)
        if gene_profiles_map.get(fg):
            desc_parts.append("profiles:" + ",".join(gene_profiles_map[fg]))
        if not desc_parts and not is_dk_fg:
            desc_parts.append("annotated_catalogue_protein")
        flank_annot_map[fg] = {
            "family_id": ffid,
            "is_dark": is_dk_fg,
            "signalp_class": str(cg_row.signalp_class),
            "functional_description": "; ".join(desc_parts) if desc_parts else "hypothetical protein",
        }

    # Build per-contig ordered gene arrays for +/-5 gene neighbourhood panels and circular permutation nulls
    contig_ordered_genes: dict[tuple[str, str], list[str]] = {}
    for (g_mag, c_id), grp in audited_contig_genes.sort_values("contig_rank").groupby(["genome", "contig"], sort=False):
        contig_ordered_genes[(str(g_mag), str(c_id))] = grp["gene"].astype(str).tolist()

    # Validate Stage 1b verified Route A1 proposals schema and map verified priority reactions & pathways per family
    r1_context_map = extract_verified_route_a1_reaction_pathway_context(prop_df)

    # 6. Evaluate each family across all 7 evidence axes
    eval_rows: list[dict[str, Any]] = []
    neighbourhood_rows: list[dict[str, Any]] = []
    competing_rows: list[dict[str, Any]] = []

    for fam_id, sel_route, is_dk, is_neg_ctrl in target_families_ordered:
        fam_mems = member_manifest_df[member_manifest_df["family_id"] == fam_id].copy()
        if fam_mems.empty:
            continue

        # Identify representative locus: prefer an actually verified complete interior member (matching fam_id if verified, else longest verified interior member)
        ver_mems = fam_mems[fam_mems["verified_complete_interior_member"].astype(bool)]
        exact_ver_rep = ver_mems[ver_mems["cat_id"] == fam_id]
        if not exact_ver_rep.empty:
            rep_row = exact_ver_rep.iloc[0]
        elif not ver_mems.empty:
            rep_row = ver_mems.sort_values(["length", "gene"], ascending=[False, True]).iloc[0]
        else:
            exact_any_rep = fam_mems[fam_mems["cat_id"] == fam_id]
            rep_row = exact_any_rep.iloc[0] if not exact_any_rep.empty else fam_mems.iloc[0]

        rep_gene = str(rep_row["gene"])
        rep_cat_id = str(rep_row["cat_id"])
        rep_genome = str(rep_row["genome"])
        rep_aa = raw_aa_map.get(rep_gene)  # NEVER fabricate a fallback AA sequence!

        # Assemble per-member structural & domain records strictly keyed to each member's own cat_id and model source
        mem_ev_records: list[dict[str, Any]] = []
        for m in fam_mems.itertuples(index=False):
            m_gene = str(m.gene)
            m_cid = str(m.cat_id)
            m_raw_aa = raw_aa_map.get(m_gene)
            m_clean = clean_candidate_aa_sequence(m_raw_aa)
            m_clean_aa = str(m_clean["aa_body"]) if m_clean["status"] == "ok" else ""
            h_hit = raw_ev["hh_hits"].get(m_cid, {})
            i_doms = raw_ev["ipr_domains"].get(m_cid, [])
            i_non = raw_ev["ipr_nondomain"].get(m_cid, [])

            emitted_any_source = False
            for src_name in ("ColabFold_AF2", "ESMFold"):
                s_mod = raw_ev["struct_models"].get((m_cid, src_name), {})
                p_hits = raw_ev["pdb_hits"].get((m_cid, src_name), [])
                a_hits = raw_ev["afdb_hits"].get((m_cid, src_name), [])
                if not s_mod and not p_hits and not a_hits and emitted_any_source:
                    continue
                best_pdb = max(p_hits, key=lambda x: x["qtmscore"]) if p_hits else {}
                best_afdb = min(a_hits, key=lambda x: x["evalue"]) if a_hits else {}
                poly_seq = s_mod.get("model_polymer_aa_seq")

                mem_ev_records.append(
                    {
                        "locus_gene": m_gene,
                        "locus_cat_id": m_cid,
                        "locus_aa_seq": m_raw_aa,
                        "claim_as_representative_own": (m_gene == rep_gene),
                        "model_source": s_mod.get("model_source", ""),
                        "model_polymer_aa_seq": poly_seq,
                        "plddt_raw": s_mod.get("plddt_0to1", np.nan),
                        "plddt_source_scale": s_mod.get("plddt_source_scale", "0_to_1"),
                        "pdb_hit_source": best_pdb.get("hit_source", src_name if best_pdb else ""),
                        "pdb_target": best_pdb.get("target", "missing"),
                        "pdb_qtmscore": best_pdb.get("qtmscore", np.nan),
                        "pdb_qcov": best_pdb.get("qcov", np.nan),
                        "afdb_hit_source": best_afdb.get("hit_source", src_name if best_afdb else ""),
                        "afdb_target": best_afdb.get("target", "missing"),
                        "afdb_evalue": best_afdb.get("evalue", np.nan),
                        "afdb_qcov": best_afdb.get("qcov", np.nan),
                        "afdb_qtmscore": best_afdb.get("qtmscore", np.nan),
                        "hhsearch_pfam": h_hit.get("pfam_id", "missing") if not emitted_any_source else "missing",
                        "hhsearch_prob": h_hit.get("pfam_prob", np.nan) if not emitted_any_source else np.nan,
                        "hhsearch_qcov": h_hit.get("qcov", np.nan) if not emitted_any_source else np.nan,
                        "interpro_domains": i_doms if not emitted_any_source else [],
                        "interpro_non_domain_flags": i_non if not emitted_any_source else [],
                    }
                )
                emitted_any_source = True

                # Record competing structural assignments per model source
                poly_ok = bool(poly_seq and m_clean_aa and poly_seq == m_clean_aa)
                for ph in sorted(p_hits, key=lambda x: -x["qtmscore"])[:3]:
                    competing_rows.append(
                        {
                            "family_id": fam_id,
                            "locus_gene": m_gene,
                            "locus_cat_id": m_cid,
                            "is_representative_locus": bool(m_gene == rep_gene),
                            "evidence_source": f"Foldseek_PDB_{src_name}",
                            "target_or_signature": str(ph["target"]),
                            "score_metric": f"qtm={ph['qtmscore']:.3f},ttm={ph['ttmscore']:.3f},poly_match={poly_ok}",
                            "query_coverage": float(ph["qcov"]),
                            "qualified_pass": bool(
                                poly_ok
                                and s_mod.get("plddt_0to1", 0.0) >= 0.70
                                and len(m_clean_aa) >= 100
                                and ph["qtmscore"] >= 0.50
                                and pd.notna(ph["qcov"])
                                and ph["qcov"] >= 0.50
                            ),
                        }
                    )

            for idom in i_doms:
                competing_rows.append(
                    {
                        "family_id": fam_id,
                        "locus_gene": m_gene,
                        "locus_cat_id": m_cid,
                        "is_representative_locus": bool(m_gene == rep_gene),
                        "evidence_source": "InterProScan_domain",
                        "target_or_signature": idom,
                        "score_metric": "domain_interval",
                        "query_coverage": np.nan,
                        "qualified_pass": bool(m_clean_aa),
                    }
                )
            for inon in i_non:
                competing_rows.append(
                    {
                        "family_id": fam_id,
                        "locus_gene": m_gene,
                        "locus_cat_id": m_cid,
                        "is_representative_locus": bool(m_gene == rep_gene),
                        "evidence_source": "InterProScan_non_domain_excluded",
                        "target_or_signature": inon,
                        "score_metric": "disorder_or_coil_not_domain",
                        "query_coverage": np.nan,
                        "qualified_pass": False,
                    }
                )

        qual_ev = qualify_structural_and_domain_evidence_strict(
            rep_gene,
            rep_aa,
            mem_ev_records,
            allow_member_summary_separate=True,
            forbid_mixed_locus_in_rep=True,
        )

        # Check curated SwissProt reference in gapseq Archaea/rev/ if Route A1 reactions are linked
        fam_r1_ctx = r1_context_map.get(
            fam_id,
            {"missing_reactions": [], "pathways": [], "reaction_names": [], "ecs": []},
        )
        curated_ref = select_curated_gapseq_reference_for_family(fam_id, r1_context_map, P_GAPSEQ_SEQ)

        # Run MAFFT and IQ-TREE on verified complete interior member AA sequences (+ curated reference if available)
        ver_interior_mems = fam_mems[fam_mems["verified_complete_interior_member"].astype(bool)].copy()
        ver_seq_map = {
            str(r.gene): str(r.audited_aa_sequence)
            for r in ver_interior_mems.itertuples(index=False)
            if str(r.audited_aa_sequence)
        }
        msa_res = run_mafft_and_iqtree_for_family(
            fam_id,
            ver_seq_map,
            curated_ref if curated_ref.get("has_reviewed_ref") else None,
        )

        # Compute query-to-curated-reference coverage if curated reference was aligned
        q_to_ref_cov = np.nan
        if msa_res["reference_id"] and msa_res["reference_id"] in msa_res["aligned_sequences"] and rep_gene in msa_res["aligned_sequences"]:
            ref_aln_s = msa_res["aligned_sequences"][msa_res["reference_id"]]
            rep_aln_s = msa_res["aligned_sequences"][rep_gene]
            ref_ungapped = sum(1 for c in ref_aln_s if c != "-")
            shared_cols = sum(1 for c1, c2 in zip(ref_aln_s, rep_aln_s) if c1 != "-" and c2 != "-")
            q_to_ref_cov = (shared_cols / float(ref_ungapped)) if ref_ungapped > 0 else np.nan

        # Classify family role strictly (forced MSA coverage + generic domain + reaction reference CANNOT establish positive control)
        role_info = classify_family_control_or_candidate_role(
            fam_id,
            is_dark=is_dk,
            selection_route=sel_route,
            is_matched_negative_exported_control=is_neg_ctrl,
            claim_annotated_positive_control=(not is_dk and not is_neg_ctrl),
            curated_reference_accession=curated_ref.get("accession") if curated_ref.get("has_reviewed_ref") else None,
            curated_function_description=curated_ref.get("description") if curated_ref.get("has_reviewed_ref") else None,
            curated_evidence_source=curated_ref.get("source") if curated_ref.get("has_reviewed_ref") else None,
            query_domain_signatures=raw_ev["ipr_domains"].get(rep_cat_id, []),
            reference_domain_signatures=None,
            independently_verified_curated_homology=False,
            query_to_curated_ref_cov=q_to_ref_cov,
            strict_reject_unverified_positive_control=False,
        )

        # Evaluate reference motif conservation (never inventing a motif without a documented reference active site)
        motif_res = evaluate_reference_motif_conservation(
            msa_res["aligned_sequences"],
            reference_id=msa_res["reference_id"],
            documented_catalytic_sites=None,  # SwissProt FASTA headers have no catalytic residue coordinates
        )

        # Build +/-5 gene neighbourhood panels and run within-contig circular permutation nulls
        # Separate SPECIFIC guild/pathway hypothesis from generic known-gene proximity
        guild_counts: collections.Counter[str] = collections.Counter()
        member_contigs_for_null: list[dict[str, Any]] = []
        guild_genus_blocks: dict[str, set[str]] = collections.defaultdict(set)

        for m in ver_interior_mems.itertuples(index=False):
            m_gene = str(m.gene)
            m_mag = str(m.genome)
            m_contig = str(m.contig)
            m_gblock = str(m.genus_singleton_block)
            c_genes = contig_ordered_genes.get((m_mag, m_contig), [])
            if m_gene not in c_genes:
                continue
            cand_idx = c_genes.index(m_gene)

            # Build per-gene guild vector along this contig
            contig_guilds: list[set[str]] = []
            for idx_g, fg in enumerate(c_genes):
                fa = flank_annot_map.get(fg, {})
                f_desc = str(fa.get("functional_description") or fa.get("interpro_desc") or fa.get("pfam_desc") or "hypothetical protein")
                f_profs = gene_profiles_map.get(fg, [])
                g_set = classify_flanking_gene_guild(f_desc, f_profs)
                contig_guilds.append(g_set)
                if abs(idx_g - cand_idx) <= 5 and idx_g != cand_idx:
                    for g_tag in g_set:
                        if g_tag != "generic_any_known_annotation":
                            guild_counts[g_tag] += 1
                            guild_genus_blocks[g_tag].add(m_gblock)
                    neighbourhood_rows.append(
                        {
                            "family_id": fam_id,
                            "candidate_role": str(role_info["candidate_role"]),
                            "member_gene": m_gene,
                            "genome": m_mag,
                            "tax_order": str(m.tax_order),
                            "genus_singleton_block": m_gblock,
                            "contig": m_contig,
                            "relative_gene_offset": int(idx_g - cand_idx),
                            "flanking_gene": fg,
                            "flanking_family_id": str(fa.get("family_id", "missing")),
                            "flanking_is_dark": bool(fa.get("is_dark", True)),
                            "flanking_signalp_class": str(fa.get("signalp_class", "missing")),
                            "flanking_functional_description": f_desc,
                            "flanking_profiles": ";".join(f_profs) if f_profs else "none",
                            "flanking_specific_guilds": ";".join(sorted(g_set - {"generic_any_known_annotation"})) or "none",
                        }
                    )

            member_contigs_for_null.append(
                {
                    "genome": m_mag,
                    "contig": m_contig,
                    "genus_block": m_gblock,
                    "cand_idx": cand_idx,
                    "contig_guilds": contig_guilds,
                }
            )

        # Determine dominant specific flanking guild (or specific gapfill pathway for Route A1)
        # Explicitly record that post-hoc selection of the winning neighbourhood guild is exploratory & selection-dependent
        if guild_counts:
            top_guild = sorted(guild_counts.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
            guild_selection_mode = "exploratory_posthoc_observed_neighbourhood_winner_selection_dependent"
        elif fam_r1_ctx["pathways"]:
            top_guild = f"pathway_window:{fam_r1_ctx['pathways'][0]}"
            guild_selection_mode = "route_a1_gapfill_pathway_hypothesis_exploratory"
        else:
            top_guild = "none_localisation_only"
            guild_selection_mode = "none_localisation_only"

        # Run within-contig circular permutation null (1,000 permutations) grouped by unique (genome, contig),
        # merging multiple target members on the same contig instead of pseudoreplicating the contig:
        # (a) Specific hypothesised guild (`top_guild`)
        # (b) Generic any-known-gene adjacency (`generic_any_known_annotation`)
        null_bundle = run_grouped_contig_circular_nulls(
            member_contigs_for_null,
            specific_guild=top_guild,
            specific_guild_selection_mode=guild_selection_mode,
            window=5,
            n_perm=1000,
            seed=20261108,
        )
        spec_null = null_bundle["specific_null"]
        gen_null = null_bundle["generic_null"]
        n_spec_genus_blocks = len(guild_genus_blocks.get(top_guild, set()))

        # Recurrence and paralogy strictly from verified complete interior member loci
        n_total_loci = int(len(fam_mems))
        n_ver_loci = int(len(ver_interior_mems))
        n_ver_genomes = int(ver_interior_mems["genome"].nunique()) if n_ver_loci > 0 else 0
        n_ver_orders = int(ver_interior_mems["tax_order"].nunique()) if n_ver_loci > 0 else 0
        n_ver_gblocks = int(ver_interior_mems["genus_singleton_block"].nunique()) if n_ver_loci > 0 else 0
        max_copies = int(ver_interior_mems.groupby("genome")["gene"].count().max()) if n_ver_loci > 0 else 0
        frac_exported = (
            float((ver_interior_mems["sp_any"].astype(bool) | ver_interior_mems["tmbed_S"].astype(bool)).mean())
            if n_ver_loci > 0
            else 0.0
        )

        # Summarise competing / alternative hypotheses
        rep_c6 = ";".join(raw_ev["c6_hits"].get(rep_gene, [])) or "none"
        rep_deepfri = ";".join(raw_ev["deepfri_hits"].get(rep_cat_id, [])[:3]) or "none"
        rc_info = rnacode_map.get(fam_id, {})
        pnps_info = pnps_map.get(fam_id, {})

        competing_parts: list[str] = []
        if qual_ev["rep_own_pdb_target"] != "missing":
            competing_parts.append(
                f"rep_pdb={qual_ev['rep_own_pdb_target']}(qtm={qual_ev['rep_own_pdb_qtmscore']:.2f},qcov={qual_ev['rep_own_pdb_qcov']:.2f})"
            )
        if qual_ev["rep_own_afdb_target"] != "missing":
            competing_parts.append(
                f"rep_afdb={qual_ev['rep_own_afdb_target']}(ev={qual_ev['rep_own_afdb_evalue']:.1e},{qual_ev['rep_own_afdb_support_status']})"
            )
        if rep_c6 != "none":
            competing_parts.append(f"rep_c6={rep_c6}")
        if rep_deepfri != "none":
            competing_parts.append(f"rep_deepfri={rep_deepfri}")
        if not competing_parts:
            competing_parts.append("no_competing_structural_or_domain_hit_on_representative")

        axes = evaluate_seven_evidence_axes(
            candidate_role=str(role_info["candidate_role"]),
            rep_final_promotion_gate_pass=bool(rep_row["verified_complete_interior_member"]),
            n_verified_complete_loci=n_ver_loci,
            n_total_loci=n_total_loci,
            rep_signalp_class=str(rep_row["signalp_class"]),
            rep_tmbed_s=bool(rep_row["tmbed_S"]),
            rep_tm_helices=int(rep_row["tmbed_H"]),
            rep_tm_beta=int(rep_row["tmbed_B"]),
            frac_family_exported=frac_exported,
            rep_has_own_domain_or_profile=bool(qual_ev["rep_own_has_domain_or_profile_support"]),
            rep_own_domain_summary=f"hh={qual_ev['rep_own_hhsearch_pfam']};ipr={qual_ev['rep_own_interpro_domain_signatures']}",
            rep_has_own_qualified_whole_fold=bool(qual_ev["rep_own_has_qualified_whole_fold_support"]),
            rep_own_fold_summary=f"pdb={qual_ev['rep_own_pdb_target']}(qtm={qual_ev['rep_own_pdb_qtmscore']},qcov={qual_ev['rep_own_pdb_qcov']})",
            other_member_evidence_summary=f"struct={qual_ev['other_member_qualified_structural_evidence']};prof={qual_ev['other_member_qualified_profile_evidence']}",
            motif_conservation_status=str(motif_res["status"]),
            median_alignment_coverage=float(msa_res["median_ungapped_coverage"]),
            n_verified_orders=n_ver_orders,
            n_verified_independent_genus_blocks=n_ver_gblocks,
            max_copies_per_genome=max_copies,
            specific_hypothesis_name=top_guild,
            specific_context_obs_adj=int(spec_null["obs_adjacent"]),
            specific_context_p_emp=float(spec_null["p_empirical_one_sided"]),
            n_verified_genus_blocks_with_specific_context=n_spec_genus_blocks,
            competing_assignments_summary=";".join(competing_parts),
            positive_control_evidence_status=str(role_info["positive_control_evidence_status"]),
            specific_guild_selection_mode=guild_selection_mode,
        )

        eval_rows.append(
            {
                "family_id": fam_id,
                "selection_route": sel_route,
                "candidate_role": str(role_info["candidate_role"]),
                "positive_control_evidence_status": str(role_info["positive_control_evidence_status"]),
                "has_curated_positive_control_evidence": bool(role_info["has_curated_positive_control_evidence"]),
                "route_a1_proposed_missing_reactions": ";".join(fam_r1_ctx["missing_reactions"]) or "none",
                "route_a1_proposed_pathways": ";".join(fam_r1_ctx["pathways"]) or "none",
                "route_a1_proposed_reaction_names": ";".join(fam_r1_ctx["reaction_names"]) or "none",
                "route_a1_proposed_ecs": ";".join(fam_r1_ctx["ecs"]) or "none",
                "curated_reference_accession": curated_ref.get("accession", "none") if curated_ref.get("has_reviewed_ref") else "none",
                "curated_reference_description": curated_ref.get("description", "none") if curated_ref.get("has_reviewed_ref") else "none",
                "query_to_curated_ref_coverage": q_to_ref_cov,
                "representative_gene": rep_gene,
                "representative_cat_id": rep_cat_id,
                "representative_genome": rep_genome,
                "rep_verified_complete_interior_cds": bool(rep_row["verified_complete_interior_member"]),
                "n_total_member_loci": n_total_loci,
                "n_verified_complete_interior_loci": n_ver_loci,
                "n_verified_complete_interior_genomes": n_ver_genomes,
                "n_verified_orders": n_ver_orders,
                "n_verified_independent_genus_blocks": n_ver_gblocks,
                "max_verified_copies_per_genome": max_copies,
                "frac_verified_members_exported": frac_exported,
                **qual_ev,
                "rep_own_c6_hits": rep_c6,
                "rep_own_deepfri_hits": rep_deepfri,
                "rnacode_p_best": rc_info.get("p_best", np.nan),
                "rnacode_status": rc_info.get("status", "not_assayed"),
                "family_median_pnps": pnps_info.get("median_pnps", np.nan),
                "msa_status": msa_res["msa_status"],
                "n_aligned_verified_sequences": msa_res["n_aligned_sequences"],
                "alignment_cols": msa_res["alignment_cols"],
                "median_ungapped_alignment_coverage": msa_res["median_ungapped_coverage"],
                "min_pairwise_aa_identity": msa_res["min_pairwise_aa_identity"],
                "median_pairwise_aa_identity": msa_res["median_pairwise_aa_identity"],
                "max_pairwise_aa_identity": msa_res["max_pairwise_aa_identity"],
                "phylogenetic_tree_status": msa_res["tree_status"],
                "phylogenetic_tree_file": msa_res["tree_file"],
                "motif_conservation_status": motif_res["status"],
                "motif_site_audit_details": motif_res["site_audit_details"],
                "specific_flanking_hypothesis_guild": top_guild,
                "specific_flanking_guild_selection_status": guild_selection_mode,
                "specific_context_n_unique_contigs": int(null_bundle["n_unique_contigs"]),
                "specific_context_n_target_loci": int(spec_null["n_target_loci"]),
                "specific_context_observed_adjacent": int(spec_null["obs_adjacent"]),
                "specific_context_null_mean_adjacent": float(spec_null["null_mean"]),
                "specific_context_p_empirical": float(spec_null["p_empirical_one_sided"]),
                "n_verified_genus_blocks_with_specific_context": n_spec_genus_blocks,
                "generic_known_flank_observed_adjacent": int(gen_null["obs_adjacent"]),
                "generic_known_flank_null_mean_adjacent": float(gen_null["null_mean"]),
                "generic_known_flank_p_empirical": float(gen_null["p_empirical_one_sided"]),
                **axes,
            }
        )

    eval_df = pd.DataFrame(eval_rows)
    validate_primary_key(eval_df, ["family_id"], "task_c_targeted_family_functional_evaluation.tsv")
    eval_df.to_csv(P_RES / "task_c_targeted_family_functional_evaluation.tsv", sep="\t", index=False)

    pd.DataFrame(neighbourhood_rows).to_csv(P_RES / "task_c_candidate_gene_neighbourhood_panels.tsv", sep="\t", index=False)
    pd.DataFrame(competing_rows).to_csv(P_RES / "task_c_competing_structural_and_domain_assignments.tsv", sep="\t", index=False)

    # 7. Honest Leave-Order-Out Threshold Calibration on masked_context_benchmark.tsv
    loo_df = run_leave_order_out_benchmark_evaluation(tax_order_map, genus_block_map)
    loo_df.to_csv(P_RES / "task_c_masked_benchmark_leave_order_out_threshold_calibration.tsv", sep="\t", index=False)
    loo_df.to_csv(P_RES / "task_c_masked_benchmark_leave_order_out.tsv", sep="\t", index=False)

    # 8. Order-stratified metabolic/gap-filling summary and 3-genome GSPA pilot descriptive audit
    order_met_df, pilot_gspa_df = build_order_stratified_metabolic_and_gspa_tables(meta_df, prop_df)
    order_met_df.to_csv(P_RES / "task_c_order_stratified_metabolic_and_gspa_descriptive_summary.tsv", sep="\t", index=False)
    pilot_gspa_df.to_csv(P_RES / "task_c_gspa_pilot_3genome_descriptive_audit.tsv", sep="\t", index=False)

    # 9. Write stage3_summary.json with exact input hashes and outcome status
    summary_payload = {
        "stage": "stage3_candidate_functional",
        "status": "completed",
        **sel_outcome,
        "n_evaluated_total_families": int(len(eval_df)),
        "n_evaluated_new_dark_candidates": int((eval_df["candidate_role"] == "de_novo_dark_candidate_family").sum()),
        "n_evaluated_prior_heat_stress_examples": int((eval_df["candidate_role"] == "prior_heat_stress_candidate_example").sum()),
        "n_evaluated_annotated_positive_controls": int((eval_df["candidate_role"] == "annotated_positive_control").sum()),
        "n_evaluated_sequence_annotated_comparators": int((eval_df["candidate_role"] == "sequence_annotated_comparator").sum()),
        "n_evaluated_matched_negative_exported_controls": int((eval_df["candidate_role"] == "matched_annotation_negative_exported_control").sum()),
        "n_verified_member_loci_checked": int(len(member_manifest_df)),
        "n_verified_complete_interior_member_loci": int(member_manifest_df["verified_complete_interior_member"].sum()),
        "input_sha256": {
            "task_a_stage1b_route1_gapfill_family_ranking_verified.tsv": sha256_file(r1_ver_path),
            "task_a_stage1b_route2_exported_dark_family_ranking_verified.tsv": sha256_file(r2_ver_path),
            "task_a_stage1b_all_candidate_member_loci_verified.tsv.gz": sha256_file(mem_ver_path),
            "genome_manifest.tsv": sha256_file(P_IN / "genome_manifest.tsv"),
            "masked_context_benchmark.tsv": sha256_file(P_ARCH / "results/masked_context_benchmark.tsv"),
        },
    }
    (P_RES / "stage3_summary.json").write_text(json.dumps(summary_payload, indent=2) + "\n")
    print("Stage 3 completed successfully:", json.dumps(summary_payload, indent=2))


if __name__ == "__main__":
    if "--summary-only" in sys.argv:
        refresh_stage3_metabolic_and_gspa_summaries_only()
    else:
        main()
