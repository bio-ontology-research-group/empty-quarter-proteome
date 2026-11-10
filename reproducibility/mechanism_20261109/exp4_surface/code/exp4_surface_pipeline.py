#!/usr/bin/env python3
"""Experiment 4 (`exp4_surface`): Exported & Cell-Surface Dark Protein Functional Discrimination Pipeline.

Executes on Shaheen III (`k10018`, `workq`) under `/scratch/project/k10018/eq-mechanism-20261109/exp4_surface`.
Implements all requirements from:
- `mechanism-20261109/EXP4-CONTRACT.txt`
- `mechanism-20261109/COMMON-CONTRACT.txt`
- `mechanism-20261109/EXP4-CURATED-SOURCE-NOTE.txt`
- `mechanism-20261109/EXP4-EVIDENCE-REVIEW.txt`
- `mechanism-20261109/EXP4-CODE-REVIEW.txt` (Items 1-9)
"""
from __future__ import annotations

import argparse
import collections
import glob
import gzip
import hashlib
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import tarfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_curated_reference_panel import PANEL_SPECS, PANEL_VERSION
from common_science import (
    PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES,
    clean_candidate_aa_sequence,
    compute_canonical_top_genome_per_genus_block,
    normalize_genus_groups,
    parse_and_verify_refined_fasta_bytes,
    require_dir,
    require_file,
    scale_plddt_explicit_source,
    sha256_file,
    validate_primary_key,
)
from exp4_surface_helpers import (
    BOUNDED_GUILD_PATTERNS,
    KNOWN_SIGNALP_CLASSES,
    PANEL_SPECS_BY_ID,
    ROLE_CLASSES,
    SPANY_CLASSES,
    apply_bh_fdr,
    classify_export_and_surface_topology,
    classify_flanking_machinery_guilds,
    clean_aa_sequence,
    compute_sequence_integrity_and_compositional_metrics as compute_sequence_composition_and_motifs,
    evaluate_competing_roles_for_family,
    evaluate_gspa_heldout_recovery_and_leakage_guard,
    evaluate_same_role_evidence_and_convergence,
    extract_tmbed_segments_1based,
    qualify_structural_hit_for_member,
    resolve_curated_reference_id_from_header,
    resolve_pdb_chain_and_architecture_guards,
    verify_prediction_sequence_hash_match,
)
from stage3_helpers import (
    NON_DOMAIN_INTERPRO_ANALYSES,
    audit_gene_table_syntax_and_fasta,
    compute_alignment_summary_metrics,
    run_grouped_contig_circular_nulls,
    verify_member_locus_strict,
)

# Authoritative read-only Shaheen input roots
P_ANNOT = Path(os.environ.get("EQ_ANNOT_ROOT", "/project/k10018/work/eq-annot-20261030"))
P_ARCH = Path(os.environ.get("EQ_ARCH_ROOT", "/scratch/project/k10018/eq-archaeal-function-20261106"))
P_PROT_INT = Path(os.environ.get("EQ_PROT_INT_ROOT", "/scratch/project/k10018/eq-proteome-integration-20261106"))
P_CANON_SCI = Path(
    os.environ.get("EQ_CANON_SCI_ROOT", "/scratch/project/k10018/eq-science-20261109-canonical-taxonomy")
)
P_LEG_SCI = Path(os.environ.get("EQ_LEG_SCI_ROOT", "/scratch/project/k10018/eq-science-20261108"))
P_GSPA_PILOT = Path(
    os.environ.get("EQ_GSPA_PILOT_ROOT", "/scratch/project/k10018/eq-proteome-gspa-pilot-20261105")
)

# Containers and binaries on Shaheen
P_SIF_FOLDSEEK = Path(os.environ.get("EQ_SIF_FOLDSEEK", "/project/k10018/sw/images/eq-annot/foldseek.sif"))
P_SIF_DIAMOND = Path(os.environ.get("EQ_SIF_DIAMOND", "/project/k10018/sw/images/eq-annot/diamond.sif"))
P_SIF_MAFFT = Path(os.environ.get("EQ_SIF_MAFFT", "/project/k10018/sw/images/eq-annot/mafft.sif"))
P_SIF_IQTREE = Path(os.environ.get("EQ_SIF_IQTREE", "/project/k10018/sw/images/eq-annot/iqtree.sif"))
P_BIN_DIAMOND = Path(os.environ.get("EQ_BIN_DIAMOND", "/project/k10018/sw/envs/diamond-2.1.11/bin/diamond"))

# Branch output root on Shaheen
P_EXP4_ROOT = Path(
    os.environ.get("EQ_EXP4_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp4_surface")
)
P_IN = P_EXP4_ROOT / "inputs"
P_PROV = P_EXP4_ROOT / "provenance"
P_RES = P_EXP4_ROOT / "results"
P_LOGS = P_EXP4_ROOT / "logs"
P_TMP = P_EXP4_ROOT / "tmp"
P_SHARD_IN = P_RES / "shard_inputs"
P_SHARD_OUT = P_RES / "shards"
P_MSA_OUT = P_RES / "family_msas_and_trees"

N_ARRAY_SHARDS = 16  # Shards 1..16 for full array; Shard 0 is the Pilot shard

STAGE3_FOCAL_DARK_CANDIDATES = (
    "51Dr3::51Dr3_552805_4",
    "58PRr3::58PRr3_78238_7",
    "V54PRr1::V54PRr1_1398418_22",
    "59PRr1::59PRr1_330165_10",
    "T1PRr2::T1PRr2_2909672_25",
    "60PRr1::60PRr1_861859_6",
    "F41Dr3::F41Dr3_256695_46",
    "F22PRr2::F22PRr2_1552766_13",
    "F41Dr3::F41Dr3_256695_58",
    "58PRr3::58PRr3_784227_12",
)

# 27 exploratory dark-family candidates from completed v1 baseline (`EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt` Item C)
EXPLORATORY_DARK_CANDIDATES_V1: set[str] = {
    "19Dr1::19Dr1_436937_5",
    "59PRr1::59PRr1_33956_8",
    "34PRr1::34PRr1_571293_2",
    "F4PRr2::F4PRr2_215129_17",
    "V39PRr2::V39PRr2_1502258_21",
    "S22PRr1::S22PRr1_528275_6",
    "S46Dr1::S46Dr1_193475_151",
    "27Dr3::27Dr3_1657614_16",
    "F4PRr3::F4PRr3_995500_4",
    "V60PRr2::V60PRr2_46430_4",
    "51Dr1::51Dr1_699885_22",
    "60PRr1::60PRr1_285908_17",
    "F22PRr2::F22PRr2_1141396_50",
    "V45Dr1::V45Dr1_339889_128",
    "34PRr3::34PRr3_94396_4",
    "57PRr2::57PRr2_346723_2",
    "F12PRr2::F12PRr2_1087883_2",
    "24Dr1::24Dr1_1510502_4",
    "F26PRr1::F26PRr1_2434251_8",
    "44Dr2::44Dr2_1410473_3",
    "F34PRr2::F34PRr2_511316_37",
    "62PRr2::62PRr2_1114592_4",
    "36Sr2::36Sr2_1322350_1",
    "27Dr3::27Dr3_608182_78",
    "S54PRr1::S54PRr1_594119_11",
    "V45Dr1::V45Dr1_276274_5",
    "44Dr2::44Dr2_2152055_9",
}

# Pre-declared context leads with v1 empirical p <= 0.01 (followup selection, not independent discovery)
PREDECLARED_SMALL_P_FAMILIES_V1: set[str] = {
    "59PRr1::59PRr1_464316_7",
    "V45Dr1::V45Dr1_323118_21",
    "F35Dr1::F35Dr1_990114_1",
    "S25PRr1::S25PRr1_803204_4",
    "F3PRr1::F3PRr1_1124057_24",
    "F38Dr3::F38Dr3_679630_1",
    "50PRr3::50PRr3_907761_1",
    "22PRr3::22PRr3_756560_4",
    "V5PRr1::V5PRr1_2276374_7",
    "V60PRr2::V60PRr2_631919_1",
    "F26PRr1::F26PRr1_482741_1",
    "22PRr1::22PRr1_2531955_3",
}


def parse_gtdb_lineage(tax_str: str) -> dict[str, str]:
    ranks = {"d": "", "p": "", "c": "", "o": "", "f": "", "g": "", "s": ""}
    if not isinstance(tax_str, str):
        return ranks
    for part in tax_str.split(";"):
        part = part.strip()
        if "__" in part:
            k, v = part.split("__", 1)
            if k in ranks:
                ranks[k] = f"{k}__{v}" if v else f"{k}__"
    return ranks


def _norm_sp_tm_header(h: str) -> tuple[str, str]:
    tok = h.split()[0]
    if "|" in tok:
        return ("T", tok.split("|", 1)[1])
    return ("D", tok)


def _count_segments(lab: str, chars: str) -> int:
    n = 0
    prev = False
    for c in lab:
        cur = c in chars
        if cur and not prev:
            n += 1
        prev = cur
    return n


def _parse_tmbed_chunk_with_spans(
    path_str: str,
) -> dict[tuple[str, str], tuple[int, bool, int, int, str, list[tuple[int, int]]]]:
    """Parse TMbed `.pred.gz` chunk preserving exact 1-based `(h_start, h_end)` helix spans."""
    out: dict[tuple[str, str], tuple[int, bool, int, int, str, list[tuple[int, int]]]] = {}
    with gzip.open(path_str, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline().strip()
            lab = fh.readline().strip()
            if not h.startswith(">") or len(seq) != len(lab):
                raise ValueError(f"Malformed TMbed record in {path_str}: {h[:60]}")
            k = _norm_sp_tm_header(h[1:])
            has_s = "S" in lab
            h_spans = extract_tmbed_segments_1based(lab, "Hh")
            n_h = len(h_spans)
            n_b = _count_segments(lab, "Bb")
            topo = f"S:{int(has_s)}_H:{n_h}_B:{n_b}"
            out[k] = (len(seq), has_s, n_h, n_b, topo, h_spans)
    return out


def _parse_signalp_chunk_with_cs(path_str: str) -> dict[tuple[str, str], tuple[str, float, str, int | None]]:
    """Parse SignalP 6.0h prediction results preserving parsed integer cleavage site position when reported."""
    cs_re = re.compile(r"CS pos:\s*(\d+)-(\d+)")
    out: dict[tuple[str, str], tuple[str, float, str, int | None]] = {}
    with gzip.open(path_str, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            k = _norm_sp_tm_header(parts[0])
            pred = parts[1].strip()
            other_prob = float(parts[2])
            sp_prob = max(0.0, 1.0 - other_prob)
            cs_raw = parts[8].strip() if len(parts) > 8 else ""
            m_cs = cs_re.search(cs_raw)
            cs_pos = int(m_cs.group(1)) if m_cs else None
            out[k] = (pred, sp_prob, cs_raw, cs_pos)
    return out


def load_all_sp_tm_with_spans(
    n_workers: int = 32,
) -> tuple[
    dict[tuple[str, str], tuple[int, bool, int, int, str, list[tuple[int, int]]]],
    dict[tuple[str, str], tuple[str, float, str, int | None]],
]:
    tm_files = sorted(glob.glob(str(P_ANNOT / "results/sp_tm_analysis/tmbed/chunk_*.pred.gz")))
    sp_files = sorted(glob.glob(str(P_ANNOT / "results/sp_tm_analysis/signalp/chunk_*.prediction_results.txt.gz")))
    if not tm_files or not sp_files:
        raise FileNotFoundError("Missing SignalP or TMbed chunk files in sp_tm_analysis")

    tm_all: dict[tuple[str, str], tuple[int, bool, int, int, str, list[tuple[int, int]]]] = {}
    sp_all: dict[tuple[str, str], tuple[str, float, str, int | None]] = {}
    with ProcessPoolExecutor(max_workers=min(n_workers, 40)) as ex:
        for d in ex.map(_parse_tmbed_chunk_with_spans, tm_files):
            tm_all.update(d)
        for d in ex.map(_parse_signalp_chunk_with_cs, sp_files):
            sp_all.update(d)
    return tm_all, sp_all


def load_c6_thermo_gene_hits(thermo_genomes: set[str]) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    cazy_gene: dict[str, str] = {}
    mer_gene: dict[str, str] = {}
    tc_gene: dict[str, str] = {}
    for target_dict, fname in (
        (tc_gene, "tcdb_thermo.tsv.gz"),
        (mer_gene, "merops_thermo.tsv.gz"),
        (cazy_gene, "dbcan_thermo_overview.tsv.gz"),
    ):
        p = P_ANNOT / "results/c6_cazy_merops_tcdb" / fname
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                q_raw = parts[0]
                if q_raw.startswith("Gene_ID") or q_raw.startswith("query"):
                    continue
                q_gene = (
                    q_raw.split("::", 1)[1]
                    if "::" in q_raw
                    else (q_raw.split("|", 1)[1] if "|" in q_raw else q_raw)
                )
                if q_gene not in target_dict:
                    target_dict[q_gene] = parts[1]
    return cazy_gene, mer_gene, tc_gene


def load_raw_archaeal_aa_map(target_genes: set[str] | None = None) -> dict[str, str]:
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


def load_verified_refined_genomes_subset(
    target_genomes: set[str],
    sci_in_dir: Path,
) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    manifest_path = require_file(sci_in_dir / "genome_manifest.tsv", "genome_manifest.tsv")
    man_df = pd.read_csv(manifest_path, sep="\t")
    validate_primary_key(man_df, ["mag"], "genome_manifest.tsv")
    man_lookup = man_df.set_index("mag")

    tar_path = require_file(sci_in_dir / "thermo237_refined_mags.tar.gz", "thermo237_refined_mags.tar.gz")
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


def _load_foldseek_db_polymer_sequences(db_prefix: Path, target_snames: set[str]) -> dict[str, str]:
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
    struct_models: dict[tuple[str, str], dict[str, Any]] = {}
    sname_to_pids: dict[tuple[str, str], list[str]] = collections.defaultdict(list)
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
                sname_to_pids[(src_name, sname)].append(norm_pid)

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

    hh_hits: dict[str, dict[str, Any]] = {}
    c4_pq = require_file(P_ANNOT / "results/_agent_analysis/c4_c5_c7/c4/c4_per_query.tsv.gz", "c4_per_query.tsv.gz")
    with gzip.open(c4_pq, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            q_raw = row.get("query", "")
            q = q_raw.replace("__", "::", 1) if "::" not in q_raw and "__" in q_raw else q_raw
            if q not in target_cat_ids:
                continue
            hh_hits[q] = {
                "pfam_id": row.get("pfam_id") or row.get("best_pfam") or "missing",
                "pfam_prob": float(row["pfam_prob"]) if row.get("pfam_prob") not in ("", "NA", None) else np.nan,
                "pfam_desc": row.get("pfam_desc") or "missing",
                "pdb_id": row.get("pdb_id") or "missing",
                "pdb_prob": float(row["pdb_prob"]) if row.get("pdb_prob") not in ("", "NA", None) else np.nan,
                "pdb_desc": row.get("pdb_desc") or "missing",
                "qcov": np.nan,
            }

    ipr_domains: dict[str, list[str]] = collections.defaultdict(list)
    ipr_nondomain: dict[str, list[str]] = collections.defaultdict(list)
    for p_str in sorted(glob.glob(str(P_ANNOT / "results/c5_interpro/families/chunk_*.tsv.gz"))):
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
                ipr_desc = parts[12] if len(parts) > 12 else "-"
                tag = f"{analysis}:{sig_acc}({sig_desc if sig_desc != '-' else ipr_desc})[{start_s}-{end_s}]"
                if analysis in NON_DOMAIN_INTERPRO_ANALYSES:
                    if tag not in ipr_nondomain[pid]:
                        ipr_nondomain[pid].append(tag)
                else:
                    if tag not in ipr_domains[pid]:
                        ipr_domains[pid].append(tag)

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

    return {
        "struct_models": struct_models,
        "sname_to_pids": sname_to_pids,
        "pdb_hits": pdb_hits,
        "hh_hits": hh_hits,
        "ipr_domains": ipr_domains,
        "ipr_nondomain": ipr_nondomain,
        "c6_hits": c6_hits,
    }


def run_curated_panel_diamond_and_foldseek(
    universe_faa: Path,
    sname_to_pids: dict[tuple[str, str], list[str]],
    threads: int = 64,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, list[dict[str, Any]]]]:
    """Run DIAMOND blastp and Foldseek against the v1.1-20261109 curated reference panel,
    including explicit positive self-query controls and confusable negative control verification
    (`EXP4-PIPELINE-REVIEW.txt` Items 1 & 3).
    """
    P_TMP.mkdir(parents=True, exist_ok=True)
    panel_faa = require_file(P_IN / "curated_role_reference_panel.faa", "curated_role_reference_panel.faa")
    pdb_chains_dir = require_dir(P_IN / "pdb_chains", "inputs/pdb_chains directory")

    # Prepend positive self-query controls (`POSCTRL_SELF__{ref_id}`) from curated_role_reference_panel.faa
    combined_query_faa = P_TMP / "diamond_query_with_posctrl.faa"
    expected_dmnd_posctrls: set[str] = set()
    with open(combined_query_faa, "wt") as q_out:
        with open(panel_faa, "rt") as p_in:
            for ln in p_in:
                if ln.startswith(">"):
                    rid = resolve_curated_reference_id_from_header(ln[1:].strip(), require_known=True)
                    expected_dmnd_posctrls.add(rid)
                    q_out.write(f">POSCTRL_SELF__{rid}\n")
                else:
                    q_out.write(ln)
        with open(universe_faa, "rt") as u_in:
            shutil.copyfileobj(u_in, q_out)

    # 1. DIAMOND blastp against curated_role_reference_panel.faa
    dmnd_db = P_TMP / "curated_panel_db"
    dmnd_out = P_RES / "curated_panel_diamond_hits.tsv"
    if P_BIN_DIAMOND.exists():
        dmnd_cmd_prefix = [str(P_BIN_DIAMOND)]
    else:
        require_file(P_SIF_DIAMOND, "diamond.sif")
        dmnd_cmd_prefix = ["singularity", "exec", "-B", "/project,/scratch,/dev/shm", str(P_SIF_DIAMOND), "diamond"]

    subprocess.run(
        dmnd_cmd_prefix + ["makedb", "--in", str(panel_faa), "-d", str(dmnd_db), "--threads", "1", "--quiet"],
        check=True,
    )
    subprocess.run(
        dmnd_cmd_prefix
        + [
            "blastp",
            "-q",
            str(combined_query_faa),
            "-d",
            f"{dmnd_db}.dmnd",
            "-o",
            str(dmnd_out),
            "--ultra-sensitive",
            "--max-target-seqs",
            "25",
            "--evalue",
            "10.0",
            "--threads",
            str(threads),
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen",
            "--quiet",
        ],
        check=True,
    )

    seq_hits_by_gene: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    dmnd_posctrl_best: dict[str, dict[str, Any]] = {}
    if dmnd_out.exists():
        with open(dmnd_out, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 14:
                    continue
                qid, sid, pident, alnlen, _, _, qst, qen, sst, sen, ev, bits, qlen, slen = parts[:14]
                ref_id = resolve_curated_reference_id_from_header(sid, require_known=True)
                spec = PANEL_SPECS_BY_ID[ref_id]
                ql = float(qlen)
                sl = float(slen)
                qc = (abs(float(qen) - float(qst)) + 1.0) / ql if ql > 0 else np.nan
                tc = (abs(float(sen) - float(sst)) + 1.0) / sl if sl > 0 else np.nan
                rec_d = {
                    "query_gene": qid,
                    "ref_id": ref_id,
                    "target_header": sid,
                    "pident": float(pident),
                    "alnlen": int(alnlen),
                    "qcov": qc,
                    "tcov": tc,
                    "evalue": float(ev),
                    "bits": float(bits),
                    "experimentally_supported_roles": spec.get("experimentally_supported_roles", "none"),
                    "hypothesized_roles": spec.get("hypothesized_roles", "none"),
                    "control_type": spec.get("control_type", ""),
                }
                if qid.startswith("POSCTRL_SELF__"):
                    expected_rid = qid.split("POSCTRL_SELF__", 1)[1]
                    if ref_id == expected_rid:
                        prev_b = dmnd_posctrl_best.get(expected_rid)
                        if prev_b is None or float(bits) > float(prev_b["bits"]):
                            dmnd_posctrl_best[expected_rid] = rec_d
                else:
                    seq_hits_by_gene[qid].append(rec_d)

    missing_dmnd_pos = sorted(expected_dmnd_posctrls - set(dmnd_posctrl_best.keys()))
    if missing_dmnd_pos:
        raise RuntimeError(
            f"DIAMOND positive self-query control failed to recover reference IDs: {missing_dmnd_pos}"
        )

    # 2. Foldseek search of ColabFold and ESMFold DBs against single-chain PDB files in inputs/pdb_chains/*.pdb
    require_file(P_SIF_FOLDSEEK, "foldseek.sif")
    fs_tmp = P_TMP / "foldseek_curated"
    if fs_tmp.exists():
        shutil.rmtree(fs_tmp)
    fs_tmp.mkdir(parents=True, exist_ok=True)

    pdb_only_dir = fs_tmp / "pdb_only"
    pdb_only_dir.mkdir(parents=True, exist_ok=True)
    expected_fs_posctrls: set[str] = set()
    for p_pdb in sorted(pdb_chains_dir.glob("REF_*.pdb")):
        shutil.copy2(p_pdb, pdb_only_dir / p_pdb.name)
        expected_fs_posctrls.add(p_pdb.stem)

    target_db = fs_tmp / "curated_pdb_db"
    subprocess.run(
        [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch,/dev/shm",
            str(P_SIF_FOLDSEEK),
            "foldseek",
            "createdb",
            str(pdb_only_dir),
            str(target_db),
            "--threads",
            str(min(threads, 16)),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )

    # Parse target_db.lookup and polymer sequences for unambiguous Foldseek target resolution
    target_lookup_map: dict[str, str] = {}
    t_lookup_path = Path(f"{target_db}.lookup")
    if t_lookup_path.exists():
        with open(t_lookup_path, "rt") as fh:
            for ln in fh:
                pts = ln.rstrip("\n").split("\t")
                if len(pts) >= 2:
                    target_lookup_map[pts[0]] = pts[1]
                    target_lookup_map[pts[1]] = pts[1]
    target_poly_map = _load_foldseek_db_polymer_sequences(target_db, set(target_lookup_map.values()))

    # Run Foldseek positive self-control (`target_db` vs `target_db`) and verify all 18 PDB references
    self_aln = fs_tmp / "aln_self_posctrl"
    self_m8 = fs_tmp / "fs_self_posctrl.m8"
    self_work = fs_tmp / "work_self_posctrl"
    self_work.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch,/dev/shm",
            str(P_SIF_FOLDSEEK),
            "foldseek",
            "search",
            str(target_db),
            str(target_db),
            str(self_aln),
            str(self_work),
            "--alignment-type",
            "1",
            "-a",
            "-e",
            "10.0",
            "--threads",
            str(min(threads, 16)),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    subprocess.run(
        [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch,/dev/shm",
            str(P_SIF_FOLDSEEK),
            "foldseek",
            "convertalis",
            str(target_db),
            str(target_db),
            str(self_aln),
            str(self_m8),
            "--threads",
            str(min(threads, 16)),
            "--format-output",
            "query,target,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,qtmscore,ttmscore,alntmscore,evalue,bits",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    fs_posctrl_best: dict[str, dict[str, Any]] = {}
    if self_m8.exists():
        with open(self_m8, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 17:
                    continue
                q_raw, t_raw, fident, alnlen, _, _, _, _, _, _, qcov, tcov, qtm, ttm, _, ev, bits = parts[:17]
                q_rid = resolve_curated_reference_id_from_header(
                    q_raw,
                    target_polymer_seq=target_poly_map.get(q_raw),
                    db_lookup_map=target_lookup_map,
                    require_known=True,
                )
                t_rid = resolve_curated_reference_id_from_header(
                    t_raw,
                    target_polymer_seq=target_poly_map.get(t_raw),
                    db_lookup_map=target_lookup_map,
                    require_known=True,
                )
                if q_rid == t_rid:
                    prev_f = fs_posctrl_best.get(q_rid)
                    if prev_f is None or float(qtm) > float(prev_f["qtmscore"]):
                        fs_posctrl_best[q_rid] = {
                            "ref_id": q_rid,
                            "query_raw": q_raw,
                            "target_raw": t_raw,
                            "fident": float(fident),
                            "qcov": float(qcov),
                            "tcov": float(tcov),
                            "qtmscore": float(qtm),
                            "ttmscore": float(ttm),
                            "evalue": float(ev),
                            "bits": float(bits),
                        }
    missing_fs_pos = sorted(expected_fs_posctrls - set(fs_posctrl_best.keys()))
    if missing_fs_pos:
        raise RuntimeError(
            f"Foldseek positive self-query control failed to recover PDB chain reference IDs: {missing_fs_pos}"
        )

    # Record positive self-control and confusable negative control audit table
    ctrl_audit_rows: list[dict[str, Any]] = []
    for rid in sorted(expected_dmnd_posctrls):
        d_hit = dmnd_posctrl_best[rid]
        spec = PANEL_SPECS_BY_ID[rid]
        pass_gate = bool(
            d_hit["pident"] >= 99.0
            and d_hit["qcov"] >= 0.85
            and d_hit["tcov"] >= 0.85
            and d_hit["evalue"] <= 1e-10
            and bool(spec.get("ref_id"))
        )
        if not pass_gate:
            raise RuntimeError(f"DIAMOND positive control failed coverage/identity gate for {rid}: {d_hit}")
        ctrl_audit_rows.append(
            {
                "control_category": "positive_self_control_diamond",
                "ref_id": rid,
                "raw_target_header": d_hit["target_header"],
                "resolved_ref_id": d_hit["ref_id"],
                "metric_primary": d_hit["pident"],
                "qcov": d_hit["qcov"],
                "tcov": d_hit["tcov"],
                "experimentally_supported_roles": d_hit["experimentally_supported_roles"],
                "gate_passed": pass_gate,
                "audit_note": "Verified exact `|`-delimited header mapping and PANEL_SPECS_BY_ID lookup",
            }
        )
    for rid in sorted(expected_fs_posctrls):
        f_hit = fs_posctrl_best[rid]
        spec = PANEL_SPECS_BY_ID[rid]
        pass_gate = bool(
            f_hit["qtmscore"] >= 0.98
            and f_hit["qcov"] >= 0.98
            and f_hit["tcov"] >= 0.98
            and bool(spec.get("ref_id"))
        )
        if not pass_gate:
            raise RuntimeError(f"Foldseek positive control failed TM/coverage gate for {rid}: {f_hit}")
        ctrl_audit_rows.append(
            {
                "control_category": "positive_self_control_foldseek",
                "ref_id": rid,
                "raw_target_header": f_hit["target_raw"],
                "resolved_ref_id": f_hit["ref_id"],
                "metric_primary": f_hit["qtmscore"],
                "qcov": f_hit["qcov"],
                "tcov": f_hit["tcov"],
                "experimentally_supported_roles": spec.get("experimentally_supported_roles", "none"),
                "gate_passed": pass_gate,
                "audit_note": "Verified Foldseek target lookup + chain suffix mapping and PANEL_SPECS_BY_ID lookup",
            }
        )
    for neg_rid in (
        "REF_CTRL_9UES_ChainA_Propeptide_NegativeCtrl",
        "REF_CTRL_8CE8_Entity5_AuthE_LabelI_CcmE",
        "REF_CTRL_8CE8_Entity1_LabelE_Autha_CcmA",
        "REF_SBP_SSO3043_SugarConfusableCtrl",
    ):
        spec = PANEL_SPECS_BY_ID[neg_rid]
        g_check = resolve_pdb_chain_and_architecture_guards(
            neg_rid, "M" + "A" * 160, 0.72, 0.72, 0.85, 0.85, verified_n_terminal_helix=True
        )
        blocked_ok = not bool(g_check["has_full_mechanism_support"])
        if not blocked_ok:
            raise RuntimeError(f"Confusable negative control {neg_rid} failed guard check: {g_check}")
        ctrl_audit_rows.append(
            {
                "control_category": "confusable_negative_architecture_guard",
                "ref_id": neg_rid,
                "raw_target_header": neg_rid,
                "resolved_ref_id": g_check["resolved_ref_id"],
                "metric_primary": 0.72,
                "qcov": 0.85,
                "tcov": 0.85,
                "experimentally_supported_roles": spec.get("experimentally_supported_roles", "none"),
                "gate_passed": blocked_ok,
                "audit_note": g_check["pdb_guard_status"],
            }
        )
    pd.DataFrame(ctrl_audit_rows).to_csv(P_RES / "positive_and_negative_search_controls_audit.tsv", sep="\t", index=False)

    struct_hits_by_cat: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    all_fs_rows: list[dict[str, Any]] = []

    for src_name, db_rel in (
        ("ColabFold_AF2", "results/c10_structclust/colabfold/db/all"),
        ("ESMFold", "results/c10_structclust/esmfold/db/all"),
    ):
        q_db = P_ANNOT / db_rel
        if not Path(f"{q_db}.dbtype").exists():
            continue
        aln_db = fs_tmp / f"aln_{src_name}"
        out_m8 = fs_tmp / f"fs_{src_name}.m8"
        tmp_work = fs_tmp / f"work_{src_name}"
        tmp_work.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                "singularity",
                "exec",
                "-B",
                "/project,/scratch,/dev/shm",
                str(P_SIF_FOLDSEEK),
                "foldseek",
                "search",
                str(q_db),
                str(target_db),
                str(aln_db),
                str(tmp_work),
                "--alignment-type",
                "1",
                "-a",
                "-e",
                "10.0",
                "--threads",
                str(threads),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        subprocess.run(
            [
                "singularity",
                "exec",
                "-B",
                "/project,/scratch,/dev/shm",
                str(P_SIF_FOLDSEEK),
                "foldseek",
                "convertalis",
                str(q_db),
                str(target_db),
                str(aln_db),
                str(out_m8),
                "--threads",
                str(min(threads, 32)),
                "--format-output",
                "query,target,fident,alnlen,qstart,qend,qlen,tstart,tend,tlen,qcov,tcov,qtmscore,ttmscore,alntmscore,evalue,bits",
            ],
            check=True,
            stdout=subprocess.DEVNULL,
        )
        if out_m8.exists():
            with open(out_m8, "rt") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 17:
                        continue
                    q_sname, t_raw, fident, alnlen, qst, qen, qlen, tst, ten, tlen, qcov, tcov, qtm, ttm, alntm, ev, bits = parts[:17]
                    # Strictly map Foldseek query via members.tsv.gz lookup; NEVER rewrite unmapped sname (`EXP4-PIPELINE-REVIEW.txt` Item 1c)
                    matched_pids = sname_to_pids.get((src_name, q_sname), [])
                    if not matched_pids:
                        continue
                    ref_id = resolve_curated_reference_id_from_header(
                        t_raw,
                        target_polymer_seq=target_poly_map.get(t_raw),
                        db_lookup_map=target_lookup_map,
                        require_known=True,
                    )
                    spec = PANEL_SPECS_BY_ID[ref_id]
                    for m_pid in matched_pids:
                        rec = {
                            "cat_id": m_pid,
                            "model_source": src_name,
                            "struct_cluster": q_sname,
                            "ref_id": ref_id,
                            "fident": float(fident),
                            "alnlen": int(alnlen),
                            "qlen": int(qlen),
                            "tlen": int(tlen),
                            "qcov": float(qcov),
                            "tcov": float(tcov),
                            "qtmscore": float(qtm),
                            "ttmscore": float(ttm),
                            "alntmscore": float(alntm),
                            "evalue": float(ev),
                            "bits": float(bits),
                            "experimentally_supported_roles": spec.get("experimentally_supported_roles", "none"),
                            "hypothesized_roles": spec.get("hypothesized_roles", "none"),
                            "control_type": spec.get("control_type", ""),
                        }
                        struct_hits_by_cat[m_pid].append(rec)
                        all_fs_rows.append(rec)

    pd.DataFrame(all_fs_rows).to_csv(P_RES / "curated_panel_foldseek_hits.tsv", sep="\t", index=False)
    return seq_hits_by_gene, struct_hits_by_cat


def run_mafft_and_iqtree_for_family(
    family_id: str,
    member_seqs: dict[str, str],
    run_iqtree: bool = False,
) -> dict[str, Any]:
    """Run MAFFT (`mafft.sif`) and optionally IQ-TREE 3 (`iqtree.sif`, `LG+G4`) on verified member AA sequences."""
    P_MSA_OUT.mkdir(parents=True, exist_ok=True)
    safe_fid = family_id.replace("::", "__").replace("/", "_")
    faa_path = P_MSA_OUT / f"{safe_fid}.faa"
    aln_path = P_MSA_OUT / f"{safe_fid}.aln.faa"
    tree_prefix = P_MSA_OUT / f"{safe_fid}.iqtree"

    all_seqs = {mid: mseq for mid, mseq in member_seqs.items() if mseq}
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
        }

    if len(all_seqs) == 1:
        metrics = compute_alignment_summary_metrics(all_seqs)
        return {
            "msa_status": "single_sequence_no_msa_needed",
            **metrics,
            "tree_status": "single_sequence_tree_not_applicable",
            "tree_file": "none",
        }

    with open(faa_path, "wt") as fh:
        for sid, sseq in all_seqs.items():
            fh.write(f">{sid}\n{sseq}\n")

    aligned_map: dict[str, str] = {}
    if P_SIF_MAFFT.exists():
        cmd = [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch,/dev/shm",
            str(P_SIF_MAFFT),
            "mafft",
            "--quiet",
            "--retree",
            "2",
            "--thread",
            "4",
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

    metrics = compute_alignment_summary_metrics(aligned_map if aligned_map else all_seqs)

    unique_by_seq: dict[str, str] = {}
    for sid, aseq in aligned_map.items():
        if aseq not in unique_by_seq.values():
            unique_by_seq[sid] = aseq

    tree_status = "iqtree_not_requested_for_non_focal_family"
    tree_file_str = "none"
    if run_iqtree:
        if len(unique_by_seq) < 4:
            tree_status = f"fewer_than_4_nonredundant_sequences_tree_skipped(n_nr={len(unique_by_seq)})"
        elif P_SIF_IQTREE.exists() and aln_path.exists():
            nr_aln_path = P_MSA_OUT / f"{safe_fid}.nr.aln.faa"
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
                tree_file_str = str(tf.relative_to(P_EXP4_ROOT)) if tf.is_relative_to(P_EXP4_ROOT) else str(tf)
            else:
                tree_status = f"iqtree3_failed_or_constant_alignment(n_nr={len(unique_by_seq)})"

    return {
        "msa_status": msa_status,
        **metrics,
        "tree_status": tree_status,
        "tree_file": tree_file_str,
    }


def evaluate_shard_bundle(shard_id: int) -> dict[str, Any]:
    """Evaluate all families assigned to `shard_id` (`0..16`) using pre-staged shard input bundle."""
    P_SHARD_OUT.mkdir(parents=True, exist_ok=True)
    bundle_path = require_file(P_SHARD_IN / f"shard_{shard_id}.pkl", f"shard_{shard_id}.pkl")
    with open(bundle_path, "rb") as fh:
        bundle = pickle.load(fh)

    families_meta: list[dict[str, Any]] = bundle["families_meta"]
    members_by_fam: dict[str, list[dict[str, Any]]] = bundle["members_by_fam"]
    contig_ordered_genes: dict[tuple[str, str], list[str]] = bundle["contig_ordered_genes"]
    flank_annot_map: dict[str, dict[str, Any]] = bundle["flank_annot_map"]
    gene_profiles_map: dict[str, list[str]] = bundle["gene_profiles_map"]
    raw_ev: dict[str, Any] = bundle["raw_ev"]
    seq_hits_by_gene: dict[str, list[dict[str, Any]]] = bundle["seq_hits_by_gene"]
    struct_hits_by_cat: dict[str, list[dict[str, Any]]] = bundle["struct_hits_by_cat"]

    eval_rows: list[dict[str, Any]] = []
    competing_rows: list[dict[str, Any]] = []
    neighbourhood_rows: list[dict[str, Any]] = []
    struct_hit_audit_rows: list[dict[str, Any]] = []

    null_guilds_to_test = [
        "nutrient_acquisition_peptidase",
        "substrate_binding_or_uptake",
        "surface_assembly_pilin_archaellin",
        "cell_coat_or_envelope",
        "host_defense_antiviral",
        "viral_or_mobile_element_origin",
    ]

    top_mag_by_genus_block = bundle.get("top_mag_by_genus_block")
    if top_mag_by_genus_block is None:
        canon_rep_path = P_LEG_SCI / "inputs/refined_representatives.tsv"
        if canon_rep_path.exists():
            c_df = pd.read_csv(canon_rep_path, sep="\t")
            top_mag_by_genus_block = compute_canonical_top_genome_per_genus_block(c_df)
        else:
            top_mag_by_genus_block = {}

    for fam_rec in families_meta:
        fid = str(fam_rec["family_id"])
        cohort_group = str(fam_rec["cohort_group"])
        cand_role = str(fam_rec["candidate_role"])
        is_dk = bool(fam_rec["is_dark"])
        fam_mems = members_by_fam.get(fid, [])
        if not fam_mems:
            continue

        ver_mems = [m for m in fam_mems if m["verified_complete_interior_member"]]
        active_mems = ver_mems if ver_mems else fam_mems

        # Representative locus selection
        exact_rep = [m for m in active_mems if m["cat_id"] == fid]
        if exact_rep:
            rep_m = exact_rep[0]
        else:
            rep_m = sorted(active_mems, key=lambda x: (-int(x["length"]), str(x["gene"])))[0]

        rep_gene = str(rep_m["gene"])
        rep_cat = str(rep_m["cat_id"])
        rep_genome = str(rep_m["genome"])
        rep_seq = str(rep_m["raw_aa_sequence"] or "")

        # Sequence composition, N-X-S/T sequons, prepilin/archaeosortase motifs on representative
        seq_comp = compute_sequence_composition_and_motifs(rep_seq)

        # Surface topology classification on representative with real TMbed 1-based helix segment coordinates
        rep_topo = classify_export_and_surface_topology(
            rep_m["signalp_class"],
            rep_m["tmbed_S"],
            rep_m["tmbed_H"],
            rep_m["tmbed_B"],
            aa_length=seq_comp["aa_length"] if seq_comp["aa_length"] > 0 else None,
            signalp_cs_pos=rep_m["signalp_cs_pos"],
            tm_helix_spans_1based=rep_m["tmbed_h_spans"],
            has_prepilin_motif=bool(seq_comp["has_prepilin_classIII_motif_pred"]),
            has_archaeosortase_cterm_motif=bool(seq_comp["has_archaeosortase_cterm_motif_pred"]),
            sequence_valid=bool(seq_comp["aa_length"] > 0 and seq_comp["n_internal_stops"] == 0),
        )

        # Family-wide topology & export fractions across verified complete interior members
        member_topos = [
            classify_export_and_surface_topology(
                m["signalp_class"],
                m["tmbed_S"],
                m["tmbed_H"],
                m["tmbed_B"],
                aa_length=int(m["length"]),
                signalp_cs_pos=m["signalp_cs_pos"],
                tm_helix_spans_1based=m["tmbed_h_spans"],
                sequence_valid=True,
            )
            for m in active_mems
        ]
        known_exp_vals = [t["is_exported_any"] for t in member_topos if t["is_exported_any"] is not None]
        known_sec_vals = [t["is_freely_secreted_candidate"] for t in member_topos if t["is_freely_secreted_candidate"] is not None]
        known_anc_vals = [t["is_surface_or_membrane_anchored"] for t in member_topos if t["is_surface_or_membrane_anchored"] is not None]
        n_missing_tmbed = sum(1 for m in active_mems if m["tmbed_S"] is None)

        frac_exported = float(np.mean(known_exp_vals)) if known_exp_vals else np.nan
        frac_freely_secreted = float(np.mean(known_sec_vals)) if known_sec_vals else np.nan
        frac_surface_anchored = float(np.mean(known_anc_vals)) if known_anc_vals else np.nan

        # Gather curated sequence hits and qualify every structural hit per-member and per-source (`EXP4-PIPELINE-REVIEW.txt` Item 2)
        fam_curated_seq_hits: list[dict[str, Any]] = []
        fam_curated_struct_hits: list[dict[str, Any]] = []
        rep_qualified_struct_hits: list[dict[str, Any]] = []
        for m in active_mems:
            m_gene_id = str(m["gene"])
            m_cid = str(m["cat_id"])
            m_aa = str(m.get("raw_aa_sequence") or "")
            fam_curated_seq_hits.extend(seq_hits_by_gene.get(m_gene_id, []))
            for raw_h in struct_hits_by_cat.get(m_cid, []):
                m_src = str(raw_h.get("model_source", ""))
                mod_info = raw_ev["struct_models"].get((m_cid, m_src), {})
                q_hit = qualify_structural_hit_for_member(
                    raw_h,
                    member_gene=m_gene_id,
                    member_cat_id=m_cid,
                    member_aa_seq=m_aa,
                    model_info=mod_info,
                )
                fam_curated_struct_hits.append(q_hit)
                struct_hit_audit_rows.append({"family_id": fid, "cohort_group": cohort_group, **q_hit})
                if m_cid == rep_cat and bool(q_hit.get("qualified_for_role_evidence", False)):
                    rep_qualified_struct_hits.append(q_hit)

        # Best raw PDB hit (qualified by exact sequence & pLDDT on representative) and qualified curated PDB hit on representative
        rep_pdb_hits: list[dict[str, Any]] = []
        for h_src in ("ColabFold_AF2", "ESMFold"):
            mod_rep = raw_ev["struct_models"].get((rep_cat, h_src), {})
            p_chk_rep = verify_prediction_sequence_hash_match(rep_seq, mod_rep.get("model_polymer_aa_seq"))
            if p_chk_rep["is_valid_exact_match"] and float(mod_rep.get("plddt_0to1") or 0.0) >= 0.70 and len(rep_seq) >= 100:
                rep_pdb_hits.extend(raw_ev["pdb_hits"].get((rep_cat, h_src), []))
        rep_curated_struct = rep_qualified_struct_hits

        # Verify prediction polymer sequence hash match against candidate AA sequence
        cf_mod = raw_ev["struct_models"].get((rep_cat, "ColabFold_AF2"), {})
        esm_mod = raw_ev["struct_models"].get((rep_cat, "ESMFold"), {})
        cf_chk = verify_prediction_sequence_hash_match(rep_seq, cf_mod.get("model_polymer_aa_seq"))
        esm_chk = verify_prediction_sequence_hash_match(rep_seq, esm_mod.get("model_polymer_aa_seq"))
        if cf_chk["is_valid_exact_match"]:
            best_mod = cf_mod
        elif esm_chk["is_valid_exact_match"]:
            best_mod = esm_mod
        else:
            best_mod = cf_mod if cf_mod else esm_mod
        poly_check = verify_prediction_sequence_hash_match(rep_seq, best_mod.get("model_polymer_aa_seq"))

        # Build +/-5 gene neighbourhood panels and run grouped within-contig circular permutation nulls per role guild
        member_contigs_for_null: list[dict[str, Any]] = []
        guild_genus_blocks: dict[str, set[str]] = collections.defaultdict(set)

        for m in active_mems:
            m_gene = str(m["gene"])
            m_mag = str(m["genome"])
            m_contig = str(m["contig"])
            m_gblock = str(m["genus_singleton_block"])
            c_genes = contig_ordered_genes.get((m_mag, m_contig), [])
            if m_gene not in c_genes:
                continue
            cand_idx = c_genes.index(m_gene)

            contig_guilds: list[set[str]] = []
            for idx_g, fg in enumerate(c_genes):
                fa = flank_annot_map.get(fg, {})
                f_desc = str(fa.get("functional_description", "hypothetical protein"))
                f_profs = gene_profiles_map.get(fg, [])
                g_set = classify_flanking_machinery_guilds(f_desc, f_profs)
                contig_guilds.append(g_set)
                if abs(idx_g - cand_idx) <= 5 and idx_g != cand_idx:
                    for g_tag in g_set:
                        if g_tag != "generic_any_known_annotation":
                            guild_genus_blocks[g_tag].add(m_gblock)
                    if (
                        shard_id == 0
                        or fid in STAGE3_FOCAL_DARK_CANDIDATES
                        or fid in PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES
                        or fid in EXPLORATORY_DARK_CANDIDATES_V1
                        or fid in PREDECLARED_SMALL_P_FAMILIES_V1
                    ):
                        neighbourhood_rows.append(
                            {
                                "family_id": fid,
                                "cohort_group": cohort_group,
                                "member_gene": m_gene,
                                "genome": m_mag,
                                "tax_class": str(m["tax_class"]),
                                "tax_order": str(m["tax_order"]),
                                "genus_singleton_block": m_gblock,
                                "contig": m_contig,
                                "relative_gene_offset": int(idx_g - cand_idx),
                                "flanking_gene": fg,
                                "flanking_family_id": str(fa.get("family_id", "missing")),
                                "flanking_is_dark": bool(fa.get("is_dark", True)),
                                "flanking_signalp_class": str(fa.get("signalp_class", "missing")),
                                "flanking_functional_description": f_desc,
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

        role_null_results: dict[str, dict[str, Any]] = {}
        generic_null_res: dict[str, Any] = {"obs_adjacent": 0, "null_mean": 0.0, "p_empirical_one_sided": 1.0}
        n_unique_contigs = 0
        for g_idx, g_name in enumerate(null_guilds_to_test):
            nb = run_grouped_contig_circular_nulls(
                member_contigs_for_null,
                specific_guild=g_name,
                specific_guild_selection_mode="pre_specified_competing_role_guild",
                window=5,
                n_perm=1000,
                seed=20261109 + g_idx,
                genus_dereplicate_map=top_mag_by_genus_block,
            )
            n_unique_contigs = int(nb["n_unique_contigs"])
            s_null = dict(nb["specific_null"])
            s_null["n_genus_blocks_adjacent"] = len(guild_genus_blocks.get(g_name, set()))
            gd_null = dict(nb.get("genus_dereplicated_null", {}))
            s_null["genus_derep_obs_adjacent"] = int(gd_null.get("obs_adjacent", 0))
            s_null["genus_derep_null_mean"] = float(gd_null.get("null_mean", 0.0))
            s_null["genus_derep_p_exact"] = float(
                gd_null.get("p_exact_one_sided", gd_null.get("p_empirical_one_sided", np.nan))
            )
            s_null["genus_derep_status"] = str(gd_null.get("status", "unknown"))
            s_null["genus_derep_n_genomes_retained"] = int(gd_null.get("n_genomes_retained", 0))
            s_null["genus_derep_n_genus_blocks_retained"] = int(gd_null.get("n_genus_blocks_retained", 0))
            role_null_results[g_name] = s_null
            if g_idx == 0:
                generic_null_res = dict(nb["generic_null"])

        # Populate legacy alias so coat lookups under old key continue to work
        role_null_results["s_layer_or_charted_coat"] = dict(role_null_results["cell_coat_or_envelope"])

        # Apply within-family BH FDR across the 6 tested specific guilds
        guild_pvals = [float(role_null_results[g]["p_empirical_one_sided"]) for g in null_guilds_to_test]
        guild_qvals = apply_bh_fdr(guild_pvals)
        for g_name, qv in zip(null_guilds_to_test, guild_qvals):
            role_null_results[g_name]["q_bh"] = float(qv)

        # Evaluate PDB chain & architecture guards on best structural hit
        has_ver_n_helix = rep_topo["tm_anchor_position_status"] in (
            "verified_n_terminal_prepilin_alpha1_helix",
            "verified_n_terminal_tm_anchor",
        ) or str(rep_topo["tm_anchor_position_status"]).startswith("verified_n_terminal_tm_anchor")

        if rep_curated_struct:
            best_cs = sorted(rep_curated_struct, key=lambda x: (-float(x["qtmscore"]), -float(x["tcov"])))[0]
            pdb_guard = resolve_pdb_chain_and_architecture_guards(
                best_cs["ref_id"],
                rep_seq,
                best_cs["qtmscore"],
                best_cs["ttmscore"],
                best_cs["qcov"],
                best_cs["tcov"],
                verified_n_terminal_helix=has_ver_n_helix,
                verified_alignment_active_site_conserved=False,
                has_verified_operon_context=any(
                    int(role_null_results[g]["obs_adjacent"]) >= 2 and float(role_null_results[g]["p_empirical_one_sided"]) <= 0.05
                    for g in null_guilds_to_test
                ),
            )
        elif rep_pdb_hits:
            best_raw_pdb = sorted(rep_pdb_hits, key=lambda x: -float(x["qtmscore"]))[0]
            pdb_guard = resolve_pdb_chain_and_architecture_guards(
                best_raw_pdb["target"],
                rep_seq,
                best_raw_pdb["qtmscore"],
                best_raw_pdb["ttmscore"],
                best_raw_pdb["qcov"],
                best_raw_pdb["tcov"],
                verified_n_terminal_helix=has_ver_n_helix,
            )
        else:
            pdb_guard = resolve_pdb_chain_and_architecture_guards(
                "missing", rep_seq, np.nan, np.nan, np.nan, np.nan
            )

        own_dom_descs = list(dict.fromkeys(raw_ev["ipr_domains"].get(rep_cat, []) + raw_ev["ipr_domains"].get(fid, [])))
        hh_info = raw_ev["hh_hits"].get(rep_cat) or raw_ev["hh_hits"].get(fid) or {}
        if hh_info.get("pfam_id") not in (None, "", "missing"):
            own_dom_descs.append(f"HHsearch:{hh_info['pfam_id']}({hh_info.get('pfam_desc', '')})")
        own_c6 = list(raw_ev["c6_hits"].get(rep_gene, []))

        tax_classes_present = {str(m["tax_class"]) for m in active_mems if m.get("tax_class")}
        tax_orders_present = {str(m["tax_order"]) for m in active_mems if m.get("tax_order")}
        genus_blocks_present = {str(m["genus_singleton_block"]) for m in active_mems}

        comp_res = evaluate_competing_roles_for_family(
            family_id=fid,
            candidate_role=cand_role,
            is_dark=is_dk,
            surface_localization_class=str(rep_topo["surface_localization_class"]),
            frac_exported=frac_exported,
            curated_seq_hits=fam_curated_seq_hits,
            curated_struct_hits=fam_curated_struct_hits,
            own_domain_descriptions=own_dom_descs,
            own_c6_hits=own_c6,
            pdb_guard_result=pdb_guard,
            role_null_results=role_null_results,
            tax_classes_present=tax_classes_present,
        )

        # Run MAFFT and IQ-TREE on Pilot shard (shard 0), focal Stage 3 candidates, or multi-order convergent candidates
        run_msa = bool(
            shard_id == 0
            or fid in STAGE3_FOCAL_DARK_CANDIDATES
            or fid in PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES
            or comp_res["max_role_evidence_score"] >= 0.40
        )
        if run_msa:
            ver_seq_map = {str(m["gene"]): str(m["raw_aa_sequence"]) for m in active_mems if m.get("raw_aa_sequence")}
            msa_res = run_mafft_and_iqtree_for_family(
                fid,
                ver_seq_map,
                run_iqtree=(shard_id == 0 or fid in STAGE3_FOCAL_DARK_CANDIDATES or comp_res["max_role_evidence_score"] >= 0.65),
            )
        else:
            msa_res = {
                "msa_status": "deferred_non_focal_no_role_hit",
                "n_aligned_sequences": len(active_mems),
                "alignment_cols": int(seq_comp["aa_length"]),
                "median_ungapped_coverage": np.nan,
                "min_pairwise_aa_identity": np.nan,
                "median_pairwise_aa_identity": np.nan,
                "max_pairwise_aa_identity": np.nan,
                "tree_status": "deferred_non_focal_no_role_hit",
                "tree_file": "none",
            }

        best_seq_hit = sorted(fam_curated_seq_hits, key=lambda x: -float(x["bits"]))[0] if fam_curated_seq_hits else {}
        qualified_fam_struct_hits = [
            h for h in fam_curated_struct_hits if bool(h.get("qualified_for_role_evidence", False))
        ]
        best_struct_hit = (
            sorted(qualified_fam_struct_hits, key=lambda x: -float(x["qtmscore"]))[0]
            if qualified_fam_struct_hits
            else {}
        )

        competing_rows.append(
            {
                "shard_id": shard_id,
                "cohort_group": cohort_group,
                **comp_res,
            }
        )

        eval_rows.append(
            {
                "shard_id": shard_id,
                "family_id": fid,
                "cohort_group": cohort_group,
                "candidate_role": cand_role,
                "is_dark": is_dk,
                "representative_gene": rep_gene,
                "representative_cat_id": rep_cat,
                "representative_genome": rep_genome,
                "n_total_member_loci": len(fam_mems),
                "n_verified_complete_interior_loci": len(ver_mems),
                "n_verified_classes": len(tax_classes_present),
                "verified_classes": ";".join(sorted(tax_classes_present)),
                "n_verified_orders": len(tax_orders_present),
                "verified_orders": ";".join(sorted(tax_orders_present)),
                "n_verified_independent_genus_blocks": len(genus_blocks_present),
                "rep_signalp_class": str(rep_m["signalp_class"]),
                "rep_signalp_cs_pos": rep_m["signalp_cs_pos"] if rep_m["signalp_cs_pos"] is not None else np.nan,
                "rep_tmbed_S": rep_m["tmbed_S"],
                "rep_tmbed_H": rep_m["tmbed_H"] if rep_m["tmbed_H"] is not None else np.nan,
                "rep_tmbed_B": rep_m["tmbed_B"] if rep_m["tmbed_B"] is not None else np.nan,
                "rep_tmbed_h_spans_1based": (
                    ";".join(f"{s}-{e}" for s, e in rep_m["tmbed_h_spans"])
                    if rep_m["tmbed_h_spans"]
                    else ("missing_tmbed" if rep_m["tmbed_S"] is None else "none")
                ),
                "n_members_missing_tmbed": n_missing_tmbed,
                **rep_topo,
                "frac_family_exported_known": frac_exported,
                "frac_family_freely_secreted_known": frac_freely_secreted,
                "frac_family_surface_or_membrane_anchored_known": frac_surface_anchored,
                **seq_comp,
                "model_source": best_mod.get("model_source", "none"),
                "model_plddt_0to1": best_mod.get("plddt_0to1", np.nan),
                "model_polymer_hash_status": poly_check["sequence_match_status"],
                "model_polymer_exact_match": poly_check["is_valid_exact_match"],
                "n_struct_hits_qualified": int(comp_res.get("n_struct_hits_qualified", len(qualified_fam_struct_hits))),
                "n_struct_hits_rejected_gate": int(comp_res.get("n_struct_hits_rejected_gate", len(fam_curated_struct_hits) - len(qualified_fam_struct_hits))),
                "best_curated_seq_ref_id": best_seq_hit.get("ref_id", "none"),
                "best_curated_seq_pident": best_seq_hit.get("pident", np.nan),
                "best_curated_seq_qcov": best_seq_hit.get("qcov", np.nan),
                "best_curated_seq_tcov": best_seq_hit.get("tcov", np.nan),
                "best_curated_seq_evalue": best_seq_hit.get("evalue", np.nan),
                "best_curated_seq_bits": best_seq_hit.get("bits", np.nan),
                "best_curated_struct_ref_id": best_struct_hit.get("ref_id", "none"),
                "best_curated_struct_qtmscore": best_struct_hit.get("qtmscore", np.nan),
                "best_curated_struct_ttmscore": best_struct_hit.get("ttmscore", np.nan),
                "best_curated_struct_qcov": best_struct_hit.get("qcov", np.nan),
                "best_curated_struct_tcov": best_struct_hit.get("tcov", np.nan),
                "pdb_guard_status": pdb_guard["pdb_guard_status"],
                "pdb_resolved_ref_id": pdb_guard["resolved_ref_id"],
                "pdb_structural_fold_compatible": pdb_guard["structural_fold_compatible"],
                "pdb_has_full_mechanism_support": pdb_guard["has_full_mechanism_support"],
                "rep_own_domains": ";".join(own_dom_descs) if own_dom_descs else "none",
                "rep_own_c6_hits": ";".join(own_c6) if own_c6 else "none",
                "n_unique_contigs_evaluated": n_unique_contigs,
                "generic_known_flank_obs_adj": int(generic_null_res["obs_adjacent"]),
                "generic_known_flank_p_emp": float(generic_null_res["p_empirical_one_sided"]),
                **{
                    f"null_obs_adj_{g}": int(role_null_results[g]["obs_adjacent"])
                    for g in null_guilds_to_test
                },
                **{
                    f"null_n_gblocks_{g}": int(role_null_results[g]["n_genus_blocks_adjacent"])
                    for g in null_guilds_to_test
                },
                **{
                    f"null_p_emp_{g}": float(role_null_results[g]["p_empirical_one_sided"])
                    for g in null_guilds_to_test
                },
                **{
                    f"null_p_exact_{g}": float(role_null_results[g].get("p_exact_one_sided", role_null_results[g]["p_empirical_one_sided"]))
                    for g in null_guilds_to_test
                },
                **{
                    f"null_q_bh_{g}": float(role_null_results[g]["q_bh"])
                    for g in null_guilds_to_test
                },
                **{
                    f"genus_derep_obs_adj_{g}": int(role_null_results[g].get("genus_derep_obs_adjacent", 0))
                    for g in null_guilds_to_test
                },
                **{
                    f"genus_derep_p_exact_{g}": float(role_null_results[g].get("genus_derep_p_exact", np.nan))
                    for g in null_guilds_to_test
                },
                **{
                    f"genus_derep_n_gblocks_{g}": int(role_null_results[g].get("genus_derep_n_genus_blocks_retained", 0))
                    for g in null_guilds_to_test
                },
                "null_obs_adj_s_layer_or_charted_coat": int(role_null_results["cell_coat_or_envelope"]["obs_adjacent"]),
                "null_n_gblocks_s_layer_or_charted_coat": int(role_null_results["cell_coat_or_envelope"]["n_genus_blocks_adjacent"]),
                "null_p_emp_s_layer_or_charted_coat": float(role_null_results["cell_coat_or_envelope"]["p_empirical_one_sided"]),
                "null_p_exact_s_layer_or_charted_coat": float(role_null_results["cell_coat_or_envelope"].get("p_exact_one_sided", role_null_results["cell_coat_or_envelope"]["p_empirical_one_sided"])),
                "null_q_bh_s_layer_or_charted_coat": float(role_null_results["cell_coat_or_envelope"]["q_bh"]),
                "sbp_substrate_specificity_status": comp_res["sbp_substrate_specificity_status"],
                "viral_origin_audit_status": comp_res["viral_origin_audit_status"],
                "supported_overlapping_roles": comp_res["supported_overlapping_roles"],
                "n_supported_overlapping_roles": comp_res["n_supported_overlapping_roles"],
                "top_competing_roles": comp_res["top_competing_roles"],
                "max_role_evidence_score": comp_res["max_role_evidence_score"],
                "role_discrimination_outcome": comp_res["role_discrimination_outcome"],
                **{f"tier_{r}": comp_res[f"tier_{r}"] for r in ROLE_CLASSES},
                **{f"score_{r}": comp_res[f"score_{r}"] for r in ROLE_CLASSES},
                **{f"detail_{r}": comp_res[f"detail_{r}"] for r in ROLE_CLASSES},
                "msa_status": msa_res["msa_status"],
                "n_aligned_sequences": msa_res["n_aligned_sequences"],
                "alignment_cols": msa_res["alignment_cols"],
                "median_ungapped_coverage": msa_res["median_ungapped_coverage"],
                "median_pairwise_aa_identity": msa_res["median_pairwise_aa_identity"],
                "tree_status": msa_res["tree_status"],
                "tree_file": msa_res["tree_file"],
            }
        )

    eval_df = pd.DataFrame(eval_rows)
    comp_df = pd.DataFrame(competing_rows)
    neigh_df = pd.DataFrame(neighbourhood_rows)
    struct_audit_df = pd.DataFrame(struct_hit_audit_rows)

    eval_out = P_SHARD_OUT / f"shard_{shard_id}_family_evaluation.tsv"
    comp_out = P_SHARD_OUT / f"shard_{shard_id}_competing_roles.tsv"
    neigh_out = P_SHARD_OUT / f"shard_{shard_id}_neighbourhoods.tsv"
    struct_audit_out = P_SHARD_OUT / f"shard_{shard_id}_struct_hits_audit.tsv"
    done_out = P_SHARD_OUT / f"shard_{shard_id}.done"

    eval_df.to_csv(eval_out, sep="\t", index=False)
    comp_df.to_csv(comp_out, sep="\t", index=False)
    neigh_df.to_csv(neigh_out, sep="\t", index=False)
    struct_audit_df.to_csv(struct_audit_out, sep="\t", index=False)
    done_out.write_text(f"completed_shard_{shard_id}_n_families={len(eval_df)}\n", encoding="utf-8")

    return {
        "shard_id": shard_id,
        "n_families_evaluated": len(eval_df),
        "n_neighbourhood_rows": len(neigh_df),
        "n_struct_hits_audited": len(struct_audit_df),
        "eval_tsv": str(eval_out),
    }


def run_freeze_and_pilot(threads: int = 64) -> None:
    """Freeze the exact candidate + control family universe, run curated panel DIAMOND/Foldseek,
    build per-shard input bundles (`0..16`), and execute + validate the Pilot shard (`shard_id = 0`).
    """
    for d in (P_PROV, P_RES, P_LOGS, P_TMP, P_SHARD_IN, P_SHARD_OUT, P_MSA_OUT):
        d.mkdir(parents=True, exist_ok=True)

    # 1. Load canonical Stage 1b verified family rankings & controls
    r2_path = require_file(
        P_CANON_SCI / "results/task_a_stage1b_route2_exported_dark_family_ranking_verified.tsv",
        "Canonical Stage 1b Route A2 exported dark family ranking",
    )
    r1_path = require_file(
        P_CANON_SCI / "results/task_a_stage1b_route1_gapfill_family_ranking_verified.tsv",
        "Canonical Stage 1b Route A1 gapfill family ranking",
    )
    neg_ctrl_path = require_file(
        P_CANON_SCI / "results/task_a_route2_annotation_negative_exported_controls.tsv",
        "Route A2 annotation-negative exported controls",
    )
    canon_rep_path = require_file(
        P_LEG_SCI / "inputs/refined_representatives.tsv",
        "Canonical refined_representatives.tsv",
    )

    r2_df = pd.read_csv(r2_path, sep="\t")
    r1_df = pd.read_csv(r1_path, sep="\t")
    neg_ctrl_df = pd.read_csv(neg_ctrl_path, sep="\t")
    canon_rep_df = pd.read_csv(canon_rep_path, sep="\t")

    # 2. Load canonical taxonomy for all 237 EQ Thermoplasmatota MAGs
    meta_path = require_file(P_ARCH / "inputs/thermoplasmatota.tsv", "thermoplasmatota.tsv")
    meta_df = pd.read_csv(meta_path, sep="\t")
    thermo_set = set(meta_df["genome"].astype(str))
    canon_sub = canon_rep_df[canon_rep_df["mag"].astype(str).isin(thermo_set)].copy()
    if len(canon_sub) != 237:
        raise ValueError(f"Expected 237 EQ Thermoplasmatota MAGs in refined_representatives.tsv, got {len(canon_sub)}")

    tax_col = "gtdb_r232_classification" if "gtdb_r232_classification" in canon_sub.columns else "gtdb_taxonomy"
    tax_parsed = [parse_gtdb_lineage(x) for x in canon_sub[tax_col]]
    for rk, col_name in [("c", "tax_class"), ("o", "tax_order"), ("f", "tax_family"), ("g", "tax_genus")]:
        canon_sub[col_name] = [d[rk] for d in tax_parsed]
    g_blocks, _ = normalize_genus_groups(canon_sub["mag"].astype(str).tolist(), canon_sub["tax_genus"].astype(str).tolist())
    canon_sub["genus_singleton_block"] = g_blocks

    tax_class_map = dict(zip(canon_sub["mag"].astype(str), canon_sub["tax_class"].astype(str)))
    tax_order_map = dict(zip(canon_sub["mag"].astype(str), canon_sub["tax_order"].astype(str)))
    tax_genus_map = dict(zip(canon_sub["mag"].astype(str), canon_sub["tax_genus"].astype(str)))
    genus_block_map = dict(zip(canon_sub["mag"].astype(str), canon_sub["genus_singleton_block"].astype(str)))

    # 3. Load archaeal_genes.tsv and family assignments across all 237 Thermoplasmatota genomes
    genes_path = require_file(P_ARCH / "inputs/archaeal_genes.tsv", "archaeal_genes.tsv")
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
    genes_thermo["tax_class"] = genes_thermo["genome"].astype(str).map(lambda g: tax_class_map.get(g, "unknown"))
    genes_thermo["tax_order"] = genes_thermo["genome"].astype(str).map(lambda g: tax_order_map.get(g, "unknown"))
    genes_thermo["genus_singleton_block"] = genes_thermo["genome"].astype(str).map(
        lambda g: genus_block_map.get(g, f"UNCLASSIFIED_GENUS_SINGLETON::{g}")
    )

    # Load SignalP (with CS pos), TMbed (with exact 1-based helix spans), and C6 hits
    tm_all, sp_all = load_all_sp_tm_with_spans(n_workers=min(threads, 32))
    cazy_gene, mer_gene, tc_gene = load_c6_thermo_gene_hits(thermo_set)

    sp_cls_l: list[str | None] = []
    sp_cs_l: list[int | None] = []
    sp_any_l: list[bool] = []
    tm_s_l: list[bool | None] = []
    tm_h_l: list[int | None] = []
    tm_b_l: list[int | None] = []
    tm_spans_l: list[list[tuple[int, int]] | None] = []
    missing_sp_or_tm_l: list[bool] = []
    verified_soluble_non_exp_l: list[bool] = []

    for g in genes_thermo["gene"].astype(str):
        sp_e = sp_all.get(("T", g)) or sp_all.get(("D", g))
        tm_e = tm_all.get(("T", g)) or tm_all.get(("D", g))
        sc = sp_e[0] if sp_e is not None else None
        cs_p = sp_e[3] if sp_e is not None else None
        sp_cls_l.append(sc)
        sp_cs_l.append(cs_p)
        sp_any_l.append(bool(sc in SPANY_CLASSES) if sc is not None else False)
        if tm_e is not None:
            t_s = bool(tm_e[1])
            t_h = int(tm_e[2])
            t_b = int(tm_e[3])
            tm_s_l.append(t_s)
            tm_h_l.append(t_h)
            tm_b_l.append(t_b)
            tm_spans_l.append(tm_e[5])
            missing_sp_or_tm_l.append(sc is None)
            verified_soluble_non_exp_l.append(bool(sc == "OTHER" and (not t_s) and t_h == 0 and t_b == 0))
        else:
            # Never coerce missing TMbed to False/0 (`EXP4-CODE-REVIEW.txt` item 2)!
            tm_s_l.append(None)
            tm_h_l.append(None)
            tm_b_l.append(None)
            tm_spans_l.append(None)
            missing_sp_or_tm_l.append(True)
            verified_soluble_non_exp_l.append(False)

    genes_thermo["signalp_class"] = sp_cls_l
    genes_thermo["signalp_cs_pos"] = sp_cs_l
    genes_thermo["sp_any"] = sp_any_l
    genes_thermo["tmbed_S"] = tm_s_l
    genes_thermo["tmbed_H"] = tm_h_l
    genes_thermo["tmbed_B"] = tm_b_l
    genes_thermo["tmbed_h_spans"] = tm_spans_l
    genes_thermo["missing_sp_or_tm"] = missing_sp_or_tm_l
    genes_thermo["verified_soluble_non_exported"] = verified_soluble_non_exp_l
    genes_thermo["has_c6_hit"] = [
        (g in tc_gene) or (g in mer_gene) or (g in cazy_gene) for g in genes_thermo["gene"].astype(str)
    ]

    # 4. Construct the frozen candidate + control family universe (`EXP4-PIPELINE-REVIEW.txt` Item 4)
    r2_elig = r2_df[r2_df["final_candidate_promotion_eligible"].astype(bool)].copy()
    r1_elig_dark = r1_df[
        r1_df["inherited_catalogue_dark"].astype(bool) & r1_df["final_candidate_promotion_eligible"].astype(bool)
    ].copy()
    r1_elig_nondark = r1_df[
        (~r1_df["inherited_catalogue_dark"].astype(bool)) & r1_df["rep_final_promotion_gate_pass"].astype(bool)
    ].copy()

    r2_elig_ids = [str(x) for x in r2_elig["family_id"].tolist()]
    r1_elig_dark_ids = [str(x) for x in r1_elig_dark["family_id"].tolist()]
    prior_ids = sorted(PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES)
    neg_ctrl_ids = [str(x) for x in neg_ctrl_df["family_id"].astype(str).unique().tolist()]

    # Build family summary across genes_thermo with strict missingness tracking
    fam_grp = genes_thermo.groupby("family_id", sort=False).agg(
        n_loci=("gene", "count"),
        n_gblocks=("genus_singleton_block", "nunique"),
        is_dark=("is_dark", "all"),
        frac_sp=("sp_any", "mean"),
        n_missing_sp_or_tm=("missing_sp_or_tm", "sum"),
        all_verified_soluble=("verified_soluble_non_exported", "all"),
        any_c6=("has_c6_hit", "any"),
        median_len=("length", "median"),
    ).reset_index()

    already_chosen = set(r2_elig_ids) | set(r1_elig_dark_ids) | set(prior_ids) | set(neg_ctrl_ids)

    # Known-function exported background panel (`n_missing_sp_or_tm == 0`, non-dark, exported frac >= 0.50, >= 2 genus blocks)
    known_exp_pool = fam_grp[
        (~fam_grp["is_dark"])
        & (fam_grp["n_missing_sp_or_tm"] == 0)
        & (fam_grp["frac_sp"] >= 0.50)
        & (fam_grp["n_gblocks"] >= 2)
        & (~fam_grp["family_id"].isin(already_chosen))
    ].sort_values(["any_c6", "n_gblocks", "n_loci", "family_id"], ascending=[False, False, False, True])
    known_exp_ctrl_ids = [str(x) for x in r1_elig_nondark["family_id"].tolist() if str(x) not in already_chosen][:20]
    for fid in known_exp_pool["family_id"].astype(str):
        if fid not in known_exp_ctrl_ids and fid not in already_chosen:
            known_exp_ctrl_ids.append(fid)
        if len(known_exp_ctrl_ids) >= 120:
            break
    already_chosen.update(known_exp_ctrl_ids)

    # Verified soluble non-exported dark background panel (`EXP4-PIPELINE-REVIEW.txt` Item 4a:
    # strictly requires `n_missing_sp_or_tm == 0` AND `all_verified_soluble == True`, matched across length & genus-block strata)
    non_exp_dark_pool = fam_grp[
        fam_grp["is_dark"]
        & (fam_grp["n_missing_sp_or_tm"] == 0)
        & fam_grp["all_verified_soluble"]
        & (fam_grp["frac_sp"] == 0.0)
        & (fam_grp["n_gblocks"] >= 2)
        & (fam_grp["median_len"] >= 90)
        & (fam_grp["median_len"] <= 600)
        & (~fam_grp["family_id"].isin(already_chosen))
    ].copy()
    # Stratify selection across length quartiles (90-180, 180-280, 280-400, 400-600) and genus-block counts (2, 3, >=4)
    non_exp_dark_ctrl_ids: list[str] = []
    for l_lo, l_hi in ((90, 180), (180, 280), (280, 400), (400, 601)):
        for gb_lo, gb_hi in ((2, 2), (3, 3), (4, 999)):
            sub_p = non_exp_dark_pool[
                (non_exp_dark_pool["median_len"] >= l_lo)
                & (non_exp_dark_pool["median_len"] < l_hi)
                & (non_exp_dark_pool["n_gblocks"] >= gb_lo)
                & (non_exp_dark_pool["n_gblocks"] <= gb_hi)
            ].sort_values(["n_gblocks", "n_loci", "family_id"], ascending=[False, False, True])
            for fid in sub_p["family_id"].astype(str).head(10):
                if fid not in non_exp_dark_ctrl_ids and fid not in already_chosen:
                    non_exp_dark_ctrl_ids.append(fid)
                    already_chosen.add(fid)
    if len(non_exp_dark_ctrl_ids) < 120:
        for fid in non_exp_dark_pool.sort_values(["n_gblocks", "n_loci", "family_id"], ascending=[False, False, True])["family_id"].astype(str):
            if fid not in non_exp_dark_ctrl_ids and fid not in already_chosen:
                non_exp_dark_ctrl_ids.append(fid)
                already_chosen.add(fid)
            if len(non_exp_dark_ctrl_ids) >= 120:
                break

    # Assemble frozen family universe records and assign Pilot shard (0) vs Full Array shards (1..16)
    pilot_priority_ids: list[str] = []
    for fid in list(STAGE3_FOCAL_DARK_CANDIDATES) + prior_ids:
        if fid not in pilot_priority_ids:
            pilot_priority_ids.append(fid)
    for fid in r2_elig_ids[:15]:
        if fid not in pilot_priority_ids:
            pilot_priority_ids.append(fid)
    for fid in known_exp_ctrl_ids[:8] + neg_ctrl_ids[:6] + non_exp_dark_ctrl_ids[:6]:
        if fid not in pilot_priority_ids:
            pilot_priority_ids.append(fid)
    pilot_set = set(pilot_priority_ids)

    universe_records: list[dict[str, Any]] = []
    seen_fids: set[str] = set()

    def _add_fam(fid: str, cohort_group: str, cand_role: str, is_dk: bool) -> None:
        if fid in seen_fids:
            return
        seen_fids.add(fid)
        universe_records.append(
            {
                "family_id": fid,
                "cohort_group": cohort_group,
                "candidate_role": cand_role,
                "is_dark": is_dk,
                "in_route_a2_eligible": fid in set(r2_elig_ids),
                "in_route_a1_eligible_dark": fid in set(r1_elig_dark_ids),
            }
        )

    for fid in r2_elig_ids:
        _add_fam(fid, "candidate_route_a2_exported_dark_eligible", "de_novo_exported_dark_candidate", True)
    for fid in r1_elig_dark_ids:
        _add_fam(fid, "candidate_route_a1_gapfill_dark_eligible", "de_novo_gapfill_dark_candidate", True)
    for fid in prior_ids:
        _add_fam(fid, "prior_heat_stress_candidate_example", "prior_heat_stress_candidate_example", True)
    for fid in known_exp_ctrl_ids:
        _add_fam(fid, "control_known_function_exported_background_panel", "annotated_exported_background_comparator_panel", False)
    for fid in neg_ctrl_ids:
        _add_fam(fid, "control_annotation_negative_exported_matched", "matched_annotation_negative_exported_control", True)
    for fid in non_exp_dark_ctrl_ids:
        _add_fam(fid, "control_non_exported_dark_verified_soluble_background_panel", "verified_soluble_non_exported_dark_background_panel", True)

    non_pilot_idx = 0
    for rec in universe_records:
        if rec["family_id"] in pilot_set:
            rec["shard_id"] = 0
        else:
            rec["shard_id"] = 1 + (non_pilot_idx % N_ARRAY_SHARDS)
            non_pilot_idx += 1

    universe_df = pd.DataFrame(universe_records)
    validate_primary_key(universe_df, ["family_id"], "frozen_family_universe.tsv")
    universe_df.to_csv(P_PROV / "frozen_family_universe.tsv", sep="\t", index=False)

    # Write control cohort covariate balance & missingness audit (`EXP4-PIPELINE-REVIEW.txt` Item 4b)
    u_with_stats = universe_df.merge(fam_grp, on="family_id", how="left")
    cov_balance_rows: list[dict[str, Any]] = []
    for c_grp, g_sub in u_with_stats.groupby("cohort_group", sort=True):
        q25, q50, q75 = np.nanpercentile(g_sub["median_len"].to_numpy(dtype=float), [25, 50, 75])
        cov_balance_rows.append(
            {
                "cohort_group": str(c_grp),
                "n_families": int(len(g_sub)),
                "median_aa_length": float(q50),
                "iqr_aa_length": float(q75 - q25),
                "median_genus_blocks": float(g_sub["n_gblocks"].median()),
                "median_member_loci": float(g_sub["n_loci"].median()),
                "n_families_with_any_missing_sp_or_tm": int((g_sub["n_missing_sp_or_tm"] > 0).sum()),
                "frac_families_all_verified_soluble": float(g_sub["all_verified_soluble"].mean()),
                "cohort_comparison_scope": "descriptive_background_comparator_no_causal_enrichment_p_value_claimed",
            }
        )
    pd.DataFrame(cov_balance_rows).to_csv(P_RES / "control_cohort_covariate_balance.tsv", sep="\t", index=False)

    # 5. Audit and verify member loci across all families in the frozen universe
    universe_fam_set = set(universe_df["family_id"].astype(str))
    target_member_rows = genes_thermo[genes_thermo["family_id"].astype(str).isin(universe_fam_set)].copy()
    target_contigs_by_genome = set(
        zip(target_member_rows["genome"].astype(str), target_member_rows["contig"].astype(str))
    )
    contig_genes_df = genes_thermo[
        [
            (str(g), str(c)) in target_contigs_by_genome
            for g, c in zip(genes_thermo["genome"], genes_thermo["contig"])
        ]
    ].copy()

    raw_aa_map = load_raw_archaeal_aa_map(set(contig_genes_df["gene"].astype(str)))
    audited_contig_genes = audit_gene_table_syntax_and_fasta(contig_genes_df, raw_aa_map)
    audited_gene_lookup = audited_contig_genes.set_index("gene")

    member_genomes_set = set(target_member_rows["genome"].astype(str))
    prov_status_map, refined_contigs_map = load_verified_refined_genomes_subset(
        member_genomes_set, P_LEG_SCI / "inputs"
    )

    members_by_fam: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    target_cat_ids: set[str] = set()
    target_gene_ids: set[str] = set()
    universe_faa_path = P_RES / "frozen_universe_members.faa"

    with open(universe_faa_path, "wt") as faa_out:
        for a_row in target_member_rows.itertuples(index=False):
            fid = str(a_row.family_id)
            gid = str(a_row.gene)
            g_mag = str(a_row.genome)
            cat_id = (
                str(a_row.cat_id)
                if pd.notna(a_row.cat_id) and str(a_row.cat_id).strip()
                else f"{g_mag.split('__', 1)[0]}::{gid}"
            )
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
            clean_seq = clean_aa_sequence(raw_seq)
            if clean_seq:
                faa_out.write(f">{gid}\n{clean_seq}\n")

            members_by_fam[fid].append(
                {
                    "family_id": fid,
                    "gene": gid,
                    "cat_id": cat_id,
                    "genome": g_mag,
                    "tax_class": str(a_row.tax_class),
                    "tax_order": str(a_row.tax_order),
                    "genus_singleton_block": str(a_row.genus_singleton_block),
                    "contig": str(g_info["contig"]),
                    "length": int(g_info["length"]),
                    "contig_edge_status": edge_st,
                    "verified_complete_interior_member": verified_interior,
                    "signalp_class": a_row.signalp_class,
                    "signalp_cs_pos": a_row.signalp_cs_pos,
                    "tmbed_S": a_row.tmbed_S,
                    "tmbed_H": a_row.tmbed_H,
                    "tmbed_B": a_row.tmbed_B,
                    "tmbed_h_spans": a_row.tmbed_h_spans,
                    "raw_aa_sequence": clean_seq,
                }
            )

    # 6. Load per-protein raw evidence and run DIAMOND + Foldseek against the curated reference panel
    raw_ev = load_per_protein_raw_evidence(target_cat_ids | universe_fam_set, target_gene_ids)
    seq_hits_by_gene, struct_hits_by_cat = run_curated_panel_diamond_and_foldseek(
        universe_faa_path, raw_ev["sname_to_pids"], threads=threads
    )

    # Write explicit 8CE8 chain resolution audit table for 51Dr3::51Dr3_552805_4
    audit_8ce8_rows: list[dict[str, Any]] = []
    for m in members_by_fam.get("51Dr3::51Dr3_552805_4", []):
        m_cid = str(m["cat_id"])
        m_seq = str(m["raw_aa_sequence"])
        for h in struct_hits_by_cat.get(m_cid, []):
            if "8CE8" in str(h.get("ref_id", "")):
                g_res = resolve_pdb_chain_and_architecture_guards(
                    h["ref_id"], m_seq, h["qtmscore"], h["ttmscore"], h["qcov"], h["tcov"]
                )
                audit_8ce8_rows.append({**h, **g_res})
        for h in raw_ev["pdb_hits"].get((m_cid, "ColabFold_AF2"), []) + raw_ev["pdb_hits"].get((m_cid, "ESMFold"), []):
            if "8ce8" in str(h.get("target", "")).lower():
                g_res = resolve_pdb_chain_and_architecture_guards(
                    h["target"], m_seq, h["qtmscore"], h["ttmscore"], h["qcov"], h["tcov"]
                )
                audit_8ce8_rows.append({"cat_id": m_cid, **h, **g_res})
    pd.DataFrame(audit_8ce8_rows).to_csv(P_RES / "pdb_8ce8_chain_resolution_audit.tsv", sep="\t", index=False)

    # 7. Build flanking annotation lookups for all contigs containing frozen universe members
    prof_hits_path = require_file(P_PROT_INT / "results/profile_gene_hits.tsv.gz", "profile_gene_hits.tsv.gz")
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
                bp = r_c4.get("best_pfam") or r_c4.get("pfam_id") or ""
                pd_d = r_c4.get("pdb_desc") or r_c4.get("pfam_desc") or ""
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
            "signalp_class": cg_row.signalp_class,
            "functional_description": "; ".join(desc_parts) if desc_parts else "hypothetical protein",
        }

    contig_ordered_genes: dict[tuple[str, str], list[str]] = {}
    for (g_mag, c_id), grp in audited_contig_genes.sort_values("contig_rank").groupby(["genome", "contig"], sort=False):
        contig_ordered_genes[(str(g_mag), str(c_id))] = grp["gene"].astype(str).tolist()

    # 8. Write self-contained shard input bundles (`results/shard_inputs/shard_{0..16}.pkl`)
    for s_id in range(0, N_ARRAY_SHARDS + 1):
        s_fams = [r for r in universe_records if int(r["shard_id"]) == s_id]
        s_fids = {str(r["family_id"]) for r in s_fams}
        s_mems = {fid: members_by_fam.get(fid, []) for fid in s_fids}
        s_contigs = {
            (str(m["genome"]), str(m["contig"]))
            for fid in s_fids
            for m in s_mems.get(fid, [])
        }
        s_c_ordered = {k: contig_ordered_genes[k] for k in s_contigs if k in contig_ordered_genes}
        s_c_genes = {g for gl in s_c_ordered.values() for g in gl}
        s_flank_annot = {g: flank_annot_map[g] for g in s_c_genes if g in flank_annot_map}
        s_gene_profs = {g: gene_profiles_map[g] for g in s_c_genes if g in gene_profiles_map}
        s_cats = {str(m["cat_id"]) for fid in s_fids for m in s_mems.get(fid, [])} | s_fids
        s_genes = {str(m["gene"]) for fid in s_fids for m in s_mems.get(fid, [])}

        s_raw_ev = {
            "struct_models": {k: v for k, v in raw_ev["struct_models"].items() if k[0] in s_cats},
            "pdb_hits": {k: v for k, v in raw_ev["pdb_hits"].items() if k[0] in s_cats},
            "hh_hits": {k: v for k, v in raw_ev["hh_hits"].items() if k in s_cats},
            "ipr_domains": {k: v for k, v in raw_ev["ipr_domains"].items() if k in s_cats},
            "ipr_nondomain": {k: v for k, v in raw_ev["ipr_nondomain"].items() if k in s_cats},
            "c6_hits": {k: v for k, v in raw_ev["c6_hits"].items() if k in s_genes},
        }
        s_seq_hits = {g: seq_hits_by_gene[g] for g in s_genes if g in seq_hits_by_gene}
        s_struct_hits = {c: struct_hits_by_cat[c] for c in s_cats if c in struct_hits_by_cat}

        shard_payload = {
            "shard_id": s_id,
            "families_meta": s_fams,
            "members_by_fam": s_mems,
            "contig_ordered_genes": s_c_ordered,
            "flank_annot_map": s_flank_annot,
            "gene_profiles_map": s_gene_profs,
            "raw_ev": s_raw_ev,
            "seq_hits_by_gene": s_seq_hits,
            "struct_hits_by_cat": s_struct_hits,
        }
        with open(P_SHARD_IN / f"shard_{s_id}.pkl", "wb") as fh:
            pickle.dump(shard_payload, fh, protocol=pickle.HIGHEST_PROTOCOL)

    # 9. Run and scientifically validate Pilot shard (`shard_id = 0`)
    pilot_res = evaluate_shard_bundle(0)
    pilot_df = pd.read_csv(P_SHARD_OUT / "shard_0_family_evaluation.tsv", sep="\t")
    if pilot_df.empty or "51Dr3::51Dr3_552805_4" not in set(pilot_df["family_id"].astype(str)):
        raise RuntimeError("Pilot validation failed: focal family 51Dr3::51Dr3_552805_4 missing from shard_0")

    # Verify that 51Dr3::51Dr3_552805_4 does NOT receive automatic CcmE functional transfer
    row_51 = pilot_df[pilot_df["family_id"] == "51Dr3::51Dr3_552805_4"].iloc[0]
    if bool(row_51["pdb_has_full_mechanism_support"]):
        raise RuntimeError("Pilot scientific gate failed: 51Dr3::51Dr3_552805_4 falsely granted full mechanism support!")

    freeze_summary = {
        "panel_version": PANEL_VERSION,
        "n_curated_reference_entries": len(PANEL_SPECS),
        "n_frozen_universe_families": int(len(universe_df)),
        "cohort_counts": {
            str(k): int(v) for k, v in universe_df["cohort_group"].value_counts().to_dict().items()
        },
        "n_pilot_shard0_families": int((universe_df["shard_id"] == 0).sum()),
        "n_array_shards": N_ARRAY_SHARDS,
        "n_member_loci_extracted": int(sum(len(v) for v in members_by_fam.values())),
        "pilot_validation_status": "PASSED_scientific_and_schema_gates",
        "pilot_result": pilot_res,
        "input_sha256": {
            "curated_role_reference_panel.tsv": sha256_file(P_IN / "curated_role_reference_panel.tsv"),
            "curated_role_reference_panel.faa": sha256_file(P_IN / "curated_role_reference_panel.faa"),
            "task_a_stage1b_route2_exported_dark_family_ranking_verified.tsv": sha256_file(r2_path),
            "task_a_stage1b_route1_gapfill_family_ranking_verified.tsv": sha256_file(r1_path),
            "refined_representatives.tsv": sha256_file(canon_rep_path),
        },
    }
    (P_PROV / "frozen_family_universe_summary.json").write_text(json.dumps(freeze_summary, indent=2) + "\n")
    (P_RES / "pilot_validation_summary.json").write_text(json.dumps(freeze_summary, indent=2) + "\n")
    print("[exp4 freeze_and_pilot] Completed successfully:", json.dumps(freeze_summary, indent=2))


def _load_go_ancestors_from_obo(obo_path: Path) -> tuple[dict[str, set[str]], dict[str, str], dict[str, str]]:
    parents: dict[str, set[str]] = collections.defaultdict(set)
    names: dict[str, str] = {}
    namespaces: dict[str, str] = {}
    cur_term: str | None = None
    with open(obo_path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if line == "[Term]":
                cur_term = None
            elif line.startswith("id: GO:"):
                cur_term = line[4:].strip()
            elif cur_term and line.startswith("name: "):
                names[cur_term] = line[6:].strip()
            elif cur_term and line.startswith("namespace: "):
                namespaces[cur_term] = line[11:].strip()
            elif cur_term and line.startswith("is_a: "):
                p_go = line[6:].split()[0].strip()
                parents[cur_term].add(p_go)
            elif cur_term and line.startswith("relationship: part_of GO:"):
                p_go = line.split("part_of", 1)[1].strip().split()[0]
                parents[cur_term].add(p_go)

    memo: dict[str, set[str]] = {}

    def _get_anc(t: str, stack: frozenset[str] = frozenset()) -> set[str]:
        if t in memo:
            return memo[t]
        out: set[str] = set()
        for p in parents.get(t, set()):
            if p not in stack:
                out.add(p)
                out.update(_get_anc(p, stack | {t}))
        memo[t] = out
        return out

    for t in list(names.keys()):
        _get_anc(t)
    return memo, names, namespaces


def run_gspa_evaluation() -> None:
    """Evaluate frozen GSPA 3-genome pilot and check train-time holdout manifest status (`EXP4-CODE-REVIEW.txt` item 9)."""
    P_RES.mkdir(parents=True, exist_ok=True)
    meta_path = require_file(P_GSPA_PILOT / "inputs/gene_metadata.tsv", "GSPA pilot gene_metadata.tsv")
    obo_path = require_file(P_GSPA_PILOT / "assets/go.obo", "GSPA pilot go.obo")
    sum_path = require_file(P_GSPA_PILOT / "results/summary.tsv", "GSPA pilot summary.tsv")
    ctx_path = require_file(P_GSPA_PILOT / "results/context-changes.tsv", "GSPA pilot context-changes.tsv")

    go_ancestors, go_names, go_ns = _load_go_ancestors_from_obo(obo_path)
    gspa_meta = pd.read_csv(meta_path, sep="\t")
    pilot_genomes = sorted(gspa_meta["genome"].astype(str).unique().tolist())

    # Map all identifiers (`gene`, `cat_id`, `seq_cluster_rep`, `cat_cluster_rep`, `example_id`) to pilot `gene` IDs
    id_to_pilot_genes: dict[str, list[str]] = collections.defaultdict(list)
    gene_to_meta: dict[str, dict[str, Any]] = {}
    for r in gspa_meta.to_dict(orient="records"):
        gid = str(r["gene"])
        gene_to_meta[gid] = r
        for col_k in ("gene", "cat_id", "seq_cluster_rep", "cat_cluster_rep", "example_id"):
            val = str(r.get(col_k, "")).strip()
            if val and val != "nan":
                id_to_pilot_genes[val].append(gid)
                if "::" in val:
                    id_to_pilot_genes[val.split("::", 1)[1]].append(gid)

    # Check whether a disjoint train-time holdout manifest exists for pretrained DeepGO-PlusPlus-Light assets
    disjoint_manifest_path = P_GSPA_PILOT / "assets/disjoint_archaeal_train_holdout_manifest.tsv"
    has_disjoint_manifest = disjoint_manifest_path.exists()
    disjoint_sha = sha256_file(disjoint_manifest_path) if has_disjoint_manifest else None

    # Extract independent InterProScan GO annotations (`parts[13]`) on non-dark proteins of the 3 pilot genomes
    go_re = re.compile(r"GO:\d{7}")
    ipr_go_by_gene: dict[str, set[str]] = collections.defaultdict(set)
    for p_str in sorted(glob.glob(str(P_ANNOT / "results/c5_interpro/families/chunk_*.tsv.gz"))):
        with gzip.open(p_str, "rt") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 14:
                    continue
                raw_pid = parts[0]
                matched_gids = id_to_pilot_genes.get(raw_pid, [])
                if not matched_gids and "::" in raw_pid:
                    matched_gids = id_to_pilot_genes.get(raw_pid.split("::", 1)[1], [])
                if not matched_gids:
                    continue
                go_field = parts[13]
                gos_found = [
                    m_go
                    for m_go in go_re.findall(go_field)
                    if m_go not in ("GO:0003674", "GO:0008150", "GO:0005575")
                ]
                if gos_found:
                    for mg in matched_gids:
                        ipr_go_by_gene[mg].update(gos_found)

    # Build ground-truth DataFrame strictly from non-dark proteins with independent InterProScan GO terms (NEVER gapfill demands)
    gt_rows: list[dict[str, Any]] = []
    for gid, r in gene_to_meta.items():
        is_dk = str(r.get("dark", "1.0")) == "1.0"
        true_gos = sorted(ipr_go_by_gene.get(gid, set()))
        if (not is_dk) and true_gos:
            gt_rows.append(
                {
                    "protein_id": gid,
                    "genome": str(r["genome"]),
                    "family_id": str(r.get("seq_cluster_rep") or gid),
                    "tax_order": str(r["genome"]).split("__", 1)[0],
                    "true_go_terms": ";".join(true_gos),
                    "is_synthetic_gapfill_demand": False,
                    "is_hypothesized_untested_role": False,
                }
            )
    gt_df = pd.DataFrame(gt_rows)
    gt_protein_set = set(gt_df["protein_id"].astype(str)) if not gt_df.empty else set()

    # Load bounded integrated claims across baseline, context, and shuffled modes for both fast and cnn variants
    claim_rows: list[dict[str, Any]] = []
    dark_shift_rows: list[dict[str, Any]] = []

    for g in pilot_genomes:
        for variant in ("fast", "cnn"):
            sub_dir = P_GSPA_PILOT / "results" / (g if variant == "fast" else f"{g}-cnn")
            base_map: dict[tuple[str, str], float] = {}
            ctx_map: dict[tuple[str, str], float] = {}
            shuf_map: dict[tuple[str, str], float] = {}

            for mode, target_map in (("baseline", base_map), ("context", ctx_map), ("shuffled", shuf_map)):
                p_tsv = sub_dir / f"bounded-{mode}.tsv"
                if not p_tsv.exists():
                    continue
                with open(p_tsv, "rt") as fh:
                    hdr = fh.readline().rstrip("\n").split("\t")
                    for line in fh:
                        row = dict(zip(hdr, line.rstrip("\n").split("\t")))
                        pid = str(row.get("protein_id", ""))
                        fid_go = str(row.get("function_id") or row.get("term") or "")
                        sc = float(row.get("posterior_prob") or row.get("score") or 0.0)
                        if not pid or not fid_go or fid_go in ("GO:0003674", "GO:0008150", "GO:0005575"):
                            continue
                        target_map[(pid, fid_go)] = sc
                        if pid in gt_protein_set:
                            claim_rows.append(
                                {
                                    "protein_id": pid,
                                    "arm": f"{variant}_{mode}",
                                    "predicted_go_term": fid_go,
                                    "score": sc,
                                    "is_leaked_self_or_same_family_hit": False,
                                }
                            )

            # Audit context vs baseline vs shuffled shifts on dark exported proteins
            all_keys = set(base_map.keys()) | set(ctx_map.keys()) | set(shuf_map.keys())
            for pid, fid_go in all_keys:
                m_info = gene_to_meta.get(pid, {})
                if str(m_info.get("dark", "0.0")) != "1.0":
                    continue
                b_sc = base_map.get((pid, fid_go), 0.0)
                c_sc = ctx_map.get((pid, fid_go), 0.0)
                s_sc = shuf_map.get((pid, fid_go), 0.0)
                if abs(c_sc - b_sc) > 1e-6 or abs(s_sc - b_sc) > 1e-6:
                    dark_shift_rows.append(
                        {
                            "genome": g,
                            "variant": variant,
                            "protein_id": pid,
                            "family_id": str(m_info.get("seq_cluster_rep", "")),
                            "go_term": fid_go,
                            "go_name": go_names.get(fid_go, ""),
                            "go_namespace": go_ns.get(fid_go, ""),
                            "score_baseline_sequence_only": b_sc,
                            "score_real_context": c_sc,
                            "score_shuffled_context": s_sc,
                            "delta_real_minus_baseline": c_sc - b_sc,
                            "delta_real_minus_shuffled": c_sc - s_sc,
                            "is_calibrated_heldout": False,
                            "calibration_status": "BLOCKED_unverified_pretrained_train_split_and_graph_leakage_exploratory_only",
                        }
                    )

    claims_df = pd.DataFrame(claim_rows)
    if not gt_df.empty and not claims_df.empty:
        gt_df["pooled_all"] = "ALL_3_PILOT_GENOMES_ESTABLISHED_NON_DARK"
        eval_pooled = evaluate_gspa_heldout_recovery_and_leakage_guard(
            claims_df,
            gt_df,
            holdout_col="pooled_all",
            score_threshold=0.50,
            training_split_manifest_verified=has_disjoint_manifest,
            disjoint_train_manifest_sha256=disjoint_sha,
            go_ancestors_map=go_ancestors,
        )
        eval_by_genome = evaluate_gspa_heldout_recovery_and_leakage_guard(
            claims_df,
            gt_df,
            holdout_col="genome",
            score_threshold=0.50,
            training_split_manifest_verified=has_disjoint_manifest,
            disjoint_train_manifest_sha256=disjoint_sha,
            go_ancestors_map=go_ancestors,
        )
        cal_df = pd.concat([eval_pooled, eval_by_genome], ignore_index=True)
    else:
        cal_df = pd.DataFrame()

    cal_df.to_csv(P_RES / "gspa_calibration_and_context_ablation_audit.tsv", sep="\t", index=False)
    pd.DataFrame(dark_shift_rows).to_csv(P_RES / "gspa_exported_dark_family_context_shifts.tsv", sep="\t", index=False)

    gspa_summary = {
        "stage": "gspa_eval",
        "gspa_execution_scope": "re_evaluation_of_frozen_3_genome_pilot_outputs_not_new_model_training_or_full_cohort_inference",
        "is_calibrated_heldout": False,
        "calibration_status": "BLOCKED_unverified_pretrained_train_split_and_graph_leakage_exploratory_only",
        "blocker_explanation": (
            "Frozen DeepGO-PlusPlus-Light assets (cnn_mcm.pt, train_db.dmnd, train_net_index.tsv, train_terms.tsv) "
            "in /scratch/project/k10018/eq-proteome-gspa-pilot-20261105/assets are pretrained models without a disjoint "
            "family/order train-time holdout manifest for the EQ archaeal families. Per EXP4-CODE-REVIEW.txt Item 9, "
            "post-hoc grouping and output filtering cannot eliminate train-time or STRING graph leakage; therefore "
            "strict held-out calibration is marked BLOCKED and all GSPA context rankings are reported strictly as "
            "uncalibrated exploratory contextual priors."
        ),
        "n_pilot_genomes": len(pilot_genomes),
        "n_established_non_dark_proteins_with_interpro_go": int(len(gt_df)),
        "n_calibration_audit_rows": int(len(cal_df)),
        "n_dark_context_shift_rows": int(len(dark_shift_rows)),
        "pilot_summary_tsv_sha256": sha256_file(sum_path),
        "pilot_context_changes_tsv_sha256": sha256_file(ctx_path),
    }
    (P_RES / "gspa_evaluation_summary.json").write_text(json.dumps(gspa_summary, indent=2) + "\n")
    print("[exp4 gspa_eval] Completed:", json.dumps(gspa_summary, indent=2))


def run_collect() -> None:
    """Merge all 17 shards (`0..16`), apply global BH FDR across all tested `(family_id, guild)` pairs FIRST,
    then compute post-global-FDR convergence vs exploratory ranking (`EXP4-PIPELINE-REVIEW.txt` Items 5 & 6).
    """
    eval_dfs: list[pd.DataFrame] = []
    comp_dfs: list[pd.DataFrame] = []
    neigh_dfs: list[pd.DataFrame] = []
    struct_audit_dfs: list[pd.DataFrame] = []

    for s_id in range(0, N_ARRAY_SHARDS + 1):
        ep = require_file(P_SHARD_OUT / f"shard_{s_id}_family_evaluation.tsv", f"shard_{s_id}_family_evaluation.tsv")
        cp = require_file(P_SHARD_OUT / f"shard_{s_id}_competing_roles.tsv", f"shard_{s_id}_competing_roles.tsv")
        np_p = P_SHARD_OUT / f"shard_{s_id}_neighbourhoods.tsv"
        sa_p = P_SHARD_OUT / f"shard_{s_id}_struct_hits_audit.tsv"
        eval_dfs.append(pd.read_csv(ep, sep="\t"))
        comp_dfs.append(pd.read_csv(cp, sep="\t"))
        if np_p.exists() and np_p.stat().st_size > 10:
            try:
                neigh_dfs.append(pd.read_csv(np_p, sep="\t"))
            except pd.errors.EmptyDataError:
                pass
        if sa_p.exists() and sa_p.stat().st_size > 10:
            try:
                struct_audit_dfs.append(pd.read_csv(sa_p, sep="\t"))
            except pd.errors.EmptyDataError:
                pass

    all_eval = pd.concat(eval_dfs, ignore_index=True)
    all_comp = pd.concat(comp_dfs, ignore_index=True)
    all_neigh = pd.concat(neigh_dfs, ignore_index=True) if neigh_dfs else pd.DataFrame()
    all_struct_audit = pd.concat(struct_audit_dfs, ignore_index=True) if struct_audit_dfs else pd.DataFrame()

    validate_primary_key(all_eval, ["family_id"], "exp4_surface_family_evaluation_all.tsv")

    # Apply global Benjamini-Hochberg FDR across all tested families and across the full (family_id x 6 guilds) matrix
    # Per EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt Item 3: primary multiplicity family is all 1460 x 6 = 8760 joint tests
    null_guilds = [
        "nutrient_acquisition_peptidase",
        "substrate_binding_or_uptake",
        "surface_assembly_pilin_archaellin",
        "cell_coat_or_envelope",
        "host_defense_antiviral",
        "viral_or_mobile_element_origin",
    ]

    # Per-guild sensitivity FDR across families for each guild
    for g_name in null_guilds:
        p_col = f"null_p_exact_{g_name}" if f"null_p_exact_{g_name}" in all_eval.columns else f"null_p_emp_{g_name}"
        if p_col in all_eval.columns:
            all_eval[f"global_q_bh_{g_name}"] = apply_bh_fdr(all_eval[p_col].to_numpy(dtype=float))

    # Primary joint matrix FDR across all 1460 x 6 tests
    p_cols = [
        (f"null_p_exact_{g_name}" if f"null_p_exact_{g_name}" in all_eval.columns else f"null_p_emp_{g_name}")
        for g_name in null_guilds
        if (f"null_p_exact_{g_name}" in all_eval.columns or f"null_p_emp_{g_name}" in all_eval.columns)
    ]
    if len(p_cols) == len(null_guilds):
        p_matrix = all_eval[p_cols].to_numpy(dtype=float)
        q_matrix_joint = apply_bh_fdr(p_matrix.ravel()).reshape(p_matrix.shape)
        for idx_g, g_name in enumerate(null_guilds):
            all_eval[f"joint_matrix_q_bh_{g_name}"] = q_matrix_joint[:, idx_g]

    # Populate legacy aliases for s_layer_or_charted_coat
    all_eval["global_q_bh_s_layer_or_charted_coat"] = all_eval["global_q_bh_cell_coat_or_envelope"]
    all_eval["joint_matrix_q_bh_s_layer_or_charted_coat"] = all_eval["joint_matrix_q_bh_cell_coat_or_envelope"]

    # Re-evaluate same-role convergence and evidence routes (EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt Items 4 & 5)
    post_fdr_guilds_list: list[str] = []
    n_post_fdr_guilds_list: list[int] = []
    strict_convergent_list: list[bool] = []
    post_fdr_status_list: list[str] = []
    supported_routes_list: list[str] = []
    hypotheses_list: list[str] = []
    deep_dive_rows: list[dict[str, Any]] = []

    for row in all_eval.to_dict(orient="records"):
        fid = str(row["family_id"])
        is_dk = bool(row.get("is_dark", False))
        sig_guilds: list[str] = []
        for g_name in null_guilds:
            obs_adj = int(row.get(f"null_obs_adj_{g_name}", 0))
            n_gb = int(row.get(f"null_n_gblocks_{g_name}", 0))
            joint_q = float(row.get(f"joint_matrix_q_bh_{g_name}", 1.0))
            if obs_adj >= 2 and n_gb >= 2 and joint_q <= 0.10:
                sig_guilds.append(g_name)
        post_fdr_guilds_list.append(";".join(sig_guilds) if sig_guilds else "none")
        n_post_fdr_guilds_list.append(len(sig_guilds))

        # Evaluate role-specific convergence across sequence, structure, architecture, and context for the SAME role
        conv_res = evaluate_same_role_evidence_and_convergence(row, fdr_threshold=0.10, use_joint_matrix_q=True)
        supp_route = conv_res["supported_evidence_route"]
        cand_hyp = conv_res["candidate_computational_role_hypothesis"]
        supported_routes_list.append(supp_route)
        hypotheses_list.append(cand_hyp)

        is_strict_conv = bool(
            is_dk
            and supp_route in (
                "convergent_sequence_structure_and_same_role_context",
                "qualified_structure_and_same_role_context_with_architecture",
                "sequence_and_same_role_context_with_architecture",
            )
        )
        strict_convergent_list.append(is_strict_conv)

        if is_strict_conv:
            post_fdr_status_list.append(f"{supp_route}:{cand_hyp}")
        elif supp_route == "domain_profile_and_same_role_context_with_architecture":
            post_fdr_status_list.append(f"domain_profile_and_same_role_context_with_architecture:{cand_hyp}")
        elif supp_route == "exploratory_fold_or_motif_candidate_no_operon_convergence":
            post_fdr_status_list.append("exploratory_fold_or_motif_candidate_no_global_fdr_operon_convergence")
        else:
            post_fdr_status_list.append("uncertain_insufficient_evidence")

        # Check if this family belongs in the comprehensive deep dive table (Item C)
        has_small_p = any(
            float(row.get(f"null_p_exact_{g}", row.get(f"null_p_emp_{g}", 1.0))) <= 0.01
            for g in null_guilds
        )
        is_exploratory_dark = fid in EXPLORATORY_DARK_CANDIDATES_V1
        is_predeclared_p = fid in PREDECLARED_SMALL_P_FAMILIES_V1
        is_domain_conv = (supp_route == "domain_profile_and_same_role_context_with_architecture")

        if is_strict_conv or is_domain_conv or is_exploratory_dark or is_predeclared_p or has_small_p:
            selection_bases: list[str] = []
            if is_strict_conv:
                selection_bases.append("strict_convergent_candidate_v2")
            if is_domain_conv:
                selection_bases.append("domain_profile_context_candidate_v2")
            if is_exploratory_dark:
                selection_bases.append("v1_exploratory_dark_candidate")
            if is_predeclared_p:
                selection_bases.append("v1_predeclared_context_p_le_0.01")
            if has_small_p:
                selection_bases.append("v2_exact_context_p_le_0.01")

            deep_dive_rows.append(
                {
                    "family_id": fid,
                    "cohort_group": str(row["cohort_group"]),
                    "is_dark": is_dk,
                    "selection_basis": ";".join(selection_bases),
                    "representative_gene": str(row["representative_gene"]),
                    "representative_genome": str(row["representative_genome"]),
                    "representative_cat_id": str(row["representative_cat_id"]),
                    "aa_length": row["aa_length"],
                    "surface_localization_class": str(row["surface_localization_class"]),
                    "topology_status": str(row["topology_status"]),
                    "model_source": str(row.get("model_source", "none")),
                    "model_polymer_hash_status": str(row.get("model_polymer_hash_status", "missing_model")),
                    "model_plddt_0to1": row.get("model_plddt_0to1", np.nan),
                    "best_curated_seq_ref_id": str(row.get("best_curated_seq_ref_id", "none")),
                    "best_curated_seq_pident": row.get("best_curated_seq_pident", np.nan),
                    "best_curated_seq_qcov": row.get("best_curated_seq_qcov", np.nan),
                    "best_curated_seq_tcov": row.get("best_curated_seq_tcov", np.nan),
                    "best_curated_seq_evalue": row.get("best_curated_seq_evalue", np.nan),
                    "best_curated_seq_bits": row.get("best_curated_seq_bits", np.nan),
                    "best_curated_struct_ref_id": str(row.get("best_curated_struct_ref_id", "none")),
                    "best_curated_struct_qtmscore": row.get("best_curated_struct_qtmscore", np.nan),
                    "best_curated_struct_ttmscore": row.get("best_curated_struct_ttmscore", np.nan),
                    "best_curated_struct_qcov": row.get("best_curated_struct_qcov", np.nan),
                    "best_curated_struct_tcov": row.get("best_curated_struct_tcov", np.nan),
                    "pdb_guard_status": str(row.get("pdb_guard_status", "none")),
                    "pdb_resolved_ref_id": str(row.get("pdb_resolved_ref_id", "none")),
                    **{f"null_obs_adj_{g}": row.get(f"null_obs_adj_{g}", 0) for g in null_guilds},
                    **{f"null_p_exact_{g}": row.get(f"null_p_exact_{g}", row.get(f"null_p_emp_{g}", 1.0)) for g in null_guilds},
                    **{f"genus_derep_obs_adj_{g}": row.get(f"genus_derep_obs_adj_{g}", 0) for g in null_guilds},
                    **{f"genus_derep_p_exact_{g}": row.get(f"genus_derep_p_exact_{g}", np.nan) for g in null_guilds},
                    **{f"joint_matrix_q_bh_{g}": row.get(f"joint_matrix_q_bh_{g}", 1.0) for g in null_guilds},
                    **{f"global_q_bh_{g}": row.get(f"global_q_bh_{g}", 1.0) for g in null_guilds},
                    "msa_status": str(row.get("msa_status", "none")),
                    "n_aligned_sequences": row.get("n_aligned_sequences", 0),
                    "median_pairwise_aa_identity": row.get("median_pairwise_aa_identity", np.nan),
                    "tree_status": str(row.get("tree_status", "none")),
                    "supported_evidence_route": supp_route,
                    "candidate_computational_role_hypothesis": cand_hyp,
                    "is_computational_role_hypothesis": True,
                    "is_biochemical_assignment": False,
                    "scientific_qualification_note": (
                        "Computational role hypothesis; domain annotations and structural matches are predictors, "
                        "not experimental biochemical validation."
                    ),
                }
            )

    all_eval["post_global_fdr_significant_guilds"] = post_fdr_guilds_list
    all_eval["n_post_global_fdr_significant_guilds"] = n_post_fdr_guilds_list
    all_eval["strict_convergent_candidate_awaiting_review"] = strict_convergent_list
    all_eval["post_global_fdr_discrimination_status"] = post_fdr_status_list
    all_eval["supported_evidence_route"] = supported_routes_list
    all_eval["candidate_computational_role_hypothesis"] = hypotheses_list
    all_eval["is_computational_role_hypothesis"] = True
    all_eval["is_biochemical_assignment"] = False

    # Write versioned v2 tables and maintain primary _all tables
    all_eval.to_csv(P_RES / "exp4_surface_family_evaluation_v2.tsv", sep="\t", index=False)
    all_eval.to_csv(P_RES / "exp4_surface_family_evaluation_all.tsv", sep="\t", index=False)
    all_comp.to_csv(P_RES / "exp4_surface_competing_role_matrix_v2.tsv", sep="\t", index=False)
    all_comp.to_csv(P_RES / "exp4_surface_competing_role_matrix_all.tsv", sep="\t", index=False)
    all_neigh.to_csv(P_RES / "exp4_surface_neighbourhood_panels_v2.tsv", sep="\t", index=False)
    all_neigh.to_csv(P_RES / "exp4_surface_neighbourhood_panels_all.tsv", sep="\t", index=False)
    all_struct_audit.to_csv(P_RES / "exp4_surface_struct_hits_audit_v2.tsv", sep="\t", index=False)
    all_struct_audit.to_csv(P_RES / "exp4_surface_struct_hits_audit_all.tsv", sep="\t", index=False)

    deep_dive_df = pd.DataFrame(deep_dive_rows)
    deep_dive_df.to_csv(P_RES / "exp4_exploratory_and_candidate_deep_dive_v2.tsv", sep="\t", index=False)

    # Summarise surface localization classes, evidence routes, and model availability across cohort groups
    cohort_rows: list[dict[str, Any]] = []
    for c_grp, grp_df in all_eval.groupby("cohort_group", sort=True):
        n_fams = len(grp_df)
        loc_counts = grp_df["surface_localization_class"].value_counts().to_dict()
        disc_counts = grp_df["post_global_fdr_discrimination_status"].value_counts().to_dict()
        route_counts = grp_df["supported_evidence_route"].value_counts().to_dict()
        model_hash_counts = grp_df["model_polymer_hash_status"].value_counts().to_dict()
        cohort_rows.append(
            {
                "cohort_group": str(c_grp),
                "n_families": int(n_fams),
                "n_models_exact_match": int((grp_df["model_polymer_hash_status"].isin(["exact_match", "exact_sequence_hash_match"])).sum()),
                "n_models_mismatch": int((grp_df["model_polymer_hash_status"].isin(["mismatch", "sequence_mismatch"])).sum()),
                "n_models_missing": int((grp_df["model_polymer_hash_status"] == "missing_model").sum()),
                "n_freely_secreted_rep": int((grp_df["surface_localization_class"] == "freely_secreted_extracellular_candidate").sum()),
                "n_lipoprotein_spII_rep": int((grp_df["surface_localization_class"] == "lipoprotein_membrane_anchored_SPII").sum()),
                "n_prepilin_archaellin_rep": int((grp_df["surface_localization_class"] == "prepilin_or_archaellin_filament_subunit_candidate").sum()),
                "n_single_pass_anchored_rep": int(
                    grp_df["surface_localization_class"].isin(
                        ["n_terminal_anchored_single_pass_surface", "c_terminal_anchored_surface_coat_candidate", "internal_single_pass_membrane_with_export_signal"]
                    ).sum()
                ),
                "n_multipass_or_beta_barrel_rep": int(
                    (grp_df["surface_localization_class"] == "multipass_or_beta_barrel_membrane_with_export_signal").sum()
                ),
                "n_missing_tmbed_unknown_topology": int((grp_df["topology_status"] == "unknown_missing_tmbed").sum()),
                "median_aa_length": float(grp_df["aa_length"].median()),
                "median_predicted_nxst_per_100aa": float(grp_df["predicted_nxst_sequons_per_100aa"].median()),
                "n_strict_post_global_fdr_convergent_candidates_awaiting_review": int(
                    grp_df["strict_convergent_candidate_awaiting_review"].astype(bool).sum()
                ),
                "n_route_convergent_sequence_structure_and_context": int(
                    (grp_df["supported_evidence_route"] == "convergent_sequence_structure_and_same_role_context").sum()
                ),
                "n_route_qualified_structure_and_context": int(
                    (grp_df["supported_evidence_route"] == "qualified_structure_and_same_role_context_with_architecture").sum()
                ),
                "n_route_sequence_and_context": int(
                    (grp_df["supported_evidence_route"] == "sequence_and_same_role_context_with_architecture").sum()
                ),
                "n_route_exploratory_fold_or_motif_no_context": int(
                    (grp_df["supported_evidence_route"] == "exploratory_fold_or_motif_candidate_no_operon_convergence").sum()
                ),
                "n_route_uncertain_insufficient_evidence": int(
                    (grp_df["supported_evidence_route"] == "uncertain_insufficient_evidence").sum()
                ),
                "cohort_comparison_scope": "descriptive_background_comparator_no_causal_enrichment_p_value_claimed",
                "localization_breakdown_json": json.dumps(loc_counts),
                "post_global_fdr_discrimination_breakdown_json": json.dumps(disc_counts),
                "supported_evidence_routes_breakdown_json": json.dumps(route_counts),
                "model_polymer_hash_breakdown_json": json.dumps(model_hash_counts),
            }
        )
    pd.DataFrame(cohort_rows).to_csv(P_RES / "exp4_cohort_architecture_and_role_comparison_v2.tsv", sep="\t", index=False)
    pd.DataFrame(cohort_rows).to_csv(P_RES / "exp4_cohort_architecture_and_role_comparison.tsv", sep="\t", index=False)

    # Evaluate Stage 5 conditional detailed modeling gate (strictly requires same-role post-global-FDR convergence)
    stage5_gate_df = all_eval[all_eval["strict_convergent_candidate_awaiting_review"].astype(bool)].copy()
    stage5_gate_df.to_csv(P_RES / "exp4_convergent_candidates_v2.tsv", sep="\t", index=False)
    stage5_gate_df.to_csv(P_RES / "exp4_convergent_candidates_for_stage5_gate.tsv", sep="\t", index=False)

    final_summary = {
        "experiment": "exp4_surface",
        "pipeline_version": "v2",
        "status": "COMPLETED",
        "panel_version": PANEL_VERSION,
        "n_total_families_evaluated": int(len(all_eval)),
        "n_exploratory_and_candidate_deep_dive_families": int(len(deep_dive_df)),
        "cohort_summary": cohort_rows,
        "n_dark_families_passing_strict_post_global_fdr_stage5_gate": int(len(stage5_gate_df)),
        "stage5_detailed_modeling_status": "awaiting_review_of_convergent_evidence_no_docking_launched",
        "stage5_detailed_modeling_disposition": (
            f"{len(stage5_gate_df)} dark families passed strict post-global-FDR multi-axis same-role convergence gate; "
            "stage5_detailed_modeling_status = awaiting_review_of_convergent_evidence_no_docking_launched."
        ),
        "primary_multiplicity_correction": "joint_matrix_q_bh_across_all_1460x6_tests",
        "sensitivity_multiplicity_correction": "global_q_bh_per_guild_across_1460_tests",
        "null_calculation_method": "exact_conditional_circular_shift_discrete_convolution_no_monte_carlo_floor",
        "genus_replication_sensitivity": "deterministic_selection_highest_completeness_lowest_contamination_per_genus_block",
        "gspa_execution_scope": "re_evaluation_of_frozen_3_genome_pilot_outputs_not_new_model_training_or_full_cohort_inference",
        "gspa_calibration_status": "BLOCKED_unverified_pretrained_train_split_and_graph_leakage_exploratory_only",
    }
    (P_RES / "exp4_surface_final_summary_v2.json").write_text(json.dumps(final_summary, indent=2) + "\n")
    (P_RES / "exp4_surface_final_summary.json").write_text(json.dumps(final_summary, indent=2) + "\n")
    print("[exp4 collect v2] Completed:", json.dumps(final_summary, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser(description="Experiment 4 Surface Dark Protein Pipeline")
    ap.add_argument(
        "--mode",
        required=True,
        choices=["freeze_and_pilot", "run_shard", "gspa_eval", "collect"],
    )
    ap.add_argument("--shard-id", type=int, default=0)
    ap.add_argument("--threads", type=int, default=64)
    args = ap.parse_args()

    if args.mode == "freeze_and_pilot":
        run_freeze_and_pilot(threads=args.threads)
    elif args.mode == "run_shard":
        res = evaluate_shard_bundle(args.shard_id)
        print(f"[exp4 run_shard {args.shard_id}] Completed:", json.dumps(res))
    elif args.mode == "gspa_eval":
        run_gspa_evaluation()
    elif args.mode == "collect":
        run_collect()


if __name__ == "__main__":
    main()
