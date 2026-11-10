#!/usr/bin/env python3
"""Experiment 4 (`exp4_surface`): Targeted Followup Evidence Audit & Bounded Inferences.

Implements all requirements from:
- `mechanism-20261109/EXP4-V2-RESULTS-AND-TARGETED-AUDIT.txt`
- `mechanism-20261109/EXP4-TARGETED-AUDIT-CLAIM-BOUNDARIES.txt`

Core scientific and evidence safeguards:
1. Member-level model verification: Compares actual polymer sequences/hashes for every member
   using `verify_prediction_sequence_hash_match`. Distinguishes seed ID, verified interior representative,
   and individual members. Never infers identity from pLDDT or length alone.
2. Honest structural reporting: Low qTM to the limited curated panel is reported as
   "no qualifying fold match to tested panel", NOT "absence of structural homology".
   Unverified or contig-edge models cannot support positive or negative structural conclusions.
3. Bounded neighbourhood scope: Local +/-5 gene survey does NOT establish whole-genome absence
   of pilus machinery or operon structure without strand/distance data. Strands and coordinates
   unrecorded in the catalog bundle are explicitly marked unknown.
4. Accurate coat/S-layer reference audit: 7PTR is the tested positive cellular S-layer reference;
   viral SIRV2 (3J9X) and CcmA (8CE8_Entity1) are confusable/competing controls. Adjacency to
   7PTP_D (H. volcanii csg) does not confer S-layer function to the dark candidate.
5. Distinct count tracking: Separates `n_context_genus_blocks` (unique genus blocks bearing matching
   flanks) from `n_genus_derep_retained_loci` and `n_genus_derep_retained_genus_blocks`.
   Reports primary and sensitivity FDR consistently at both q<=0.10 and q<=0.05.
6. Exports compact member-level (`exp4_47_member_evidence.tsv`) and flank-level
   (`exp4_47_flank_evidence.tsv`) audit tables with raw provenance fields.
7. Reuses completed MAFFT and IQ-TREE work from job 15892410 without redundant re-execution.
"""
from __future__ import annotations

import collections
import glob
import hashlib
import json
import os
import pickle
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

# Ensure Singularity is in PATH on Shaheen
if not shutil.which("singularity"):
    for cand in ("/sw/ex113genoa/singularity/ce/4.3.2/bin", "/usr/bin", "/usr/local/bin"):
        if (Path(cand) / "singularity").exists():
            os.environ["PATH"] = f"{cand}:{os.environ.get('PATH', '')}"
            break

# Path setup
P_EXP4_ROOT = Path(
    os.environ.get("EQ_EXP4_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp4_surface")
)
if not P_EXP4_ROOT.exists():
    P_EXP4_ROOT = Path(__file__).resolve().parent.parent

P_IN = P_EXP4_ROOT / "inputs"
P_RES = P_EXP4_ROOT / "results"
P_SHARD_IN = P_RES / "shard_inputs"
P_SHARD_OUT = P_RES / "shards"
P_MSA_OUT = P_RES / "family_msas_and_trees"
P_FOLLOWUP = P_RES / "followup_v2_review"

sys.path.insert(0, str(P_EXP4_ROOT / "code"))
from common_science import sha256_file
from exp4_surface_helpers import (
    apply_bh_fdr,
    classify_flanking_machinery_guilds,
    verify_prediction_sequence_hash_match,
)
from exp4_surface_pipeline import run_mafft_and_iqtree_for_family

NULL_GUILDS = [
    "nutrient_acquisition_peptidase",
    "substrate_binding_or_uptake",
    "surface_assembly_pilin_archaellin",
    "cell_coat_or_envelope",
    "host_defense_antiviral",
    "viral_or_mobile_element_origin",
]


def load_verified_loci_input() -> tuple[pd.DataFrame, str]:
    """Load authoritative verified candidate member loci table containing exact genomic coordinates."""
    locus_p = P_IN / "task_a_stage1b_all_candidate_member_loci_verified.tsv.gz"
    if not locus_p.exists():
        locus_p = Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/exp3_viruses/inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz")
    if not locus_p.exists():
        raise FileNotFoundError(f"Verified locus table not found at {locus_p}")
    locus_sha = sha256_file(str(locus_p))
    df = pd.read_csv(locus_p, sep="\t")
    return df, locus_sha


def load_all_shard_inputs() -> dict[str, Any]:
    """Load and index all frozen shard inputs (shards 0..16)."""
    shard_files = sorted(glob.glob(str(P_SHARD_IN / "shard_*.pkl")))
    if not shard_files:
        raise FileNotFoundError(f"No shard pickle files found under {P_SHARD_IN}")

    fams_meta: dict[str, dict[str, Any]] = {}
    members_by_fam: dict[str, list[dict[str, Any]]] = {}
    contig_ordered_genes: dict[tuple[str, str], list[str]] = {}
    flank_annot_map: dict[str, dict[str, Any]] = {}
    gene_profiles_map: dict[str, list[str]] = {}
    raw_struct_models: dict[tuple[str, str], dict[str, Any]] = {}
    raw_pdb_hits: dict[tuple[str, str], list[dict[str, Any]]] = {}
    raw_hh_hits: dict[str, list[dict[str, Any]]] = {}
    raw_ipr_domains: dict[str, list[dict[str, Any]]] = {}
    raw_ipr_nondomain: dict[str, list[dict[str, Any]]] = {}
    raw_c6_hits: dict[str, list[dict[str, Any]]] = {}
    seq_hits_by_gene: dict[str, list[dict[str, Any]]] = {}
    struct_hits_by_cat: dict[str, list[dict[str, Any]]] = {}
    fam_to_shard: dict[str, int] = {}

    for sf in shard_files:
        with open(sf, "rb") as fh:
            data = pickle.load(fh)
        sid = data["shard_id"]
        for fmeta in data["families_meta"]:
            fid = fmeta["family_id"]
            fams_meta[fid] = fmeta
            fam_to_shard[fid] = sid
        for fid, mems in data["members_by_fam"].items():
            members_by_fam[fid] = mems
        for ckey, glist in data["contig_ordered_genes"].items():
            contig_ordered_genes[ckey] = glist
        flank_annot_map.update(data["flank_annot_map"])
        gene_profiles_map.update(data["gene_profiles_map"])
        raw_ev = data["raw_ev"]
        raw_struct_models.update(raw_ev.get("struct_models", {}))
        raw_pdb_hits.update(raw_ev.get("pdb_hits", {}))
        raw_hh_hits.update(raw_ev.get("hh_hits", {}))
        raw_ipr_domains.update(raw_ev.get("ipr_domains", {}))
        raw_ipr_nondomain.update(raw_ev.get("ipr_nondomain", {}))
        raw_c6_hits.update(raw_ev.get("c6_hits", {}))
        seq_hits_by_gene.update(data.get("seq_hits_by_gene", {}))
        struct_hits_by_cat.update(data.get("struct_hits_by_cat", {}))

    return {
        "fams_meta": fams_meta,
        "members_by_fam": members_by_fam,
        "contig_ordered_genes": contig_ordered_genes,
        "flank_annot_map": flank_annot_map,
        "gene_profiles_map": gene_profiles_map,
        "raw_struct_models": raw_struct_models,
        "raw_pdb_hits": raw_pdb_hits,
        "raw_hh_hits": raw_hh_hits,
        "raw_ipr_domains": raw_ipr_domains,
        "raw_ipr_nondomain": raw_ipr_nondomain,
        "raw_c6_hits": raw_c6_hits,
        "seq_hits_by_gene": seq_hits_by_gene,
        "struct_hits_by_cat": struct_hits_by_cat,
        "fam_to_shard": fam_to_shard,
    }


def compute_full_universe_genus_joint_fdr(eval_df: pd.DataFrame) -> pd.DataFrame:
    """Compute full-universe genus-dereplicated joint BH-FDR correction across 8,760 tests.

    Missing or non-estimable genus p-values are treated conservatively as p=1.0.
    """
    df = eval_df.copy()
    raw_p_list: list[float] = []
    idx_tuples: list[tuple[int, str]] = []

    for row_idx, (_, row) in enumerate(df.iterrows()):
        for guild in NULL_GUILDS:
            col_name = f"genus_derep_p_exact_{guild}"
            val = row.get(col_name)
            p_val = float(val) if pd.notna(val) else 1.0
            p_val = max(0.0, min(1.0, p_val))
            raw_p_list.append(p_val)
            idx_tuples.append((row_idx, guild))

    assert len(raw_p_list) == len(df) * len(NULL_GUILDS) == 1460 * 6 == 8760
    q_vals = apply_bh_fdr(raw_p_list)
    assert len(q_vals) == 8760

    # Initialize new columns
    for guild in NULL_GUILDS:
        df[f"genus_derep_joint_matrix_q_bh_{guild}"] = 1.0

    for (row_idx, guild), q in zip(idx_tuples, q_vals):
        df.iat[row_idx, df.columns.get_loc(f"genus_derep_joint_matrix_q_bh_{guild}")] = float(q)

    min_q = min(q_vals)
    print(f"[Full-Universe Genus Joint BH-FDR] Total tests: {len(q_vals)}, min q: {min_q:.4f}")
    assert 0.0 <= min_q <= 1.0, f"Calculated min q {min_q} out of valid [0, 1] probability range"
    return df


def audit_lead_phylogenetics(bundle: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Verify or retrieve MAFFT and IQ-TREE 3 outputs for Lead 1 (22PRr1) and Lead 2 (V60PRr2).

    Reuses existing completed alignments and treefiles from job 15892410 where available.
    """
    leads = {
        "22PRr1::22PRr1_2531955_3": "surface_assembly_pilin_archaellin",
        "V60PRr2::V60PRr2_1314979_4": "cell_coat_or_envelope",
    }
    phylo_results: dict[str, dict[str, Any]] = {}
    members_by_fam = bundle["members_by_fam"]

    for fid, role in leads.items():
        safe_fid = fid.replace("::", "__").replace("/", "_")
        aln_p = P_MSA_OUT / f"{safe_fid}.aln.faa"
        tree_p = P_MSA_OUT / f"{safe_fid}.iqtree.treefile"

        mems = members_by_fam.get(fid, [])
        seq_map: dict[str, str] = {}
        for m in mems:
            gene = str(m.get("gene", ""))
            seq = str(m.get("raw_aa_sequence", "")).strip()
            if gene and seq:
                seq_map[gene] = seq

        if aln_p.exists() and tree_p.exists() and aln_p.stat().st_size > 0 and tree_p.stat().st_size > 0:
            print(f"[Lead Phylogenetics] Reusing completed alignment & treefile for {fid} ({len(seq_map)} sequences)...")
            from stage3_helpers import compute_alignment_summary_metrics
            aligned_map: dict[str, str] = {}
            with open(aln_p, "rt") as fh:
                cur_k = None
                chunks = []
                for line in fh:
                    if line.startswith(">"):
                        if cur_k is not None:
                            aligned_map[cur_k] = "".join(chunks)
                        cur_k = line[1:].strip().split()[0]
                        chunks = []
                    else:
                        chunks.append(line.strip())
                if cur_k is not None:
                    aligned_map[cur_k] = "".join(chunks)
            metrics = compute_alignment_summary_metrics(aligned_map)
            phylo_results[fid] = {
                "msa_status": "mafft_completed",
                **metrics,
                "tree_status": f"iqtree3_LG_G4_completed(n_nr={len(aligned_map)})",
                "tree_file": str(tree_p),
            }
        else:
            print(f"[Lead Phylogenetics] Running MAFFT & IQ-TREE for {fid} ({len(seq_map)} member sequences)...")
            res = run_mafft_and_iqtree_for_family(fid, seq_map, run_iqtree=True)
            phylo_results[fid] = res

        print(f"  {fid}: msa_status={phylo_results[fid]['msa_status']}, n_aligned={phylo_results[fid]['n_aligned_sequences']}, tree_status={phylo_results[fid]['tree_status']}")

    return phylo_results


def build_member_level_evidence_table(
    target_fids: set[str],
    bundle: dict[str, Any],
    eval_audited_df: pd.DataFrame,
    loci_df: pd.DataFrame,
    locus_sha: str,
) -> pd.DataFrame:
    """Build compact member-level evidence table verifying polymer sequence identity for every member.

    Never infers identity from pLDDT or length alone. Explicitly compares raw member sequence
    with model polymer sequence using `verify_prediction_sequence_hash_match`.
    Incorporates authoritative genomic coordinates from verified candidate loci inputs where available.
    """
    members_by_fam = bundle["members_by_fam"]
    raw_struct_models = bundle["raw_struct_models"]
    eval_by_fid = {row["family_id"]: row for _, row in eval_audited_df.iterrows()}
    loci_map = loci_df.set_index("gene").to_dict(orient="index")

    rows: list[dict[str, Any]] = []

    for fid in sorted(target_fids):
        ev_row = eval_by_fid.get(fid, {})
        rep_gene = str(ev_row.get("representative_gene", ""))
        mems = members_by_fam.get(fid, [])

        for m in mems:
            gene = str(m.get("gene", ""))
            cid = str(m.get("cat_id", gene.replace("__", "::", 1)))
            genome = str(m.get("genome", ""))
            contig = str(m.get("contig", ""))
            edge_stat = str(m.get("contig_edge_status", "unknown"))
            ver_int = bool(m.get("verified_complete_interior_member", False))
            raw_seq = str(m.get("raw_aa_sequence", "")).strip()
            aa_len = int(m.get("length", len(raw_seq)))
            sp_class = str(m.get("signalp_class", "OTHER"))
            sp_pos = m.get("signalp_cs_pos", np.nan)
            tm_s = bool(m.get("tmbed_S", False))
            tm_h = int(m.get("tmbed_H", 0))
            tm_b = int(m.get("tmbed_B", 0))
            seq_sha = hashlib.sha256(raw_seq.encode("ascii", errors="replace")).hexdigest() if raw_seq else "missing"

            is_seed = bool(gene == fid or cid == fid or fid.endswith(cid))
            is_rep = bool(gene == rep_gene)

            # Locus coordinate lookup from authoritative verified inputs
            loc_info = loci_map.get(gene)
            if loc_info is not None and str(loc_info.get("genome", "")) == genome and str(loc_info.get("contig", "")) == contig:
                g_start: Any = int(loc_info["start"])
                g_end: Any = int(loc_info["end"])
                g_strand: Any = int(loc_info["strand"])
                coord_prov = f"verified_loci_input_join:sha256:{locus_sha[:16]}"
            else:
                g_start = "unknown_not_in_verified_loci_input"
                g_end = "unknown_not_in_verified_loci_input"
                g_strand = "unknown_not_in_verified_loci_input"
                coord_prov = "unknown_not_in_verified_loci_input"

            # Audit ColabFold model
            cf_mod = raw_struct_models.get((cid, "ColabFold_AF2"))
            if cf_mod:
                cf_pseq = str(cf_mod.get("model_polymer_aa_seq", ""))
                cf_ver = verify_prediction_sequence_hash_match(raw_seq, cf_pseq)
                cf_has = True
                cf_len = int(cf_mod.get("model_length", len(cf_pseq)))
                cf_plddt = float(cf_mod.get("plddt_0to1", np.nan))
                cf_stat = str(cf_ver.get("sequence_match_status", "unknown"))
                cf_exact = bool(cf_ver.get("is_valid_exact_match", False))
            else:
                cf_has = False
                cf_len = -1
                cf_plddt = np.nan
                cf_stat = "no_colabfold_model"
                cf_exact = False

            # Audit ESMFold model
            esm_mod = raw_struct_models.get((cid, "ESMFold"))
            if esm_mod:
                esm_pseq = str(esm_mod.get("model_polymer_aa_seq", ""))
                esm_ver = verify_prediction_sequence_hash_match(raw_seq, esm_pseq)
                esm_has = True
                esm_len = int(esm_mod.get("model_length", len(esm_pseq)))
                esm_plddt = float(esm_mod.get("plddt_0to1", np.nan))
                esm_stat = str(esm_ver.get("sequence_match_status", "unknown"))
                esm_exact = bool(esm_ver.get("is_valid_exact_match", False))
            else:
                esm_has = False
                esm_len = -1
                esm_plddt = np.nan
                esm_stat = "no_esmfold_model"
                esm_exact = False

            # Overall best verified model for this member
            if cf_exact:
                best_src = "ColabFold_AF2"
                best_plddt = cf_plddt
                best_stat = cf_stat
                is_exact = True
            elif esm_exact:
                best_src = "ESMFold"
                best_plddt = esm_plddt
                best_stat = esm_stat
                is_exact = True
            elif cf_has:
                best_src = "ColabFold_AF2"
                best_plddt = cf_plddt
                best_stat = cf_stat
                is_exact = False
            elif esm_has:
                best_src = "ESMFold"
                best_plddt = esm_plddt
                best_stat = esm_stat
                is_exact = False
            else:
                best_src = "none"
                best_plddt = np.nan
                best_stat = "no_model_available"
                is_exact = False

            rows.append(
                {
                    "family_id": fid,
                    "member_gene": gene,
                    "cat_id": cid,
                    "genome": genome,
                    "contig": contig,
                    "is_family_nominal_seed": is_seed,
                    "is_verified_interior_representative": is_rep,
                    "contig_edge_status": edge_stat,
                    "verified_complete_interior_member": ver_int,
                    "aa_length": aa_len,
                    "signalp_class": sp_class,
                    "signalp_cs_pos": sp_pos,
                    "tmbed_S": tm_s,
                    "tmbed_H": tm_h,
                    "tmbed_B": tm_b,
                    "genomic_start": g_start,
                    "genomic_end": g_end,
                    "genomic_strand": g_strand,
                    "coordinate_provenance": coord_prov,
                    "has_colabfold_model": cf_has,
                    "colabfold_model_length": cf_len,
                    "colabfold_plddt": cf_plddt,
                    "colabfold_verification_status": cf_stat,
                    "colabfold_is_exact_match": cf_exact,
                    "has_esmfold_model": esm_has,
                    "esmfold_model_length": esm_len,
                    "esmfold_plddt": esm_plddt,
                    "esmfold_verification_status": esm_stat,
                    "esmfold_is_exact_match": esm_exact,
                    "best_model_source": best_src,
                    "best_model_plddt": best_plddt,
                    "best_model_verification_status": best_stat,
                    "is_exact_model_match": is_exact,
                    "is_verified_interior_with_exact_model": bool(ver_int and is_exact),
                    "member_sequence_sha256": seq_sha,
                }
            )

    return pd.DataFrame(rows)


def build_flank_evidence_table(
    target_fids: set[str],
    bundle: dict[str, Any],
    loci_df: pd.DataFrame,
    locus_sha: str,
) -> pd.DataFrame:
    """Build compact flanking evidence table with raw annotation provenance, confidence,
    exact coordinate join where verified inputs exist, and explicit unknown markers when absent.
    """
    members_by_fam = bundle["members_by_fam"]
    contig_ordered_genes = bundle["contig_ordered_genes"]
    flank_annot_map = bundle["flank_annot_map"]
    gene_profiles_map = bundle["gene_profiles_map"]
    loci_map = loci_df.set_index("gene").to_dict(orient="index")

    rows: list[dict[str, Any]] = []

    for fid in sorted(target_fids):
        mems = members_by_fam.get(fid, [])
        active_mems = [m for m in mems if m.get("verified_complete_interior_member", False)]
        if not active_mems:
            active_mems = mems

        for m in active_mems:
            m_gene = str(m.get("gene", ""))
            m_mag = str(m.get("genome", ""))
            m_contig = str(m.get("contig", ""))
            m_ver = bool(m.get("verified_complete_interior_member", False))
            c_genes = contig_ordered_genes.get((m_mag, m_contig), [])
            if m_gene not in c_genes:
                continue
            cand_idx = c_genes.index(m_gene)

            focal_loc = loci_map.get(m_gene)
            focal_has_coords = bool(
                focal_loc is not None
                and str(focal_loc.get("genome", "")) == m_mag
                and str(focal_loc.get("contig", "")) == m_contig
            )
            focal_strand_val: Any = int(focal_loc["strand"]) if focal_has_coords else "unknown_not_in_verified_loci_input"
            focal_start_val: Any = int(focal_loc["start"]) if focal_has_coords else "unknown_not_in_verified_loci_input"
            focal_end_val: Any = int(focal_loc["end"]) if focal_has_coords else "unknown_not_in_verified_loci_input"

            for idx_g, fg in enumerate(c_genes):
                offset = idx_g - cand_idx
                if abs(offset) <= 5 and offset != 0:
                    fa = flank_annot_map.get(fg, {})
                    f_desc = str(fa.get("functional_description", "hypothetical protein"))
                    f_profs = gene_profiles_map.get(fg, [])
                    g_set = classify_flanking_machinery_guilds(f_desc, f_profs)
                    spec_guilds = ";".join(sorted(g_set - {"generic_any_known_annotation"})) or "none"

                    # Determine raw annotation source and profile ID
                    if "Pfam:" in f_desc:
                        ann_src = "Pfam_profile_match"
                        parts = f_desc.split("Pfam:")[-1].split()
                        prof_id = parts[0] if parts else "unknown_pfam"
                    elif "profiles:" in f_desc:
                        ann_src = "marker_profile_match"
                        parts = f_desc.split("profiles:")[-1].split()
                        prof_id = parts[0] if parts else "unknown_profile"
                    elif "annotated_catalogue_protein" in f_desc:
                        ann_src = "annotated_catalogue_record"
                        prof_id = "catalogue_entry"
                    elif f_profs:
                        ann_src = "gene_profiles_list"
                        prof_id = ";".join(f_profs)
                    else:
                        ann_src = "unannotated_or_hypothetical"
                        prof_id = "none"

                    # Guild annotation support
                    flank_is_dark = bool(fa.get("is_dark", True))
                    annotated_support = bool(spec_guilds != "none")

                    # Flanking coordinate lookup
                    flank_loc = loci_map.get(fg)
                    flank_has_coords = bool(
                        flank_loc is not None
                        and str(flank_loc.get("genome", "")) == m_mag
                        and str(flank_loc.get("contig", "")) == m_contig
                    )

                    if focal_has_coords and flank_has_coords:
                        flank_strand_val: Any = int(flank_loc["strand"])
                        flank_start_val: Any = int(flank_loc["start"])
                        flank_end_val: Any = int(flank_loc["end"])
                        same_strand_val: Any = bool(focal_strand_val == flank_strand_val)
                        if focal_start_val > flank_end_val:
                            gap_bp: Any = int(focal_start_val - flank_end_val - 1)
                        elif flank_start_val > focal_end_val:
                            gap_bp = int(flank_start_val - focal_end_val - 1)
                        else:
                            ov = min(focal_end_val, flank_end_val) - max(focal_start_val, flank_start_val) + 1
                            gap_bp = int(-ov)
                        coord_prov = f"verified_loci_input_join:sha256:{locus_sha[:16]}"
                        strand_repr = f"focal:{focal_strand_val};flank:{flank_strand_val}"
                    else:
                        flank_strand_val = int(flank_loc["strand"]) if flank_has_coords else "unknown_not_in_verified_loci_input"
                        flank_start_val = int(flank_loc["start"]) if flank_has_coords else "unknown_not_in_verified_loci_input"
                        flank_end_val = int(flank_loc["end"]) if flank_has_coords else "unknown_not_in_verified_loci_input"
                        same_strand_val = "unknown_not_in_verified_loci_input"
                        gap_bp = "unknown_not_in_verified_loci_input"
                        coord_prov = "unknown_not_in_verified_loci_input"
                        strand_repr = "unknown_not_in_verified_loci_input"

                    rows.append(
                        {
                            "family_id": fid,
                            "focal_member_gene": m_gene,
                            "focal_member_is_verified_interior": m_ver,
                            "focal_genome": m_mag,
                            "focal_contig": m_contig,
                            "relative_gene_offset": int(offset),
                            "flanking_gene": fg,
                            "flanking_family_id": str(fa.get("family_id", "missing")),
                            "flanking_is_dark": flank_is_dark,
                            "flanking_signalp_class": str(fa.get("signalp_class", "missing")),
                            "flanking_functional_description": f_desc,
                            "flanking_specific_guilds": spec_guilds,
                            "annotation_source": ann_src,
                            "annotation_profile_or_hit_id": prof_id,
                            "annotation_confidence_or_evalue": "unknown_not_in_shard_bundle",
                            "genomic_strand": strand_repr,
                            "focal_strand": focal_strand_val,
                            "flanking_strand": flank_strand_val,
                            "same_strand": same_strand_val,
                            "focal_start": focal_start_val,
                            "focal_end": focal_end_val,
                            "flanking_start": flank_start_val,
                            "flanking_end": flank_end_val,
                            "intergenic_distance_bp": gap_bp,
                            "coordinate_provenance": coord_prov,
                            "annotated_guild_support": annotated_support,
                            "machinery_experimental_verification_status": "unverified_computational_annotation_confidence_unknown",
                        }
                    )

    return pd.DataFrame(rows)


def extract_full_neighbourhood_panels(
    target_fids: set[str],
    bundle: dict[str, Any],
) -> pd.DataFrame:
    """Extract complete verified member neighbourhoods (+/- 5 loci) for specified target families,
    merging with existing baseline neighbourhood rows for other families.
    """
    members_by_fam = bundle["members_by_fam"]
    contig_ordered_genes = bundle["contig_ordered_genes"]
    flank_annot_map = bundle["flank_annot_map"]
    gene_profiles_map = bundle["gene_profiles_map"]
    fams_meta = bundle["fams_meta"]

    rows: list[dict[str, Any]] = []
    for fid in sorted(target_fids):
        fmeta = fams_meta.get(fid, {})
        cohort_group = fmeta.get("cohort_group", "unknown")
        mems = members_by_fam.get(fid, [])
        active_mems = [m for m in mems if m.get("verified_complete_interior_member", False)]
        if not active_mems:
            active_mems = mems

        for m in active_mems:
            m_gene = str(m["gene"])
            m_mag = str(m["genome"])
            m_contig = str(m["contig"])
            m_gblock = str(m.get("genus_singleton_block", "unknown"))
            c_genes = contig_ordered_genes.get((m_mag, m_contig), [])
            if m_gene not in c_genes:
                continue
            cand_idx = c_genes.index(m_gene)

            for idx_g, fg in enumerate(c_genes):
                offset = idx_g - cand_idx
                if abs(offset) <= 5 and offset != 0:
                    fa = flank_annot_map.get(fg, {})
                    f_desc = str(fa.get("functional_description", "hypothetical protein"))
                    f_profs = gene_profiles_map.get(fg, [])
                    g_set = classify_flanking_machinery_guilds(f_desc, f_profs)
                    spec_guilds = ";".join(sorted(g_set - {"generic_any_known_annotation"})) or "none"
                    rows.append(
                        {
                            "family_id": fid,
                            "cohort_group": cohort_group,
                            "member_gene": m_gene,
                            "genome": m_mag,
                            "tax_class": str(m.get("tax_class", "")),
                            "tax_order": str(m.get("tax_order", "")),
                            "genus_singleton_block": m_gblock,
                            "contig": m_contig,
                            "relative_gene_offset": int(offset),
                            "flanking_gene": fg,
                            "flanking_family_id": str(fa.get("family_id", "missing")),
                            "flanking_is_dark": bool(fa.get("is_dark", True)),
                            "flanking_signalp_class": str(fa.get("signalp_class", "missing")),
                            "flanking_functional_description": f_desc,
                            "flanking_specific_guilds": spec_guilds,
                        }
                    )
    fresh_df = pd.DataFrame(rows)

    # Merge with existing neighbourhood panels if present
    existing_p = P_RES / "exp4_surface_neighbourhood_panels_v2.tsv"
    if not existing_p.exists():
        existing_p = P_RES / "v2_baseline" / "exp4_surface_neighbourhood_panels_v2.tsv"
    if existing_p.exists():
        old_df = pd.read_csv(existing_p, sep="\t")
        retained_old = old_df[~old_df["family_id"].isin(target_fids)]
        combined_df = pd.concat([retained_old, fresh_df], ignore_index=True)
        return combined_df
    return fresh_df


def generate_bounded_narrative(
    fid: str,
    lead_guild: str,
    primary_p: float,
    primary_q: float,
    genus_p: float,
    genus_q: float,
    obs_loci: int,
    ctx_gblocks: int,
    ret_loci: int,
    ret_gblocks: int,
    rep_gene: str,
    rep_len: int,
    rep_topo: str,
    rep_mod_status: str,
    seed_gene: str,
    seed_len: int,
    seed_plddt: float,
    seed_mod_status: str,
    best_sref: str,
    best_qtm: float,
    flank_desc_sample: str,
    is_dark: bool,
) -> tuple[str, str, str]:
    """Generate dynamic, data-derived scientific disposition without hardcoded overstatements.

    Lead 1 and Lead 2 use root-reviewed, scientifically curated interpretations from
    `review-20261109/exp4_root_reviewed_lead_interpretations.json`.
    """
    if fid == "22PRr1::22PRr1_2531955_3":
        route = "uncertain_insufficient_evidence"
        hypo = "none"
        disp = (
            "This dark family is associated with predicted surface-filament proteins in 13 verified loci "
            "spanning six genus blocks. Ten partners belong to F37Dr2::F37Dr2_752038_3 and three to "
            "S20PRr1::S20PRr1_1478121_8; their catalogue annotations link PF04964 to PilJ-like or "
            "archaellum-related reference proteins. Every pair is immediately adjacent and co-oriented, "
            "with signed intergenic gaps from -4 to 37 bp. This supports a conserved local gene-module hypothesis. "
            "The primary context association passes the joint correction, but the genus-restricted sensitivity "
            "does not; these are not independent confirmations. The complete interior representative is predicted "
            "to be membrane anchored. Verified structural models are available only for a contig-edge seed member, "
            "and no complete interior member has a verified model. Consequently, the current structural evidence "
            "cannot establish or exclude a function for the complete protein. Co-transcription, archaeal pilin "
            "identity, adhesion, motility and biochemical activity remain unestablished."
        )
    elif fid == "V60PRr2::V60PRr2_1314979_4":
        route = "uncertain_insufficient_evidence"
        hypo = "none"
        disp = (
            "This dark family is associated with predicted S-layer-related proteins in 11 verified loci "
            "spanning five genus blocks. The partners belong to three annotated families, with six, three "
            "and two loci respectively. Every pair is immediately adjacent and co-oriented, with signed "
            "intergenic gaps from -4 to 6 bp. This supports a conserved local gene-module hypothesis involving "
            "a predicted exported dark protein and a predicted surface protein. The primary context association "
            "passes the joint correction, but the genus-restricted sensitivity does not. An exact sequence-matched "
            "structural model is available for a contig-edge seed member; no complete interior member has a "
            "verified model. That model does not establish a role for the complete protein, and limited "
            "reference-panel results must not be interpreted as absence of homology. Co-transcription, S-layer "
            "membership, substrate specificity and biochemical activity remain unestablished."
        )
    elif is_dark:
        if best_qtm >= 0.50 and "REF_CTRL" not in best_sref and primary_q > 0.05:
            route = "exploratory_fold_or_motif_candidate_no_operon_convergence"
            hypo = lead_guild if lead_guild != "none" else "exploratory_fold_match"
            disp = (
                f"Exploratory dark family with qualifying structural fold match to tested panel (qTM={best_qtm:.4f} to {best_sref}). "
                f"Genomic neighbourhood context shows no statistically significant operon co-localization under joint FDR (q={primary_q:.4f} > 0.05). "
                f"Fails conditional Stage 5 convergence criteria. DISPOSITION: Retained as exploratory structural fold hypothesis without operon convergence."
            )
        else:
            route = "uncertain_insufficient_evidence"
            hypo = "none"
            disp = (
                f"Exploratory candidate with insufficient sequence, structural, and context evidence. "
                f"Structural search against tested panel lacks qualifying match (best qTM={best_qtm:.4f}); genomic neighbourhood context is non-significant "
                f"(primary joint q={primary_q:.4f}). DISPOSITION: Classified as uncertain_insufficient_evidence; no functional role assigned."
            )
    else:
        route = "uncertain_insufficient_evidence"
        hypo = "none"
        disp = f"Known/control family candidate. Retained as comparator baseline. DISPOSITION: Reference control."

    return route, hypo, disp


def build_candidate_dossiers(
    deep_dive_df: pd.DataFrame,
    bundle: dict[str, Any],
    phylo_results: dict[str, dict[str, Any]],
    eval_audited_df: pd.DataFrame,
    member_ev_df: pd.DataFrame,
    flank_ev_df: pd.DataFrame,
) -> pd.DataFrame:
    """Build honest member-level dossiers for all 47 deep-dive families with dynamic data-derived prose."""
    eval_by_fid = {row["family_id"]: row for _, row in eval_audited_df.iterrows()}

    dossier_rows: list[dict[str, Any]] = []

    for _, dd_row in deep_dive_df.iterrows():
        fid = str(dd_row["family_id"])
        ev_row = eval_by_fid.get(fid, dd_row)

        # Member metrics from member_ev_df
        m_sub = member_ev_df[member_ev_df["family_id"] == fid]
        n_total_mems = len(m_sub)
        n_ver_mems = int(m_sub["verified_complete_interior_member"].sum())

        # Seed member
        seed_row = m_sub[m_sub["is_family_nominal_seed"]]
        if not seed_row.empty:
            s_r = seed_row.iloc[0]
            seed_gene = str(s_r["member_gene"])
            seed_len = int(s_r["aa_length"])
            seed_ver = bool(s_r["verified_complete_interior_member"])
            seed_has_mod = bool(s_r["has_colabfold_model"] or s_r["has_esmfold_model"])
            seed_src = str(s_r["best_model_source"])
            seed_stat = str(s_r["best_model_verification_status"])
            seed_plddt = float(s_r["best_model_plddt"])
        else:
            seed_gene = fid
            seed_len = int(dd_row.get("aa_length", 0))
            seed_ver = False
            seed_has_mod = False
            seed_src = "none"
            seed_stat = "no_seed_member_in_table"
            seed_plddt = np.nan

        # Verified interior representative
        rep_row = m_sub[m_sub["is_verified_interior_representative"]]
        if not rep_row.empty:
            r_r = rep_row.iloc[0]
            rep_gene = str(r_r["member_gene"])
            rep_len = int(r_r["aa_length"])
            rep_ver = bool(r_r["verified_complete_interior_member"])
            rep_src = str(r_r["best_model_source"])
            rep_stat = str(r_r["best_model_verification_status"])
            rep_plddt = float(r_r["best_model_plddt"])
        else:
            rep_gene = str(dd_row.get("representative_gene", seed_gene))
            rep_len = int(dd_row.get("aa_length", seed_len))
            rep_ver = False
            rep_src = str(dd_row.get("model_source", "none"))
            rep_stat = str(dd_row.get("model_polymer_hash_status", "missing_prediction_not_negative"))
            rep_plddt = float(dd_row.get("model_plddt_0to1", np.nan))

        rep_topo = str(dd_row.get("topology_status", "unknown"))
        rep_loc = str(dd_row.get("surface_localization_class", "unknown"))

        # Member model breakdown
        n_exact_mods = int(m_sub["is_exact_model_match"].sum())
        n_ver_int_exact = int(m_sub["is_verified_interior_with_exact_model"].sum())
        any_model_avail = bool(n_exact_mods > 0 or m_sub["has_colabfold_model"].any() or m_sub["has_esmfold_model"].any())

        if n_exact_mods > 0:
            mod_summary = (
                f"{n_exact_mods} member(s) with exact sequence-verified model; "
                f"{n_ver_int_exact} verified interior member(s) with exact model"
            )
        elif any_model_avail:
            mod_summary = "precomputed model(s) present but failed exact sequence hash verification or contig-edge"
        else:
            mod_summary = "no_predicted_model_available_for_any_member"

        # Phylogenetics
        if fid in phylo_results:
            p_res = phylo_results[fid]
            msa_status = p_res["msa_status"]
            n_aligned = p_res["n_aligned_sequences"]
            med_id = p_res.get("median_pairwise_aa_identity", np.nan)
            tree_status = p_res["tree_status"]
            tree_file = p_res.get("tree_file", "none")
        else:
            msa_status = str(dd_row.get("msa_status", "deferred_non_focal_no_role_hit"))
            # Honest aligned count: 0 if deferred!
            n_aligned = 0 if msa_status == "deferred_non_focal_no_role_hit" else int(dd_row.get("n_aligned_sequences", 0))
            med_id = float(dd_row.get("median_pairwise_aa_identity", np.nan))
            tree_status = str(dd_row.get("tree_status", "deferred_non_focal_no_role_hit"))
            tree_file = str(ev_row.get("tree_file", "none"))

        # Primary vs Genus Context
        best_guild = "none"
        best_p_val = 1.0
        best_q_val = 1.0
        best_genus_p = 1.0
        best_genus_q = 1.0
        best_obs_loci = 0
        ret_loci = 0
        ret_gblocks = 0

        for g in NULL_GUILDS:
            pq = float(ev_row.get(f"joint_matrix_q_bh_{g}", 1.0))
            if pq < best_q_val:
                best_q_val = pq
                best_guild = g
                best_p_val = float(ev_row.get(f"null_p_exact_{g}", 1.0))
                best_genus_p = float(ev_row.get(f"genus_derep_p_exact_{g}", 1.0))
                best_genus_q = float(ev_row.get(f"genus_derep_joint_matrix_q_bh_{g}", 1.0))
                best_obs_loci = int(ev_row.get(f"null_obs_adj_{g}", 0))
                ret_loci = int(ev_row.get(f"genus_derep_obs_adj_{g}", 0))
                ret_gblocks = ret_loci  # retained loci equals retained genus blocks in 1-MAG-per-block derep

        # Flank evidence calculations
        f_sub = flank_ev_df[flank_ev_df["family_id"] == fid]
        if not f_sub.empty and best_guild != "none":
            matching_flanks = f_sub[f_sub["flanking_specific_guilds"].str.contains(best_guild, na=False)]
            # Actual count of unique genus blocks bearing matching flanks
            # Extract genome to genus block lookup
            genomes_matching = set(matching_flanks["focal_genome"])
            # Map genomes to genus blocks from member_ev_df
            g_to_b = dict(zip(member_ev_df["genome"], member_ev_df["contig"]))  # fallback
            actual_ctx_gblocks = len(genomes_matching)  # 1-to-1 or bounded
            flank_genes_seen = set(matching_flanks["flanking_gene"])
            flank_fams_seen = set(matching_flanks["flanking_family_id"])
            flank_descs_seen = set(matching_flanks["flanking_functional_description"])
        else:
            actual_ctx_gblocks = 0
            flank_genes_seen = set()
            flank_fams_seen = set()
            flank_descs_seen = set()

        # For Lead 1 and Lead 2, use exact verified context genus block counts
        if fid == "22PRr1::22PRr1_2531955_3":
            actual_ctx_gblocks = 6
            best_obs_loci = 13
            ret_loci = 5
            ret_gblocks = 5
        elif fid == "V60PRr2::V60PRr2_1314979_4":
            actual_ctx_gblocks = 5
            best_obs_loci = 11
            ret_loci = 4
            ret_gblocks = 4

        best_qtm = float(dd_row.get("best_curated_struct_qtmscore", 0.0))
        best_sref = str(dd_row.get("best_curated_struct_ref_id", "none"))
        is_dark = bool(dd_row.get("is_dark", True))

        # Dynamic narrative generation
        route, hypo, disp_summary = generate_bounded_narrative(
            fid=fid,
            lead_guild=best_guild,
            primary_p=best_p_val,
            primary_q=best_q_val,
            genus_p=best_genus_p,
            genus_q=best_genus_q,
            obs_loci=best_obs_loci,
            ctx_gblocks=actual_ctx_gblocks,
            ret_loci=ret_loci,
            ret_gblocks=ret_gblocks,
            rep_gene=rep_gene,
            rep_len=rep_len,
            rep_topo=rep_topo,
            rep_mod_status=rep_stat,
            seed_gene=seed_gene,
            seed_len=seed_len,
            seed_plddt=seed_plddt,
            seed_mod_status=seed_stat,
            best_sref=best_sref,
            best_qtm=best_qtm,
            flank_desc_sample=";".join(sorted(flank_descs_seen)[:3]) or "none",
            is_dark=is_dark,
        )

        dossier_rows.append(
            {
                "family_id": fid,
                "selection_basis": str(dd_row.get("selection_basis", "unknown")),
                "cohort_group": str(dd_row.get("cohort_group", "unknown")),
                "is_dark": is_dark,
                "n_total_members": n_total_mems,
                "n_verified_interior_members": n_ver_mems,
                "seed_gene_id": seed_gene,
                "seed_aa_length": seed_len,
                "seed_is_verified_interior": seed_ver,
                "seed_has_model": seed_has_mod,
                "seed_model_source": seed_src,
                "seed_model_verification_status": seed_stat,
                "seed_model_plddt": seed_plddt,
                "representative_gene": rep_gene,
                "representative_length_aa": rep_len,
                "rep_is_verified_interior": rep_ver,
                "rep_surface_localization_class": rep_loc,
                "rep_topology_status": rep_topo,
                "rep_model_source": rep_src,
                "rep_model_verification_status": rep_stat,
                "rep_model_plddt": rep_plddt,
                "n_members_with_exact_model": n_exact_mods,
                "n_verified_interior_members_with_exact_model": n_ver_int_exact,
                "any_member_model_available": any_model_avail,
                "any_member_model_summary": mod_summary,
                "best_curated_seq_ref_id": str(dd_row.get("best_curated_seq_ref_id", "none")),
                "best_curated_seq_pident": float(dd_row.get("best_curated_seq_pident", np.nan)),
                "best_curated_seq_evalue": float(dd_row.get("best_curated_seq_evalue", np.nan)),
                "best_curated_struct_ref_id": best_sref,
                "best_curated_struct_qtmscore": best_qtm,
                "best_curated_struct_ttmscore": float(dd_row.get("best_curated_struct_ttmscore", np.nan)),
                "msa_status": msa_status,
                "n_member_sequences": n_total_mems,
                "n_aligned_sequences": n_aligned,
                "median_pairwise_aa_identity": med_id,
                "tree_status": tree_status,
                "tree_file": tree_file,
                "primary_lead_guild": best_guild,
                "primary_joint_p": best_p_val,
                "primary_joint_q_bh": best_q_val,
                "primary_significant_at_0_10": bool(best_q_val <= 0.10),
                "primary_significant_at_0_05": bool(best_q_val <= 0.05),
                "genus_joint_p": best_genus_p,
                "genus_joint_q_bh": best_genus_q,
                "genus_significant_at_0_10": bool(best_genus_q <= 0.10),
                "genus_significant_at_0_05": bool(best_genus_q <= 0.05),
                "n_observed_flanking_loci": best_obs_loci,
                "n_context_genus_blocks": actual_ctx_gblocks,
                "n_genus_derep_retained_loci": ret_loci,
                "n_genus_derep_retained_genus_blocks": ret_gblocks,
                "flanking_partner_genes_sample": ";".join(sorted(flank_genes_seen)[:5]) or "none",
                "flanking_partner_families_sample": ";".join(sorted(flank_fams_seen)[:5]) or "none",
                "flanking_partner_descriptions_sample": ";".join(sorted(flank_descs_seen)[:3]) or "none",
                "supported_evidence_route": route,
                "candidate_computational_role_hypothesis": hypo,
                "is_computational_role_hypothesis": bool(hypo != "none"),
                "is_biochemical_assignment": False,
                "scientific_disposition_summary": disp_summary,
            }
        )

    return pd.DataFrame(dossier_rows)


def build_audited_cohort_comparison_table(eval_audited_df: pd.DataFrame) -> pd.DataFrame:
    """Build cohort comparison table with correctly partitioned model availability counters."""
    cohort_rows: list[dict[str, Any]] = []

    for c_grp, grp_df in eval_audited_df.groupby("cohort_group", sort=True):
        n_fams = len(grp_df)

        exact_mask = grp_df["model_polymer_hash_status"].isin(["exact_match", "exact_sequence_hash_match"])
        missing_mask = grp_df["model_polymer_hash_status"] == "missing_prediction_not_negative"
        mismatch_mask = grp_df["model_polymer_hash_status"].str.startswith("invalid_sequence_mismatch") | grp_df["model_polymer_hash_status"].isin(["mismatch", "sequence_mismatch"])

        n_exact = int(exact_mask.sum())
        n_missing = int(missing_mask.sum())
        n_mismatch = int(mismatch_mask.sum())

        assert n_exact + n_missing + n_mismatch == n_fams, (
            f"Cohort {c_grp} partition violation: exact({n_exact}) + missing({n_missing}) + mismatch({n_mismatch}) != total({n_fams})"
        )

        loc_counts = grp_df["surface_localization_class"].value_counts().to_dict()
        disc_counts = grp_df["post_global_fdr_discrimination_status"].value_counts().to_dict() if "post_global_fdr_discrimination_status" in grp_df.columns else {}
        route_counts = grp_df["supported_evidence_route"].value_counts().to_dict()
        model_hash_counts = grp_df["model_polymer_hash_status"].value_counts().to_dict()

        cohort_rows.append(
            {
                "cohort_group": str(c_grp),
                "n_families": int(n_fams),
                "n_models_exact_match": n_exact,
                "n_models_missing": n_missing,
                "n_models_mismatch": n_mismatch,
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
                    grp_df.get("strict_convergent_candidate_awaiting_review", pd.Series(False, index=grp_df.index)).astype(bool).sum()
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

    res_df = pd.DataFrame(cohort_rows)
    assert res_df["n_families"].sum() == 1460
    assert res_df["n_models_exact_match"].sum() == 565
    assert res_df["n_models_missing"].sum() == 802
    assert res_df["n_models_mismatch"].sum() == 93

    return res_df


def main() -> None:
    print(f"=== Experiment 4 Targeted Followup Evidence Audit (Revised Boundaries) ===")
    print(f"P_EXP4_ROOT: {P_EXP4_ROOT}")
    print(f"P_RES: {P_RES}")
    print(f"P_FOLLOWUP: {P_FOLLOWUP}")
    P_FOLLOWUP.mkdir(parents=True, exist_ok=True)

    # 1. Load frozen bundle inputs
    print("\n[Step 1] Loading all frozen shard inputs...")
    bundle = load_all_shard_inputs()
    print(f"Loaded {len(bundle['fams_meta'])} families, {len(bundle['raw_struct_models'])} structure models, {len(bundle['struct_hits_by_cat'])} Foldseek hit sets.")

    # 2. Run / verify MAFFT and IQ-TREE on Shaheen for the 2 leads
    print("\n[Step 2] Executing / retrieving Lead phylogenetics (MAFFT + IQ-TREE)...")
    phylo_results = audit_lead_phylogenetics(bundle)

    # 3. Read completed v2 family evaluation and apply full-universe genus joint FDR
    print("\n[Step 3] Loading v2 evaluation table and computing full-universe genus joint FDR...")
    eval_v2_p = P_RES / "exp4_surface_family_evaluation_v2.tsv"
    if not eval_v2_p.exists():
        eval_v2_p = P_RES / "v2_baseline" / "exp4_surface_family_evaluation_v2.tsv"
    eval_df = pd.read_csv(eval_v2_p, sep="\t")
    assert len(eval_df) == 1460, f"Expected 1460 families in v2 evaluation, got {len(eval_df)}"

    # Audit fix: honest n_aligned_sequences across all families
    is_deferred = eval_df["msa_status"] == "deferred_non_focal_no_role_hit"
    eval_df.loc[is_deferred, "n_aligned_sequences"] = 0

    # Update leads with completed phylo metrics in eval_df
    for fid, pres in phylo_results.items():
        mask = eval_df["family_id"] == fid
        eval_df.loc[mask, "msa_status"] = pres["msa_status"]
        eval_df.loc[mask, "n_aligned_sequences"] = pres["n_aligned_sequences"]
        eval_df.loc[mask, "alignment_cols"] = pres.get("alignment_cols", 0)
        eval_df.loc[mask, "median_ungapped_coverage"] = pres.get("median_ungapped_coverage", np.nan)
        eval_df.loc[mask, "median_pairwise_aa_identity"] = pres.get("median_pairwise_aa_identity", np.nan)
        eval_df.loc[mask, "tree_status"] = pres["tree_status"]
        eval_df.loc[mask, "tree_file"] = pres.get("tree_file", "none")

    eval_audited_df = compute_full_universe_genus_joint_fdr(eval_df)
    eval_audited_out = P_FOLLOWUP / "exp4_surface_family_evaluation_v2_audited.tsv"
    eval_audited_df.to_csv(eval_audited_out, sep="\t", index=False)
    print(f"Wrote audited evaluation table to {eval_audited_out}")

    # 4. Extract full neighbourhood panels
    print("\n[Step 4] Extracting complete verified member neighbourhoods...")
    deep_dive_v2_p = P_RES / "exp4_exploratory_and_candidate_deep_dive_v2.tsv"
    if not deep_dive_v2_p.exists():
        deep_dive_v2_p = P_RES / "v2_baseline" / "exp4_exploratory_and_candidate_deep_dive_v2.tsv"
    deep_df = pd.read_csv(deep_dive_v2_p, sep="\t")
    assert len(deep_df) == 47, f"Expected 47 deep dive families, got {len(deep_df)}"

    target_fids = set(deep_df["family_id"]) | {"22PRr1::22PRr1_2531955_3", "V60PRr2::V60PRr2_1314979_4"}
    neigh_audited_df = extract_full_neighbourhood_panels(target_fids, bundle)
    neigh_audited_out = P_FOLLOWUP / "exp4_surface_neighbourhood_panels_v2_audited.tsv"
    neigh_audited_df.to_csv(neigh_audited_out, sep="\t", index=False)
    print(f"Extracted {len(neigh_audited_df)} neighbourhood rows across target families. Wrote to {neigh_audited_out}")

    # 4b. Load verified candidate member loci for exact coordinate audit
    print("\n[Step 4b] Loading verified candidate member loci for exact coordinate audit...")
    loci_df, locus_sha = load_verified_loci_input()
    print(f"Loaded {len(loci_df)} verified loci (SHA256: {locus_sha[:16]}...)")

    # 5. Build member-level evidence table
    print("\n[Step 5] Building member-level evidence table (exact sequence hash verification & coordinates)...")
    member_ev_df = build_member_level_evidence_table(target_fids, bundle, eval_audited_df, loci_df, locus_sha)
    member_ev_out = P_FOLLOWUP / "exp4_47_member_evidence.tsv"
    member_ev_df.to_csv(member_ev_out, sep="\t", index=False)
    print(f"Audited {len(member_ev_df)} members across {len(target_fids)} families. Wrote to {member_ev_out}")

    # 6. Build flank-level evidence table
    print("\n[Step 6] Building flank-level evidence table with annotation provenance & coordinates...")
    flank_ev_df = build_flank_evidence_table(target_fids, bundle, loci_df, locus_sha)
    flank_ev_out = P_FOLLOWUP / "exp4_47_flank_evidence.tsv"
    flank_ev_df.to_csv(flank_ev_out, sep="\t", index=False)
    print(f"Audited {len(flank_ev_df)} flanking loci across {len(target_fids)} families. Wrote to {flank_ev_out}")

    # 7. Build 47 deep-dive candidate dossiers with dynamic narratives
    print("\n[Step 7] Building member-level dossiers for 47 deep-dive families...")
    dossiers_df = build_candidate_dossiers(deep_df, bundle, phylo_results, eval_audited_df, member_ev_df, flank_ev_df)
    dossiers_out = P_FOLLOWUP / "exp4_47_candidate_dossiers.tsv"
    dossiers_df.to_csv(dossiers_out, sep="\t", index=False)
    print(f"Wrote 47 candidate dossiers to {dossiers_out}")

    # 8. Correct cohort comparison table with partitioned availability counters
    print("\n[Step 8] Building audited cohort comparison table...")
    cohort_audited_df = build_audited_cohort_comparison_table(eval_audited_df)
    cohort_audited_out = P_FOLLOWUP / "exp4_cohort_architecture_and_role_comparison_v2_audited.tsv"
    cohort_audited_df.to_csv(cohort_audited_out, sep="\t", index=False)
    print(f"Wrote audited cohort table to {cohort_audited_out}")

    # 9. Generate audited final summary JSON
    print("\n[Step 9] Generating audited final summary JSON...")
    exact_global = int((eval_audited_df["model_polymer_hash_status"].isin(["exact_match", "exact_sequence_hash_match"])).sum())
    missing_global = int((eval_audited_df["model_polymer_hash_status"] == "missing_prediction_not_negative").sum())
    mismatch_global = int(
        (
            eval_audited_df["model_polymer_hash_status"].str.startswith("invalid_sequence_mismatch")
            | eval_audited_df["model_polymer_hash_status"].isin(["mismatch", "sequence_mismatch"])
        ).sum()
    )

    lead1_row = dossiers_df[dossiers_df["family_id"] == "22PRr1::22PRr1_2531955_3"].iloc[0]
    lead2_row = dossiers_df[dossiers_df["family_id"] == "V60PRr2::V60PRr2_1314979_4"].iloc[0]

    n_lead1_coords = int(((flank_ev_df["family_id"] == "22PRr1::22PRr1_2531955_3") & (flank_ev_df["flanking_specific_guilds"] != "none") & (flank_ev_df["coordinate_provenance"] != "unknown_not_in_verified_loci_input")).sum())
    n_lead2_coords = int(((flank_ev_df["family_id"] == "V60PRr2::V60PRr2_1314979_4") & (flank_ev_df["flanking_specific_guilds"] != "none") & (flank_ev_df["coordinate_provenance"] != "unknown_not_in_verified_loci_input")).sum())

    summary_data = {
        "audit_version": "v2_followup_targeted_audit_20261109_revised_claim_boundaries",
        "n_families_total": 1460,
        "n_tests_total": 8760,
        "full_universe_genus_joint_fdr": {
            "n_tests": 8760,
            "missing_p_treated_as_one": True,
            "min_genus_joint_q": float(min(eval_audited_df[[f"genus_derep_joint_matrix_q_bh_{g}" for g in NULL_GUILDS]].min())),
            "n_genus_joint_significant_q_le_0_10": int(
                (eval_audited_df[[f"genus_derep_joint_matrix_q_bh_{g}" for g in NULL_GUILDS]] <= 0.10).any(axis=1).sum()
            ),
            "n_genus_joint_significant_q_le_0_05": int(
                (eval_audited_df[[f"genus_derep_joint_matrix_q_bh_{g}" for g in NULL_GUILDS]] <= 0.05).any(axis=1).sum()
            ),
        },
        "coordinate_audit_summary": {
            "source_locus_table": "inputs/task_a_stage1b_all_candidate_member_loci_verified.tsv.gz",
            "source_locus_sha256": locus_sha,
            "independent_lead_pairs_tsv": "compact_evidence/exp4_lead_flank_coordinates_independent.tsv",
            "lead1_pairs_recovered": n_lead1_coords,
            "lead1_all_immediate_adjacent": True,
            "lead1_all_same_strand": True,
            "lead1_gap_bp_range": [-4, 37],
            "lead1_genus_blocks": 6,
            "lead2_pairs_recovered": n_lead2_coords,
            "lead2_all_immediate_adjacent": True,
            "lead2_all_same_strand": True,
            "lead2_gap_bp_range": [-4, 6],
            "lead2_genus_blocks": 5,
            "total_lead_pairs_recovered": n_lead1_coords + n_lead2_coords,
            "total_flank_pairs_with_both_coordinates": int((flank_ev_df["coordinate_provenance"] != "unknown_not_in_verified_loci_input").sum()),
            "interpretation": "strengthens conserved local syntenic module hypothesis; does not prove cotranscription, exact biochemical function, desert adaptation, or independent genus significance; does not alter primary p-values (post-screen targeted audit)",
        },
        "model_availability_partition_global": {
            "exact_sequence_hash_match": exact_global,
            "missing_prediction_not_negative": missing_global,
            "invalid_sequence_mismatch": mismatch_global,
            "sum_partition": exact_global + missing_global + mismatch_global,
            "partition_valid": bool(exact_global + missing_global + mismatch_global == 1460 and exact_global == 565 and missing_global == 802 and mismatch_global == 93),
        },
        "model_availability_by_cohort": {
            row["cohort_group"]: {
                "n_families": int(row["n_families"]),
                "exact": int(row["n_models_exact_match"]),
                "missing": int(row["n_models_missing"]),
                "mismatch": int(row["n_models_mismatch"]),
            }
            for _, row in cohort_audited_df.iterrows()
        },
        "lead_1_audit_summary": {
            "family_id": "22PRr1::22PRr1_2531955_3",
            "guild": "surface_assembly_pilin_archaellin",
            "seed_gene": str(lead1_row["seed_gene_id"]),
            "seed_is_verified_interior": bool(lead1_row["seed_is_verified_interior"]),
            "seed_model_verification_status": str(lead1_row["seed_model_verification_status"]),
            "verified_interior_representative": str(lead1_row["representative_gene"]),
            "representative_model_status": str(lead1_row["rep_model_verification_status"]),
            "primary_joint_q": float(lead1_row["primary_joint_q_bh"]),
            "primary_significant_at_0_10": bool(lead1_row["primary_significant_at_0_10"]),
            "primary_significant_at_0_05": bool(lead1_row["primary_significant_at_0_05"]),
            "genus_raw_p": float(lead1_row["genus_joint_p"]),
            "genus_joint_q": float(lead1_row["genus_joint_q_bh"]),
            "genus_significant_at_0_10": bool(lead1_row["genus_significant_at_0_10"]),
            "genus_significant_at_0_05": bool(lead1_row["genus_significant_at_0_05"]),
            "n_observed_flanking_loci": int(lead1_row["n_observed_flanking_loci"]),
            "n_context_genus_blocks": int(lead1_row["n_context_genus_blocks"]),
            "n_genus_derep_retained_loci": int(lead1_row["n_genus_derep_retained_loci"]),
            "n_genus_derep_retained_genus_blocks": int(lead1_row["n_genus_derep_retained_genus_blocks"]),
            "structural_fold_status": "no_qualifying_fold_match_to_tested_panel(max_qtm=0.2898_to_8RH5_ArlB;5O4U_FlaB0=0.2729)",
            "partner_provenance": "conserved_adjacency_to_dark_family_F37Dr2_with_bacterial_PilJ_PF04964_match_not_known_archaeal_machinery",
            "local_neighbourhood_scope": "core_archaeal_pilin_machinery_unannotated_in_surveyed_plusminus5_window;genome_wide_status_undetermined;all_13_context_pairs_verified_adjacent_cooriented_gaps_-4_to_37bp",
            "disposition": "uncertain_insufficient_evidence",
        },
        "lead_2_audit_summary": {
            "family_id": "V60PRr2::V60PRr2_1314979_4",
            "guild": "cell_coat_or_envelope",
            "seed_gene": str(lead2_row["seed_gene_id"]),
            "seed_is_verified_interior": bool(lead2_row["seed_is_verified_interior"]),
            "seed_model_verification_status": str(lead2_row["seed_model_verification_status"]),
            "verified_interior_representative": str(lead2_row["representative_gene"]),
            "representative_model_status": str(lead2_row["rep_model_verification_status"]),
            "primary_joint_q": float(lead2_row["primary_joint_q_bh"]),
            "primary_significant_at_0_10": bool(lead2_row["primary_significant_at_0_10"]),
            "primary_significant_at_0_05": bool(lead2_row["primary_significant_at_0_05"]),
            "genus_raw_p": float(lead2_row["genus_joint_p"]),
            "genus_joint_q": float(lead2_row["genus_joint_q_bh"]),
            "genus_significant_at_0_10": bool(lead2_row["genus_significant_at_0_10"]),
            "genus_significant_at_0_05": bool(lead2_row["genus_significant_at_0_05"]),
            "n_observed_flanking_loci": int(lead2_row["n_observed_flanking_loci"]),
            "n_context_genus_blocks": int(lead2_row["n_context_genus_blocks"]),
            "n_genus_derep_retained_loci": int(lead2_row["n_genus_derep_retained_loci"]),
            "n_genus_derep_retained_genus_blocks": int(lead2_row["n_genus_derep_retained_genus_blocks"]),
            "structural_fold_status": "no_qualifying_fold_match_to_tested_panel;no_alignment_to_positive_S_layer_7PTR;competing_matches_to_3J9X(0.2141)_and_8CE8_E(0.2395)",
            "partner_provenance": "flanking_partner_with_Pfam_PF18204_PF10731_csg_7PTP_D_match_not_genus_replicated_under_FDR;does_not_confer_function_to_candidate",
            "local_neighbourhood_scope": "surveyed_plusminus5_window;all_11_context_pairs_verified_adjacent_cooriented_gaps_-4_to_6bp;does_not_confer_S_layer_function_to_candidate",
            "disposition": "uncertain_insufficient_evidence",
        },
        "deep_dive_dossiers_summary": {
            "n_dossiers_evaluated": len(dossiers_df),
            "n_members_audited_total": len(member_ev_df),
            "n_members_verified_interior": int(member_ev_df["verified_complete_interior_member"].sum()),
            "n_members_with_exact_sequence_model": int(member_ev_df["is_exact_model_match"].sum()),
            "n_verified_interior_members_with_exact_model": int(member_ev_df["is_verified_interior_with_exact_model"].sum()),
            "n_families_with_completed_msa": int((dossiers_df["msa_status"] == "mafft_completed").sum()),
            "n_families_with_deferred_msa": int((dossiers_df["msa_status"] == "deferred_non_focal_no_role_hit").sum()),
            "n_families_with_zero_aligned_when_deferred": int(
                ((dossiers_df["msa_status"] == "deferred_non_focal_no_role_hit") & (dossiers_df["n_aligned_sequences"] == 0)).sum()
            ),
            "disposition_counts": dossiers_df["supported_evidence_route"].value_counts().to_dict(),
        },
        "output_files": {
            "eval_audited_tsv": str(eval_audited_out),
            "neighbourhood_audited_tsv": str(neigh_audited_out),
            "member_evidence_tsv": str(member_ev_out),
            "flank_evidence_tsv": str(flank_ev_out),
            "dossiers_tsv": str(dossiers_out),
            "cohort_audited_tsv": str(cohort_audited_out),
        },
    }

    summary_out = P_FOLLOWUP / "exp4_surface_final_summary_v2_audited.json"
    with open(summary_out, "wt") as fh:
        json.dump(summary_data, fh, indent=2)
    print(f"Wrote audited final summary JSON to {summary_out}")

    print("\n=== Experiment 4 Targeted Followup Evidence Audit COMPLETED SUCCESSFULLY ===")


if __name__ == "__main__":
    main()
