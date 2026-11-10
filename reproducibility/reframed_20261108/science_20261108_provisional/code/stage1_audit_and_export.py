#!/usr/bin/env python3
"""Stage 1: Candidate Sequence/Provenance Audit (Routes A1 & A2) and Lineage-Adjusted Export Analysis (Task B).

Executes on Shaheen III (`workq`) and writes outputs to `/scratch/project/k10018/eq-science-20261108/results`.
Follows all supervisor scientific design rules and code-review fixes (2026-11-08).
"""
from __future__ import annotations

import collections
import glob
import gzip
import hashlib
import json
import math
import os
import re
import sys
import tarfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_science import (
    PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES,
    audit_orf_integrity,
    block_bootstrap_stat,
    classify_candidate_role,
    classify_genome_quality,
    fit_clustered_binomial_glm,
    identify_mgnify_gtdb_duplicates,
    normalize_genus_groups,
    require_file,
    validate_phylum_scope,
    validate_primary_key,
    validated_merge,
)

# Paths on Shaheen
P_ANNOT = Path("/project/k10018/work/eq-annot-20261030")
P_ARCH = Path("/scratch/project/k10018/eq-archaeal-function-20261106")
P_PROT_INT = Path("/scratch/project/k10018/eq-proteome-integration-20261106")
P_OUT_ROOT = Path(os.environ.get("EQ_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261108"))
P_IN = P_OUT_ROOT / "inputs"
P_RES = P_OUT_ROOT / "results"

SPANY_CLASSES = frozenset({"SP", "LIPO", "TAT", "TATLIPO"})


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


def _parse_tmbed_chunk(path_str: str) -> dict[tuple[str, str], tuple[int, bool, int, int, str]]:
    out: dict[tuple[str, str], tuple[int, bool, int, int, str]] = {}
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
            n_h = _count_segments(lab, "Hh")
            n_b = _count_segments(lab, "Bb")
            # Compact topology summary (e.g. S + H count + B count)
            topo = f"S:{int(has_s)}_H:{n_h}_B:{n_b}"
            out[k] = (len(seq), has_s, n_h, n_b, topo)
    return out


def _parse_signalp_chunk(path_str: str) -> dict[tuple[str, str], tuple[str, float, str]]:
    out: dict[tuple[str, str], tuple[str, float, str]] = {}
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
            cs_pos = parts[8].strip() if len(parts) > 8 else ""
            out[k] = (pred, sp_prob, cs_pos)
    return out


def load_all_sp_tm(n_workers: int = 32) -> tuple[dict[tuple[str, str], tuple[int, bool, int, int, str]], dict[tuple[str, str], tuple[str, float, str]]]:
    tm_files = sorted(glob.glob(str(P_ANNOT / "results/sp_tm_analysis/tmbed/chunk_*.pred.gz")))
    sp_files = sorted(glob.glob(str(P_ANNOT / "results/sp_tm_analysis/signalp/chunk_*.prediction_results.txt.gz")))
    if not tm_files or not sp_files:
        raise FileNotFoundError("Missing SignalP or TMbed chunk files in sp_tm_analysis")

    tm_all: dict[tuple[str, str], tuple[int, bool, int, int, str]] = {}
    sp_all: dict[tuple[str, str], tuple[str, float, str]] = {}
    with ProcessPoolExecutor(max_workers=min(n_workers, 40)) as ex:
        for d in ex.map(_parse_tmbed_chunk, tm_files):
            tm_all.update(d)
        for d in ex.map(_parse_signalp_chunk, sp_files):
            sp_all.update(d)
    return tm_all, sp_all


def load_fasta_md5_and_len(fasta_path: Path, target_genes: set[str] | None = None) -> dict[str, tuple[int, str, str]]:
    """Return gene_id -> (aa_len, md5_hex, aa_seq_if_target_else_empty)."""
    require_file(fasta_path, "archaeal proteins FASTA")
    res: dict[str, tuple[int, str, str]] = {}
    cur_id: str | None = None
    chunks: list[str] = []

    def _flush(gid: str | None, seq_parts: list[str]) -> None:
        if gid is None:
            return
        raw = "".join(seq_parts).strip()
        seq = raw[:-1] if raw.endswith("*") else raw
        md5 = hashlib.md5(seq.encode("ascii")).hexdigest()
        keep_seq = seq if (target_genes is None or gid in target_genes) else ""
        res[gid] = (len(seq), md5, keep_seq)

    with open(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                _flush(cur_id, chunks)
                hdr = line[1:].strip().split()[0]
                cur_id = hdr.split("|", 1)[1] if "|" in hdr else hdr
                chunks = []
            else:
                chunks.append(line.strip())
        _flush(cur_id, chunks)
    return res


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


def load_structural_and_domain_summaries() -> dict[str, Any]:
    """Load ColabFold, ESMFold, Foldseek, HHsearch, InterPro, and DeepFRI summaries keyed by family_rep / protein_id."""
    # 1. ColabFold members
    cf_mem: dict[str, dict[str, Any]] = {}
    cf_mem_path = require_file(P_ANNOT / "results/c10_structclust/colabfold/members.tsv.gz", "ColabFold members")
    with gzip.open(cf_mem_path, "rt") as fh:
        for line in fh:
            arm, pid, length, plddt, ptm, sname = line.rstrip("\n").split("\t")
            cf_mem[pid] = {
                "cf_plddt": float(plddt),
                "cf_ptm": float(ptm),
                "cf_struct_name": sname,
            }

    # 2. ESMFold members
    esm_mem: dict[str, dict[str, Any]] = {}
    esm_mem_path = require_file(P_ANNOT / "results/c10_structclust/esmfold/members.tsv.gz", "ESMFold members")
    with gzip.open(esm_mem_path, "rt") as fh:
        for line in fh:
            arm, pid, length, plddt, ptm, sname = line.rstrip("\n").split("\t")
            esm_mem[pid] = {
                "esm_arm": arm,
                "esm_plddt": float(plddt),
                "esm_ptm": float(ptm),
                "esm_struct_name": sname,
            }

    # 3. ColabFold Foldseek PDB & AFDB best hits
    cf_pdb_best: dict[str, dict[str, Any]] = {}
    for fn in ("pdb_whole.tsv.gz", "pdb_regions.tsv.gz"):
        p = require_file(P_ANNOT / "results/c9cf_colabfold" / fn, f"c9cf {fn}")
        with gzip.open(p, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                row = dict(zip(hdr, line.rstrip("\n").split("\t")))
                pid = row["protein_id"]
                try:
                    qtm = float(row["qtmscore"])
                    lddt = float(row["lddt"])
                    ev = float(row["evalue"])
                except ValueError:
                    continue
                prev = cf_pdb_best.get(pid)
                if prev is None or qtm > prev["cf_pdb_qtmscore"]:
                    cf_pdb_best[pid] = {
                        "cf_pdb_target": row["target"],
                        "cf_pdb_region": row["region"],
                        "cf_pdb_qtmscore": qtm,
                        "cf_pdb_lddt": lddt,
                        "cf_pdb_evalue": ev,
                        "cf_pdb_fident": float(row.get("fident", "nan")),
                    }

    cf_afdb_best: dict[str, dict[str, Any]] = {}
    p_af = require_file(P_ANNOT / "results/c9cf_colabfold/afdb_whole.tsv.gz", "c9cf afdb_whole")
    with gzip.open(p_af, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            pid = row["protein_id"]
            try:
                ev = float(row["evalue"])
                bits = float(row["bits"])
            except ValueError:
                continue
            prev = cf_afdb_best.get(pid)
            if prev is None or ev < prev["cf_afdb_evalue"]:
                cf_afdb_best[pid] = {
                    "cf_afdb_target": row["target"],
                    "cf_afdb_evalue": ev,
                    "cf_afdb_bits": bits,
                }

    # 4. HHsearch c4_per_query.tsv.gz (query uses '__' instead of '::')
    c4_hits: dict[str, dict[str, Any]] = {}
    p_c4 = require_file(P_ANNOT / "results/_agent_analysis/c4_c5_c7/c4/c4_per_query.tsv.gz", "c4_per_query")
    with gzip.open(p_c4, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            qid = row["query"].replace("__", "::", 1)
            c4_hits[qid] = {
                "hh_pfam_prob": float(row["pfam_prob"]) if row["pfam_prob"] else 0.0,
                "hh_pdb_prob": float(row["pdb_prob"]) if row["pdb_prob"] else 0.0,
                "hh_best_pfam": row["best_pfam"] or "missing",
                "hh_pfam_id": row["pfam_id"] or "missing",
                "hh_best_pdb": row["best_pdb"] or "missing",
                "hh_pdb_desc": row["pdb_desc"] or "missing",
            }

    # 5. InterPro c5_family_bits.tsv.gz
    c5_hits: dict[str, dict[str, Any]] = {}
    p_c5 = require_file(P_ANNOT / "results/_agent_analysis/c4_c5_c7/c5/c5_family_bits.tsv.gz", "c5_family_bits")
    with gzip.open(p_c5, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            c5_hits[row["family_rep"]] = {
                "interpro_rep_bits": int(row["rep_bits"]),
                "interpro_anymember_bits": int(row["anymember_bits"]),
                "interpro_n_members": int(row["n_members"]),
            }

    return {
        "cf_mem": cf_mem,
        "esm_mem": esm_mem,
        "cf_pdb_best": cf_pdb_best,
        "cf_afdb_best": cf_afdb_best,
        "c4_hits": c4_hits,
        "c5_hits": c5_hits,
    }


def _parse_c6_thermo_batch(args: tuple[str, set[str]]) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    b, thermo_genomes = args
    d = P_ANNOT / "results/c6_cazy_merops_tcdb"
    cazy_gene: dict[str, str] = {}
    mer_gene: dict[str, str] = {}
    tc_gene: dict[str, str] = {}

    dbcan_tar = d / f"batch_{b}.dbcan.tar.gz"
    if dbcan_tar.exists():
        with tarfile.open(dbcan_tar, "r|gz") as t:
            for m in t:
                if m.name.endswith("overview.tsv"):
                    f = (l.decode() for l in t.extractfile(m))
                    next(f, None)
                    for line in f:
                        x = line.rstrip("\n").split("\t")
                        if "|" in x[0]:
                            g, gid = x[0].split("|", 1)
                        else:
                            g, gid = x[0].split("::", 1)[0], x[0]
                        if g not in thermo_genomes:
                            continue
                        fams = []
                        for col in x[2:5]:
                            if col not in ("-", ""):
                                fams.append(col)
                        if fams:
                            cazy_gene[gid] = ";".join(fams)

    mer_path = d / f"batch_{b}.merops.m8.gz"
    if mer_path.exists():
        with gzip.open(mer_path, "rt") as fh:
            for line in fh:
                q, s, pid, ln, qlen, slen, qs, qe, ss, se, ev, bs = line.rstrip("\n").split("\t")
                g, gid = q.split("|", 1) if "|" in q else (q.split("::", 1)[0], q)
                if g not in thermo_genomes or gid in mer_gene:
                    continue
                if float(ev) <= 1e-5 and float(pid) >= 30.0:
                    mer_gene[gid] = f"{s}(pident={pid},e={ev})"

    tc_path = d / f"batch_{b}.tcdb.m8.gz"
    if tc_path.exists():
        with gzip.open(tc_path, "rt") as fh:
            for line in fh:
                q, s, pid, ln, qlen, slen, qs, qe, ss, se, ev, bs = line.rstrip("\n").split("\t")
                g, gid = q.split("|", 1) if "|" in q else (q.split("::", 1)[0], q)
                if g not in thermo_genomes or gid in tc_gene:
                    continue
                if float(ev) <= 1e-5 and float(pid) >= 25.0:
                    tcid = s.split("|")[3] if s.count("|") >= 3 else s
                    tc_gene[gid] = f"{tcid}(pident={pid},e={ev})"

    return cazy_gene, mer_gene, tc_gene


def load_c6_thermo_gene_hits(thermo_genomes: set[str], n_workers: int = 16) -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    batches = [(f"{i:02d}", thermo_genomes) for i in range(16)]
    cazy_all: dict[str, str] = {}
    mer_all: dict[str, str] = {}
    tc_all: dict[str, str] = {}
    with ProcessPoolExecutor(max_workers=min(n_workers, 16)) as ex:
        for c, m, t in ex.map(_parse_c6_thermo_batch, batches):
            cazy_all.update(c)
            mer_all.update(m)
            tc_all.update(t)
    return cazy_all, mer_all, tc_all


def load_orthofinder_gene_to_og(thermo_genomes: set[str]) -> tuple[dict[str, str], pd.DataFrame]:
    """Extract gene -> OG mapping for the 237 EQ Thermoplasmatota genomes and load c7_orthogroups_all538.tsv.gz."""
    og_summary_path = require_file(
        P_ANNOT / "results/_agent_analysis/c4_c5_c7/c7/c7_orthogroups_all538.tsv.gz",
        "c7_orthogroups_all538",
    )
    og_df = pd.read_csv(og_summary_path, sep="\t")
    validate_primary_key(og_df, ["og"], "c7_orthogroups_all538")

    tar_path = require_file(P_ANNOT / "results/c7_thermo_evol/orthofinder_results.tar.gz", "orthofinder_results.tar.gz")
    gene_to_og: dict[str, str] = {}
    with tarfile.open(tar_path, "r:gz") as tf:
        f = tf.extractfile("Orthogroups/Orthogroups.tsv")
        if f is None:
            raise ValueError("Missing Orthogroups/Orthogroups.tsv in archive")
        hdr = f.readline().decode("utf-8").rstrip("\r\n").split("\t")
        eq_col_indices = [i for i, col in enumerate(hdr) if col in thermo_genomes]
        for raw_line in f:
            parts = raw_line.decode("utf-8").rstrip("\r\n").split("\t")
            og_id = parts[0]
            for idx in eq_col_indices:
                if idx < len(parts) and parts[idx]:
                    for g_tok in parts[idx].split(", "):
                        g_tok = g_tok.strip()
                        if g_tok:
                            gene_id = g_tok.split("|", 1)[1] if "|" in g_tok else g_tok
                            gene_to_og[gene_id] = og_id
    return gene_to_og, og_df


def classify_candidate_stratum(
    pathway: str,
    rxn_name: str,
    ec: str,
    sp_or_s: bool,
    tm_any: bool,
    tc_hit: str,
    mer_hit: str,
    cazy_hit: str,
    hh_desc: str,
) -> str:
    """Stratify candidate families into substrate_acquisition_interface, maintenance, or other."""
    txt = f"{pathway} {rxn_name} {ec} {tc_hit} {mer_hit} {cazy_hit} {hh_desc}".lower()
    if any(
        kw in txt
        for kw in (
            "trans-",
            "transport",
            "abc",
            "permease",
            "symport",
            "antiport",
            "uptake",
            "import",
            "export",
            "peptidase",
            "protease",
            "glycosyl",
            "hydrolase",
            "chitin",
            "cellulose",
            "xylan",
            "s-layer",
            "surface",
            "adhesin",
            "pilin",
            "flagell",
            "archaeaellum",
            "extracellular",
        )
    ) or tc_hit != "missing" or mer_hit != "missing" or cazy_hit != "missing" or sp_or_s:
        return "substrate_acquisition_interface"
    if any(
        kw in txt
        for kw in (
            "repair",
            "chaperon",
            "grpe",
            "dnaj",
            "dnak",
            "hsp",
            "superoxide",
            "peroxid",
            "catalase",
            "thioredoxin",
            "glutaredoxin",
            "rad50",
            "rec",
            "uvr",
            "mut",
            "trehalose",
            "betaine",
            "ectoine",
            "osmol",
            "homeostasis",
            "stress",
            "defense",
        )
    ):
        return "maintenance"
    return "other"


def classify_exported_dark_hypothesis(
    sp_any: bool,
    tmbed_s: bool,
    tm_helices: int,
    tm_beta: int,
    hh_pfam: str,
    hh_desc: str,
    interpro_bits: int,
    cazy_hit: str,
    mer_hit: str,
    tc_hit: str,
    flank_summary: str,
) -> tuple[str, str]:
    """Classify independent exported-dark families across biological alternatives without overclaiming mechanism.

    Returns (primary_biological_category, s_layer_or_enzyme_counterexample_note).
    A signal peptide is localization evidence, not evidence of substrate degradation.
    """
    combined = f"{hh_pfam} {hh_desc} {cazy_hit} {mer_hit} {tc_hit} {flank_summary}".lower()
    notes = []
    if "s-layer" in combined or "slh" in combined or "pkd" in combined or "ig-like" in combined or "fn3" in combined:
        notes.append("shows_S_layer_or_surface_repeat_relationship")
    else:
        notes.append("no_known_S_layer_domain_hit")

    if mer_hit != "missing" or "peptidase" in combined or "protease" in combined:
        notes.append("shows_peptidase_relationship")
    elif cazy_hit != "missing" or "glycosid" in combined or "hydrolase" in combined:
        notes.append("shows_CAZyme_or_hydrolase_relationship")
    elif tc_hit != "missing" or "transport" in combined or "abc" in combined or "permease" in combined:
        notes.append("shows_transporter_relationship")
    else:
        notes.append("signal_peptide_localization_only_no_catalytic_peptidase_or_CAZyme_match")

    if mer_hit != "missing" or cazy_hit != "missing" or tc_hit != "missing" or "transport" in combined or "abc" in combined:
        cat = "substrate_acquisition_candidate_conditional"
    elif "s-layer" in combined or "pilin" in combined or "adhes" in combined or "flagell" in combined or tm_beta > 0:
        cat = "surface_organisation_or_adhesion_candidate_conditional"
    elif "chaperon" in combined or "thioredoxin" in combined or "peroxid" in combined or "repair" in combined:
        cat = "exported_maintenance_or_protection_candidate_conditional"
    else:
        cat = "unresolved_exported_accessory_localisation_only"

    return cat, ";".join(notes)


def main() -> None:
    P_RES.mkdir(parents=True, exist_ok=True)
    print("[Stage 1] Loading and validating Thermoplasmatota genome table and quality...", flush=True)

    thermo_df = pd.read_csv(require_file(P_ARCH / "inputs/thermoplasmatota.tsv", "thermoplasmatota.tsv"), sep="\t")
    validate_primary_key(thermo_df, ["genome"], "thermoplasmatota.tsv")
    validate_phylum_scope(thermo_df, "phylum_group", allowed={"Thermoplasmatota"})
    assert len(thermo_df) == 237, f"Expected 237 EQ Thermoplasmatota genomes, got {len(thermo_df)}"

    # Parse GTDB taxonomy ranks
    tax_parsed = [parse_gtdb_lineage(x) for x in thermo_df["gtdb_classification"]]
    for rk, col_name in [("c", "tax_class"), ("o", "tax_order"), ("f", "tax_family"), ("g", "tax_genus")]:
        thermo_df[col_name] = [d[rk] for d in tax_parsed]

    # Classify genome quality with explicit RNA evidence unknown (no rRNA/tRNA columns in CheckM2 table)
    qual_df = classify_genome_quality(
        thermo_df["checkm_completeness"],
        thermo_df["contamination"],
        rna_complete=None,
    )
    for col in ("sensitivity_comp70_contam5_pass", "strict_comp50_contam_lt10_pass", "comp_contam_tier", "mimag_strict_tier"):
        thermo_df[col] = qual_df[col].to_numpy()

    # Join GUNC
    gunc_df = pd.read_csv(require_file(P_ARCH / "inputs/c3_gunc_per_genome.tsv", "c3_gunc_per_genome.tsv"), sep="\t")
    thermo_df = validated_merge(
        thermo_df,
        gunc_df[["genome", "pass.GUNC", "clade_separation_score", "reference_representation_score"]],
        on="genome",
        how="left",
        validate="one_to_one",
        require_all_left=True,
        left_name="thermoplasmatota",
        right_name="c3_gunc",
    )

    # Join comparator AAI/ANI pairs (`v2_thermo_relatives_pairs.tsv`)
    pairs_path = require_file(P_IN / "v2_thermo_relatives_pairs.tsv", "v2_thermo_relatives_pairs.tsv")
    pairs_df = pd.read_csv(pairs_path, sep="\t")
    thermo_df = validated_merge(
        thermo_df,
        pairs_df[["genome", "aai_nearest", "aai", "aai_source", "aai_nearest_desert", "aai_desert"]],
        on="genome",
        how="left",
        validate="one_to_one",
        require_all_left=True,
        left_name="thermoplasmatota",
        right_name="v2_thermo_relatives_pairs",
    )
    thermo_df["has_external_comparator_aai65"] = pd.to_numeric(thermo_df["aai"], errors="coerce").fillna(0.0) >= 65.0

    # Build non-pooling genus/lineage blocks (unclassified genera become unique singletons!)
    genus_blocks, is_unclass_genus = normalize_genus_groups(thermo_df["genome"].tolist(), thermo_df["tax_genus"].tolist())
    thermo_df["genus_singleton_block"] = genus_blocks
    thermo_df["is_unclassified_genus"] = is_unclass_genus

    thermo_genomes_set = set(thermo_df["genome"])
    print(f"[Stage 1] Verified 237 Thermoplasmatota genomes ({int(thermo_df['sensitivity_comp70_contam5_pass'].sum())} pass 70/5 sensitivity).", flush=True)

    # Load candidate_reaction_gene_review.tsv.gz first so we know which genes need full AA sequence strings retained
    cand_path = require_file(P_PROT_INT / "results/candidate_reaction_gene_review.tsv.gz", "candidate_reaction_gene_review")
    cand_df = pd.read_csv(cand_path, sep="\t")
    cand_genes_set = set(cand_df["candidate_gene"].astype(str))
    print(f"[Stage 1] Loaded candidate_reaction_gene_review: {len(cand_df)} total rows, {int(cand_df['review_priority'].sum())} priority rows, {len(cand_genes_set)} unique candidate genes.", flush=True)

    # Load FASTA MD5 & lengths for all archaeal proteins (retaining full AA sequence for candidate genes + dark family reps)
    print("[Stage 1] Auditing FASTA sequences and MD5 digests from all_archaea.faa...", flush=True)
    fasta_info = load_fasta_md5_and_len(P_ARCH / "proteins/all_archaea.faa", target_genes=cand_genes_set)

    # Load archaeal_genes.tsv filtered to the 237 Thermoplasmatota genomes
    print("[Stage 1] Loading archaeal_genes.tsv for 237 Thermoplasmatota genomes...", flush=True)
    arch_genes_path = require_file(P_ARCH / "inputs/archaeal_genes.tsv", "archaeal_genes.tsv")
    genes_df = pd.read_csv(arch_genes_path, sep="\t", low_memory=False)
    genes_df = genes_df[genes_df["genome"].isin(thermo_genomes_set)].copy().reset_index(drop=True)
    validate_primary_key(genes_df, ["gene"], "thermoplasmatota_archaeal_genes")
    assert len(genes_df) == 777177, f"Expected 777,177 Thermoplasmatota genes, got {len(genes_df)}"

    # Compute exact within-contig gene rank and contig_n_genes
    genes_df["start_num"] = pd.to_numeric(genes_df["start"], errors="coerce")
    genes_df["end_num"] = pd.to_numeric(genes_df["end"], errors="coerce")
    genes_df = genes_df.sort_values(["genome", "contig", "start_num", "end_num", "gene"]).reset_index(drop=True)
    genes_df["contig_rank"] = genes_df.groupby(["genome", "contig"]).cumcount() + 1
    genes_df["contig_n_genes"] = genes_df.groupby(["genome", "contig"])["gene"].transform("count")

    # Verify FASTA presence, AA length, MD5 match, and coordinate validity
    fasta_len = [fasta_info.get(g, (-1, "", ""))[0] for g in genes_df["gene"]]
    fasta_md5 = [fasta_info.get(g, (-1, "", ""))[1] for g in genes_df["gene"]]
    genes_df["fasta_aa_length"] = fasta_len
    genes_df["fasta_md5"] = fasta_md5
    genes_df["has_full_seq"] = (genes_df["fasta_aa_length"] > 0) & (genes_df["fasta_aa_length"] == genes_df["length"]) & (genes_df["fasta_md5"] == genes_df["md5"])
    genes_df["coords_valid"] = (
        genes_df["start_num"].notna()
        & genes_df["end_num"].notna()
        & (genes_df["start_num"] >= 1)
        & (genes_df["end_num"] > genes_df["start_num"])
        & (genes_df["strand"].isin([1, -1, "1", "-1", "+", "-"]))
    )

    # Run strict ORF & contig-edge integrity audit
    genes_df = audit_orf_integrity(
        genes_df,
        partial_col="partial",
        internal_stop_col="internal_stop",
        length_col="length",
        contig_rank_col="contig_rank",
        contig_n_genes_col="contig_n_genes",
        contig_in_cat_col="contig_in_catalogue",
        has_full_seq_col="has_full_seq",
        coords_valid_col="coords_valid",
        edge_gene_buffer=2,
    )

    # Load SignalP 6.0h & TMbed predictions in parallel
    print("[Stage 1] Parsing 40 SignalP 6.0h and TMbed chunks in parallel...", flush=True)
    tm_all, sp_all = load_all_sp_tm(n_workers=32)

    # Map SignalP & TMbed onto the 777,177 Thermoplasmatota genes
    g_list = genes_df["gene"].tolist()
    sp_pred_col = []
    sp_prob_col = []
    sp_cs_col = []
    tm_s_col = []
    tm_h_col = []
    tm_b_col = []
    tm_topo_col = []
    has_sp_tm = []

    for g in g_list:
        kt = ("T", g)
        kd = ("D", g)
        sp_entry = sp_all.get(kt) or sp_all.get(kd)
        tm_entry = tm_all.get(kt) or tm_all.get(kd)
        if sp_entry is not None and tm_entry is not None:
            has_sp_tm.append(True)
            sp_pred_col.append(sp_entry[0])
            sp_prob_col.append(sp_entry[1])
            sp_cs_col.append(sp_entry[2] if sp_entry[2] else "none")
            tm_s_col.append(bool(tm_entry[1]))
            tm_h_col.append(int(tm_entry[2]))
            tm_b_col.append(int(tm_entry[3]))
            tm_topo_col.append(tm_entry[4])
        else:
            has_sp_tm.append(False)
            sp_pred_col.append("missing")
            sp_prob_col.append(np.nan)
            sp_cs_col.append("missing")
            tm_s_col.append(False)
            tm_h_col.append(0)
            tm_b_col.append(0)
            tm_topo_col.append("missing")

    genes_df["has_sp_tm_prediction"] = has_sp_tm
    genes_df["signalp_class"] = sp_pred_col
    genes_df["signalp_prob"] = sp_prob_col
    genes_df["signalp_cs"] = sp_cs_col
    genes_df["sp_any"] = genes_df["signalp_class"].isin(SPANY_CLASSES)
    genes_df["sec_spi"] = genes_df["signalp_class"] == "SP"
    genes_df["sec_spii"] = genes_df["signalp_class"] == "LIPO"
    genes_df["tat"] = genes_df["signalp_class"].isin({"TAT", "TATLIPO"})
    genes_df["pilin"] = genes_df["signalp_class"] == "PILIN"
    genes_df["tmbed_S"] = tm_s_col
    genes_df["sp_or_S"] = genes_df["sp_any"] | genes_df["tmbed_S"]
    genes_df["sp_and_S"] = genes_df["sp_any"] & genes_df["tmbed_S"]
    genes_df["tm_helices"] = tm_h_col
    genes_df["tm_beta"] = tm_b_col
    genes_df["tm_any"] = (genes_df["tm_helices"] + genes_df["tm_beta"]) >= 1
    genes_df["tmbed_topology"] = tm_topo_col

    # Map family memberships (member_fam_final.tsv, ctrl_member_fam.tsv, fam_phylum.tsv)
    print("[Stage 1] Loading family memberships and OrthoFinder 538-genome orthogroups...", flush=True)
    mem2fam: dict[str, str] = {}
    with open(require_file(P_ANNOT / "inputs/fam/member_fam_final.tsv", "member_fam_final.tsv")) as fh:
        for line in fh:
            m, r = line.rstrip("\n").split("\t")[:2]
            mem2fam[m] = r
    ctrl2fam: dict[str, str] = {}
    with open(require_file(P_ANNOT / "inputs/fam/ctrl_member_fam.tsv", "ctrl_member_fam.tsv")) as fh:
        for line in fh:
            m, r = line.rstrip("\n").split("\t")[:2]
            ctrl2fam[m] = r
    fam_phylum: dict[str, tuple[int, str]] = {}
    with open(require_file(P_ANNOT / "inputs/fam/fam_phylum.tsv", "fam_phylum.tsv")) as fh:
        next(fh)
        for line in fh:
            r, n_p, phy = line.rstrip("\n").split("\t")
            fam_phylum[r] = (int(n_p), phy)

    # Note: member_fam_final keys use `cat_id` or `example_id` or `seq_cluster_rep` (format `sample::contig_gene`)
    def _resolve_dark_fam(row_cat_id: str, row_seq_rep: str, row_cat_rep: str) -> str:
        for cand_k in (row_cat_id, row_seq_rep, row_cat_rep):
            if isinstance(cand_k, str) and cand_k in mem2fam:
                return mem2fam[cand_k]
        if isinstance(row_seq_rep, str) and row_seq_rep and row_seq_rep != "nan":
            return row_seq_rep
        if isinstance(row_cat_rep, str) and row_cat_rep and row_cat_rep != "nan":
            return row_cat_rep
        return "unassigned"

    genes_df["family_id"] = [
        _resolve_dark_fam(str(c), str(s), str(cr))
        for c, s, cr in zip(genes_df["cat_id"], genes_df["seq_cluster_rep"], genes_df["cat_cluster_rep"])
    ]

    gene_to_og, og_df = load_orthofinder_gene_to_og(thermo_genomes_set)
    genes_df["orthogroup"] = [gene_to_og.get(g, "unassigned") for g in genes_df["gene"]]
    og_comp_map = dict(zip(og_df["og"], og_df["n_comparator_genomes"]))
    og_des_map = dict(zip(og_df["og"], og_df["n_desert_genomes"]))
    genes_df["og_n_comparator_genomes"] = [int(og_comp_map.get(og, 0)) for og in genes_df["orthogroup"]]
    genes_df["og_n_desert_genomes"] = [int(og_des_map.get(og, 0)) for og in genes_df["orthogroup"]]

    # Join genome-level metadata onto genes_df
    genes_df = validated_merge(
        genes_df,
        thermo_df[
            [
                "genome",
                "tax_class",
                "tax_order",
                "tax_family",
                "tax_genus",
                "genus_singleton_block",
                "is_unclassified_genus",
                "checkm_completeness",
                "contamination",
                "sensitivity_comp70_contam5_pass",
                "comp_contam_tier",
                "mimag_strict_tier",
                "pass.GUNC",
                "has_external_comparator_aai65",
            ]
        ],
        on="genome",
        how="left",
        validate="many_to_one",
        require_all_left=True,
        left_name="thermoplasmatota_genes",
        right_name="thermoplasmatota_genomes",
    )

    # Load structural, HHsearch, InterPro, and C6 (CAZy/MEROPS/TCDB) annotations
    print("[Stage 1] Loading structural (ColabFold/ESMFold/Foldseek), HHsearch, InterPro, and CAZy/MEROPS/TCDB evidence...", flush=True)
    struct_dom = load_structural_and_domain_summaries()
    cazy_gene, mer_gene, tc_gene = load_c6_thermo_gene_hits(thermo_genomes_set, n_workers=16)

    # Build fast flank lookup: (genome, contig) -> list of (contig_rank, gene, dark, family_id, og, cazy, mer, tc)
    print("[Stage 1] Building contig neighborhood index for flanking gene audit...", flush=True)
    contig_groups: dict[tuple[str, str], list[tuple[int, str, Any, str, str]]] = collections.defaultdict(list)
    for r_idx, g_id, dk, fid, og_id, g_nm, c_nm in zip(
        genes_df["contig_rank"],
        genes_df["gene"],
        genes_df["dark"],
        genes_df["family_id"],
        genes_df["orthogroup"],
        genes_df["genome"],
        genes_df["contig"],
    ):
        tag_parts = []
        if g_id in tc_gene:
            tag_parts.append(f"TC:{tc_gene[g_id]}")
        if g_id in mer_gene:
            tag_parts.append(f"MER:{mer_gene[g_id]}")
        if g_id in cazy_gene:
            tag_parts.append(f"CAZy:{cazy_gene[g_id]}")
        if fid in struct_dom["c4_hits"] and struct_dom["c4_hits"][fid]["hh_best_pfam"] != "missing":
            tag_parts.append(f"Pfam:{struct_dom['c4_hits'][fid]['hh_best_pfam']}")
        annot_tag = "|".join(tag_parts) if tag_parts else ("dark" if pd.notna(dk) and float(dk) == 1.0 else "annotated_other")
        contig_groups[(g_nm, c_nm)].append((int(r_idx), g_id, dk, fid, annot_tag))

    def _summarize_flanks(genome: str, contig: str, rank: int, window: int = 5) -> tuple[int, int, str]:
        entries = contig_groups.get((genome, contig), [])
        flanks = [e for e in entries if e[0] != rank and abs(e[0] - rank) <= window]
        n_flank = len(flanks)
        n_known_flank = sum(1 for e in flanks if pd.notna(e[2]) and float(e[2]) == 0.0)
        desc = ";".join(f"{e[0]-rank:+d}:{e[4]}" for e in flanks) if flanks else "no_flanking_genes_on_contig"
        return n_flank, n_known_flank, desc

    # =========================================================================
    # PART A1: Audit of the 3,766 filtered gap-filling proposals (and all 7,381 rows)
    # =========================================================================
    print("[Stage 1] Running Part A1: Sequence/coordinate/provenance audit of candidate_reaction_gene_review...", flush=True)
    gene_lookup = genes_df.set_index("gene")
    cand_audit_rows = []
    for row in cand_df.itertuples(index=False):
        cg = str(row.candidate_gene)
        if cg not in gene_lookup.index:
            raise ValueError(f"Candidate gene {cg} missing from audited Thermoplasmatota gene table")
        g_rec = gene_lookup.loc[cg]
        fid = str(g_rec["family_id"])
        seq_rep = str(g_rec["seq_cluster_rep"]) if pd.notna(g_rec["seq_cluster_rep"]) else "missing"

        # Structural & domain lookups (check both family_id and seq_cluster_rep)
        cf_m = struct_dom["cf_mem"].get(fid) or struct_dom["cf_mem"].get(seq_rep)
        esm_m = struct_dom["esm_mem"].get(fid) or struct_dom["esm_mem"].get(seq_rep)
        cf_p = struct_dom["cf_pdb_best"].get(fid) or struct_dom["cf_pdb_best"].get(seq_rep)
        cf_a = struct_dom["cf_afdb_best"].get(fid) or struct_dom["cf_afdb_best"].get(seq_rep)
        hh_m = struct_dom["c4_hits"].get(fid) or struct_dom["c4_hits"].get(seq_rep)
        ip_m = struct_dom["c5_hits"].get(fid) or struct_dom["c5_hits"].get(seq_rep)

        n_flk, n_known_flk, flk_desc = _summarize_flanks(
            str(g_rec["genome"]), str(g_rec["contig"]), int(g_rec["contig_rank"]), window=5
        )
        tc_val = tc_gene.get(cg, "missing")
        mer_val = mer_gene.get(cg, "missing")
        cazy_val = cazy_gene.get(cg, "missing")
        hh_desc = hh_m["hh_pdb_desc"] if hh_m else "missing"

        stratum = classify_candidate_stratum(
            str(row.pathway),
            str(row.reaction_name),
            str(row.ec),
            bool(g_rec["sp_or_S"]),
            bool(g_rec["tm_any"]),
            tc_val,
            mer_val,
            cazy_val,
            hh_desc,
        )

        aa_len, aa_md5, aa_seq = fasta_info.get(cg, (-1, "missing", "missing"))
        cand_audit_rows.append(
            {
                "genome": str(row.genome),
                "candidate_gene": cg,
                "pathway": str(row.pathway),
                "missing_reaction": str(row.missing_reaction),
                "reaction_name": str(row.reaction_name),
                "ec": str(row.ec),
                "review_priority": bool(row.review_priority),
                "inherited_catalogue_dark": bool(row.dark),
                "inherited_catalogue_dark_strict": bool(float(row.dark_strict) == 1.0) if pd.notna(row.dark_strict) else False,
                "family_id": fid,
                "seq_cluster_rep": seq_rep,
                "cat_cluster_rep": str(g_rec["cat_cluster_rep"]),
                "orthogroup": str(g_rec["orthogroup"]),
                "og_n_desert_genomes": int(g_rec["og_n_desert_genomes"]),
                "og_n_comparator_genomes": int(g_rec["og_n_comparator_genomes"]),
                "contig": str(g_rec["contig"]),
                "start": int(g_rec["start_num"]),
                "end": int(g_rec["end_num"]),
                "strand": str(g_rec["strand"]),
                "partial_flag": str(g_rec["partial"]),
                "internal_stop": int(g_rec["internal_stop"]),
                "orf_status": str(g_rec["orf_status"]),
                "coord_aa_length": int(g_rec["length"]),
                "fasta_aa_length": int(aa_len),
                "md5": str(aa_md5),
                "has_full_seq": bool(g_rec["has_full_seq"]),
                "coords_valid": bool(g_rec["coords_valid"]),
                "contig_in_catalogue": bool(g_rec["contig_in_catalogue"] == 1),
                "contig_len_match": bool(g_rec["contig_len_match"] == 1),
                "contig_rank": int(g_rec["contig_rank"]),
                "contig_n_genes": int(g_rec["contig_n_genes"]),
                "dist_from_contig_end_genes": float(g_rec["dist_from_contig_end_genes"]),
                "contig_edge_status": str(g_rec["contig_edge_status"]),
                "candidate_promotion_eligible": bool(g_rec["candidate_promotion_eligible"]),
                "n_flank_genes_pm5": int(n_flk),
                "n_known_flank_genes_pm5": int(n_known_flk),
                "flank_annotations_pm5": flk_desc,
                "context_score": float(row.context_score),
                "subthreshold_sequence_hit": bool(row.subthreshold_sequence_hit),
                "sequence_bitscore": float(row.sequence_bitscore) if pd.notna(row.sequence_bitscore) else np.nan,
                "sequence_coverage": float(row.sequence_coverage) if pd.notna(row.sequence_coverage) else np.nan,
                "signalp_class": str(g_rec["signalp_class"]),
                "signalp_prob": float(g_rec["signalp_prob"]) if pd.notna(g_rec["signalp_prob"]) else np.nan,
                "sp_any": bool(g_rec["sp_any"]),
                "tmbed_S": bool(g_rec["tmbed_S"]),
                "sp_or_S": bool(g_rec["sp_or_S"]),
                "tm_helices": int(g_rec["tm_helices"]),
                "tm_beta": int(g_rec["tm_beta"]),
                "tm_any": bool(g_rec["tm_any"]),
                "tmbed_topology": str(g_rec["tmbed_topology"]),
                "colabfold_plddt": cf_m["cf_plddt"] if cf_m else np.nan,
                "colabfold_ptm": cf_m["cf_ptm"] if cf_m else np.nan,
                "colabfold_pdb_target": cf_p["cf_pdb_target"] if cf_p else "missing",
                "colabfold_pdb_qtmscore": cf_p["cf_pdb_qtmscore"] if cf_p else np.nan,
                "colabfold_pdb_lddt": cf_p["cf_pdb_lddt"] if cf_p else np.nan,
                "colabfold_afdb_target": cf_a["cf_afdb_target"] if cf_a else "missing",
                "colabfold_afdb_evalue": cf_a["cf_afdb_evalue"] if cf_a else np.nan,
                "esmfold_plddt": esm_m["esm_plddt"] if esm_m else np.nan,
                "esmfold_ptm": esm_m["esm_ptm"] if esm_m else np.nan,
                "hhsearch_pfam_prob": hh_m["hh_pfam_prob"] if hh_m else np.nan,
                "hhsearch_best_pfam": hh_m["hh_best_pfam"] if hh_m else "missing",
                "hhsearch_pfam_id": hh_m["hh_pfam_id"] if hh_m else "missing",
                "hhsearch_pdb_prob": hh_m["hh_pdb_prob"] if hh_m else np.nan,
                "hhsearch_best_pdb": hh_m["hh_best_pdb"] if hh_m else "missing",
                "hhsearch_pdb_desc": hh_desc,
                "interpro_anymember_bits": ip_m["interpro_anymember_bits"] if ip_m else 0,
                "tcdb_hit": tc_val,
                "merops_hit": mer_val,
                "cazy_hit": cazy_val,
                "tax_class": str(g_rec["tax_class"]),
                "tax_order": str(g_rec["tax_order"]),
                "tax_genus": str(g_rec["tax_genus"]),
                "genus_singleton_block": str(g_rec["genus_singleton_block"]),
                "checkm_completeness": float(g_rec["checkm_completeness"]),
                "contamination": float(g_rec["contamination"]),
                "sensitivity_comp70_contam5_pass": bool(g_rec["sensitivity_comp70_contam5_pass"]),
                "comp_contam_tier": str(g_rec["comp_contam_tier"]),
                "mimag_strict_tier": str(g_rec["mimag_strict_tier"]),
                "pass_GUNC": bool(g_rec["pass.GUNC"]),
                "candidate_stratum": stratum,
                "candidate_role": classify_candidate_role(fid),
                "aa_sequence": aa_seq,
            }
        )

    cand_audit_df = pd.DataFrame(cand_audit_rows)
    cand_audit_df.to_csv(P_RES / "task_a_route1_gapfill_proposals_audit.tsv.gz", sep="\t", index=False, compression="gzip")

    # Filter to the 3,766 review_priority proposals (incl. 3,387 dark) and collapse to FAMILY level
    prio_df = cand_audit_df[cand_audit_df["review_priority"]].copy()
    assert len(prio_df) == 3766, f"Expected 3,766 priority proposals, got {len(prio_df)}"
    assert int(prio_df["inherited_catalogue_dark"].sum()) == 3387, f"Expected 3,387 dark priority proposals, got {int(prio_df['inherited_catalogue_dark'].sum())}"

    # Family-level aggregation for Route A1
    fam_rows_a1 = []
    for fid, grp in prio_df.groupby("family_id", sort=False):
        n_rows = len(grp)
        n_genes = grp["candidate_gene"].nunique()
        n_genomes = grp["genome"].nunique()
        n_orders = grp["tax_order"].nunique()
        n_genus_blocks = grp["genus_singleton_block"].nunique()
        frac_promo = float(grp["candidate_promotion_eligible"].mean())
        frac_interior = float((grp["contig_edge_status"] == "interior").mean())
        frac_sp_or_s = float(grp["sp_or_S"].mean())
        frac_tm = float(grp["tm_any"].mean())
        max_ctx = float(grp["context_score"].max())
        max_bit = float(grp["sequence_bitscore"].max()) if grp["sequence_bitscore"].notna().any() else np.nan
        cf_plddt = float(grp["colabfold_plddt"].max()) if grp["colabfold_plddt"].notna().any() else np.nan
        cf_qtm = float(grp["colabfold_pdb_qtmscore"].max()) if grp["colabfold_pdb_qtmscore"].notna().any() else np.nan
        esm_plddt = float(grp["esmfold_plddt"].max()) if grp["esmfold_plddt"].notna().any() else np.nan
        hh_prob = float(grp["hhsearch_pfam_prob"].max()) if grp["hhsearch_pfam_prob"].notna().any() else np.nan
        stratum = grp["candidate_stratum"].mode().iloc[0]
        is_dark = bool(grp["inherited_catalogue_dark"].any())

        # Conditional evidence count (never assigning mechanism from context alone):
        # 1. Complete interior ORF & retained contig (frac_promo == 1.0 and frac_interior >= 0.5)
        # 2. Independent lineage recurrence (n_genus_blocks >= 2)
        # 3. Structural confidence & fold match (cf_plddt >= 0.70 or esm_plddt >= 0.70, with cf_qtm >= 0.50 or hh_prob >= 80)
        # 4. Domain / subthreshold sequence / transporter / peptidase / CAZyme evidence
        has_struct_fold = bool((pd.notna(cf_qtm) and cf_qtm >= 0.50) or (pd.notna(hh_prob) and hh_prob >= 80.0))
        has_seq_or_domain = bool(
            pd.notna(max_bit)
            or (grp["tcdb_hit"] != "missing").any()
            or (grp["merops_hit"] != "missing").any()
            or (grp["cazy_hit"] != "missing").any()
            or (grp["interpro_anymember_bits"] > 0).any()
        )
        convergent_score = (
            (1.0 if frac_promo >= 0.8 and frac_interior >= 0.5 else 0.0)
            + (1.0 if n_genus_blocks >= 2 else 0.0)
            + (1.0 if has_struct_fold else 0.0)
            + (1.0 if has_seq_or_domain else 0.0)
        )

        fam_rows_a1.append(
            {
                "selection_route": "Route_A1_gapfill_selected",
                "family_id": fid,
                "candidate_role": classify_candidate_role(fid),
                "inherited_catalogue_dark": is_dark,
                "candidate_stratum": stratum,
                "n_proposal_rows": n_rows,
                "n_distinct_genes": n_genes,
                "n_distinct_genomes": n_genomes,
                "n_distinct_orders": n_orders,
                "n_independent_genus_blocks": n_genus_blocks,
                "median_aa_length": float(grp["coord_aa_length"].median()),
                "frac_candidate_promotion_eligible": frac_promo,
                "frac_contig_interior": frac_interior,
                "frac_sp_or_S": frac_sp_or_s,
                "frac_tm_any": frac_tm,
                "max_context_score": max_ctx,
                "max_subthreshold_bitscore": max_bit,
                "colabfold_plddt": cf_plddt,
                "colabfold_pdb_target": grp["colabfold_pdb_target"].iloc[0],
                "colabfold_pdb_qtmscore": cf_qtm,
                "esmfold_plddt": esm_plddt,
                "hhsearch_pfam_prob": hh_prob,
                "hhsearch_best_pfam": grp["hhsearch_best_pfam"].iloc[0],
                "hhsearch_pdb_desc": grp["hhsearch_pdb_desc"].iloc[0],
                "interpro_anymember_bits": int(grp["interpro_anymember_bits"].max()),
                "tcdb_hit": grp["tcdb_hit"].iloc[0],
                "merops_hit": grp["merops_hit"].iloc[0],
                "cazy_hit": grp["cazy_hit"].iloc[0],
                "proposed_reactions": ";".join(sorted(grp["missing_reaction"].unique())[:5]),
                "proposed_pathways": ";".join(sorted(grp["pathway"].unique())[:5]),
                "convergent_evidence_score": convergent_score,
                "evidence_decision_status": (
                    "conditional_candidate_for_targeted_followup"
                    if (convergent_score >= 3.0 and frac_promo >= 0.8)
                    else "context_or_incomplete_only_not_promoted"
                ),
            }
        )

    fam_a1_df = pd.DataFrame(fam_rows_a1).sort_values(
        ["inherited_catalogue_dark", "convergent_evidence_score", "n_independent_genus_blocks", "max_context_score"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    fam_a1_df.to_csv(P_RES / "task_a_route1_gapfill_family_ranking.tsv", sep="\t", index=False)

    # Construct matched controls for Route A1 (matched on length bin, comp_contam_tier, tax_order, and contig_edge_status)
    non_cand_pool = genes_df[
        (~genes_df["gene"].isin(cand_genes_set))
        & (pd.to_numeric(genes_df["dark"], errors="coerce") == 1.0)
        & genes_df["has_sp_tm_prediction"]
    ].copy()
    non_cand_pool["len_bin"] = np.where(
        non_cand_pool["length"] < 100,
        "<100",
        np.where(non_cand_pool["length"] < 250, "100-249", ">=250"),
    )
    prio_dark_df = prio_df[prio_df["inherited_catalogue_dark"]].copy()
    prio_dark_df["len_bin"] = np.where(
        prio_dark_df["coord_aa_length"] < 100,
        "<100",
        np.where(prio_dark_df["coord_aa_length"] < 250, "100-249", ">=250"),
    )

    rng_ctrl = np.random.default_rng(20261108)
    ctrl_rows = []
    for (ord_k, lb_k, qt_k, ed_k), c_grp in prio_dark_df.groupby(
        ["tax_order", "len_bin", "comp_contam_tier", "contig_edge_status"], sort=False
    ):
        n_cand_stratum = c_grp["candidate_gene"].nunique()
        pool_sub = non_cand_pool[
            (non_cand_pool["tax_order"] == ord_k)
            & (non_cand_pool["len_bin"] == lb_k)
            & (non_cand_pool["comp_contam_tier"] == qt_k)
            & (non_cand_pool["contig_edge_status"] == ed_k)
        ]
        if len(pool_sub) > 0:
            take_n = min(len(pool_sub), n_cand_stratum)
            sampled = pool_sub.iloc[rng_ctrl.choice(len(pool_sub), size=take_n, replace=False)]
            ctrl_sp_rate = float(sampled["sp_or_S"].mean())
            ctrl_promo_rate = float(sampled["candidate_promotion_eligible"].mean())
        else:
            take_n = 0
            ctrl_sp_rate = np.nan
            ctrl_promo_rate = np.nan
        ctrl_rows.append(
            {
                "tax_order": ord_k,
                "length_bin": lb_k,
                "comp_contam_tier": qt_k,
                "contig_edge_status": ed_k,
                "n_unique_dark_candidates": n_cand_stratum,
                "n_matched_dark_controls": take_n,
                "candidate_sp_or_S_rate": float(c_grp.drop_duplicates("candidate_gene")["sp_or_S"].mean()),
                "matched_control_sp_or_S_rate": ctrl_sp_rate,
                "candidate_promotion_eligible_rate": float(
                    c_grp.drop_duplicates("candidate_gene")["candidate_promotion_eligible"].mean()
                ),
                "matched_control_promotion_eligible_rate": ctrl_promo_rate,
            }
        )
    pd.DataFrame(ctrl_rows).to_csv(P_RES / "task_a_route1_matched_controls.tsv", sep="\t", index=False)

    # =========================================================================
    # PART A2: Independent Full Exported-Dark-Family Selection across All 237 Genomes
    # =========================================================================
    print("[Stage 1] Running Part A2: Independent full exported-dark-family selection across 237 Thermoplasmatota genomes...", flush=True)
    dark_genes_df = genes_df[
        (pd.to_numeric(genes_df["dark"], errors="coerce") == 1.0)
        & (genes_df["family_id"] != "unassigned")
    ].copy()

    fam_rows_a2 = []
    a1_fam_set = set(fam_a1_df["family_id"])
    for fid, grp in dark_genes_df.groupby("family_id", sort=False):
        n_genes = len(grp)
        n_genomes = grp["genome"].nunique()
        if n_genomes < 3:
            continue
        frac_sp_any = float(grp["sp_any"].mean())
        frac_tmbed_s = float(grp["tmbed_S"].mean())
        frac_sp_or_s = float(grp["sp_or_S"].mean())
        if frac_sp_or_s < 0.50:
            continue

        n_orders = grp["tax_order"].nunique()
        n_genus_blocks = grp["genus_singleton_block"].nunique()
        frac_promo = float(grp["candidate_promotion_eligible"].mean())
        frac_interior = float((grp["contig_edge_status"] == "interior").mean())
        med_len = float(grp["length"].median())
        mean_h = float(grp["tm_helices"].mean())
        mean_b = float(grp["tm_beta"].mean())

        seq_rep = str(grp["seq_cluster_rep"].iloc[0])
        cf_m = struct_dom["cf_mem"].get(fid) or struct_dom["cf_mem"].get(seq_rep)
        esm_m = struct_dom["esm_mem"].get(fid) or struct_dom["esm_mem"].get(seq_rep)
        cf_p = struct_dom["cf_pdb_best"].get(fid) or struct_dom["cf_pdb_best"].get(seq_rep)
        cf_a = struct_dom["cf_afdb_best"].get(fid) or struct_dom["cf_afdb_best"].get(seq_rep)
        hh_m = struct_dom["c4_hits"].get(fid) or struct_dom["c4_hits"].get(seq_rep)
        ip_m = struct_dom["c5_hits"].get(fid) or struct_dom["c5_hits"].get(seq_rep)

        # Check C6 hits across all member genes of this family
        member_genes = grp["gene"].tolist()
        tc_hits = [tc_gene[g] for g in member_genes if g in tc_gene]
        mer_hits = [mer_gene[g] for g in member_genes if g in mer_gene]
        cazy_hits = [cazy_gene[g] for g in member_genes if g in cazy_gene]
        tc_str = tc_hits[0] if tc_hits else "missing"
        mer_str = mer_hits[0] if mer_hits else "missing"
        cazy_str = cazy_hits[0] if cazy_hits else "missing"

        # Flank completeness & summary on representative locus
        rep_row = grp.sort_values(["candidate_promotion_eligible", "checkm_completeness"], ascending=[False, False]).iloc[0]
        n_flk, n_known_flk, flk_desc = _summarize_flanks(
            str(rep_row["genome"]), str(rep_row["contig"]), int(rep_row["contig_rank"]), window=5
        )

        hh_pfam = hh_m["hh_best_pfam"] if hh_m else "missing"
        hh_desc = hh_m["hh_pdb_desc"] if hh_m else "missing"
        ip_bits = ip_m["interpro_anymember_bits"] if ip_m else 0

        bio_cat, counterexample_note = classify_exported_dark_hypothesis(
            frac_sp_any >= 0.5,
            frac_tmbed_s >= 0.5,
            int(round(mean_h)),
            int(round(mean_b)),
            hh_pfam,
            hh_desc,
            ip_bits,
            cazy_str,
            mer_str,
            tc_str,
            flk_desc,
        )

        cf_plddt = cf_m["cf_plddt"] if cf_m else np.nan
        esm_plddt = esm_m["esm_plddt"] if esm_m else np.nan
        best_plddt = np.nanmax([cf_plddt, esm_plddt]) if (pd.notna(cf_plddt) or pd.notna(esm_plddt)) else np.nan
        cf_qtm = cf_p["cf_pdb_qtmscore"] if cf_p else np.nan
        hh_prob = hh_m["hh_pfam_prob"] if hh_m else np.nan

        # Independent quality/recurrence/structure priority score (not based on gap-filling context!)
        indep_score = (
            (2.0 if frac_promo >= 0.80 and frac_interior >= 0.70 else 0.0)
            + (2.0 if n_genus_blocks >= 4 else (1.0 if n_genus_blocks >= 2 else 0.0))
            + (1.5 if (frac_sp_any >= 0.50 and frac_tmbed_s >= 0.50) else 0.5)
            + (1.5 if (pd.notna(best_plddt) and best_plddt >= 0.75 and med_len >= 100) else 0.0)
            + (1.0 if n_flk >= 4 else 0.0)
        )

        fam_rows_a2.append(
            {
                "selection_route": "Route_A2_independent_exported_dark",
                "family_id": fid,
                "in_route_a1_gapfill_priority": fid in a1_fam_set,
                "primary_biological_category": bio_cat,
                "s_layer_and_enzyme_counterexample_audit": counterexample_note,
                "n_member_genes_thermo237": n_genes,
                "n_distinct_genomes": n_genomes,
                "n_distinct_orders": n_orders,
                "n_independent_genus_blocks": n_genus_blocks,
                "median_aa_length": med_len,
                "frac_candidate_promotion_eligible": frac_promo,
                "frac_contig_interior": frac_interior,
                "frac_sp_any": frac_sp_any,
                "frac_tmbed_S": frac_tmbed_s,
                "frac_sp_or_S": frac_sp_or_s,
                "mean_tm_helices": mean_h,
                "mean_tm_beta": mean_b,
                "colabfold_plddt": cf_plddt,
                "colabfold_ptm": cf_m["cf_ptm"] if cf_m else np.nan,
                "colabfold_pdb_target": cf_p["cf_pdb_target"] if cf_p else "missing",
                "colabfold_pdb_qtmscore": cf_qtm,
                "colabfold_afdb_target": cf_a["cf_afdb_target"] if cf_a else "missing",
                "colabfold_afdb_evalue": cf_a["cf_afdb_evalue"] if cf_a else np.nan,
                "esmfold_plddt": esm_plddt,
                "hhsearch_pfam_prob": hh_prob,
                "hhsearch_best_pfam": hh_pfam,
                "hhsearch_pdb_desc": hh_desc,
                "interpro_anymember_bits": ip_bits,
                "tcdb_hit": tc_str,
                "merops_hit": mer_str,
                "cazy_hit": cazy_str,
                "rep_locus_gene": str(rep_row["gene"]),
                "rep_locus_n_flank_pm5": int(n_flk),
                "rep_locus_n_known_flank_pm5": int(n_known_flk),
                "rep_locus_flank_summary": flk_desc,
                "independent_selection_score": indep_score,
            }
        )

    fam_a2_df = pd.DataFrame(fam_rows_a2).sort_values(
        ["independent_selection_score", "n_independent_genus_blocks", "n_distinct_orders", "colabfold_plddt"],
        ascending=[False, False, False, False],
    ).reset_index(drop=True)
    fam_a2_df.to_csv(P_RES / "task_a_route2_exported_dark_family_ranking.tsv", sep="\t", index=False)

    # Save annotation-negative exported controls matched by length and lineage
    annot_neg_ctrls = fam_a2_df[
        (fam_a2_df["hhsearch_pfam_prob"].isna() | (fam_a2_df["hhsearch_pfam_prob"] < 50.0))
        & (fam_a2_df["interpro_anymember_bits"] == 0)
        & (fam_a2_df["tcdb_hit"] == "missing")
        & (fam_a2_df["merops_hit"] == "missing")
        & (fam_a2_df["cazy_hit"] == "missing")
        & (fam_a2_df["colabfold_pdb_qtmscore"].isna() | (fam_a2_df["colabfold_pdb_qtmscore"] < 0.50))
        & (fam_a2_df["frac_candidate_promotion_eligible"] >= 0.80)
    ].head(25)
    annot_neg_ctrls.to_csv(P_RES / "task_a_route2_annotation_negative_exported_controls.tsv", sep="\t", index=False)

    # Compare Route A1 (gapfill) vs Route A2 (independent exported dark)
    a1_dark_fams = set(fam_a1_df[fam_a1_df["inherited_catalogue_dark"]]["family_id"])
    a2_fams = set(fam_a2_df["family_id"])
    overlap_fams = sorted(a1_dark_fams & a2_fams)
    overlap_summary = pd.DataFrame(
        [
            {
                "n_route_a1_dark_families": len(a1_dark_fams),
                "n_route_a2_recurrent_exported_dark_families": len(a2_fams),
                "n_overlap_families": len(overlap_fams),
                "jaccard_overlap": len(overlap_fams) / float(len(a1_dark_fams | a2_fams)) if (a1_dark_fams | a2_fams) else 0.0,
                "overlapping_family_examples": ";".join(overlap_fams[:15]) if overlap_fams else "none",
                "interpretation": (
                    "Gap-filling proposals (Route A1) are selected toward cytoplasmic metabolic reactions in reference models, "
                    "whereas Route A2 screens all recurrent exported sequence-dark families across the 237 genomes; "
                    "low overlap confirms why independent exported-dark selection is essential."
                ),
            }
        ]
    )
    overlap_summary.to_csv(P_RES / "task_a_route1_vs_route2_overlap.tsv", sep="\t", index=False)

    # =========================================================================
    # PART B: Dark-vs-Known Export Association in Thermoplasmatota & Lineage Bootstrap
    # =========================================================================
    print("[Stage 1] Running Part B: Dark-vs-known export association with lineage/length/completeness covariates...", flush=True)
    labelled_df = genes_df[
        genes_df["has_sp_tm_prediction"]
        & genes_df["dark"].notna()
        & genes_df["dark"].astype(str).isin({"0.0", "1.0", "0", "1"})
    ].copy()
    labelled_df["is_dark"] = pd.to_numeric(labelled_df["dark"], errors="coerce") == 1.0
    labelled_df["is_dark_strict"] = pd.to_numeric(labelled_df["dark_strict"], errors="coerce") == 1.0
    labelled_df["len_bin"] = np.where(
        labelled_df["length"] < 100,
        "<100",
        np.where(labelled_df["length"] < 250, "100-249", ">=250"),
    )

    # 1. Unadjusted and stratified Mantel-Haenszel summary across cohorts
    cohorts = {
        "all_labelled_orfs": labelled_df,
        "complete_interior_orfs_only": labelled_df[
            (labelled_df["orf_status"] == "complete") & (labelled_df["contig_edge_status"] == "interior")
        ],
        "comp70_contam5_complete_interior": labelled_df[
            labelled_df["sensitivity_comp70_contam5_pass"]
            & (labelled_df["orf_status"] == "complete")
            & (labelled_df["contig_edge_status"] == "interior")
        ],
        "ge100aa_complete_interior": labelled_df[
            (labelled_df["length"] >= 100)
            & (labelled_df["orf_status"] == "complete")
            & (labelled_df["contig_edge_status"] == "interior")
        ],
    }

    mh_rows = []
    for c_name, sub in cohorts.items():
        for endpoint in ("sp_any", "tmbed_S", "sp_or_S", "sp_and_S", "tm_any", "sec_spi", "sec_spii", "tat", "pilin"):
            d_mask = sub["is_dark"].to_numpy(dtype=bool)
            y_mask = sub[endpoint].to_numpy(dtype=bool)
            n_d = int(d_mask.sum())
            n_k = int((~d_mask).sum())
            a = int((d_mask & y_mask).sum())
            b = n_d - a
            c = int(((~d_mask) & y_mask).sum())
            d = n_k - c
            pct_d = 100.0 * a / n_d if n_d else np.nan
            pct_k = 100.0 * c / n_k if n_k else np.nan
            diff_pp = pct_d - pct_k
            or_unadj = (a * d) / float(b * c) if (b > 0 and c > 0) else np.nan

            # Mantel-Haenszel OR stratified by (tax_order, len_bin, comp_contam_tier)
            num_mh = 0.0
            den_mh = 0.0
            for _, s_grp in sub.groupby(["tax_order", "len_bin", "comp_contam_tier"], sort=False):
                sd = s_grp["is_dark"].to_numpy(dtype=bool)
                sy = s_grp[endpoint].to_numpy(dtype=bool)
                n_tot = len(s_grp)
                if n_tot < 2:
                    continue
                sa = float((sd & sy).sum())
                sb = float((sd & (~sy)).sum())
                sc = float(((~sd) & sy).sum())
                sd_ = float(((~sd) & (~sy)).sum())
                num_mh += (sa * sd_) / n_tot
                den_mh += (sb * sc) / n_tot
            or_mh_order_len = num_mh / den_mh if den_mh > 0 else np.nan

            # Mantel-Haenszel OR stratified by (genus_singleton_block, len_bin)
            num_g = 0.0
            den_g = 0.0
            for _, s_grp in sub.groupby(["genus_singleton_block", "len_bin"], sort=False):
                sd = s_grp["is_dark"].to_numpy(dtype=bool)
                sy = s_grp[endpoint].to_numpy(dtype=bool)
                n_tot = len(s_grp)
                if n_tot < 2:
                    continue
                sa = float((sd & sy).sum())
                sb = float((sd & (~sy)).sum())
                sc = float(((~sd) & sy).sum())
                sd_ = float(((~sd) & (~sy)).sum())
                num_g += (sa * sd_) / n_tot
                den_g += (sb * sc) / n_tot
            or_mh_genus_len = num_g / den_g if den_g > 0 else np.nan

            mh_rows.append(
                {
                    "cohort": c_name,
                    "endpoint": endpoint,
                    "n_dark_genes": n_d,
                    "n_known_genes": n_k,
                    "dark_positive": a,
                    "known_positive": c,
                    "pct_dark": pct_d,
                    "pct_known": pct_k,
                    "diff_percentage_points": diff_pp,
                    "or_unadjusted": or_unadj,
                    "or_mh_order_x_lenbin_x_qual": or_mh_order_len,
                    "or_mh_genusblock_x_lenbin": or_mh_genus_len,
                }
            )

    pd.DataFrame(mh_rows).to_csv(P_RES / "task_b_export_unadjusted_and_stratified.tsv", sep="\t", index=False)

    # 2. Lineage-blocked bootstrap (Order blocks, Genus-singleton blocks, Genome blocks)
    # Pre-aggregate counts per genome x len_bin so each bootstrap resample over blocks runs in <1ms!
    print("[Stage 1] Running 1,000-replicate lineage-blocked bootstrap (order, genus-singleton, genome blocks)...", flush=True)
    boot_rows = []
    for c_name, sub in [
        ("all_labelled_orfs", cohorts["all_labelled_orfs"]),
        ("complete_interior_orfs_only", cohorts["complete_interior_orfs_only"]),
        ("comp70_contam5_complete_interior", cohorts["comp70_contam5_complete_interior"]),
    ]:
        # Pre-aggregate per genome for fast exact block bootstrap
        g_agg = (
            sub.groupby(["genome", "tax_order", "genus_singleton_block"], sort=False)
            .agg(
                n_dark=("is_dark", lambda s: int(s.sum())),
                n_known=("is_dark", lambda s: int((~s).sum())),
                sp_dark=("sp_any", lambda s: int((s & sub.loc[s.index, "is_dark"]).sum())),
                sp_known=("sp_any", lambda s: int((s & (~sub.loc[s.index, "is_dark"])).sum())),
                tms_dark=("tmbed_S", lambda s: int((s & sub.loc[s.index, "is_dark"]).sum())),
                tms_known=("tmbed_S", lambda s: int((s & (~sub.loc[s.index, "is_dark"])).sum())),
                spor_dark=("sp_or_S", lambda s: int((s & sub.loc[s.index, "is_dark"]).sum())),
                spor_known=("sp_or_S", lambda s: int((s & (~sub.loc[s.index, "is_dark"])).sum())),
                tm_dark=("tm_any", lambda s: int((s & sub.loc[s.index, "is_dark"]).sum())),
                tm_known=("tm_any", lambda s: int((s & (~sub.loc[s.index, "is_dark"])).sum())),
            )
            .reset_index()
        )

        nd_arr = g_agg["n_dark"].to_numpy(dtype="float64")
        nk_arr = g_agg["n_known"].to_numpy(dtype="float64")

        for ep_label, d_col, k_col in [
            ("sp_any", "sp_dark", "sp_known"),
            ("tmbed_S", "tms_dark", "tms_known"),
            ("sp_or_S", "spor_dark", "spor_known"),
            ("tm_any", "tm_dark", "tm_known"),
        ]:
            yd_arr = g_agg[d_col].to_numpy(dtype="float64")
            yk_arr = g_agg[k_col].to_numpy(dtype="float64")

            def _stat_diff_and_logor(idx: np.ndarray) -> np.ndarray:
                nd = nd_arr[idx].sum()
                nk = nk_arr[idx].sum()
                yd = yd_arr[idx].sum()
                yk = yk_arr[idx].sum()
                pd_ = 100.0 * yd / nd
                pk_ = 100.0 * yk / nk
                diff = pd_ - pk_
                log_or = math.log((yd * (nk - yk)) / ((nd - yd) * yk))
                return np.array([diff, log_or], dtype="float64")

            for block_scheme, block_col in [
                ("tax_order_blocked", "tax_order"),
                ("genus_singleton_blocked", "genus_singleton_block"),
                ("genome_blocked", "genome"),
            ]:
                b_res = block_bootstrap_stat(
                    g_agg[block_col].tolist(),
                    _stat_diff_and_logor,
                    n_boot=1000,
                    seed=20261108,
                )
                est = b_res["estimate"]
                se = b_res["se"]
                lo95 = b_res["ci95_low"]
                hi95 = b_res["ci95_high"]
                boot_rows.append(
                    {
                        "cohort": c_name,
                        "endpoint": ep_label,
                        "block_scheme": block_scheme,
                        "n_independent_blocks": b_res["n_blocks"],
                        "n_genomes": b_res["n_obs"],
                        "diff_pp_estimate": float(est[0]),
                        "diff_pp_se": float(se[0]),
                        "diff_pp_ci95_low": float(lo95[0]),
                        "diff_pp_ci95_high": float(hi95[0]),
                        "or_estimate": float(math.exp(est[1])),
                        "or_ci95_low": float(math.exp(lo95[1])),
                        "or_ci95_high": float(math.exp(hi95[1])),
                    }
                )

    pd.DataFrame(boot_rows).to_csv(P_RES / "task_b_export_lineage_blocked_bootstrap.tsv", sep="\t", index=False)

    # 3. Multivariable logistic regression with order fixed effects, length, completeness, contamination, and contig-edge covariates
    print("[Stage 1] Fitting multivariable logistic regression with order-clustered standard errors...", flush=True)
    glm_rows = []
    for c_name, sub in [
        ("all_labelled_orfs", cohorts["all_labelled_orfs"]),
        ("complete_interior_orfs_only", cohorts["complete_interior_orfs_only"]),
    ]:
        # Aggregate into sufficient-statistic cells: (genome, tax_order, genus_singleton_block, is_dark, len_bin, is_edge)
        sub_mod = sub.copy()
        sub_mod["is_edge_num"] = (sub_mod["contig_edge_status"] == "edge").astype(float)
        sub_mod["log_len"] = np.log(np.maximum(sub_mod["length"].to_numpy(dtype="float64"), 1.0))
        cell_df = (
            sub_mod.groupby(["genome", "tax_order", "genus_singleton_block", "is_dark", "len_bin", "is_edge_num"], sort=False)
            .agg(
                n_total=("is_dark", "size"),
                n_sp_any=("sp_any", "sum"),
                n_tmbed_S=("tmbed_S", "sum"),
                n_sp_or_S=("sp_or_S", "sum"),
                mean_log_len=("log_len", "mean"),
                comp=("checkm_completeness", "first"),
                contam=("contamination", "first"),
            )
            .reset_index()
        )
        cell_df["is_dark_num"] = cell_df["is_dark"].astype(float)
        order_dummies = pd.get_dummies(cell_df["tax_order"], prefix="ord", drop_first=True, dtype=float)

        for ep_col, ep_name in [("n_sp_any", "sp_any"), ("n_tmbed_S", "tmbed_S"), ("n_sp_or_S", "sp_or_S")]:
            x_cols = ["is_dark_num", "mean_log_len", "comp", "contam"]
            if c_name == "all_labelled_orfs":
                x_cols.append("is_edge_num")
            X_df = pd.concat(
                [pd.Series(1.0, index=cell_df.index, name="const"), cell_df[x_cols], order_dummies],
                axis=1,
            )
            col_names = list(X_df.columns)
            dark_idx = col_names.index("is_dark_num")
            res_cov = fit_clustered_binomial_glm(
                cell_df[ep_col].to_numpy(dtype=float),
                cell_df["n_total"].to_numpy(dtype=float),
                X_df.to_numpy(dtype=float),
                cell_df["genus_singleton_block"].tolist(),
            )
            beta = float(res_cov["params"][dark_idx])
            se = float(res_cov["bse"][dark_idx])
            pval = float(res_cov["pvalues"][dark_idx])
            glm_rows.append(
                {
                    "cohort": c_name,
                    "endpoint": ep_name,
                    "covariates": "+".join(x_cols) + "+order_FE",
                    "cluster_se_group": "genus_singleton_block",
                    "n_clusters": int(res_cov["n_clusters"]),
                    "log_or_dark": beta,
                    "se_log_or": se,
                    "adjusted_or_dark": math.exp(beta),
                    "adjusted_or_ci95_low": math.exp(beta - 1.96 * se),
                    "adjusted_or_ci95_high": math.exp(beta + 1.96 * se),
                    "p_value_clustered": pval,
                }
            )

    pd.DataFrame(glm_rows).to_csv(P_RES / "task_b_export_adjusted_models.tsv", sep="\t", index=False)

    # 4. Independent lineage counts & External comparator harmonisation
    print("[Stage 1] Computing independent lineage counts and external comparator harmonisation...", flush=True)
    comp_path = require_file(P_IN / "v2_thermo_relatives_comparators.tsv", "v2_thermo_relatives_comparators.tsv")
    comp_df = pd.read_csv(comp_path, sep="\t")
    validate_primary_key(comp_df, ["genome_id"], "v2_thermo_relatives_comparators")
    assert len(comp_df) == 301, f"Expected 301 external comparators, got {len(comp_df)}"

    # Identify the 54 MGnify duplicates of GTDB reps (41 same assembly + 13 skani >= 99.9%)
    skani_path = require_file(P_IN / "v2_thermo_relatives_skani.tsv.gz", "v2_thermo_relatives_skani.tsv.gz")
    skani_df = pd.read_csv(skani_path, sep="\t")
    dup_ids = identify_mgnify_gtdb_duplicates(comp_df, skani_df, ani_threshold=99.9)
    comp_df["is_duplicate_of_gtdb_rep"] = comp_df["genome_id"].astype(str).isin(dup_ids)
    comp_dedup_df = comp_df[~comp_df["is_duplicate_of_gtdb_rep"]].copy().reset_index(drop=True)
    assert len(comp_dedup_df) == 247, f"Expected 247 deduplicated external comparators, got {len(comp_dedup_df)}"

    comp_tax = [parse_gtdb_lineage(x) for x in comp_dedup_df["lineage"]]
    for rk, col_name in [("c", "tax_class"), ("o", "tax_order"), ("f", "tax_family"), ("g", "tax_genus")]:
        comp_dedup_df[col_name] = [d[rk] for d in comp_tax]
    comp_genus_blocks, _ = normalize_genus_groups(comp_dedup_df["genome_id"].tolist(), comp_dedup_df["tax_genus"].tolist())
    comp_dedup_df["genus_singleton_block"] = comp_genus_blocks

    lineage_counts_df = pd.DataFrame(
        [
            {
                "cohort": "EQ_Thermoplasmatota_237",
                "n_genomes_total": len(thermo_df),
                "n_genomes_deduplicated": len(thermo_df),
                "n_classes": int(thermo_df["tax_class"].nunique()),
                "n_orders": int(thermo_df["tax_order"].nunique()),
                "n_families": int(thermo_df["tax_family"].nunique()),
                "n_named_genera": int(thermo_df[~thermo_df["is_unclassified_genus"]]["tax_genus"].nunique()),
                "n_unclassified_genus_genomes": int(thermo_df["is_unclassified_genus"].sum()),
                "n_genus_singleton_blocks": int(thermo_df["genus_singleton_block"].nunique()),
                "n_genomes_with_external_aai_ge65": int(thermo_df["has_external_comparator_aai65"].sum()),
            },
            {
                "cohort": "External_Comparators_Deduplicated_247",
                "n_genomes_total": len(comp_df),
                "n_genomes_deduplicated": len(comp_dedup_df),
                "n_classes": int(comp_dedup_df["tax_class"].nunique()),
                "n_orders": int(comp_dedup_df["tax_order"].nunique()),
                "n_families": int(comp_dedup_df["tax_family"].nunique()),
                "n_named_genera": int(
                    comp_dedup_df[~comp_dedup_df["tax_genus"].isin(["", "g__"])]["tax_genus"].nunique()
                ),
                "n_unclassified_genus_genomes": int(comp_dedup_df["tax_genus"].isin(["", "g__"]).sum()),
                "n_genus_singleton_blocks": int(comp_dedup_df["genus_singleton_block"].nunique()),
                "n_genomes_with_external_aai_ge65": len(comp_dedup_df),
            },
        ]
    )
    lineage_counts_df.to_csv(P_RES / "task_b_independent_lineage_counts.tsv", sep="\t", index=False)

    # 5. Presence/absence matched detectability in 538-genome OrthoFinder orthogroups
    # Partition EQ Thermoplasmatota complete interior genes into:
    #   - shared_with_external_comparators (`og_n_comparator_genomes > 0`)
    #   - eq_restricted_orthogroup (`og_n_comparator_genomes == 0` and `orthogroup != 'unassigned'`)
    #   - unassigned_singleton (`orthogroup == 'unassigned'`)
    comp_int = cohorts["complete_interior_orfs_only"].copy()
    comp_int["detectability_stratum"] = np.where(
        comp_int["og_n_comparator_genomes"] > 0,
        "shared_with_external_comparators",
        np.where(comp_int["orthogroup"] != "unassigned", "eq_restricted_orthogroup", "unassigned_singleton"),
    )

    det_rows = []
    for strat, s_grp in comp_int.groupby("detectability_stratum", sort=False):
        for lb in ["all", "<100", "100-249", ">=250"]:
            sub_lb = s_grp if lb == "all" else s_grp[s_grp["len_bin"] == lb]
            if sub_lb.empty:
                continue
            d_m = sub_lb["is_dark"].to_numpy(dtype=bool)
            nd = int(d_m.sum())
            nk = int((~d_m).sum())
            sp_d = int((d_m & sub_lb["sp_or_S"].to_numpy(dtype=bool)).sum())
            sp_k = int(((~d_m) & sub_lb["sp_or_S"].to_numpy(dtype=bool)).sum())
            det_rows.append(
                {
                    "detectability_stratum": strat,
                    "length_bin": lb,
                    "n_genes_total": len(sub_lb),
                    "n_inherited_dark": nd,
                    "n_inherited_known": nk,
                    "pct_inherited_dark": 100.0 * nd / len(sub_lb),
                    "pct_sp_or_S_all": 100.0 * float(sub_lb["sp_or_S"].mean()),
                    "pct_sp_or_S_in_dark": 100.0 * sp_d / nd if nd > 0 else np.nan,
                    "pct_sp_or_S_in_known": 100.0 * sp_k / nk if nk > 0 else np.nan,
                    "note": (
                        "Direct SignalP/TMbed export prediction on EQ genes stratified by external orthogroup presence/absence; "
                        "direct SignalP/TMbed on external comparator proteomes is blocked on Shaheen by IBEX-only signalp/6.0h license."
                    ),
                }
            )
    pd.DataFrame(det_rows).to_csv(P_RES / "task_b_orthogroup_detectability_export.tsv", sep="\t", index=False)

    # Save external comparator harmonisation status & gapseq metabolic comparison
    gapseq_sum_path = require_file(P_ANNOT / "results/c11_gapseq/tables/c12_summary.tsv.gz", "c12_summary.tsv.gz")
    gs_df = pd.read_csv(gapseq_sum_path, sep="\t")
    gs_df["is_eq_desert"] = gs_df["genome"].isin(thermo_genomes_set)
    dedup_ext_set = set(comp_dedup_df["genome_id"])
    gs_dedup = gs_df[gs_df["is_eq_desert"] | gs_df["genome"].isin(dedup_ext_set)].copy()
    gs_dedup["cohort"] = np.where(gs_dedup["is_eq_desert"], "EQ_Thermoplasmatota_237", "External_Comparators_Deduplicated_247")
    gs_dedup.to_csv(P_RES / "task_b_external_comparator_gapseq_harmonised.tsv", sep="\t", index=False)

    summary_payload = {
        "stage": "stage1_audit_and_export",
        "status": "completed",
        "n_eq_thermoplasmatota_genomes": int(len(thermo_df)),
        "n_eq_thermoplasmatota_genes": int(len(genes_df)),
        "n_gapfill_proposals_total": int(len(cand_audit_df)),
        "n_gapfill_proposals_priority": int(len(prio_df)),
        "n_gapfill_proposals_priority_dark": int(prio_df["inherited_catalogue_dark"].sum()),
        "n_route_a1_families": int(len(fam_a1_df)),
        "n_route_a2_independent_exported_dark_families": int(len(fam_a2_df)),
        "n_route_a1_a2_overlap_families": int(len(overlap_fams)),
        "external_comparators_total": int(len(comp_df)),
        "external_comparators_deduplicated": int(len(comp_dedup_df)),
        "external_signalp_tmbed_status": "BLOCKED_ON_SHAHEEN_BY_IBEX_ONLY_SIGNALP6_LICENSE_AND_MISSING_TMBED_CONTAINER",
        "external_uniref50_direct_search_status": "LAUNCHED_IN_STAGE2_MMSEQS2_UNIREF50_2026_01_S7.5_E1E-3_MAX10",
    }
    (P_RES / "stage1_summary.json").write_text(json.dumps(summary_payload, indent=2))
    print("[Stage 1] Completed successfully!", json.dumps(summary_payload), flush=True)


if __name__ == "__main__":
    main()
