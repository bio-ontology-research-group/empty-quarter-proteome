#!/usr/bin/env python3
"""Stage 1b: Authoritative Refined-FASTA (`genome_manifest.tsv`) Contig Retention & CDS Translation Audit,
Strict Representative-Own vs Family-Member Evidence Separation, and Within-Genome/Within-Genus Export Sensitivity.

Implements all requirements from:
- `coordinated-revision-20261108/stage1-scientific-review.txt` (Items 1-5)
- `coordinated-revision-20261108/stage1b-final-review.txt` (Items 1-5)
- `coordinated-revision-20261108/structural-source-location.txt`
"""
from __future__ import annotations

import collections
import glob
import gzip
import json
import math
import os
import shutil
import sys
import tarfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common_science import (
    audit_orf_integrity,
    classify_candidate_role,
    fit_clustered_binomial_glm,
    normalize_genus_groups,
    parse_and_verify_refined_fasta_bytes,
    require_file,
    scale_plddt_explicit_source,
    validate_primary_key,
    validate_stage1_to_stage1b_schema,
    validated_merge,
    verify_locus_against_refined_contig,
)

P_ANNOT = Path("/project/k10018/work/eq-annot-20261030")
P_ARCH = Path("/scratch/project/k10018/eq-archaeal-function-20261106")
P_PROT_INT = Path("/scratch/project/k10018/eq-proteome-integration-20261106")
P_SHAHEEN_DEREP = Path("/project/k10018/empty-quarter/data/catalog/mags/dereplicated")
P_OUT_ROOT = Path(os.environ.get("EQ_SCIENCE_OUT", "/scratch/project/k10018/eq-science-20261108"))
P_IN = P_OUT_ROOT / "inputs"
P_RES = P_OUT_ROOT / "results"


def load_authoritative_237_refined_genomes(
    thermo_genomes: list[str],
) -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    """Load and verify all 237 authoritative refined Thermoplasmatota FASTAs against `genome_manifest.tsv`.

    Reads the authoritative refined dRep FASTAs from `inputs/thermo237_refined_mags.tar.gz`
    (staged from `/ibex/scratch/projects/c2014/rob/metagenomics/05_refine/_drep/out/dereplicated_genomes/<genome>.fa`),
    verifies every genome's exact SHA-256, contig count, and base count against `inputs/genome_manifest.tsv`,
    rejects any duplicate contig IDs, and also audits the Shaheen `/project/k10018/empty-quarter/data/catalog/mags/dereplicated`
    links to document where `_sub` unrefined per-sample bins differed.
    """
    manifest_path = require_file(P_IN / "genome_manifest.tsv", "authoritative genome_manifest.tsv")
    man_df = pd.read_csv(manifest_path, sep="\t")
    validate_primary_key(man_df, ["mag"], "genome_manifest.tsv")
    man_lookup = man_df.set_index("mag")

    tar_path = require_file(P_IN / "thermo237_refined_mags.tar.gz", "thermo237_refined_mags.tar.gz")
    target_set = set(thermo_genomes)
    genome_contigs: dict[str, dict[str, str]] = {}
    rows: list[dict[str, Any]] = []

    with tarfile.open(tar_path, "r:gz") as tf:
        for member in tf:
            if not member.isfile():
                continue
            bname = Path(member.name).name
            if not bname.endswith(".fa"):
                continue
            gid = bname[:-3]
            if gid not in target_set:
                continue
            fobj = tf.extractfile(member)
            if fobj is None:
                raise ValueError(f"Could not extract {member.name} from {tar_path}")
            raw_bytes = fobj.read()
            if gid not in man_lookup.index:
                raise ValueError(f"Genome {gid} missing from authoritative genome_manifest.tsv")
            m_row = man_lookup.loc[gid]
            ver = parse_and_verify_refined_fasta_bytes(
                gid,
                raw_bytes,
                expected_sha256=str(m_row["sha256"]),
                expected_contigs=int(m_row["contigs"]),
                expected_bases=int(m_row["bases"]),
            )
            if ver["provenance_status"] != "verified_exact_match_genome_manifest":
                raise ValueError(
                    f"Authoritative refined FASTA mismatch for {gid}: {ver['provenance_status']} "
                    f"(computed_sha256={ver['computed_sha256']}, expected={m_row['sha256']})"
                )

            # Also inspect Shaheen dereplicated directory file/symlink for provenance transparency
            sh_gz = P_SHAHEEN_DEREP / f"{gid}.fa.gz"
            sh_fa = P_SHAHEEN_DEREP / f"{gid}.fa"
            sh_n_contigs = -1
            sh_compression = "missing"
            if sh_gz.exists():
                sh_compression = "gzip_per_sample_bin_link"
                with gzip.open(sh_gz, "rt") as fh:
                    sh_n_contigs = sum(1 for ln in fh if ln.startswith(">"))
            elif sh_fa.exists():
                sh_compression = "uncompressed_fa"
                with open(sh_fa, "rt") as fh:
                    sh_n_contigs = sum(1 for ln in fh if ln.startswith(">"))

            c_dict = ver.pop("contigs_dict")
            genome_contigs[gid] = c_dict  # type: ignore[assignment]
            ver["refined_archive_member"] = member.name
            ver["shaheen_catalog_storage"] = sh_compression
            ver["shaheen_catalog_n_contigs"] = sh_n_contigs
            ver["shaheen_catalog_contigs_match_refined"] = bool(sh_n_contigs == int(ver["n_refined_contigs"]))
            rows.append(ver)

    if set(genome_contigs.keys()) != target_set:
        missing = sorted(target_set - set(genome_contigs.keys()))
        raise ValueError(f"Missing {len(missing)} Thermoplasmatota genomes in thermo237_refined_mags.tar.gz: {missing[:5]}")

    out_df = pd.DataFrame(rows).sort_values("genome").reset_index(drop=True)
    return out_df, genome_contigs


def load_archaeal_aa_sequences(target_genes: set[str]) -> dict[str, str]:
    """Load raw AA sequences from `all_archaea.faa`, stripping at most ONE terminal '*' and preserving any internal '*'."""
    fasta_path = require_file(P_ARCH / "proteins/all_archaea.faa", "all_archaea.faa")
    seqs: dict[str, str] = {}
    cur_id: str | None = None
    chunks: list[str] = []

    def _store(gid: str | None, parts: list[str]) -> None:
        if gid is None or gid not in target_genes:
            return
        raw = "".join(parts).strip()
        seqs[gid] = raw[:-1] if raw.endswith("*") else raw

    with open(fasta_path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                _store(cur_id, chunks)
                hdr = line[1:].strip().split()[0]
                cur_id = hdr.split("|", 1)[1] if "|" in hdr else hdr
                chunks = []
            else:
                if cur_id in target_genes:
                    chunks.append(line.strip())
        _store(cur_id, chunks)
    return seqs


def load_raw_per_protein_structural_and_domain_evidence() -> dict[str, Any]:
    """Load raw per-protein structural and profile evidence keyed by exact protein ID (`sample::contig_gene`).

    Enforces explicit `source_scale='0_to_1'` via `scale_plddt_explicit_source`, captures alignment
    coverage (`qcov`, `tcov`) and model length (`qlen`), marks AFDB qTM-score as unavailable in the
    C9 AFDB tables (`afdb_qtmscore_status='unavailable_in_c9_afdb_table'`), and never conflates
    family-member evidence with representative-own evidence.
    """
    # 1. ColabFold members (`c10_structclust/colabfold/members.tsv.gz`, pLDDT on 0..1 scale)
    cf_mem: dict[str, dict[str, Any]] = {}
    with gzip.open(require_file(P_ANNOT / "results/c10_structclust/colabfold/members.tsv.gz", "cf members"), "rt") as fh:
        for line in fh:
            arm, pid, length_s, plddt_s, ptm_s, sname = line.rstrip("\n").split("\t")
            norm_pid = pid.replace("__", "::", 1) if "::" not in pid and "__" in pid else pid
            cf_mem[norm_pid] = {
                "model_protein_id": norm_pid,
                "model_length": int(length_s) if length_s not in ("NA", "", "nan") else -1,
                "plddt_0to1": scale_plddt_explicit_source(float(plddt_s), "0_to_1"),
                "ptm_0to1": scale_plddt_explicit_source(float(ptm_s), "0_to_1"),
                "plddt_source_scale": "0_to_1",
                "struct_name": sname,
            }

    # 2. ESMFold members (`c10_structclust/esmfold/members.tsv.gz`, pLDDT on 0..1 scale)
    esm_mem: dict[str, dict[str, Any]] = {}
    with gzip.open(require_file(P_ANNOT / "results/c10_structclust/esmfold/members.tsv.gz", "esm members"), "rt") as fh:
        for line in fh:
            arm, pid, length_s, plddt_s, ptm_s, sname = line.rstrip("\n").split("\t")
            norm_pid = pid.replace("__", "::", 1) if "::" not in pid and "__" in pid else pid
            esm_mem[norm_pid] = {
                "model_protein_id": norm_pid,
                "esm_arm": arm,
                "model_length": int(length_s) if length_s not in ("NA", "", "nan") else -1,
                "plddt_0to1": scale_plddt_explicit_source(float(plddt_s), "0_to_1"),
                "ptm_0to1": scale_plddt_explicit_source(float(ptm_s), "0_to_1"),
                "plddt_source_scale": "0_to_1",
                "struct_name": sname,
            }

    # 3. Foldseek PDB hits (both ColabFold `c9cf_colabfold` and ESMFold `c9_foldseek/c9pdb_*`)
    pdb_best: dict[str, dict[str, Any]] = {}
    pdb_files = [
        str(P_ANNOT / "results/c9cf_colabfold/pdb_whole.tsv.gz"),
        str(P_ANNOT / "results/c9cf_colabfold/pdb_regions.tsv.gz"),
    ] + sorted(glob.glob(str(P_ANNOT / "results/c9_foldseek/c9pdb_*/chunk_*.pdb_*.tsv.gz")))

    for p_str in pdb_files:
        p = Path(p_str)
        if not p.exists():
            continue
        src_tag = "ColabFold_AF2" if "c9cf_colabfold" in p_str else "ESMFold"
        with gzip.open(p, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                row = dict(zip(hdr, line.rstrip("\n").split("\t")))
                raw_pid = row.get("protein_id", "")
                pid = raw_pid.replace("__", "::", 1) if "::" not in raw_pid and "__" in raw_pid else raw_pid
                try:
                    qtm = float(row["qtmscore"])
                    ttm = float(row.get("ttmscore", "nan"))
                    lddt = float(row.get("lddt", "nan"))
                    qcov = float(row.get("qcov", "nan"))
                    tcov = float(row.get("tcov", "nan"))
                    qlen = int(float(row["qlen"])) if row.get("qlen") else -1
                except ValueError:
                    continue
                prev = pdb_best.get(pid)
                if prev is None or qtm > prev["pdb_qtmscore"]:
                    pdb_best[pid] = {
                        "pdb_model_source": src_tag,
                        "pdb_model_protein_id": pid,
                        "pdb_target": row.get("target", "missing"),
                        "pdb_region": row.get("region", "whole"),
                        "pdb_qtmscore": qtm,
                        "pdb_ttmscore": ttm,
                        "pdb_lddt": lddt,
                        "pdb_qcov": qcov,
                        "pdb_tcov": tcov,
                        "pdb_qlen": qlen,
                    }

    # 4. Foldseek AFDB hits (note: `afdb_whole.tsv.gz` has evalue, bits, qstart, qend, qlen, tstart, tend, tlen; NO qtmscore)
    afdb_best: dict[str, dict[str, Any]] = {}
    afdb_files = [
        str(P_ANNOT / "results/c9cf_colabfold/afdb_whole.tsv.gz"),
    ] + sorted(glob.glob(str(P_ANNOT / "results/c9_foldseek/c9afdb_*/chunk_*.afdb_whole.tsv.gz")))

    for p_str in afdb_files:
        p = Path(p_str)
        if not p.exists():
            continue
        with gzip.open(p, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                row = dict(zip(hdr, line.rstrip("\n").split("\t")))
                raw_pid = row.get("protein_id", "")
                pid = raw_pid.replace("__", "::", 1) if "::" not in raw_pid and "__" in raw_pid else raw_pid
                try:
                    ev = float(row["evalue"])
                    qstart = float(row["qstart"])
                    qend = float(row["qend"])
                    qlen = float(row["qlen"])
                    tstart = float(row["tstart"])
                    tend = float(row["tend"])
                    tlen = float(row["tlen"])
                except (KeyError, ValueError):
                    continue
                qcov = (abs(qend - qstart) + 1.0) / qlen if qlen > 0 else np.nan
                tcov = (abs(tend - tstart) + 1.0) / tlen if tlen > 0 else np.nan
                prev = afdb_best.get(pid)
                if prev is None or ev < prev["afdb_evalue"]:
                    afdb_best[pid] = {
                        "afdb_model_protein_id": pid,
                        "afdb_target": row.get("target", "missing"),
                        "afdb_evalue": ev,
                        "afdb_qcov": float(qcov),
                        "afdb_tcov": float(tcov),
                        "afdb_qlen": int(qlen),
                        "afdb_qtmscore_status": "unavailable_in_c9_afdb_table",
                    }

    # 5. HHsearch per-query (`c4_per_query.tsv.gz`, keyed by exact searched protein ID)
    c4_hits: dict[str, dict[str, Any]] = {}
    with gzip.open(require_file(P_ANNOT / "results/_agent_analysis/c4_c5_c7/c4/c4_per_query.tsv.gz", "c4_per_query"), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            qid = row["query"].replace("__", "::", 1)
            c4_hits[qid] = {
                "hh_query_protein_id": qid,
                "hh_pfam_prob": float(row["pfam_prob"]) if row["pfam_prob"] else np.nan,
                "hh_pdb_prob": float(row["pdb_prob"]) if row["pdb_prob"] else np.nan,
                "hh_best_pfam": row["best_pfam"] or "missing",
                "hh_pfam_id": row["pfam_id"] or "missing",
                "hh_best_pdb": row["best_pdb"] or "missing",
                "hh_pdb_desc": row["pdb_desc"] or "missing",
            }

    # 6. InterPro (`c5_family_bits.tsv.gz`: `rep_bits` applies to `family_rep`; `anymember_bits` is family-aggregate)
    c5_hits: dict[str, dict[str, Any]] = {}
    with gzip.open(require_file(P_ANNOT / "results/_agent_analysis/c4_c5_c7/c5/c5_family_bits.tsv.gz", "c5_family_bits"), "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(hdr, line.rstrip("\n").split("\t")))
            frep = row["family_rep"]
            c5_hits[frep] = {
                "interpro_family_rep_id": frep,
                "interpro_rep_own_bits": int(row["rep_bits"]),
                "family_interpro_anymember_bits": int(row["anymember_bits"]),
                "family_interpro_n_members": int(row["n_members"]),
            }

    return {
        "cf_mem": cf_mem,
        "esm_mem": esm_mem,
        "pdb_best": pdb_best,
        "afdb_best": afdb_best,
        "c4_hits": c4_hits,
        "c5_hits": c5_hits,
    }


def len_bin_7(length_aa: pd.Series | np.ndarray) -> pd.Series:
    """7-category fine amino-acid length stratification."""
    s = pd.Series(length_aa, dtype="float64")
    return pd.Series(
        np.select(
            [s < 100, s < 150, s < 200, s < 250, s < 350, s < 500],
            ["<100", "100-149", "150-199", "200-249", "250-349", "350-499"],
            default=">=500",
        ),
        index=s.index if isinstance(length_aa, pd.Series) else None,
    )


def mantel_haenszel_with_rbg_ci(
    df: pd.DataFrame,
    strata_cols: list[str],
    exposure_col: str,
    outcome_col: str,
) -> dict[str, float | int]:
    """Compute Mantel-Haenszel pooled odds ratio and Robins-Breslow-Greenland 95% CI across strata."""
    R_sum = 0.0
    S_sum = 0.0
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0
    n_informative_strata = 0
    n_genes_in_informative = 0

    for _, grp in df.groupby(strata_cols, sort=False):
        n_i = float(len(grp))
        if n_i < 2.0:
            continue
        x = grp[exposure_col].to_numpy(dtype=bool)
        y = grp[outcome_col].to_numpy(dtype=bool)
        a = float((x & y).sum())
        b = float((x & (~y)).sum())
        c = float(((~x) & y).sum())
        d = float(((~x) & (~y)).sum())
        if (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
            continue
        n_informative_strata += 1
        n_genes_in_informative += int(n_i)
        r_i = (a * d) / n_i
        s_i = (b * c) / n_i
        p_i = (a + d) / n_i
        q_i = (b + c) / n_i
        R_sum += r_i
        S_sum += s_i
        term1 += p_i * r_i
        term2 += p_i * s_i + q_i * r_i
        term3 += q_i * s_i

    if R_sum <= 0.0 or S_sum <= 0.0:
        return {
            "n_informative_strata": n_informative_strata,
            "n_genes_informative": n_genes_in_informative,
            "or_mh": np.nan,
            "se_log_or_rbg": np.nan,
            "or_ci95_low": np.nan,
            "or_ci95_high": np.nan,
            "p_value": np.nan,
        }

    or_mh = R_sum / S_sum
    var_log_or = term1 / (2.0 * R_sum * R_sum) + term2 / (2.0 * R_sum * S_sum) + term3 / (2.0 * S_sum * S_sum)
    se_log = math.sqrt(max(var_log_or, 1e-16))
    log_or = math.log(or_mh)
    from scipy import stats as sp_stats
    z = log_or / se_log
    pval = float(2.0 * sp_stats.norm.sf(abs(z)))
    return {
        "n_informative_strata": int(n_informative_strata),
        "n_genes_informative": int(n_genes_in_informative),
        "or_mh": float(or_mh),
        "se_log_or_rbg": float(se_log),
        "or_ci95_low": float(math.exp(log_or - 1.96 * se_log)),
        "or_ci95_high": float(math.exp(log_or + 1.96 * se_log)),
        "p_value": pval,
    }


def _evaluate_rep_and_family_evidence(
    rep_gene: str,
    rep_cat_id: str,
    rep_aa_len: int,
    family_id: str,
    member_cat_ids: list[str],
    member_genes: list[str],
    raw_ev: dict[str, Any],
    cazy_gene: dict[str, str],
    mer_gene: dict[str, str],
    tc_gene: dict[str, str],
) -> dict[str, Any]:
    """Strictly separate representative-own (`rep_own_*`) from family-member (`family_member_*`) evidence."""
    cf_mem = raw_ev["cf_mem"]
    esm_mem = raw_ev["esm_mem"]
    pdb_best = raw_ev["pdb_best"]
    afdb_best = raw_ev["afdb_best"]
    c4_hits = raw_ev["c4_hits"]
    c5_hits = raw_ev["c5_hits"]

    # 1. Representative-own structural model (`rep_cat_id` or `rep_gene`)
    rep_cf = cf_mem.get(rep_cat_id) or cf_mem.get(rep_gene)
    rep_esm = esm_mem.get(rep_cat_id) or esm_mem.get(rep_gene)
    rep_pdb = pdb_best.get(rep_cat_id) or pdb_best.get(rep_gene)
    rep_afdb = afdb_best.get(rep_cat_id) or afdb_best.get(rep_gene)

    rep_cf_plddt = float(rep_cf["plddt_0to1"]) if rep_cf is not None else np.nan
    rep_esm_plddt = float(rep_esm["plddt_0to1"]) if rep_esm is not None else np.nan
    rep_best_plddt = (
        max(v for v in (rep_cf_plddt, rep_esm_plddt) if pd.notna(v))
        if (pd.notna(rep_cf_plddt) or pd.notna(rep_esm_plddt))
        else np.nan
    )
    rep_model_len = (
        int(rep_cf["model_length"])
        if (rep_cf is not None and int(rep_cf["model_length"]) > 0)
        else (int(rep_esm["model_length"]) if (rep_esm is not None and int(rep_esm["model_length"]) > 0) else -1)
    )
    if rep_model_len <= 0 and rep_pdb is not None and int(rep_pdb["pdb_qlen"]) > 0:
        rep_model_len = int(rep_pdb["pdb_qlen"])

    rep_model_len_matches = bool(rep_model_len == rep_aa_len) if rep_model_len > 0 else False
    rep_own_model_qualified = bool(
        (rep_cf is not None or rep_esm is not None)
        and rep_aa_len >= 100
        and pd.notna(rep_best_plddt)
        and rep_best_plddt >= 0.70
        and (rep_model_len <= 0 or rep_model_len_matches)
    )

    rep_pdb_qtm = float(rep_pdb["pdb_qtmscore"]) if rep_pdb is not None else np.nan
    rep_pdb_qcov = float(rep_pdb["pdb_qcov"]) if rep_pdb is not None else np.nan
    rep_pdb_qlen_matches = (
        bool(int(rep_pdb["pdb_qlen"]) == rep_aa_len)
        if (rep_pdb is not None and int(rep_pdb["pdb_qlen"]) > 0)
        else False
    )
    rep_has_qualified_pdb_support = bool(
        rep_own_model_qualified
        and rep_pdb is not None
        and rep_pdb_qlen_matches
        and pd.notna(rep_pdb_qtm)
        and rep_pdb_qtm >= 0.50
        and pd.notna(rep_pdb_qcov)
        and rep_pdb_qcov >= 0.50
    )

    # 2. Representative-own profile/domain evidence (strictly on `rep_cat_id` / `rep_gene`)
    rep_hh = c4_hits.get(rep_cat_id) or c4_hits.get(rep_gene)
    rep_hh_prob = float(rep_hh["hh_pfam_prob"]) if rep_hh is not None else np.nan
    rep_hh_pfam = str(rep_hh["hh_best_pfam"]) if rep_hh is not None else "not_searched_on_rep"
    rep_hh_pdb_desc = str(rep_hh["hh_pdb_desc"]) if rep_hh is not None else "not_searched_on_rep"

    c5_entry = c5_hits.get(family_id)
    rep_is_c5_family_rep = bool(c5_entry is not None and rep_cat_id == family_id)
    rep_own_interpro_bits = int(c5_entry["interpro_rep_own_bits"]) if rep_is_c5_family_rep and c5_entry is not None else -1

    rep_tc = tc_gene.get(rep_gene, "missing")
    rep_mer = mer_gene.get(rep_gene, "missing")
    rep_cazy = cazy_gene.get(rep_gene, "missing")

    rep_has_own_profile_domain = bool(
        (pd.notna(rep_hh_prob) and rep_hh_prob >= 80.0)
        or (rep_own_interpro_bits > 0)
        or (rep_tc != "missing")
        or (rep_mer != "missing")
        or (rep_cazy != "missing")
    )

    # 3. Family-member evidence (explicitly preserving the contributing member ID when != rep_cat_id)
    search_ids = [family_id] + [cid for cid in member_cat_ids if cid != family_id]
    fam_cf_id = next((cid for cid in search_ids if cid in cf_mem), "none")
    fam_esm_id = next((cid for cid in search_ids if cid in esm_mem), "none")
    fam_pdb_id = next((cid for cid in search_ids if cid in pdb_best), "none")
    fam_hh_id = next((cid for cid in search_ids if cid in c4_hits), "none")

    fam_pdb_entry = pdb_best.get(fam_pdb_id) if fam_pdb_id != "none" else None
    fam_hh_entry = c4_hits.get(fam_hh_id) if fam_hh_id != "none" else None
    fam_c6_genes = [
        g for g in member_genes if (g in tc_gene or g in mer_gene or g in cazy_gene)
    ]

    return {
        "rep_cat_id": rep_cat_id,
        "rep_own_colabfold_plddt_0to1": rep_cf_plddt,
        "rep_own_esmfold_plddt_0to1": rep_esm_plddt,
        "rep_own_plddt_source_scale": "0_to_1",
        "rep_own_model_length": rep_model_len,
        "rep_own_model_length_matches_cds": rep_model_len_matches,
        "rep_own_model_qualified_ge100_plddt070": rep_own_model_qualified,
        "rep_own_pdb_target": str(rep_pdb["pdb_target"]) if rep_pdb is not None else "missing",
        "rep_own_pdb_qtmscore": rep_pdb_qtm,
        "rep_own_pdb_qcov": rep_pdb_qcov,
        "rep_has_qualified_pdb_structural_support": rep_has_qualified_pdb_support,
        "rep_own_afdb_target": str(rep_afdb["afdb_target"]) if rep_afdb is not None else "missing",
        "rep_own_afdb_evalue": float(rep_afdb["afdb_evalue"]) if rep_afdb is not None else np.nan,
        "rep_own_afdb_qcov": float(rep_afdb["afdb_qcov"]) if rep_afdb is not None else np.nan,
        "rep_own_afdb_qtmscore_status": str(rep_afdb["afdb_qtmscore_status"]) if rep_afdb is not None else "no_afdb_hit",
        "rep_own_hhsearch_pfam": rep_hh_pfam,
        "rep_own_hhsearch_pfam_prob": rep_hh_prob,
        "rep_own_hhsearch_pdb_desc": rep_hh_pdb_desc,
        "rep_own_interpro_bits": rep_own_interpro_bits,
        "rep_own_tcdb_hit": rep_tc,
        "rep_own_merops_hit": rep_mer,
        "rep_own_cazy_hit": rep_cazy,
        "rep_has_own_profile_domain_evidence": rep_has_own_profile_domain,
        # Separate family-aggregate / other-member columns with explicit member IDs
        "family_struct_model_member_id": fam_cf_id if fam_cf_id != "none" else fam_esm_id,
        "family_pdb_match_member_id": fam_pdb_id,
        "family_pdb_match_target": str(fam_pdb_entry["pdb_target"]) if fam_pdb_entry is not None else "missing",
        "family_pdb_match_qtmscore": float(fam_pdb_entry["pdb_qtmscore"]) if fam_pdb_entry is not None else np.nan,
        "family_pdb_match_qcov": float(fam_pdb_entry["pdb_qcov"]) if fam_pdb_entry is not None else np.nan,
        "family_hhsearch_member_id": fam_hh_id,
        "family_hhsearch_pfam": str(fam_hh_entry["hh_best_pfam"]) if fam_hh_entry is not None else "missing",
        "family_hhsearch_pfam_prob": float(fam_hh_entry["hh_pfam_prob"]) if fam_hh_entry is not None else np.nan,
        "family_interpro_rep_member_id": str(c5_entry["interpro_family_rep_id"]) if c5_entry is not None else "missing",
        "family_interpro_anymember_bits": int(c5_entry["family_interpro_anymember_bits"]) if c5_entry is not None else 0,
        "family_c6_contributing_genes": ";".join(sorted(fam_c6_genes)[:10]) if fam_c6_genes else "none",
    }


def main() -> None:
    P_RES.mkdir(parents=True, exist_ok=True)
    print("[Stage 1b] Loading and verifying 237 authoritative refined MAG FASTAs against genome_manifest.tsv...", flush=True)

    # Copy all-family structural summaries from `c17_struct_summary` per `structural-source-location.txt`
    for fn in ("esmfold_darkonly_summary.txt", "esmfold_summary.txt", "colabfold_summary.txt"):
        src_p = P_ANNOT / "results/c17_struct_summary" / fn
        if src_p.exists():
            shutil.copy2(src_p, P_RES / f"c17_{fn}")

    thermo_df = pd.read_csv(require_file(P_ARCH / "inputs/thermoplasmatota.tsv", "thermoplasmatota.tsv"), sep="\t")
    from stage1_audit_and_export import load_all_sp_tm, load_c6_thermo_gene_hits, parse_gtdb_lineage

    tax_parsed = [parse_gtdb_lineage(x) for x in thermo_df["gtdb_classification"]]
    for rk, col_name in [("c", "tax_class"), ("o", "tax_order"), ("f", "tax_family"), ("g", "tax_genus")]:
        thermo_df[col_name] = [d[rk] for d in tax_parsed]
    g_blocks, is_unclass = normalize_genus_groups(thermo_df["genome"].tolist(), thermo_df["tax_genus"].tolist())
    thermo_df["genus_singleton_block"] = g_blocks
    thermo_df["is_unclassified_genus"] = is_unclass

    thermo_genomes = sorted(thermo_df["genome"].tolist())
    fasta_manifest_df, genome_contigs = load_authoritative_237_refined_genomes(thermo_genomes)
    fasta_manifest_df.to_csv(P_RES / "task_a_stage1b_refined_genome_fasta_hashes.tsv", sep="\t", index=False)
    n_exact_man = int((fasta_manifest_df["provenance_status"] == "verified_exact_match_genome_manifest").sum())
    n_sh_diff = int((~fasta_manifest_df["shaheen_catalog_contigs_match_refined"]).sum())
    print(
        f"[Stage 1b] Verified {n_exact_man}/237 refined MAG FASTAs against genome_manifest.tsv "
        f"({n_sh_diff} _sub genomes on Shaheen had unrefined per-sample bin links).",
        flush=True,
    )

    # Load Stage 1 provisional tables and validate producer/consumer schema BEFORE the 6-min locus audit
    a1_props = pd.read_csv(require_file(P_RES / "task_a_route1_gapfill_proposals_audit.tsv.gz", "Route A1 proposals"), sep="\t")
    a2_fams = pd.read_csv(require_file(P_RES / "task_a_route2_exported_dark_family_ranking.tsv", "Route A2 ranking"), sep="\t")
    validate_stage1_to_stage1b_schema(a1_props, a2_fams)

    # Rename syntax/inherited columns explicitly
    a1_props = a1_props.rename(
        columns={
            "coords_valid": "coordinate_syntax_valid",
            "contig_in_catalogue": "inherited_contig_in_catalogue_flag",
            "candidate_promotion_eligible": "provisional_syntax_promotion_flag",
        }
    )

    # Load full Thermoplasmatota gene table & family mapping so we can verify EVERY member locus of all Route A1 & Route A2 families
    print("[Stage 1b] Loading Thermoplasmatota gene table, raw AA sequences, and raw per-protein evidence...", flush=True)
    arch_genes_df = pd.read_csv(require_file(P_ARCH / "inputs/archaeal_genes.tsv", "archaeal_genes.tsv"), sep="\t", low_memory=False)
    genes_df = arch_genes_df[arch_genes_df["genome"].isin(set(thermo_genomes))].copy().reset_index(drop=True)
    genes_df["start_num"] = pd.to_numeric(genes_df["start"], errors="coerce")
    genes_df["end_num"] = pd.to_numeric(genes_df["end"], errors="coerce")
    genes_df = genes_df.sort_values(["genome", "contig", "start_num", "end_num", "gene"]).reset_index(drop=True)
    genes_df["contig_rank"] = genes_df.groupby(["genome", "contig"]).cumcount() + 1
    genes_df["contig_n_genes"] = genes_df.groupby(["genome", "contig"])["gene"].transform("count")
    genes_df["has_full_seq"] = True
    genes_df["coords_valid"] = True
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

    mem2fam: dict[str, str] = {}
    with open(require_file(P_ANNOT / "inputs/fam/member_fam_final.tsv", "member_fam_final.tsv")) as fh:
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

    genes_df["family_id"] = [
        _resolve_fam(str(c), str(s), str(cr))
        for c, s, cr in zip(genes_df["cat_id"], genes_df["seq_cluster_rep"], genes_df["cat_cluster_rep"])
    ]
    genes_df = validated_merge(
        genes_df,
        thermo_df[
            [
                "genome",
                "tax_order",
                "tax_genus",
                "genus_singleton_block",
                "is_unclassified_genus",
                "checkm_completeness",
                "contamination",
            ]
        ],
        on="genome",
        how="left",
        validate="many_to_one",
        require_all_left=True,
        left_name="genes_df",
        right_name="thermo_df",
    )

    # Identify all candidate families across Route A1 and Route A2 and verify ALL their member loci against refined FASTAs
    cand_fams_all = set(a1_props["family_id"].astype(str)) | set(a2_fams["family_id"].astype(str))
    cand_members_df = genes_df[genes_df["family_id"].isin(cand_fams_all)].copy().reset_index(drop=True)
    target_aa_genes = set(cand_members_df["gene"].astype(str)) | set(a1_props["candidate_gene"].astype(str))
    aa_seqs = load_archaeal_aa_sequences(target_aa_genes)

    print(
        f"[Stage 1b] Verifying all {len(cand_members_df)} member loci across {len(cand_fams_all)} Route A1+A2 families "
        f"against authoritative refined genome FASTAs...",
        flush=True,
    )
    mem_v_cols: dict[str, list[Any]] = collections.defaultdict(list)
    for row in cand_members_df.itertuples(index=False):
        vres = verify_locus_against_refined_contig(
            genome_contigs.get(str(row.genome), {}),
            str(row.contig),
            row.start,
            row.end,
            str(row.strand),
            str(row.partial),
            aa_seqs.get(str(row.gene)),
        )
        for k, v in vres.items():
            mem_v_cols[k].append(v)
    for k, vals in mem_v_cols.items():
        cand_members_df[k] = vals

    cand_members_df["final_promotion_gate_pass"] = (
        cand_members_df["final_promotion_gate_pass"]
        & cand_members_df["contig_edge_status"].isin(["interior", "edge"])
    )
    cand_members_df["verified_complete_interior_member"] = (
        cand_members_df["final_promotion_gate_pass"]
        & (cand_members_df["contig_edge_status"] == "interior")
    )

    member_audit_cols = [
        "family_id",
        "gene",
        "cat_id",
        "genome",
        "contig",
        "start",
        "end",
        "strand",
        "partial",
        "length",
        "contig_edge_status",
        "tax_order",
        "tax_genus",
        "genus_singleton_block",
        "checkm_completeness",
        "contamination",
        "refined_contig_present",
        "refined_contig_bp_len",
        "locus_within_contig_bounds",
        "has_admissible_start_codon",
        "cds_translation_exact_match",
        "n_internal_stops_translated",
        "n_internal_stops_in_candidate",
        "has_terminal_stop_codon",
        "verification_status",
        "final_promotion_gate_pass",
        "verified_complete_interior_member",
    ]
    cand_members_df[member_audit_cols].to_csv(
        P_RES / "task_a_stage1b_all_candidate_member_loci_verified.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    # Index member verification by gene and by family_id
    gene_ver_lookup = cand_members_df.set_index("gene")
    fam_members_grouped = {fid: grp for fid, grp in cand_members_df.groupby("family_id", sort=False)}

    # Also attach row-level verification to all 7,381 Route A1 proposal rows
    a1_v_cols: dict[str, list[Any]] = collections.defaultdict(list)
    for row in a1_props.itertuples(index=False):
        vres = verify_locus_against_refined_contig(
            genome_contigs.get(str(row.genome), {}),
            str(row.contig),
            row.start,
            row.end,
            str(row.strand),
            str(row.partial_flag),
            aa_seqs.get(str(row.candidate_gene)),
        )
        for k, v in vres.items():
            a1_v_cols[k].append(v)
    for k, vals in a1_v_cols.items():
        a1_props[k] = vals
    a1_props["final_promotion_gate_pass"] = (
        a1_props["final_promotion_gate_pass"]
        & a1_props["contig_edge_status"].isin(["interior", "edge"])
        & a1_props["has_full_seq"]
    )
    a1_props.to_csv(P_RES / "task_a_stage1b_route1_gapfill_proposals_verified.tsv.gz", sep="\t", index=False, compression="gzip")

    # Load raw per-protein structural and domain evidence
    raw_ev = load_raw_per_protein_structural_and_domain_evidence()
    cazy_gene, mer_gene, tc_gene = load_c6_thermo_gene_hits(set(thermo_genomes), n_workers=16)

    # =========================================================================
    # Build Verified Route A1 Family Ranking (`task_a_stage1b_route1_gapfill_family_ranking_verified.tsv`)
    # =========================================================================
    prio_df = a1_props[a1_props["review_priority"]].copy()
    a1_fam_rows = []
    for fid, grp in prio_df.groupby("family_id", sort=False):
        # Prefer representative gene that passes final_promotion_gate_pass, is contig-interior, and has own model/profile
        grp_sorted = grp.assign(
            _is_int=(grp["contig_edge_status"] == "interior").astype(int),
            _gate=grp["final_promotion_gate_pass"].astype(int),
            _is_frep=(grp["candidate_gene"].map(lambda g: str(gene_ver_lookup.loc[g, "cat_id"]) if g in gene_ver_lookup.index else "") == fid).astype(int),
        ).sort_values(
            ["_gate", "_is_int", "_is_frep", "checkm_completeness", "context_score"],
            ascending=[False, False, False, False, False],
        )
        rep = grp_sorted.iloc[0]
        rep_gene = str(rep["candidate_gene"])
        rep_cat_id = str(gene_ver_lookup.loc[rep_gene, "cat_id"]) if rep_gene in gene_ver_lookup.index else rep_gene
        rep_len = int(rep["coord_aa_length"])

        fam_mem_df = fam_members_grouped.get(fid, pd.DataFrame())
        mem_cat_ids = fam_mem_df["cat_id"].astype(str).tolist() if not fam_mem_df.empty else [rep_cat_id]
        mem_genes = fam_mem_df["gene"].astype(str).tolist() if not fam_mem_df.empty else [rep_gene]

        ev_dict = _evaluate_rep_and_family_evidence(
            rep_gene,
            rep_cat_id,
            rep_len,
            fid,
            mem_cat_ids,
            mem_genes,
            raw_ev,
            cazy_gene,
            mer_gene,
            tc_gene,
        )

        # Verified member counts across all supporting member loci of this family
        ver_int_mem = fam_mem_df[fam_mem_df["verified_complete_interior_member"]] if not fam_mem_df.empty else pd.DataFrame()
        n_ver_int_loci = int(len(ver_int_mem))
        n_ver_int_genomes = int(ver_int_mem["genome"].nunique()) if not ver_int_mem.empty else 0
        n_ver_genus_blocks = int(ver_int_mem["genus_singleton_block"].nunique()) if not ver_int_mem.empty else 0

        rep_sub_bits = float(rep["sequence_bitscore"]) if pd.notna(rep["sequence_bitscore"]) else np.nan
        rep_own_profile = bool(ev_dict["rep_has_own_profile_domain_evidence"] or pd.notna(rep_sub_bits))
        rep_qual_pdb = bool(ev_dict["rep_has_qualified_pdb_structural_support"])

        row_out: dict[str, Any] = {
            "selection_route": "Route_A1_gapfill_selected",
            "family_id": fid,
            "candidate_role": classify_candidate_role(fid),
            "inherited_catalogue_dark": bool(grp["inherited_catalogue_dark"].any()),
            "actual_representative_gene": rep_gene,
            "actual_representative_cat_id": rep_cat_id,
            "actual_representative_genome": str(rep["genome"]),
            "rep_coord_aa_length": rep_len,
            "rep_verification_status": str(rep["verification_status"]),
            "rep_has_admissible_start_codon": bool(rep["has_admissible_start_codon"]),
            "rep_has_terminal_stop_codon": bool(rep["has_terminal_stop_codon"]),
            "rep_final_promotion_gate_pass": bool(rep["final_promotion_gate_pass"]),
            "rep_contig_edge_status": str(rep["contig_edge_status"]),
            "n_proposal_rows": int(len(grp)),
            "n_distinct_proposal_genes": int(grp["candidate_gene"].nunique()),
            "n_total_family_member_loci": int(len(fam_mem_df)),
            "n_verified_complete_interior_member_loci": n_ver_int_loci,
            "n_verified_complete_interior_genomes": n_ver_int_genomes,
            "n_verified_independent_genus_blocks": n_ver_genus_blocks,
            "rep_own_signalp_class": str(rep["signalp_class"]),
            "rep_own_sp_any": bool(rep["sp_any"]),
            "rep_own_tmbed_S": bool(rep["tmbed_S"]),
            "rep_own_tm_helices": int(rep["tm_helices"]),
            "rep_own_tm_beta": int(rep["tm_beta"]),
            "rep_own_subthreshold_bitscore": rep_sub_bits,
        }
        row_out.update(ev_dict)
        row_out["rep_has_own_profile_domain_evidence"] = rep_own_profile
        row_out["flanking_context_proposed_reactions"] = ";".join(sorted(grp["missing_reaction"].unique())[:5])
        row_out["flanking_context_proposed_pathways"] = ";".join(sorted(grp["pathway"].unique())[:5])
        row_out["rep_flanking_context_max_score"] = float(rep["context_score"])
        row_out["rep_flanking_genes_pm5_summary"] = str(rep["flank_annotations_pm5"])
        row_out["final_candidate_promotion_eligible"] = bool(
            rep["final_promotion_gate_pass"]
            and (rep["contig_edge_status"] == "interior")
            and n_ver_genus_blocks >= 2
            and (rep_qual_pdb or rep_own_profile)
        )
        a1_fam_rows.append(row_out)

    a1_fam_verified_df = pd.DataFrame(a1_fam_rows).sort_values(
        [
            "inherited_catalogue_dark",
            "final_candidate_promotion_eligible",
            "rep_has_qualified_pdb_structural_support",
            "rep_has_own_profile_domain_evidence",
            "n_verified_independent_genus_blocks",
        ],
        ascending=[False, False, False, False, False],
    ).reset_index(drop=True)
    a1_fam_verified_df.to_csv(P_RES / "task_a_stage1b_route1_gapfill_family_ranking_verified.tsv", sep="\t", index=False)

    # =========================================================================
    # Build Verified Route A2 Family Ranking (`task_a_stage1b_route2_exported_dark_family_ranking_verified.tsv`)
    # =========================================================================
    print("[Stage 1b] Building verified Route A2 independent exported-dark family ranking across all member loci...", flush=True)
    a2_verified_rows = []
    for row in a2_fams.itertuples(index=False):
        fid = str(row.family_id)
        fam_mem_df = fam_members_grouped.get(fid, pd.DataFrame())
        if fam_mem_df.empty:
            continue

        # Select concrete representative protein: prioritize verified_complete_interior_member, cat_id == fid, highest completeness
        fam_sorted = fam_mem_df.assign(
            _v_int=fam_mem_df["verified_complete_interior_member"].astype(int),
            _is_frep=(fam_mem_df["cat_id"].astype(str) == fid).astype(int),
        ).sort_values(
            ["_v_int", "_is_frep", "checkm_completeness", "length"],
            ascending=[False, False, False, False],
        )
        rep_m = fam_sorted.iloc[0]
        rep_gene = str(rep_m["gene"])
        rep_cat_id = str(rep_m["cat_id"])
        rep_len = int(rep_m["length"])

        mem_cat_ids = fam_mem_df["cat_id"].astype(str).tolist()
        mem_genes = fam_mem_df["gene"].astype(str).tolist()
        ev_dict = _evaluate_rep_and_family_evidence(
            rep_gene,
            rep_cat_id,
            rep_len,
            fid,
            mem_cat_ids,
            mem_genes,
            raw_ev,
            cazy_gene,
            mer_gene,
            tc_gene,
        )

        ver_int_mem = fam_mem_df[fam_mem_df["verified_complete_interior_member"]]
        n_ver_int_loci = int(len(ver_int_mem))
        n_ver_int_genomes = int(ver_int_mem["genome"].nunique())
        n_ver_genus_blocks = int(ver_int_mem["genus_singleton_block"].nunique())
        frac_ver_int = float(len(ver_int_mem) / len(fam_mem_df))

        # Own-protein domain classification (strictly from `rep_own_*`, never using flanking text or unverified other-member hits)
        own_txt = (
            f"{ev_dict['rep_own_hhsearch_pfam']} {ev_dict['rep_own_hhsearch_pdb_desc']} "
            f"{ev_dict['rep_own_tcdb_hit']} {ev_dict['rep_own_merops_hit']} {ev_dict['rep_own_cazy_hit']}"
        ).lower()
        if any(w in own_txt for w in ("pkd", "ig-like", "fn3", "fibronectin", "immunoglobulin")):
            s_layer_note = "broad_repeat_domain_in_own_protein_not_diagnostic_of_S_layer"
        elif "s-layer" in own_txt or "slh" in own_txt:
            s_layer_note = "putative_S_layer_homology_in_own_protein"
        else:
            s_layer_note = "no_S_layer_or_repeat_domain_in_own_protein"

        if ev_dict["rep_own_merops_hit"] != "missing":
            own_cat = "own_protein_peptidase_similarity_candidate"
        elif ev_dict["rep_own_cazy_hit"] != "missing":
            own_cat = "own_protein_CAZyme_similarity_candidate"
        elif ev_dict["rep_own_tcdb_hit"] != "missing":
            own_cat = "own_protein_transporter_similarity_candidate"
        elif ev_dict["rep_has_qualified_pdb_structural_support"] or ev_dict["rep_has_own_profile_domain_evidence"]:
            own_cat = "own_protein_divergent_fold_or_profile_match_candidate"
        else:
            own_cat = "exported_localization_only_no_own_domain_or_fold_match"

        d: dict[str, Any] = {
            "selection_route": "Route_A2_independent_exported_dark",
            "family_id": fid,
            "candidate_role": classify_candidate_role(fid),
            "actual_representative_gene": rep_gene,
            "actual_representative_cat_id": rep_cat_id,
            "actual_representative_genome": str(rep_m["genome"]),
            "rep_coord_aa_length": rep_len,
            "rep_verification_status": str(rep_m["verification_status"]),
            "rep_has_admissible_start_codon": bool(rep_m["has_admissible_start_codon"]),
            "rep_has_terminal_stop_codon": bool(rep_m["has_terminal_stop_codon"]),
            "rep_final_promotion_gate_pass": bool(rep_m["final_promotion_gate_pass"]),
            "rep_contig_edge_status": str(rep_m["contig_edge_status"]),
            "n_total_family_member_loci": int(len(fam_mem_df)),
            "n_verified_complete_interior_member_loci": n_ver_int_loci,
            "n_verified_complete_interior_genomes": n_ver_int_genomes,
            "n_verified_independent_genus_blocks": n_ver_genus_blocks,
            "frac_verified_complete_interior_loci": frac_ver_int,
            "frac_exported_sp_or_S": float(row.frac_sp_or_S),
            "frac_sp_any": float(row.frac_sp_any),
            "frac_tmbed_S": float(row.frac_tmbed_S),
        }
        d.update(ev_dict)
        d["own_protein_evidence_category"] = own_cat
        d["own_protein_s_layer_and_repeat_audit"] = s_layer_note
        d["flanking_context_hypothesis_only"] = str(row.rep_locus_flank_summary)
        d["final_candidate_promotion_eligible"] = bool(
            rep_m["verified_complete_interior_member"]
            and n_ver_int_genomes >= 3
            and n_ver_genus_blocks >= 2
            and frac_ver_int >= 0.70
        )
        a2_verified_rows.append(d)

    a2_verified_df = pd.DataFrame(a2_verified_rows).sort_values(
        [
            "final_candidate_promotion_eligible",
            "rep_has_qualified_pdb_structural_support",
            "rep_has_own_profile_domain_evidence",
            "n_verified_independent_genus_blocks",
            "rep_own_colabfold_plddt_0to1",
        ],
        ascending=[False, False, False, False, False],
    ).reset_index(drop=True)
    a2_verified_df.to_csv(P_RES / "task_a_stage1b_route2_exported_dark_family_ranking_verified.tsv", sep="\t", index=False)

    # Neutral overlap summary (Review Item 4)
    a1_set = set(a1_fam_verified_df[a1_fam_verified_df["inherited_catalogue_dark"]]["family_id"])
    a2_set = set(a2_verified_df["family_id"])
    ov = sorted(a1_set & a2_set)
    union_n = len(a1_set | a2_set)
    pd.DataFrame(
        [
            {
                "route_a1_definition": "Dark families among the 3,766 review-priority gap-filling reaction proposals",
                "route_a2_definition": "Dark families present in >=3 EQ Thermoplasmatota genomes with >=50% predicted export (sp_any or tmbed_S)",
                "n_route_a1_dark_families": len(a1_set),
                "n_route_a2_exported_dark_families": len(a2_set),
                "n_overlapping_families": len(ov),
                "jaccard_index": len(ov) / float(union_n) if union_n else 0.0,
                "overlapping_family_ids": ";".join(ov[:25]) if ov else "none",
            }
        ]
    ).to_csv(P_RES / "task_a_stage1b_route1_vs_route2_overlap.tsv", sep="\t", index=False)

    # =========================================================================
    # Review Item 3: Categorical-Length GLM & Within-Genome / Within-Named-Genus Mantel-Haenszel Sensitivity
    # =========================================================================
    print("[Stage 1b] Computing categorical-length GLM and within-genome / within-named-genus Mantel-Haenszel export associations...", flush=True)
    tm_all, sp_all = load_all_sp_tm(n_workers=32)
    sp_any_col = []
    tm_s_col = []
    has_pred = []
    for g in genes_df["gene"].tolist():
        kt = ("T", g)
        kd = ("D", g)
        s_e = sp_all.get(kt) or sp_all.get(kd)
        t_e = tm_all.get(kt) or tm_all.get(kd)
        if s_e is not None and t_e is not None:
            has_pred.append(True)
            sp_any_col.append(s_e[0] in {"SP", "LIPO", "TAT", "TATLIPO"})
            tm_s_col.append(bool(t_e[1]))
        else:
            has_pred.append(False)
            sp_any_col.append(False)
            tm_s_col.append(False)

    genes_df["has_sp_tm"] = has_pred
    genes_df["sp_any"] = sp_any_col
    genes_df["tmbed_S"] = tm_s_col
    genes_df["sp_or_S"] = genes_df["sp_any"] | genes_df["tmbed_S"]

    lab_df = genes_df[
        genes_df["has_sp_tm"]
        & genes_df["dark"].notna()
        & genes_df["dark"].astype(str).isin({"0.0", "1.0", "0", "1"})
    ].copy()
    lab_df["is_dark"] = pd.to_numeric(lab_df["dark"], errors="coerce") == 1.0
    lab_df["len_bin7"] = len_bin_7(lab_df["length"]).to_numpy()

    sens_rows = []
    for cohort_name, sub_df in [
        ("all_labelled_orfs", lab_df),
        (
            "complete_interior_orfs_only",
            lab_df[(lab_df["orf_status"] == "complete") & (lab_df["contig_edge_status"] == "interior")],
        ),
    ]:
        cell_cat = (
            sub_df.groupby(
                ["genome", "tax_order", "genus_singleton_block", "is_dark", "len_bin7", "contig_edge_status"],
                sort=False,
            )
            .agg(
                n_total=("is_dark", "size"),
                n_sp_any=("sp_any", "sum"),
                n_tmbed_S=("tmbed_S", "sum"),
                n_sp_or_S=("sp_or_S", "sum"),
                comp=("checkm_completeness", "first"),
                contam=("contamination", "first"),
            )
            .reset_index()
        )
        cell_cat["is_dark_num"] = cell_cat["is_dark"].astype(float)
        cell_cat["is_edge_num"] = (cell_cat["contig_edge_status"] == "edge").astype(float)
        len_dummies = pd.get_dummies(cell_cat["len_bin7"], prefix="len7", drop_first=True, dtype=float)
        ord_dummies = pd.get_dummies(cell_cat["tax_order"], prefix="ord", drop_first=True, dtype=float)

        for ep in ("sp_any", "tmbed_S", "sp_or_S"):
            mh_genome = mantel_haenszel_with_rbg_ci(sub_df, ["genome", "len_bin7"], "is_dark", ep)
            named_sub = sub_df[~sub_df["is_unclassified_genus"]]
            mh_named_genus = mantel_haenszel_with_rbg_ci(named_sub, ["tax_genus", "len_bin7"], "is_dark", ep)
            mh_genus_block = mantel_haenszel_with_rbg_ci(sub_df, ["genus_singleton_block", "len_bin7"], "is_dark", ep)

            base_cols = ["is_dark_num", "comp", "contam"]
            if cohort_name == "all_labelled_orfs":
                base_cols.append("is_edge_num")
            X_df = pd.concat(
                [pd.Series(1.0, index=cell_cat.index, name="const"), cell_cat[base_cols], len_dummies, ord_dummies],
                axis=1,
            )
            d_idx = list(X_df.columns).index("is_dark_num")
            glm_res = fit_clustered_binomial_glm(
                cell_cat[f"n_{ep}"].to_numpy(dtype=float),
                cell_cat["n_total"].to_numpy(dtype=float),
                X_df.to_numpy(dtype=float),
                cell_cat["genus_singleton_block"].tolist(),
            )
            b_glm = float(glm_res["params"][d_idx])
            se_glm = float(glm_res["bse"][d_idx])

            sens_rows.extend(
                [
                    {
                        "cohort": cohort_name,
                        "endpoint": ep,
                        "estimator": "within_genome_x_len7_Mantel_Haenszel",
                        "length_adjustment": "categorical_7bin_exact_strata",
                        "lineage_conditioning": "exact_within_genome_strata (237 genomes)",
                        "n_independent_lineage_units": int(sub_df["genome"].nunique()),
                        "n_informative_strata": mh_genome["n_informative_strata"],
                        "or_estimate": mh_genome["or_mh"],
                        "or_ci95_low": mh_genome["or_ci95_low"],
                        "or_ci95_high": mh_genome["or_ci95_high"],
                        "p_value": mh_genome["p_value"],
                    },
                    {
                        "cohort": cohort_name,
                        "endpoint": ep,
                        "estimator": "within_named_genus_x_len7_Mantel_Haenszel",
                        "length_adjustment": "categorical_7bin_exact_strata",
                        "lineage_conditioning": "exact_within_named_genus_strata (classified genera only)",
                        "n_independent_lineage_units": int(named_sub["tax_genus"].nunique()),
                        "n_informative_strata": mh_named_genus["n_informative_strata"],
                        "or_estimate": mh_named_genus["or_mh"],
                        "or_ci95_low": mh_named_genus["or_ci95_low"],
                        "or_ci95_high": mh_named_genus["or_ci95_high"],
                        "p_value": mh_named_genus["p_value"],
                    },
                    {
                        "cohort": cohort_name,
                        "endpoint": ep,
                        "estimator": "within_genus_singleton_block_x_len7_Mantel_Haenszel",
                        "length_adjustment": "categorical_7bin_exact_strata",
                        "lineage_conditioning": "exact_within_genus_or_unclassified_singleton_strata",
                        "n_independent_lineage_units": int(sub_df["genus_singleton_block"].nunique()),
                        "n_informative_strata": mh_genus_block["n_informative_strata"],
                        "or_estimate": mh_genus_block["or_mh"],
                        "or_ci95_low": mh_genus_block["or_ci95_low"],
                        "or_ci95_high": mh_genus_block["or_ci95_high"],
                        "p_value": mh_genus_block["p_value"],
                    },
                    {
                        "cohort": cohort_name,
                        "endpoint": ep,
                        "estimator": "order_FE_categorical_len7_Binomial_GLM_genus_clustered_SE",
                        "length_adjustment": "categorical_7bin_exact_cell_counts",
                        "lineage_conditioning": f"order_fixed_effects ({int(sub_df['tax_order'].nunique())} orders) + genus_singleton_clustered_SE",
                        "n_independent_lineage_units": int(glm_res["n_clusters"]),
                        "n_informative_strata": len(cell_cat),
                        "or_estimate": math.exp(b_glm),
                        "or_ci95_low": math.exp(b_glm - 1.96 * se_glm),
                        "or_ci95_high": math.exp(b_glm + 1.96 * se_glm),
                        "p_value": float(glm_res["pvalues"][d_idx]),
                    },
                ]
            )

    pd.DataFrame(sens_rows).to_csv(P_RES / "task_b_export_stage1b_within_genome_and_genus_sensitivity.tsv", sep="\t", index=False)

    # Update external comparator provenance & IBEX SignalP handoff table per external-comparison-scope-review.txt
    ext_handoff_df = pd.DataFrame(
        [
            {
                "comparison_dimension": "UniRef50_2026_01_direct_MMseqs2_search",
                "eq_237_source": "c7_thermo_evol/thermo_proteomes.tar.gz (from eq-annot-20261030/inputs/proteomes.faa)",
                "external_247_dedup_source": "c7_thermo_evol/thermo_proteomes.tar.gz (from eq-annot-20261030/inputs/thermo_comparators.tar)",
                "eq_gene_caller_and_mode": "Pyrodigal/Prodigal metagenomic mode (-p meta) from representative MAG catalogue",
                "external_gene_caller_and_mode": "Prodigal v2.6.3 (/project/k10018/sw/envs/mapping) single-genome mode (-p single, fallback -p meta) in scripts/c7b_orthofinder.sbatch",
                "database_and_search_harmonised": True,
                "gene_caller_mode_harmonised": False,
                "execution_status": "RUNNING_ON_SHAHEEN_STAGE2_ARRAY_15867464_AND_COLLECTOR_15867465",
                "scientific_scope_note": "Harmonises UniRef50 2026_01 database and MMseqs2 search rule (-s 7.5, E<=1e-3, max-accept 10, no extra id/cov gate), but does NOT harmonise gene-calling mode (-p meta vs -p single).",
            },
            {
                "comparison_dimension": "SignalP_6_0h_direct_export_prediction",
                "eq_237_source": "IBEX /ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/gpu_annot_20261030/out/signalp/ (40 chunks, 777,177 proteins)",
                "external_247_dedup_source": "IBEX /ibex/scratch/projects/c2014/rob/metagenomics/eq-science-20261109-external-signalp/ (60 chunks, 479,004 proteins across 247 deduplicated external comparators)",
                "eq_gene_caller_and_mode": "Pyrodigal/Prodigal metagenomic mode (-p meta)",
                "external_gene_caller_and_mode": "Prodigal v2.6.3 single-genome mode (-p single, fallback -p meta) from scripts/c7b_orthofinder.sbatch",
                "database_and_search_harmonised": True,
                "gene_caller_mode_harmonised": False,
                "execution_status": "RUNNING_ON_IBEX_BATCH_CPU_SIGNALP_6_0H_ORGANISM_OTHER_MODE_FAST",
                "scientific_scope_note": "Runs licensed SignalP 6.0h directly on IBEX CPU nodes with identical parameters (--organism other --mode fast) and joins with Stage 2 direct UniRef50 darkness via stage2b_join_external_signalp_and_direct_uniref50.py. Gene-caller mode difference (-p meta vs -p single) remains an explicit caveat unless external genomes are recalled in -p meta.",
            },
            {
                "comparison_dimension": "TMbed_transmembrane_and_signal_peptide_prediction",
                "eq_237_source": "IBEX/Shaheen gpu_annot_20261030/out/tmbed/ (40 chunks, 777,177 proteins)",
                "external_247_dedup_source": "not_run_for_external_comparators",
                "eq_gene_caller_and_mode": "Pyrodigal/Prodigal metagenomic mode (-p meta)",
                "external_gene_caller_and_mode": "Prodigal v2.6.3 single-genome mode (-p single, fallback -p meta)",
                "database_and_search_harmonised": False,
                "gene_caller_mode_harmonised": False,
                "execution_status": "UNAVAILABLE_EXTERNAL_TMBED_NOT_RUN",
                "scientific_scope_note": "External comparator export comparison uses direct SignalP 6.0h on IBEX; TMbed is used only within the 237 EQ MAG cohort.",
            },
        ]
    )
    ext_handoff_df.to_csv(P_RES / "task_b_stage1b_external_comparator_provenance_and_ibex_handoff.tsv", sep="\t", index=False)

    summary_1b = {
        "stage": "stage1b_refined_fasta_and_corrections",
        "status": "completed",
        "n_refined_genome_fastas_verified_against_manifest": n_exact_man,
        "n_shaheen_sub_genomes_with_unrefined_bin_links_documented": n_sh_diff,
        "n_candidate_member_loci_audited": int(len(cand_members_df)),
        "n_candidate_member_loci_verified_complete_interior": int(cand_members_df["verified_complete_interior_member"].sum()),
        "route_a1_priority_rows_total": int(len(prio_df)),
        "route_a1_priority_rows_verified_exact_complete_cds": int(prio_df["final_promotion_gate_pass"].sum()),
        "route_a1_dark_families_passing_final_promotion": int(
            a1_fam_verified_df[a1_fam_verified_df["inherited_catalogue_dark"]]["final_candidate_promotion_eligible"].sum()
        ),
        "route_a2_exported_dark_families_total": int(len(a2_verified_df)),
        "route_a2_exported_dark_families_passing_final_promotion": int(
            a2_verified_df["final_candidate_promotion_eligible"].sum()
        ),
    }
    (P_RES / "stage1b_summary.json").write_text(json.dumps(summary_1b, indent=2))
    print("[Stage 1b] Completed successfully!", json.dumps(summary_1b), flush=True)


if __name__ == "__main__":
    main()
