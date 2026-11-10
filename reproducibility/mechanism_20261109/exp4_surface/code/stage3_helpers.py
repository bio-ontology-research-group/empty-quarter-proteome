#!/usr/bin/env python3
"""Stage 3 specific scientific helpers and validation gates.

Implements all supervisor scientific review requirements from:
- `coordinated-revision-20261108/stage3-scientific-review.txt` (Items 1-6)
- `coordinated-revision-20261108/stage1b-final-review.txt` (Items 1-5 as applied to Stage 3)
- `coordinated-revision-20261108/scientific-design.txt`
- `coordinated-revision-20261108/compute-followup-design.txt`

Owned exclusively by the Stage 3 worker. Does not modify `common_science.py` or Stage 1/1b/2/4 files.
"""
from __future__ import annotations

import hashlib
import math
import subprocess
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from common_science import (
    PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES,
    TABLE_11_START_CODONS,
    audit_orf_integrity,
    clean_candidate_aa_sequence,
    parse_and_verify_refined_fasta_bytes,
    scale_plddt_explicit_source,
    translate_cds_table11,
    verify_locus_against_refined_contig,
    within_contig_circular_null,
)

NON_DOMAIN_INTERPRO_ANALYSES = frozenset({"MobiDBLite", "Coils", "AntiFam"})
BROAD_REPEAT_KEYWORDS = ("pkd", "ig-like", "fn3", "fibronectin", "immunoglobulin", "wd40", "kelch", "tpr", "ankyrin")
SLAYER_KEYWORDS = ("s-layer", "slh", "surface layer", "csaB")
PEPTIDASE_KEYWORDS = ("peptidase", "protease", "proteinase", "merops", "subtilisin", "metallopeptidase", "serine protease")
CAZY_KEYWORDS = ("glycosid", "glycosyl", "hydrolase", "cellulase", "chitinase", "xylanase", "esterase", "pectate", "cbm", "gh", "gt", "pl", "ce")
TRANSPORTER_KEYWORDS = ("transport", "permease", "symport", "antiport", "abc", "efflux", "importer", "exporter", "channel", "porin")
CHAPERONE_REPAIR_KEYWORDS = ("grpe", "dnaj", "dnak", "hsp", "chaperon", "thioredoxin", "peroxiredoxin", "ahpc", "rad50", "smc", "uvr", "muts", "rec")


def is_exact_integer_value(val: object) -> bool:
    """Return True iff `val` represents an exact finite integer (never truncating non-integer floats)."""
    if val is None or isinstance(val, bool) or pd.isna(val):
        return False
    if isinstance(val, (int, np.integer)):
        return True
    if isinstance(val, (float, np.floating)):
        v = float(val)
        return bool(np.isfinite(v) and v == math.floor(v))
    if isinstance(val, str):
        s = val.strip()
        if not s:
            return False
        if s.startswith(("+", "-")):
            s = s[1:]
        return s.isdigit()
    return False


def audit_gene_table_syntax_and_fasta(
    genes_df: pd.DataFrame,
    fasta_aa_map: Mapping[str, str],
    *,
    edge_gene_buffer: int = 2,
) -> pd.DataFrame:
    """Audit `genes_df` without ever setting `has_full_seq=True` or `coordinate_syntax_valid=True` hardcoded.

    - Checks exact integer coordinates (`start >= 1`, `end > start`, `(end - start + 1) % 3 == 0`) and valid strand (`1`, `-1`, `+`, `-`).
    - Checks actual sequence presence in `fasta_aa_map`, permitting removal of at most ONE terminal `*` and
      rejecting any internal `*` stop codon, length mismatch, or MD5 mismatch.
    - Records missing sequence as `'unknown_missing_aa_sequence'` (and `has_full_seq = False`), never as `True`.
    - Computes `contig_rank`, `contig_n_genes`, and runs `audit_orf_integrity`.
    - Note: `provisional_syntax_promotion_eligible` is syntax/FASTA-table level only; biological promotion
      strictly requires `member_final_promotion_gate_pass` against the authoritative refined genome FASTA.
    """
    req = {"genome", "gene", "contig", "start", "end", "strand", "partial", "length", "internal_stop", "md5", "contig_in_catalogue"}
    missing_cols = sorted(req - set(genes_df.columns))
    if missing_cols:
        raise ValueError(f"audit_gene_table_syntax_and_fasta: missing required columns {missing_cols}")

    out = genes_df.copy()
    out["start_num"] = pd.to_numeric(out["start"], errors="coerce")
    out["end_num"] = pd.to_numeric(out["end"], errors="coerce")
    out = out.sort_values(["genome", "contig", "start_num", "end_num", "gene"]).reset_index(drop=True)
    out["contig_rank"] = out.groupby(["genome", "contig"]).cumcount() + 1
    out["contig_n_genes"] = out.groupby(["genome", "contig"])["gene"].transform("count")

    coord_syntax_ok: list[bool] = []
    has_full_seq_list: list[bool] = []
    seq_audit_status_list: list[str] = []
    clean_aa_list: list[str] = []

    for row in out.itertuples(index=False):
        s_ok = is_exact_integer_value(row.start)
        e_ok = is_exact_integer_value(row.end)
        strd_ok = str(row.strand).strip() in {"1", "-1", "+1", "+", "-"}
        if s_ok and e_ok and strd_ok:
            s_i = int(float(row.start))
            e_i = int(float(row.end))
            nt_len = e_i - s_i + 1
            c_ok = (s_i >= 1) and (e_i > s_i) and ((nt_len % 3) == 0)
        else:
            c_ok = False
        coord_syntax_ok.append(bool(c_ok))

        gid = str(row.gene)
        raw_seq = fasta_aa_map.get(gid)
        cand_clean = clean_candidate_aa_sequence(raw_seq)
        if cand_clean["status"] == "missing_candidate_sequence":
            has_full_seq_list.append(False)
            seq_audit_status_list.append("unknown_missing_aa_sequence")
            clean_aa_list.append("")
        elif cand_clean["status"] != "ok":
            has_full_seq_list.append(False)
            seq_audit_status_list.append(str(cand_clean["status"]))
            clean_aa_list.append("")
        else:
            aa_body = str(cand_clean["aa_body"])
            md5_hex = hashlib.md5(aa_body.encode("ascii")).hexdigest()
            len_ok = is_exact_integer_value(row.length) and (len(aa_body) == int(float(row.length)))
            md5_ok = (str(row.md5).strip().lower() == md5_hex)
            if len_ok and md5_ok:
                has_full_seq_list.append(True)
                seq_audit_status_list.append("fasta_aa_matches_table_len_and_md5")
                clean_aa_list.append(aa_body)
            else:
                has_full_seq_list.append(False)
                seq_audit_status_list.append("fasta_aa_len_or_md5_mismatch")
                clean_aa_list.append(aa_body)

    out["coordinate_syntax_valid"] = coord_syntax_ok
    out["has_full_seq"] = has_full_seq_list
    out["fasta_sequence_audit_status"] = seq_audit_status_list
    out["audited_aa_sequence"] = clean_aa_list
    out["inherited_contig_in_catalogue_flag"] = out["contig_in_catalogue"]

    out = audit_orf_integrity(
        out,
        partial_col="partial",
        internal_stop_col="internal_stop",
        length_col="length",
        contig_rank_col="contig_rank",
        contig_n_genes_col="contig_n_genes",
        contig_in_cat_col="inherited_contig_in_catalogue_flag",
        has_full_seq_col="has_full_seq",
        coords_valid_col="coordinate_syntax_valid",
        edge_gene_buffer=edge_gene_buffer,
    )
    out = out.rename(columns={"candidate_promotion_eligible": "provisional_syntax_promotion_flag"})
    return out


def verify_member_locus_strict(
    genome_manifest_Status: str,
    contigs_dict: Mapping[str, str],
    contig_id: str,
    start: object,
    end: object,
    strand: object,
    partial_flag: object,
    candidate_aa_seq: object,
    contig_edge_status: str,
) -> dict[str, Any]:
    """Verify a single member locus against the authoritative refined genome FASTA and `genome_manifest.tsv`.

    A failed or unknown sequence audit or unverified genome FASTA manifest hash prevents biological promotion.
    """
    vres = verify_locus_against_refined_contig(
        dict(contigs_dict),
        str(contig_id),
        start,
        end,
        strand,
        str(partial_flag),
        candidate_aa_seq,
    )
    manifest_ok = (genome_manifest_Status == "verified_exact_match_genome_manifest")
    edge_known = contig_edge_status in ("interior", "edge")
    gate_pass = bool(
        manifest_ok
        and vres["final_promotion_gate_pass"]
        and edge_known
    )
    out = dict(vres)
    out["genome_manifest_provenance_status"] = genome_manifest_Status
    out["contig_edge_status"] = contig_edge_status
    out["member_final_promotion_gate_pass"] = gate_pass
    if not manifest_ok and out["verification_status"] == "verified_exact_complete_cds":
        out["verification_status"] = f"unverified_genome_manifest_{genome_manifest_Status}"
    elif not edge_known and out["verification_status"] == "verified_exact_complete_cds":
        out["verification_status"] = "unknown_contig_edge_status"
    return out


def classify_family_control_or_candidate_role(
    family_id: str,
    *,
    is_dark: bool,
    selection_route: str,
    is_matched_negative_exported_control: bool = False,
    claim_annotated_positive_control: bool = False,
    curated_reference_accession: str | None = None,
    curated_function_description: str | None = None,
    curated_evidence_source: str | None = None,
    query_domain_signatures: Sequence[str] | None = None,
    reference_domain_signatures: Sequence[str] | None = None,
    independently_verified_curated_homology: bool = False,
    query_to_curated_ref_cov: float | None = None,
    strict_reject_unverified_positive_control: bool = False,
) -> dict[str, str | bool]:
    """Classify a family's role while preventing non-dark families from becoming positive controls without curated evidence.

    Implements `stage3-scientific-review.txt` Item 2 & Item 6 and `stage3-final-integrity-review.txt` Item 2:
    - Prior heat-stress candidate examples (`F41Dr3::F41Dr3_256695_46`, `49Dr1::49Dr1_177214_2`, `24Dr1::24Dr1_1575644_3`)
      are strictly `'prior_heat_stress_candidate_example'`; claiming them as positive controls raises `ValueError`.
    - Forced MSA alignment coverage (`>= 0.70`) + a generic InterPro domain on the query + a reference sequence
      from a proposed reaction CANNOT establish an `'annotated_positive_control'`.
    - To be classified as `'annotated_positive_control'`, a family must have:
        * `is_dark == False`
        * `independently_verified_curated_homology == True`
        * non-empty `curated_reference_accession` and `curated_function_description`
        * `curated_evidence_source` in `{'SwissProt_Archaea_rev', 'SwissProt_curated', 'experimentally_characterised_archaeal_reference'}`
        * non-empty `query_domain_signatures` and `reference_domain_signatures` with exact matching domain architecture
          (`set(query_domain_signatures) == set(reference_domain_signatures)`)
        * finite `query_to_curated_ref_cov >= 0.70`
      Otherwise it is labelled `'sequence_annotated_comparator'` with
      `positive_control_evidence_status = 'positive_control_evidence_unavailable'` (or raises `ValueError` when
      `strict_reject_unverified_positive_control=True`).
    """
    fid = str(family_id).strip()
    if not fid:
        raise ValueError("Empty family_id in classify_family_control_or_candidate_role")

    if fid in PRIOR_HEAT_STRESS_CANDIDATE_EXAMPLES:
        if claim_annotated_positive_control:
            raise ValueError(
                f"{fid} is a prior heat-stress dark candidate example with unresolved genus-aggregated effects, "
                f"never an annotated positive control."
            )
        return {
            "family_id": fid,
            "candidate_role": "prior_heat_stress_candidate_example",
            "positive_control_evidence_status": "not_applicable_prior_dark_candidate_example",
            "has_curated_positive_control_evidence": False,
        }

    if is_matched_negative_exported_control:
        if claim_annotated_positive_control:
            raise ValueError(f"{fid}: matched annotation-negative exported control cannot be a positive control")
        return {
            "family_id": fid,
            "candidate_role": "matched_annotation_negative_exported_control",
            "positive_control_evidence_status": "not_applicable_negative_exported_control_arm",
            "has_curated_positive_control_evidence": False,
        }

    if is_dark:
        if claim_annotated_positive_control:
            raise ValueError(f"{fid}: sequence-dark family cannot be classified as an annotated positive control")
        return {
            "family_id": fid,
            "candidate_role": "de_novo_dark_candidate_family",
            "positive_control_evidence_status": "not_applicable_dark_candidate",
            "has_curated_positive_control_evidence": False,
        }

    acc = "" if curated_reference_accession is None else str(curated_reference_accession).strip()
    desc = "" if curated_function_description is None else str(curated_function_description).strip()
    src = "" if curated_evidence_source is None else str(curated_evidence_source).strip()
    q_doms = {str(x).strip() for x in (query_domain_signatures or []) if str(x).strip()}
    r_doms = {str(x).strip() for x in (reference_domain_signatures or []) if str(x).strip()}
    arch_match = bool(q_doms and r_doms and q_doms == r_doms)
    cov_ok = (
        query_to_curated_ref_cov is not None
        and not pd.isna(query_to_curated_ref_cov)
        and np.isfinite(float(query_to_curated_ref_cov))
        and float(query_to_curated_ref_cov) >= 0.70
    )
    valid_curated = (
        bool(independently_verified_curated_homology)
        and bool(acc)
        and acc not in {"missing", "none", "nan", "NA"}
        and bool(desc)
        and desc not in {"missing", "none", "nan", "NA", "hypothetical protein", "uncharacterized protein"}
        and src in {"SwissProt_Archaea_rev", "SwissProt_curated", "experimentally_characterised_archaeal_reference"}
        and arch_match
        and cov_ok
    )

    if claim_annotated_positive_control and not valid_curated:
        if strict_reject_unverified_positive_control:
            raise ValueError(
                f"Cannot promote non-dark family '{fid}' to 'annotated_positive_control' without independently verified "
                f"curated homology and matching reference domain architecture "
                f"(got acc={acc!r}, src={src!r}, ind_homology={independently_verified_curated_homology}, "
                f"arch_match={arch_match}, cov={query_to_curated_ref_cov})."
            )
        return {
            "family_id": fid,
            "candidate_role": "sequence_annotated_comparator",
            "positive_control_evidence_status": "positive_control_evidence_unavailable",
            "has_curated_positive_control_evidence": False,
        }

    if valid_curated and claim_annotated_positive_control:
        return {
            "family_id": fid,
            "candidate_role": "annotated_positive_control",
            "positive_control_evidence_status": f"verified_curated_reference:{acc}({src})",
            "has_curated_positive_control_evidence": True,
        }

    return {
        "family_id": fid,
        "candidate_role": "sequence_annotated_comparator",
        "positive_control_evidence_status": "positive_control_evidence_unavailable",
        "has_curated_positive_control_evidence": False,
    }


def qualify_structural_and_domain_evidence_strict(
    representative_gene: str,
    representative_aa_seq: str | None,
    member_evidence_records: Sequence[Mapping[str, Any]],
    *,
    allow_member_summary_separate: bool = True,
    forbid_mixed_locus_in_rep: bool = True,
) -> dict[str, Any]:
    """Bind structural and domain evidence to `representative_gene` without ever taking cross-locus maxima or fabricating sequences.

    Implements `stage3-scientific-review.txt` Items 3 & 6, `stage1b-final-review.txt` Items 3 & 4,
    and `stage3-final-integrity-review.txt` Item 1:
    - NEVER fabricates a fallback AA sequence if `representative_aa_seq` is missing or invalid; instead records
      `rep_sequence_audit_status = rep_clean['status']`, `rep_aa_length = 0`, and sets all `rep_own_*` qualified
      structural/domain flags to `False`.
    """
    rep_clean = clean_candidate_aa_sequence(representative_aa_seq)
    rep_seq_valid = (rep_clean["status"] == "ok")
    rep_seq = str(rep_clean["aa_body"]) if rep_seq_valid else ""
    rep_len = len(rep_seq) if rep_seq_valid else 0

    rep_plddt = np.nan
    rep_model_source = "none"
    rep_model_seq_prov_status = "no_structural_model" if rep_seq_valid else f"unavailable_invalid_rep_sequence_{rep_clean['status']}"
    rep_model_qualified = False
    rep_pdb_target = "missing"
    rep_pdb_qtm = np.nan
    rep_pdb_qcov = np.nan
    rep_qualified_pdb_match = False

    rep_afdb_target = "missing"
    rep_afdb_evalue = np.nan
    rep_afdb_qcov = np.nan
    rep_afdb_qtm = np.nan
    rep_qualified_afdb_match = False
    rep_afdb_status = "no_afdb_hit" if rep_seq_valid else f"unavailable_invalid_rep_sequence_{rep_clean['status']}"

    rep_hh_pfam = "missing"
    rep_hh_prob = np.nan
    rep_hh_qcov = np.nan
    rep_qualified_hh_profile = False

    rep_ipr_domains: list[str] = []
    rep_ipr_nondomain: list[str] = []

    other_qualified_struct_members: list[str] = []
    other_qualified_profile_members: list[str] = []

    for rec in member_evidence_records:
        loc_gene = str(rec.get("locus_gene", "")).strip()
        loc_seq_raw = rec.get("locus_aa_seq")
        claim_rep = bool(rec.get("claim_as_representative_own", False))

        c_info = clean_candidate_aa_sequence(loc_seq_raw)
        loc_seq_ok = (c_info["status"] == "ok")
        loc_seq = str(c_info["aa_body"]) if loc_seq_ok else ""

        is_exact_rep = (loc_gene == representative_gene) or (rep_seq_valid and loc_seq_ok and loc_seq == rep_seq)
        if claim_rep and not is_exact_rep and forbid_mixed_locus_in_rep:
            raise ValueError(
                f"Mixed-locus evidence aggregation rejected: record from locus '{loc_gene}' cannot be claimed "
                f"as representative-own evidence for '{representative_gene}'."
            )

        scale_str = str(rec.get("plddt_source_scale", "0_to_1"))
        plddt_01 = scale_plddt_explicit_source(rec.get("plddt_raw"), scale_str)
        target_audited_seq = loc_seq if loc_seq_ok else (rep_seq if (loc_gene == representative_gene and rep_seq_valid) else "")
        seq_len_for_gate = len(target_audited_seq) if target_audited_seq else -1
        has_known_valid_seq = bool(target_audited_seq)

        # Item 3 of stage3-final-integrity-review.txt:
        # Verify actual model polymer sequence against the audited locus AA sequence (length equality alone is NOT proof),
        # and verify model_source matches hit_source (never pair AF2 pLDDT with an ESMFold Foldseek hit or vice versa).
        mod_src = str(rec.get("model_source") or "").strip()
        mod_poly_raw = rec.get("model_polymer_aa_seq")
        mod_poly_info = clean_candidate_aa_sequence(mod_poly_raw)
        if not pd.notna(plddt_01) and not mod_src:
            mod_seq_status = "no_structural_model"
            mod_seq_matches_locus = False
        elif mod_poly_info["status"] != "ok":
            mod_seq_status = "unverified_missing_model_polymer_sequence"
            mod_seq_matches_locus = False
        elif not has_known_valid_seq:
            mod_seq_status = "unavailable_invalid_locus_sequence"
            mod_seq_matches_locus = False
        elif str(mod_poly_info["aa_body"]) != target_audited_seq:
            mod_seq_status = (
                "mismatch_model_polymer_sequence_same_length"
                if len(str(mod_poly_info["aa_body"])) == len(target_audited_seq)
                else "mismatch_model_polymer_sequence_and_length"
            )
            mod_seq_matches_locus = False
        else:
            mod_seq_status = "verified_exact_model_polymer_sequence_match"
            mod_seq_matches_locus = True

        model_qual = bool(
            has_known_valid_seq
            and mod_seq_matches_locus
            and bool(mod_src)
            and seq_len_for_gate >= 100
            and pd.notna(plddt_01)
            and float(plddt_01) >= 0.70
        )

        pdb_hit_src = str(rec.get("pdb_hit_source") or "").strip()
        pdb_src_agrees = bool(mod_src and pdb_hit_src and mod_src == pdb_hit_src)
        qtm_raw = rec.get("pdb_qtmscore")
        qcov_raw = rec.get("pdb_qcov")
        qtm_val = float(qtm_raw) if (qtm_raw is not None and not pd.isna(qtm_raw)) else np.nan
        qcov_val = float(qcov_raw) if (qcov_raw is not None and not pd.isna(qcov_raw)) else np.nan
        pdb_qual = bool(
            model_qual
            and pdb_src_agrees
            and pd.notna(qtm_val)
            and qtm_val >= 0.50
            and pd.notna(qcov_val)
            and qcov_val >= 0.50
        )

        afdb_hit_src = str(rec.get("afdb_hit_source") or "").strip()
        afdb_src_agrees = bool(mod_src and afdb_hit_src and mod_src == afdb_hit_src)
        af_ev_raw = rec.get("afdb_evalue")
        af_qcov_raw = rec.get("afdb_qcov")
        af_qtm_raw = rec.get("afdb_qtmscore")
        af_ev = float(af_ev_raw) if (af_ev_raw is not None and not pd.isna(af_ev_raw)) else np.nan
        af_qcov = float(af_qcov_raw) if (af_qcov_raw is not None and not pd.isna(af_qcov_raw)) else np.nan
        af_qtm = float(af_qtm_raw) if (af_qtm_raw is not None and not pd.isna(af_qtm_raw)) else np.nan
        afdb_qual = bool(
            model_qual
            and afdb_src_agrees
            and pd.notna(af_ev)
            and af_ev < 1e-3
            and pd.notna(af_qtm)
            and af_qtm >= 0.50
            and pd.notna(af_qcov)
            and af_qcov >= 0.50
        )

        hh_prob_raw = rec.get("hhsearch_prob")
        hh_qcov_raw = rec.get("hhsearch_qcov")
        hh_prob = float(hh_prob_raw) if (hh_prob_raw is not None and not pd.isna(hh_prob_raw)) else np.nan
        hh_qcov = float(hh_qcov_raw) if (hh_qcov_raw is not None and not pd.isna(hh_qcov_raw)) else np.nan
        hh_pf = str(rec.get("hhsearch_pfam") or "missing")
        hh_qual = bool(
            has_known_valid_seq
            and hh_pf not in ("missing", "", "None", "nan")
            and pd.notna(hh_prob)
            and hh_prob >= 80.0
            and (pd.isna(hh_qcov) or hh_qcov >= 0.40)
        )

        ipr_doms = [str(x) for x in (rec.get("interpro_domains") or []) if str(x).strip()] if has_known_valid_seq else []
        ipr_non = [str(x) for x in (rec.get("interpro_non_domain_flags") or []) if str(x).strip()]

        if is_exact_rep:
            if rep_seq_valid and pd.notna(plddt_01):
                if pd.isna(rep_plddt) or (model_qual and not rep_model_qualified) or (model_qual == rep_model_qualified and float(plddt_01) > rep_plddt):
                    rep_plddt = float(plddt_01)
                    rep_model_source = mod_src or "unknown_source"
                    rep_model_seq_prov_status = mod_seq_status
            elif rep_seq_valid and pd.notna(qtm_val) and not pdb_src_agrees and rep_model_seq_prov_status == "no_structural_model":
                rep_model_seq_prov_status = "mismatched_or_missing_model_source_for_hit"
            elif rep_seq_valid and mod_seq_status != "no_structural_model" and rep_model_seq_prov_status == "no_structural_model":
                rep_model_seq_prov_status = mod_seq_status

            rep_model_qualified = rep_model_qualified or model_qual
            if rep_seq_valid and pd.notna(qtm_val) and (pd.isna(rep_pdb_qtm) or (pdb_qual and not rep_qualified_pdb_match) or (pdb_qual == rep_qualified_pdb_match and qtm_val > rep_pdb_qtm)):
                rep_pdb_qtm = qtm_val
                rep_pdb_qcov = qcov_val
                rep_pdb_target = str(rec.get("pdb_target") or "missing")
                if not pdb_src_agrees and mod_seq_matches_locus:
                    rep_model_seq_prov_status = f"mismatched_model_vs_hit_source({mod_src}_vs_{pdb_hit_src})"
            rep_qualified_pdb_match = rep_qualified_pdb_match or pdb_qual

            if rep_seq_valid and pd.notna(af_ev) and (pd.isna(rep_afdb_evalue) or af_ev < rep_afdb_evalue):
                rep_afdb_evalue = af_ev
                rep_afdb_qcov = af_qcov
                rep_afdb_qtm = af_qtm
                rep_afdb_target = str(rec.get("afdb_target") or "missing")
                if afdb_qual:
                    rep_afdb_status = "qualified_afdb_match_with_qtm_and_coverage"
                elif pd.isna(af_qtm):
                    rep_afdb_status = "afdb_3di_aa_alignment_only_qtm_unavailable"
                else:
                    rep_afdb_status = "subthreshold_afdb_match_or_unqualified_model"
            rep_qualified_afdb_match = rep_qualified_afdb_match or afdb_qual

            if rep_seq_valid and pd.notna(hh_prob) and (pd.isna(rep_hh_prob) or hh_prob > rep_hh_prob):
                rep_hh_prob = hh_prob
                rep_hh_qcov = hh_qcov
                rep_hh_pfam = hh_pf
            rep_qualified_hh_profile = rep_qualified_hh_profile or hh_qual

            for d_item in ipr_doms:
                if d_item not in rep_ipr_domains:
                    rep_ipr_domains.append(d_item)
            for nd_item in ipr_non:
                if nd_item not in rep_ipr_nondomain:
                    rep_ipr_nondomain.append(nd_item)
        elif allow_member_summary_separate and has_known_valid_seq:
            if pdb_qual or afdb_qual:
                other_qualified_struct_members.append(f"{loc_gene}[{mod_src}](plddt={plddt_01:.2f},qtm={qtm_val:.2f},qcov={qcov_val:.2f})")
            if hh_qual or ipr_doms:
                dom_tag = hh_pf if hh_qual else ";".join(ipr_doms[:2])
                other_qualified_profile_members.append(f"{loc_gene}({dom_tag})")

    rep_has_own_domain = bool(rep_seq_valid and (rep_qualified_hh_profile or len(rep_ipr_domains) > 0))
    rep_has_own_struct = bool(rep_seq_valid and (rep_qualified_pdb_match or rep_qualified_afdb_match))

    return {
        "representative_gene": representative_gene,
        "rep_sequence_audit_status": str(rep_clean["status"]),
        "rep_aa_length": rep_len,
        "rep_own_model_source": rep_model_source,
        "rep_own_model_seq_provenance_status": rep_model_seq_prov_status,
        "rep_own_plddt_0to1": rep_plddt,
        "rep_own_model_qualified": rep_model_qualified,
        "rep_own_pdb_target": rep_pdb_target,
        "rep_own_pdb_qtmscore": rep_pdb_qtm,
        "rep_own_pdb_qcov": rep_pdb_qcov,
        "rep_own_qualified_pdb_fold_match": rep_qualified_pdb_match,
        "rep_own_afdb_target": rep_afdb_target,
        "rep_own_afdb_evalue": rep_afdb_evalue,
        "rep_own_afdb_qcov": rep_afdb_qcov,
        "rep_own_afdb_qtmscore": rep_afdb_qtm,
        "rep_own_qualified_afdb_fold_match": rep_qualified_afdb_match,
        "rep_own_afdb_support_status": rep_afdb_status,
        "rep_own_has_qualified_whole_fold_support": rep_has_own_struct,
        "rep_own_hhsearch_pfam": rep_hh_pfam,
        "rep_own_hhsearch_prob": rep_hh_prob,
        "rep_own_hhsearch_qcov": rep_hh_qcov,
        "rep_own_interpro_domain_signatures": ";".join(rep_ipr_domains) if rep_ipr_domains else "none",
        "rep_own_interpro_nondomain_disorder_coils": ";".join(rep_ipr_nondomain) if rep_ipr_nondomain else "none",
        "rep_own_has_domain_or_profile_support": rep_has_own_domain,
        "other_member_qualified_structural_evidence": ";".join(other_qualified_struct_members) if other_qualified_struct_members else "none",
        "other_member_qualified_profile_evidence": ";".join(other_qualified_profile_members) if other_qualified_profile_members else "none",
    }


def evaluate_reference_motif_conservation(
    aligned_sequences: Mapping[str, str],
    *,
    reference_id: str | None = None,
    documented_catalytic_sites: Sequence[tuple[int, str, str]] | None = None,
) -> dict[str, Any]:
    """Check motif/catalytic residue conservation ONLY against a stated reference with documented functional residues.

    Implements `stage3-scientific-review.txt` Item 3:
    - Never invents a motif from a family name or generic annotation.
    - `documented_catalytic_sites`: list of `(ref_1based_residue_pos, expected_aa_set_str, site_description)`.
    - If `reference_id` is None/missing from `aligned_sequences` or `documented_catalytic_sites` is empty/None,
      explicitly returns `status = 'curated_reference_or_active_site_residues_absent_functional_specificity_pending'`.
    """
    if (
        not reference_id
        or reference_id not in aligned_sequences
        or not documented_catalytic_sites
    ):
        return {
            "status": "curated_reference_or_active_site_residues_absent_functional_specificity_pending",
            "reference_id": reference_id or "none",
            "n_documented_sites_checked": 0,
            "n_sites_conserved_in_all_members": 0,
            "site_audit_details": "no_curated_active_site_reference_available",
        }

    ref_aln = aligned_sequences[reference_id]
    member_ids = [k for k in aligned_sequences if k != reference_id]
    if not member_ids:
        raise ValueError("evaluate_reference_motif_conservation requires at least 1 family member sequence alongside reference_id")
    aln_len = len(ref_aln)
    for k, s in aligned_sequences.items():
        if len(s) != aln_len:
            raise ValueError(f"Aligned sequence '{k}' length ({len(s)}) != reference alignment length ({aln_len})")

    # Map 1-based ungapped residue position in reference_id to 0-based alignment column
    ref_pos_to_col: dict[int, int] = {}
    res_counter = 0
    for col_idx, ch in enumerate(ref_aln):
        if ch != "-":
            res_counter += 1
            ref_pos_to_col[res_counter] = col_idx

    n_checked = 0
    n_conserved_all = 0
    details: list[str] = []
    for pos_1b, allowed_aas, desc in documented_catalytic_sites:
        n_checked += 1
        if pos_1b not in ref_pos_to_col:
            details.append(f"pos{pos_1b}({desc}):out_of_reference_range")
            continue
        col_i = ref_pos_to_col[pos_1b]
        ref_aa = ref_aln[col_i].upper()
        allowed_set = set(allowed_aas.upper())
        member_aas = [aligned_sequences[m][col_i].upper() for m in member_ids]
        n_match = sum(1 for a in member_aas if a in allowed_set)
        if n_match == len(member_ids) and ref_aa in allowed_set:
            n_conserved_all += 1
        obs_str = "".join(sorted(set(member_aas)))
        details.append(f"pos{pos_1b}_{desc}[ref={ref_aa},allowed={allowed_aas},obs={obs_str},match={n_match}/{len(member_ids)}]")

    if n_conserved_all == n_checked and n_checked > 0:
        status = f"conserved_at_all_{n_checked}_documented_reference_sites"
    elif n_conserved_all > 0:
        status = f"partially_conserved_{n_conserved_all}_of_{n_checked}_documented_reference_sites"
    else:
        status = f"divergent_at_documented_reference_sites_0_of_{n_checked}"

    return {
        "status": status,
        "reference_id": reference_id,
        "n_documented_sites_checked": n_checked,
        "n_sites_conserved_in_all_members": n_conserved_all,
        "site_audit_details": ";".join(details),
    }


def compute_alignment_summary_metrics(aligned_seqs: Mapping[str, str]) -> dict[str, float | int]:
    """Compute exact pairwise amino-acid identity (over ungapped shared positions) and alignment coverage."""
    ids = list(aligned_seqs.keys())
    n_seqs = len(ids)
    if n_seqs == 0:
        return {
            "n_aligned_sequences": 0,
            "alignment_cols": 0,
            "median_ungapped_coverage": np.nan,
            "min_pairwise_aa_identity": np.nan,
            "median_pairwise_aa_identity": np.nan,
            "max_pairwise_aa_identity": np.nan,
        }
    aln_len = len(aligned_seqs[ids[0]])
    for k in ids:
        if len(aligned_seqs[k]) != aln_len:
            raise ValueError(f"Inconsistent alignment length for {k}")
    if aln_len == 0:
        raise ValueError("Empty alignment columns")

    covs = [sum(1 for c in aligned_seqs[k] if c != "-") / float(aln_len) for k in ids]
    if n_seqs == 1:
        return {
            "n_aligned_sequences": 1,
            "alignment_cols": aln_len,
            "median_ungapped_coverage": float(covs[0]),
            "min_pairwise_aa_identity": 1.0,
            "median_pairwise_aa_identity": 1.0,
            "max_pairwise_aa_identity": 1.0,
        }

    pids: list[float] = []
    for i in range(n_seqs):
        s1 = aligned_seqs[ids[i]].upper()
        for j in range(i + 1, n_seqs):
            s2 = aligned_seqs[ids[j]].upper()
            shared = 0
            ident = 0
            for c1, c2 in zip(s1, s2):
                if c1 != "-" and c2 != "-":
                    shared += 1
                    if c1 == c2:
                        ident += 1
            pids.append((ident / float(shared)) if shared > 0 else 0.0)

    arr = np.asarray(pids, dtype="float64")
    return {
        "n_aligned_sequences": n_seqs,
        "alignment_cols": aln_len,
        "median_ungapped_coverage": float(np.median(covs)),
        "min_pairwise_aa_identity": float(np.min(arr)),
        "median_pairwise_aa_identity": float(np.median(arr)),
        "max_pairwise_aa_identity": float(np.max(arr)),
    }


def evaluate_seven_evidence_axes(
    *,
    candidate_role: str,
    rep_final_promotion_gate_pass: bool,
    n_verified_complete_loci: int,
    n_total_loci: int,
    rep_signalp_class: str,
    rep_tmbed_s: bool,
    rep_tm_helices: int,
    rep_tm_beta: int,
    frac_family_exported: float,
    rep_has_own_domain_or_profile: bool,
    rep_own_domain_summary: str,
    rep_has_own_qualified_whole_fold: bool,
    rep_own_fold_summary: str,
    other_member_evidence_summary: str,
    motif_conservation_status: str,
    median_alignment_coverage: float,
    n_verified_orders: int,
    n_verified_independent_genus_blocks: int,
    max_copies_per_genome: int,
    specific_hypothesis_name: str,
    specific_context_obs_adj: int,
    specific_context_p_emp: float,
    n_verified_genus_blocks_with_specific_context: int,
    competing_assignments_summary: str,
    positive_control_evidence_status: str,
    specific_guild_selection_mode: str = "exploratory_posthoc_observed_neighbourhood_winner",
) -> dict[str, str]:
    """Populate the 7 separate evidence axes required by `stage3-scientific-review.txt` Item 5.

    Never promotes repeated observations across >=3 genera to an established mechanism.
    Leaves final biological interpretation to supervising Codex.
    """
    # Axis 1: Sequence audit
    if not rep_final_promotion_gate_pass or n_verified_complete_loci == 0:
        ax1 = f"failed_or_unknown_sequence_audit(verified_loci={n_verified_complete_loci}/{n_total_loci})"
    elif n_verified_complete_loci == n_total_loci:
        ax1 = f"verified_exact_complete_cds_all_members({n_verified_complete_loci}/{n_total_loci})"
    else:
        ax1 = f"verified_rep_and_subset_members({n_verified_complete_loci}/{n_total_loci})"

    # Axis 2: Localisation (explicitly noting signal peptide is localisation evidence only)
    sp_cls = str(rep_signalp_class).strip()
    if sp_cls in ("SP", "LIPO", "TAT", "TATLIPO") or rep_tmbed_s:
        ax2 = (
            f"predicted_export_localisation_only(rep_signalp={sp_cls},rep_tmbed_S={rep_tmbed_s},"
            f"rep_tmH={rep_tm_helices},rep_tmB={rep_tm_beta},fam_export_frac={frac_family_exported:.2f};"
            f"localisation_not_catalytic_or_S_layer_proof)"
        )
    elif (rep_tm_helices + rep_tm_beta) > 0:
        ax2 = f"predicted_transmembrane_no_cleaved_sp(rep_tmH={rep_tm_helices},rep_tmB={rep_tm_beta})"
    elif sp_cls in ("missing", "", "unknown"):
        ax2 = "unknown_missing_localisation_prediction"
    else:
        ax2 = f"predicted_cytoplasmic_or_non_exported(rep_signalp={sp_cls},fam_export_frac={frac_family_exported:.2f})"

    # Axis 3: Own-domain and structure (strictly separating representative-own from other family members)
    ax3 = (
        f"rep_own_domain={rep_has_own_domain_or_profile}({rep_own_domain_summary});"
        f"rep_own_qualified_whole_fold={rep_has_own_qualified_whole_fold}({rep_own_fold_summary});"
        f"other_member_evidence_separate=({other_member_evidence_summary})"
    )

    # Axis 4: Motif and architecture
    cov_str = f"{median_alignment_coverage:.2f}" if pd.notna(median_alignment_coverage) else "NA"
    ax4 = f"msa_median_cov={cov_str};motif_status={motif_conservation_status}"

    # Axis 5: Phylogenetic breadth & paralogy
    if n_verified_independent_genus_blocks < 2:
        breadth_tag = "single_lineage_or_genus_block_restricted"
    elif n_verified_orders < 2:
        breadth_tag = "within_single_order_multi_genus_block"
    else:
        breadth_tag = "multi_order_and_multi_genus_block_recurrent"
    paralog_tag = "paralogous_multicopy_present_avoid_transfer" if max_copies_per_genome > 1 else "single_copy_per_verified_genome"
    ax5 = (
        f"{breadth_tag}(verified_orders={n_verified_orders},"
        f"verified_genus_blocks={n_verified_independent_genus_blocks},{paralog_tag})"
    )

    # Axis 6: Specific-context support (never a function label; separate from generic known-gene proximity;
    # explicitly records exploratory post-hoc guild-selection status per stage3-real-schema-followup.txt)
    if not specific_hypothesis_name or specific_hypothesis_name in ("none", "none_localisation_only"):
        ax6 = "no_specific_pathway_or_guild_context_hypothesis(localisation_screen_only)"
    else:
        p_str = f"{specific_context_p_emp:.4f}" if pd.notna(specific_context_p_emp) else "NA"
        ax6 = (
            f"specific_context_supporting_only[hypothesis={specific_hypothesis_name},"
            f"selection_mode={specific_guild_selection_mode},"
            f"obs_adj={specific_context_obs_adj},p_emp_exploratory={p_str},"
            f"verified_genus_blocks_with_specific_flank={n_verified_genus_blocks_with_specific_context};"
            f"posthoc_guild_selection_dependent_exploratory_never_establishes_mechanism]"
        )

    # Axis 7: Competing / alternative assignments and explicit missing evidence
    missing_items: list[str] = []
    if not rep_final_promotion_gate_pass:
        missing_items.append("verified_complete_representative_cds_missing")
    if not rep_has_own_domain_or_profile:
        missing_items.append("rep_own_sequence_or_profile_domain_absent")
    if not rep_has_own_qualified_whole_fold:
        missing_items.append("rep_own_qualified_whole_fold_pdb_match_absent")
    if "absent" in motif_conservation_status or "pending" in motif_conservation_status:
        missing_items.append("curated_active_site_reference_and_residue_proof_absent")
    if candidate_role == "sequence_annotated_comparator" and positive_control_evidence_status == "positive_control_evidence_unavailable":
        missing_items.append("experimental_or_curated_positive_control_proof_unavailable")
    missing_items.append("no_transcriptomic_proteomic_or_biochemical_activity_assay")

    ax7 = f"competing_assignments=({competing_assignments_summary})"
    missing_str = ";".join(missing_items)

    return {
        "axis1_sequence_audit_status": ax1,
        "axis2_localisation_status": ax2,
        "axis3_own_domain_and_structure_status": ax3,
        "axis4_motif_and_architecture_status": ax4,
        "axis5_phylogenetic_breadth_status": ax5,
        "axis6_specific_context_support_status": ax6,
        "axis7_alternative_and_competing_assignments_status": ax7,
        "explicit_unavailable_or_pending_evidence": missing_str,
        "supervisor_interpretation_handoff": (
            "seven_evidence_axes_reported_no_automated_mechanistic_promotion_final_interpretation_by_supervising_codex"
        ),
    }


def calibrate_leave_order_out_threshold(
    train_scores: np.ndarray | Sequence[float],
    train_top1: np.ndarray | Sequence[float],
    eval_scores: np.ndarray | Sequence[float],
    eval_top1: np.ndarray | Sequence[float],
    *,
    candidate_thresholds: Sequence[float] = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0),
    documented_target_top1_recovery: float = 0.50,
    min_train_passing: int = 20,
    is_held_out_split: bool = True,
) -> dict[str, float | int | str | bool]:
    """Calibrate a `context_score` threshold on `train` and evaluate conditional candidate recovery on `eval`.

    Implements `stage3-scientific-review.txt` Item 4 & Item 6:
    - Uses a single documented target (`documented_target_top1_recovery = 0.50` by default).
    - Never falls back to an arbitrary threshold (like `2.0`) when no threshold achieves the target on `train`;
      instead returns `calibrated_threshold = np.nan` and `threshold_calibration_status = 'no_qualifying_threshold_attained'`.
    - Explicitly distinguishes `leave_order_out_threshold_calibration` (`is_held_out_split=True`) from
      `in_sample_pooled_descriptive_not_held_out` (`is_held_out_split=False`).
    """
    tr_s = np.asarray(train_scores, dtype="float64")
    tr_y = np.asarray(train_top1, dtype="float64")
    ev_s = np.asarray(eval_scores, dtype="float64")
    ev_y = np.asarray(eval_top1, dtype="float64")
    if len(tr_s) != len(tr_y) or len(ev_s) != len(ev_y):
        raise ValueError("Score and top1 arrays must have matching length")

    chosen_t: float = np.nan
    train_pass_n: int = 0
    train_pass_rec: float = np.nan

    for cand_t in sorted(candidate_thresholds):
        m = tr_s >= float(cand_t)
        n_m = int(m.sum())
        if n_m >= min_train_passing:
            rec = float(tr_y[m].mean())
            if rec >= documented_target_top1_recovery:
                chosen_t = float(cand_t)
                train_pass_n = n_m
                train_pass_rec = rec
                break

    eval_mode = (
        "leave_order_out_threshold_calibration"
        if is_held_out_split
        else "in_sample_pooled_descriptive_not_held_out"
    )

    if pd.isna(chosen_t):
        return {
            "evaluation_mode": eval_mode,
            "is_held_out_evaluation": bool(is_held_out_split),
            "documented_target_top1_recovery": float(documented_target_top1_recovery),
            "min_train_passing_required": int(min_train_passing),
            "threshold_calibration_status": "no_qualifying_threshold_attained",
            "calibrated_score_threshold": np.nan,
            "train_n_passing_threshold": 0,
            "train_top1_recovery_at_threshold": np.nan,
            "eval_n_passing_threshold": 0,
            "eval_fraction_passing_threshold": np.nan,
            "eval_conditional_top1_recovery_at_threshold": np.nan,
        }

    ev_m = ev_s >= chosen_t
    ev_pass_n = int(ev_m.sum())
    ev_frac = ev_pass_n / float(len(ev_s)) if len(ev_s) > 0 else np.nan
    ev_rec = float(ev_y[ev_m].mean()) if ev_pass_n > 0 else np.nan

    return {
        "evaluation_mode": eval_mode,
        "is_held_out_evaluation": bool(is_held_out_split),
        "documented_target_top1_recovery": float(documented_target_top1_recovery),
        "min_train_passing_required": int(min_train_passing),
        "threshold_calibration_status": "qualifying_threshold_attained",
        "calibrated_score_threshold": chosen_t,
        "train_n_passing_threshold": train_pass_n,
        "train_top1_recovery_at_threshold": train_pass_rec,
        "eval_n_passing_threshold": ev_pass_n,
        "eval_fraction_passing_threshold": ev_frac,
        "eval_conditional_top1_recovery_at_threshold": ev_rec,
    }


def summarize_candidate_selection_outcome(
    route1_df: pd.DataFrame | None,
    route2_df: pd.DataFrame | None,
    *,
    max_route1_dark: int = 5,
    max_route2_dark: int = 5,
) -> dict[str, Any]:
    """Select up to 10 verified dark candidate families across Route A1 and Route A2.

    Implements `stage3-negative-outcome-review.txt`:
    - Distinguishes missing/empty/unverified upstream tables (which raise `ValueError`) from a
      legitimate scientific negative outcome where `0 <= n_selected_dark <= 10` (including `n = 0`
      families surviving all strict Stage 1b/Stage 3 member-locus verification gates).
    """
    if route1_df is None or route2_df is None or route1_df.empty or route2_df.empty:
        raise ValueError("Missing or empty Stage 1b verified family ranking inputs; cannot run Stage 3 selection.")
    for name, df in (("route1_df", route1_df), ("route2_df", route2_df)):
        for col in ("family_id", "final_candidate_promotion_eligible", "n_verified_complete_interior_member_loci"):
            if col not in df.columns:
                raise ValueError(f"{name} missing required Stage 1b verified column '{col}'")

    r1_dark = route1_df[route1_df.get("inherited_catalogue_dark", pd.Series(True, index=route1_df.index)).astype(bool)].copy()
    r1_elig = r1_dark[r1_dark["final_candidate_promotion_eligible"].astype(bool)].head(max_route1_dark)

    r2_elig = route2_df[
        route2_df["final_candidate_promotion_eligible"].astype(bool)
        & (~route2_df["family_id"].isin(set(r1_elig["family_id"])))
    ].head(max_route2_dark)

    selected_ids = list(r1_elig["family_id"].astype(str)) + list(r2_elig["family_id"].astype(str))
    n_sel = len(selected_ids)
    if n_sel > (max_route1_dark + max_route2_dark):
        raise ValueError(f"Selected {n_sel} dark families, exceeding cap of {max_route1_dark + max_route2_dark}")

    return {
        "selection_outcome_status": (
            "no_eligible_dark_family_after_strict_audit"
            if n_sel == 0
            else "eligible_dark_families_selected"
        ),
        "no_eligible_family": bool(n_sel == 0),
        "n_route1_screened_dark": int(len(r1_dark)),
        "n_route2_screened_dark": int(len(route2_df)),
        "n_route1_selected_dark": int(len(r1_elig)),
        "n_route2_selected_dark": int(len(r2_elig)),
        "n_total_selected_new_dark": int(n_sel),
        "selected_route1_family_ids": list(r1_elig["family_id"].astype(str)),
        "selected_route2_family_ids": list(r2_elig["family_id"].astype(str)),
    }


STAGE1B_VERIFIED_PROPOSALS_REQUIRED_COLS: tuple[str, ...] = (
    "genome",
    "candidate_gene",
    "pathway",
    "missing_reaction",
    "reaction_name",
    "ec",
    "review_priority",
    "inherited_catalogue_dark",
    "family_id",
    "contig_edge_status",
    "final_promotion_gate_pass",
)


def extract_verified_route_a1_reaction_pathway_context(
    prop_df: pd.DataFrame,
) -> dict[str, dict[str, list[str]]]:
    """Validate Stage 1b verified Route A1 proposals schema and extract reaction/pathway context.

    Implements `stage3-real-schema-followup.txt`:
    - Strictly requires the real Stage 1b columns (`missing_reaction`, `pathway`, `reaction_name`, `ec`,
      `review_priority`, `final_promotion_gate_pass`, `contig_edge_status`, `family_id`, `genome`,
      `candidate_gene`, `inherited_catalogue_dark`).
    - Never uses `getattr` defaults for non-existent `rxn_id`/`pathway_id` columns; raises `ValueError`
      if required columns are missing.
    - Filters to verified complete interior priority supporting proposals (`review_priority == True`,
      `final_promotion_gate_pass == True`, `contig_edge_status == 'interior'`) consistently with
      Stage 1b/3 candidate promotion.
    """
    if prop_df is None or prop_df.empty:
        raise ValueError("Stage 1b Route A1 verified proposals table is missing or empty")
    missing_cols = [c for c in STAGE1B_VERIFIED_PROPOSALS_REQUIRED_COLS if c not in prop_df.columns]
    if missing_cols:
        raise ValueError(
            f"Stage 1b Route A1 verified proposals table is missing required columns: {missing_cols}"
        )

    ver_mask = (
        prop_df["review_priority"].astype(bool)
        & prop_df["final_promotion_gate_pass"].astype(bool)
        & (prop_df["contig_edge_status"].astype(str) == "interior")
    )
    ver_df = prop_df.loc[ver_mask]

    out: dict[str, dict[str, list[str]]] = {}

    def _clean_tok(val: object) -> str:
        if val is None or pd.isna(val):
            return ""
        s = str(val).strip()
        if s in ("", "nan", "None", "missing", "<NA>"):
            return ""
        return s

    for row in ver_df.itertuples(index=False):
        fid = _clean_tok(row.family_id)
        if not fid:
            continue
        entry = out.setdefault(
            fid,
            {
                "missing_reactions": [],
                "pathways": [],
                "reaction_names": [],
                "ecs": [],
            },
        )
        for col_attr, target_key in (
            ("missing_reaction", "missing_reactions"),
            ("pathway", "pathways"),
            ("reaction_name", "reaction_names"),
            ("ec", "ecs"),
        ):
            v = _clean_tok(getattr(row, col_attr))
            if v and v not in entry[target_key]:
                entry[target_key].append(v)

    return out


def load_curated_swissprot_gapseq_reference_from_root(
    rxn_id: str,
    gapseq_seq_root: Path,
) -> dict[str, Any]:
    """Inspect `{gapseq_seq_root}/rev/{rxn_id}.fasta` for a reviewed SwissProt reference sequence."""
    clean_rxn = str(rxn_id or "").strip()
    if not clean_rxn:
        return {"has_reviewed_ref": False}
    rev_fa = Path(gapseq_seq_root) / "rev" / f"{clean_rxn}.fasta"
    if not rev_fa.exists() or rev_fa.stat().st_size == 0:
        return {"has_reviewed_ref": False}
    with open(rev_fa, "rt") as fh:
        cur_hdr: str | None = None
        chunks: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if cur_hdr is not None and chunks:
                    break
                cur_hdr = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if not cur_hdr or not chunks:
        return {"has_reviewed_ref": False}
    seq = "".join(chunks)
    parts = cur_hdr.split(None, 1)
    token = parts[0]
    desc = parts[1] if len(parts) > 1 else token
    if "|" in token:
        sp_fields = token.split("|")
        acc = sp_fields[1] if len(sp_fields) >= 2 else token
    else:
        acc = token
    return {
        "has_reviewed_ref": True,
        "rxn_id": clean_rxn,
        "accession": acc,
        "description": desc,
        "sequence": seq,
        "source": "SwissProt_Archaea_rev",
    }


def select_curated_gapseq_reference_for_family(
    family_id: str,
    r1_context_map: Mapping[str, Mapping[str, Sequence[str]]],
    gapseq_seq_root: Path,
) -> dict[str, Any]:
    """Select the first available reviewed SwissProt gapseq reference for `family_id` from verified `missing_reactions`."""
    ctx = r1_context_map.get(str(family_id), {})
    for rxn_id in ctx.get("missing_reactions", []):
        cref = load_curated_swissprot_gapseq_reference_from_root(str(rxn_id), gapseq_seq_root)
        if cref.get("has_reviewed_ref"):
            return cref
    return {"has_reviewed_ref": False}


def run_grouped_contig_circular_nulls(
    member_contig_records: Sequence[Mapping[str, Any]],
    *,
    specific_guild: str,
    specific_guild_selection_mode: str = "exploratory_posthoc_observed_neighbourhood_winner",
    window: int = 5,
    n_perm: int = 1000,
    seed: int = 20261108,
    genus_dereplicate_map: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Group verified member loci by unique real `(genome, contig)`, merge same-contig target members,
    and run `common_science.within_contig_circular_null` separately for the specific guild and generic known flanks.

    Prevents contig pseudoreplication when multiple target members of a family reside on the same contig,
    and explicitly labels post-hoc neighbourhood guild selection as exploratory and selection-dependent.
    When `genus_dereplicate_map` is provided (mapping genus_singleton_block -> canonical top genome per
    `EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt` Item B), also computes the exact null restricted to the deterministic
    top-completeness genome per genus block.
    """
    if not member_contig_records:
        empty_null = {
            "n_target_loci": 0,
            "obs_adjacent": 0,
            "null_mean": 0.0,
            "null_sd": 0.0,
            "enrichment_ratio": np.nan,
            "p_empirical_one_sided": np.nan,
            "p_exact_one_sided": np.nan,
        }
        return {
            "n_unique_contigs": 0,
            "n_member_records": 0,
            "n_merged_same_contig_target_loci": 0,
            "specific_guild": specific_guild,
            "specific_guild_selection_mode": specific_guild_selection_mode,
            "specific_null_inference_scope": "exploratory_posthoc_selection_dependent_not_independent_confirmation",
            "specific_null": dict(empty_null),
            "genus_dereplicated_null": dict(empty_null),
            "generic_null": dict(empty_null),
        }

    grouped: dict[str, dict[str, Any]] = {}
    for rec in member_contig_records:
        genome = str(rec.get("genome", "")).strip()
        contig = str(rec.get("contig", "")).strip()
        if not genome or not contig:
            raise ValueError("Each member_contig_record must have non-empty 'genome' and 'contig'")
        c_guilds = list(rec.get("contig_guilds", []))
        if not c_guilds:
            raise ValueError(f"Empty contig_guilds for ({genome}, {contig})")
        cand_idx = int(rec["cand_idx"])
        if cand_idx < 0 or cand_idx >= len(c_guilds):
            raise ValueError(f"cand_idx={cand_idx} out of bounds for ({genome}, {contig}) of length {len(c_guilds)}")

        cid_key = f"{genome}::{contig}"
        if cid_key not in grouped:
            grouped[cid_key] = {
                "contig_guilds": c_guilds,
                "target_indices": {cand_idx},
            }
        else:
            existing = grouped[cid_key]["contig_guilds"]
            if len(existing) != len(c_guilds) or [set(s) for s in existing] != [set(s) for s in c_guilds]:
                raise ValueError(f"Inconsistent contig_guilds on duplicate contig record for {cid_key}")
            grouped[cid_key]["target_indices"].add(cand_idx)

    contig_ids_flat: list[str] = []
    target_mask_flat: list[bool] = []
    spec_flank_flat: list[bool] = []
    gen_flank_flat: list[bool] = []

    has_spec_guild = bool(specific_guild and specific_guild not in ("none", "none_localisation_only"))
    for cid_key in sorted(grouped.keys()):
        info = grouped[cid_key]
        c_guilds = info["contig_guilds"]
        t_idxs: set[int] = info["target_indices"]
        n_c = len(c_guilds)
        contig_ids_flat.extend([cid_key] * n_c)
        for i in range(n_c):
            is_tgt = i in t_idxs
            target_mask_flat.append(is_tgt)
            # Exclude target member loci themselves from flank_target_mask so multi-member contigs
            # test flanking guild adjacency rather than self/paralog guild labels
            spec_flank_flat.append((not is_tgt) and has_spec_guild and (specific_guild in c_guilds[i]))
            gen_flank_flat.append((not is_tgt) and ("generic_any_known_annotation" in c_guilds[i]))

    spec_res = within_contig_circular_null(
        contig_ids_flat,
        np.asarray(target_mask_flat, dtype=bool),
        np.asarray(spec_flank_flat, dtype=bool),
        window=window,
        n_perm=n_perm,
        seed=seed,
    )
    gen_res = within_contig_circular_null(
        contig_ids_flat,
        np.asarray(target_mask_flat, dtype=bool),
        np.asarray(gen_flank_flat, dtype=bool),
        window=window,
        n_perm=n_perm,
        seed=seed,
    )

    # Genus-dereplicated sensitivity null (EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt Item B & EXP4-NULL-ACCEPTANCE.txt)
    if genus_dereplicate_map is not None:
        derep_contig_ids: list[str] = []
        derep_tgt_mask: list[bool] = []
        derep_spec_flank: list[bool] = []
        derep_genomes: set[str] = set()
        derep_gblocks: set[str] = set()

        for cid_key in sorted(grouped.keys()):
            genome_id = cid_key.split("::", 1)[0]
            gblock = None
            for rec in member_contig_records:
                if str(rec.get("genome")) == genome_id and str(rec.get("contig")) == cid_key.split("::", 1)[1]:
                    gblock = rec.get("genus_block")
                    break

            if gblock is None or str(gblock).strip() == "" or str(gblock) not in genus_dereplicate_map:
                # Missing genus block or unmapped in genus_dereplicate_map -> excluded from dereplicated sensitivity (no silent fallback!)
                continue

            canonical_top_genome = genus_dereplicate_map[str(gblock)]
            if genome_id != canonical_top_genome:
                # Non-top genome in this genus block -> excluded from dereplicated sensitivity
                continue

            info = grouped[cid_key]
            c_guilds = info["contig_guilds"]
            t_idxs = info["target_indices"]
            n_c = len(c_guilds)
            derep_contig_ids.extend([cid_key] * n_c)
            derep_genomes.add(genome_id)
            derep_gblocks.add(str(gblock))
            for i in range(n_c):
                is_tgt = i in t_idxs
                derep_tgt_mask.append(is_tgt)
                derep_spec_flank.append((not is_tgt) and has_spec_guild and (specific_guild in c_guilds[i]))

        if derep_contig_ids and any(derep_tgt_mask):
            genus_derep_res = within_contig_circular_null(
                derep_contig_ids,
                np.asarray(derep_tgt_mask, dtype=bool),
                np.asarray(derep_spec_flank, dtype=bool),
                window=window,
                n_perm=n_perm,
                seed=seed,
            )
            genus_derep_res["status"] = "evaluated"
            genus_derep_res["n_genomes_retained"] = len(derep_genomes)
            genus_derep_res["n_genus_blocks_retained"] = len(derep_gblocks)
        else:
            # Per EXP4-NULL-ACCEPTANCE.txt: if dereplication eliminates all family loci, return n_target_loci=0, p=NaN, status nonestimable_no_retained_loci
            genus_derep_res = {
                "n_target_loci": 0,
                "obs_adjacent": 0,
                "null_mean": 0.0,
                "null_sd": 0.0,
                "enrichment_ratio": np.nan,
                "p_empirical_one_sided": np.nan,
                "p_exact_one_sided": np.nan,
                "status": "nonestimable_no_retained_loci",
                "n_genomes_retained": len(derep_genomes),
                "n_genus_blocks_retained": len(derep_gblocks),
            }
    else:
        # Per EXP4-NULL-ACCEPTANCE.txt: If genus map absent mark sensitivity not run (don't copy full result)
        genus_derep_res = {
            "n_target_loci": 0,
            "obs_adjacent": 0,
            "null_mean": 0.0,
            "null_sd": 0.0,
            "enrichment_ratio": np.nan,
            "p_empirical_one_sided": np.nan,
            "p_exact_one_sided": np.nan,
            "status": "sensitivity_not_run_missing_genus_map",
            "n_genomes_retained": 0,
            "n_genus_blocks_retained": 0,
        }

    n_unique = len(grouped)
    n_total_t = int(sum(len(v["target_indices"]) for v in grouped.values()))
    return {
        "n_unique_contigs": n_unique,
        "n_member_records": len(member_contig_records),
        "n_merged_same_contig_target_loci": int(n_total_t - n_unique),
        "specific_guild": specific_guild,
        "specific_guild_selection_mode": specific_guild_selection_mode,
        "specific_null_inference_scope": "exploratory_posthoc_selection_dependent_not_independent_confirmation",
        "specific_null": spec_res,
        "genus_dereplicated_null": genus_derep_res,
        "generic_null": gen_res,
    }


GSPA_SUMMARY_REQUIRED_COLS: tuple[str, ...] = (
    "genome",
    "mode",
    "proteins",
    "dark",
    "strict_dark",
    "claims",
    "nonroot_0.5_proteins",
    "nonroot_0.5_dark",
    "nonroot_0.5_strict_dark",
    "mf_bp_depth4_0.5_proteins",
    "mf_bp_depth4_0.5_dark",
    "mf_bp_depth4_0.5_strict_dark",
    "nonroot_0.8_proteins",
    "nonroot_0.8_dark",
    "nonroot_0.8_strict_dark",
    "mf_bp_depth4_0.8_proteins",
    "mf_bp_depth4_0.8_dark",
    "mf_bp_depth4_0.8_strict_dark",
)

GSPA_CONTEXT_CHANGES_REQUIRED_COLS: tuple[str, ...] = (
    "genome",
    "variant",
    "mode",
    "increased_claims",
    "increased_proteins",
    "crossed_05",
    "dark_crossed_05",
    "max_increase",
)

METABOLIC_MODEL_SUMMARY_REQUIRED_COLS: tuple[str, ...] = (
    "genome",
    "eq_thermoplasmatota",
    "pathways",
    "predicted_pathways",
    "good_blast_reactions",
    "transporter_hits",
)

C12_SUMMARY_REQUIRED_COLS: tuple[str, ...] = (
    "genome",
    "medium",
    "filled",
    "n_rxn",
    "n_draft",
    "origin_counts",
    "n_genes",
    "growth",
    "fba_status",
)

C18_WINDOWS_REQUIRED_COLS: tuple[str, ...] = (
    "genome",
    "pathway",
    "has_hole",
    "n_window_genes",
    "n_dark",
    "dark_families",
)


def build_manifest_genome_alias_resolver(
    thermo_genomes: Sequence[str],
    manifest_df: pd.DataFrame | None = None,
    canonical_rep_df: pd.DataFrame | None = None,
) -> dict[str, str]:
    """Build a collision-checked 1:1 genome alias resolver backed by `genome_manifest.tsv` / `refined_representatives.tsv`.

    Implements `stage3-gspa-summary-schema-review.txt` and `stage3-final-summary-review.txt`:
    - Registers only literal manifest `mag` and actual listed `drep_fasta_basename` / `path` stem aliases.
    - Never infers or strips `_sub` to create an unsuffixed parent alias (`mag_id[:-4]`), because refined `_sub`
      genomes and unrefined parent bins are different genome sequences.
    - Raises `ValueError` if any alias maps to more than one canonical target genome.
    """
    thermo_list = [str(g).strip() for g in thermo_genomes if str(g).strip()]
    if not thermo_list:
        raise ValueError("thermo_genomes cannot be empty in build_manifest_genome_alias_resolver")
    if len(set(thermo_list)) != len(thermo_list):
        raise ValueError("Duplicate genome IDs in thermo_genomes")

    thermo_set = set(thermo_list)
    alias_to_target: dict[str, str] = {}

    def _register_alias(alias_id: str, target_id: str) -> None:
        a = str(alias_id).strip()
        t = str(target_id).strip()
        if not a or not t:
            return
        if a in alias_to_target and alias_to_target[a] != t:
            raise ValueError(
                f"Ambiguous genome alias '{a}' maps to both '{alias_to_target[a]}' and '{t}'"
            )
        alias_to_target[a] = t

    for g in thermo_list:
        _register_alias(g, g)

    # Register only literal manifest mag and actual listed path stem (never strip _sub to parent)
    if manifest_df is not None and not manifest_df.empty:
        if "mag" not in manifest_df.columns:
            raise ValueError("genome_manifest.tsv missing required column 'mag'")
        for row in manifest_df.to_dict(orient="records"):
            mag_id = str(row["mag"]).strip()
            path_stem = Path(str(row.get("path", f"{mag_id}.fa"))).stem
            candidates = {mag_id, path_stem}
            matched_targets = [c for c in candidates if c in thermo_set]
            if len(matched_targets) > 1:
                raise ValueError(f"Manifest entry '{mag_id}' matches multiple thermo_genomes: {matched_targets}")
            if len(matched_targets) == 1:
                tgt = matched_targets[0]
                for c in candidates:
                    _register_alias(c, tgt)

    if canonical_rep_df is not None and not canonical_rep_df.empty:
        for req_c in ("mag", "drep_fasta_basename"):
            if req_c not in canonical_rep_df.columns:
                raise ValueError(f"refined_representatives.tsv missing required column '{req_c}'")
        for row in canonical_rep_df.to_dict(orient="records"):
            mag_id = str(row["mag"]).strip()
            drep_stem = Path(str(row["drep_fasta_basename"])).stem
            candidates = {mag_id, drep_stem}
            matched_targets = [c for c in candidates if c in thermo_set]
            if len(matched_targets) > 1:
                raise ValueError(f"Canonical entry '{mag_id}' matches multiple thermo_genomes: {matched_targets}")
            if len(matched_targets) == 1:
                tgt = matched_targets[0]
                for c in candidates:
                    _register_alias(c, tgt)

    return alias_to_target


def reconcile_and_summarize_order_metabolic_tables(
    thermo_meta: pd.DataFrame,
    mod_df: pd.DataFrame,
    c12_df: pd.DataFrame,
    c18_df: pd.DataFrame,
    proposals_df: pd.DataFrame,
    alias_map: Mapping[str, str],
    *,
    primary_medium: str = "native",
) -> pd.DataFrame:
    """Validate producer/consumer schemas, reconcile `_sub` aliases via `alias_map`, enforce 1:1 medium cardinality,
    and build the order-stratified metabolic and gap-filling summary without fabricating zeros for missing models.
    """
    for req_c in ("genome", "tax_order", "genus_singleton_block", "checkm_completeness"):
        if req_c not in thermo_meta.columns:
            raise ValueError(f"thermo_meta missing required column '{req_c}'")
    if thermo_meta["genome"].duplicated().any():
        raise ValueError("thermo_meta contains duplicate genome IDs")

    for name, df_in, req_cols in (
        ("metabolic_model_summary.tsv", mod_df, METABOLIC_MODEL_SUMMARY_REQUIRED_COLS),
        ("c12_summary.tsv.gz", c12_df, C12_SUMMARY_REQUIRED_COLS),
        ("c18_pathway_holes/windows.tsv", c18_df, C18_WINDOWS_REQUIRED_COLS),
        ("Stage 1b verified proposals", proposals_df, STAGE1B_VERIFIED_PROPOSALS_REQUIRED_COLS),
    ):
        missing = [c for c in req_cols if c not in df_in.columns]
        if missing:
            raise ValueError(f"{name} missing required columns: {missing}")

    thermo_set = set(thermo_meta["genome"].astype(str))

    # 1. Reconcile metabolic_model_summary.tsv (restrict to eq_thermoplasmatota == True if column has non-EQ rows)
    mod_work = mod_df.copy()
    if "eq_thermoplasmatota" in mod_work.columns:
        mod_work = mod_work[mod_work["eq_thermoplasmatota"].astype(bool)].copy()
    mod_work["reconciled_genome"] = [alias_map.get(str(g).strip(), "") for g in mod_work["genome"]]
    unmapped_mod = sorted(set(mod_work.loc[mod_work["reconciled_genome"] == "", "genome"].astype(str)))
    if unmapped_mod:
        raise ValueError(f"Unmapped EQ Thermoplasmatota genome IDs in metabolic_model_summary.tsv: {unmapped_mod}")
    if mod_work["reconciled_genome"].duplicated().any():
        dups = sorted(set(mod_work.loc[mod_work["reconciled_genome"].duplicated(), "reconciled_genome"]))
        raise ValueError(f"Duplicate reconciled genome IDs in metabolic_model_summary.tsv: {dups}")

    # 2. Reconcile c12_summary.tsv.gz and enforce 1:1 medium selection (never count multiple media as multiple genomes)
    c12_work = c12_df.copy()
    c12_work["reconciled_genome"] = [alias_map.get(str(g).strip(), "") for g in c12_work["genome"]]
    c12_unmapped_all = sorted(set(c12_work.loc[c12_work["reconciled_genome"] == "", "genome"].astype(str)))
    # Separate external comparator accessions (GCA_*, GCF_*, MGYG*) from unexpected unmapped EQ bin IDs
    c12_unmapped_eq = [
        g for g in c12_unmapped_all if not (g.startswith("GCA_") or g.startswith("GCF_") or g.startswith("MGYG"))
    ]
    if c12_unmapped_eq:
        raise ValueError(f"Unmapped EQ genome IDs in c12_summary.tsv.gz: {c12_unmapped_eq}")
    n_c12_external_excluded = len(c12_unmapped_all)

    c12_eq = c12_work[c12_work["reconciled_genome"].isin(thermo_set)].copy()
    c12_eq["medium_clean"] = c12_eq["medium"].astype(str).str.strip()
    if c12_eq.duplicated(subset=["reconciled_genome", "medium_clean"]).any():
        raise ValueError("Duplicate (reconciled_genome, medium) rows in c12_summary.tsv.gz")
    avail_media = set(c12_eq["medium_clean"].unique())
    if primary_medium not in avail_media:
        raise ValueError(f"Requested primary_medium '{primary_medium}' not found in c12_summary.tsv.gz media: {sorted(avail_media)}")

    c12_pri = c12_eq[c12_eq["medium_clean"] == primary_medium].copy()
    if c12_pri["reconciled_genome"].duplicated().any():
        raise ValueError(f"Duplicate genomes in c12_summary.tsv.gz for medium='{primary_medium}'")

    c12_pri_keep = pd.DataFrame(
        {
            "genome": c12_pri["reconciled_genome"].astype(str),
            "has_c12_model": True,
            "c12_primary_filled": c12_pri["filled"].astype(str).str.upper() == "TRUE",
            "n_rxn": pd.to_numeric(c12_pri["n_rxn"], errors="coerce"),
            "n_draft": pd.to_numeric(c12_pri["n_draft"], errors="coerce"),
            "c12_n_genes": pd.to_numeric(c12_pri["n_genes"], errors="coerce"),
            "c12_growth": pd.to_numeric(c12_pri["growth"], errors="coerce"),
        }
    )

    c12_des = c12_eq[c12_eq["medium_clean"] == "desert"].copy()
    if not c12_des.empty:
        c12_des_keep = pd.DataFrame(
            {
                "genome": c12_des["reconciled_genome"].astype(str),
                "c12_desert_filled": c12_des["filled"].astype(str).str.upper() == "TRUE",
                "desert_n_rxn": pd.to_numeric(c12_des["n_rxn"], errors="coerce"),
                "desert_n_draft": pd.to_numeric(c12_des["n_draft"], errors="coerce"),
            }
        )
    else:
        c12_des_keep = pd.DataFrame(columns=["genome", "c12_desert_filled", "desert_n_rxn", "desert_n_draft"])

    # 3. Reconcile c18_pathway_holes/windows.tsv (never fill missing genomes with 0!)
    c18_work = c18_df.copy()
    c18_work["reconciled_genome"] = [alias_map.get(str(g).strip(), "") for g in c18_work["genome"]]
    c18_unmapped = sorted(set(c18_work.loc[c18_work["reconciled_genome"] == "", "genome"].astype(str)))
    c18_unmapped_eq = [
        g for g in c18_unmapped if not (g.startswith("GCA_") or g.startswith("GCF_") or g.startswith("MGYG"))
    ]
    if c18_unmapped_eq:
        raise ValueError(f"Unmapped EQ genome IDs in c18_pathway_holes/windows.tsv: {c18_unmapped_eq}")
    c18_eq = c18_work[c18_work["reconciled_genome"].isin(thermo_set)].copy()
    c18_eq["has_hole_num"] = pd.to_numeric(c18_eq["has_hole"], errors="coerce")
    c18_eq["n_dark_num"] = pd.to_numeric(c18_eq["n_dark"], errors="coerce")
    if c18_eq["has_hole_num"].isna().any() or c18_eq["n_dark_num"].isna().any():
        raise ValueError("c18_pathway_holes/windows.tsv contains non-numeric has_hole or n_dark values")

    c18_g = (
        c18_eq.groupby("reconciled_genome", as_index=False)
        .agg(
            n_c18_pathway_windows=("pathway", "count"),
            n_c18_windows_with_hole=("has_hole_num", lambda s: int((s > 0).sum())),
            n_c18_dark_genes_in_windows=("n_dark_num", lambda s: int(s.sum())),
        )
        .rename(columns={"reconciled_genome": "genome"})
    )
    c18_g["has_c18_windows"] = True

    # 4. Reconcile Stage 1b verified Route A1 proposals
    prop_work = proposals_df.copy()
    prop_work["reconciled_genome"] = [alias_map.get(str(g).strip(), "") for g in prop_work["genome"]]
    unmapped_prop = sorted(set(prop_work.loc[prop_work["reconciled_genome"] == "", "genome"].astype(str)))
    if unmapped_prop:
        raise ValueError(f"Unmapped genome IDs in Stage 1b verified proposals: {unmapped_prop}")

    ver_prop = prop_work[
        prop_work["review_priority"].astype(bool)
        & prop_work["final_promotion_gate_pass"].astype(bool)
        & (prop_work["contig_edge_status"].astype(str) == "interior")
    ].copy()
    if not ver_prop.empty:
        prop_g = (
            ver_prop.groupby("reconciled_genome", as_index=False)
            .agg(
                n_verified_gapfill_proposals=("candidate_gene", "count"),
                n_verified_gapfill_dark_proposals=("inherited_catalogue_dark", lambda s: int(s.astype(bool).sum())),
            )
            .rename(columns={"reconciled_genome": "genome"})
        )
    else:
        prop_g = pd.DataFrame(columns=["genome", "n_verified_gapfill_proposals", "n_verified_gapfill_dark_proposals"])

    # 5. Merge strictly 1:1 onto thermo_meta and assert cardinality invariance
    base = thermo_meta[["genome", "tax_order", "genus_singleton_block", "checkm_completeness"]].copy()
    mod_keep = mod_work[["reconciled_genome", "pathways", "predicted_pathways", "good_blast_reactions", "transporter_hits"]].rename(
        columns={"reconciled_genome": "genome"}
    )
    merged = (
        base.merge(mod_keep, on="genome", how="left", validate="one_to_one")
        .merge(c12_pri_keep, on="genome", how="left", validate="one_to_one")
        .merge(c12_des_keep, on="genome", how="left", validate="one_to_one")
        .merge(c18_g, on="genome", how="left", validate="one_to_one")
        .merge(prop_g, on="genome", how="left", validate="one_to_one")
    )
    if len(merged) != len(base):
        raise ValueError(f"Join cardinality violation: expected {len(base)} genomes, got {len(merged)}")

    merged["has_c12_model"] = merged["has_c12_model"].eq(True)
    merged["has_c18_windows"] = merged["has_c18_windows"].eq(True)
    merged["n_gapfilled_reactions"] = merged["n_rxn"] - merged["n_draft"]
    merged["frac_gapfilled_of_total_rxn"] = merged["n_gapfilled_reactions"] / merged["n_rxn"].replace(0, np.nan)
    merged["desert_n_gapfilled_reactions"] = merged["desert_n_rxn"] - merged["desert_n_draft"]

    # Proposals are evaluated across all genomes with gapfill models; for genomes with a model and 0 proposals, count is 0,
    # whereas for genomes missing a gapseq model altogether, keep proposal counts NaN so missing models are never conflated with 0
    for col in ("n_verified_gapfill_proposals", "n_verified_gapfill_dark_proposals"):
        merged[col] = np.where(
            merged["has_c12_model"],
            pd.to_numeric(merged[col], errors="coerce").fillna(0.0),
            np.nan,
        )

    order_rows: list[dict[str, Any]] = []
    for ord_name, grp in merged.groupby("tax_order", sort=True):
        missing_model_g = sorted(grp.loc[~grp["has_c12_model"], "genome"].astype(str).tolist())
        c18_sub = grp[grp["has_c18_windows"]]
        mod_sub = grp[grp["has_c12_model"]]
        order_rows.append(
            {
                "tax_order": str(ord_name),
                "n_genomes": int(len(grp)),
                "n_independent_genus_blocks": int(grp["genus_singleton_block"].nunique()),
                "n_genomes_with_gapseq_model": int(grp["has_c12_model"].sum()),
                "n_genomes_missing_gapseq_model": int((~grp["has_c12_model"]).sum()),
                "missing_gapseq_model_genomes": ";".join(missing_model_g) if missing_model_g else "none",
                "selected_c12_medium": primary_medium,
                "n_genomes_filled_native": int(grp["c12_primary_filled"].eq(True).sum()),
                "n_genomes_filled_desert": int(grp["c12_desert_filled"].eq(True).sum()),
                "median_checkm_completeness": float(pd.to_numeric(grp["checkm_completeness"], errors="coerce").median()),
                "median_sequence_supported_draft_rxn": float(mod_sub["n_draft"].median()) if not mod_sub.empty else np.nan,
                "median_gapfilled_rxn_added": float(mod_sub["n_gapfilled_reactions"].median()) if not mod_sub.empty else np.nan,
                "median_total_model_rxn": float(mod_sub["n_rxn"].median()) if not mod_sub.empty else np.nan,
                "median_fraction_gapfilled_rxn": float(mod_sub["frac_gapfilled_of_total_rxn"].median()) if not mod_sub.empty else np.nan,
                "median_desert_gapfilled_rxn_added": float(mod_sub["desert_n_gapfilled_reactions"].median()) if not mod_sub.empty else np.nan,
                "median_desert_total_model_rxn": float(mod_sub["desert_n_rxn"].median()) if not mod_sub.empty else np.nan,
                "median_good_blast_reactions": float(pd.to_numeric(grp["good_blast_reactions"], errors="coerce").median()),
                "median_predicted_pathways": float(pd.to_numeric(grp["predicted_pathways"], errors="coerce").median()),
                "n_genomes_with_c18_windows": int(len(c18_sub)),
                "total_c18_pathway_windows": int(c18_sub["n_c18_pathway_windows"].sum()) if not c18_sub.empty else 0,
                "total_c18_windows_with_hole": int(c18_sub["n_c18_windows_with_hole"].sum()) if not c18_sub.empty else 0,
                "total_c18_dark_genes_in_windows": int(c18_sub["n_c18_dark_genes_in_windows"].sum()) if not c18_sub.empty else 0,
                "total_verified_complete_gapfill_proposals": int(mod_sub["n_verified_gapfill_proposals"].sum()) if not mod_sub.empty else 0,
                "total_verified_complete_dark_gapfill_proposals": int(mod_sub["n_verified_gapfill_dark_proposals"].sum()) if not mod_sub.empty else 0,
                "c12_external_comparator_genomes_excluded_n": int(n_c12_external_excluded),
                "unmapped_eq_source_ids": "none",
                "evidence_separation_note": (
                    "sequence_supported_draft_reactions_kept_strictly_separate_from_stoichiometric_gapfilled_additions_and_context_proposals;"
                    "missing_gapseq_models_recorded_explicitly_never_filled_as_zero"
                ),
            }
        )
    out_df = pd.DataFrame(order_rows)
    if int(out_df["n_genomes"].sum()) != int(len(thermo_meta)):
        raise ValueError(
            f"Order summary genome count mismatch: sum(n_genomes)={int(out_df['n_genomes'].sum())} != {len(thermo_meta)}"
        )
    return out_df


def build_gspa_pilot_3genome_audit_table(
    gspa_sum_df: pd.DataFrame,
    gspa_ctx_df: pd.DataFrame,
    thermo_meta: pd.DataFrame,
    alias_map: Mapping[str, str],
    gspa_gap_pilot_root: Path | None = None,
) -> pd.DataFrame:
    """Build the 3-genome GSPA pilot descriptive audit using exact column access (never `itertuples`/`getattr`).

    Implements `stage3-gspa-summary-schema-review.txt`:
    - Strictly validates `GSPA_SUMMARY_REQUIRED_COLS` (including dotted columns `nonroot_0.5_proteins`,
      `nonroot_0.5_dark`, `nonroot_0.5_strict_dark`, `mf_bp_depth4_0.5_proteins`, `mf_bp_depth4_0.5_dark`,
      `nonroot_0.8_proteins`, `nonroot_0.8_dark`, etc.) and `GSPA_CONTEXT_CHANGES_REQUIRED_COLS`.
    - Raises `ValueError` if any required column is missing or contains `NaN`/`None` (no zero fallback).
    - Preserves `(genome, mode)` rows and separates `context` vs `shuffled` modes explicitly.
    - Never equates prediction counts to validated functional precision.
    """
    if gspa_sum_df is None or gspa_sum_df.empty:
        raise ValueError("gspa summary.tsv is missing or empty")
    if gspa_ctx_df is None or gspa_ctx_df.empty:
        raise ValueError("gspa context-changes.tsv is missing or empty")

    missing_sum = [c for c in GSPA_SUMMARY_REQUIRED_COLS if c not in gspa_sum_df.columns]
    if missing_sum:
        raise ValueError(f"gspa summary.tsv missing required columns: {missing_sum}")
    missing_ctx = [c for c in GSPA_CONTEXT_CHANGES_REQUIRED_COLS if c not in gspa_ctx_df.columns]
    if missing_ctx:
        raise ValueError(f"gspa context-changes.tsv missing required columns: {missing_ctx}")

    if gspa_sum_df[list(GSPA_SUMMARY_REQUIRED_COLS)].isna().any().any():
        bad_cols = [c for c in GSPA_SUMMARY_REQUIRED_COLS if gspa_sum_df[c].isna().any()]
        raise ValueError(f"gspa summary.tsv contains missing/NaN values in required columns: {bad_cols}")
    if gspa_ctx_df[list(GSPA_CONTEXT_CHANGES_REQUIRED_COLS)].isna().any().any():
        bad_cols = [c for c in GSPA_CONTEXT_CHANGES_REQUIRED_COLS if gspa_ctx_df[c].isna().any()]
        raise ValueError(f"gspa context-changes.tsv contains missing/NaN values in required columns: {bad_cols}")

    if gspa_sum_df.duplicated(subset=["genome", "mode"]).any():
        raise ValueError("Duplicate (genome, mode) rows in gspa summary.tsv")
    if gspa_ctx_df.duplicated(subset=["genome", "variant", "mode"]).any():
        raise ValueError("Duplicate (genome, variant, mode) rows in gspa context-changes.tsv")

    tax_lookup = dict(zip(thermo_meta["genome"].astype(str), thermo_meta["tax_order"].astype(str)))

    # Index context-changes.tsv by (reconciled_genome, variant, mode)
    ctx_lookup: dict[tuple[str, str, str], dict[str, Any]] = {}
    for cr in gspa_ctx_df.to_dict(orient="records"):
        raw_g = str(cr["genome"]).strip()
        rec_g = alias_map.get(raw_g, raw_g)
        var = str(cr["variant"]).strip()
        cmode = str(cr["mode"]).strip()
        ctx_lookup[(rec_g, var, cmode)] = {
            "increased_claims": int(cr["increased_claims"]),
            "increased_proteins": int(cr["increased_proteins"]),
            "crossed_05": int(cr["crossed_05"]),
            "dark_crossed_05": int(cr["dark_crossed_05"]),
            "max_increase": float(cr["max_increase"]),
        }

    mode_to_ctx_key: dict[str, tuple[str, str]] = {
        "context": ("fast", "context"),
        "shuffled": ("fast", "shuffled"),
        "cnn_context": ("cnn", "context"),
        "cnn_shuffled": ("cnn", "shuffled"),
    }

    pilot_rows: list[dict[str, Any]] = []
    for r in gspa_sum_df.to_dict(orient="records"):
        raw_gid = str(r["genome"]).strip()
        gid = alias_map.get(raw_gid, "")
        if not gid or gid not in tax_lookup:
            raise ValueError(f"Unmapped GSPA pilot genome '{raw_gid}' not found in Thermoplasmatota metadata")
        mode_str = str(r["mode"]).strip()
        var_family = "cnn" if mode_str.startswith("cnn") else "fast"

        ctx_pair_obs = ctx_lookup.get((gid, var_family, "context"))
        ctx_pair_shuf = ctx_lookup.get((gid, var_family, "shuffled"))
        if ctx_pair_obs is None or ctx_pair_shuf is None:
            raise ValueError(
                f"Missing paired context/shuffled rows in context-changes.tsv for genome='{gid}', variant='{var_family}'"
            )

        if mode_str in mode_to_ctx_key:
            v_k, m_k = mode_to_ctx_key[mode_str]
            row_ctx = ctx_lookup[(gid, v_k, m_k)]
            cond_type = "observed_contig_context" if m_k == "context" else "shuffled_contig_order_control"
            m_inc_claims: float | int = int(row_ctx["increased_claims"])
            m_inc_prots: float | int = int(row_ctx["increased_proteins"])
            m_cross_05: float | int = int(row_ctx["crossed_05"])
            m_dk_cross_05: float | int = int(row_ctx["dark_crossed_05"])
            m_max_inc: float = float(row_ctx["max_increase"])
        else:
            cond_type = "unconditioned_or_base_mode"
            m_inc_claims = np.nan
            m_inc_prots = np.nan
            m_cross_05 = np.nan
            m_dk_cross_05 = np.nan
            m_max_inc = np.nan

        obs_n: float | int = np.nan
        shuf_n: float | int = np.nan
        obs_top: float = np.nan
        shuf_top: float = np.nan
        obs_minus_shuf: float | int = np.nan
        if gspa_gap_pilot_root is not None:
            obs_p = Path(gspa_gap_pilot_root) / raw_gid / "observed_suggestions.tsv"
            shuf_p = Path(gspa_gap_pilot_root) / raw_gid / "shuffled_suggestions.tsv"
            if not obs_p.exists() and gid != raw_gid:
                obs_p = Path(gspa_gap_pilot_root) / gid / "observed_suggestions.tsv"
                shuf_p = Path(gspa_gap_pilot_root) / gid / "shuffled_suggestions.tsv"
            if obs_p.exists() and shuf_p.exists():
                obs_df = pd.read_csv(obs_p, sep="\t")
                shuf_df = pd.read_csv(shuf_p, sep="\t")
                obs_sc_col = "suggestion_score" if "suggestion_score" in obs_df.columns else "score"
                shuf_sc_col = "suggestion_score" if "suggestion_score" in shuf_df.columns else "score"
                obs_n = int(len(obs_df))
                shuf_n = int(len(shuf_df))
                obs_top = (
                    float(pd.to_numeric(obs_df[obs_sc_col], errors="coerce").max())
                    if (not obs_df.empty and obs_sc_col in obs_df.columns)
                    else np.nan
                )
                shuf_top = (
                    float(pd.to_numeric(shuf_df[shuf_sc_col], errors="coerce").max())
                    if (not shuf_df.empty and shuf_sc_col in shuf_df.columns)
                    else np.nan
                )
                obs_minus_shuf = int(len(obs_df) - len(shuf_df))

        pilot_rows.append(
            {
                "pilot_genome": gid,
                "gspa_mode": mode_str,
                "variant_family": var_family,
                "context_conditioning_type": cond_type,
                "tax_order": tax_lookup[gid],
                "gspa_total_proteins": int(r["proteins"]),
                "gspa_dark_proteins": int(r["dark"]),
                "gspa_strict_dark_proteins": int(r["strict_dark"]),
                "gspa_claims": int(r["claims"]),
                "gspa_nonroot_05_proteins": int(r["nonroot_0.5_proteins"]),
                "gspa_nonroot_05_dark_proteins": int(r["nonroot_0.5_dark"]),
                "gspa_nonroot_05_strict_dark_proteins": int(r["nonroot_0.5_strict_dark"]),
                "gspa_mf_bp_depth4_05_proteins": int(r["mf_bp_depth4_0.5_proteins"]),
                "gspa_mf_bp_depth4_05_dark_proteins": int(r["mf_bp_depth4_0.5_dark"]),
                "gspa_mf_bp_depth4_05_strict_dark_proteins": int(r["mf_bp_depth4_0.5_strict_dark"]),
                "gspa_nonroot_08_proteins": int(r["nonroot_0.8_proteins"]),
                "gspa_nonroot_08_dark_proteins": int(r["nonroot_0.8_dark"]),
                "gspa_nonroot_08_strict_dark_proteins": int(r["nonroot_0.8_strict_dark"]),
                "gspa_mf_bp_depth4_08_proteins": int(r["mf_bp_depth4_0.8_proteins"]),
                "gspa_mf_bp_depth4_08_dark_proteins": int(r["mf_bp_depth4_0.8_dark"]),
                "gspa_mf_bp_depth4_08_strict_dark_proteins": int(r["mf_bp_depth4_0.8_strict_dark"]),
                "mode_increased_claims": m_inc_claims,
                "mode_increased_proteins": m_inc_prots,
                "mode_crossed_05": m_cross_05,
                "mode_dark_crossed_05": m_dk_cross_05,
                "mode_max_increase": m_max_inc,
                "variant_context_increased_proteins": int(ctx_pair_obs["increased_proteins"]),
                "variant_shuffled_increased_proteins": int(ctx_pair_shuf["increased_proteins"]),
                "variant_context_crossed_05": int(ctx_pair_obs["crossed_05"]),
                "variant_shuffled_crossed_05": int(ctx_pair_shuf["crossed_05"]),
                "variant_context_dark_crossed_05": int(ctx_pair_obs["dark_crossed_05"]),
                "variant_shuffled_dark_crossed_05": int(ctx_pair_shuf["dark_crossed_05"]),
                "variant_context_minus_shuffled_crossed_05": int(ctx_pair_obs["crossed_05"] - ctx_pair_shuf["crossed_05"]),
                "variant_context_minus_shuffled_dark_crossed_05": int(
                    ctx_pair_obs["dark_crossed_05"] - ctx_pair_shuf["dark_crossed_05"]
                ),
                "gap_pilot_observed_suggestions_n": obs_n,
                "gap_pilot_shuffled_suggestions_n": shuf_n,
                "gap_pilot_observed_top_score": obs_top,
                "gap_pilot_shuffled_top_score": shuf_top,
                "gap_pilot_observed_minus_shuffled_count": obs_minus_shuf,
                "prediction_status_note": (
                    "uncalibrated_gspa_prediction_counts_only_not_validated_functional_precision;"
                    "dotted_column_schema_verified_without_fallback"
                ),
                "scope_and_limitation_note": (
                    "three-genome exploratory pilot; observed and shuffled counts reported separately; "
                    "prediction counts do not establish functional precision or statistical separation"
                ),
            }
        )
    return pd.DataFrame(pilot_rows)



