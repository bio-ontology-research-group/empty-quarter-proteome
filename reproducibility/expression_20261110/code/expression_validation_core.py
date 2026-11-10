#!/usr/bin/env python3
"""
Core validation, sequence canonicalization, metadata-verified MS-GF+ PSM & peptide Q-value
filtering, multipattern background I/L specificity, unique spectrum-union quantification,
and n=3 paired-plot statistical functions for Analysis 4 (Karaoz et al. 2026 external
metaproteomic expression validation).

Complies with:
- DESIGN.txt
- ROOT-ANALYSIS-CONTRACT.txt
- ROOT-STAGE1-REVIEW.txt
- ROOT-DOWNLOAD-REVIEW.txt
- ROOT-TRANSLATION-CORRECTION.txt
- ROOT-PSM-SCORES.txt
- ROOT-BACKGROUND-REVIEW.txt
"""
import hashlib
import itertools
import math
import re
from collections import Counter, defaultdict
from typing import Any, Callable, Dict, Iterable, List, Optional, Set, Tuple

VALID_AA_CHARS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def canonicalize_and_validate_protein_seq(raw_seq: str) -> Dict[str, Any]:
    """
    Canonicalize and validate a protein sequence according to ROOT-TRANSLATION-CORRECTION.txt:
    - Preserves exact raw_sha256 and raw_md5 on the input string.
    - Strips ONLY whitespace (spaces, tabs, CR, LF) and uppercases.
    - Strips at most ONE terminal stop codon ('*').
    - NEVER deletes internal stop codons ('*'), gaps ('-', '.'), digits, or invalid characters.
    - Flags any sequence with internal '*' or non-AA characters as invalid (is_valid=False)
      so it can never be silently repaired or used in valid hash joins.
    """
    if raw_seq is None:
        raise ValueError("Protein sequence cannot be None")
    raw_bytes = raw_seq.encode("utf-8", errors="replace")
    raw_sha256 = hashlib.sha256(raw_bytes).hexdigest()
    raw_md5 = hashlib.md5(raw_bytes).hexdigest()

    no_ws = re.sub(r"\s+", "", raw_seq).upper()
    had_terminal_stop = no_ws.endswith("*")
    canonical = no_ws[:-1] if had_terminal_stop else no_ws

    n_internal_stops = canonical.count("*")
    invalid_chars = sorted({ch for ch in canonical if ch not in VALID_AA_CHARS and ch != "*"})

    if not canonical:
        validity_status = "invalid_empty_sequence"
        is_valid = False
    elif n_internal_stops > 0:
        validity_status = "invalid_internal_stop_codon"
        is_valid = False
    elif invalid_chars:
        validity_status = f"invalid_non_aa_characters:{''.join(invalid_chars)}"
        is_valid = False
    else:
        validity_status = "valid_protein_sequence"
        is_valid = True

    canon_bytes = canonical.encode("ascii", errors="replace")
    canon_sha256 = hashlib.sha256(canon_bytes).hexdigest()
    canon_md5 = hashlib.md5(canon_bytes).hexdigest()

    return {
        "canonical_seq": canonical,
        "aa_length": len(canonical),
        "is_valid": is_valid,
        "validity_status": validity_status,
        "had_terminal_stop": had_terminal_stop,
        "n_internal_stops": n_internal_stops,
        "invalid_chars": "".join(invalid_chars),
        "canonical_sha256": canon_sha256,
        "canonical_md5": canon_md5,
        "raw_sha256": raw_sha256,
        "raw_md5": raw_md5,
    }


def translate_pyrodigal_gene_consistent(
    gene: Any,
    requested_translation_table: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Translate a Pyrodigal predicted gene using the EXACT genetic code selected during gene calling
    (ROOT-TRANSLATION-CORRECTION.txt & ROOT-PSM-SCORES.txt).
    - Requires gene.translation_table attribute explicitly (no silent fallback to 11).
    - Raises ValueError if caller attempts to force a different translation table than gene.translation_table.
    - Calls gene.translate() with default None so Pyrodigal uses its calling table.
    """
    if not hasattr(gene, "translation_table") or getattr(gene, "translation_table") is None:
        raise AttributeError("Gene object is missing required 'translation_table' provenance attribute")
    called_table = int(gene.translation_table)
    if requested_translation_table is not None and int(requested_translation_table) != called_table:
        raise ValueError(
            f"Genetic code mismatch: gene was called with translation_table={called_table}, "
            f"cannot force translation with translation_table={requested_translation_table}."
        )
    raw_aa = gene.translate()
    res = canonicalize_and_validate_protein_seq(raw_aa)
    res["called_translation_table"] = called_table
    return res


def build_unique_alias_map(
    records: Iterable[Dict[str, Any]],
    id_field: str,
    normalizer: Callable[[str], str],
    strict: bool = True,
) -> Tuple[Dict[str, Dict[str, Any]], List[Dict[str, str]]]:
    """
    Build a normalized ID lookup dictionary while detecting and recording any alias collisions.
    Raises ValueError on collision when strict=True instead of silently overwriting entries.
    """
    mapped: Dict[str, Dict[str, Any]] = {}
    collisions: List[Dict[str, str]] = []
    for rec in records:
        if id_field not in rec or not str(rec[id_field]).strip():
            raise ValueError(f"Record missing required id_field '{id_field}': {rec}")
        raw_id = str(rec[id_field]).strip()
        norm_id = normalizer(raw_id)
        if not norm_id:
            raise ValueError(f"Normalizer produced empty key for raw_id='{raw_id}'")
        if norm_id in mapped:
            existing_raw = str(mapped[norm_id][id_field]).strip()
            if existing_raw != raw_id:
                col = {
                    "normalized_id": norm_id,
                    "first_raw_id": existing_raw,
                    "colliding_raw_id": raw_id,
                }
                collisions.append(col)
                if strict:
                    raise ValueError(
                        f"Normalized ID collision for '{norm_id}': '{existing_raw}' vs '{raw_id}'"
                    )
        else:
            mapped[norm_id] = rec
    return mapped, collisions


def lookup_strict(mapping: Dict[str, Any], key: str, context: str) -> Any:
    """
    Strict dictionary lookup that raises KeyError with context rather than inventing 'unknown' or {}.
    """
    if key not in mapping:
        raise KeyError(f"Unexpected/unmapped identifier '{key}' in {context}")
    return mapping[key]


def verify_sequence_hash_join(
    source_id: str,
    source_seq_or_sha256: str,
    target_id: str,
    target_seq_or_sha256: str,
    is_raw_sequence: bool = False,
) -> bool:
    """
    Verify that a join between a query/genome ORF and a search-database protein has an identical
    valid canonical SHA-256 sequence hash. Rejects replaced sequences even if IDs or counts match.
    """
    if is_raw_sequence:
        s_info = canonicalize_and_validate_protein_seq(source_seq_or_sha256)
        t_info = canonicalize_and_validate_protein_seq(target_seq_or_sha256)
        if not s_info["is_valid"] or not t_info["is_valid"]:
            return False
        s_sha = s_info["canonical_sha256"]
        t_sha = t_info["canonical_sha256"]
    else:
        s_sha = source_seq_or_sha256.strip().lower()
        t_sha = target_seq_or_sha256.strip().lower()
    if not s_sha or not t_sha or len(s_sha) != 64 or len(t_sha) != 64:
        raise ValueError(f"Invalid SHA-256 hash in join verification between {source_id} and {target_id}")
    return s_sha == t_sha


def evaluate_homology_alignment(
    pident: float,
    qstart: int,
    qend: int,
    qlen: int,
    sstart: int,
    send: int,
    slen: int,
    evalue: float,
    bitscore: float,
) -> Dict[str, Any]:
    """
    Validate alignment coordinate ranges, finite e-value/bitscore, bidirectional alignment span
    coverage, and short-domain rejection gates (ROOT-STAGE1-REVIEW.txt items 3 & 4).
    Note: Full-length sequence similarity establishes homology eligibility, not proof of orthology.
    """
    for name, val in [("pident", pident), ("evalue", evalue), ("bitscore", bitscore)]:
        if not math.isfinite(val):
            raise ValueError(f"Non-finite {name}={val}")
    if pident < 0.0 or pident > 100.0:
        raise ValueError(f"pident out of [0, 100] range: {pident}")
    if evalue < 0.0:
        raise ValueError(f"Negative evalue: {evalue}")
    if bitscore < 0.0:
        raise ValueError(f"Negative bitscore: {bitscore}")
    if qlen <= 0 or slen <= 0:
        raise ValueError(f"Invalid sequence lengths: qlen={qlen}, slen={slen}")
    if qstart < 1 or qend < qstart or qend > qlen:
        raise ValueError(f"Invalid query alignment coordinates: {qstart}..{qend} (qlen={qlen})")
    if sstart < 1 or send < sstart or send > slen:
        raise ValueError(f"Invalid target alignment coordinates: {sstart}..{send} (slen={slen})")

    q_span = qend - qstart + 1
    t_span = send - sstart + 1
    qcov = q_span / float(qlen)
    tcov = t_span / float(slen)

    if not (0.0 < qcov <= 1.0 and 0.0 < tcov <= 1.0):
        raise ValueError(f"Computed coverage out of (0, 1] range: qcov={qcov}, tcov={tcov}")

    min_span = min(q_span, t_span)
    is_short_protein_full_match = (min(qlen, slen) < 70) and (qcov >= 0.80) and (tcov >= 0.80)
    if min_span < 50 and not is_short_protein_full_match:
        gate_status = "rejected_short_or_generic_domain_only"
    elif evalue <= 1e-10 and pident >= 30.0 and qcov >= 0.70 and tcov >= 0.70:
        gate_status = "strict_full_length_homolog"
    elif evalue <= 1e-5 and pident >= 25.0 and qcov >= 0.50 and tcov >= 0.50:
        gate_status = "moderate_bidirectional_homolog"
    elif evalue <= 1e-5 and (qcov < 0.50 or tcov < 0.50):
        gate_status = "rejected_partial_or_generic_domain_only"
    else:
        gate_status = "rejected_weak_similarity"

    return {
        "q_span": q_span,
        "t_span": t_span,
        "qcov": round(qcov, 6),
        "tcov": round(tcov, 6),
        "pident": round(pident, 3),
        "evalue": evalue,
        "bitscore": bitscore,
        "homology_gate_status": gate_status,
        "passes_strict_homology_gate": gate_status == "strict_full_length_homolog",
        "passes_moderate_homology_gate": gate_status in (
            "strict_full_length_homolog",
            "moderate_bidirectional_homolog",
        ),
    }


def clean_peptide_sequence(raw_peptide: str) -> str:
    """
    Extract clean uppercase amino-acid sequence from an MS-GF+ / mzTab peptide string,
    stripping flanking cleavage notation (e.g. 'K.PEPTIDER.A') and mass-modification tags
    (e.g. '+15.995', '(Oxidation)', '[UNIMOD:35]'), while validating that remaining letters
    are valid amino-acid codes.
    """
    if not raw_peptide or not raw_peptide.strip():
        raise ValueError("Empty peptide sequence")
    s = raw_peptide.strip()
    m = re.match(r"^[A-Z\-]\.(.+)\.[A-Z\-]$", s)
    if m:
        s = m.group(1)
    s = re.sub(r"\[[^\]]*\]|\([^\)]*\)|[+-]\d+(?:\.\d+)?", "", s)
    s = s.upper()
    invalid = [ch for ch in s if ch not in VALID_AA_CHARS]
    if invalid or not s:
        raise ValueError(f"Invalid characters in peptide sequence '{raw_peptide}' -> '{s}'")
    return s


def normalize_peptide_il(peptide_seq: str) -> str:
    """Normalize isobaric Isoleucine (I) and Leucine (L) to 'L' for MS specificity and deduplication."""
    clean = clean_peptide_sequence(peptide_seq)
    return clean.replace("I", "L")


def stream_multipattern_peptide_search_records(
    records: Iterable[Tuple[str, str]],
    candidate_clean_peptides: Iterable[str],
    batch_size: int = 4096,
) -> Dict[str, Dict[str, Set[str]]]:
    """
    Efficient multipattern streaming search of candidate peptides (both exact and I/L-normalized)
    against an arbitrary stream of (protein_id, protein_sequence) records (ROOT-BACKGROUND-REVIEW.txt).
    Uses chunked C-engine Boyer-Moore-Horspool substring filtering (il_pat in joined_chunk_il) so
    chunks of 4,096 background sequences are filtered in microseconds in C without regex backtracking.
    """
    clean_peps = sorted({clean_peptide_sequence(p) for p in candidate_clean_peptides})
    il_to_cleans: Dict[str, List[str]] = defaultdict(list)
    for cp in clean_peps:
        il_to_cleans[cp.replace("I", "L")].append(cp)

    il_patterns = sorted(il_to_cleans.keys(), key=len, reverse=True)
    results: Dict[str, Dict[str, Set[str]]] = {
        cp: {"exact_matches": set(), "il_matches": set()} for cp in clean_peps
    }
    if not il_patterns:
        return results

    def _process_chunk(chunk_items: List[Tuple[str, str, str]]) -> None:
        joined_il = "\n".join(x[2] for x in chunk_items)
        for il_pat in il_patterns:
            if il_pat in joined_il:
                c_list = il_to_cleans[il_pat]
                for pid, seq_up, seq_il in chunk_items:
                    if il_pat in seq_il:
                        for cp in c_list:
                            results[cp]["il_matches"].add(pid)
                            if cp in seq_up:
                                results[cp]["exact_matches"].add(pid)

    chunk: List[Tuple[str, str, str]] = []
    for pid, raw_seq in records:
        seq_up = raw_seq.strip().upper()
        if seq_up.endswith("*"):
            seq_up = seq_up[:-1]
        seq_il = seq_up.replace("I", "L")
        chunk.append((pid, seq_up, seq_il))
        if len(chunk) >= batch_size:
            _process_chunk(chunk)
            chunk = []
    if chunk:
        _process_chunk(chunk)

    return results


def classify_peptide_specificity_against_background(
    raw_peptide: str,
    target_protein_ids: Set[str],
    background_proteins: Dict[str, str],
    protein_to_mags: Dict[str, Set[str]],
    protein_to_families: Optional[Dict[str, Set[str]]] = None,
    precomputed_exact_matches: Optional[Set[str]] = None,
    precomputed_il_matches: Optional[Set[str]] = None,
) -> Dict[str, Any]:
    """
    Evaluate a peptide's exact and I/L-isobaric specificity against a protein background
    (ROOT-ANALYSIS-CONTRACT.txt Rule 3 & 7, ROOT-PSM-SCORES.txt, ROOT-BACKGROUND-REVIEW.txt):
    - Maps both exact sequence and I/L-normalized sequence across all background proteins.
    - Detects isobaric I/L collisions (where a background protein matches after I->L normalization
      but differs in exact I/L residues).
    - Separately evaluates single-MAG specificity and single-family specificity:
      within-MAG specificity does NOT imply family specificity if multiple paralog families
      (or an unassigned paralog) in the same MAG share the peptide!
    """
    clean_pep = clean_peptide_sequence(raw_peptide)
    il_pep = clean_pep.replace("I", "L")

    exact_hits: Set[str] = set(precomputed_exact_matches or set())
    il_hits: Set[str] = set(precomputed_il_matches or set())
    for pid, seq in background_proteins.items():
        seq_up = seq.upper()
        if clean_pep in seq_up:
            exact_hits.add(pid)
        if il_pep in seq_up.replace("I", "L"):
            il_hits.add(pid)

    il_only_collisions = il_hits - exact_hits
    in_target = bool(target_protein_ids & exact_hits)

    matched_mags: Set[str] = set()
    for pid in il_hits:
        matched_mags.update(protein_to_mags.get(pid, {f"unbinned_or_external::{pid}"}))

    matched_families: Set[str] = set()
    if protein_to_families is not None:
        for pid in il_hits:
            fset = protein_to_families.get(pid)
            if fset:
                matched_families.update(fset)
            else:
                matched_families.add(f"unassigned_paralog_or_background::{pid}")

    single_family_ok = (len(matched_families) == 1) if (protein_to_families is not None) else (len(il_hits) == 1)

    if not in_target:
        status = "rejected_not_in_target_protein"
        unambiguous_for_mag = False
        unambiguous_for_family = False
    elif il_only_collisions:
        collision_mags: Set[str] = set()
        for pid in il_only_collisions:
            collision_mags.update(protein_to_mags.get(pid, {f"unbinned_or_external::{pid}"}))
        target_mags: Set[str] = set()
        for pid in (target_protein_ids & exact_hits):
            target_mags.update(protein_to_mags.get(pid, {f"unbinned_or_external::{pid}"}))
        if collision_mags - target_mags or (il_only_collisions - target_protein_ids):
            status = "ambiguous_isobaric_IL_collision_with_background"
            unambiguous_for_mag = False
            unambiguous_for_family = False
        else:
            unambiguous_for_mag = len(matched_mags) == 1
            unambiguous_for_family = single_family_ok
            status = (
                "unambiguous_unique_single_mag"
                if (unambiguous_for_mag and unambiguous_for_family)
                else "ambiguous_shared_across_multiple_paralog_families_within_mag"
            )
    elif len(il_hits) == 1 and len(matched_mags) == 1:
        status = "unambiguous_unique_single_orf_single_mag"
        unambiguous_for_mag = True
        unambiguous_for_family = single_family_ok
    elif len(matched_mags) == 1:
        unambiguous_for_mag = True
        unambiguous_for_family = single_family_ok
        if single_family_ok:
            status = "unambiguous_unique_single_mag_duplicate_or_same_family_paralog_orf"
        else:
            status = "ambiguous_shared_across_multiple_paralog_families_within_mag"
    elif len(matched_mags) > 1 and single_family_ok:
        status = "shared_within_single_homolog_family_multiple_mags_cannot_assign_single_mag"
        unambiguous_for_mag = False
        unambiguous_for_family = True
    else:
        status = "ambiguous_shared_across_multiple_mags_or_families"
        unambiguous_for_mag = False
        unambiguous_for_family = False

    return {
        "clean_peptide": clean_pep,
        "il_normalized_peptide": il_pep,
        "exact_matched_proteins": sorted(exact_hits),
        "il_matched_proteins": sorted(il_hits),
        "isobaric_il_collision_proteins": sorted(il_only_collisions),
        "has_isobaric_il_collision": len(il_only_collisions) > 0,
        "matched_mags": sorted(matched_mags),
        "matched_families": sorted(matched_families),
        "specificity_status": status,
        "is_unambiguous_for_single_mag": unambiguous_for_mag,
        "is_unambiguous_for_family": unambiguous_for_family,
    }


def parse_mztab_score_metadata(mtd_lines: Iterable[str]) -> Dict[str, str]:
    """
    Parse MTD lines of an MS-GF+ mzTab file to map CV terms to exact search_engine_score[k] columns
    (ROOT-PSM-SCORES.txt):
      - MS:1002055 / MS-GF:PepQValue -> e.g. 'search_engine_score[3]'
      - MS:1002054 / MS-GF:QValue    -> e.g. 'search_engine_score[6]'
      - MS:1002052 / MS-GF:SpecEValue -> e.g. 'search_engine_score[5]'
    Raises ValueError if any required score column is absent from MTD metadata.
    """
    score_cols: Dict[str, str] = {}
    for raw in mtd_lines:
        line = raw.strip()
        if not line.startswith("MTD\t"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        key = parts[1].strip()
        val = parts[2].strip()
        m = re.match(r"^psm_(search_engine_score\[\d+\])$", key)
        if not m:
            continue
        col_name = m.group(1)
        if "MS:1002055" in val or "MS-GF:PepQValue" in val:
            score_cols["pep_q_col"] = col_name
        elif "MS:1002054" in val or "MS-GF:QValue" in val:
            score_cols["psm_q_col"] = col_name
        elif "MS:1002052" in val or "MS-GF:SpecEValue" in val:
            score_cols["spec_evalue_col"] = col_name

    missing = [k for k in ("pep_q_col", "psm_q_col", "spec_evalue_col") if k not in score_cols]
    if missing:
        raise ValueError(
            f"mzTab MTD metadata missing required MS-GF+ score mappings {missing}; found={score_cols}. "
            "Cannot enter primary peptide-q-filtered analysis without upstream MS-GF:PepQValue and MS-GF:QValue."
        )
    return score_cols


def _parse_bounded_probability(raw_val: str, field_name: str, row_idx: int) -> float:
    if not raw_val or raw_val == "null":
        raise ValueError(f"PSM row {row_idx} missing required '{field_name}'")
    try:
        val = float(raw_val)
    except ValueError as exc:
        raise ValueError(f"PSM row {row_idx} has non-numeric {field_name}='{raw_val}'") from exc
    if not math.isfinite(val) or val < 0.0 or val > 1.0:
        raise ValueError(f"PSM row {row_idx} has out-of-bounds [0,1] {field_name}={val}")
    return val


def filter_and_audit_mztab_psms(
    psm_rows: List[Dict[str, str]],
    mtd_lines: Optional[Iterable[str]] = None,
    max_psm_q: float = 0.01,
    max_peptide_q: float = 0.01,
    q_discordance_tol: float = 1e-5,
) -> Dict[str, Any]:
    """
    Parse and filter MS-GF+ mzTab PSM records using metadata-verified upstream MS-GF+ per-hit
    scores (ROOT-PSM-SCORES.txt & ROOT-ANALYSIS-CONTRACT.txt Rule 3).
    """
    if not psm_rows:
        raise ValueError("Empty PSM row list provided to filter_and_audit_mztab_psms")

    if mtd_lines is not None:
        score_cols = parse_mztab_score_metadata(mtd_lines)
    else:
        score_cols = {
            "pep_q_col": "search_engine_score[3]",
            "spec_evalue_col": "search_engine_score[5]",
            "psm_q_col": "search_engine_score[6]",
        }

    pep_q_col = score_cols["pep_q_col"]
    psm_q_col = score_cols["psm_q_col"]
    spec_col = score_cols["spec_evalue_col"]

    qualifying_target_psms = []
    excluded_invalid_rows = []
    n_total = 0
    n_target = 0
    n_decoy = 0
    n_discordant_q = 0

    for idx, row in enumerate(psm_rows):
        n_total += 1
        acc = (row.get("accession") or "").strip()
        if not acc or acc == "null":
            raise ValueError(f"PSM row {idx} missing accession")

        opt_valid = (row.get("opt_global_valid") or "VALID").strip()
        if opt_valid != "VALID":
            excluded_invalid_rows.append({
                "row_idx": idx,
                "psm_id": row.get("PSM_ID", str(idx)),
                "accession": acc,
                "reason": f"opt_global_valid={opt_valid}:{row.get('opt_global_invalid_reason', '')}",
            })
            continue

        raw_opt_q = (row.get("opt_global_cv_MS:1002354_PSM-level_q-value") or "").strip()
        raw_eng_q = (row.get(psm_q_col) or "").strip()
        if not raw_opt_q or raw_opt_q == "null":
            raise ValueError(
                f"PSM row {idx} (accession={acc}) missing per-hit 'opt_global_cv_MS:1002354_PSM-level_q-value'"
            )
        if not raw_eng_q or raw_eng_q == "null":
            raise ValueError(
                f"PSM row {idx} (accession={acc}) missing per-hit '{psm_q_col}' (MS-GF:QValue)"
            )

        opt_psm_q = _parse_bounded_probability(raw_opt_q, "opt_global_cv_MS:1002354_PSM-level_q-value", idx)
        eng_psm_q = _parse_bounded_probability(raw_eng_q, f"{psm_q_col} (MS-GF:QValue)", idx)

        if abs(opt_psm_q - eng_psm_q) > q_discordance_tol:
            n_discordant_q += 1
            raise ValueError(
                f"PSM row {idx} (accession={acc}) has discordant PSM q-values: "
                f"opt_global_cv_MS:1002354_PSM-level_q-value={opt_psm_q} vs {psm_q_col}={eng_psm_q}"
            )

        raw_pep_q = (row.get(pep_q_col) or "").strip()
        if not raw_pep_q or raw_pep_q == "null":
            raise ValueError(
                f"PSM row {idx} (accession={acc}) missing upstream '{pep_q_col}' (MS-GF:PepQValue); "
                "cannot substitute reconstructed peptide q-value or dataset-global FDR."
            )
        pep_q = _parse_bounded_probability(raw_pep_q, f"{pep_q_col} (MS-GF:PepQValue)", idx)

        raw_spec = (row.get(spec_col) or "").strip()
        if not raw_spec or raw_spec == "null":
            raise ValueError(f"PSM row {idx} (accession={acc}) missing '{spec_col}' (MS-GF:SpecEValue)")
        try:
            spec_eval = float(raw_spec)
        except ValueError as exc:
            raise ValueError(f"PSM row {idx} has malformed {spec_col}='{raw_spec}'") from exc
        if not math.isfinite(spec_eval) or spec_eval <= 0.0:
            raise ValueError(f"PSM row {idx} has non-positive or non-finite {spec_col}={spec_eval}")

        decoy_flag = (row.get("opt_global_cv_MS:1002217_decoy_peptide") or "").strip()
        is_xxx = acc.startswith("XXX_")
        if decoy_flag not in ("0", "1"):
            raise ValueError(f"PSM row {idx} (accession={acc}) has invalid decoy flag '{decoy_flag}'")
        if is_xxx and decoy_flag == "0":
            raise ValueError(
                f"Inconsistent decoy annotation at row {idx}: accession '{acc}' has XXX_ prefix but decoy_peptide=0"
            )
        is_decoy = is_xxx or (decoy_flag == "1")
        if is_decoy:
            n_decoy += 1
        else:
            n_target += 1

        raw_seq = (row.get("sequence") or "").strip()
        clean_pep = clean_peptide_sequence(raw_seq)
        il_pep = clean_pep.replace("I", "L")
        psm_id_str = (row.get("PSM_ID") or str(idx)).strip()
        spectra_ref = (row.get("spectra_ref") or "").strip()
        if not spectra_ref or spectra_ref == "null":
            spectra_ref = f"psm_id::{psm_id_str}"

        rec = {
            "psm_id": psm_id_str,
            "spectra_ref": spectra_ref,
            "accession": acc,
            "clean_peptide": clean_pep,
            "il_peptide": il_pep,
            "is_decoy": is_decoy,
            "psm_q_value": eng_psm_q,
            "opt_psm_q_value": opt_psm_q,
            "peptide_q_value": pep_q,
            "spec_evalue": spec_eval,
            "opt_global_valid": opt_valid,
            "opt_global_pass_threshold": (row.get("opt_global_pass_threshold") or "").strip(),
        }

        if (not is_decoy) and eng_psm_q <= max_psm_q and pep_q <= max_peptide_q:
            qualifying_target_psms.append(rec)

    return {
        "score_columns_verified": score_cols,
        "n_total_psms": n_total,
        "n_excluded_invalid_rows": len(excluded_invalid_rows),
        "excluded_invalid_rows": excluded_invalid_rows,
        "n_target_psms": n_target,
        "n_decoy_psms": n_decoy,
        "n_discordant_psm_q": n_discordant_q,
        "n_qualifying_target_psms": len(qualifying_target_psms),
        "qualifying_target_psms": qualifying_target_psms,
    }


def aggregate_unique_spectra_for_entity(
    qualifying_psms: Iterable[Dict[str, Any]],
    target_accessions: Set[str],
    allowed_clean_peptides: Set[str],
) -> Set[str]:
    """
    Union unique spectrum identifiers (spectra_ref) across qualifying PSM rows matching
    target_accessions AND whose clean_peptide is in allowed_clean_peptides
    (ROOT-BACKGROUND-REVIEW.txt line 9).
    - Prevents double-counting when one spectrum matches multiple duplicate aliases or paralogs.
    - Prevents ambiguous peptides (not in allowed_clean_peptides for the entity's specificity level)
      from driving a family- or clade-specific warming abundance.
    """
    unique_spectra: Set[str] = set()
    for psm in qualifying_psms:
        if psm["accession"] in target_accessions and psm["clean_peptide"] in allowed_clean_peptides:
            spec_key = psm.get("spectra_ref") or f"psm_id::{psm['psm_id']}"
            unique_spectra.add(spec_key)
    return unique_spectra


def classify_protein_and_family_detection(
    has_karaoz_homolog: bool,
    in_search_database: bool,
    n_unambiguous_il_peptides_q01: int,
    n_family_shared_il_peptides_q01: int,
    reported_spectral_count: int,
    full_background_verified: bool = False,
    search_database_exclusion_verified: Optional[bool] = None,
) -> str:
    """
    Classify protein/family metaproteomic detection tier according to ROOT-ANALYSIS-CONTRACT.txt
    Rule 3, ROOT-BACKGROUND-REVIEW.txt line 7, and ROOT-FINAL-CLAIM-REVIEW.txt line 2:
    - Absent from study MAGs -> 'nonestimable_no_homolog_in_study_genomes'
    - Not matched in observed EMSL50964 subset:
      * If search_database_exclusion_verified is False (default when full_background_verified=False,
        because the complete 5%-FDR composite search FASTA is not publicly deposited):
        returns 'nonestimable_search_database_inclusion_unverified' (absence from observed EMSL50964
        does NOT prove absence from the searched FASTA).
      * Only when search_database_exclusion_verified is True (explicit complete searched FASTA proves
        exclusion) does it return 'nonestimable_absent_from_search_database'.
    - If full_background_verified is False:
      * NEVER returns 'confirmed_detected_ge2_unambiguous_peptides_q01'!
      * Instead returns 'provisional_observed_subset_ge2_peptides_q01_full_background_unverified_primary_uniqueness_nonestimable'
        when n_unambiguous_il_peptides_q01 >= 2.
    - Only when full_background_verified is True AND n_unambiguous_il_peptides_q01 >= 2 does it
      return 'confirmed_detected_ge2_unambiguous_peptides_q01'.
    """
    if search_database_exclusion_verified is None:
        search_database_exclusion_verified = False
    if not has_karaoz_homolog:
        return "nonestimable_no_homolog_in_study_genomes"
    if not in_search_database:
        if search_database_exclusion_verified:
            return "nonestimable_absent_from_search_database"
        return "nonestimable_search_database_inclusion_unverified"
    if n_unambiguous_il_peptides_q01 >= 2:
        if full_background_verified:
            return "confirmed_detected_ge2_unambiguous_peptides_q01"
        return "provisional_observed_subset_ge2_peptides_q01_full_background_unverified_primary_uniqueness_nonestimable"
    if n_unambiguous_il_peptides_q01 == 1:
        if full_background_verified:
            return "provisional_single_unambiguous_peptide_q01"
        return "provisional_observed_subset_single_peptide_q01_full_background_unverified"
    if n_family_shared_il_peptides_q01 >= 2:
        if full_background_verified:
            return "clade_or_family_level_detected_shared_peptides_cannot_assign_single_mag"
        return "provisional_clade_or_family_observed_subset_ge2_shared_peptides_q01_full_background_unverified"
    if n_family_shared_il_peptides_q01 == 1:
        return "provisional_clade_or_family_single_shared_peptide"
    if reported_spectral_count > 0:
        return "rejected_spectral_count_without_qualifying_q01_unambiguous_peptide"
    return "undetected_in_soluble_metaproteome_not_proof_of_absence"


def classify_primary_warming_estimability(
    has_karaoz_homolog: bool,
    in_search_database: bool,
    has_qualifying_unambiguous_peptides_in_observed_bg: bool,
    reported_spectral_count: int,
    full_background_verified: bool = False,
    search_database_exclusion_verified: Optional[bool] = None,
) -> str:
    """
    Classify primary family/clade warming estimability status (ROOT-FINAL-CLAIM-REVIEW.txt line 3):
    - When full_background_verified is False, even if a family or clade has observed-subset
      unambiguous peptides, primary family/clade warming remains
      'nonestimable_full_search_background_unavailable', and any warming statistics computed
      on the observed subset are retained strictly as 'observed_subset_provisional_sensitivity'.
    """
    if search_database_exclusion_verified is None:
        search_database_exclusion_verified = False
    if not has_karaoz_homolog:
        return "nonestimable_no_homolog_in_study_genomes"
    if not in_search_database:
        if search_database_exclusion_verified:
            return "nonestimable_absent_from_search_database"
        return "nonestimable_search_database_inclusion_unverified"
    if not has_qualifying_unambiguous_peptides_in_observed_bg:
        if reported_spectral_count > 0:
            return "nonestimable_rejected_spectral_count_without_qualifying_q01_unambiguous_peptide"
        return "nonestimable_undetected_in_soluble_metaproteome_not_proof_of_absence"
    if not full_background_verified:
        return "nonestimable_full_search_background_unavailable"
    return "estimable_confirmed_full_search_background"


def validate_karaoz_experimental_design(
    sample_records: List[Dict[str, Any]],
    claimed_independent_n: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Validate the Karaoz et al. 2026 soil warming experimental design (ROOT-ANALYSIS-CONTRACT.txt Rules 4 & 5).
    """
    if not sample_records:
        raise ValueError("Empty sample_records in validate_karaoz_experimental_design")

    valid_plots = {"1", "2", "3"}
    valid_treatments = {"control", "heated"}
    valid_depths = {"0-10", "10-30", "30-45", "45-60", "60-80"}

    seen_cells: Set[Tuple[str, str, str]] = set()
    seen_sample_ids: Set[str] = set()
    plot_pairs_seen: Set[str] = set()

    for idx, rec in enumerate(sample_records):
        sid = str(rec.get("sample_id", f"row_{idx}")).strip()
        if sid in seen_sample_ids:
            raise ValueError(f"Duplicate sample_id '{sid}' in experimental design")
        seen_sample_ids.add(sid)

        year = int(rec.get("sampling_year", 0))
        dna_year = int(rec.get("paired_dna_year", year))
        if year != 2018 or dna_year != 2018:
            raise ValueError(
                f"Sample '{sid}' has unmatched or non-2018 sampling year (proteome_year={year}, dna_year={dna_year}); "
                "2014/2016 metagenomes must not be pooled with 2018 metaproteomes."
            )

        plot_pair = str(rec.get("plot_pair_id", "")).strip()
        if plot_pair not in valid_plots:
            raise ValueError(f"Invalid plot_pair_id '{plot_pair}' for sample '{sid}'; expected one of {sorted(valid_plots)}")
        plot_pairs_seen.add(plot_pair)

        trt = str(rec.get("treatment", "")).strip().lower()
        if trt not in valid_treatments:
            raise ValueError(f"Invalid treatment '{trt}' for sample '{sid}'")

        depth = str(rec.get("depth_cm", "")).strip()
        if depth not in valid_depths:
            raise ValueError(f"Invalid depth_cm '{depth}' for sample '{sid}'")

        cell = (plot_pair, trt, depth)
        if cell in seen_cells:
            raise ValueError(
                f"Duplicate biological/technical cell {cell} for sample '{sid}'; "
                "cannot duplicate biological replicates."
            )
        seen_cells.add(cell)

    n_true_plot_pairs = len(plot_pairs_seen)
    if claimed_independent_n is not None and claimed_independent_n > n_true_plot_pairs:
        raise ValueError(
            f"Pseudoreplication error: claimed_independent_n={claimed_independent_n} exceeds true "
            f"biological plot-pair replication n={n_true_plot_pairs} (depths are nested repeated measures)."
        )

    return {
        "n_samples": len(sample_records),
        "n_independent_plot_pairs": n_true_plot_pairs,
        "n_nested_depth_strata": len({r["depth_cm"] for r in sample_records}),
        "max_whole_plot_sign_flips": 2 ** n_true_plot_pairs,
        "min_possible_two_sided_sign_flip_p": 2.0 / (2 ** n_true_plot_pairs) if n_true_plot_pairs > 0 else 1.0,
    }


def compute_paired_warming_effect_n3_plots(
    sample_records: List[Dict[str, Any]],
    value_field: str = "value",
    claimed_independent_n: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Compute paired heated-vs-control warming differences respecting the n=3 paired-plot block structure
    (ROOT-ANALYSIS-CONTRACT.txt Rule 4 and ROOT-FINAL-CLAIM-REVIEW.txt line 4):
    - Uses paired differences (heated - control) as the primary paired quantity.
    - NEVER uses arbitrary 1e-6 pseudocounts for log2 fold changes at zero values (which produce
      meaningless ~24-30 artifacts).
    - Descriptive log2(heated/control) is reported ONLY for individual pairs where both heated > 0
      and control > 0; zero-involved pairs are None ('NA') with an explicit detection_state
      ('both_positive', 'heated_only_zero_control', 'control_only_zero_heated', 'both_zero'),
      and no mean log2FC over changing non-zero subsets is computed or presented as a full-plot estimator.
    """
    design = validate_karaoz_experimental_design(sample_records, claimed_independent_n=claimed_independent_n)

    grid: Dict[Tuple[str, str], Dict[str, float]] = defaultdict(dict)
    for rec in sample_records:
        p = str(rec["plot_pair_id"]).strip()
        d = str(rec["depth_cm"]).strip()
        t = str(rec["treatment"]).strip().lower()
        val = float(rec[value_field])
        if not math.isfinite(val) or val < 0.0:
            raise ValueError(f"Invalid abundance value {val} in sample {rec.get('sample_id')}")
        grid[(p, d)][t] = val

    depth_deltas_by_plot: Dict[str, List[float]] = defaultdict(list)
    depth_profile: Dict[str, Dict[str, Any]] = defaultdict(dict)
    n_heated_pos = 0
    n_control_pos = 0
    n_zero_samples = 0
    pair_state_counts: Counter = Counter()

    for (p, d), trt_map in sorted(grid.items()):
        if "heated" not in trt_map or "control" not in trt_map:
            raise ValueError(f"Unpaired plot-depth cell ({p}, {d}): present={sorted(trt_map.keys())}")
        y_h = trt_map["heated"]
        y_c = trt_map["control"]
        if y_h > 0.0:
            n_heated_pos += 1
        else:
            n_zero_samples += 1
        if y_c > 0.0:
            n_control_pos += 1
        else:
            n_zero_samples += 1

        if y_h > 0.0 and y_c > 0.0:
            det_state = "both_positive"
            desc_l2r: Optional[float] = math.log2(y_h / y_c)
        elif y_h > 0.0 and y_c == 0.0:
            det_state = "heated_only_zero_control"
            desc_l2r = None
        elif y_h == 0.0 and y_c > 0.0:
            det_state = "control_only_zero_heated"
            desc_l2r = None
        else:
            det_state = "both_zero"
            desc_l2r = None

        pair_state_counts[det_state] += 1
        delta = y_h - y_c
        depth_deltas_by_plot[p].append(delta)
        depth_profile[d][f"plot{p}_heated"] = y_h
        depth_profile[d][f"plot{p}_control"] = y_c
        depth_profile[d][f"plot{p}_delta"] = delta
        depth_profile[d][f"plot{p}_detection_state"] = det_state
        depth_profile[d][f"plot{p}_descriptive_log2_ratio_both_positive_only"] = desc_l2r

    plot_ids = sorted(depth_deltas_by_plot.keys())
    plot_mean_deltas = {p: sum(depth_deltas_by_plot[p]) / len(depth_deltas_by_plot[p]) for p in plot_ids}

    deltas_vec = [plot_mean_deltas[p] for p in plot_ids]
    n_pairs = len(deltas_vec)
    obs_mean_delta = sum(deltas_vec) / float(n_pairs)

    n_pos = sum(1 for x in deltas_vec if x > 0)
    n_neg = sum(1 for x in deltas_vec if x < 0)
    n_zero = sum(1 for x in deltas_vec if x == 0)

    if all(x == 0.0 for x in deltas_vec):
        exact_p = 1.0
    else:
        n_extreme = 0
        total_flips = 0
        for signs in itertools.product((-1.0, 1.0), repeat=n_pairs):
            total_flips += 1
            perm_mean = sum(s * x for s, x in zip(signs, deltas_vec)) / float(n_pairs)
            if abs(perm_mean) >= abs(obs_mean_delta) - 1e-12:
                n_extreme += 1
        exact_p = n_extreme / float(total_flips)

    return {
        "n_independent_plot_pairs": n_pairs,
        "n_nested_depth_pairs": len(grid),
        "n_heated_positive_samples_of_15": n_heated_pos,
        "n_control_positive_samples_of_15": n_control_pos,
        "n_zero_samples_of_30": n_zero_samples,
        "n_depth_pairs_both_positive": pair_state_counts["both_positive"],
        "n_depth_pairs_heated_only_zero_control": pair_state_counts["heated_only_zero_control"],
        "n_depth_pairs_control_only_zero_heated": pair_state_counts["control_only_zero_heated"],
        "n_depth_pairs_both_zero": pair_state_counts["both_zero"],
        "plot_mean_deltas": plot_mean_deltas,
        "overall_mean_delta_heated_minus_control": obs_mean_delta,
        "n_plot_pairs_positive": n_pos,
        "n_plot_pairs_negative": n_neg,
        "n_plot_pairs_zero": n_zero,
        "replicate_sign_consistency": f"{max(n_pos, n_neg)}/{n_pairs}",
        "exact_whole_plot_sign_flip_two_sided_p": exact_p,
        "min_achievable_two_sided_p_for_n3": design["min_possible_two_sided_sign_flip_p"],
        "depth_profile": dict(depth_profile),
    }


def benjamini_hochberg_fdr(p_values: Dict[str, Optional[float]]) -> Dict[str, Optional[float]]:
    """
    Compute Benjamini-Hochberg adjusted q-values across a cohort of tested items,
    preserving None for non-estimable items (so non-estimable families never inflate or distort m).
    """
    valid = [(k, float(v)) for k, v in p_values.items() if v is not None and math.isfinite(float(v))]
    m = len(valid)
    out: Dict[str, Optional[float]] = {k: None for k in p_values}
    if m == 0:
        return out
    valid.sort(key=lambda kv: kv[1])
    raw_q = []
    for rank_idx, (k, pval) in enumerate(valid, start=1):
        raw_q.append((k, min(1.0, (pval * m) / float(rank_idx))))
    running_min = 1.0
    for k, qval in reversed(raw_q):
        running_min = min(running_min, qval)
        out[k] = running_min
    return out


def parse_and_validate_flank_evidence_rows(
    rows: Iterable[Dict[str, str]],
    primary_leads: Dict[str, str],
    partner_families: Dict[str, str],
    require_all_leads_and_partners: bool = True,
) -> Dict[str, Any]:
    """
    Parse and validate rows from exp4_47_flank_evidence.tsv using its exact real headers:
    - family_id
    - focal_member_gene
    - focal_genome
    - flanking_gene
    - flanking_family_id
    - annotated_guild_support
    (ROOT-RUNTIME-REPAIR.txt).

    Raises KeyError if any required column is absent, and ValueError if require_all_leads_and_partners=True
    and any primary lead or partner family has 0 mapped member genes after join.
    """
    required_cols = (
        "family_id",
        "focal_member_gene",
        "focal_genome",
        "flanking_gene",
        "flanking_family_id",
    )
    gene_to_fam_role: Dict[str, Tuple[str, str]] = {}
    needed_from_genome_faa: Dict[str, List[Tuple[str, str, str]]] = defaultdict(list)
    found_leads: Set[str] = set()
    found_partners: Set[str] = set()

    for idx, row in enumerate(rows):
        for col in required_cols:
            if col not in row or row[col] is None:
                raise KeyError(
                    f"Missing required column '{col}' in exp4_47_flank_evidence.tsv row {idx}; "
                    f"available columns={sorted(row.keys())}"
                )
        focal_fid = row["family_id"].strip()
        focal_gid = row["focal_member_gene"].strip()
        gname = row["focal_genome"].strip()
        flank_gid = row["flanking_gene"].strip()
        flank_fid = row["flanking_family_id"].strip()

        if not focal_fid or not focal_gid or not gname or not flank_gid or not flank_fid:
            raise ValueError(f"Empty required field in exp4_47_flank_evidence.tsv row {idx}: {row}")

        if flank_fid in partner_families:
            rtag = partner_families[flank_fid]
            gene_to_fam_role[flank_gid] = (flank_fid, rtag)
            needed_from_genome_faa[gname].append((flank_gid, flank_fid, rtag))
            found_partners.add(flank_fid)

        if focal_fid in primary_leads:
            rtag = primary_leads[focal_fid]
            gene_to_fam_role[focal_gid] = (focal_fid, rtag)
            needed_from_genome_faa[gname].append((focal_gid, focal_fid, rtag))
            found_leads.add(focal_fid)

    if require_all_leads_and_partners:
        missing_leads = sorted(set(primary_leads.keys()) - found_leads)
        missing_partners = sorted(set(partner_families.keys()) - found_partners)
        if missing_leads or missing_partners:
            raise ValueError(
                f"Flank evidence join incomplete: missing_leads={missing_leads}, "
                f"missing_partners={missing_partners}"
            )

    return {
        "gene_to_fam_role": gene_to_fam_role,
        "needed_from_genome_faa": dict(needed_from_genome_faa),
        "found_leads": sorted(found_leads),
        "found_partners": sorted(found_partners),
    }

