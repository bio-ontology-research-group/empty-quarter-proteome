#!/usr/bin/env python3
"""Core pure scientific functions and provenance validators for exp3_viruses.

Implements strict rules from:
- `mechanism-20261109/EXP3-CONTRACT.txt`
- `mechanism-20261109/EXP3-STAGING-NOTE.txt`
- `mechanism-20261109/EXP3-EVIDENCE-REVIEW.txt`
- `mechanism-20261109/EXP3-CODE-REVIEW.txt`
- `coordinated-revision-20261108/mechanism-experiments-20261109.txt`

Key scientific principles enforced:
1. Non-circular canonical 591 archaeal MAG FASTA verification:
   - Primary canonical archive verification strictly enforces exact SHA-256,
     contig count, and base count against `genome_manifest.tsv` with
     `verified_canonical_reference_bytes=None`.
   - Only after a canonical reference FASTA is independently verified against
     `expected_sha256` may its bytes be passed as `verified_canonical_reference_bytes`
     to evaluate an alternative product or registered header transform.
2. Exact namespaced `sample::contig` identifiers across MAGs, assemblies,
   proviruses, and cellular contig indexes so identical bare contig suffixes
   in different samples never collide.
3. Strict host ID reconciliation (no heuristic suffix-stripping joins such as
   stripping `_sub` or guessing `self_`/`single_`).
4. Keep CRISPR match evidence strictly separate from viral identity:
   - Mask spacer-array self-matches (`Technical_self_crispr_array_masked`) on
     the same namespaced `sample::contig` locus as technical matches.
   - Treat `gapopen` as gap EVENT count (not gap residue count): require total
     `gaps` (gap residues) for numerical gapped mismatch calculations or mark
     gapped alignments without residue counts as unsupported; reject gapped
     alignments from strict/exploratory CRISPR spacer gates.
5. Independent viral identity vs unknown vs cellular-like vs curated controls:
   - Missing annotations (`None`/`NA`) are never treated as zero (`0`).
   - Overlapping evidence across tools (e.g. Foldseek structural hits, geNomad
     hallmarks, archaeal HMMs, and CheckV viral genes hitting the same ORF) is
     never summed as independent gene counts; distinct ORF union count or
     single-source maximum is used instead.
   - `CheckV Not-determined` + no recognised hallmark is `Unknown_unverified_viral_identity`
     (NEVER `Verified_nonviral` and NEVER turned into a negative control).
   - External reference viruses from other habitats (`decoys_derep.fna`) are STILL
     viral sequences (`Verified_viral_external_habitat_reference`), not negative
     viral-identity controls; their role as sample/mapping specificity controls
     is recorded separately (`habitat_specificity_role`), and hosts are never
     inferred wrong simply because a reference virus is not from EQ.
   - Only deliberately selected verified non-viral plasmid or cellular chromosomal
     reference controls receive `Decoy_curated_*` non-viral labels.
"""
from __future__ import annotations

import hashlib
import math
import re
from pathlib import Path
from typing import Any

VALID_DNA_SET = set("ACGT")
VALID_IUPAC_SET = set("ACGTN")
HOMOPOLYMER_RE = re.compile(r"A{7,}|C{7,}|G{7,}|T{7,}")


def sha256_bytes(data: bytes) -> str:
    """Return hex SHA-256 digest of raw bytes."""
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path | str) -> str:
    """Return hex SHA-256 digest of a file on disk."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _opt_int(val: Any) -> int | None:
    """Parse optional integer without conflating missing/NA with 0."""
    if val is None or isinstance(val, bool):
        return None
    if isinstance(val, int):
        return val
    if isinstance(val, float):
        if math.isnan(val):
            return None
        return int(val)
    s = str(val).strip()
    if s in ("", "NA", "NaN", "nan", "None", "null", "unscored", "missing"):
        return None
    return int(float(s))


def _opt_float(val: Any) -> float | None:
    """Parse optional float without conflating missing/NA with 0.0."""
    if val is None or isinstance(val, bool):
        return None
    if isinstance(val, (int, float)):
        f = float(val)
        return None if math.isnan(f) else f
    s = str(val).strip()
    if s in ("", "NA", "NaN", "nan", "None", "null", "unscored", "missing"):
        return None
    return float(s)


def parse_fasta_bytes(raw_bytes: bytes) -> tuple[list[tuple[str, str]], list[str]]:
    """Parse FASTA bytes into list of (header, uppercase_seq) and list of duplicate raw IDs."""
    records: list[tuple[str, str]] = []
    seen_ids: set[str] = set()
    dup_ids: list[str] = []
    cur_header: str | None = None
    cur_seq: list[str] = []

    text = raw_bytes.decode("utf-8", errors="replace")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_header is not None:
                records.append((cur_header, "".join(cur_seq).upper()))
            cur_header = line[1:].strip()
            raw_id = cur_header.split()[0] if cur_header else ""
            if raw_id in seen_ids:
                dup_ids.append(raw_id)
            seen_ids.add(raw_id)
            cur_seq = []
        else:
            cur_seq.append(line)
    if cur_header is not None:
        records.append((cur_header, "".join(cur_seq).upper()))
    return records, dup_ids


def normalize_contig_header(header: str, transform: str, source_sample: str = "") -> str:
    """Apply a registered deterministic header transform to a FASTA header.

    When `source_sample` is provided and `transform == "strip_sample_double_colon"`,
    only strips the exact `<source_sample>::` prefix (or preserves explicit sample
    namespace so distinct samples never collapse).
    """
    raw_id = header.split()[0] if header else ""
    if transform == "raw_id":
        return raw_id
    if transform == "strip_sample_double_colon":
        if "::" in raw_id:
            samp, rest = raw_id.split("::", 1)
            if source_sample and samp != source_sample:
                return raw_id
            return rest
        return raw_id
    raise ValueError(f"Unregistered header transform: {transform}")


def verify_and_reconcile_fasta(
    mag_id: str,
    raw_bytes: bytes,
    expected_sha256: str,
    expected_contigs: int,
    expected_bases: int,
    verified_canonical_reference_bytes: bytes | None = None,
) -> dict[str, Any]:
    """Verify a genome FASTA against canonical manifest expectations without circularity.

    1. Primary canonical archive verification (`verified_canonical_reference_bytes=None`):
       Requires exact byte SHA-256 match (`comp_sha256 == expected_sha256`), exact
       `n_contigs == expected_contigs`, and exact `n_bases == expected_bases`.
       Any sequence mutation or SHA-256 mismatch is strictly rejected.
    2. Secondary/alternative product verification (`verified_canonical_reference_bytes` provided):
       First independently verifies that `verified_canonical_reference_bytes` itself
       matches `expected_sha256`, `expected_contigs`, and `expected_bases` (raising
       `ValueError` if not). Only then checks whether `raw_bytes` has identical
       contig sequences and headers under a registered header transform.
    """
    comp_sha256 = sha256_bytes(raw_bytes)
    records, dup_ids = parse_fasta_bytes(raw_bytes)
    n_contigs = len(records)
    n_bases = sum(len(seq) for _, seq in records)

    if dup_ids:
        return {
            "mag": mag_id,
            "computed_sha256": comp_sha256,
            "expected_sha256": expected_sha256,
            "n_contigs": n_contigs,
            "expected_contigs": int(expected_contigs),
            "n_bases": n_bases,
            "expected_bases": int(expected_bases),
            "header_transform_applied": "none",
            "status": "rejected_duplicate_contig_ids",
            "is_verified_canonical": False,
        }

    if (
        comp_sha256 == expected_sha256
        and n_contigs == int(expected_contigs)
        and n_bases == int(expected_bases)
    ):
        return {
            "mag": mag_id,
            "computed_sha256": comp_sha256,
            "expected_sha256": expected_sha256,
            "n_contigs": n_contigs,
            "expected_contigs": int(expected_contigs),
            "n_bases": n_bases,
            "expected_bases": int(expected_bases),
            "header_transform_applied": "exact_bytes",
            "status": "verified_exact_match_genome_manifest",
            "is_verified_canonical": True,
        }

    if verified_canonical_reference_bytes is not None:
        ref_sha = sha256_bytes(verified_canonical_reference_bytes)
        ref_records, ref_dups = parse_fasta_bytes(verified_canonical_reference_bytes)
        ref_n_contigs = len(ref_records)
        ref_n_bases = sum(len(s) for _, s in ref_records)
        if (
            ref_dups
            or ref_sha != expected_sha256
            or ref_n_contigs != int(expected_contigs)
            or ref_n_bases != int(expected_bases)
        ):
            raise ValueError(
                f"Cannot reconcile {mag_id}: verified_canonical_reference_bytes failed "
                f"independent manifest SHA-256/count verification "
                f"(ref_sha={ref_sha}, expected_sha256={expected_sha256})"
            )

        if n_contigs == ref_n_contigs and n_bases == ref_n_bases:
            source_sample = mag_id.split("__")[0] if "__" in mag_id else ""
            for transform in ("raw_id", "strip_sample_double_colon"):
                cand_set = sorted(
                    (
                        normalize_contig_header(h, transform, source_sample),
                        sha256_bytes(s.encode("ascii")),
                        len(s),
                    )
                    for h, s in records
                )
                ref_set = sorted(
                    (
                        normalize_contig_header(h, transform, source_sample),
                        sha256_bytes(s.encode("ascii")),
                        len(s),
                    )
                    for h, s in ref_records
                )
                if cand_set == ref_set:
                    return {
                        "mag": mag_id,
                        "computed_sha256": comp_sha256,
                        "expected_sha256": expected_sha256,
                        "n_contigs": n_contigs,
                        "expected_contigs": int(expected_contigs),
                        "n_bases": n_bases,
                        "expected_bases": int(expected_bases),
                        "header_transform_applied": transform,
                        "status": "verified_identical_contigs_registered_header_transform",
                        "is_verified_canonical": True,
                    }

    status = (
        "rejected_canonical_sha256_sequence_mismatch"
        if (n_contigs == int(expected_contigs) and n_bases == int(expected_bases))
        else "rejected_wrong_product_contig_membership_differs"
    )
    return {
        "mag": mag_id,
        "computed_sha256": comp_sha256,
        "expected_sha256": expected_sha256,
        "n_contigs": n_contigs,
        "expected_contigs": int(expected_contigs),
        "n_bases": n_bases,
        "expected_bases": int(expected_bases),
        "header_transform_applied": "none",
        "status": status,
        "is_verified_canonical": False,
    }


def reconcile_host_alias_strict(
    host_id: str,
    canonical_591_set: set[str],
    nonrep_to_species_map: dict[str, str] | None = None,
) -> tuple[str, str | None]:
    """Reconcile a candidate host MAG identifier without heuristic suffix stripping.

    Never strips `_sub`, never guesses `self_` or `single_` insertions, and never
    performs substring matching.
    Returns `(reconciliation_status, canonical_591_mag_id_or_None)`.
    """
    if not host_id or not isinstance(host_id, str):
        return ("rejected_empty_host_id", None)
    clean_id = host_id.strip()
    if clean_id in canonical_591_set:
        return ("canonical_591_direct", clean_id)
    if nonrep_to_species_map and clean_id in nonrep_to_species_map:
        sp_rep = nonrep_to_species_map[clean_id]
        if sp_rep in canonical_591_set:
            return ("nonrep_refined_bin_mapped_to_species_rep", sp_rep)
    return ("rejected_unreconciled_host_alias", None)


def shannon_entropy_dna(seq: str) -> float:
    """Compute Shannon entropy in bits per base over A/C/G/T."""
    s = seq.upper()
    n = len(s)
    if n == 0:
        return 0.0
    counts = [s.count(b) for b in "ACGT"]
    tot = sum(counts)
    if tot == 0:
        return 0.0
    return -sum((c / tot) * math.log2(c / tot) for c in counts if c > 0)


def is_valid_crispr_repeat_and_spacer(
    repeat_seq: str,
    spacer_seq: str,
    min_repeat_len: int = 23,
    max_repeat_len: int = 55,
    min_spacer_len: int = 18,
    max_spacer_len: int = 60,
    min_entropy: float = 1.6,
) -> tuple[bool, str]:
    """Validate a candidate CRISPR repeat and spacer against low-complexity & length gates."""
    rep = (repeat_seq or "").strip().upper()
    spc = (spacer_seq or "").strip().upper()
    if not (min_repeat_len <= len(rep) <= max_repeat_len):
        return (False, "rejected_repeat_length")
    if not (min_spacer_len <= len(spc) <= max_spacer_len):
        return (False, "rejected_spacer_length")
    if not (set(rep) <= VALID_IUPAC_SET and set(spc) <= VALID_IUPAC_SET):
        return (False, "rejected_invalid_characters")
    if spc.count("N") > 0 or rep.count("N") > 1:
        return (False, "rejected_ambiguous_bases")
    if HOMOPOLYMER_RE.search(rep) or HOMOPOLYMER_RE.search(spc):
        return (False, "rejected_homopolymer_run")
    if shannon_entropy_dna(rep) < min_entropy or shannon_entropy_dna(spc) < min_entropy:
        return (False, "rejected_low_complexity_entropy")
    for di in ("AT", "TA", "CG", "GC", "AC", "CA", "AG", "GA", "CT", "TC", "GT", "TG"):
        if spc.count(di) * 2 >= 0.75 * len(spc):
            return (False, "rejected_dinucleotide_tandem_repeat")
    return (True, "valid_crispr_spacer")


def score_spacer_alignment(
    qlen: int,
    qstart: int,
    qend: int,
    sstart: int,
    send: int,
    length: int,
    mismatch: int,
    gapopen: int,
    pident: float,
    gaps: int | None = None,
) -> dict[str, Any]:
    """Explicitly score a BLASTn-short spacer-to-target alignment (Axis 1: CRISPR match quality).

    Accounts for internal mismatches, gaps, and unaligned 5'/3' spacer ends.
    Treats `gapopen` strictly as gap EVENT count (never gap residue count):
    - If an alignment is gapped (`gapopen > 0`, `gaps > 0`, or span != alignment length):
      - When `gaps` (total gap residues) is provided, `effective_mismatches` uses `gaps`.
      - When `gaps` is not provided, `effective_mismatches` and `full_spacer_identity`
        are set to `None` (`unsupported_gapped_alignment_without_gap_residue_count`)
        so inaccurate identity is never claimed.
      - Gapped alignments are rejected from strict and exploratory spacer gates.
    """
    if qlen <= 0:
        raise ValueError(f"Invalid spacer qlen: {qlen}")
    q_lo = min(int(qstart), int(qend))
    q_hi = max(int(qstart), int(qend))
    s_lo = min(int(sstart), int(send))
    s_hi = max(int(sstart), int(send))

    aligned_q_span = max(0, q_hi - q_lo + 1)
    aligned_s_span = max(0, s_hi - s_lo + 1)
    uncovered_5p = max(0, q_lo - 1)
    uncovered_3p = max(0, int(qlen) - q_hi)
    max_end_uncovered = max(uncovered_5p, uncovered_3p)
    full_spacer_coverage = aligned_q_span / float(qlen)

    parsed_gaps = _opt_int(gaps)
    has_gap = (
        int(gapopen) > 0
        or (parsed_gaps is not None and parsed_gaps > 0)
        or aligned_q_span != int(length)
        or aligned_s_span != int(length)
    )

    if has_gap:
        if parsed_gaps is not None:
            total_gap_residues: int | None = parsed_gaps
            eff_mm: int | None = int(mismatch) + parsed_gaps + uncovered_5p + uncovered_3p
            full_id: float | None = round(max(0.0, (float(qlen) - float(eff_mm)) / float(qlen)), 6)
            match_tier = "rejected_gapped_alignment"
            gap_status = "gapped_rejected_using_total_gap_residues"
        else:
            total_gap_residues = None
            eff_mm = None
            full_id = None
            match_tier = "rejected_gapped_alignment_unsupported"
            gap_status = "unsupported_gapped_alignment_gapopen_is_event_count_not_residue_count"
        passes_strict = False
        passes_exploratory = False
    else:
        total_gap_residues = 0
        eff_mm = int(mismatch) + uncovered_5p + uncovered_3p
        full_id = round(max(0.0, (float(qlen) - float(eff_mm)) / float(qlen)), 6)
        gap_status = "ungapped_exact_calculation"

        if (
            eff_mm == 0
            and int(mismatch) == 0
            and aligned_q_span == int(qlen)
            and float(pident) >= 99.9
        ):
            match_tier = "exact_0mm_full"
            passes_strict = True
            passes_exploratory = True
        elif (
            int(mismatch) <= 1
            and uncovered_5p <= 2
            and uncovered_3p <= 2
            and full_spacer_coverage >= 0.95
            and float(pident) >= 95.0
            and eff_mm <= 2
        ):
            match_tier = "strict_1mm_end2"
            passes_strict = True
            passes_exploratory = True
        elif (
            eff_mm <= 2
            and max_end_uncovered <= 3
            and full_spacer_coverage >= 0.90
            and float(pident) >= 90.0
        ):
            match_tier = "near_2mm_exploratory"
            passes_strict = False
            passes_exploratory = True
        else:
            match_tier = "rejected_mismatch_or_end_coverage"
            passes_strict = False
            passes_exploratory = False

    return {
        "qlen": int(qlen),
        "qstart": int(qstart),
        "qend": int(qend),
        "sstart": int(sstart),
        "send": int(send),
        "target_hit_lo": s_lo,
        "target_hit_hi": s_hi,
        "aligned_q_span": aligned_q_span,
        "aligned_s_span": aligned_s_span,
        "uncovered_5p": uncovered_5p,
        "uncovered_3p": uncovered_3p,
        "max_end_uncovered": max_end_uncovered,
        "internal_mismatch": int(mismatch),
        "gapopen": int(gapopen),
        "gap_residues": total_gap_residues,
        "gap_handling_status": gap_status,
        "effective_mismatches": eff_mm,
        "blast_pident": float(pident),
        "full_spacer_coverage": round(full_spacer_coverage, 6),
        "full_spacer_identity": full_id,
        "match_tier": match_tier,
        "passes_strict_spacer_gate": passes_strict,
        "passes_exploratory_spacer_gate": passes_exploratory,
    }


DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
HEADER_SAMPLE_PREFIX_RE = re.compile(r"^([0-9A-Za-z]+r\d+|coasm|N\d+[A-Za-z0-9]*)_\d+")


def revcomp_dna(seq: str) -> str:
    """Return reverse-complement of a DNA sequence."""
    return (seq or "").strip().upper().translate(DNA_COMPLEMENT)[::-1]


def canonical_spacer_seq_rc(seq: str) -> tuple[str, str]:
    """Return `(canonical_seq, orientation)` where `canonical_seq = min(seq, revcomp(seq))`.

    Audits reverse-complement duplication so identical spacers assembled or called
    on opposite strands are not counted twice as independent spacers.
    """
    fwd = (seq or "").strip().upper()
    rc = revcomp_dna(fwd)
    if fwd <= rc:
        return (fwd, "fwd")
    return (rc, "rc")


def parse_target_provirus_coords(target_id: str) -> tuple[str, str, int | None, int | None]:
    """Parse sample component, base contig ID, and excised provirus coordinates from a target ID."""
    sample_part = ""
    rest = target_id.strip()
    if "::" in rest:
        sample_part, rest = rest.split("::", 1)

    m_prov = re.match(r"^(.+?)\|provirus_(\d+)_(\d+)(?:_\d+)?$", rest)
    if m_prov:
        return (sample_part, m_prov.group(1), int(m_prov.group(2)), int(m_prov.group(3)))

    return (sample_part, rest, None, None)


def qualify_namespaced_contig_id(
    raw_id: str,
    fallback_sample: str = "",
) -> tuple[str, str, str, int | None, int | None]:
    """Return `(sample_id, base_contig_id, qualified_sample_contig_id, prov_start, prov_end)`.

    Priority for sample namespace:
    1. Explicit `<sample>::<contig>` in header.
    2. Embedded sample prefix in contig header (`<sample>_<contig_num>`, e.g. `41PRr2_12345`),
       so co-assembled or cross-sample contigs retain their true header sample rather than
       being overwritten by `fallback_sample`.
    3. `fallback_sample` (e.g. from MAG name) ONLY when the contig header has no explicit or
       embedded sample prefix.
    """
    sample_part, base_contig, prov_start, prov_end = parse_target_provirus_coords(raw_id)
    sample_clean = sample_part.strip()
    if not sample_clean and base_contig:
        m_emb = HEADER_SAMPLE_PREFIX_RE.match(base_contig)
        if m_emb:
            sample_clean = m_emb.group(1)
    if not sample_clean:
        sample_clean = (fallback_sample or "").strip()

    if sample_clean:
        qualified = f"{sample_clean}::{base_contig}"
    else:
        qualified = base_contig
    return (sample_clean, base_contig, qualified, prov_start, prov_end)


def are_same_namespaced_contig(
    target_id: str,
    host_contig_id: str,
    host_mag: str = "",
    target_sample_hint: str = "",
    target_contig_sha256: str = "",
    host_contig_sha256: str = "",
) -> tuple[bool, str, str]:
    """Determine whether `target_id` and `host_contig_id` refer to the exact same namespaced contig.

    - Never discards the sample component when comparing `sample::contig` IDs.
    - If both have sample namespaces and `target_sample != host_sample`, returns `False`
      even if `target_base_contig == host_base_contig`.
    - If only one side has a sample namespace, requires either matching sequence SHA-256
      or that the unprefixed contig ID explicitly starts with `<sample>_` matching the
      other side's sample namespace.
    """
    host_sample_from_mag = host_mag.split("__")[0].strip() if "__" in (host_mag or "") else ""
    t_samp, t_base, t_qual, _, _ = qualify_namespaced_contig_id(target_id, fallback_sample=target_sample_hint)
    h_samp, h_base, h_qual, _, _ = qualify_namespaced_contig_id(host_contig_id, fallback_sample=host_sample_from_mag)

    if not t_base or not h_base:
        return (False, t_qual, h_qual)

    if target_contig_sha256 and host_contig_sha256 and target_contig_sha256 not in ("NA", "") and host_contig_sha256 not in ("NA", ""):
        if target_contig_sha256 != host_contig_sha256:
            return (False, t_qual, h_qual)
        if t_samp and h_samp and t_samp != h_samp:
            return (False, t_qual, h_qual)
        if t_base == h_base:
            return (True, t_qual, h_qual)

    if t_samp and h_samp:
        return (t_samp == h_samp and t_base == h_base, t_qual, h_qual)

    if not t_samp and not h_samp:
        return (t_base == h_base, t_qual, h_qual)

    # Exactly one side has an explicit sample namespace: only match if the bare contig
    # carries the exact `<sample>_` prefix of that sample AND base contigs match.
    known_samp = t_samp or h_samp
    if t_base == h_base and t_base.startswith(f"{known_samp}_"):
        return (True, f"{known_samp}::{t_base}", f"{known_samp}::{h_base}")

    return (False, t_qual, h_qual)


def evaluate_viral_identity_and_decoys(
    target_id: str,
    target_hit_lo: int,
    target_hit_hi: int,
    host_mag: str,
    host_contig_id: str,
    host_array_lo: int,
    host_array_hi: int,
    target_meta: dict[str, Any],
    all_crispr_arrays_on_contig: list[tuple[int, int]] | None = None,
) -> dict[str, Any]:
    """Evaluate Axis 2 (viral identity vs unknown vs cellular-like vs curated controls) and technical self-array masks.

    Strictly follows `EXP3-EVIDENCE-REVIEW.txt` and `EXP3-CODE-REVIEW.txt`:
    1. Preserves exact `sample::contig` namespace: never conflates same bare contig names
       from different samples.
    2. Technical CRISPR array self-matches on the same namespaced contig or known CRISPR locus
       on the target contig are masked as `Technical_self_crispr_array_masked`.
    3. External reference viruses from other habitats (`is_external_habitat_viral_control` /
       `is_refseq_non_eq_decoy`) are STILL viral sequences (`Verified_viral_external_habitat_reference`,
       `is_verified_viral=True`, `is_decoy=False`), recorded separately under
       `habitat_specificity_role="External_habitat_viral_specificity_control"`.
    4. Only deliberately selected verified non-viral plasmid or cellular reference controls
       receive `Decoy_curated_*` non-viral labels (`is_decoy=True`).
    5. Missing annotations (`None`/`NA`) are distinguished from scored zero (`0`), and
       overlapping counts across Foldseek, geNomad, archaeal HMMs, and CheckV are never
       summed uncritically as independent genes.
    6. `CheckV Not-determined` + no recognised hallmark is `Unknown_unverified_viral_identity`
       (NEVER `Verified_nonviral` and NEVER turned into a negative control).
    """
    _, _, _, prov_start, prov_end = qualify_namespaced_contig_id(
        target_id,
        fallback_sample=str(target_meta.get("source_sample", "") or ""),
    )
    is_same_base_contig, target_qualified_contig, host_qualified_contig = are_same_namespaced_contig(
        target_id=target_id,
        host_contig_id=host_contig_id,
        host_mag=host_mag,
        target_sample_hint=str(target_meta.get("source_sample", "") or ""),
        target_contig_sha256=str(target_meta.get("target_contig_sha256", "") or ""),
        host_contig_sha256=str(target_meta.get("host_contig_sha256", "") or ""),
    )

    if prov_start is not None:
        parent_hit_lo = prov_start + int(target_hit_lo) - 1
        parent_hit_hi = prov_start + int(target_hit_hi) - 1
    else:
        parent_hit_lo = int(target_hit_lo)
        parent_hit_hi = int(target_hit_hi)

    target_in_host_mag = bool(target_meta.get("contig_in_same_host_mag", False)) or is_same_base_contig

    # 1. Technical self-CRISPR-array match mask (only on the exact same namespaced contig or target's own CRISPR loci)
    arrays_to_check: list[tuple[int, int]] = []
    if is_same_base_contig and host_array_lo > 0 and host_array_hi > 0:
        arrays_to_check.append((min(host_array_lo, host_array_hi), max(host_array_lo, host_array_hi)))
    if all_crispr_arrays_on_contig:
        for alo, ahi in all_crispr_arrays_on_contig:
            arrays_to_check.append((min(alo, ahi), max(alo, ahi)))

    for alo, ahi in arrays_to_check:
        if parent_hit_lo <= (ahi + 500) and parent_hit_hi >= (alo - 500):
            return {
                "target_qualified_contig": target_qualified_contig,
                "host_qualified_contig": host_qualified_contig,
                "is_same_namespaced_contig": is_same_base_contig,
                "viral_identity_status": "Technical_self_crispr_array_masked",
                "habitat_specificity_role": "Technical_crispr_array_locus",
                "is_verified_viral": False,
                "is_provisional_viral": False,
                "is_unknown_viral_identity": False,
                "is_cellular_like": False,
                "is_technical_crispr_array_self_match": True,
                "is_decoy": False,
                "is_external_viral_specificity_control": False,
                "min_independent_viral_orfs": 0,
                "status_detail": "technical_spacer_match_to_crispr_array_locus_masked",
                "provirus_boundary_status": "masked_crispr_array_locus",
            }

    # Parse optional annotation fields (distinguishing None/missing from 0!)
    is_curated_plasmid_ctrl = bool(
        target_meta.get("is_curated_plasmid_control", False) or target_meta.get("in_plasmid_pool", False)
    )
    is_curated_cellular_ctrl = bool(target_meta.get("is_curated_cellular_control", False))
    is_external_viral_ctrl = bool(
        target_meta.get("is_external_habitat_viral_control", False)
        or target_meta.get("is_refseq_non_eq_decoy", False)
    )

    plasmid_hallmarks = _opt_int(target_meta.get("plasmid_hallmarks"))
    virus_hallmarks = _opt_int(target_meta.get("virus_hallmarks"))
    arch_hallmarks = _opt_int(target_meta.get("archaeal_viral_hallmarks"))
    foldseek_viral_hallmarks = _opt_int(target_meta.get("foldseek_viral_hallmarks"))
    checkv_viral_genes = _opt_int(target_meta.get("checkv_viral_genes"))
    checkv_host_genes = _opt_int(target_meta.get("checkv_host_genes"))
    distinct_viral_orf_count = _opt_int(target_meta.get("distinct_viral_orf_count"))
    genomad_virus_score = _opt_float(target_meta.get("genomad_virus_score"))

    raw_cq = target_meta.get("checkv_quality")
    if raw_cq is None or str(raw_cq).strip() in ("", "NA", "NaN", "None", "null", "unscored", "missing"):
        checkv_quality = "Missing_unscored"
    else:
        checkv_quality = str(raw_cq).strip()

    # Filename `|provirus_start_end` is a predicted excision coordinate only (used for parent-contig
    # coordinate mapping), NEVER verified host flanks or a verified viral region.
    has_host_flank_support = bool(target_meta.get("has_verified_host_flanks", False))
    read_junction_status = str(target_meta.get("read_junction_support_status", "untested_missing") or "untested_missing")
    provirus_viral_lo = _opt_int(target_meta.get("verified_viral_region_lo"))
    provirus_viral_hi = _opt_int(target_meta.get("verified_viral_region_hi"))

    # Avoid summing overlapping tool counts (Foldseek/geNomad/CheckV/archaeal HMMs can hit the same ORF).
    # Use `distinct_viral_orf_count` (union of viral ORF IDs) if available; otherwise use the
    # single-source maximum as a conservative non-additive lower bound on distinct viral genes.
    known_viral_counts = [
        c for c in (virus_hallmarks, arch_hallmarks, foldseek_viral_hallmarks, checkv_viral_genes)
        if c is not None
    ]
    if distinct_viral_orf_count is not None:
        min_independent_viral_orfs: int | None = distinct_viral_orf_count
    elif known_viral_counts:
        min_independent_viral_orfs = max(known_viral_counts)
    else:
        min_independent_viral_orfs = None

    has_seq_hallmark = any((c is not None and c > 0) for c in (virus_hallmarks, arch_hallmarks))
    has_struct_hallmark = (foldseek_viral_hallmarks is not None and foldseek_viral_hallmarks > 0)
    has_any_hallmark = has_seq_hallmark or has_struct_hallmark
    has_checkv_viral = (checkv_viral_genes is not None and checkv_viral_genes > 0)
    has_any_viral_marker = (
        has_any_hallmark
        or has_checkv_viral
        or (min_independent_viral_orfs is not None and min_independent_viral_orfs > 0)
    )
    has_positive_host_genes = (checkv_host_genes is not None and checkv_host_genes > 0)
    has_positive_plasmid_hallmarks = (plasmid_hallmarks is not None and plasmid_hallmarks > 0)
    all_viral_annotations_missing = (
        min_independent_viral_orfs is None and checkv_quality == "Missing_unscored" and genomad_virus_score is None
    )

    # 2. External habitat reference viruses (e.g. RefSeq phages from non-desert hosts):
    # These are STILL viral sequences (`is_verified_viral=True`, `is_decoy=False`), NOT negative
    # viral-identity controls. Their role as sample/mapping specificity controls is recorded separately.
    if is_external_viral_ctrl:
        return {
            "target_qualified_contig": target_qualified_contig,
            "host_qualified_contig": host_qualified_contig,
            "is_same_namespaced_contig": is_same_base_contig,
            "viral_identity_status": "Verified_viral_external_habitat_reference",
            "habitat_specificity_role": "External_habitat_viral_specificity_control",
            "is_verified_viral": True,
            "is_provisional_viral": False,
            "is_unknown_viral_identity": False,
            "is_cellular_like": False,
            "is_technical_crispr_array_self_match": False,
            "is_decoy": False,
            "is_external_viral_specificity_control": True,
            "min_independent_viral_orfs": min_independent_viral_orfs if min_independent_viral_orfs is not None else 1,
            "status_detail": "external_habitat_reference_virus_is_viral_sequence_used_for_mapping_specificity_not_identity_decoy",
            "provirus_boundary_status": "external_reference_virus",
        }

    # 3. Deliberately selected curated non-viral controls ONLY get `Decoy_curated_*` labels
    if is_curated_plasmid_ctrl and not has_any_viral_marker:
        return {
            "target_qualified_contig": target_qualified_contig,
            "host_qualified_contig": host_qualified_contig,
            "is_same_namespaced_contig": is_same_base_contig,
            "viral_identity_status": "Decoy_curated_plasmid_control",
            "habitat_specificity_role": "Curated_nonviral_plasmid_control",
            "is_verified_viral": False,
            "is_provisional_viral": False,
            "is_unknown_viral_identity": False,
            "is_cellular_like": False,
            "is_technical_crispr_array_self_match": False,
            "is_decoy": True,
            "is_external_viral_specificity_control": False,
            "min_independent_viral_orfs": 0,
            "status_detail": "deliberately_selected_verified_plasmid_reference_control",
            "provirus_boundary_status": "not_applicable",
        }

    if is_curated_cellular_ctrl and not has_any_viral_marker:
        return {
            "target_qualified_contig": target_qualified_contig,
            "host_qualified_contig": host_qualified_contig,
            "is_same_namespaced_contig": is_same_base_contig,
            "viral_identity_status": "Decoy_curated_cellular_reference_control",
            "habitat_specificity_role": "Curated_nonviral_cellular_control",
            "is_verified_viral": False,
            "is_provisional_viral": False,
            "is_unknown_viral_identity": False,
            "is_cellular_like": True,
            "is_technical_crispr_array_self_match": False,
            "is_decoy": True,
            "is_external_viral_specificity_control": False,
            "min_independent_viral_orfs": 0,
            "status_detail": "deliberately_selected_verified_cellular_chromosomal_reference_control",
            "provirus_boundary_status": "not_applicable",
        }

    # 4. Same chromosomal contig / provirus boundary evaluation
    if target_in_host_mag or prov_start is not None or str(target_meta.get("checkv_provirus", "No")) == "Yes":
        if provirus_viral_lo is not None and provirus_viral_hi is not None:
            v_lo = int(provirus_viral_lo)
            v_hi = int(provirus_viral_hi)
            hit_inside_viral = (parent_hit_lo >= v_lo and parent_hit_hi <= v_hi)
            if not hit_inside_viral:
                return {
                    "target_qualified_contig": target_qualified_contig,
                    "host_qualified_contig": host_qualified_contig,
                    "is_same_namespaced_contig": is_same_base_contig,
                    "viral_identity_status": "Cellular_like_flanking_host_region_outside_provirus",
                    "habitat_specificity_role": "EQ_native_sequence",
                    "is_verified_viral": False,
                    "is_provisional_viral": False,
                    "is_unknown_viral_identity": False,
                    "is_cellular_like": True,
                    "is_technical_crispr_array_self_match": False,
                    "is_decoy": False,
                    "is_external_viral_specificity_control": False,
                    "min_independent_viral_orfs": min_independent_viral_orfs,
                    "status_detail": "spacer_hit_falls_in_flanking_host_dna_outside_provirus_region",
                    "provirus_boundary_status": "rejected_hit_outside_viral_interval",
                }
            if target_in_host_mag and (not has_host_flank_support or not has_any_viral_marker):
                return {
                    "target_qualified_contig": target_qualified_contig,
                    "host_qualified_contig": host_qualified_contig,
                    "is_same_namespaced_contig": is_same_base_contig,
                    "viral_identity_status": "Unresolved_unexcised_provirus_or_unverified_host_contig",
                    "habitat_specificity_role": "EQ_native_sequence",
                    "is_verified_viral": False,
                    "is_provisional_viral": has_any_viral_marker,
                    "is_unknown_viral_identity": not has_any_viral_marker,
                    "is_cellular_like": has_positive_host_genes and not has_any_viral_marker,
                    "is_technical_crispr_array_self_match": False,
                    "is_decoy": False,
                    "is_external_viral_specificity_control": False,
                    "min_independent_viral_orfs": min_independent_viral_orfs,
                    "status_detail": "same_host_contig_lacks_verified_host_flanks_or_viral_hallmarks",
                    "provirus_boundary_status": "unresolved_provirus_boundary_on_host_contig",
                }
        elif target_in_host_mag:
            if not has_any_viral_marker and not all_viral_annotations_missing:
                return {
                    "target_qualified_contig": target_qualified_contig,
                    "host_qualified_contig": host_qualified_contig,
                    "is_same_namespaced_contig": is_same_base_contig,
                    "viral_identity_status": "Cellular_like_chromosomal_contig_not_provirus",
                    "habitat_specificity_role": "EQ_native_sequence",
                    "is_verified_viral": False,
                    "is_provisional_viral": False,
                    "is_unknown_viral_identity": False,
                    "is_cellular_like": True,
                    "is_technical_crispr_array_self_match": False,
                    "is_decoy": False,
                    "is_external_viral_specificity_control": False,
                    "min_independent_viral_orfs": 0,
                    "read_junction_support_status": read_junction_status,
                    "status_detail": "same_host_mag_chromosomal_contig_without_viral_genes_or_provirus_boundaries",
                    "provirus_boundary_status": "rejected_not_provirus_chromosomal_contig",
                }
            if has_any_viral_marker:
                return {
                    "target_qualified_contig": target_qualified_contig,
                    "host_qualified_contig": host_qualified_contig,
                    "is_same_namespaced_contig": is_same_base_contig,
                    "viral_identity_status": "Unresolved_unexcised_provirus_candidate",
                    "habitat_specificity_role": "EQ_native_sequence",
                    "is_verified_viral": False,
                    "is_provisional_viral": True,
                    "is_unknown_viral_identity": False,
                    "is_cellular_like": False,
                    "is_technical_crispr_array_self_match": False,
                    "is_decoy": False,
                    "is_external_viral_specificity_control": False,
                    "min_independent_viral_orfs": min_independent_viral_orfs,
                    "read_junction_support_status": read_junction_status,
                    "status_detail": "unexcised_host_contig_with_viral_markers_requires_boundary_resolution",
                    "provirus_boundary_status": "unresolved_unexcised_provirus_candidate",
                }

    # 5. Plasmid-like candidate (carries plasmid hallmarks without viral markers)
    if has_positive_plasmid_hallmarks and not has_any_viral_marker:
        return {
            "target_qualified_contig": target_qualified_contig,
            "host_qualified_contig": host_qualified_contig,
            "is_same_namespaced_contig": is_same_base_contig,
            "viral_identity_status": "Plasmid_like_nonviral_element",
            "habitat_specificity_role": "EQ_native_sequence",
            "is_verified_viral": False,
            "is_provisional_viral": False,
            "is_unknown_viral_identity": False,
            "is_cellular_like": False,
            "is_technical_crispr_array_self_match": False,
            "is_decoy": False,
            "is_external_viral_specificity_control": False,
            "min_independent_viral_orfs": 0,
            "status_detail": "carries_plasmid_hallmarks_with_zero_viral_markers",
            "provirus_boundary_status": "not_applicable",
        }

    # 6. Cellular-like contig based on positive chromosomal/host gene evidence (not missing/absent annotation!)
    in_cellular_bin = bool(target_meta.get("in_cellular_mag_contig_index", False))
    if not has_any_viral_marker and not all_viral_annotations_missing and (has_positive_host_genes or in_cellular_bin):
        return {
            "target_qualified_contig": target_qualified_contig,
            "host_qualified_contig": host_qualified_contig,
            "is_same_namespaced_contig": is_same_base_contig,
            "viral_identity_status": "Cellular_like_host_gene_supported",
            "habitat_specificity_role": "EQ_native_sequence",
            "is_verified_viral": False,
            "is_provisional_viral": False,
            "is_unknown_viral_identity": False,
            "is_cellular_like": True,
            "is_technical_crispr_array_self_match": False,
            "is_decoy": False,
            "is_external_viral_specificity_control": False,
            "min_independent_viral_orfs": 0,
            "status_detail": "positive_host_gene_or_cellular_bin_evidence_with_zero_viral_markers",
            "provirus_boundary_status": "not_provirus",
        }

    # 7. Unknown / Unverified viral identity:
    # Either annotations are missing (`all_viral_annotations_missing`), or CheckV is Not-determined/Low-quality
    # with 0 recognised viral markers and 0 host genes. NEVER classify as Verified_nonviral or Decoy!
    if not has_any_viral_marker:
        detail = (
            "missing_unscored_viral_and_checkv_annotations_not_zero"
            if all_viral_annotations_missing
            else f"scored_checkv_{checkv_quality}_zero_recognised_hallmarks_and_zero_host_genes_is_unknown_not_nonviral"
        )
        return {
            "target_qualified_contig": target_qualified_contig,
            "host_qualified_contig": host_qualified_contig,
            "is_same_namespaced_contig": is_same_base_contig,
            "viral_identity_status": "Unknown_unverified_viral_identity",
            "habitat_specificity_role": "EQ_native_sequence",
            "is_verified_viral": False,
            "is_provisional_viral": False,
            "is_unknown_viral_identity": True,
            "is_cellular_like": False,
            "is_technical_crispr_array_self_match": False,
            "is_decoy": False,
            "is_external_viral_specificity_control": False,
            "min_independent_viral_orfs": min_independent_viral_orfs,
            "status_detail": detail,
            "provirus_boundary_status": "unknown_unverified",
        }

    # 8. Unexcised cellular bin contig that ALSO carries viral markers -> unresolved unexcised provirus candidate
    has_checkv_hq_mq = checkv_quality in ("Complete", "High-quality", "Medium-quality")
    if in_cellular_bin and prov_start is None and not has_host_flank_support and has_positive_host_genes:
        return {
            "target_qualified_contig": target_qualified_contig,
            "host_qualified_contig": host_qualified_contig,
            "is_same_namespaced_contig": is_same_base_contig,
            "viral_identity_status": "Unresolved_unexcised_provirus_candidate",
            "habitat_specificity_role": "EQ_native_sequence",
            "is_verified_viral": False,
            "is_provisional_viral": True,
            "is_unknown_viral_identity": False,
            "is_cellular_like": False,
            "is_technical_crispr_array_self_match": False,
            "is_decoy": False,
            "is_external_viral_specificity_control": False,
            "min_independent_viral_orfs": min_independent_viral_orfs,
            "status_detail": "cellular_bin_contig_with_both_host_and_viral_genes_unexcised_provirus_candidate",
            "provirus_boundary_status": "unresolved_unexcised_provirus_candidate",
        }

    # 9. Independent viral verification tiers (using non-additive `min_independent_viral_orfs` lower bound so
    # overlapping Foldseek/CheckV/geNomad hits on a single ORF are never double-counted as >=2 independent viral genes)
    has_verified_proviral_interval_and_flanks = (
        has_host_flank_support
        and provirus_viral_lo is not None
        and provirus_viral_hi is not None
    )
    prov_status = (
        "verified_provirus_with_host_flanks"
        if has_verified_proviral_interval_and_flanks
        else (
            "predicted_excised_provirus_interval_unverified_flanks"
            if prov_start is not None
            else "unintegrated_viral_contig"
        )
    )
    n_indep_orfs = min_independent_viral_orfs or 0
    g_score = genomad_virus_score if genomad_virus_score is not None else 0.0

    if (has_checkv_hq_mq and (has_checkv_viral or has_any_hallmark)) or (
        has_any_hallmark and n_indep_orfs >= 2 and g_score >= 0.70
    ):
        return {
            "target_qualified_contig": target_qualified_contig,
            "host_qualified_contig": host_qualified_contig,
            "is_same_namespaced_contig": is_same_base_contig,
            "viral_identity_status": "Verified_viral_independent_hallmark_and_checkv",
            "habitat_specificity_role": "EQ_native_sequence",
            "is_verified_viral": True,
            "is_provisional_viral": False,
            "is_unknown_viral_identity": False,
            "is_cellular_like": False,
            "is_technical_crispr_array_self_match": False,
            "is_decoy": False,
            "is_external_viral_specificity_control": False,
            "min_independent_viral_orfs": n_indep_orfs,
            "read_junction_support_status": read_junction_status,
            "status_detail": "verified_independent_viral_hallmark_and_checkv_support",
            "provirus_boundary_status": prov_status,
        }

    return {
        "target_qualified_contig": target_qualified_contig,
        "host_qualified_contig": host_qualified_contig,
        "is_same_namespaced_contig": is_same_base_contig,
        "viral_identity_status": "Provisional_viral_low_completeness_fragment",
        "habitat_specificity_role": "EQ_native_sequence",
        "is_verified_viral": False,
        "is_provisional_viral": True,
        "is_unknown_viral_identity": False,
        "is_cellular_like": False,
        "is_technical_crispr_array_self_match": False,
        "is_decoy": False,
        "is_external_viral_specificity_control": False,
        "min_independent_viral_orfs": n_indep_orfs,
        "read_junction_support_status": read_junction_status,
        "status_detail": "low_completeness_viral_fragment_with_partial_marker_support",
        "provirus_boundary_status": prov_status,
    }


def classify_host_link_evidence_tier(
    host_reconciliation_status: str,
    viral_eval: dict[str, Any],
    n_strict_spacers: int,
    n_exact_0mm_spacers: int,
    n_exploratory_spacers: int,
    has_cas_operon_in_mag: bool,
    has_verified_provirus_junction_in_mag: bool,
    has_classifier_corroboration: bool,
    n_distinct_host_genera_for_virus: int,
    n_distinct_host_phyla_for_virus: int,
) -> dict[str, Any]:
    """Assign an explicit evidence tier combining Axis 1 (CRISPR/junction) and Axis 2 (viral identity).

    Keeps CRISPR match evidence, viral identity, and habitat specificity strictly separate:
    - Technical self-CRISPR matches -> `Technical_self_crispr_array_masked`
    - External habitat reference virus matches -> `External_habitat_viral_control_match_retained`
      (recorded as a mapping/habitat specificity control match; NEVER treated as non-viral
       or proof that the host is wrong)
    - Curated non-viral plasmid/cellular decoys -> `Decoy_control_rejected`
    - Cellular-like / chromosomal / host-flank hits -> `Cellular_or_plasmid_excluded`
    - Unknown viral identity (`CheckV Not-determined` + no hallmark, or missing annotation) ->
      `Tier3_unverified_viral_identity_retained`
    - Unexcised provirus candidates -> `Tier3_unresolved_unexcised_provirus_retained`
    - Ambiguous multi-host across genera/phyla -> `Tier3_ambiguous_multihost_retained`
    - Verified EQ viral + strong direct host link -> `Tier1_verified_direct` or
      `Tier2_corroborated_single_spacer_or_provirus`
    """
    if host_reconciliation_status.startswith("rejected_"):
        return {
            "evidence_tier": "Rejected_unreconciled_host_alias",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": host_reconciliation_status,
        }

    if viral_eval.get("is_technical_crispr_array_self_match", False):
        return {
            "evidence_tier": "Technical_self_crispr_array_masked",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": str(viral_eval.get("status_detail", "technical_crispr_self_match")),
        }

    if viral_eval.get("is_external_viral_specificity_control", False):
        return {
            "evidence_tier": "External_habitat_viral_control_match_retained",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": (
                "external_non_eq_reference_virus_recorded_for_habitat_mapping_specificity_"
                "not_nonviral_and_not_wrong_host"
            ),
        }

    if viral_eval.get("is_decoy", False):
        return {
            "evidence_tier": "Decoy_control_rejected",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": str(viral_eval.get("viral_identity_status", "Decoy_control_rejected")),
        }

    v_status = str(viral_eval.get("viral_identity_status", ""))
    if viral_eval.get("is_cellular_like", False) or v_status == "Plasmid_like_nonviral_element":
        return {
            "evidence_tier": "Cellular_or_plasmid_excluded",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": v_status,
        }

    if viral_eval.get("is_unknown_viral_identity", False):
        return {
            "evidence_tier": "Tier3_unverified_viral_identity_retained",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": "unknown_viral_identity_checkv_not_determined_or_missing_no_hallmark_retained_not_negative_control",
        }

    if "Unresolved_unexcised_provirus" in v_status:
        return {
            "evidence_tier": "Tier3_unresolved_unexcised_provirus_retained",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": "unexcised_provirus_candidate_retained_pending_boundary_resolution",
        }

    is_direct_591 = (host_reconciliation_status == "canonical_591_direct")
    is_verified_viral = bool(viral_eval.get("is_verified_viral", False))
    is_provisional_viral = bool(viral_eval.get("is_provisional_viral", False))
    is_multihost_cross_genus = (int(n_distinct_host_genera_for_virus) > 1 or int(n_distinct_host_phyla_for_virus) > 1)

    if n_strict_spacers == 0 and n_exploratory_spacers == 0 and not has_verified_provirus_junction_in_mag:
        if has_classifier_corroboration:
            return {
                "evidence_tier": "Tier3_classifier_only_no_direct_sequence_link",
                "eligible_stage4_primary": False,
                "eligible_stage4_sensitivity": False,
                "exclusion_reason": "classifier_prediction_alone_is_not_direct_species_proof",
            }
        return {
            "evidence_tier": "Rejected_no_qualifying_evidence",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": "no_spacer_or_provirus_evidence",
        }

    if is_multihost_cross_genus and n_strict_spacers >= 1:
        return {
            "evidence_tier": "Tier3_ambiguous_multihost_retained",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": False,
            "exclusion_reason": (
                f"ambiguous_multihost_target_matches_{n_distinct_host_genera_for_virus}_genera_"
                f"{n_distinct_host_phyla_for_virus}_phyla"
            ),
        }

    if is_direct_591 and is_verified_viral and not is_multihost_cross_genus:
        if (
            n_strict_spacers >= 2
            or (n_exact_0mm_spacers >= 1 and (has_cas_operon_in_mag or has_verified_provirus_junction_in_mag or has_classifier_corroboration))
            or (has_verified_provirus_junction_in_mag and n_strict_spacers >= 1)
        ):
            return {
                "evidence_tier": "Tier1_verified_direct",
                "eligible_stage4_primary": True,
                "eligible_stage4_sensitivity": True,
                "exclusion_reason": "none",
            }

        if n_strict_spacers == 1 or has_verified_provirus_junction_in_mag:
            return {
                "evidence_tier": "Tier2_corroborated_single_spacer_or_provirus",
                "eligible_stage4_primary": False,
                "eligible_stage4_sensitivity": True,
                "exclusion_reason": "single_spacer_or_provirus_only_included_in_sensitivity",
            }

    if is_direct_591 and is_provisional_viral and n_strict_spacers >= 2 and not is_multihost_cross_genus:
        return {
            "evidence_tier": "Tier2_corroborated_single_spacer_or_provirus",
            "eligible_stage4_primary": False,
            "eligible_stage4_sensitivity": True,
            "exclusion_reason": "low_completeness_viral_fragment_with_multi_spacer_support",
        }

    return {
        "evidence_tier": "Tier3_exploratory_or_species_projected",
        "eligible_stage4_primary": False,
        "eligible_stage4_sensitivity": False,
        "exclusion_reason": "weak_or_projected_link_retained_in_report_excluded_from_mechanism_tests",
    }


def summarize_mag_host_disposition(
    mag_id: str,
    fasta_verified: bool,
    crispr_screened: bool,
    n_raw_arrays: int,
    n_verified_spacers: int,
    pair_tiers: list[str],
    caller_failed: bool = False,
) -> str:
    """Return the exact mutually exclusive disposition for each of the 591 archaeal MAGs.

    Never conflates zero detected spacers, spacers with no viral hit, unverified viral
    candidates (`Tier3_unverified_viral_identity_retained`), technical self-matches,
    external habitat viral control matches, curated non-viral decoy hits, failed caller
    execution (`Failed_caller_execution`), and unknown/unscreened status.
    """
    if caller_failed:
        return "Failed_caller_execution"
    if not fasta_verified or not crispr_screened:
        return "Unknown_unscreened"
    if "Tier1_verified_direct" in pair_tiers:
        return "Linked_Tier1_verified_direct"
    if "Tier2_corroborated_single_spacer_or_provirus" in pair_tiers:
        return "Linked_Tier2_corroborated"
    if "Tier3_unverified_viral_identity_retained" in pair_tiers:
        return "Candidate_link_unknown_viral_identity_retained"
    if "Tier3_unresolved_unexcised_provirus_retained" in pair_tiers:
        return "Candidate_link_unresolved_unexcised_provirus_retained"
    if any(t.startswith("Tier3_") for t in pair_tiers):
        return "Linked_Tier3_ambiguous_or_exploratory_only"
    if "External_habitat_viral_control_match_retained" in pair_tiers:
        return "External_habitat_viral_reference_match_only"
    if "Cellular_or_plasmid_excluded" in pair_tiers:
        return "Cellular_or_plasmid_target_excluded_only"
    if "Decoy_control_rejected" in pair_tiers:
        return "Curated_decoy_control_hits_only"
    if "Technical_self_crispr_array_masked" in pair_tiers and n_verified_spacers > 0:
        return "Spacers_detected_only_technical_self_array_masked"
    if n_verified_spacers > 0:
        return "Spacers_detected_no_viral_hit"
    return "Zero_spacers_detected"
