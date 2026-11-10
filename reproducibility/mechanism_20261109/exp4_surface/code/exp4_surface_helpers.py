#!/usr/bin/env python3
"""Core scientific helpers for Experiment 4 (`exp4_surface`):
Functional discrimination of dark exported Thermoplasmatota families across
nutrient acquisition, substrate binding, surface assembly, cell coat, and virus interaction.

Implements all scientific rules and code-review fixes from:
- `mechanism-20261109/EXP4-CONTRACT.txt` (Items 1-6)
- `mechanism-20261109/EXP4-CURATED-SOURCE-NOTE.txt`
- `mechanism-20261109/EXP4-EVIDENCE-REVIEW.txt`
- `mechanism-20261109/EXP4-CODE-REVIEW.txt` (Items 1-7)
"""
from __future__ import annotations

import collections
import hashlib
import math
import re
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd

from build_curated_reference_panel import PANEL_SPECS

SPANY_CLASSES = frozenset({"SP", "LIPO", "TAT", "TATLIPO"})
KNOWN_SIGNALP_CLASSES = frozenset({"SP", "LIPO", "TAT", "TATLIPO", "OTHER"})

ROLE_CLASSES = (
    "nutrient_acquisition_peptidase",
    "substrate_binding_or_uptake",
    "surface_assembly_pilin_archaellin",
    "cell_coat_or_envelope",
    "virus_interaction_or_entry",
)

PANEL_SPECS_BY_ID: dict[str, dict[str, str]] = {s["ref_id"]: dict(s) for s in PANEL_SPECS}

# Bounded token/profile patterns per neighbourhood guild (no raw substring collisions like 'opp' in 'copper' or 'cas1' in unrelated words)
BOUNDED_GUILD_PATTERNS: dict[str, re.Pattern[str]] = {
    "nutrient_acquisition_peptidase": re.compile(
        r"(?:\bmer:[a-z0-9_.-]+|\bmerops:[a-z0-9_.-]+|\bpeptidase\b|\bprotease\b|\bsubtilisin\b|\bthermopsin\b|\btricorn\b|\bmetalloprotease\b|\baminopeptidase\b|\bcarboxypeptidase\b|\bdipeptidase\b|\boligopeptidase\b|\bpf00082\b|\bpf00026\b)",
        re.IGNORECASE,
    ),
    "substrate_binding_or_uptake": re.compile(
        r"(?:\btc:3\.a\.1\b|\btcdb:3\.a\.1\b|\babc\s+transporter\b|\bsolute-binding\b|\bsubstrate-binding\b|\bperiplasmic\s+binding\b|\bsbp_bac_[0-9]+\b|\bpf00496\b|\bpf01547\b|\bpf13416\b|\boppa\b|\bdppa\b|\bappa\b|\bappb\b|\bappc\b|\boppb\b|\boppc\b|\boppd\b|\boppf\b|\bdppb\b|\bdppc\b|\blivk\b|\blivh\b|\blivm\b)",
        re.IGNORECASE,
    ),
    "surface_assembly_pilin_archaellin": re.compile(
        r"(?:\barchaellin\b|\barchaellum\b|\bflagellin\b|\bpilin\b|\bpilus\b|\bflai\b|\bflaj\b|\bflah\b|\bflaf\b|\bflag\b|\bflak\b|\bpibd\b|\barlb\b|\barli\b|\barlj\b|\barlh\b|\baapb\b|\baape\b|\baapf\b|\bcdip[abef]\b|\btype\s+ii\s+secretion\s+system\b|\btype\s+iv\s+pilus\b|\bpilt\b|\bpilb\b|\bpilc\b|\bprepilin\s+peptidase\b|\bpf00437\b|\bpf01478\b)",
        re.IGNORECASE,
    ),
    "cell_coat_or_envelope": re.compile(
        r"(?:\bs-layer\b|\bsurface\s+layer\b|\bduf1608\b|\bpf07752\b|\bpf05123\b|\barchaeosortase\b|\barta\b|\bpgf-cterm\b|\bpef-cterm\b|\boligosaccharyltransferase\b|\baglb\b|\bslaa\b|\bslab\b|\bcsxa\b)",
        re.IGNORECASE,
    ),
    # Host antiviral defense is strictly separated from viral/mobile element structural hallmarks
    "host_defense_antiviral": re.compile(
        r"(?:\bcrispr\b|\bcas1\b|\bcas2\b|\bcas3\b|\bcas4\b|\bcas5\b|\bcas6\b|\bcas7\b|\bcas8\b|\bcas10\b|\bcsa[0-9]+\b|\bcsc[0-9]+\b|\bcmr[0-9]+\b|\bcsm[0-9]+\b|\brestriction\s+endonuclease\b|\babortive\s+infection\b)",
        re.IGNORECASE,
    ),
    "viral_or_mobile_element_origin": re.compile(
        r"(?:\bviral\s+capsid\b|\bmajor\s+capsid\b|\bminor\s+capsid\b|\bportal\s+protein\b|\bterminase\b|\btail\s+fiber\b|\bturret\s+protein\b|\bvirion\b|\bprophage\b|\bprovirus\b)",
        re.IGNORECASE,
    ),
}

PREPILIN_MOTIF_RE = re.compile(r"^[A-Z]{2,18}[KR][GA][ALIFVM][A-Z]{10,28}")
# Overlapping lookahead regex so e.g. 'NNST' yields two sequons ('NNS' at pos 0 and 'NST' at pos 1)
NXST_SEQUON_LOOKAHEAD_RE = re.compile(r"(?=(N[^P][ST]))")
ARCHAEOSORTASE_CTERM_RE = re.compile(r"(?:P[GSA]F|PEF)[A-Z]{2,12}[FILVWM]{12,26}[KRH]{1,8}$")


def clean_aa_sequence(raw_seq: str | None) -> str:
    """Strip whitespace and terminal stop '*' while preserving internal characters for auditing."""
    if not raw_seq or not isinstance(raw_seq, str):
        return ""
    s = "".join(raw_seq.split()).upper()
    if s.endswith("*"):
        s = s[:-1]
    return s


def extract_tmbed_segment_spans(label_seq: str | None) -> dict[str, list[tuple[int, int]]]:
    """Extract 1-based inclusive `(start, end)` coordinate spans for `S`, `H` (`H`/`h`), and `B` (`B`/`b`)
    from a per-residue TMbed prediction label string.
    """
    out: dict[str, list[tuple[int, int]]] = {"S": [], "H": [], "B": []}
    if not label_seq or not isinstance(label_seq, str):
        return out
    s = label_seq.strip()
    for key, chars in (("S", "S"), ("H", "Hh"), ("B", "Bb")):
        in_seg = False
        seg_start = 0
        for idx, ch in enumerate(s, start=1):
            if ch in chars:
                if not in_seg:
                    in_seg = True
                    seg_start = idx
            else:
                if in_seg:
                    out[key].append((seg_start, idx - 1))
                    in_seg = False
        if in_seg:
            out[key].append((seg_start, len(s)))
    return out


def parse_signalp_cleavage_site(cs_pos_str: str | None) -> int | None:
    """Parse SignalP cleavage site coordinate string (e.g. 'CS pos: 24-25. Pr: 0.91') into 1-based last SP residue (24)."""
    if not cs_pos_str or not isinstance(cs_pos_str, str) or pd.isna(cs_pos_str):
        return None
    m = re.search(r"(\d+)\s*-\s*(\d+)", cs_pos_str)
    if not m:
        return None
    return int(m.group(1))


def compute_sequence_integrity_and_compositional_metrics(raw_seq: str | None) -> dict[str, Any]:
    """Compute exact sequence hashes (MD5, SHA256), ORF integrity, compositional bias, and predicted motif counts.

    Uses lookahead regex `(?=(N[^P][ST]))` so overlapping sequons (e.g. `NNST` -> `NNS` and `NST`)
    are all counted. Glycosylation sequons remain explicitly labelled as predicted sequence motifs,
    never measured modifications.
    """
    seq = clean_aa_sequence(raw_seq)
    n = len(seq)
    if n == 0:
        return {
            "aa_length": 0,
            "aa_md5": "missing",
            "aa_sha256": "missing",
            "n_internal_stops": 0,
            "n_ambiguous_residues": 0,
            "shannon_entropy_bits": np.nan,
            "max_homopolymer_run": 0,
            "low_complexity_window20_frac": np.nan,
            "acidic_fraction_DE": np.nan,
            "basic_fraction_KR": np.nan,
            "net_charge_per_residue": np.nan,
            "hydrophobic_fraction_FILVWM": np.nan,
            "ser_thr_fraction_ST": np.nan,
            "predicted_nxst_sequon_count": 0,
            "predicted_nxst_sequon_positions_1based": "none",
            "predicted_nxst_sequons_per_100aa": np.nan,
            "glycosylation_evidence_status": "missing_sequence",
            "has_prepilin_classIII_motif_pred": False,
            "has_archaeosortase_cterm_motif_pred": False,
        }

    md5_hex = hashlib.md5(seq.encode("ascii", errors="replace")).hexdigest()
    sha256_hex = hashlib.sha256(seq.encode("ascii", errors="replace")).hexdigest()
    n_stops = seq.count("*")
    std_aa = frozenset("ACDEFGHIKLMNPQRSTVWY")
    n_ambig = sum(1 for c in seq if c not in std_aa and c != "*")

    counts = collections.Counter(seq)
    entropy = -sum((cnt / n) * math.log2(cnt / n) for cnt in counts.values() if cnt > 0)

    max_run = 1
    cur_run = 1
    for i in range(1, n):
        if seq[i] == seq[i - 1]:
            cur_run += 1
            if cur_run > max_run:
                max_run = cur_run
        else:
            cur_run = 1

    if n >= 20:
        n_win = n - 19
        lc_win = sum(1 for i in range(n_win) if len(set(seq[i : i + 20])) <= 5)
        lc_frac = lc_win / float(n_win)
    else:
        lc_frac = 1.0 if len(set(seq)) <= 5 else 0.0

    n_de = counts.get("D", 0) + counts.get("E", 0)
    n_kr = counts.get("K", 0) + counts.get("R", 0)
    n_hydro = sum(counts.get(c, 0) for c in "FILVWM")
    n_st = counts.get("S", 0) + counts.get("T", 0)

    nxst_matches = list(NXST_SEQUON_LOOKAHEAD_RE.finditer(seq))
    n_nxst = len(nxst_matches)
    nxst_pos = [str(m.start() + 1) for m in nxst_matches]
    nxst_per_100 = (100.0 * n_nxst) / float(n)

    has_prepilin = bool(PREPILIN_MOTIF_RE.search(seq[:65]))
    has_art_cterm = bool(ARCHAEOSORTASE_CTERM_RE.search(seq[-55:])) if n >= 55 else False

    return {
        "aa_length": n,
        "aa_md5": md5_hex,
        "aa_sha256": sha256_hex,
        "n_internal_stops": n_stops,
        "n_ambiguous_residues": n_ambig,
        "shannon_entropy_bits": float(entropy),
        "max_homopolymer_run": int(max_run),
        "low_complexity_window20_frac": float(lc_frac),
        "acidic_fraction_DE": float(n_de / n),
        "basic_fraction_KR": float(n_kr / n),
        "net_charge_per_residue": float((n_kr - n_de) / n),
        "hydrophobic_fraction_FILVWM": float(n_hydro / n),
        "ser_thr_fraction_ST": float(n_st / n),
        "predicted_nxst_sequon_count": int(n_nxst),
        "predicted_nxst_sequon_positions_1based": ";".join(nxst_pos) if nxst_pos else "none",
        "predicted_nxst_sequons_per_100aa": float(nxst_per_100),
        "glycosylation_evidence_status": "predicted_sequence_motif_only_not_observed_modification",
        "has_prepilin_classIII_motif_pred": has_prepilin,
        "has_archaeosortase_cterm_motif_pred": has_art_cterm,
    }


def verify_prediction_sequence_hash_match(
    candidate_raw_seq: str | None,
    prediction_polymer_seq: str | None,
) -> dict[str, Any]:
    """Verify exact sequence identity/hash between a gene product and a structural/functional model input.

    Invalid or missing predictions are explicitly marked and NEVER treated as negative evidence.
    """
    cand = clean_aa_sequence(candidate_raw_seq)
    pred = clean_aa_sequence(prediction_polymer_seq)
    if not cand:
        return {
            "sequence_match_status": "invalid_missing_candidate_sequence",
            "is_valid_exact_match": False,
            "is_missing_or_invalid_not_negative": True,
            "candidate_sha256": "missing",
            "prediction_sha256": "missing",
        }
    cand_sha = hashlib.sha256(cand.encode("ascii", errors="replace")).hexdigest()
    if not pred:
        return {
            "sequence_match_status": "missing_prediction_not_negative",
            "is_valid_exact_match": False,
            "is_missing_or_invalid_not_negative": True,
            "candidate_sha256": cand_sha,
            "prediction_sha256": "missing",
        }
    pred_sha = hashlib.sha256(pred.encode("ascii", errors="replace")).hexdigest()
    if cand == pred and cand_sha == pred_sha:
        return {
            "sequence_match_status": "exact_sequence_hash_match",
            "is_valid_exact_match": True,
            "is_missing_or_invalid_not_negative": False,
            "candidate_sha256": cand_sha,
            "prediction_sha256": pred_sha,
        }
    return {
        "sequence_match_status": f"invalid_sequence_mismatch(cand_len={len(cand)},pred_len={len(pred)})",
        "is_valid_exact_match": False,
        "is_missing_or_invalid_not_negative": True,
        "candidate_sha256": cand_sha,
        "prediction_sha256": pred_sha,
    }


def extract_tmbed_segments_1based(lab: str | None, chars: str = "Hh") -> list[tuple[int, int]]:
    """Extract 1-based inclusive (start, end) residue coordinate spans for contiguous segments of `chars` in a TMbed label string."""
    if not lab:
        return []
    spans: list[tuple[int, int]] = []
    in_seg = False
    seg_start = 0
    for idx, ch in enumerate(lab, start=1):
        if ch in chars:
            if not in_seg:
                in_seg = True
                seg_start = idx
        else:
            if in_seg:
                spans.append((seg_start, idx - 1))
                in_seg = False
    if in_seg:
        spans.append((seg_start, len(lab)))
    return spans


def classify_export_and_surface_topology(
    signalp_class: Any,
    tmbed_s: Any,
    tmbed_h: Any,
    tmbed_b: Any,
    *,
    aa_length: int | None = None,
    signalp_cs_pos: int | None = None,
    tm_helix_spans_1based: list[tuple[int, int]] | None = None,
    has_prepilin_motif: bool = False,
    has_archaeosortase_cterm_motif: bool = False,
    sequence_valid: bool = True,
) -> dict[str, Any]:
    """Explicitly distinguish `exported != freely_secreted` across Sec/Tat/lipoprotein/pilin and membrane anchors.

    Enforces `EXP4-CODE-REVIEW.txt` items 2 and 3:
    - Missing `tmbed_s`, `tmbed_h`, or `tmbed_b` (`None` or `pd.isna(...)`) is NEVER coerced to `False`/`0`.
      When TMbed is missing, `topology_status = "unknown_missing_tmbed"` and `is_freely_secreted_candidate = None`.
    - Unknown/missing `signalp_class` (anything not in `KNOWN_SIGNALP_CLASSES`) is NEVER treated as non-secretory.
    - `h_cnt == 1` is NEVER assumed to be an N-terminal anchor without verified segment coordinates
      (`tm_helix_spans_1based`) and cleavage-site overlap checking (`signalp_cs_pos`).
    """
    if not sequence_valid:
        return {
            "is_exported_any": None,
            "is_freely_secreted_candidate": None,
            "is_surface_or_membrane_anchored": None,
            "topology_status": "invalid_sequence",
            "tm_anchor_position_status": "invalid_sequence",
            "secretion_pathway_class": "missing_or_invalid_not_negative",
            "surface_localization_class": "prediction_missing_or_invalid_not_negative",
            "exported_vs_secreted_note": "sequence invalid; not counted as negative",
        }

    sp_missing = signalp_class is None or pd.isna(signalp_class)
    sp = "" if sp_missing else str(signalp_class).strip().upper()
    sp_known = sp in KNOWN_SIGNALP_CLASSES

    tm_missing = (
        tmbed_s is None
        or pd.isna(tmbed_s)
        or tmbed_h is None
        or pd.isna(tmbed_h)
        or tmbed_b is None
        or pd.isna(tmbed_b)
    )

    if not sp_known and tm_missing:
        return {
            "is_exported_any": None,
            "is_freely_secreted_candidate": None,
            "is_surface_or_membrane_anchored": None,
            "topology_status": "unknown_missing_signalp_and_tmbed",
            "tm_anchor_position_status": "unknown_missing_tmbed",
            "secretion_pathway_class": "unknown_missing_or_unrecognized_signalp",
            "surface_localization_class": "prediction_missing_or_invalid_not_negative",
            "exported_vs_secreted_note": "SignalP and TMbed missing or unrecognized; not counted as negative",
        }

    if tm_missing:
        # SignalP is known, but TMbed topology is missing -> anchor/freely-secreted status is UNKNOWN
        if sp in ("LIPO", "TATLIPO"):
            pathway = "Tat_Lipoprotein_SPII" if sp == "TATLIPO" else "Sec_Lipoprotein_SPII"
            return {
                "is_exported_any": True,
                "is_freely_secreted_candidate": False,
                "is_surface_or_membrane_anchored": True,
                "topology_status": "unknown_missing_tmbed",
                "tm_anchor_position_status": "lipobox_n_terminal_tmbed_missing",
                "secretion_pathway_class": pathway,
                "surface_localization_class": "lipoprotein_membrane_anchored_SPII_tmbed_unknown",
                "exported_vs_secreted_note": "SignalP lipoprotein SPII known; TMbed topology missing",
            }
        if sp in ("SP", "TAT"):
            pathway = "Tat_SPI" if sp == "TAT" else "Sec_SPI"
            return {
                "is_exported_any": True,
                "is_freely_secreted_candidate": None,
                "is_surface_or_membrane_anchored": None,
                "topology_status": "unknown_missing_tmbed",
                "tm_anchor_position_status": "unknown_missing_tmbed",
                "secretion_pathway_class": pathway,
                "surface_localization_class": "exported_signalp_known_topology_unknown_missing_tmbed",
                "exported_vs_secreted_note": "SignalP export known, but TMbed topology missing so freely secreted vs membrane-anchored is unknown",
            }
        return {
            "is_exported_any": None,
            "is_freely_secreted_candidate": None,
            "is_surface_or_membrane_anchored": None,
            "topology_status": "unknown_missing_tmbed",
            "tm_anchor_position_status": "unknown_missing_tmbed",
            "secretion_pathway_class": "SignalP_OTHER_tmbed_missing",
            "surface_localization_class": "signalp_other_topology_unknown_missing_tmbed",
            "exported_vs_secreted_note": "SignalP OTHER with missing TMbed topology; export/membrane status unknown",
        }

    # Here TMbed values are present and non-NA; check boolean/integer values safely after pd.isna guard
    if isinstance(tmbed_s, str):
        s_flag = tmbed_s.strip().lower() in ("true", "1", "yes")
    else:
        s_flag = bool(tmbed_s)
    h_cnt = int(tmbed_h)
    b_cnt = int(tmbed_b)

    if not sp_known:
        # SignalP is missing/unrecognized while TMbed is present: cannot claim negative SignalP
        exported_by_tm = True if s_flag else None
        return {
            "is_exported_any": exported_by_tm,
            "is_freely_secreted_candidate": None,
            "is_surface_or_membrane_anchored": True if (h_cnt >= 1 or b_cnt >= 1) else None,
            "topology_status": f"tmbed_present_signalp_unknown(S={int(s_flag)},H={h_cnt},B={b_cnt})",
            "tm_anchor_position_status": "signalp_missing_cleavage_unknown",
            "secretion_pathway_class": "unknown_missing_or_unrecognized_signalp",
            "surface_localization_class": (
                "exported_tmbed_S_only_signalp_missing" if s_flag else "signalp_missing_tmbed_only_status_unresolved"
            ),
            "exported_vs_secreted_note": "SignalP missing/unrecognized; cannot classify SPII/Tat or freely secreted status",
        }

    sp_any = sp in SPANY_CLASSES

    # Inspect verified TM helix spans relative to SignalP cleavage site if segment coordinates are provided
    post_cleavage_spans: list[tuple[int, int]] | None = None
    overlaps_signal_peptide_only = False
    if tm_helix_spans_1based is not None:
        if len(tm_helix_spans_1based) != h_cnt:
            # Segment list count does not match h_cnt -> treat coordinates as unverified
            post_cleavage_spans = None
        else:
            cs_cut = signalp_cs_pos if (signalp_cs_pos is not None and signalp_cs_pos > 0) else (35 if (sp_any or s_flag) else 0)
            kept: list[tuple[int, int]] = []
            for st, en in tm_helix_spans_1based:
                if cs_cut > 0 and en <= cs_cut:
                    continue
                kept.append((st, en))
            post_cleavage_spans = kept
            if h_cnt >= 1 and len(kept) == 0 and (sp_any or s_flag):
                overlaps_signal_peptide_only = True

    # Check if N-terminal class III prepilin motif has a verified N-terminal helix
    has_verified_n_term_helix = bool(
        tm_helix_spans_1based is not None
        and len(tm_helix_spans_1based) == 1
        and tm_helix_spans_1based[0][0] <= 35
        and tm_helix_spans_1based[0][1] <= 65
    )

    exported_any = bool(sp_any or s_flag or (has_prepilin_motif and has_verified_n_term_helix and b_cnt == 0))

    if sp in ("LIPO", "TATLIPO"):
        pathway = "Tat_Lipoprotein_SPII" if sp == "TATLIPO" else "Sec_Lipoprotein_SPII"
    elif sp == "TAT":
        pathway = "Tat_SPI"
    elif sp == "SP":
        pathway = "Sec_SPI"
    elif has_prepilin_motif and has_verified_n_term_helix and b_cnt == 0:
        pathway = "ClassIII_Prepilin_Archaellin_PibD_candidate"
    elif s_flag:
        pathway = "TMbed_SignalPeptide_only"
    else:
        pathway = "Non_secretory"

    if not exported_any:
        if h_cnt >= 1 or b_cnt >= 1:
            return {
                "is_exported_any": False,
                "is_freely_secreted_candidate": False,
                "is_surface_or_membrane_anchored": True,
                "topology_status": f"verified_tmbed(S={int(s_flag)},H={h_cnt},B={b_cnt})",
                "tm_anchor_position_status": "integral_membrane_non_sp",
                "secretion_pathway_class": pathway,
                "surface_localization_class": "integral_membrane_non_sp",
                "exported_vs_secreted_note": f"non-secretory integral membrane protein (TM_H={h_cnt}, TM_B={b_cnt})",
            }
        return {
            "is_exported_any": False,
            "is_freely_secreted_candidate": False,
            "is_surface_or_membrane_anchored": False,
            "topology_status": "verified_tmbed(S=0,H=0,B=0)",
            "tm_anchor_position_status": "none_soluble",
            "secretion_pathway_class": pathway,
            "surface_localization_class": "intracellular_soluble",
            "exported_vs_secreted_note": "intracellular soluble protein (SignalP=OTHER, TMbed S=0, H=0, B=0)",
        }

    # Exported proteins
    if sp in ("LIPO", "TATLIPO"):
        return {
            "is_exported_any": True,
            "is_freely_secreted_candidate": False,
            "is_surface_or_membrane_anchored": True,
            "topology_status": f"verified_tmbed(S={int(s_flag)},H={h_cnt},B={b_cnt})",
            "tm_anchor_position_status": "n_terminal_lipobox_SPII",
            "secretion_pathway_class": pathway,
            "surface_localization_class": "lipoprotein_membrane_anchored_SPII",
            "exported_vs_secreted_note": "exported N-terminally lipid-anchored surface/membrane protein (not freely secreted)",
        }

    if b_cnt >= 1 or h_cnt >= 2:
        return {
            "is_exported_any": True,
            "is_freely_secreted_candidate": False,
            "is_surface_or_membrane_anchored": True,
            "topology_status": f"verified_tmbed(S={int(s_flag)},H={h_cnt},B={b_cnt})",
            "tm_anchor_position_status": "multipass_or_beta_barrel",
            "secretion_pathway_class": pathway,
            "surface_localization_class": "multipass_or_beta_barrel_membrane_with_export_signal",
            "exported_vs_secreted_note": f"exported polytopic/beta-barrel membrane protein (TM_H={h_cnt}, TM_B={b_cnt}; not freely secreted)",
        }

    if has_prepilin_motif and has_verified_n_term_helix and b_cnt == 0:
        return {
            "is_exported_any": True,
            "is_freely_secreted_candidate": False,
            "is_surface_or_membrane_anchored": True,
            "topology_status": f"verified_tmbed(S={int(s_flag)},H={h_cnt},B={b_cnt})",
            "tm_anchor_position_status": "verified_n_terminal_prepilin_alpha1_helix",
            "secretion_pathway_class": pathway,
            "surface_localization_class": "prepilin_or_archaellin_filament_subunit_candidate",
            "exported_vs_secreted_note": "predicted class III prepilin/archaellin subunit with verified N-terminal hydrophobic helix (not freely secreted)",
        }

    if h_cnt == 1:
        # Must check actual segment coordinates before claiming N-terminal vs C-terminal vs cleaved signal helix!
        if post_cleavage_spans is None:
            return {
                "is_exported_any": True,
                "is_freely_secreted_candidate": None,
                "is_surface_or_membrane_anchored": None,
                "topology_status": f"tmbed_summary_only_no_segment_coords(S={int(s_flag)},H=1,B=0)",
                "tm_anchor_position_status": "unresolved_without_segment_coordinates",
                "secretion_pathway_class": pathway,
                "surface_localization_class": "exported_single_tm_segment_position_unresolved_signal_vs_anchor",
                "exported_vs_secreted_note": "H=1 without segment coordinates: cannot resolve N/C/internal anchor vs uncleaved signal helix",
            }
        if overlaps_signal_peptide_only and len(post_cleavage_spans) == 0:
            return {
                "is_exported_any": True,
                "is_freely_secreted_candidate": True,
                "is_surface_or_membrane_anchored": False,
                "topology_status": f"verified_tmbed_H1_within_cleaved_SP(S={int(s_flag)},H=1,B=0)",
                "tm_anchor_position_status": "h1_overlaps_cleaved_signal_peptide_post_cleavage_0_tm",
                "secretion_pathway_class": pathway,
                "surface_localization_class": "freely_secreted_extracellular_candidate",
                "exported_vs_secreted_note": "single TMbed H segment lies entirely within cleaved N-terminal signal peptide; mature protein has 0 TM segments",
            }
        st_pos, en_pos = post_cleavage_spans[0]
        if aa_length is not None and aa_length >= 90 and st_pos >= aa_length - 45:
            return {
                "is_exported_any": True,
                "is_freely_secreted_candidate": False,
                "is_surface_or_membrane_anchored": True,
                "topology_status": f"verified_tmbed(S={int(s_flag)},H=1[{st_pos}-{en_pos}],B=0)",
                "tm_anchor_position_status": f"verified_c_terminal_tm_anchor[{st_pos}-{en_pos}]",
                "secretion_pathway_class": pathway,
                "surface_localization_class": "c_terminal_anchored_surface_coat_candidate",
                "exported_vs_secreted_note": f"exported protein with verified C-terminal transmembrane anchor [{st_pos}-{en_pos}]",
            }
        if en_pos <= 65:
            return {
                "is_exported_any": True,
                "is_freely_secreted_candidate": False,
                "is_surface_or_membrane_anchored": True,
                "topology_status": f"verified_tmbed(S={int(s_flag)},H=1[{st_pos}-{en_pos}],B=0)",
                "tm_anchor_position_status": f"verified_n_terminal_tm_anchor[{st_pos}-{en_pos}]",
                "secretion_pathway_class": pathway,
                "surface_localization_class": "n_terminal_anchored_single_pass_surface",
                "exported_vs_secreted_note": f"exported protein with verified post-cleavage N-terminal transmembrane anchor [{st_pos}-{en_pos}]",
            }
        return {
            "is_exported_any": True,
            "is_freely_secreted_candidate": False,
            "is_surface_or_membrane_anchored": True,
            "topology_status": f"verified_tmbed(S={int(s_flag)},H=1[{st_pos}-{en_pos}],B=0)",
            "tm_anchor_position_status": f"verified_internal_tm_helix[{st_pos}-{en_pos}]",
            "secretion_pathway_class": pathway,
            "surface_localization_class": "internal_single_pass_membrane_with_export_signal",
            "exported_vs_secreted_note": f"exported protein with verified internal transmembrane helix [{st_pos}-{en_pos}]",
        }

    # h_cnt == 0 and b_cnt == 0
    return {
        "is_exported_any": True,
        "is_freely_secreted_candidate": True,
        "is_surface_or_membrane_anchored": False,
        "topology_status": f"verified_tmbed(S={int(s_flag)},H=0,B=0)",
        "tm_anchor_position_status": "none_0_tm_segments",
        "secretion_pathway_class": pathway,
        "surface_localization_class": "freely_secreted_extracellular_candidate",
        "exported_vs_secreted_note": "signal peptide with verified 0 transmembrane segments (TMbed H=0, B=0); candidate freely secreted or non-covalently S-layer-associated protein",
    }


def resolve_curated_reference_id_from_header(
    raw_target_id: str | None,
    *,
    target_polymer_seq: str | None = None,
    db_lookup_map: dict[str, str] | None = None,
    curated_panel_by_id: dict[str, dict[str, Any]] | None = None,
    require_known: bool = True,
) -> str:
    """ deterministically map a DIAMOND or Foldseek target header/name to a canonical `REF_*` key in `PANEL_SPECS_BY_ID`.

    Enforces `EXP4-PIPELINE-REVIEW.txt` item 1:
    - Parses exact `REF_ID|uniprot:...|role:...` headers from `curated_role_reference_panel.faa` (delimited by `|`).
    - Parses `REF_ID__PDB_...` headers from `curated_pdb_chain_sequences.faa`.
    - Resolves Foldseek target names (including `.pdb` and chain suffixes like `REF_PIL_8REY_CdiPB.pdb_A`)
      via `db_lookup_map` and exact prefix matching against `PANEL_SPECS_BY_ID` or exact `target_polymer_seq`.
    - Raises `ValueError` when `require_known=True` and the target does not map to a registered `REF_*` ID.
    """
    panel_map = curated_panel_by_id if curated_panel_by_id is not None else PANEL_SPECS_BY_ID
    raw = str(raw_target_id or "").strip()
    if not raw:
        if require_known:
            raise ValueError("Empty curated reference target header")
        return "unresolved_reference"

    if db_lookup_map and raw in db_lookup_map:
        raw = str(db_lookup_map[raw]).strip()

    # Strip path prefix if present
    tok = Path(raw).name
    # Strip | metadata suffix (`REF_ID|uniprot:...|role:...`)
    tok = tok.split("|", 1)[0]
    # Strip __PDB_ metadata suffix (`REF_ID__PDB_...`)
    tok = tok.split("__PDB_", 1)[0]

    if tok in panel_map:
        return tok
    if tok.endswith(".pdb") and tok[:-4] in panel_map:
        return tok[:-4]

    # Handle Foldseek chain suffix appended after .pdb (e.g. REF_PIL_8REY_CdiPB.pdb_A or REF_PIL_8REY_CdiPB_A)
    if ".pdb_" in tok:
        stem = tok.split(".pdb_", 1)[0]
        if stem in panel_map:
            return stem
    if "_" in tok:
        stem = tok.rsplit("_", 1)[0]
        if stem in panel_map:
            return stem

    # Check exact polymer sequence match if provided
    t_poly = clean_aa_sequence(target_polymer_seq)
    if t_poly:
        for cand_id, cand_spec in panel_map.items():
            p_s = clean_aa_sequence(str(cand_spec.get("pdb_entity_sequence", "")))
            u_s = clean_aa_sequence(str(cand_spec.get("uniprot_sequence", "")))
            if (p_s and p_s != "NONE" and t_poly == p_s) or (u_s and t_poly == u_s):
                return str(cand_id)

    if require_known:
        raise ValueError(f"Curated reference target header {raw_target_id!r} could not be mapped to PANEL_SPECS_BY_ID")
    return "unresolved_reference"


def qualify_structural_hit_for_member(
    hit: dict[str, Any],
    *,
    member_gene: str,
    member_cat_id: str,
    member_aa_seq: str | None,
    model_info: dict[str, Any] | None = None,
    model_source: str | None = None,
    model_polymer_aa_seq: str | None = None,
    model_plddt_0to1: float | None = None,
    min_plddt: float = 0.70,
    min_aa_len: int = 100,
    min_qtm: float = 0.50,
    min_cov: float = 0.50,
) -> dict[str, Any]:
    """Qualify an individual structural hit by exact member-to-model sequence identity and source-specific pLDDT (`EXP4-PIPELINE-REVIEW.txt` item 2).

    A matching AF2 model on one locus CANNOT validate a mismatched ESMFold hit or a cluster-propagated
    hit on another family member whose sequence differs from the structural model polymer.
    """
    if model_info is not None:
        if model_source is None:
            model_source = str(model_info.get("model_source") or hit.get("model_source") or "unknown")
        if model_polymer_aa_seq is None:
            model_polymer_aa_seq = model_info.get("model_polymer_aa_seq")
        if model_plddt_0to1 is None:
            model_plddt_0to1 = model_info.get("plddt_0to1")
    if model_source is None:
        model_source = str(hit.get("model_source") or "unknown")

    seq_clean = clean_aa_sequence(member_aa_seq)
    poly_check = verify_prediction_sequence_hash_match(seq_clean, model_polymer_aa_seq)
    exact_match = bool(poly_check["is_valid_exact_match"])
    plddt_val = float(model_plddt_0to1) if model_plddt_0to1 is not None and not pd.isna(model_plddt_0to1) else np.nan
    plddt_ok = bool(np.isfinite(plddt_val) and plddt_val >= min_plddt)
    len_ok = len(seq_clean) >= min_aa_len

    qtm = float(hit.get("qtmscore", np.nan))
    qc = float(hit.get("qcov", np.nan))
    tc = float(hit.get("tcov", np.nan))
    cov_tm_ok = bool(
        np.isfinite(qtm) and qtm >= min_qtm and np.isfinite(qc) and qc >= min_cov and np.isfinite(tc) and tc >= min_cov
    )

    reasons: list[str] = []
    if not exact_match:
        reasons.append(f"rejected_model_sequence_mismatch_or_missing({poly_check['sequence_match_status']})")
    if not plddt_ok:
        reasons.append(f"rejected_low_or_missing_model_plddt_below_0.70({plddt_val})")
    if not len_ok:
        reasons.append(f"rejected_short_protein_lt{min_aa_len}({len(seq_clean)})")
    if not cov_tm_ok:
        reasons.append(f"rejected_low_or_missing_tm_or_coverage(qtm={qtm},qcov={qc},tcov={tc})")

    qualified = len(reasons) == 0
    out = dict(hit)
    out.update(
        {
            "member_gene": member_gene,
            "member_cat_id": member_cat_id,
            "model_source": model_source,
            "model_exact_seq_match": exact_match,
            "model_polymer_hash_status": poly_check["sequence_match_status"],
            "model_plddt_0to1": plddt_val,
            "member_aa_length": len(seq_clean),
            "qualified_for_role_evidence": qualified,
            "rejection_reason": "none" if qualified else ";".join(reasons),
        }
    )
    return out


def resolve_pdb_chain_and_architecture_guards(
    pdb_target: str | None,
    query_aa_seq: str | None,
    qtmscore: float | None,
    ttmscore: float | None,
    qcov: float | None,
    tcov: float | None,
    *,
    target_polymer_aa_seq: str | None = None,
    curated_panel_by_id: dict[str, dict[str, Any]] | None = None,
    verified_n_terminal_helix: bool = False,
    verified_alignment_active_site_conserved: bool = False,
    has_verified_operon_context: bool = False,
) -> dict[str, Any]:
    """Enforce strict chain-level and architecture guards on structural hits (`EXP4-CODE-REVIEW.txt` items 1, 4, 7):
    - NEVER infers `8CE8` or any PDB chain identity from `qtmscore / ttmscore` ratios!
    - Resolves curated `REF_*` identifiers directly through `curated_panel_by_id` (`PANEL_SPECS_BY_ID`),
      or resolves raw PDB targets ONLY when `target_polymer_aa_seq` matches the exact sequence/SHA256
      of a curated chain entity. Absent exact sequence/hash match or `REF_*` ID, `8CE8` chain E
      remains `unresolved_8CE8_chain_E_label_vs_auth_ambiguity` with `structural_fold_compatible = False`
      and `is_functional_role_positive = False`.
    - Requires BOTH `qcov` AND `tcov` to be finite and above threshold (`np.isfinite(tc) and tc >= min_cov`);
      missing `tcov` (`NaN` or `None`) NEVER passes a full-architecture gate!
    - Strictly separates `structural_fold_compatible` (whole-chain fold/coverage compatibility) from
      `has_full_mechanism_support` / `is_functional_role_positive` (which requires verified active-site/topology/operon
      evidence, never coverage alone).
    """
    panel_map = curated_panel_by_id if curated_panel_by_id is not None else PANEL_SPECS_BY_ID
    t_raw = str(pdb_target or "missing").strip()
    seq = clean_aa_sequence(query_aa_seq)
    qlen = len(seq)
    qtm = float(qtmscore) if qtmscore is not None and not pd.isna(qtmscore) else np.nan
    ttm = float(ttmscore) if ttmscore is not None and not pd.isna(ttmscore) else np.nan
    qc = float(qcov) if qcov is not None and not pd.isna(qcov) else np.nan
    tc = float(tcov) if tcov is not None and not pd.isna(tcov) else np.nan
    t_poly = clean_aa_sequence(target_polymer_aa_seq)
    t_poly_sha = hashlib.sha256(t_poly.encode("ascii", errors="replace")).hexdigest() if t_poly else "missing"

    # Match against canonical REF_* ID (including Foldseek target names like REF_PIL_8REY_CdiPB.pdb or REF_PIL_8REY_CdiPB__PDB_...)
    matched_ref_spec: dict[str, Any] | None = None
    base_tok = t_raw.split("__PDB_")[0]
    if base_tok.endswith(".pdb"):
        base_tok = base_tok[:-4]
    if base_tok in panel_map:
        matched_ref_spec = panel_map[base_tok]
    elif t_poly:
        # Check exact polymer sequence / SHA256 against curated panel entries
        for cand_spec in panel_map.values():
            u_s = clean_aa_sequence(str(cand_spec.get("uniprot_sequence", "")))
            p_s = clean_aa_sequence(str(cand_spec.get("pdb_entity_sequence", "")))
            if (p_s and p_s != "NONE" and t_poly == p_s) or (u_s and t_poly == u_s):
                matched_ref_spec = cand_spec
                break

    t_low = t_raw.lower()

    # 1. Handle canonical curated REF_* or exact-sequence-matched panel entries
    if matched_ref_spec is not None:
        ref_id = str(matched_ref_spec["ref_id"])
        ctrl_type = str(matched_ref_spec.get("control_type", ""))
        pdb_id = str(matched_ref_spec.get("pdb_id", "none"))
        ent_id = str(matched_ref_spec.get("pdb_entity_id", "none"))
        l_id = str(matched_ref_spec.get("pdb_label_asym_id", "none"))
        a_id = str(matched_ref_spec.get("pdb_auth_asym_id", "none"))
        u_id = str(matched_ref_spec.get("uniprot_acc", "none"))
        exp_roles = str(matched_ref_spec.get("experimentally_supported_roles", "none"))

        if ctrl_type == "negative_propeptide_chain_control" or ref_id == "REF_CTRL_9UES_ChainA_Propeptide_NegativeCtrl":
            return {
                "pdb_guard_status": "rejected_9UES_chainA_propeptide_negative_architecture_control",
                "resolved_ref_id": ref_id,
                "resolved_pdb_id": pdb_id,
                "resolved_entity_id": ent_id,
                "resolved_label_asym_id": l_id,
                "resolved_auth_asym_id": a_id,
                "resolved_uniprot_acc": u_id,
                "structural_fold_compatible": False,
                "has_full_mechanism_support": False,
                "is_functional_role_positive": False,
                "allowed_role_transfer": "none_negative_architecture_control",
                "guard_explanation": (
                    "9UES chain A (Entity 1, 70 aa, UniProt P58502) is the cleaved Tk-subtilisin inhibitory propeptide, "
                    "NOT the catalytic protease; blocked from peptidase functional transfer."
                ),
            }

        if ctrl_type.startswith("pdb_chain_ambiguity_control"):
            both_cov_ok = bool(np.isfinite(qtm) and qtm >= 0.50 and np.isfinite(qc) and qc >= 0.50 and np.isfinite(tc) and tc >= 0.50 and qlen >= 100)
            mech_ok = bool(both_cov_ok and verified_alignment_active_site_conserved and has_verified_operon_context and ent_id == "5")
            return {
                "pdb_guard_status": (
                    "qualified_8CE8_Entity5_CcmE_with_verified_alignment_and_operon"
                    if mech_ok
                    else f"blocked_automatic_CcmE_transfer__resolved_{ref_id}_fold_only"
                ),
                "resolved_ref_id": ref_id,
                "resolved_pdb_id": pdb_id,
                "resolved_entity_id": ent_id,
                "resolved_label_asym_id": l_id,
                "resolved_auth_asym_id": a_id,
                "resolved_uniprot_acc": u_id,
                "structural_fold_compatible": both_cov_ok,
                "has_full_mechanism_support": mech_ok,
                "is_functional_role_positive": mech_ok,
                "allowed_role_transfer": (
                    "heme_chaperone_CcmE" if mech_ok else "none_generic_fold_only_no_automatic_CcmE_transfer"
                ),
                "guard_explanation": (
                    f"Resolved by explicit reference chain {ref_id} (PDB {pdb_id} Entity {ent_id}, label={l_id}, auth={a_id}, {u_id}). "
                    f"fold_compatible={both_cov_ok}, verified_active_site={verified_alignment_active_site_conserved}, "
                    f"verified_operon={has_verified_operon_context}. Automatic CcmE transfer is blocked."
                ),
            }

        if ref_id == "REF_PEP_9UES_ChainB_TK1689_Protease":
            both_cov_ok = bool(
                np.isfinite(qtm) and qtm >= 0.50 and np.isfinite(qc) and qc >= 0.50 and np.isfinite(tc) and tc >= 0.50 and qlen >= 250
            )
            mech_ok = bool(both_cov_ok and verified_alignment_active_site_conserved)
            return {
                "pdb_guard_status": (
                    "qualified_9UES_chainB_TK1689_protease_with_verified_active_site"
                    if mech_ok
                    else (
                        "structural_fold_compatible_9UES_chainB_mutant_S359A_active_site_unverified"
                        if both_cov_ok
                        else "rejected_9UES_chainB_missing_or_low_coverage"
                    )
                ),
                "resolved_ref_id": ref_id,
                "resolved_pdb_id": pdb_id,
                "resolved_entity_id": ent_id,
                "resolved_label_asym_id": l_id,
                "resolved_auth_asym_id": a_id,
                "resolved_uniprot_acc": u_id,
                "structural_fold_compatible": both_cov_ok,
                "has_full_mechanism_support": mech_ok,
                "is_functional_role_positive": mech_ok,
                "allowed_role_transfer": "nutrient_acquisition_peptidase" if mech_ok else "fold_compatibility_only_not_verified_enzyme",
                "guard_explanation": (
                    "9UES chain B (Entity 2, 427 aa, Q5JIZ5 TK1689) carries inactive S359A mutation; "
                    "requires both qcov>=0.50 and tcov>=0.50 for fold compatibility plus verified catalytic triad for mechanism support."
                ),
            }

        if ref_id in ("REF_PIL_8REY_CdiPB", "REF_PIL_8RH5_ArlB", "REF_PIL_8Q30_AapB", "REF_PIL_5O4U_FlaB0"):
            both_cov_ok = bool(
                np.isfinite(qtm) and qtm >= 0.50 and np.isfinite(qc) and qc >= 0.60 and np.isfinite(tc) and tc >= 0.60
            )
            arch_and_topo_ok = bool(both_cov_ok and verified_n_terminal_helix)
            mech_ok = bool(arch_and_topo_ok and has_verified_operon_context)
            return {
                "pdb_guard_status": (
                    f"convergent_{ref_id}_architecture_topology_and_machinery"
                    if mech_ok
                    else (
                        f"structural_and_n_term_helix_compatible_{ref_id}_machinery_context_required_outside_Thermoplasmatales"
                        if arch_and_topo_ok
                        else f"rejected_{ref_id}_missing_tcov_or_n_terminal_helix_architecture"
                    )
                ),
                "resolved_ref_id": ref_id,
                "resolved_pdb_id": pdb_id,
                "resolved_entity_id": ent_id,
                "resolved_label_asym_id": l_id,
                "resolved_auth_asym_id": a_id,
                "resolved_uniprot_acc": u_id,
                "structural_fold_compatible": arch_and_topo_ok,
                "has_full_mechanism_support": mech_ok,
                "is_functional_role_positive": mech_ok,
                "allowed_role_transfer": (
                    "surface_assembly_pilin_archaellin" if mech_ok else "architecture_candidate_only_not_confirmed_pilin"
                ),
                "guard_explanation": (
                    f"{ref_id}: requires finite qcov>=0.60 and tcov>=0.60 plus verified N-terminal helix for structural compatibility, "
                    "and assembly machinery context for mechanism support. Virus entry is untested."
                ),
            }

        if ctrl_type == "viral_source_protein_reference":
            both_cov_ok = bool(np.isfinite(qtm) and qtm >= 0.50 and np.isfinite(qc) and qc >= 0.50 and np.isfinite(tc) and tc >= 0.50)
            return {
                "pdb_guard_status": (
                    f"viral_virion_fold_match_{ref_id}_requires_proviral_origin_audit_not_host_receptor"
                    if both_cov_ok
                    else f"rejected_{ref_id}_low_or_missing_coverage"
                ),
                "resolved_ref_id": ref_id,
                "resolved_pdb_id": pdb_id,
                "resolved_entity_id": ent_id,
                "resolved_label_asym_id": l_id,
                "resolved_auth_asym_id": a_id,
                "resolved_uniprot_acc": u_id,
                "structural_fold_compatible": both_cov_ok,
                "has_full_mechanism_support": False,
                "is_functional_role_positive": False,
                "allowed_role_transfer": "none_viral_virion_protein_not_cellular_host_receptor",
                "guard_explanation": (
                    f"{ref_id} is a viral-encoded virion protein; structural match in a MAG requires proviral/mobile origin audit "
                    "and does NOT establish cellular host receptor or defense function."
                ),
            }

        both_cov_ok = bool(np.isfinite(qtm) and qtm >= 0.50 and np.isfinite(qc) and qc >= 0.50 and np.isfinite(tc) and tc >= 0.50 and qlen >= 100)
        mech_ok = bool(both_cov_ok and verified_alignment_active_site_conserved and has_verified_operon_context and exp_roles != "none")
        return {
            "pdb_guard_status": (
                f"convergent_{ref_id}_fold_and_mechanism_supported"
                if mech_ok
                else (f"structural_fold_compatible_{ref_id}_only" if both_cov_ok else f"rejected_{ref_id}_low_or_missing_coverage")
            ),
            "resolved_ref_id": ref_id,
            "resolved_pdb_id": pdb_id,
            "resolved_entity_id": ent_id,
            "resolved_label_asym_id": l_id,
            "resolved_auth_asym_id": a_id,
            "resolved_uniprot_acc": u_id,
            "structural_fold_compatible": both_cov_ok,
            "has_full_mechanism_support": mech_ok,
            "is_functional_role_positive": mech_ok,
            "allowed_role_transfer": exp_roles if mech_ok else "fold_compatibility_only",
            "guard_explanation": f"Curated reference {ref_id}: both qcov and tcov required for fold compatibility; mechanism support requires alignment and context.",
        }

    # 2. Raw non-curated PDB IDs: check 9UES and 8CE8 explicitly without qtm/ttm ratio inference!
    if t_low.startswith("9ues"):
        chain_tok = t_raw.split("_")[-1].upper() if "_" in t_raw else ""
        if chain_tok == "A":
            return resolve_pdb_chain_and_architecture_guards(
                "REF_CTRL_9UES_ChainA_Propeptide_NegativeCtrl",
                query_aa_seq,
                qtmscore,
                ttmscore,
                qcov,
                tcov,
                curated_panel_by_id=panel_map,
            )
        if chain_tok == "B":
            return resolve_pdb_chain_and_architecture_guards(
                "REF_PEP_9UES_ChainB_TK1689_Protease",
                query_aa_seq,
                qtmscore,
                ttmscore,
                qcov,
                tcov,
                curated_panel_by_id=panel_map,
                verified_alignment_active_site_conserved=verified_alignment_active_site_conserved,
                has_verified_operon_context=has_verified_operon_context,
            )

    if t_low.startswith("8ce8"):
        # Without exact target polymer sequence hash matching Entity 1 (P33931) or Entity 5 (P69490),
        # 8ce8-assembly1_E is strictly UNRESOLVED (never inferred from qtm/ttm ratio!).
        return {
            "pdb_guard_status": "unresolved_8CE8_chain_E_label_vs_auth_ambiguity_no_target_sequence_hash",
            "resolved_ref_id": "none_unresolved_8CE8_chain_E",
            "resolved_pdb_id": "8CE8",
            "resolved_entity_id": "unresolved_1_or_5",
            "resolved_label_asym_id": "unresolved_E_or_I",
            "resolved_auth_asym_id": "unresolved_a_or_E",
            "resolved_uniprot_acc": "unresolved_P33931_or_P69490",
            "structural_fold_compatible": False,
            "has_full_mechanism_support": False,
            "is_functional_role_positive": False,
            "allowed_role_transfer": "none_unresolved_chain_ambiguity_no_CcmE_transfer",
            "guard_explanation": (
                "8CE8 target E is ambiguous between PDB label_asym_id=E (Entity 1 CcmA, P33931, 218 aa) and "
                "auth_asym_id=E (Entity 5 CcmE, P69490, 159 aa). TM-score ratios cannot resolve chain identity; "
                "requires direct alignment against explicit chain coordinate files REF_CTRL_8CE8_Entity1/Entity5. "
                "No automatic CcmE functional transfer allowed."
            ),
        }

    # Generic raw PDB hit: requires finite qtm>=0.50, qcov>=0.50, AND tcov>=0.50 for fold compatibility; NEVER sets is_functional_role_positive=True
    both_cov_ok = bool(
        np.isfinite(qtm) and qtm >= 0.50 and np.isfinite(qc) and qc >= 0.50 and np.isfinite(tc) and tc >= 0.50 and qlen >= 100
    )
    return {
        "pdb_guard_status": "generic_pdb_fold_compatible_only" if both_cov_ok else "unqualified_or_missing_coverage_pdb_hit",
        "resolved_ref_id": "none",
        "resolved_pdb_id": t_raw.split("-")[0].split("_")[0].upper() if t_raw != "missing" else "none",
        "resolved_entity_id": "unknown",
        "resolved_label_asym_id": "unknown",
        "resolved_auth_asym_id": t_raw.split("_")[-1] if "_" in t_raw else "unknown",
        "resolved_uniprot_acc": "unknown",
        "structural_fold_compatible": both_cov_ok,
        "has_full_mechanism_support": False,
        "is_functional_role_positive": False,
        "allowed_role_transfer": "fold_compatibility_only_not_specific_function",
        "guard_explanation": "Whole-fold similarity alone does not establish specific biochemical function without curated chain, active-site, and context verification.",
    }


def classify_flanking_machinery_guilds_detailed(
    desc: str | None,
    profiles: list[str] | None = None,
) -> dict[str, Any]:
    """Classify a flanking gene's functional guild(s) using bounded regex tokens and explicit profile IDs (`EXP4-CODE-REVIEW.txt` item 5).

    - Prevents substring collisions (`opp` in `copper`, `cas1` in unrelated tokens, generic `modification` or `glycosyltransferase` or `integrase`).
    - Explicitly separates `host_defense_antiviral` (`CRISPR-Cas`, restriction endonuclease, abortive infection)
      from `viral_or_mobile_element_origin` (capsid, portal, terminase, virion turret).
    """
    guilds: set[str] = set()
    matched_tokens: dict[str, list[str]] = {}
    d_str = str(desc or "").strip()
    p_list = [str(p).strip() for p in (profiles or []) if str(p).strip()]
    combined = f"{d_str} {' '.join(p_list)}".strip()
    if not combined or combined.lower() in ("missing", "nan", "hypothetical protein", "uncharacterized protein"):
        return {"guilds": guilds, "matched_tokens": matched_tokens}

    guilds.add("generic_any_known_annotation")
    for g_name, pat in BOUNDED_GUILD_PATTERNS.items():
        hits = sorted({m.group(0) for m in pat.finditer(combined)})
        if hits:
            guilds.add(g_name)
            matched_tokens[g_name] = hits

    return {"guilds": guilds, "matched_tokens": matched_tokens}


def classify_flanking_machinery_guilds(desc: str | None, profiles: list[str] | None = None) -> set[str]:
    """Return the set of bounded flanking guilds for a gene."""
    return classify_flanking_machinery_guilds_detailed(desc, profiles)["guilds"]


def evaluate_competing_roles_for_family(
    *,
    family_id: str,
    candidate_role: str,
    is_dark: bool,
    surface_localization_class: str,
    frac_exported: float,
    curated_seq_hits: list[dict[str, Any]],
    curated_struct_hits: list[dict[str, Any]],
    own_domain_descriptions: list[str],
    own_c6_hits: list[str],
    pdb_guard_result: dict[str, Any],
    role_null_results: dict[str, dict[str, Any]],
    tax_classes_present: set[str],
) -> dict[str, Any]:
    """Score competing biological roles simultaneously without forcing mutually exclusive single-class labels.

    Enforces `EXP4-EVIDENCE-REVIEW.txt` and `EXP4-CODE-REVIEW.txt`:
    - Only `experimentally_supported_roles` on curated positive references can contribute to positive reference
      homolog/structural evidence; `hypothesized_roles` (`untested`) NEVER count as positive reference hits.
    - Viral-source virion proteins (`viral_source_entry_or_capsid_protein`: `4KFR`, `2W77`, `3J9X`) are flagged
      separately (`viral_virion_homolog_requires_proviral_origin_audit`) and NEVER establish cellular host receptor/defense!
    - Requires BOTH finite `qcov >= 0.50` AND `tcov >= 0.50` on structural and sequence reference hits.
    - Distinguishes structural/domain compatibility (`tier3_structural_or_domain_compatible_only`) from
      convergent mechanism support (`tier1`/`tier2`).
    """
    dom_text = " ".join(own_domain_descriptions + own_c6_hits)

    seq_by_role: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    struct_by_role: dict[str, list[dict[str, Any]]] = collections.defaultdict(list)
    seq_by_ref: dict[str, dict[str, Any]] = {}
    struct_by_ref: dict[str, dict[str, Any]] = {}
    viral_source_hits: list[str] = []
    assembly_machinery_hits: list[str] = []

    for h in curated_seq_hits:
        rid = str(h.get("ref_id", ""))
        spec = PANEL_SPECS_BY_ID.get(rid, {})
        ctrl_type = str(h.get("control_type") or spec.get("control_type", ""))
        exp_roles_str = str(h.get("experimentally_supported_roles") or spec.get("experimentally_supported_roles", "none"))
        qc = float(h.get("qcov", np.nan))
        tc = float(h.get("tcov", np.nan))
        ev = float(h.get("evalue", np.nan))
        if not (np.isfinite(qc) and qc >= 0.50 and np.isfinite(tc) and tc >= 0.50 and np.isfinite(ev) and ev <= 1e-5):
            continue
        seq_by_ref[rid] = h
        if ctrl_type in (
            "negative_propeptide_chain_control",
            "pdb_chain_ambiguity_control_authE_labelI",
            "pdb_chain_ambiguity_control_labelE_autha",
        ):
            continue
        if ctrl_type == "assembly_machinery_reference":
            # Per EXP4-AUDIT-FINAL.txt Item 4: ArlI/J are assembly machinery context ONLY, excluded from direct pilin positive
            assembly_machinery_hits.append(f"seq:{rid}")
            continue
        if ctrl_type == "viral_source_protein_reference" or exp_roles_str == "viral_source_entry_or_capsid_protein":
            viral_source_hits.append(f"seq:{rid}")
            continue
        for er in exp_roles_str.split(";"):
            er = er.strip()
            if er in ROLE_CLASSES:
                seq_by_role[er].append(h)

    n_struct_hits_qualified = 0
    n_struct_hits_rejected_gate = 0

    for h in curated_struct_hits:
        if h.get("qualified_for_role_evidence") is False or h.get("model_exact_seq_match") is False:
            n_struct_hits_rejected_gate += 1
            continue
        if "model_plddt_0to1" in h:
            pl_v = float(h.get("model_plddt_0to1", np.nan))
            if not (np.isfinite(pl_v) and pl_v >= 0.70):
                n_struct_hits_rejected_gate += 1
                continue
        rid = str(h.get("ref_id", ""))
        spec = PANEL_SPECS_BY_ID.get(rid, {})
        ctrl_type = str(h.get("control_type") or spec.get("control_type", ""))
        exp_roles_str = str(h.get("experimentally_supported_roles") or spec.get("experimentally_supported_roles", "none"))
        qc = float(h.get("qcov", np.nan))
        tc = float(h.get("tcov", np.nan))
        qtm = float(h.get("qtmscore", np.nan))
        if not (np.isfinite(qc) and qc >= 0.50 and np.isfinite(tc) and tc >= 0.50 and np.isfinite(qtm) and qtm >= 0.50):
            n_struct_hits_rejected_gate += 1
            continue
        n_struct_hits_qualified += 1
        struct_by_ref[rid] = h
        if ctrl_type in (
            "negative_propeptide_chain_control",
            "pdb_chain_ambiguity_control_authE_labelI",
            "pdb_chain_ambiguity_control_labelE_autha",
        ):
            continue
        if ctrl_type == "assembly_machinery_reference":
            # Per EXP4-AUDIT-FINAL.txt Item 4: ArlI/J are assembly machinery context ONLY, excluded from direct pilin positive
            assembly_machinery_hits.append(f"struct:{rid}")
            continue
        if ctrl_type == "viral_source_protein_reference" or exp_roles_str == "viral_source_entry_or_capsid_protein":
            viral_source_hits.append(f"struct:{rid}")
            continue
        for er in exp_roles_str.split(";"):
            er = er.strip()
            if er in ROLE_CLASSES:
                struct_by_role[er].append(h)

    # Check confusable OppA vs sugar-binder specificity when substrate_binding_or_uptake is triggered
    true_oppa_score = max(
        [float(seq_by_ref[r]["bits"]) for r in ("REF_SBP_SSO1288_TrueOppA", "REF_SBP_APE1583_TrueOppA") if r in seq_by_ref]
        + [0.0]
    )
    sugar_ctrl_score = max(
        [
            float(seq_by_ref[r]["bits"])
            for r in ("REF_SBP_SSO3043_SugarConfusableCtrl", "REF_SBP_1EU8_TrehaloseMaltoseCtrl", "REF_SBP_2B3B_GlucoseBinderCtrl")
            if r in seq_by_ref
        ]
        + [
            float(struct_by_ref[r]["qtmscore"]) * 100.0
            for r in ("REF_SBP_1EU8_TrehaloseMaltoseCtrl", "REF_SBP_2B3B_GlucoseBinderCtrl")
            if r in struct_by_ref
        ]
        + [0.0]
    )
    has_sbp_dom = bool(BOUNDED_GUILD_PATTERNS["substrate_binding_or_uptake"].search(dom_text))
    if true_oppa_score > 0 and true_oppa_score >= max(sugar_ctrl_score * 1.15, 50.0):
        sbp_specificity_status = "favours_true_oligopeptide_binder_over_sugar_control"
    elif true_oppa_score > 0 or sugar_ctrl_score > 0 or has_sbp_dom:
        sbp_specificity_status = "generic_sbp_fold_substrate_unresolved_peptide_vs_sugar"
    else:
        sbp_specificity_status = "no_sbp_match"

    role_tiers: dict[str, str] = {}
    role_scores: dict[str, float] = {}
    role_details: dict[str, str] = {}
    supported_overlapping_roles: list[str] = []

    for r_name in ROLE_CLASSES:
        s_hits = seq_by_role.get(r_name, [])
        st_hits = struct_by_role.get(r_name, [])
        null_key = "host_defense_antiviral" if r_name == "virus_interaction_or_entry" else r_name
        null_info = role_null_results.get(
            null_key,
            role_null_results.get("s_layer_or_charted_coat" if null_key == "cell_coat_or_envelope" else null_key, {}),
        )
        obs_adj = int(null_info.get("obs_adjacent", 0))
        p_emp = float(null_info.get("p_exact_one_sided", null_info.get("p_empirical_one_sided", 1.0)))
        q_bh = float(null_info.get("q_bh", 1.0))
        ctx_sig = bool(obs_adj >= 2 and p_emp <= 0.05)

        pat_key = "host_defense_antiviral" if r_name == "virus_interaction_or_entry" else r_name
        has_dom_kw = bool(BOUNDED_GUILD_PATTERNS[pat_key].search(dom_text))

        if r_name == "nutrient_acquisition_peptidase" and pdb_guard_result.get("pdb_guard_status", "").startswith("rejected_9UES_chainA"):
            st_hits = [h for h in st_hits if "9UES" not in str(h.get("ref_id", ""))]

        if r_name == "surface_assembly_pilin_archaellin":
            only_non_thermoplasmatales = "c__Thermoplasmata" not in tax_classes_present
            has_prepilin_topo = surface_localization_class == "prepilin_or_archaellin_filament_subunit_candidate"
            if s_hits and (ctx_sig or (not only_non_thermoplasmatales and has_prepilin_topo)):
                tier = "tier1_experimentally_validated_reference_and_context"
                score = 0.85
            elif (st_hits or has_prepilin_topo) and ctx_sig:
                tier = "tier2_qualified_architecture_and_machinery_context"
                score = 0.70
            elif s_hits or st_hits or has_dom_kw:
                tier = "tier3_structural_or_domain_compatible_only_ambiguous"
                score = 0.45
            elif has_prepilin_topo or ctx_sig:
                tier = "tier4_localization_or_context_only"
                score = 0.25
            else:
                tier = "none"
                score = 0.0
        elif r_name == "substrate_binding_or_uptake":
            # Exclude sugar-control-only hits from counting as positive peptide-binding support
            pos_s_hits = [h for h in s_hits if "Ctrl" not in str(h.get("ref_id", ""))]
            if pos_s_hits and sbp_specificity_status == "favours_true_oligopeptide_binder_over_sugar_control" and ctx_sig:
                tier = "tier1_experimentally_validated_reference_and_context"
                score = 0.85
            elif (pos_s_hits or st_hits or has_dom_kw) and ctx_sig:
                tier = "tier2_qualified_architecture_and_machinery_context"
                score = 0.65
            elif s_hits or st_hits or has_dom_kw:
                tier = "tier3_structural_or_domain_compatible_only_ambiguous"
                score = 0.40
            elif ctx_sig:
                tier = "tier4_localization_or_context_only"
                score = 0.25
            else:
                tier = "none"
                score = 0.0
        elif r_name == "virus_interaction_or_entry":
            # Viral-source virion protein hits DO NOT establish cellular host receptor/defense; only report them as proviral-audit leads
            if has_dom_kw and ctx_sig:
                tier = "tier2_qualified_architecture_and_machinery_context"
                score = 0.65
            elif has_dom_kw:
                tier = "tier3_structural_or_domain_compatible_only_ambiguous"
                score = 0.40
            elif ctx_sig:
                tier = "tier4_localization_or_context_only"
                score = 0.25
            else:
                tier = "none"
                score = 0.0
        else:
            if s_hits and (ctx_sig or has_dom_kw):
                tier = "tier1_experimentally_validated_reference_and_context"
                score = 0.85
            elif (st_hits or has_dom_kw) and ctx_sig:
                tier = "tier2_qualified_architecture_and_machinery_context"
                score = 0.70
            elif s_hits or st_hits or has_dom_kw:
                tier = "tier3_structural_or_domain_compatible_only_ambiguous"
                score = 0.45
            elif ctx_sig:
                tier = "tier4_localization_or_context_only"
                score = 0.25
            else:
                tier = "none"
                score = 0.0

        role_tiers[r_name] = tier
        role_scores[r_name] = score
        best_seq_id = s_hits[0]["ref_id"] if s_hits else "none"
        best_st_id = st_hits[0]["ref_id"] if st_hits else "none"
        role_details[r_name] = f"tier={tier};seq_ref={best_seq_id};struct_ref={best_st_id};ctx_adj={obs_adj}(p={p_emp:.3g},q={q_bh:.3g})"
        if score >= 0.40:
            supported_overlapping_roles.append(r_name)

    max_score = max(role_scores.values()) if role_scores else 0.0
    top_roles = [r for r, sc in role_scores.items() if sc == max_score and sc > 0.0]

    viral_origin_null = role_null_results.get("viral_or_mobile_element_origin", {})
    viral_flank_adj = int(viral_origin_null.get("obs_adjacent", 0))
    if viral_source_hits:
        viral_origin_audit_status = (
            f"viral_virion_protein_similarity({'+'.join(viral_source_hits)};viral_flank_adj={viral_flank_adj})_requires_proviral_contig_origin_audit_not_host_receptor"
        )
    elif viral_flank_adj > 0:
        viral_origin_audit_status = f"flanking_viral_or_mobile_hallmarks_present(adj={viral_flank_adj})"
    else:
        viral_origin_audit_status = "no_viral_virion_similarity_or_hallmark_flank"

    if max_score >= 0.65 and len(supported_overlapping_roles) == 1:
        disc_status = f"convergent_role_candidate:{supported_overlapping_roles[0]}"
    elif max_score >= 0.65 and len(supported_overlapping_roles) > 1:
        disc_status = f"overlapping_multi_role_candidate:{'+'.join(supported_overlapping_roles)}"
    elif max_score >= 0.40:
        disc_status = f"ambiguous_fold_or_domain_level_candidate:{'+'.join(supported_overlapping_roles)}"
    elif viral_source_hits:
        disc_status = "putative_viral_or_proviral_structural_protein_not_cellular_surface_receptor"
    elif pd.notna(frac_exported) and frac_exported >= 0.50:
        disc_status = f"uncertain_exported_{surface_localization_class}_insufficient_functional_evidence"
    else:
        disc_status = "uncertain_non_exported_or_insufficient_evidence"

    return {
        "family_id": family_id,
        "candidate_role": candidate_role,
        "is_dark": is_dark,
        "surface_localization_class": surface_localization_class,
        "sbp_substrate_specificity_status": sbp_specificity_status,
        "viral_origin_audit_status": viral_origin_audit_status,
        "assembly_machinery_hits": ";".join(assembly_machinery_hits) if assembly_machinery_hits else "none",
        "has_assembly_machinery_context": bool(assembly_machinery_hits),
        "assembly_machinery_context_status": (
            f"assembly_machinery_detected({'+'.join(assembly_machinery_hits)})"
            if assembly_machinery_hits
            else "none"
        ),
        "pdb_guard_status": pdb_guard_result.get("pdb_guard_status", "none"),
        "pdb_guard_explanation": pdb_guard_result.get("guard_explanation", "none"),
        "n_struct_hits_qualified": int(n_struct_hits_qualified),
        "n_struct_hits_rejected_gate": int(n_struct_hits_rejected_gate),
        "supported_overlapping_roles": ";".join(supported_overlapping_roles) if supported_overlapping_roles else "none",
        "n_supported_overlapping_roles": len(supported_overlapping_roles),
        "top_competing_roles": ";".join(top_roles) if top_roles else "none",
        "max_role_evidence_score": float(max_score),
        "exploratory_role_ranking_score": float(max_score),
        "role_discrimination_outcome": disc_status,
        **{f"tier_{r}": role_tiers[r] for r in ROLE_CLASSES},
        **{f"score_{r}": role_scores[r] for r in ROLE_CLASSES},
        **{f"detail_{r}": role_details[r] for r in ROLE_CLASSES},
    }


def apply_bh_fdr(p_values: list[float] | np.ndarray) -> np.ndarray:
    """Compute Benjamini-Hochberg FDR q-values preserving NaN for missing p-values."""
    arr = np.asarray(p_values, dtype="float64")
    q_out = np.full_like(arr, np.nan, dtype="float64")
    valid_mask = np.isfinite(arr)
    if not np.any(valid_mask):
        return q_out
    p_valid = np.clip(arr[valid_mask], 0.0, 1.0)
    m = len(p_valid)
    order = np.argsort(p_valid)
    ranked_p = p_valid[order]
    ranks = np.arange(1, m + 1, dtype="float64")
    q_raw = ranked_p * (m / ranks)
    q_mono = np.minimum.accumulate(q_raw[::-1])[::-1]
    q_clipped = np.clip(q_mono, 0.0, 1.0)
    q_valid = np.empty_like(p_valid)
    q_valid[order] = q_clipped
    q_out[valid_mask] = q_valid
    return q_out


def evaluate_gspa_heldout_recovery_and_leakage_guard(
    claims_df: pd.DataFrame,
    ground_truth_df: pd.DataFrame,
    *,
    holdout_col: str = "family_id",
    score_threshold: float = 0.50,
    training_split_manifest_verified: bool = False,
    disjoint_train_manifest_sha256: str | None = None,
    go_ancestors_map: dict[str, set[str]] | None = None,
) -> pd.DataFrame:
    """Evaluate GSPA / DeepGO-PlusPlus-Light predictions on established archaeal functions (`EXP4-CODE-REVIEW.txt` item 9).

    Enforces:
    1. Post-hoc grouping by `family_id`/`tax_order` and output-hit filtering CANNOT remove train-time
       leakage from pretrained weights (`cnn_mcm.pt`), `train_db.dmnd`, or STRING graph message passing.
       Unless `training_split_manifest_verified=True` AND a non-empty `disjoint_train_manifest_sha256`
       are provided, `calibration_status` is set to
       `"BLOCKED_unverified_pretrained_train_split_and_graph_leakage_exploratory_only"` and
       `is_calibrated_heldout = False`.
    2. Computes true term-level precision `TP / (TP + FP)` and term-level recall `TP / (TP + FN)`
       (both per-protein macro-average and micro-average across predicted/true GO terms after ancestor
       closure if `go_ancestors_map` is supplied), so an adversarial spam predictor emitting 1,000 false
       terms + 1 true term gets `term_precision_macro = 1 / 1001`, NEVER `1.0`!
    3. Explicitly labels open-world annotation proxy limits (`absence_of_annotation != confirmed_negative`).
    4. Rejects synthetic gap-filling demands (`is_synthetic_gapfill_demand == True`) and untested
       hypothesized roles (`is_hypothesized_untested_role == True`).
    """
    required_gt = {"protein_id", "family_id", "tax_order", "true_go_terms", "is_synthetic_gapfill_demand"}
    missing_gt = sorted(required_gt - set(ground_truth_df.columns))
    if missing_gt:
        raise ValueError(f"ground_truth_df missing required columns: {missing_gt}")

    if ground_truth_df["is_synthetic_gapfill_demand"].astype(bool).any():
        raise ValueError(
            "Synthetic gap-filling demands MUST NOT be used as ground-truth labels in GSPA evaluation!"
        )
    if (
        "is_hypothesized_untested_role" in ground_truth_df.columns
        and ground_truth_df["is_hypothesized_untested_role"].astype(bool).any()
    ):
        raise ValueError(
            "Untested hypothesized roles MUST NOT be used as ground-truth evaluation labels!"
        )

    required_cl = {"protein_id", "arm", "predicted_go_term", "score", "is_leaked_self_or_same_family_hit"}
    missing_cl = sorted(required_cl - set(claims_df.columns))
    if missing_cl:
        raise ValueError(f"claims_df missing required columns: {missing_cl}")

    is_truly_heldout_calibrated = bool(
        training_split_manifest_verified
        and disjoint_train_manifest_sha256
        and str(disjoint_train_manifest_sha256).strip() not in ("", "none", "missing")
    )
    cal_status = (
        f"verified_train_time_disjoint_holdout(sha256={disjoint_train_manifest_sha256})"
        if is_truly_heldout_calibrated
        else "BLOCKED_unverified_pretrained_train_split_and_graph_leakage_exploratory_only"
    )

    n_leaked_removed = int(claims_df["is_leaked_self_or_same_family_hit"].astype(bool).sum())
    clean_claims = claims_df[~claims_df["is_leaked_self_or_same_family_hit"].astype(bool)].copy()

    def _propagate(terms: set[str]) -> set[str]:
        if not go_ancestors_map:
            return set(terms)
        closed = set(terms)
        for t in list(terms):
            closed.update(go_ancestors_map.get(t, set()))
        closed -= {"GO:0003674", "GO:0008150", "GO:0005575"}
        return closed

    gt_lookup: dict[str, set[str]] = {}
    for r in ground_truth_df.to_dict(orient="records"):
        pid = str(r["protein_id"])
        raw_terms = {
            t.strip()
            for t in str(r["true_go_terms"]).split(";")
            if t.strip() and t.strip() not in ("none", "GO:0003674", "GO:0008150", "GO:0005575")
        }
        closed_terms = _propagate(raw_terms)
        if closed_terms:
            gt_lookup[pid] = closed_terms

    groups = sorted(ground_truth_df[holdout_col].astype(str).unique())
    arms = (
        sorted(clean_claims["arm"].astype(str).unique())
        if not clean_claims.empty
        else ["baseline", "context", "shuffled"]
    )

    rows: list[dict[str, Any]] = []
    for grp in groups:
        grp_pids = [
            str(p)
            for p in ground_truth_df.loc[ground_truth_df[holdout_col].astype(str) == grp, "protein_id"]
            if str(p) in gt_lookup
        ]
        if not grp_pids:
            continue
        for arm in arms:
            arm_sub = clean_claims[
                (clean_claims["arm"].astype(str) == arm)
                & (clean_claims["protein_id"].astype(str).isin(set(grp_pids)))
            ]
            top1_hits = 0
            any_hit_proteins = 0
            thresh_pred_proteins = 0
            per_prot_precisions: list[float] = []
            per_prot_recalls: list[float] = []
            total_tp = 0
            total_fp = 0
            total_fn = 0

            for pid in grp_pids:
                true_set = gt_lookup[pid]
                p_rows = arm_sub[arm_sub["protein_id"].astype(str) == pid].sort_values(
                    ["score", "predicted_go_term"], ascending=[False, True]
                )
                if not p_rows.empty:
                    top_term = str(p_rows.iloc[0]["predicted_go_term"])
                    top_closed = _propagate({top_term})
                    if top_term in true_set or (top_closed and top_closed.issubset(true_set)):
                        top1_hits += 1

                passing = p_rows[p_rows["score"].astype(float) >= score_threshold] if not p_rows.empty else p_rows
                raw_pred_set = {
                    str(t).strip()
                    for t in passing["predicted_go_term"]
                    if str(t).strip() not in ("", "none", "GO:0003674", "GO:0008150", "GO:0005575")
                }
                pred_set = _propagate(raw_pred_set)

                tp = len(pred_set & true_set)
                fp = len(pred_set - true_set)
                fn = len(true_set - pred_set)
                total_tp += tp
                total_fp += fp
                total_fn += fn

                if pred_set:
                    thresh_pred_proteins += 1
                    per_prot_precisions.append(tp / float(tp + fp))
                    if tp > 0:
                        any_hit_proteins += 1
                per_prot_recalls.append(tp / float(tp + fn) if (tp + fn) > 0 else np.nan)

            n_eval = len(grp_pids)
            top1_rec = top1_hits / float(n_eval)
            any_hit_rate = any_hit_proteins / float(n_eval)
            macro_prec = float(np.mean(per_prot_precisions)) if per_prot_precisions else np.nan
            macro_rec = float(np.nanmean(per_prot_recalls)) if per_prot_recalls else np.nan
            micro_prec = (total_tp / float(total_tp + total_fp)) if (total_tp + total_fp) > 0 else np.nan
            micro_rec = (total_tp / float(total_tp + total_fn)) if (total_tp + total_fn) > 0 else np.nan

            rows.append(
                {
                    "grouping_scheme": f"grouped_by_{holdout_col}",
                    "group_id": grp,
                    "arm": arm,
                    "is_calibrated_heldout": is_truly_heldout_calibrated,
                    "calibration_status": cal_status,
                    "n_eval_established_proteins": n_eval,
                    "score_threshold": float(score_threshold),
                    "top1_term_recovery": float(top1_rec),
                    "term_precision_macro_proxy": macro_prec,
                    "term_recall_macro": macro_rec,
                    "term_precision_micro_proxy": micro_prec,
                    "term_recall_micro": micro_rec,
                    "protein_any_true_term_hit_rate": float(any_hit_rate),
                    "n_proteins_above_threshold": int(thresh_pred_proteins),
                    "total_tp_terms": int(total_tp),
                    "total_fp_terms_open_world_proxy": int(total_fp),
                    "total_fn_terms": int(total_fn),
                    "n_output_leaked_hits_masked": n_leaked_removed,
                    "go_ancestor_propagation_applied": bool(go_ancestors_map is not None),
                    "synthetic_gapfill_labels_excluded": True,
                    "annotation_based_precision_proxy_note": (
                        "open-world GO annotation proxy: unannotated terms counted as FP in term_precision_*_proxy, "
                        "though absence of annotation != confirmed experimental negative; "
                        "pretrained model evaluation without disjoint train-time fold manifest is exploratory and NOT calibrated"
                    ),
                }
            )

    return pd.DataFrame(rows)


def evaluate_same_role_evidence_and_convergence(
    row: Mapping[str, Any],
    *,
    fdr_threshold: float = 0.10,
    use_joint_matrix_q: bool = True,
) -> dict[str, Any]:
    """Evaluate role-specific convergence across sequence, structure, architecture, and context for the SAME biological role.

    Enforces `EXP4-POSTRUN-SCIENTIFIC-REVIEW.txt` Items 4 & 5 and `EXP4-SCHEMA-REGRESSION-REVIEW.txt`:
    - Requires sequence/structure and context to support the SAME biological role (preventing negative-control
      folds or unrelated structures like CcmE from claiming peptidase/pilin convergence).
    - Differentiates transparent evidence routes without dropping high-confidence sequence+context leads:
        1. `convergent_sequence_structure_and_same_role_context`
        2. `qualified_structure_and_same_role_context_with_architecture`
        3. `sequence_and_same_role_context_with_architecture`
        4. `domain_profile_and_same_role_context_with_architecture`
        5. `exploratory_fold_or_motif_candidate_no_operon_convergence`
        6. `uncertain_insufficient_evidence`
    - Explicitly checks actual pipeline status strings:
        * Accepts `exact_sequence_hash_match` and `exact_match` for `model_polymer_hash_status`
        * Validates structured topology (`topology_status.startswith("verified_tmbed")`); missing stays unknown
        * Distinguishes domain profile match alone from direct positive reference sequence match
        * Audits role-specific hits from `detail_{r}` so positive second-best hits are not masked by confusable controls
    - Explicitly labels every outcome as a computational role hypothesis, NEVER a biochemical assignment.
    """
    model_source = str(row.get("model_source", "none"))
    model_hash_status = str(row.get("model_polymer_hash_status", "missing_model"))
    has_exact_model = (
        bool(row.get("model_polymer_exact_match", False))
        or model_hash_status in ("exact_match", "exact_sequence_hash_match")
    )
    plddt_val = float(row.get("model_plddt_0to1", np.nan))
    model_plddt_qualified = bool(has_exact_model and np.isfinite(plddt_val) and plddt_val >= 0.70)

    # Topology and export status
    topo_status = str(row.get("topology_status", ""))
    has_verified_topology = topo_status.startswith("verified_tmbed")
    loc_class = str(row.get("surface_localization_class", ""))
    is_exported = bool(row.get("is_exported_any", False))
    tm_status = str(row.get("tm_anchor_position_status", ""))
    has_prepilin = (
        bool(row.get("has_prepilin_classIII_motif_pred", False))
        or tm_status.startswith("verified_n_terminal")
        or loc_class == "prepilin_or_archaellin_filament_subunit_candidate"
    )
    dom_text = str(row.get("rep_own_domains", "")) + " " + str(row.get("rep_own_c6_hits", ""))

    # Guards and exclusions
    pdb_guard = str(row.get("pdb_guard_status", ""))
    pdb_guard_ok = not (pdb_guard.startswith("rejected_") or pdb_guard.startswith("blocked_"))
    sbp_status = str(row.get("sbp_substrate_specificity_status", ""))

    # Overall best references
    best_seq_ref = str(row.get("best_curated_seq_ref_id", "none"))
    best_seq_spec = PANEL_SPECS_BY_ID.get(best_seq_ref, {})
    best_seq_ev = float(row.get("best_curated_seq_evalue", np.nan))
    best_seq_qc = float(row.get("best_curated_seq_qcov", 0.0))
    best_seq_tc = float(row.get("best_curated_seq_tcov", 0.0))

    best_struct_ref = str(row.get("best_curated_struct_ref_id", "none"))
    best_struct_spec = PANEL_SPECS_BY_ID.get(best_struct_ref, {})
    best_st_qtm = float(row.get("best_curated_struct_qtmscore", 0.0))
    best_st_qc = float(row.get("best_curated_struct_qcov", 0.0))
    best_st_tc = float(row.get("best_curated_struct_tcov", 0.0))

    role_candidates: list[dict[str, Any]] = []

    for r_name in ROLE_CLASSES:
        # 1. Operon context for this specific role
        g_name = "cell_coat_or_envelope" if r_name == "cell_coat_or_envelope" else (
            "host_defense_antiviral" if r_name == "virus_interaction_or_entry" else r_name
        )
        obs_adj = int(row.get(f"null_obs_adj_{g_name}", 0))
        n_gb = int(row.get(f"null_n_gblocks_{g_name}", 0))
        joint_q = float(row.get(f"joint_matrix_q_bh_{g_name}", 1.0))
        glob_q = float(row.get(f"global_q_bh_{g_name}", 1.0))
        q_val = joint_q if use_joint_matrix_q else glob_q
        has_ctx = bool(obs_adj >= 2 and n_gb >= 2 and q_val <= fdr_threshold)

        # 2. Check role-specific hits from detail_{r_name} or overall best
        detail_str = str(row.get(f"detail_{r_name}", ""))
        detail_seq_ref = "none"
        detail_struct_ref = "none"
        if detail_str and detail_str != "none":
            for part in detail_str.split(";"):
                if part.startswith("seq_ref="):
                    detail_seq_ref = part.split("=", 1)[1]
                elif part.startswith("struct_ref="):
                    detail_struct_ref = part.split("=", 1)[1]

        # Check positive sequence hit for r_name
        has_valid_seq = False
        seq_ref_match = "none"
        if best_seq_ref != "none" and r_name in best_seq_spec.get("experimentally_supported_roles", "").split(";"):
            if best_seq_spec.get("control_type", "").startswith("positive"):
                if np.isfinite(best_seq_ev) and best_seq_ev <= 1e-5 and best_seq_qc >= 0.50 and best_seq_tc >= 0.50:
                    has_valid_seq = True
                    seq_ref_match = best_seq_ref
        if not has_valid_seq and detail_seq_ref not in ("none", "") and detail_seq_ref in PANEL_SPECS_BY_ID:
            ds_spec = PANEL_SPECS_BY_ID[detail_seq_ref]
            if r_name in ds_spec.get("experimentally_supported_roles", "").split(";") and ds_spec.get("control_type", "").startswith("positive"):
                has_valid_seq = True
                seq_ref_match = detail_seq_ref

        # Check positive structure hit for r_name
        has_valid_struct = False
        struct_ref_match = "none"
        if best_struct_ref != "none" and r_name in best_struct_spec.get("experimentally_supported_roles", "").split(";"):
            if best_struct_spec.get("control_type", "").startswith("positive"):
                if model_plddt_qualified and best_st_qtm >= 0.50 and best_st_qc >= 0.50 and best_st_tc >= 0.50 and pdb_guard_ok:
                    has_valid_struct = True
                    struct_ref_match = best_struct_ref
        if not has_valid_struct and detail_struct_ref not in ("none", "") and detail_struct_ref in PANEL_SPECS_BY_ID:
            dst_spec = PANEL_SPECS_BY_ID[detail_struct_ref]
            if dst_spec.get("control_type", "").startswith("positive") and r_name in dst_spec.get("experimentally_supported_roles", "").split(";"):
                if model_plddt_qualified and pdb_guard_ok:
                    has_valid_struct = True
                    struct_ref_match = detail_struct_ref

        # Domain profile match (distinct from direct sequence reference)
        pat_key = "host_defense_antiviral" if r_name == "virus_interaction_or_entry" else r_name
        has_dom_match = bool(BOUNDED_GUILD_PATTERNS.get(pat_key, re.compile(r"$^")).search(dom_text))

        # Negative control and confusable exclusions
        if r_name == "nutrient_acquisition_peptidase":
            if "Propeptide" in seq_ref_match or "9UES_A" in seq_ref_match:
                has_valid_seq = False
            if "Propeptide" in struct_ref_match or "9UES_A" in struct_ref_match or pdb_guard.startswith("rejected_9UES_chainA"):
                has_valid_struct = False
        elif r_name == "substrate_binding_or_uptake":
            if sbp_status == "generic_sbp_fold_substrate_unresolved_peptide_vs_sugar":
                if "Sugar" in seq_ref_match or "Ctrl" in seq_ref_match:
                    has_valid_seq = False
                if "Sugar" in struct_ref_match or "Ctrl" in struct_ref_match:
                    has_valid_struct = False
        elif r_name == "surface_assembly_pilin_archaellin":
            # Assembly ATPases cannot establish pilin subunit
            if any(mach in seq_ref_match for mach in ("ArlI", "ArlJ", "FlaI", "PilB", "PilT")):
                has_valid_seq = False
            if any(mach in struct_ref_match for mach in ("ArlI", "ArlJ", "FlaI", "PilB", "PilT")):
                has_valid_struct = False
        elif r_name == "virus_interaction_or_entry":
            viral_source_hits = str(row.get("viral_origin_audit_status", ""))
            if "viral_virion_protein_similarity" in viral_source_hits:
                has_valid_seq = False
                has_valid_struct = False
                has_dom_match = False

        # Architecture compatibility: missing TMbed stays unknown, NOT default eligible
        if not has_verified_topology:
            arch_compatible = False
        elif r_name == "surface_assembly_pilin_archaellin":
            arch_compatible = has_prepilin
        elif r_name == "cell_coat_or_envelope":
            arch_compatible = loc_class in (
                "c_terminal_anchored_surface_coat_candidate",
                "lipoprotein_membrane_anchored_SPII",
                "freely_secreted_extracellular_candidate",
                "prepilin_or_archaellin_filament_subunit_candidate",
            )
        elif r_name in ("nutrient_acquisition_peptidase", "substrate_binding_or_uptake"):
            arch_compatible = is_exported
        else:
            arch_compatible = is_exported or loc_class != "intracellular_soluble"

        # Evidence routes
        if has_ctx and arch_compatible:
            if has_valid_seq and has_valid_struct:
                role_candidates.append({
                    "role": r_name,
                    "route": "convergent_sequence_structure_and_same_role_context",
                    "score": 0.95,
                    "has_valid_seq": True,
                    "has_valid_struct": True,
                    "has_dom_match": has_dom_match,
                    "has_ctx": True,
                    "seq_ref": seq_ref_match,
                    "struct_ref": struct_ref_match,
                })
            elif has_valid_struct:
                role_candidates.append({
                    "role": r_name,
                    "route": "qualified_structure_and_same_role_context_with_architecture",
                    "score": 0.80,
                    "has_valid_seq": False,
                    "has_valid_struct": True,
                    "has_dom_match": has_dom_match,
                    "has_ctx": True,
                    "seq_ref": "none",
                    "struct_ref": struct_ref_match,
                })
            elif has_valid_seq:
                role_candidates.append({
                    "role": r_name,
                    "route": "sequence_and_same_role_context_with_architecture",
                    "score": 0.75,
                    "has_valid_seq": True,
                    "has_valid_struct": False,
                    "has_dom_match": has_dom_match,
                    "has_ctx": True,
                    "seq_ref": seq_ref_match,
                    "struct_ref": "none",
                })
            elif has_dom_match:
                role_candidates.append({
                    "role": r_name,
                    "route": "domain_profile_and_same_role_context_with_architecture",
                    "score": 0.60,
                    "has_valid_seq": False,
                    "has_valid_struct": False,
                    "has_dom_match": True,
                    "has_ctx": True,
                    "seq_ref": "none",
                    "struct_ref": "none",
                })

    if role_candidates:
        role_candidates.sort(key=lambda x: -x["score"])
        best = role_candidates[0]
        supported_route = best["route"]
        best_role = best["role"]
        has_pos_seq = best["has_valid_seq"]
        has_pos_struct = best["has_valid_struct"]
        has_dom = best["has_dom_match"]
        matched_seq_ref = best["seq_ref"]
        matched_struct_ref = best["struct_ref"]
    elif float(row.get("max_role_evidence_score", 0.0)) >= 0.40:
        supported_route = "exploratory_fold_or_motif_candidate_no_operon_convergence"
        best_role = str(row.get("top_competing_roles", "none")).split(";")[0]
        has_pos_seq = False
        has_pos_struct = False
        has_dom = False
        matched_seq_ref = best_seq_ref
        matched_struct_ref = best_struct_ref
    else:
        supported_route = "uncertain_insufficient_evidence"
        best_role = "none"
        has_pos_seq = False
        has_pos_struct = False
        has_dom = False
        matched_seq_ref = "none"
        matched_struct_ref = "none"

    return {
        "supported_evidence_route": supported_route,
        "candidate_computational_role_hypothesis": best_role,
        "is_computational_role_hypothesis": True,
        "is_biochemical_assignment": False,
        "has_direct_positive_seq_ref_match": has_pos_seq,
        "has_direct_positive_struct_ref_match": has_pos_struct,
        "has_domain_profile_match": has_dom,
        "matched_positive_seq_ref_id": matched_seq_ref,
        "matched_positive_struct_ref_id": matched_struct_ref,
        "model_availability_status": (
            f"model_source={model_source};hash={model_hash_status};plddt={plddt_val:.2f}"
            if np.isfinite(plddt_val)
            else f"model_source={model_source};hash={model_hash_status};plddt=nan"
        ),
        "scientific_disposition_note": (
            "Computational candidate role hypothesis; domain annotations and structural matches "
            "are predictors, not experimental biochemical assignments."
        ),
    }
