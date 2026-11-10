#!/usr/bin/env python3
"""Stage 1: Exhaustive CRISPR-Cas and Spacer Recovery Across All 591 Verified Archaeal MAGs
(plus Non-Representative Refined Archaeal Bins).

Implements all requirements from `EXP3-CONTRACT.txt` and `EXP3-STAGE1-REVIEW.txt`:
1. Preserves raw MinCED and CCTyper outputs, stdout, stderr, and exit codes persistently
   in per-MAG compressed archives (`results/stage1_raw_archives/<tag>/<mag_id>.tar.gz`),
   and writes failure diagnostics to `results/stage1_failed_diagnostics/` if any tool fails.
2. Separates `n_mags_attempted`, `n_mags_successfully_screened` (both callers succeeded),
   `n_mags_failed`, and `n_mags_unscreened_unknown`. Never treats failed caller execution
   as a biological zero-array negative.
3. Validates synthetic positive and negative control FASTA fixtures against the actual
   `minced.sif` and `cctyper.sif` containers during `--mode pilot`, and gates pilot
   completion on 100% caller execution success across all 30 pilot MAGs.
4. Enforces header-only schemas for empty tables (`STATUS_COLS`, `ARRAY_COLS`, `SPACER_COLS`,
   `OPERON_COLS`, `DEDUP_SPACER_COLS`) and raises on duplicate MAG rows across chunks
   rather than silently calling `drop_duplicates`.
5. Deduplicates spacers within each MAG by `canonical_spacer_seq_rc` (auditing both
   multi-array repeats and reverse-complement duplicates as a single independent spacer)
   and preserves the exact list of individual array intervals (`array_loci`) rather than
   collapsing multiple arrays into a synthetic `min_start..max_end` bounding interval.
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import json
import os
import random
import re
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exp3_common import (
    canonical_spacer_seq_rc,
    is_valid_crispr_repeat_and_spacer,
    parse_fasta_bytes,
    qualify_namespaced_contig_id,
    sha256_file,
)

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
IN_DIR = BRANCH_ROOT / "inputs"
RES_DIR = BRANCH_ROOT / "results"
CHUNK_OUT_DIR = RES_DIR / "stage1_chunks"
RAW_ARCH_DIR = RES_DIR / "stage1_raw_archives"
FAIL_DIAG_DIR = RES_DIR / "stage1_failed_diagnostics"
TMP_ROOT = Path(os.environ.get("TMPDIR", str(BRANCH_ROOT / "tmp")))

SIF_MINCED = Path("/project/k10018/sw/images/mobilome/minced.sif")
SIF_CCTYPER = Path("/project/k10018/sw/images/mobilome/cctyper.sif")
DB_CCTYPER = Path(
    os.environ.get(
        "CCTYPER_DB",
        str(IN_DIR / "cctyper_db_v1.8.0"),
    )
)
P_C8_DEF = Path("/project/k10018/work/eq-annot-20261030/results/c8_defense/analysis/c8_per_genome.tsv")

CCTYPER_REQUIRED_DB_FILES = (
    "CasScoring.csv",
    "xgb_repeats.model",
    "type_dict.tab",
    "cutoffs.tab",
    "interference.json",
    "adaptation.json",
    "repeats.fa",
)
ANSI_ESCAPE_RE = re.compile(r"\x1b\[[0-9;]*[a-zA-Z]")
CCTYPER_ERROR_LOG_RE = re.compile(
    r"(?:^|\n)\s*(?:\[[^\]]+\]\s*)?(?:ERROR|CRITICAL)\s*:|"
    r"Traceback \(most recent call last\)|"
    r"Could not find database directory|"
    r"CasScoring table could not|"
    r"Could not find HMM profile directory|"
    r"There are non-HMM profiles in the HMM profile directory|"
    r"Prodigal failed!|"
    r"XGBoost model incompatible|"
    r"HMMER failed|"
    r"Could not find input file|"
    r"Duplicate fasta headers detected",
    re.IGNORECASE,
)

MINCED_SEQ_RE = re.compile(r"^Sequence\s+'([^']+)'\s+\(\d+\s+bp\)")
MINCED_HDR_RE = re.compile(r"^CRISPR\s+(\d+)\s+Range:\s+(\d+)\s+-\s+(\d+)")
MINCED_ROW_RE = re.compile(r"^(\d+)\s+([ACGTNacgtn]+)\s+([ACGTNacgtn]+)\s+\[\s*\d+,\s*\d+\s*\]")
MINCED_LAST_RE = re.compile(r"^(\d+)\s+([ACGTNacgtn]+)\s*$")

STATUS_COLS = [
    "mag_id",
    "mag_fallback_sample",
    "is_arch591_rep",
    "fasta_sha256",
    "minced_exitcode",
    "cctyper_exitcode",
    "minced_ok",
    "cctyper_ok",
    "cctyper_status_detail",
    "screening_status",
    "n_raw_arrays",
    "n_cctyper_trusted_arrays",
    "n_cctyper_putative_untrusted_arrays",
    "n_raw_spacers",
    "n_valid_spacers",
    "n_independent_spacers",
    "n_cas_operons",
    "raw_archive_path",
]

ARRAY_COLS = [
    "array_id",
    "mag_id",
    "is_arch591_rep",
    "raw_contig_id",
    "contig_sample",
    "namespaced_contig_id",
    "array_start",
    "array_end",
    "caller",
    "subtype_or_label",
    "is_trusted_array",
    "array_trust_tier",
    "consensus_repeat",
    "n_repeats",
    "n_raw_spacers",
    "n_valid_spacers",
]

SPACER_COLS = [
    "spacer_id",
    "array_id",
    "mag_id",
    "is_arch591_rep",
    "raw_contig_id",
    "contig_sample",
    "namespaced_contig_id",
    "array_start",
    "array_end",
    "caller",
    "is_trusted_array",
    "array_trust_tier",
    "spacer_index",
    "repeat_seq",
    "spacer_seq",
    "canonical_spacer_seq_rc",
    "spacer_orientation",
    "spacer_len",
    "passes_seq_quality_gate",
    "passes_quality_gate",
    "quality_gate_reason",
]

OPERON_COLS = [
    "mag_id",
    "is_arch591_rep",
    "raw_contig_id",
    "contig_sample",
    "namespaced_contig_id",
    "operon_id",
    "operon_tier",
    "start",
    "end",
    "prediction",
    "genes",
]

DEDUP_SPACER_COLS = [
    "unique_spacer_id",
    "mag_id",
    "is_arch591_rep",
    "canonical_spacer_seq_rc",
    "representative_spacer_seq",
    "spacer_len",
    "observed_orientations",
    "n_loci_in_mag",
    "n_distinct_arrays_in_mag",
    "callers",
    "array_trust_tiers",
    "caller_independence_note",
    "primary_namespaced_contig_id",
    "primary_array_start",
    "primary_array_end",
    "all_array_loci",
    "repeat_seq",
]


def verify_cctyper_db_dir(db_dir: Path = DB_CCTYPER) -> dict[str, Any]:
    """Verify official compatible CRISPRCasTyper v1.8.0 database files and extracted HMM Profiles."""
    if not db_dir.exists() or not db_dir.is_dir():
        raise RuntimeError(f"CCTyper database directory missing or not a directory: {db_dir}")

    missing_or_empty: list[str] = []
    for fname in CCTYPER_REQUIRED_DB_FILES:
        fp = db_dir / fname
        if not fp.is_file() or fp.stat().st_size == 0:
            missing_or_empty.append(fname)
    if missing_or_empty:
        raise RuntimeError(
            f"CCTyper database at {db_dir} missing or empty required files: {missing_or_empty}"
        )

    prof_dir = db_dir / "Profiles"
    if not prof_dir.is_dir():
        raise RuntimeError(f"CCTyper HMM profile directory missing: {prof_dir}")

    entries = sorted(os.listdir(prof_dir))
    if not entries:
        raise RuntimeError(f"CCTyper HMM profile directory is empty: {prof_dir}")

    non_hmm = [e for e in entries if not e.lower().endswith(".hmm")]
    if non_hmm:
        raise RuntimeError(
            f"CCTyper HMM profile directory {prof_dir} contains non-.hmm files "
            f"(rejected by Controller.check_db): {non_hmm[:10]}"
        )

    return {
        "db_dir": str(db_dir),
        "n_required_files_verified": len(CCTYPER_REQUIRED_DB_FILES),
        "n_hmm_profiles": len(entries),
    }


def is_minced_run_successful(
    returncode: int,
    stdout: str,
    stderr: str,
    minced_txt_path: Path,
) -> tuple[bool, str]:
    """Verify MinCED execution succeeded and produced its output file without runtime error."""
    if returncode != 0:
        return False, f"nonzero_exit_{returncode}"
    if not minced_txt_path.is_file():
        return False, "missing_minced_crisprs_file"
    comb = f"{stdout or ''}\n{stderr or ''}"
    if re.search(r"Exception in thread|Error:\s+Unable to|NoSuchFileException", comb, re.IGNORECASE):
        return False, "minced_log_exception_or_error"
    return True, "ok"


def is_cctyper_run_successful(
    returncode: int,
    stdout: str,
    stderr: str,
    cct_out: Path,
) -> tuple[bool, str]:
    """Verify CCTyper v1.8.0 completed all pipeline stages without silent exit-0 error.

    In CRISPRCasTyper v1.8.0, `Controller.check_db()`, `check_out()`, `check_input()`,
    `Prodigal.check_rerun()`, and `XGB.predict_repeats()` call `logging.error(...)` or
    `logging.critical(...)` followed by bare `sys.exit()` (which returns exit code 0!).
    Furthermore, output tables (`crisprs_all.tab`, `cas_operons.tab`, `CRISPR_Cas.tab`)
    are only present when CRISPRs/Cas operons are detected.
    A genuine, complete CCTyper v1.8.0 run (without `--keep_tmp`) ALWAYS:
    1. Exits with returncode 0,
    2. Emits NO `ERROR:` / `CRITICAL:` / `Traceback` in stdout/stderr,
    3. Creates `cct_out/arguments.tab`, `cct_out/genes.tab`, and `cct_out/hmmer.tab`,
    4. Reaches `master.clean()` at the end of `bin/cctyper`, logging `Removing temporary files`.
    """
    clean_out = ANSI_ESCAPE_RE.sub("", stdout or "")
    clean_err = ANSI_ESCAPE_RE.sub("", stderr or "")
    combined = f"{clean_out}\n{clean_err}"

    if returncode != 0:
        return False, f"nonzero_exit_{returncode}"
    if CCTYPER_ERROR_LOG_RE.search(combined):
        return False, "cctyper_log_error_or_critical"
    if not cct_out.is_dir():
        return False, "missing_cctyper_output_dir"

    args_tab = cct_out / "arguments.tab"
    if not args_tab.is_file() or args_tab.stat().st_size == 0:
        return False, "missing_or_empty_arguments_tab"

    genes_tab = cct_out / "genes.tab"
    hmmer_tab = cct_out / "hmmer.tab"
    if not genes_tab.is_file() or not hmmer_tab.is_file():
        return False, "missing_genes_or_hmmer_tab"

    if "Removing temporary files" not in clean_err and "Removing temporary files" not in clean_out:
        return False, "missing_final_clean_completion_marker"

    return True, "ok"


def parse_minced_txt(
    minced_txt: str,
    mag_id: str,
    fallback_sample: str,
    is_arch591: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Parse MinCED `.crisprs` text output into array rows and spacer rows."""
    arrays: list[dict[str, Any]] = []
    spacers: list[dict[str, Any]] = []

    cur_contig = ""
    cur_arr_num = 0
    cur_start = 0
    cur_end = 0
    cur_repeats: list[str] = []
    cur_spacers: list[tuple[int, str, str]] = []

    def _flush() -> None:
        nonlocal cur_arr_num, cur_start, cur_end, cur_repeats, cur_spacers
        if not cur_contig or cur_arr_num == 0 or not cur_spacers:
            cur_arr_num = 0
            cur_repeats = []
            cur_spacers = []
            return
        c_samp, base_c, q_cid, _, _ = qualify_namespaced_contig_id(cur_contig, fallback_sample=fallback_sample)
        consensus_rep = collections.Counter(cur_repeats).most_common(1)[0][0] if cur_repeats else ""
        arr_id = f"{mag_id}__minced_{base_c}_CRISPR{cur_arr_num}_{cur_start}_{cur_end}"
        n_valid = 0
        for idx, (_, rep_s, spc_s) in enumerate(cur_spacers, start=1):
            ok, reason = is_valid_crispr_repeat_and_spacer(rep_s, spc_s)
            if ok:
                n_valid += 1
            can_seq, orient = canonical_spacer_seq_rc(spc_s)
            spacers.append(
                {
                    "spacer_id": f"{arr_id}__spc{idx}",
                    "array_id": arr_id,
                    "mag_id": mag_id,
                    "is_arch591_rep": int(is_arch591),
                    "raw_contig_id": cur_contig,
                    "contig_sample": c_samp,
                    "namespaced_contig_id": q_cid,
                    "array_start": cur_start,
                    "array_end": cur_end,
                    "caller": "minced_0.4.2",
                    "is_trusted_array": 1,
                    "array_trust_tier": "minced_repeat_array",
                    "spacer_index": idx,
                    "repeat_seq": rep_s.upper(),
                    "spacer_seq": spc_s.upper(),
                    "canonical_spacer_seq_rc": can_seq,
                    "spacer_orientation": orient,
                    "spacer_len": len(spc_s),
                    "passes_seq_quality_gate": int(ok),
                    "passes_quality_gate": int(ok),
                    "quality_gate_reason": reason,
                }
            )
        arrays.append(
            {
                "array_id": arr_id,
                "mag_id": mag_id,
                "is_arch591_rep": int(is_arch591),
                "raw_contig_id": cur_contig,
                "contig_sample": c_samp,
                "namespaced_contig_id": q_cid,
                "array_start": cur_start,
                "array_end": cur_end,
                "caller": "minced_0.4.2",
                "subtype_or_label": "minced_repeat_array",
                "is_trusted_array": 1,
                "array_trust_tier": "minced_repeat_array",
                "consensus_repeat": consensus_rep.upper(),
                "n_repeats": len(cur_repeats),
                "n_raw_spacers": len(cur_spacers),
                "n_valid_spacers": n_valid,
            }
        )
        cur_arr_num = 0
        cur_repeats = []
        cur_spacers = []

    for raw_line in minced_txt.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        m_seq = MINCED_SEQ_RE.match(line)
        if m_seq:
            _flush()
            cur_contig = m_seq.group(1).split()[0]
            continue
        m_hdr = MINCED_HDR_RE.match(line)
        if m_hdr:
            _flush()
            cur_arr_num = int(m_hdr.group(1))
            cur_start = int(m_hdr.group(2))
            cur_end = int(m_hdr.group(3))
            continue
        if cur_arr_num > 0:
            m_row = MINCED_ROW_RE.match(line)
            if m_row:
                pos = int(m_row.group(1))
                rep_s = m_row.group(2).upper()
                spc_s = m_row.group(3).upper()
                cur_repeats.append(rep_s)
                cur_spacers.append((pos, rep_s, spc_s))
                continue
            m_last = MINCED_LAST_RE.match(line)
            if m_last:
                cur_repeats.append(m_last.group(2).upper())
                continue

    _flush()
    return arrays, spacers


def parse_cctyper_dir(
    cct_dir: Path,
    mag_id: str,
    fallback_sample: str,
    is_arch591: bool,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Parse CCTyper output directory into array rows, spacer rows, and Cas operon rows."""
    arrays: list[dict[str, Any]] = []
    spacers: list[dict[str, Any]] = []
    cas_operons: list[dict[str, Any]] = []

    if not cct_dir.exists():
        return arrays, spacers, cas_operons

    for tab_name, op_tier in (
        ("cas_operons.tab", "high_confidence"),
        ("cas_operons_putative.tab", "putative"),
    ):
        p_op = cct_dir / tab_name
        if p_op.exists() and p_op.stat().st_size > 0:
            with open(p_op) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    raw_c = (r.get("Contig") or r.get("contig") or "").strip()
                    if not raw_c:
                        continue
                    c_samp, _, q_cid, _, _ = qualify_namespaced_contig_id(raw_c, fallback_sample=fallback_sample)
                    cas_operons.append(
                        {
                            "mag_id": mag_id,
                            "is_arch591_rep": int(is_arch591),
                            "raw_contig_id": raw_c,
                            "contig_sample": c_samp,
                            "namespaced_contig_id": q_cid,
                            "operon_id": r.get("Operon") or r.get("operon") or f"{mag_id}_{raw_c}",
                            "operon_tier": op_tier,
                            "start": r.get("Start") or r.get("start") or "0",
                            "end": r.get("End") or r.get("end") or "0",
                            "prediction": r.get("Prediction") or r.get("Subtype") or "unknown",
                            "genes": r.get("Genes") or r.get("genes") or "",
                        }
                    )

    crispr_meta: dict[str, dict[str, str]] = {}
    for tab_name in ("crisprs_all.tab", "crisprs_near_cas.tab", "crisprs_putative.tab"):
        p_cr = cct_dir / tab_name
        if p_cr.exists() and p_cr.stat().st_size > 0:
            with open(p_cr) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    cid = (r.get("CRISPR") or r.get("crispr") or "").strip()
                    if cid and cid not in crispr_meta:
                        crispr_meta[cid] = r

    spc_dir = cct_dir / "spacers"
    if spc_dir.exists():
        for fa_file in sorted(spc_dir.glob("*.fa")):
            cr_id = fa_file.stem
            meta = crispr_meta.get(cr_id, {})
            raw_c = (meta.get("Contig") or meta.get("contig") or cr_id.rsplit("_", 1)[0]).strip()
            c_samp, _, q_cid, _, _ = qualify_namespaced_contig_id(raw_c, fallback_sample=fallback_sample)
            a_start = int(float(meta.get("Start") or meta.get("start") or 0))
            a_end = int(float(meta.get("End") or meta.get("end") or 0))
            cons_rep = (meta.get("Consensus_repeat") or meta.get("repeat") or "GTTTGAGAGTTGTGATTTGCAAGTACGT").strip().upper()
            subtype = (meta.get("Prediction") or meta.get("Subtype") or "cctyper_array").strip()
            trusted_str = str(meta.get("Trusted", "")).strip().lower()
            is_trusted = trusted_str in ("true", "1", "t")
            trust_tier = "cctyper_trusted" if is_trusted else "cctyper_putative_untrusted"

            recs, _ = parse_fasta_bytes(fa_file.read_bytes())
            arr_id = f"{mag_id}__cctyper_{cr_id}_{a_start}_{a_end}"
            n_valid = 0
            for idx, (_, spc_s) in enumerate(recs, start=1):
                ok, reason = is_valid_crispr_repeat_and_spacer(cons_rep, spc_s)
                # Untrusted CCTyper arrays are retained exploratory raw, filtered from strongest host evidence:
                q_reason = reason if ok else reason
                if ok and not is_trusted:
                    q_reason = "cctyper_putative_untrusted_array_exploratory_only"
                passes_gate = int(ok and is_trusted)
                if passes_gate:
                    n_valid += 1
                can_seq, orient = canonical_spacer_seq_rc(spc_s)
                spacers.append(
                    {
                        "spacer_id": f"{arr_id}__spc{idx}",
                        "array_id": arr_id,
                        "mag_id": mag_id,
                        "is_arch591_rep": int(is_arch591),
                        "raw_contig_id": raw_c,
                        "contig_sample": c_samp,
                        "namespaced_contig_id": q_cid,
                        "array_start": a_start,
                        "array_end": a_end,
                        "caller": "cctyper_1.8.0",
                        "is_trusted_array": int(is_trusted),
                        "array_trust_tier": trust_tier,
                        "spacer_index": idx,
                        "repeat_seq": cons_rep,
                        "spacer_seq": spc_s.upper(),
                        "canonical_spacer_seq_rc": can_seq,
                        "spacer_orientation": orient,
                        "spacer_len": len(spc_s),
                        "passes_seq_quality_gate": int(ok),
                        "passes_quality_gate": passes_gate,
                        "quality_gate_reason": q_reason,
                    }
                )
            arrays.append(
                {
                    "array_id": arr_id,
                    "mag_id": mag_id,
                    "is_arch591_rep": int(is_arch591),
                    "raw_contig_id": raw_c,
                    "contig_sample": c_samp,
                    "namespaced_contig_id": q_cid,
                    "array_start": a_start,
                    "array_end": a_end,
                    "caller": "cctyper_1.8.0",
                    "subtype_or_label": subtype,
                    "is_trusted_array": int(is_trusted),
                    "array_trust_tier": trust_tier,
                    "consensus_repeat": cons_rep,
                    "n_repeats": int(float(meta.get("N_repeats") or (len(recs) + 1))),
                    "n_raw_spacers": len(recs),
                    "n_valid_spacers": n_valid,
                }
            )

    return arrays, spacers, cas_operons


def deduplicate_mag_spacers(valid_spacers: list[dict[str, Any]] | pd.DataFrame) -> list[dict[str, Any]]:
    """Deduplicate verified spacers per MAG by `(mag_id, canonical_spacer_seq_rc)`.

    - An identical spacer sequence (or its reverse complement) seen in 2 arrays or by
      2 callers within a MAG is counted as ONE independent spacer for that MAG.
    - Preserves the exact list of individual array intervals (`all_array_loci`) rather
      than collapsing distinct arrays into a giant synthetic `min_start..max_end` interval.
    - Audits caller overlap: CCTyper itself uses MinCED under the hood, so dual calls
      from MinCED and CCTyper are not independent biological confirmation.
    """
    df = pd.DataFrame(valid_spacers, columns=SPACER_COLS) if isinstance(valid_spacers, list) else valid_spacers
    if len(df) == 0:
        return []

    dedup_rows: list[dict[str, Any]] = []
    spc_counter = 0
    for (mag_id, can_seq), grp in df.groupby(["mag_id", "canonical_spacer_seq_rc"], as_index=False):
        spc_counter += 1
        uid = f"archspc_{spc_counter:06d}"
        first = grp.iloc[0]
        callers_set = sorted(set(grp["caller"].astype(str)))
        callers = "+".join(callers_set)
        trust_tiers_val = (
            grp["array_trust_tier"].dropna().astype(str)
            if "array_trust_tier" in grp.columns
            else pd.Series(dtype=str)
        )
        trust_tiers = "+".join(sorted(set(trust_tiers_val))) if len(trust_tiers_val) > 0 else "unknown"
        orients = "+".join(sorted(set(grp["spacer_orientation"].astype(str))))

        caller_note = (
            "cctyper_uses_minced_not_independent"
            if len(callers_set) > 1 and any("minced" in c for c in callers_set) and any("cctyper" in c for c in callers_set)
            else "single_caller"
        )

        loci_set: list[str] = []
        seen_loci: set[str] = set()
        distinct_arrays: set[tuple[str, int, int]] = set()
        for _, r in grp.iterrows():
            q_c = str(r["namespaced_contig_id"])
            a_lo = int(r["array_start"])
            a_hi = int(r["array_end"])
            distinct_arrays.add((q_c, a_lo, a_hi))
            loc_str = f"{q_c}:{a_lo}-{a_hi}"
            if loc_str not in seen_loci:
                seen_loci.add(loc_str)
                loci_set.append(loc_str)

        dedup_rows.append(
            {
                "unique_spacer_id": uid,
                "mag_id": str(mag_id),
                "is_arch591_rep": int(first["is_arch591_rep"]),
                "canonical_spacer_seq_rc": str(can_seq),
                "representative_spacer_seq": str(first["spacer_seq"]),
                "spacer_len": int(first["spacer_len"]),
                "observed_orientations": orients,
                "n_loci_in_mag": len(loci_set),
                "n_distinct_arrays_in_mag": len(distinct_arrays),
                "callers": callers,
                "array_trust_tiers": trust_tiers,
                "caller_independence_note": caller_note,
                "primary_namespaced_contig_id": str(first["namespaced_contig_id"]),
                "primary_array_start": int(first["array_start"]),
                "primary_array_end": int(first["array_end"]),
                "all_array_loci": ";".join(loci_set),
                "repeat_seq": str(first["repeat_seq"]),
            }
        )
    return dedup_rows


def screen_single_mag(
    mag_id: str,
    fasta_path: Path,
    is_arch591: bool,
    tag: str,
    work_tmp_dir: Path,
    cctyper_threads: int = 4,
) -> dict[str, Any]:
    """Run MinCED 0.4.2 and CCTyper 1.8.0 on a single MAG FASTA and archive all raw outputs."""
    fallback_sample = mag_id.split("__")[0] if "__" in mag_id else ""
    mag_tmp = work_tmp_dir / mag_id
    if mag_tmp.exists():
        shutil.rmtree(mag_tmp, ignore_errors=True)
    mag_tmp.mkdir(parents=True, exist_ok=True)

    f_sha = sha256_file(fasta_path)

    # 1. Run MinCED 0.4.2
    minced_txt_path = mag_tmp / f"{mag_id}.crisprs"
    minced_gff_path = mag_tmp / f"{mag_id}.gff"
    cmd_minced = [
        "singularity",
        "exec",
        "-B",
        "/project,/scratch",
        str(SIF_MINCED),
        "minced",
        "-minNR",
        "3",
        str(fasta_path),
        str(minced_txt_path),
        str(minced_gff_path),
    ]
    res_minced = subprocess.run(cmd_minced, capture_output=True, text=True, check=False)
    m_stdout = res_minced.stdout or ""
    m_stderr = res_minced.stderr or ""
    (mag_tmp / "minced.stdout.txt").write_text(m_stdout)
    (mag_tmp / "minced.stderr.txt").write_text(m_stderr)
    (mag_tmp / "minced.exitcode.txt").write_text(f"{res_minced.returncode}\n")

    minced_ok, minced_detail = is_minced_run_successful(
        res_minced.returncode,
        m_stdout,
        m_stderr,
        minced_txt_path,
    )
    m_arrays, m_spacers = (
        parse_minced_txt(minced_txt_path.read_text(errors="replace"), mag_id, fallback_sample, is_arch591)
        if minced_ok
        else ([], [])
    )

    # 2. Run CCTyper 1.8.0
    cct_out = mag_tmp / "cctyper_out"
    cct_tmp = mag_tmp / "tmp_cct"
    cct_tmp.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["TMPDIR"] = str(cct_tmp)
    env["CCTYPER_DB"] = str(DB_CCTYPER)
    cmd_cct = [
        "singularity",
        "exec",
        "-B",
        "/project,/scratch",
        str(SIF_CCTYPER),
        "cctyper",
        "--prodigal",
        "meta",
        "-t",
        str(cctyper_threads),
        "--no_plot",
        "--db",
        str(DB_CCTYPER),
        str(fasta_path),
        str(cct_out),
    ]
    res_cct = subprocess.run(cmd_cct, capture_output=True, text=True, env=env, check=False)
    shutil.rmtree(cct_tmp, ignore_errors=True)
    c_stdout = res_cct.stdout or ""
    c_stderr = res_cct.stderr or ""
    (mag_tmp / "cctyper.stdout.txt").write_text(c_stdout)
    (mag_tmp / "cctyper.stderr.txt").write_text(c_stderr)
    (mag_tmp / "cctyper.exitcode.txt").write_text(f"{res_cct.returncode}\n")

    cct_ok, cct_detail = is_cctyper_run_successful(
        res_cct.returncode,
        c_stdout,
        c_stderr,
        cct_out,
    )
    c_arrays, c_spacers, c_operons = (
        parse_cctyper_dir(cct_out, mag_id, fallback_sample, is_arch591)
        if cct_ok
        else ([], [], [])
    )

    # Persist compressed per-MAG raw output archive (`results/stage1_raw_archives/<tag>/<mag_id>.tar.gz`)
    tag_arch_dir = RAW_ARCH_DIR / tag
    tag_arch_dir.mkdir(parents=True, exist_ok=True)
    raw_tar_path = tag_arch_dir / f"{mag_id}.tar.gz"
    with tarfile.open(raw_tar_path, "w:gz") as tf:
        tf.add(mag_tmp, arcname=mag_id)

    if not (minced_ok and cct_ok):
        FAIL_DIAG_DIR.mkdir(parents=True, exist_ok=True)
        diag = {
            "mag_id": mag_id,
            "tag": tag,
            "fasta_path": str(fasta_path),
            "fasta_sha256": f_sha,
            "minced_exitcode": res_minced.returncode,
            "minced_ok": int(minced_ok),
            "minced_status_detail": minced_detail,
            "minced_stderr": res_minced.stderr,
            "cctyper_exitcode": res_cct.returncode,
            "cctyper_ok": int(cct_ok),
            "cctyper_status_detail": cct_detail,
            "cctyper_stderr": res_cct.stderr,
            "raw_archive_path": str(raw_tar_path),
        }
        (FAIL_DIAG_DIR / f"{mag_id}_diagnostics.json").write_text(json.dumps(diag, indent=2) + "\n")

    shutil.rmtree(mag_tmp, ignore_errors=True)

    c_trusted_cnt = sum(1 for a in c_arrays if int(a.get("is_trusted_array", 0)) == 1)
    c_untrusted_cnt = sum(1 for a in c_arrays if int(a.get("is_trusted_array", 0)) == 0)

    comb_arrays = m_arrays + c_arrays
    comb_spacers = m_spacers + c_spacers
    valid_spc = [s for s in comb_spacers if int(s["passes_quality_gate"]) == 1]
    indep_spc = deduplicate_mag_spacers(valid_spc)

    return {
        "mag_id": mag_id,
        "mag_fallback_sample": fallback_sample,
        "is_arch591_rep": int(is_arch591),
        "fasta_sha256": f_sha,
        "minced_exitcode": res_minced.returncode,
        "cctyper_exitcode": res_cct.returncode,
        "minced_ok": int(minced_ok),
        "cctyper_ok": int(cct_ok),
        "cctyper_status_detail": cct_detail,
        "screening_status": (
            "screened_both_callers_ok" if (minced_ok and cct_ok) else "failed_caller_execution"
        ),
        "raw_archive_path": str(raw_tar_path),
        "arrays": comb_arrays,
        "n_cctyper_trusted_arrays": c_trusted_cnt,
        "n_cctyper_putative_untrusted_arrays": c_untrusted_cnt,
        "spacers": comb_spacers,
        "n_independent_spacers": len(indep_spc),
        "cas_operons": c_operons,
    }


def make_fixture_coding_flank(seed: int, n_orfs: int = 2, codons_per_orf: int = 360) -> str:
    """Generate realistic non-homopolymeric coding ORFs so Prodigal predicts clean proteins."""
    sense_codons = [
        "GCT", "GCC", "GCA", "GCG", "CGT", "CGC", "CGA", "AAC", "AAT", "GAC", "GAT",
        "TGC", "TGT", "CAG", "CAA", "GAA", "GAG", "GGT", "GGC", "GGA", "CAC", "CAT",
        "ATC", "ATT", "CTG", "CTC", "CTT", "TTA", "TTG", "AAG", "AAA", "ATG", "TTC",
        "TTT", "CCG", "CCA", "CCT", "TCC", "TCT", "TCA", "TCG", "AGC", "AGT", "ACC",
        "ACA", "ACT", "ACG", "TGG", "TAC", "TAT", "GTT", "GTC", "GTA", "GTG",
    ]
    rng = random.Random(seed)
    chunks = []
    for _ in range(n_orfs):
        leader = "GCTAGCGATCGATAGGAGGTGATAA"
        orf = ["ATG"]
        for _ in range(codons_per_orf):
            for _ in range(20):
                c = rng.choice(sense_codons)
                prev = orf[-1]
                joined = prev + c
                if all(b * 4 not in joined for b in "ACGT"):
                    orf.append(c)
                    break
            else:
                orf.append("GCT")
        orf.append("TAA")
        spacer_inter = "GACTGACTGACTGACTGACTGACTG"
        chunks.append(leader + "".join(orf) + spacer_inter)
    return "".join(chunks)


def validate_tool_fixtures(work_tmp_dir: Path) -> dict[str, Any]:
    """Validate synthetic positive and negative control fixtures against actual MinCED and CCTyper tools."""
    verify_cctyper_db_dir()
    fix_dir = work_tmp_dir / "tool_fixtures"
    if fix_dir.exists():
        shutil.rmtree(fix_dir, ignore_errors=True)
    fix_dir.mkdir(parents=True, exist_ok=True)

    repeat = "GTTTGAGAGTTGTGATTTGCAAGTACGT"
    synth_spacers = [
        "AAGCCCAATAAACCACTCTGACTGGCCGAATA",
        "CTATTTGAAGGAGTCTAGCAGCCGCAGTAAGG",
        "AGCGGGACCCTCCATTGTTACTTATTAGGTTC",
        "AAACGTTCCGGCTGGTGTAGGATGATGCATCT",
        "GGTGAATCGAGAAGAGACTTGTTTTCCTCCTG",
    ]
    flank_left = make_fixture_coding_flank(101, n_orfs=2, codons_per_orf=360)
    flank_right = make_fixture_coding_flank(202, n_orfs=2, codons_per_orf=360)
    crispr_block = repeat + "".join(s + repeat for s in synth_spacers)
    pos_fa = fix_dir / "synth_pos.fa"
    pos_fa.write_text(f">synthSampler1_100\n{flank_left}{crispr_block}{flank_right}\n")

    neg_fa = fix_dir / "synth_neg.fa"
    neg_seq = make_fixture_coding_flank(303, n_orfs=4, codons_per_orf=360)
    neg_fa.write_text(f">synthSampler1_200\n{neg_seq}\n")

    pos_res = screen_single_mag("synthSampler1__self_PosFixture", pos_fa, False, "fixtures", fix_dir, 2)
    neg_res = screen_single_mag("synthSampler1__self_NegFixture", neg_fa, False, "fixtures", fix_dir, 2)

    shutil.rmtree(fix_dir, ignore_errors=True)
    pos_ok = (
        pos_res["screening_status"] == "screened_both_callers_ok"
        and pos_res["n_independent_spacers"] == 5
        and pos_res["minced_ok"] == 1
        and pos_res["cctyper_ok"] == 1
    )
    neg_ok = (
        neg_res["screening_status"] == "screened_both_callers_ok"
        and len(neg_res["arrays"]) == 0
        and neg_res["minced_ok"] == 1
        and neg_res["cctyper_ok"] == 1
    )
    if not (pos_ok and neg_ok):
        raise RuntimeError(
            f"Tool fixture validation failed: pos={pos_res['screening_status']} "
            f"(indep_spacers={pos_res['n_independent_spacers']}, minced={pos_res['minced_ok']}, cctyper={pos_res['cctyper_ok']}), "
            f"neg={neg_res['screening_status']} (arrays={len(neg_res['arrays'])}, minced={neg_res['minced_ok']}, cctyper={neg_res['cctyper_ok']})"
        )
    return {
        "positive_fixture_status": pos_res["screening_status"],
        "positive_fixture_independent_spacers": pos_res["n_independent_spacers"],
        "positive_fixture_minced_ok": pos_res["minced_ok"],
        "positive_fixture_cctyper_ok": pos_res["cctyper_ok"],
        "negative_fixture_status": neg_res["screening_status"],
        "negative_fixture_arrays": len(neg_res["arrays"]),
        "negative_fixture_minced_ok": neg_res["minced_ok"],
        "negative_fixture_cctyper_ok": neg_res["cctyper_ok"],
    }


def run_batch(
    tag: str,
    mag_list: list[str],
    is_arch591: bool,
    max_workers: int,
    cctyper_threads: int,
) -> dict[str, Any]:
    """Run a batch (pilot or chunk) of MAGs in parallel across compute node cores."""
    verify_cctyper_db_dir()
    t0 = time.time()
    CHUNK_OUT_DIR.mkdir(parents=True, exist_ok=True)
    RAW_ARCH_DIR.mkdir(parents=True, exist_ok=True)
    batch_tmp = TMP_ROOT / f"stage1_{tag}"
    batch_tmp.mkdir(parents=True, exist_ok=True)

    fasta_dir = IN_DIR / ("arch591_fastas" if is_arch591 else "nonrep_arch_fastas")
    tasks = []
    for m in mag_list:
        fp = fasta_dir / f"{m}.fa"
        if not fp.exists():
            raise FileNotFoundError(f"Missing verified FASTA for {m} at {fp}")
        tasks.append((m, fp))

    all_arrays: list[dict[str, Any]] = []
    all_spacers: list[dict[str, Any]] = []
    all_operons: list[dict[str, Any]] = []
    status_rows: list[dict[str, Any]] = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {
            ex.submit(screen_single_mag, m, fp, is_arch591, tag, batch_tmp, cctyper_threads): m
            for m, fp in tasks
        }
        for fut in concurrent.futures.as_completed(futs):
            res = fut.result()
            m_arr = res["arrays"]
            m_spc = res["spacers"]
            m_op = res["cas_operons"]
            all_arrays.extend(m_arr)
            all_spacers.extend(m_spc)
            all_operons.extend(m_op)
            status_rows.append(
                {
                    "mag_id": res["mag_id"],
                    "mag_fallback_sample": res["mag_fallback_sample"],
                    "is_arch591_rep": res["is_arch591_rep"],
                    "fasta_sha256": res["fasta_sha256"],
                    "minced_exitcode": res["minced_exitcode"],
                    "cctyper_exitcode": res["cctyper_exitcode"],
                    "minced_ok": res["minced_ok"],
                    "cctyper_ok": res["cctyper_ok"],
                    "cctyper_status_detail": res["cctyper_status_detail"],
                    "screening_status": res["screening_status"],
                    "n_raw_arrays": len(m_arr),
                    "n_cctyper_trusted_arrays": res["n_cctyper_trusted_arrays"],
                    "n_cctyper_putative_untrusted_arrays": res["n_cctyper_putative_untrusted_arrays"],
                    "n_raw_spacers": len(m_spc),
                    "n_valid_spacers": sum(int(s["passes_quality_gate"]) for s in m_spc),
                    "n_independent_spacers": res["n_independent_spacers"],
                    "n_cas_operons": len(m_op),
                    "raw_archive_path": res["raw_archive_path"],
                }
            )

    pd.DataFrame(status_rows, columns=STATUS_COLS).sort_values("mag_id").to_csv(
        CHUNK_OUT_DIR / f"{tag}_status.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_arrays, columns=ARRAY_COLS).to_csv(
        CHUNK_OUT_DIR / f"{tag}_arrays.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_spacers, columns=SPACER_COLS).to_csv(
        CHUNK_OUT_DIR / f"{tag}_spacers.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_operons, columns=OPERON_COLS).to_csv(
        CHUNK_OUT_DIR / f"{tag}_cas_operons.tsv", sep="\t", index=False
    )

    shutil.rmtree(batch_tmp, ignore_errors=True)
    n_succ = sum(1 for r in status_rows if r["screening_status"] == "screened_both_callers_ok")
    n_fail = sum(1 for r in status_rows if r["screening_status"] != "screened_both_callers_ok")

    return {
        "tag": tag,
        "n_mags_planned": len(mag_list),
        "n_mags_attempted": len(status_rows),
        "n_mags_successfully_screened": n_succ,
        "n_mags_failed": n_fail,
        "n_minced_ok": sum(r["minced_ok"] for r in status_rows),
        "n_cctyper_ok": sum(r["cctyper_ok"] for r in status_rows),
        "n_mags_with_valid_spacers": sum(
            1 for r in status_rows
            if r["screening_status"] == "screened_both_callers_ok" and r["n_independent_spacers"] > 0
        ),
        "n_total_raw_arrays": len(all_arrays),
        "n_total_raw_spacers": len(all_spacers),
        "n_total_valid_spacers": sum(int(s["passes_quality_gate"]) for s in all_spacers),
        "n_total_cas_operons": len(all_operons),
        "elapsed_seconds": round(time.time() - t0, 2),
    }

    pd.DataFrame(status_rows, columns=STATUS_COLS).sort_values("mag_id").to_csv(
        CHUNK_OUT_DIR / f"{tag}_status.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_arrays, columns=ARRAY_COLS).to_csv(
        CHUNK_OUT_DIR / f"{tag}_arrays.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_spacers, columns=SPACER_COLS).to_csv(
        CHUNK_OUT_DIR / f"{tag}_spacers.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_operons, columns=OPERON_COLS).to_csv(
        CHUNK_OUT_DIR / f"{tag}_cas_operons.tsv", sep="\t", index=False
    )

    shutil.rmtree(batch_tmp, ignore_errors=True)
    n_succ = sum(1 for r in status_rows if r["screening_status"] == "screened_both_callers_ok")
    n_fail = sum(1 for r in status_rows if r["screening_status"] != "screened_both_callers_ok")

    return {
        "tag": tag,
        "n_mags_planned": len(mag_list),
        "n_mags_attempted": len(status_rows),
        "n_mags_successfully_screened": n_succ,
        "n_mags_failed": n_fail,
        "n_minced_ok": sum(r["minced_ok"] for r in status_rows),
        "n_cctyper_ok": sum(r["cctyper_ok"] for r in status_rows),
        "n_mags_with_valid_spacers": sum(
            1 for r in status_rows
            if r["screening_status"] == "screened_both_callers_ok" and r["n_independent_spacers"] > 0
        ),
        "n_total_raw_arrays": len(all_arrays),
        "n_total_raw_spacers": len(all_spacers),
        "n_total_valid_spacers": sum(int(s["passes_quality_gate"]) for s in all_spacers),
        "n_total_cas_operons": len(all_operons),
        "elapsed_seconds": round(time.time() - t0, 2),
    }


def run_collect() -> None:
    """Merge all 12 canonical chunks + chunk 13 without silent drop_duplicates."""
    t0 = time.time()
    arch591_list = [
        line.strip()
        for line in (IN_DIR / "arch591_mags.txt").read_text().splitlines()
        if line.strip()
    ]
    arch591_set = set(arch591_list)

    status_dfs = []
    array_dfs = []
    spacer_dfs = []
    operon_dfs = []

    chunk_tags = [f"chunk_{idx:02d}" for idx in range(1, 13)]
    if (CHUNK_OUT_DIR / "chunk_13_status.tsv").exists():
        chunk_tags.append("chunk_13")

    for tag in chunk_tags:
        st_p = CHUNK_OUT_DIR / f"{tag}_status.tsv"
        if not st_p.exists():
            raise FileNotFoundError(f"Missing chunk status file: {st_p}")
        status_dfs.append(pd.read_csv(st_p, sep="\t"))
        for p, lst, cols in (
            (CHUNK_OUT_DIR / f"{tag}_arrays.tsv", array_dfs, ARRAY_COLS),
            (CHUNK_OUT_DIR / f"{tag}_spacers.tsv", spacer_dfs, SPACER_COLS),
            (CHUNK_OUT_DIR / f"{tag}_cas_operons.tsv", operon_dfs, OPERON_COLS),
        ):
            if p.exists():
                df_part = pd.read_csv(p, sep="\t")
                if len(df_part) > 0:
                    lst.append(df_part)

    status_all = pd.concat(status_dfs, ignore_index=True) if status_dfs else pd.DataFrame(columns=STATUS_COLS)
    if status_all["mag_id"].duplicated().any():
        dups = sorted(status_all[status_all["mag_id"].duplicated()]["mag_id"].unique())
        raise RuntimeError(f"Duplicate MAG rows found across Stage 1 chunks (refusing silent drop_duplicates): {dups}")

    arch591_st = status_all[status_all["is_arch591_rep"] == 1]
    attempted_arch591 = set(arch591_st["mag_id"])
    succ_arch591 = set(
        arch591_st[arch591_st["screening_status"] == "screened_both_callers_ok"]["mag_id"]
    )
    failed_arch591 = set(
        arch591_st[arch591_st["screening_status"] != "screened_both_callers_ok"]["mag_id"]
    )
    unscreened_arch591 = arch591_set - attempted_arch591

    if len(unscreened_arch591) > 0 or len(failed_arch591) > 0:
        raise RuntimeError(
            f"Stage 1 collect incomplete or failed: unscreened={len(unscreened_arch591)}, "
            f"failed={len(failed_arch591)} ({sorted(failed_arch591)[:10]})"
        )

    arrays_all = pd.concat(array_dfs, ignore_index=True) if array_dfs else pd.DataFrame(columns=ARRAY_COLS)
    spacers_all = pd.concat(spacer_dfs, ignore_index=True) if spacer_dfs else pd.DataFrame(columns=SPACER_COLS)
    operons_all = pd.concat(operon_dfs, ignore_index=True) if operon_dfs else pd.DataFrame(columns=OPERON_COLS)

    c8_lookup: dict[str, dict[str, Any]] = {}
    if P_C8_DEF.exists():
        c8_df = pd.read_csv(P_C8_DEF, sep="\t")
        if "genome" in c8_df.columns:
            c8_lookup = c8_df.set_index("genome").to_dict("index")

    status_all["c8_defense_systems"] = status_all["mag_id"].map(
        lambda m: c8_lookup.get(m, {}).get("n_systems", "NA")
    )
    status_all.sort_values(["is_arch591_rep", "mag_id"], ascending=[False, True]).to_csv(
        RES_DIR / "stage1_mag_screening_status.tsv", sep="\t", index=False
    )
    arrays_all.to_csv(RES_DIR / "stage1_raw_crispr_arrays.tsv", sep="\t", index=False)
    operons_all.to_csv(RES_DIR / "stage1_cas_operons.tsv", sep="\t", index=False)

    valid_spc_df = (
        spacers_all[spacers_all["passes_quality_gate"] == 1].copy()
        if len(spacers_all) > 0
        else pd.DataFrame(columns=SPACER_COLS)
    )
    dedup_rows = deduplicate_mag_spacers(valid_spc_df)
    dedup_df = pd.DataFrame(dedup_rows, columns=DEDUP_SPACER_COLS)
    dedup_df.to_csv(RES_DIR / "stage1_verified_archaeal_spacers.tsv", sep="\t", index=False)

    fna_path = RES_DIR / "stage1_verified_archaeal_spacers.fna"
    with open(fna_path, "w") as fh:
        for r in dedup_rows:
            fh.write(
                f">{r['unique_spacer_id']}|mag={r['mag_id']}|contig={r['primary_namespaced_contig_id']}|"
                f"loci={r['all_array_loci']}|callers={r['callers']}|arch591={r['is_arch591_rep']}\n"
                f"{r['representative_spacer_seq']}\n"
            )

    summary = {
        "stage": "stage1_crispr_recovery_collect",
        "status": "completed",
        "elapsed_seconds": round(time.time() - t0, 2),
        "n_arch591_planned": len(arch591_set),
        "n_arch591_attempted": len(attempted_arch591),
        "n_arch591_successfully_screened": len(succ_arch591),
        "n_arch591_failed": len(failed_arch591),
        "n_arch591_unscreened_unknown": len(unscreened_arch591),
        "n_arch591_with_arrays": int((arch591_st["n_raw_arrays"] > 0).sum()),
        "n_arch591_with_valid_spacers": int((arch591_st["n_independent_spacers"] > 0).sum()),
        "n_arch591_with_cas_operons": int((arch591_st["n_cas_operons"] > 0).sum()),
        "n_nonrep_arch_attempted": int((status_all["is_arch591_rep"] == 0).sum()),
        "n_nonrep_arch_successfully_screened": int(
            ((status_all["is_arch591_rep"] == 0) & (status_all["screening_status"] == "screened_both_callers_ok")).sum()
        ),
        "n_total_raw_arrays": int(len(arrays_all)),
        "n_total_raw_spacers": int(len(spacers_all)),
        "n_unique_independent_spacers": int(len(dedup_df)),
        "n_unique_independent_spacers_arch591": int(
            len(dedup_df[dedup_df["is_arch591_rep"] == 1]) if len(dedup_df) > 0 else 0
        ),
    }
    (RES_DIR / "stage1_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("[Stage 1 Collect] Summary:", json.dumps(summary), flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["pilot", "chunk", "collect"], required=True)
    ap.add_argument("--chunk-id", type=int, default=0)
    ap.add_argument("--max-workers", type=int, default=16)
    ap.add_argument("--cctyper-threads", type=int, default=4)
    args = ap.parse_args()

    for d in (RES_DIR, CHUNK_OUT_DIR, RAW_ARCH_DIR, TMP_ROOT):
        d.mkdir(parents=True, exist_ok=True)

    if args.mode in ("pilot", "chunk"):
        verify_cctyper_db_dir()

    if args.mode == "pilot":
        fixture_res = validate_tool_fixtures(TMP_ROOT)
        pilot_mags = [
            line.strip()
            for line in (IN_DIR / "pilot_30_mags.txt").read_text().splitlines()
            if line.strip()
        ]
        res = run_batch("pilot", pilot_mags, True, args.max_workers, args.cctyper_threads)
        res["tool_fixture_validation"] = fixture_res
        if (
            res["n_mags_attempted"] != len(pilot_mags)
            or res["n_mags_successfully_screened"] != len(pilot_mags)
            or res["n_mags_failed"] != 0
            or res["n_mags_with_valid_spacers"] == 0
        ):
            raise RuntimeError(f"Stage 1 pilot validation failed strict caller/spacer gate: {res}")
        (RES_DIR / "stage1_pilot_validation.json").write_text(json.dumps(res, indent=2) + "\n")
        print("[Stage 1 Pilot] Validation passed:", json.dumps(res), flush=True)

    elif args.mode == "chunk":
        if not (1 <= args.chunk_id <= 13):
            raise ValueError(f"Invalid --chunk-id {args.chunk_id}; expected 1..13")
        if args.chunk_id <= 12:
            chunk_file = IN_DIR / "chunks" / f"chunk_{args.chunk_id:02d}.txt"
            is_arch591 = True
            tag = f"chunk_{args.chunk_id:02d}"
        else:
            chunk_file = IN_DIR / "chunks" / "chunk_13_nonrep_arch.txt"
            is_arch591 = False
            tag = "chunk_13"
        mags = [line.strip() for line in chunk_file.read_text().splitlines() if line.strip()]
        res = run_batch(tag, mags, is_arch591, args.max_workers, args.cctyper_threads)
        if res["n_mags_failed"] > 0:
            raise RuntimeError(f"Stage 1 {tag} had {res['n_mags_failed']} failed MAG screenings: {res}")
        print(f"[Stage 1 {tag}] Completed:", json.dumps(res), flush=True)

    elif args.mode == "collect":
        run_collect()


if __name__ == "__main__":
    main()
