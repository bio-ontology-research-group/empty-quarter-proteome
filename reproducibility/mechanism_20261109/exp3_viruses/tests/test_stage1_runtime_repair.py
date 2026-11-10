#!/usr/bin/env python3
"""Stage 1 Runtime Repair Verification Tests.

Verifies:
1. CCTyper DB integrity checking (missing files / non-HMM files rejected).
2. MinCED and CCTyper exit code vs log error vs completion marker validation.
3. CCTyper Trusted vs putative/untrusted array parsing and spacer quality gating.
4. MinCED + CCTyper non-independence audit annotation in spacer deduplication.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

CODE_DIR = Path(__file__).resolve().parent.parent / "code"
sys.path.insert(0, str(CODE_DIR))

from stage1_crispr_recovery import (
    deduplicate_mag_spacers,
    is_cctyper_run_successful,
    is_minced_run_successful,
    parse_cctyper_dir,
    verify_cctyper_db_dir,
)


def test_verify_cctyper_db_dir_checks(tmp_path: Path):
    """Verify that verify_cctyper_db_dir rejects missing files and non-HMM profile entries."""
    db_dir = tmp_path / "mock_db"
    db_dir.mkdir()

    # Missing required files -> raise
    with pytest.raises(RuntimeError, match="missing or empty required files"):
        verify_cctyper_db_dir(db_dir)

    # Create required files
    for fname in (
        "CasScoring.csv",
        "xgb_repeats.model",
        "type_dict.tab",
        "cutoffs.tab",
        "interference.json",
        "adaptation.json",
        "repeats.fa",
    ):
        (db_dir / fname).write_text("data\n")

    # Missing Profiles dir -> raise
    with pytest.raises(RuntimeError, match="HMM profile directory missing"):
        verify_cctyper_db_dir(db_dir)

    prof_dir = db_dir / "Profiles"
    prof_dir.mkdir()

    # Empty Profiles dir -> raise
    with pytest.raises(RuntimeError, match="HMM profile directory is empty"):
        verify_cctyper_db_dir(db_dir)

    # Non-.hmm files in Profiles -> raise
    (prof_dir / "test.txt").write_text("not an hmm\n")
    with pytest.raises(RuntimeError, match="contains non-.hmm files"):
        verify_cctyper_db_dir(db_dir)

    # With only .hmm files -> pass
    (prof_dir / "test.txt").unlink()
    (prof_dir / "Cas1.hmm").write_text("HMMER3/f\n")
    res = verify_cctyper_db_dir(db_dir)
    assert res["n_required_files_verified"] == 7
    assert res["n_hmm_profiles"] == 1


def test_is_minced_run_successful(tmp_path: Path):
    """Verify MinCED success validation catches exit codes, missing files, and log exceptions."""
    out_txt = tmp_path / "sample.crisprs"

    # Nonzero exit code
    ok, reason = is_minced_run_successful(1, "", "", out_txt)
    assert not ok
    assert reason == "nonzero_exit_1"

    # Missing file
    ok, reason = is_minced_run_successful(0, "", "", out_txt)
    assert not ok
    assert reason == "missing_minced_crisprs_file"

    # Exception in log
    out_txt.write_text("CRISPR 1\n")
    ok, reason = is_minced_run_successful(0, "", "Exception in thread java.lang.NullPointerException", out_txt)
    assert not ok
    assert reason == "minced_log_exception_or_error"

    # Clean run
    ok, reason = is_minced_run_successful(0, "", "Finished", out_txt)
    assert ok
    assert reason == "ok"


def test_is_cctyper_run_successful(tmp_path: Path):
    """Verify CCTyper exit-0 error detection catches missing DB and missing completion markers."""
    cct_out = tmp_path / "cct_out"

    # Missing directory
    ok, reason = is_cctyper_run_successful(0, "", "", cct_out)
    assert not ok
    assert reason == "missing_cctyper_output_dir"

    cct_out.mkdir()

    # Exit code 0 but stderr has DB error (exact failure seen in pilot 15871683)
    stderr_fail = (
        "[2026-11-09 10:44:47] INFO: Running CRISPRCasTyper version 1.8.0\n"
        "[2026-11-09 10:44:47] ERROR: Could not find database directory\n"
    )
    ok, reason = is_cctyper_run_successful(0, "", stderr_fail, cct_out)
    assert not ok
    assert reason == "cctyper_log_error_or_critical"

    # Missing required tabs
    ok, reason = is_cctyper_run_successful(0, "", "Some log", cct_out)
    assert not ok
    assert reason == "missing_or_empty_arguments_tab"

    (cct_out / "arguments.tab").write_text("input:\ttest.fa\n")
    (cct_out / "genes.tab").write_text("Contig\tStart\tEnd\tStrand\n")
    (cct_out / "hmmer.tab").write_text("Hmm\tORF\n")

    # Missing completion clean marker
    ok, reason = is_cctyper_run_successful(0, "Done", "", cct_out)
    assert not ok
    assert reason == "missing_final_clean_completion_marker"

    # Genuine complete run with master.clean() marker
    stdout_ok = "Running CRISPRCasTyper\nRemoving temporary files\n"
    ok, reason = is_cctyper_run_successful(0, stdout_ok, "", cct_out)
    assert ok
    assert reason == "ok"


def test_cctyper_trusted_vs_untrusted_array_parsing(tmp_path: Path):
    """Verify parse_cctyper_dir preserves Trusted vs Putative untrusted arrays."""
    cct_out = tmp_path / "cct_out_test"
    cct_out.mkdir()
    spc_dir = cct_out / "spacers"
    spc_dir.mkdir()

    # Write crisprs_all.tab with 1 Trusted and 1 Untrusted array
    crisprs_tab = cct_out / "crisprs_all.tab"
    crisprs_tab.write_text(
        "Contig\tCRISPR\tStart\tEnd\tConsensus_repeat\tN_repeats\tRepeat_len\tSpacer_len_avg\tRepeat_identity\tSpacer_identity\tSpacer_len_sem\tTrusted\tPrediction\n"
        "contig_1\tcontig_1_1\t100\t300\tGTTTGAGAGTTGTGATTTGCAAGTACGT\t3\t28\t32.0\t100.0\t35.0\t0.0\tTrue\tI-B\n"
        "contig_1\tcontig_1_2\t5000\t5200\tGTTTGAGAGTTGTGATTTGCAAGTACGT\t3\t28\t32.0\t80.0\t70.0\t4.0\tFalse\tUnknown\n"
    )

    # Valid spacer sequences
    spc1 = "AAGCCCAATAAACCACTCTGACTGGCCGAATA"
    spc2 = "CTATTTGAAGGAGTCTAGCAGCCGCAGTAAGG"
    (spc_dir / "contig_1_1.fa").write_text(f">contig_1_1:1\n{spc1}\n>contig_1_1:2\n{spc2}\n")
    (spc_dir / "contig_1_2.fa").write_text(f">contig_1_2:1\n{spc1}\n>contig_1_2:2\n{spc2}\n")

    arrays, spacers, _ = parse_cctyper_dir(cct_out, "magA", "sample1", True)

    assert len(arrays) == 2
    trusted_arr = [a for a in arrays if a["array_id"].endswith("_100_300")][0]
    untrusted_arr = [a for a in arrays if a["array_id"].endswith("_5000_5200")][0]

    assert trusted_arr["is_trusted_array"] == 1
    assert trusted_arr["array_trust_tier"] == "cctyper_trusted"
    assert trusted_arr["n_valid_spacers"] == 2

    assert untrusted_arr["is_trusted_array"] == 0
    assert untrusted_arr["array_trust_tier"] == "cctyper_putative_untrusted"
    assert untrusted_arr["n_valid_spacers"] == 0

    assert len(spacers) == 4
    trusted_spcs = [s for s in spacers if s["array_id"].endswith("_100_300")]
    untrusted_spcs = [s for s in spacers if s["array_id"].endswith("_5000_5200")]

    for s in trusted_spcs:
        assert s["passes_quality_gate"] == 1
        assert s["is_trusted_array"] == 1

    for s in untrusted_spcs:
        assert s["passes_seq_quality_gate"] == 1
        assert s["passes_quality_gate"] == 0
        assert s["quality_gate_reason"] == "cctyper_putative_untrusted_array_exploratory_only"


def test_minced_and_cctyper_non_independence_dedup():
    """Verify MinCED and CCTyper calling same spacer marks caller_independence_note and dedups once."""
    spc = "AAGCCCAATAAACCACTCTGACTGGCCGAATA"
    spacer_rows = [
        {
            "spacer_id": "spc_m1",
            "array_id": "arr_m",
            "mag_id": "mag1",
            "is_arch591_rep": 1,
            "raw_contig_id": "c1",
            "contig_sample": "samp",
            "namespaced_contig_id": "samp::c1",
            "array_start": 100,
            "array_end": 300,
            "caller": "minced_0.4.2",
            "is_trusted_array": 1,
            "array_trust_tier": "minced_repeat_array",
            "spacer_index": 1,
            "repeat_seq": "GTTTGAGAGTTGTGATTTGCAAGTACGT",
            "spacer_seq": spc,
            "canonical_spacer_seq_rc": spc,
            "spacer_orientation": "fwd",
            "spacer_len": len(spc),
            "passes_seq_quality_gate": 1,
            "passes_quality_gate": 1,
            "quality_gate_reason": "valid_crispr_spacer",
        },
        {
            "spacer_id": "spc_c1",
            "array_id": "arr_c",
            "mag_id": "mag1",
            "is_arch591_rep": 1,
            "raw_contig_id": "c1",
            "contig_sample": "samp",
            "namespaced_contig_id": "samp::c1",
            "array_start": 100,
            "array_end": 300,
            "caller": "cctyper_1.8.0",
            "is_trusted_array": 1,
            "array_trust_tier": "cctyper_trusted",
            "spacer_index": 1,
            "repeat_seq": "GTTTGAGAGTTGTGATTTGCAAGTACGT",
            "spacer_seq": spc,
            "canonical_spacer_seq_rc": spc,
            "spacer_orientation": "fwd",
            "spacer_len": len(spc),
            "passes_seq_quality_gate": 1,
            "passes_quality_gate": 1,
            "quality_gate_reason": "valid_crispr_spacer",
        },
    ]

    dedup = deduplicate_mag_spacers(spacer_rows)
    assert len(dedup) == 1
    assert dedup[0]["callers"] == "cctyper_1.8.0+minced_0.4.2"
    assert dedup[0]["caller_independence_note"] == "cctyper_uses_minced_not_independent"
    assert dedup[0]["n_loci_in_mag"] == 1
    assert dedup[0]["n_distinct_arrays_in_mag"] == 1
