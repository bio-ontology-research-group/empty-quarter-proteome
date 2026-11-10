from __future__ import annotations

import shutil
import sys
import uuid
from pathlib import Path
from typing import Any, Generator

import pytest

EXP3_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = EXP3_ROOT / "code"
TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE_DIR))

import stage4_host_array_provenance_audit as host_prov
from stage4_host_array_provenance_audit import (
    _strict_int_from_tsv,
    evaluate_host_array_contig_coverage_from_mapper,
)


@pytest.fixture
def project_scratch() -> Generator[Path, None, None]:
    """Isolated project-local scratch directory under tests/scratch (never /tmp)."""
    scratch = TESTS_DIR / "scratch" / f"host_cov_{uuid.uuid4().hex[:8]}"
    scratch.mkdir(parents=True, exist_ok=True)
    try:
        yield scratch
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _build_representative_header_fixture(
    scratch: Path,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Reuses the representative-header fixture pattern from test_viral_mapping_repair.py."""
    cov_tsv = scratch / "41PRr2_filtered_coverage.tsv"
    cov_tsv.write_text(
        "#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
        "HOST__41PRr2__single_SemiBin_236__41PRr2_1588631\t1\t651471\t12000\t600000\t92.099265\t2.75\t35.1\t41.2\n"
        "HOST__41PRr2__single_SemiBin_236__41PRr2_other1\t1\t348529\t6000\t300000\t86.076051\t2.50\t35.0\t40.8\n",
        encoding="utf-8",
    )
    parsed = host_prov.parse_samtools_filtered_coverage_tsv(cov_tsv)

    host_cw = [
        {
            "source_category": "host_mag",
            "entity_id": "41PRr2__single_SemiBin_236",
            "contig_id": "41PRr2_1588631",
            "representative_header": "HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
            "seq_length": "651471",
        },
        {
            "source_category": "host_mag",
            "entity_id": "41PRr2__single_SemiBin_236",
            "contig_id": "41PRr2_other1",
            "representative_header": "HOST__41PRr2__single_SemiBin_236__41PRr2_other1",
            "seq_length": "348529",
        },
    ]
    return parsed, host_cw


def test_evaluate_host_array_contig_coverage_accepts_integer_valued_decimal_strings(project_scratch: Path):
    """Regression test for job 15910427 failure: int('17.0').
    Validates that evaluate_host_array_contig_coverage_from_mapper correctly parses
    integer-valued decimal strings (e.g. '17.0', '348529.0', '6000.0') exported by pandas
    mixed-entity TSVs into true Python integers.
    """
    parsed, host_cw = _build_representative_header_fixture(project_scratch)

    # Job 15910427 failure pattern: decimal string '17.0' in unique_mag_contigs
    eval_row = {
        "unique_mag_contigs": "17.0",
        "unique_mag_bases": "348529.0",
        "unique_filtered_reads_mapped": "6000.0",
        "unique_horizontal_breadth_pct": "30.0",
        "unique_mean_depth": "0.871323",
    }

    res = evaluate_host_array_contig_coverage_from_mapper(
        host_mag="41PRr2__single_SemiBin_236",
        array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
        array_contig_length_bp=651471,
        mag_total_contigs=2,
        mag_total_bases=1000000,
        host_crosswalk_rows=host_cw,
        coverage_by_rname=parsed,
        focal_eval_host_row=eval_row,
    )

    # Assert exact integer values
    assert res["mapper_raw_conservative_unique_mag_contigs"] == 17
    assert res["mapper_raw_conservative_unique_mag_bases"] == 348529
    assert res["mapper_raw_conservative_unique_mag_reads"] == 6000

    # Assert exact Python int type (not float and not bool)
    assert type(res["mapper_raw_conservative_unique_mag_contigs"]) is int
    assert type(res["mapper_raw_conservative_unique_mag_bases"]) is int
    assert type(res["mapper_raw_conservative_unique_mag_reads"]) is int

    # Also test zero and standard integer representations
    eval_row_zero = {
        "unique_mag_contigs": "0.0",
        "unique_mag_bases": "0.0",
        "unique_filtered_reads_mapped": "0",
        "unique_horizontal_breadth_pct": "0.0",
        "unique_mean_depth": "0.0",
    }
    res_zero = evaluate_host_array_contig_coverage_from_mapper(
        host_mag="41PRr2__single_SemiBin_236",
        array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
        array_contig_length_bp=651471,
        mag_total_contigs=2,
        mag_total_bases=1000000,
        host_crosswalk_rows=host_cw,
        coverage_by_rname=parsed,
        focal_eval_host_row=eval_row_zero,
    )
    assert res_zero["mapper_raw_conservative_unique_mag_contigs"] == 0
    assert res_zero["mapper_raw_conservative_unique_mag_bases"] == 0
    assert res_zero["mapper_raw_conservative_unique_mag_reads"] == 0
    assert type(res_zero["mapper_raw_conservative_unique_mag_contigs"]) is int


@pytest.mark.parametrize("col_name", ["unique_mag_contigs", "unique_mag_bases", "unique_filtered_reads_mapped"])
@pytest.mark.parametrize("fractional_val", ["17.5", "17.0001", "0.5", 17.5, 0.1])
def test_evaluate_host_array_contig_coverage_rejects_fractional_floats(
    project_scratch: Path,
    col_name: str,
    fractional_val: Any,
):
    """Rejects fractional count values without silent truncation across all conservative metrics."""
    parsed, host_cw = _build_representative_header_fixture(project_scratch)

    base_eval = {
        "unique_mag_contigs": "17.0",
        "unique_mag_bases": "348529.0",
        "unique_filtered_reads_mapped": "6000.0",
        "unique_horizontal_breadth_pct": "30.0",
        "unique_mean_depth": "0.871323",
    }
    bad_eval = {**base_eval, col_name: fractional_val}

    with pytest.raises(ValueError, match="non-integer float"):
        evaluate_host_array_contig_coverage_from_mapper(
            host_mag="41PRr2__single_SemiBin_236",
            array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
            array_contig_length_bp=651471,
            mag_total_contigs=2,
            mag_total_bases=1000000,
            host_crosswalk_rows=host_cw,
            coverage_by_rname=parsed,
            focal_eval_host_row=bad_eval,
        )


@pytest.mark.parametrize("col_name", ["unique_mag_contigs", "unique_mag_bases", "unique_filtered_reads_mapped"])
@pytest.mark.parametrize("nonfinite_val", ["nan", "NaN", "inf", "-inf", float("nan"), float("inf")])
def test_evaluate_host_array_contig_coverage_rejects_nonfinite_values(
    project_scratch: Path,
    col_name: str,
    nonfinite_val: Any,
):
    """Rejects non-finite values (NaN, Inf) across conservative count metrics."""
    parsed, host_cw = _build_representative_header_fixture(project_scratch)

    base_eval = {
        "unique_mag_contigs": "17.0",
        "unique_mag_bases": "348529.0",
        "unique_filtered_reads_mapped": "6000.0",
        "unique_horizontal_breadth_pct": "30.0",
        "unique_mean_depth": "0.871323",
    }
    bad_eval = {**base_eval, col_name: nonfinite_val}

    with pytest.raises(ValueError, match="empty/NaN|non-integer|non-finite"):
        evaluate_host_array_contig_coverage_from_mapper(
            host_mag="41PRr2__single_SemiBin_236",
            array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
            array_contig_length_bp=651471,
            mag_total_contigs=2,
            mag_total_bases=1000000,
            host_crosswalk_rows=host_cw,
            coverage_by_rname=parsed,
            focal_eval_host_row=bad_eval,
        )


@pytest.mark.parametrize("col_name", ["unique_mag_contigs", "unique_mag_bases", "unique_filtered_reads_mapped"])
@pytest.mark.parametrize("bool_val", [True, False])
def test_evaluate_host_array_contig_coverage_rejects_booleans(
    project_scratch: Path,
    col_name: str,
    bool_val: bool,
):
    """Rejects boolean values across conservative count metrics."""
    parsed, host_cw = _build_representative_header_fixture(project_scratch)

    base_eval = {
        "unique_mag_contigs": "17.0",
        "unique_mag_bases": "348529.0",
        "unique_filtered_reads_mapped": "6000.0",
        "unique_horizontal_breadth_pct": "30.0",
        "unique_mean_depth": "0.871323",
    }
    bad_eval = {**base_eval, col_name: bool_val}

    with pytest.raises(ValueError, match="missing or boolean value"):
        evaluate_host_array_contig_coverage_from_mapper(
            host_mag="41PRr2__single_SemiBin_236",
            array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
            array_contig_length_bp=651471,
            mag_total_contigs=2,
            mag_total_bases=1000000,
            host_crosswalk_rows=host_cw,
            coverage_by_rname=parsed,
            focal_eval_host_row=bad_eval,
        )


@pytest.mark.parametrize("col_name", ["unique_mag_contigs", "unique_mag_bases", "unique_filtered_reads_mapped"])
@pytest.mark.parametrize("missing_val", [None, "", "   "])
def test_evaluate_host_array_contig_coverage_rejects_missing_values(
    project_scratch: Path,
    col_name: str,
    missing_val: Any,
):
    """Rejects missing or empty count values across conservative metrics."""
    parsed, host_cw = _build_representative_header_fixture(project_scratch)

    base_eval = {
        "unique_mag_contigs": "17.0",
        "unique_mag_bases": "348529.0",
        "unique_filtered_reads_mapped": "6000.0",
        "unique_horizontal_breadth_pct": "30.0",
        "unique_mean_depth": "0.871323",
    }
    bad_eval = {**base_eval, col_name: missing_val}

    with pytest.raises(ValueError, match="missing|empty"):
        evaluate_host_array_contig_coverage_from_mapper(
            host_mag="41PRr2__single_SemiBin_236",
            array_rep_header="HOST__41PRr2__single_SemiBin_236__41PRr2_1588631",
            array_contig_length_bp=651471,
            mag_total_contigs=2,
            mag_total_bases=1000000,
            host_crosswalk_rows=host_cw,
            coverage_by_rname=parsed,
            focal_eval_host_row=bad_eval,
        )


def test_strict_int_from_tsv_comprehensive_boundaries():
    """Unit-level boundary verification for _strict_int_from_tsv."""
    # Integer-valued decimal floats and strings
    assert _strict_int_from_tsv("17.0", "c") == 17
    assert _strict_int_from_tsv("17", "c") == 17
    assert _strict_int_from_tsv(17, "c") == 17
    assert _strict_int_from_tsv(17.0, "c") == 17
    assert _strict_int_from_tsv("0.0", "c") == 0
    assert _strict_int_from_tsv("100.000", "c") == 100
    assert _strict_int_from_tsv("-5.0", "c") == -5

    # Fractional floats must raise without silent truncation
    for bad in ["17.5", 17.5, "17.0001", "0.1", 0.5, "-0.1"]:
        with pytest.raises(ValueError, match="non-integer float"):
            _strict_int_from_tsv(bad, "c")

    # Non-finite values
    for bad in ["nan", "NaN", "inf", "-inf", float("nan"), float("inf")]:
        with pytest.raises(ValueError, match="empty/NaN|non-integer"):
            _strict_int_from_tsv(bad, "c")

    # Booleans
    for bad in [True, False]:
        with pytest.raises(ValueError, match="missing or boolean value"):
            _strict_int_from_tsv(bad, "c")

    # Missing / empty
    for bad in [None, "", "   "]:
        with pytest.raises(ValueError, match="missing|empty"):
            _strict_int_from_tsv(bad, "c")

    # Non-numeric strings
    for bad in ["abc", "17a", "17.0.0"]:
        with pytest.raises(ValueError, match="non-integer value"):
            _strict_int_from_tsv(bad, "c")
