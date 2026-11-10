from __future__ import annotations

import gzip
import shutil
import sys
import uuid
from pathlib import Path
from typing import Generator

import pandas as pd
import pytest

EXP3_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = EXP3_ROOT / "code"
sys.path.insert(0, str(CODE_DIR))

from stage4_posthoc_instrain_and_mapping_followup import (
    INSTRAIN_PROFILE_FINISHED_BANNER,
    qualify_instrain_profile_output,
)

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures" / "root_review_v60_variants"
ACTUAL_V60_PROFILE = FIXTURES_DIR / "actual_profile"
ACTUAL_V60_REPRESENTATIVE = "VIRUS_TARGET__V60PRr2::V60PRr2_715976"


@pytest.fixture
def isolated_scratch() -> Generator[Path, None, None]:
    """Provide an isolated temporary directory in ~/tmp (never /tmp) for mutating test profiles."""
    home_tmp = Path.home() / "tmp"
    home_tmp.mkdir(parents=True, exist_ok=True)
    scratch = home_tmp / f"root_review_variant_test_{uuid.uuid4().hex[:10]}"
    scratch.mkdir(parents=True, exist_ok=True)
    try:
        yield scratch
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def test_actual_v60_profile_variant_contract_baseline():
    """Ground truth test on compact real V60 inStrain 1.9.1 profile fixture:
    Assert:
      - 24 summary variant sites in scaffold_info.tsv (24 SNV, 0 SNS, 24 divergent_site_count)
      - 22 exported non-cryptic SNVs in SNVs.tsv (all allele_count=2, class=SNV, cryptic=False)
      - 0 exported SNS
      - 2 excluded cryptic positions (zero-based 16646 and 28678)
      - raw cumulative reconciliation matches exactly and explains 24 summary vs 22 exported.
    """
    assert ACTUAL_V60_PROFILE.is_dir(), f"Actual profile fixture missing at {ACTUAL_V60_PROFILE}"

    res = qualify_instrain_profile_output(
        ACTUAL_V60_PROFILE,
        ACTUAL_V60_REPRESENTATIVE,
        min_cov=5,
    )

    # 24 summary variants in scaffold_info.tsv (unfiltered cumulative mm)
    assert res["snv_count_scaffold_info"] == 24
    assert res["snv_count_scaffold_info_unfiltered_cumulative"] == 24
    assert res["sns_count_scaffold_info"] == 0
    assert res["sns_count_scaffold_info_unfiltered_cumulative"] == 0
    assert res["divergent_site_count_scaffold_info"] == 24
    assert res["divergent_site_count_scaffold_info_unfiltered_cumulative"] == 24
    assert res["zero_allele_divergent_site_count_scaffold_info_unfiltered_cumulative"] == 0

    # 22 exported non-cryptic SNVs, 0 SNS
    assert res["exported_non_cryptic_snv_count"] == 22
    assert res["exported_non_cryptic_sns_count"] == 0
    assert res["exported_non_cryptic_zero_allele_divergent_site_count"] == 0
    assert res["snv_target_row_count"] == 22
    assert res["snv_target_non_cryptic_row_count"] == 22

    # 2 excluded cryptic positions
    assert res["cryptic_variant_site_count_excluded_from_exported_snv_table"] == 2
    excluded = res["excluded_cryptic_positions_zero_based"]
    assert len(excluded) == 2
    excluded_positions = {item["position"] for item in excluded}
    assert excluded_positions == {16646, 28678}
    for item in excluded:
        assert item["cryptic"] is True
        assert item["allele_count"] == 2

    # Overall verified status
    assert res["snv_table_status"] == "present_verified_target_rows_reconciled_with_raw_cumulative_cryptic_sites"

    # Raw cumulative reconciliation detail
    raw_rec = res["raw_cumulative_variant_reconciliation"]
    assert raw_rec is not None
    assert raw_rec["exported_noncryptic_rows_exactly_replayed"] == 22
    assert raw_rec["summary_at_max_mm"]["SNV_count"] == 24
    assert raw_rec["summary_at_max_mm"]["SNS_count"] == 0
    assert raw_rec["summary_at_max_mm"]["divergent_site_count"] == 24
    assert raw_rec["summary_at_max_mm"]["zero_allele_divergent_site_count"] == 0


def test_mutate_actual_export_var_base_requires_failure(isolated_scratch: Path):
    """Mutating an exported row's var_base while retaining row count and position must fail."""
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    snv_file = next(profile_dir.glob("output/*_SNVs.tsv"))
    lines = snv_file.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 23, "Expected header + 22 rows in actual SNVs.tsv"

    # Row 1 is position 10548 with var_base 'C'
    cols = lines[1].split("\t")
    original_var_base = cols[6]
    cols[6] = "T" if original_var_base != "T" else "A"
    lines[1] = "\t".join(cols)
    snv_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="mismatch at .* column 'var_base'"):
        qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)


def test_mutate_actual_export_ref_base_requires_failure(isolated_scratch: Path):
    """Mutating an exported row's ref_base while retaining row count and position must fail."""
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    snv_file = next(profile_dir.glob("output/*_SNVs.tsv"))
    lines = snv_file.read_text(encoding="utf-8").splitlines()
    cols = lines[1].split("\t")
    original_ref_base = cols[4]
    cols[4] = "A" if original_ref_base != "A" else "T"
    lines[1] = "\t".join(cols)
    snv_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="mismatch at .* column 'ref_base'"):
        qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)


def test_mutate_actual_export_var_frequency_requires_failure(isolated_scratch: Path):
    """Mutating an exported row's var_freq while retaining row count and position must fail."""
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    snv_file = next(profile_dir.glob("output/*_SNVs.tsv"))
    lines = snv_file.read_text(encoding="utf-8").splitlines()
    cols = lines[1].split("\t")
    # Mutate var_freq
    cols[9] = "0.2500000000000000"
    lines[1] = "\t".join(cols)
    snv_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="mismatch at .* column 'var_freq'"):
        qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)


def test_reject_malformed_cryptic_flags_in_exported_snv_table(isolated_scratch: Path):
    """SNVs.tsv must reject invalid cryptic values such as arbitrary strings or NaNs."""
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    snv_file = next(profile_dir.glob("output/*_SNVs.tsv"))
    lines = snv_file.read_text(encoding="utf-8").splitlines()
    cols = lines[1].split("\t")
    cols[14] = "maybe"  # cryptic column index
    lines[1] = "\t".join(cols)
    snv_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="invalid boolean value 'maybe' in column 'cryptic'"):
        qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)


def test_reject_malformed_cryptic_flags_in_raw_cumulative_table(isolated_scratch: Path):
    """cumulative_snv_table.csv.gz must reject invalid cryptic values rather than treating them as False."""
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    raw_path = profile_dir / "raw_data" / "cumulative_snv_table.csv.gz"
    raw_df = pd.read_csv(raw_path)
    raw_df["cryptic"] = raw_df["cryptic"].astype(object)
    raw_df.loc[0, "cryptic"] = "arbitrary_non_boolean"
    with gzip.open(raw_path, "wt", encoding="utf-8") as fh:
        raw_df.to_csv(fh, index=False)

    with pytest.raises(ValueError, match="invalid boolean value 'arbitrary_non_boolean' in column 'cryptic'"):
        qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)


def test_reject_fractional_allele_counts_in_exported_snv_table(isolated_scratch: Path):
    """Fractional allele counts in SNVs.tsv (e.g. 1.5) must be rejected rather than silently truncated."""
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    snv_file = next(profile_dir.glob("output/*_SNVs.tsv"))
    lines = snv_file.read_text(encoding="utf-8").splitlines()
    cols = lines[1].split("\t")
    cols[3] = "1.5"  # allele_count column index
    lines[1] = "\t".join(cols)
    snv_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="non-integer float|allele_count"):
        qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)


def test_reject_fractional_allele_counts_in_raw_cumulative_table(isolated_scratch: Path):
    """Fractional allele counts in raw cumulative table must be rejected rather than silently truncated."""
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    raw_path = profile_dir / "raw_data" / "cumulative_snv_table.csv.gz"
    raw_df = pd.read_csv(raw_path)
    raw_df["allele_count"] = raw_df["allele_count"].astype(object)
    raw_df.loc[0, "allele_count"] = 2.5
    with gzip.open(raw_path, "wt", encoding="utf-8") as fh:
        raw_df.to_csv(fh, index=False)

    with pytest.raises(ValueError, match="non-integer float|allele_count"):
        qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)


def test_cryptic_only_fixture_all_raw_flags_true_exported_newline(isolated_scratch: Path):
    """A profile where all raw cumulative records are flagged cryptic=True:
    inStrain SNVprofile.get_nonredundant_snv_table(cryptic=False) produces an empty DataFrame,
    exported to disk as a 1-byte newline ('\\n').
    The qualification must report:
      - 0 exported non-cryptic SNVs
      - 0 exported non-cryptic SNS
      - 0 exported rows
      - 24 summary SNVs and 24 divergent sites in scaffold_info.tsv
      - 24 cryptic sites excluded from exported table
      - status 'present_empty_dataframe_newline_cryptic_only_verified_with_raw_cumulative_table'
    """
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    # Flag all rows in cumulative_snv_table as cryptic=True
    raw_path = profile_dir / "raw_data" / "cumulative_snv_table.csv.gz"
    raw_df = pd.read_csv(raw_path)
    raw_df["cryptic"] = True
    with gzip.open(raw_path, "wt", encoding="utf-8") as fh:
        raw_df.to_csv(fh, index=False)

    # Replace SNVs.tsv with 1-byte newline
    snv_file = next(profile_dir.glob("output/*_SNVs.tsv"))
    snv_file.write_text("\n", encoding="utf-8")

    res = qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)

    assert res["snv_count_scaffold_info"] == 24
    assert res["divergent_site_count_scaffold_info"] == 24
    assert res["exported_non_cryptic_snv_count"] == 0
    assert res["exported_non_cryptic_sns_count"] == 0
    assert res["snv_target_row_count"] == 0
    assert res["snv_target_non_cryptic_row_count"] == 0
    assert res["cryptic_variant_site_count_excluded_from_exported_snv_table"] == 24
    assert res["snv_table_status"] == "present_empty_dataframe_newline_cryptic_only_verified_with_raw_cumulative_table"


def test_zero_allele_divergent_site_retains_separate_count(isolated_scratch: Path):
    """Installed call_snv_site permits morphia/allele_count==0 with class 'DivergentSite'.
    calc_snps counts len(db) as divergent_site_count, while SNS is allele_count==1 and SNV is >1.
    Therefore divergent_site_count = SNV_count + SNS_count + zero_allele_divergent_site_count.
    The qualification must NOT assume SNS+SNV==divergent, but retain the separate zero-allele count.
    """
    profile_dir = isolated_scratch / "profile"
    profile_dir.mkdir(parents=True)
    (profile_dir / "output").mkdir()
    (profile_dir / "raw_data").mkdir()
    (profile_dir / "log").mkdir()
    (profile_dir / "log" / "log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")

    rep = "VIRUS_ZERO_ALLELE_SYNTH_REF"

    # scaffold_info: 1 SNS, 1 SNV, but 3 divergent sites
    scaff_info_content = (
        "scaffold\tlength\tcoverage\tbreadth\tbreadth_minCov\tSNS_count\tSNV_count\tdivergent_site_count\n"
        f"{rep}\t1000\t15.0\t1.0\t1.0\t1\t1\t3\n"
    )
    (profile_dir / "output" / "synth_scaffold_info.tsv").write_text(scaff_info_content, encoding="utf-8")

    # cumulative_scaffold_table at mm=1
    scaff_cum = pd.DataFrame([{
        "scaffold": rep,
        "mm": 1,
        "SNS_count": 1,
        "SNV_count": 1,
        "divergent_site_count": 3,
    }])
    with gzip.open(profile_dir / "raw_data" / "cumulative_scaffold_table.csv.gz", "wt", encoding="utf-8") as fh:
        scaff_cum.to_csv(fh, index=False)

    # 3 sites: 1 SNS (ac=1), 1 SNV (ac=2), 1 DivergentSite (ac=0)
    cum_snv = pd.DataFrame([
        {
            "scaffold": rep, "position": 100, "ref_base": "A", "con_base": "T", "var_base": "T",
            "A": 0, "C": 0, "T": 15, "G": 0, "mm": 1, "allele_count": 1, "class": "SNS", "cryptic": False,
            "position_coverage": 15, "var_freq": 1.0, "con_freq": 1.0, "ref_freq": 0.0,
        },
        {
            "scaffold": rep, "position": 200, "ref_base": "C", "con_base": "C", "var_base": "G",
            "A": 0, "C": 12, "T": 0, "G": 3, "mm": 1, "allele_count": 2, "class": "SNV", "cryptic": False,
            "position_coverage": 15, "var_freq": 0.2, "con_freq": 0.8, "ref_freq": 0.8,
        },
        {
            "scaffold": rep, "position": 300, "ref_base": "G", "con_base": "G", "var_base": "G",
            "A": 0, "C": 0, "T": 0, "G": 0, "mm": 1, "allele_count": 0, "class": "DivergentSite", "cryptic": False,
            "position_coverage": 0, "var_freq": 0.0, "con_freq": 0.0, "ref_freq": 0.0,
        },
    ])
    with gzip.open(profile_dir / "raw_data" / "cumulative_snv_table.csv.gz", "wt", encoding="utf-8") as fh:
        cum_snv.to_csv(fh, index=False)

    # SNVs.tsv contains the 3 exported non-cryptic sites
    cum_snv.to_csv(profile_dir / "output" / "synth_SNVs.tsv", sep="\t", index=False)

    res = qualify_instrain_profile_output(profile_dir, rep, min_cov=5)

    # Divergent sites count (3) exceeds SNS(1) + SNV(1) = 2, because of 1 zero-allele divergent site
    assert res["snv_count_scaffold_info"] == 1
    assert res["sns_count_scaffold_info"] == 1
    assert res["divergent_site_count_scaffold_info"] == 3
    assert res["zero_allele_divergent_site_count_scaffold_info_unfiltered_cumulative"] == 1
    assert res["exported_non_cryptic_snv_count"] == 1
    assert res["exported_non_cryptic_sns_count"] == 1
    assert res["exported_non_cryptic_zero_allele_divergent_site_count"] == 1
    assert res["snv_target_row_count"] == 3


def test_missing_allele_field_in_raw_cumulative_table_fails_closed(isolated_scratch: Path):
    """Missing allele_count column in cumulative_snv_table.csv.gz must fail closed."""
    profile_dir = isolated_scratch / "profile"
    shutil.copytree(ACTUAL_V60_PROFILE, profile_dir)

    raw_path = profile_dir / "raw_data" / "cumulative_snv_table.csv.gz"
    raw_df = pd.read_csv(raw_path)
    raw_df = raw_df.drop(columns=["allele_count"])
    with gzip.open(raw_path, "wt", encoding="utf-8") as fh:
        raw_df.to_csv(fh, index=False)

    with pytest.raises(ValueError, match="missing required columns.*allele_count"):
        qualify_instrain_profile_output(profile_dir, ACTUAL_V60_REPRESENTATIVE, min_cov=5)


def test_missing_allele_field_in_snvs_tsv_must_fail_or_return_explicit_unknown(isolated_scratch: Path):
    """When SNVs.tsv lacks allele_count column (and no raw table exists to supply it),
    the parser must either fail or report explicit unknown; it must NEVER invent SNS=0.
    """
    profile_dir = isolated_scratch / "profile"
    profile_dir.mkdir(parents=True)
    (profile_dir / "output").mkdir()
    (profile_dir / "log").mkdir()
    (profile_dir / "log" / "log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")

    rep = "VIRUS_NO_AC_SYNTH_REF"
    (profile_dir / "output" / "synth_scaffold_info.tsv").write_text(
        "scaffold\tlength\tcoverage\tbreadth\tbreadth_minCov\tSNS_count\tSNV_count\tdivergent_site_count\n"
        f"{rep}\t1000\t10.0\t1.0\t1.0\t0\t2\t2\n",
        encoding="utf-8",
    )
    # SNVs.tsv with required columns but NO allele_count
    (profile_dir / "output" / "synth_SNVs.tsv").write_text(
        "scaffold\tposition\tref_base\tcon_base\n"
        f"{rep}\t10\tA\tT\n"
        f"{rep}\t20\tC\tG\n",
        encoding="utf-8",
    )

    try:
        res = qualify_instrain_profile_output(profile_dir, rep, min_cov=5)
    except ValueError:
        # Failing closed is an acceptable outcome
        return

    # If it returns a result, exported SNS count must NOT be numeric 0
    exported_sns = res.get("exported_non_cryptic_sns_count")
    assert exported_sns != 0, f"Must never invent SNS=0 when allele_count is absent; got {exported_sns!r}"
    assert "unknown" in str(exported_sns).lower(), f"Expected explicit unknown, got {exported_sns!r}"


def test_missing_sns_count_in_scaffold_info_must_fail_or_return_explicit_unknown(isolated_scratch: Path):
    """When scaffold_info.tsv lacks SNS_count, the parser must report explicit unknown or fail,
    never invent 0 SNS.
    """
    profile_dir = isolated_scratch / "profile"
    profile_dir.mkdir(parents=True)
    (profile_dir / "output").mkdir()
    (profile_dir / "log").mkdir()
    (profile_dir / "log" / "log.log").write_text(f"start\n{INSTRAIN_PROFILE_FINISHED_BANNER}\n", encoding="utf-8")

    rep = "VIRUS_NO_SNS_SYNTH_REF"
    # Lacks SNS_count column
    (profile_dir / "output" / "synth_scaffold_info.tsv").write_text(
        "scaffold\tlength\tcoverage\tbreadth\tbreadth_minCov\tSNV_count\n"
        f"{rep}\t1000\t10.0\t1.0\t1.0\t2\n",
        encoding="utf-8",
    )
    (profile_dir / "output" / "synth_SNVs.tsv").write_text(
        "scaffold\tposition\tref_base\tcon_base\n"
        f"{rep}\t10\tA\tT\n"
        f"{rep}\t20\tC\tG\n",
        encoding="utf-8",
    )

    try:
        res = qualify_instrain_profile_output(profile_dir, rep, min_cov=5)
    except ValueError:
        # Failing closed is an acceptable outcome
        return

    sns_val = res.get("sns_count_scaffold_info")
    assert sns_val != 0, f"Must never invent SNS=0 when SNS_count column is missing; got {sns_val!r}"
    assert "unknown" in str(sns_val).lower(), f"Expected explicit unknown, got {sns_val!r}"
