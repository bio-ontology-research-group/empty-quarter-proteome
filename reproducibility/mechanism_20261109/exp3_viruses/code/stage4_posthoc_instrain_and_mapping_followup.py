#!/usr/bin/env python3
"""Posthoc inStrain callable-position semantics, cross-sample ANI>0.92 sensitivity,
population comparison eligibility, and completed mapping + stats_v4 reconciliation.

Follows:
- mechanism-20261109/EXP3-AFTER-MAPPING-COMPLETION-REVIEW.txt
- mechanism-20261109/EXP3-INSTRAIN-CALLABLE-SEMANTICS-REVIEW.txt
- mechanism-20261109/EXP3-POSTHOC-COMPLETION-CONTRACT-REVIEW.txt

Never modifies running scripts (stage4_competitive_common.py, stage4_map_sample_competitive.py,
stage4_aggregate_competitive_mapping.py) or rebuilds the immutable competitive reference/index.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import math
import os
import shlex
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any

import pandas as pd

from stage4_aggregate_competitive_mapping import one_row, verify_sample_receipt
from stage4_competitive_common import (
    SHA256_HEX,
    atomic_json,
    atomic_text,
    load_analysis_pairs,
    require_columns,
    sha256_sequence,
    reverse_complement,
    stream_sha256,
)
from stage4_host_array_provenance_audit import run_host_array_provenance_audit
from stage4_map_sample_competitive import SIF_INSTRAIN

BRANCH_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))

PLANNED_14_SAMPLES: tuple[str, ...] = (
    "41PRr2",
    "V60PRr2",
    "55Dr2",
    "9Sr2",
    "55Sr3",
    "S49Dr1",
    "48Sr2",
    "F48Sr2",
    "V48Dr1",
    "S48PRr1",
    "41Sr3",
    "10Dr2",
    "43PRr3",
    "44Dr2",
)

INSTRAIN_PROFILE_FINISHED_BANNER = "..:: inStrain profile finished ::.."
INSTRAIN_COMPARE_FINISHED_BANNER = "..:: inStrain compare finished ::.."

INSTRAIN_191_ORIGINAL_SNVPROFILE_SHA256 = "a2ec6c7d3671185210e6fcac987fb2bec471b6fc6a83765e1ebe2233f129c55e"
INSTRAIN_191_PATCHED_SNVPROFILE_SHA256 = "f65ac08f1997f6a15b38902aee44afba75de7a6384b3539d66c362ec94d74bcd"
INSTRAIN_191_PATCH_DIFF_SHA256 = "963cc11e71b505b0cf024bed5311df829ae6a97f14eec8bb02a22f60a23c009a"
INSTRAIN_CONTAINER_SNVPROFILE_TARGET = "/usr/local/lib/python3.8/site-packages/inStrain/SNVprofile.py"


def verify_and_resolve_instrain_namespaced_header_overlay(
    root: Path | None = None,
    instrain_sif: Path | None = None,
    singularity_exec: str | None = None,
    verify_container_runtime: bool = False,
) -> dict[str, Any]:
    """Verify the version/hash-bound local SNVprofile.py overlay patch that fixes
    inStrain 1.9.1's `thing.split('::')` bug on namespaced `::` reference headers
    by replacing the three delimiter-unpacking lines with `rsplit('::', 1)` / `rsplit('::', 2)`.
    """
    overlay_dir: Path | None = None
    if root is not None and (root / "code/instrain_191_overlay/instrain_191_namespaced_header_overlay_manifest.json").exists():
        overlay_dir = root / "code/instrain_191_overlay"
    else:
        fallback_dir = Path(__file__).resolve().parent / "instrain_191_overlay"
        if (fallback_dir / "instrain_191_namespaced_header_overlay_manifest.json").exists():
            overlay_dir = fallback_dir
    if overlay_dir is None:
        raise FileNotFoundError("inStrain 1.9.1 namespaced-header overlay directory/manifest not found")

    manifest_path = overlay_dir / "instrain_191_namespaced_header_overlay_manifest.json"
    orig_path = overlay_dir / "SNVprofile.original.py"
    patched_path = overlay_dir / "SNVprofile.py"
    diff_path = overlay_dir / "SNVprofile_namespaced_header_rsplit.patch"
    for p in (manifest_path, orig_path, patched_path, diff_path):
        if not p.exists() or p.stat().st_size == 0:
            raise FileNotFoundError(f"inStrain 1.9.1 overlay file missing/empty: {p}")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    orig_sha = stream_sha256(orig_path)
    patched_sha = stream_sha256(patched_path)
    diff_sha = stream_sha256(diff_path)
    if (
        orig_sha != INSTRAIN_191_ORIGINAL_SNVPROFILE_SHA256
        or manifest.get("original_snvprofile_py_sha256") != INSTRAIN_191_ORIGINAL_SNVPROFILE_SHA256
    ):
        raise ValueError(f"original SNVprofile.py SHA-256 mismatch: {orig_sha}")
    if (
        patched_sha != INSTRAIN_191_PATCHED_SNVPROFILE_SHA256
        or manifest.get("patched_snvprofile_py_sha256") != INSTRAIN_191_PATCHED_SNVPROFILE_SHA256
    ):
        raise ValueError(f"patched SNVprofile.py SHA-256 mismatch: {patched_sha}")
    if (
        diff_sha != INSTRAIN_191_PATCH_DIFF_SHA256
        or manifest.get("patch_diff_sha256") != INSTRAIN_191_PATCH_DIFF_SHA256
    ):
        raise ValueError(f"SNVprofile patch diff SHA-256 mismatch: {diff_sha}")

    bind_spec = f"{patched_path.resolve()}:{INSTRAIN_CONTAINER_SNVPROFILE_TARGET}:ro"
    container_verified = False
    if verify_container_runtime and instrain_sif is not None and instrain_sif.exists() and singularity_exec:
        probe_unmounted = subprocess.run(
            [
                singularity_exec,
                "exec",
                str(instrain_sif),
                "python3",
                "-c",
                f"import hashlib, pathlib; print(hashlib.sha256(pathlib.Path({INSTRAIN_CONTAINER_SNVPROFILE_TARGET!r}).read_bytes()).hexdigest())",
            ],
            capture_output=True,
            text=True,
        )
        if probe_unmounted.returncode != 0 or probe_unmounted.stdout.strip() != INSTRAIN_191_ORIGINAL_SNVPROFILE_SHA256:
            raise ValueError(
                f"unmounted container {instrain_sif} SNVprofile.py hash mismatch: "
                f"rc={probe_unmounted.returncode}, stdout={probe_unmounted.stdout.strip()!r}, stderr={probe_unmounted.stderr.strip()!r}"
            )
        probe_mounted = subprocess.run(
            [
                singularity_exec,
                "exec",
                "-B",
                bind_spec,
                str(instrain_sif),
                "python3",
                "-c",
                f"import hashlib, pathlib; print(hashlib.sha256(pathlib.Path({INSTRAIN_CONTAINER_SNVPROFILE_TARGET!r}).read_bytes()).hexdigest())",
            ],
            capture_output=True,
            text=True,
        )
        if probe_mounted.returncode != 0 or probe_mounted.stdout.strip() != INSTRAIN_191_PATCHED_SNVPROFILE_SHA256:
            raise ValueError(
                f"bind-mounted container {instrain_sif} SNVprofile.py hash mismatch: "
                f"rc={probe_mounted.returncode}, stdout={probe_mounted.stdout.strip()!r}, stderr={probe_mounted.stderr.strip()!r}"
            )
        container_verified = True

    return {
        "tool_version": str(manifest["tool_version"]),
        "container_target_file": INSTRAIN_CONTAINER_SNVPROFILE_TARGET,
        "original_snvprofile_py_sha256": orig_sha,
        "patched_snvprofile_py_sha256": patched_sha,
        "patch_diff_sha256": diff_sha,
        "overlay_manifest_path": str(manifest_path),
        "overlay_manifest_sha256": stream_sha256(manifest_path),
        "patched_snvprofile_py_path": str(patched_path.resolve()),
        "singularity_bind_arg": bind_spec,
        "container_runtime_verified": container_verified,
    }

REQUIRED_SNV_COLUMNS = {"scaffold", "position", "ref_base", "con_base"}
REQUIRED_RAW_CUMULATIVE_SNV_COLUMNS = {
    "scaffold",
    "position",
    "ref_base",
    "A",
    "C",
    "T",
    "G",
    "con_base",
    "var_base",
    "mm",
    "allele_count",
    "class",
    "cryptic",
    "position_coverage",
    "var_freq",
    "con_freq",
    "ref_freq",
}
REQUIRED_COMPARE_COLUMNS = {
    "scaffold",
    "name1",
    "name2",
    "coverage_overlap",
    "compared_bases_count",
    "percent_genome_compared",
    "length",
    "consensus_SNPs",
    "population_SNPs",
    "popANI",
    "conANI",
}


def load_planned_samples(root: Path, pairs: pd.DataFrame) -> list[str]:
    """Load the exact 14 planned competitive-mapping sample identities from the pipeline
    submission manifest (and sbatch_stage4_map_array.sh when present) and verify that every
    direct FASTQ sample in the analysis pairs + predetermined control 10Dr2 is covered.
    """
    derived_samples = set(pairs["host_sample"].astype(str))
    for vs_list in pairs["virus_samples"]:
        if isinstance(vs_list, (list, tuple, set)):
            derived_samples.update(str(s) for s in vs_list)
        elif str(vs_list).strip():
            derived_samples.update(s.strip() for s in str(vs_list).split(";") if s.strip())
    derived_samples.add("10Dr2")
    direct_fastq_samples = {s for s in derived_samples if not s.startswith("coasm_")}

    sub_json = root / "provenance/stage4_competitive_pipeline_submission.json"
    if sub_json.exists():
        payload = json.loads(sub_json.read_text(encoding="utf-8"))
        manifest_samples = [str(s) for s in payload.get("samples", [])]
    else:
        manifest_samples = list(PLANNED_14_SAMPLES)

    if len(manifest_samples) != 14 or len(set(manifest_samples)) != 14:
        raise ValueError(f"submission manifest does not list 14 unique planned samples: {manifest_samples}")
    if set(manifest_samples) != set(PLANNED_14_SAMPLES):
        raise ValueError(
            f"planned sample mismatch between submission manifest ({manifest_samples}) "
            f"and canonical 14-sample tuple ({list(PLANNED_14_SAMPLES)})"
        )
    if not direct_fastq_samples <= set(manifest_samples):
        raise ValueError(
            f"direct FASTQ samples from analysis pairs {sorted(direct_fastq_samples)} "
            f"not covered by planned 14-sample manifest {manifest_samples}"
        )
    return manifest_samples


def _require_instrain_completion_banner(output_dir: Path, expected_banner: str, context_label: str) -> Path:
    """Require the genuine inStrain 1.9.1 terminal completion banner in log/log.log.
    Generic strings ('Finished', 'Done', 'inStrain profile', 'scaffold_info') appear before
    completion and are never sufficient to certify a completed inStrain run.
    """
    log_candidates = sorted(output_dir.glob("**/log/log.log"))
    if not log_candidates:
        log_candidates = sorted(output_dir.glob("**/log.log"))
    if len(log_candidates) != 1:
        raise ValueError(f"{context_label}: expected exactly one log/log.log in {output_dir}, found={log_candidates}")
    log_path = log_candidates[0]
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    if expected_banner not in log_text:
        raise ValueError(
            f"{context_label}: inStrain log {log_path} lacks terminal completion banner {expected_banner!r}; "
            "startup/intermediate log lines do not establish completion"
        )
    return log_path


def _extract_min_cov_provenance(commands: list[dict[str, Any]], target_sha: str, default_min_cov: int = 5) -> tuple[int, str]:
    """Determine the active min_cov setting from recorded inStrain profile command arguments
    or the verified inStrain 1.9.1 argumentParser default (-c/--min_cov default=5).
    """
    for entry in commands:
        cmd = entry.get("command", [])
        if not isinstance(cmd, list):
            continue
        cmd_str = " ".join(str(x) for x in cmd)
        if "inStrain" in cmd_str and "profile" in cmd_str and target_sha[:16] in cmd_str:
            for idx, tok in enumerate(cmd):
                if str(tok) in ("-c", "--min_cov") and idx + 1 < len(cmd):
                    return int(cmd[idx + 1]), f"explicit_cli_flag_min_cov_{int(cmd[idx + 1])}"
            return (
                int(default_min_cov),
                f"instrain_1.9.1_verified_default_min_cov_{int(default_min_cov)}_from_argumentParser",
            )
    return (
        int(default_min_cov),
        f"instrain_1.9.1_verified_default_min_cov_{int(default_min_cov)}_from_argumentParser",
    )


def _parse_strict_bool(val: Any, col_name: str, context: str) -> bool:
    """Strictly parse a boolean value without silently converting arbitrary strings, numbers, or NaN to False."""
    if pd.isna(val):
        raise ValueError(f"{context}: missing/NaN value in boolean column {col_name!r}")
    if isinstance(val, bool):
        return val
    if isinstance(val, (int, float)):
        raise ValueError(f"{context}: numeric value {val!r} not permitted in boolean column {col_name!r}")
    s = str(val).strip()
    if s in {"True", "true"}:
        return True
    if s in {"False", "false"}:
        return False
    raise ValueError(f"{context}: invalid boolean value {val!r} in column {col_name!r}")


def _parse_strict_bool_series(series: pd.Series, col_name: str, context: str) -> pd.Series:
    return pd.Series(
        [_parse_strict_bool(v, col_name, context) for v in series.tolist()],
        index=series.index,
        dtype=bool,
    )


def _parse_strict_int(val: Any, col_name: str, context: str) -> int:
    """Strictly parse an integer value without silently truncating floats or accepting booleans/NaN."""
    if pd.isna(val) or isinstance(val, bool):
        raise ValueError(f"{context}: missing or boolean value {val!r} in integer column {col_name!r}")
    if isinstance(val, int):
        return int(val)
    if isinstance(val, float):
        if not math.isfinite(val) or not val.is_integer():
            raise ValueError(f"{context}: non-integer float {val!r} in integer column {col_name!r}")
        return int(val)
    s = str(val).strip()
    try:
        return int(s)
    except ValueError as exc:
        raise ValueError(f"{context}: non-integer value {val!r} in integer column {col_name!r}") from exc


def _parse_strict_int_series(series: pd.Series, col_name: str, context: str) -> pd.Series:
    return pd.Series(
        [_parse_strict_int(v, col_name, context) for v in series.tolist()],
        index=series.index,
        dtype="int64",
    )


def _parse_strict_float(val: Any, col_name: str, context: str) -> float:
    if pd.isna(val) or isinstance(val, bool):
        raise ValueError(f"{context}: missing or boolean value {val!r} in float column {col_name!r}")
    try:
        fval = float(val)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"{context}: non-numeric value {val!r} in float column {col_name!r}") from exc
    if not math.isfinite(fval):
        raise ValueError(f"{context}: non-finite float {val!r} in float column {col_name!r}")
    return fval


def _verify_exported_snv_rows_match_raw_replayed(
    exported_df: pd.DataFrame,
    replayed_rows: list[dict[str, Any]],
    output_dir: Path,
    representative: str,
) -> None:
    """Verify that exported `<prefix>_SNVs.tsv` target rows match the replayed raw cumulative
    non-cryptic rows on all common actual export fields keyed by `(scaffold, position)`:
    `ref_base`, `con_base`, `var_base`, `allele_count`, `class`, `position_coverage`,
    `A`, `C`, `T`, `G`, `cryptic`, and frequencies `ref_freq`, `con_freq`, `var_freq`.
    """
    ctx = f"SNVs.tsv ({output_dir}, {representative})"
    if len(exported_df) != len(replayed_rows):
        raise ValueError(
            f"{ctx}: exported row count ({len(exported_df)}) != replayed raw non-cryptic row count ({len(replayed_rows)})"
        )
    exp_positions = [_parse_strict_int(v, "position", ctx) for v in exported_df["position"].tolist()]
    raw_positions = [int(r["position"]) for r in replayed_rows]
    if exp_positions != raw_positions:
        raise ValueError(f"{ctx}: exported positions {exp_positions} != replayed raw positions {raw_positions}")

    replayed_by_key = {(str(r["scaffold"]), int(r["position"])): r for r in replayed_rows}
    for erow in exported_df.to_dict(orient="records"):
        key = (str(erow["scaffold"]), _parse_strict_int(erow["position"], "position", ctx))
        if key not in replayed_by_key:
            raise ValueError(f"{ctx}: exported key {key} absent from replayed raw non-cryptic rows")
        rrow = replayed_by_key[key]

        for str_col in ("ref_base", "con_base", "var_base", "class"):
            if str_col in exported_df.columns:
                if pd.isna(erow[str_col]) or str(erow[str_col]).strip() != str(rrow[str_col]).strip():
                    raise ValueError(
                        f"{ctx}: mismatch at {key} column {str_col!r}: exported={erow[str_col]!r} != raw_replayed={rrow[str_col]!r}"
                    )
        for int_col in ("position_coverage", "allele_count", "A", "C", "T", "G"):
            if int_col in exported_df.columns:
                exp_int = _parse_strict_int(erow[int_col], int_col, ctx)
                raw_int = int(rrow[int_col])
                if exp_int != raw_int:
                    raise ValueError(
                        f"{ctx}: mismatch at {key} column {int_col!r}: exported={exp_int} != raw_replayed={raw_int}"
                    )
        if "cryptic" in exported_df.columns:
            exp_cryptic = _parse_strict_bool(erow["cryptic"], "cryptic", ctx)
            if exp_cryptic is not False or bool(rrow["cryptic"]) is not False:
                raise ValueError(
                    f"{ctx}: mismatch at {key} column 'cryptic': exported={exp_cryptic} vs raw_replayed={rrow['cryptic']}"
                )
        for float_col in ("ref_freq", "con_freq", "var_freq"):
            if float_col in exported_df.columns:
                exp_flt = _parse_strict_float(erow[float_col], float_col, ctx)
                raw_flt = float(rrow[float_col])
                if not math.isclose(exp_flt, raw_flt, rel_tol=1e-9, abs_tol=1e-12):
                    raise ValueError(
                        f"{ctx}: mismatch at {key} column {float_col!r}: exported={exp_flt} != raw_replayed={raw_flt}"
                    )


def reconcile_instrain_raw_cumulative_variants(
    output_dir: Path,
    representative: str,
) -> dict[str, Any] | None:
    """Replay installed inStrain 1.9.1 variant summary and exported SNV table algorithms from
    `raw_data/cumulative_snv_table.csv.gz` and `raw_data/cumulative_scaffold_table.csv.gz` when
    present in `output_dir`:
    - `inStrain.profile.snv_utilities.calc_snps(Odb, mm)` (lines 249-272) deduplicates
      `Odb[Odb['mm'] <= mm].sort_values('mm').drop_duplicates(subset=['position'], keep='last')`
      WITHOUT filtering `cryptic == True`, so `scaffold_info.tsv` `SNV_count` (`allele_count > 1`),
      `SNS_count` (`allele_count == 1`), and `divergent_site_count` (`len(db)`, which also includes any
      `allele_count == 0` `DivergentSite` records) count all cumulative variant sites up to `max_mm`
      including cryptic sites.
    - `SNVprofile.get_nonredundant_snv_table(cryptic=False)` (lines 495-510) filters
      `scdb[scdb['cryptic'] == False]` BEFORE `.sort_values('mm').drop_duplicates(subset=['scaffold', 'position'], keep='last')`,
      so exported `<prefix>_SNVs.tsv` excludes `cryptic == True` positions.
    """
    raw_snv_path = output_dir / "raw_data/cumulative_snv_table.csv.gz"
    raw_scaff_path = output_dir / "raw_data/cumulative_scaffold_table.csv.gz"
    if not raw_snv_path.exists() and not raw_scaff_path.exists():
        return None
    if raw_scaff_path.exists() and not raw_snv_path.exists():
        raise ValueError(
            f"raw_data/cumulative_scaffold_table.csv.gz exists without raw_data/cumulative_snv_table.csv.gz in {output_dir}"
        )

    ctx_cum = f"raw_data/cumulative_snv_table.csv.gz ({output_dir})"
    with gzip.open(raw_snv_path, "rt", encoding="utf-8") as fh:
        raw_csv_text = fh.read()

    if raw_csv_text.strip() in {"", '""'}:
        raw_snv_df = pd.DataFrame()
        target_cum = pd.DataFrame()
    else:
        raw_snv_df = pd.read_csv(io.StringIO(raw_csv_text))
        missing_cum_cols = sorted(REQUIRED_RAW_CUMULATIVE_SNV_COLUMNS - set(raw_snv_df.columns))
        if missing_cum_cols:
            raise ValueError(
                f"{ctx_cum} missing required columns {missing_cum_cols}; actual={list(raw_snv_df.columns)}"
            )
        target_cum = raw_snv_df[raw_snv_df["scaffold"].astype(str) == representative].copy()
        if not target_cum.empty:
            for int_col in ("position", "mm", "allele_count", "position_coverage", "A", "C", "T", "G"):
                target_cum[int_col] = _parse_strict_int_series(target_cum[int_col], int_col, ctx_cum)
            for flt_col in ("var_freq", "con_freq", "ref_freq"):
                target_cum[flt_col] = [_parse_strict_float(v, flt_col, ctx_cum) for v in target_cum[flt_col].tolist()]
            target_cum["_cryptic_bool"] = _parse_strict_bool_series(target_cum["cryptic"], "cryptic", ctx_cum)

            if target_cum.duplicated(subset=["scaffold", "position", "mm"]).any():
                raise ValueError(
                    f"malformed cumulative_snv_table in {output_dir}: duplicate (scaffold, position, mm) rows for {representative}"
                )
            for pos_val, pos_grp in target_cum.groupby("position"):
                if pos_grp["_cryptic_bool"].nunique() > 1:
                    raise ValueError(
                        f"malformed cumulative_snv_table in {output_dir}: position {pos_val} has inconsistent cryptic flags across mm levels"
                    )
            for _, r_row in target_cum.iterrows():
                ac = int(r_row["allele_count"])
                cls = str(r_row["class"]).strip()
                if ac < 0:
                    raise ValueError(f"malformed cumulative_snv_table row: negative allele_count={ac}")
                if ac == 0 and cls != "DivergentSite":
                    raise ValueError(f"malformed cumulative_snv_table row: allele_count=0 with class={cls!r}")
                if ac == 1 and (cls != "SNS" or str(r_row["con_base"]).strip() == str(r_row["ref_base"]).strip()):
                    raise ValueError(
                        f"malformed cumulative_snv_table row: allele_count=1 with class={cls!r}, "
                        f"ref_base={r_row['ref_base']!r}, con_base={r_row['con_base']!r}"
                    )
                if ac > 1 and cls not in {"SNV", "con_SNV", "pop_SNV"}:
                    raise ValueError(f"malformed cumulative_snv_table row: allele_count={ac} with class={cls!r}")

    summary_counts_all_mm: list[dict[str, int]] = []
    if raw_scaff_path.exists():
        ctx_scaff = f"raw_data/cumulative_scaffold_table.csv.gz ({output_dir})"
        raw_scaff_df = pd.read_csv(raw_scaff_path)
        require_columns(
            raw_scaff_df,
            ["scaffold", "mm", "SNS_count", "SNV_count", "divergent_site_count"],
            ctx_scaff,
        )
        target_scaff_cum = raw_scaff_df[raw_scaff_df["scaffold"].astype(str) == representative].copy()
        if target_scaff_cum.empty:
            raise ValueError(f"{ctx_scaff} lacks rows for {representative}")
        for int_col in ("mm", "SNS_count", "SNV_count", "divergent_site_count"):
            target_scaff_cum[int_col] = _parse_strict_int_series(target_scaff_cum[int_col], int_col, ctx_scaff)
        target_scaff_cum = target_scaff_cum.sort_values("mm", kind="mergesort")
        for _, srow in target_scaff_cum.iterrows():
            mm_val = int(srow["mm"])
            if target_cum.empty:
                replayed_sns = 0
                replayed_snv = 0
                replayed_ac0 = 0
                replayed_div = 0
            else:
                sub_mm = (
                    target_cum[target_cum["mm"] <= mm_val]
                    .sort_values("mm", kind="mergesort")
                    .drop_duplicates(subset=["position"], keep="last")
                )
                replayed_sns = int((sub_mm["allele_count"] == 1).sum())
                replayed_snv = int((sub_mm["allele_count"] > 1).sum())
                replayed_ac0 = int((sub_mm["allele_count"] == 0).sum())
                replayed_div = int(len(sub_mm))
            obs_sns = int(srow["SNS_count"])
            obs_snv = int(srow["SNV_count"])
            obs_div = int(srow["divergent_site_count"])
            if (replayed_sns, replayed_snv, replayed_div) != (obs_sns, obs_snv, obs_div):
                raise ValueError(
                    f"cumulative_scaffold_table mismatch at mm={mm_val} for {representative}: "
                    f"observed (SNS={obs_sns}, SNV={obs_snv}, div={obs_div}) != "
                    f"replayed from cumulative_snv_table (SNS={replayed_sns}, SNV={replayed_snv}, div={replayed_div})"
                )
            summary_counts_all_mm.append(
                {
                    "mm": mm_val,
                    "SNS_count": replayed_sns,
                    "SNV_count": replayed_snv,
                    "zero_allele_divergent_site_count": replayed_ac0,
                    "divergent_site_count": replayed_div,
                }
            )
        max_mm = int(target_scaff_cum["mm"].max())
        if not target_cum.empty and int(target_cum["mm"].max()) > max_mm:
            raise ValueError(
                f"cumulative_snv_table max mm ({int(target_cum['mm'].max())}) exceeds "
                f"cumulative_scaffold_table max mm ({max_mm}) for {representative}"
            )
    else:
        max_mm = int(target_cum["mm"].max()) if not target_cum.empty else 0

    if target_cum.empty:
        summary_sns_unfiltered = 0
        summary_snv_unfiltered = 0
        summary_ac0_unfiltered = 0
        summary_divergent_unfiltered = 0
        exported_non_cryptic_rows = 0
        exported_non_cryptic_snv_count = 0
        exported_non_cryptic_sns_count = 0
        exported_non_cryptic_ac0_count = 0
        exported_positions: list[int] = []
        exported_allele_counts: list[int] = []
        exported_rows_replayed: list[dict[str, Any]] = []
        excluded_cryptic_positions: list[dict[str, Any]] = []
    else:
        dedup_all_mm = (
            target_cum[target_cum["mm"] <= max_mm]
            .sort_values("mm", kind="mergesort")
            .drop_duplicates(subset=["position"], keep="last")
        )
        summary_sns_unfiltered = int((dedup_all_mm["allele_count"] == 1).sum())
        summary_snv_unfiltered = int((dedup_all_mm["allele_count"] > 1).sum())
        summary_ac0_unfiltered = int((dedup_all_mm["allele_count"] == 0).sum())
        summary_divergent_unfiltered = int(len(dedup_all_mm))

        non_cryptic_cum = target_cum[~target_cum["_cryptic_bool"]]
        if non_cryptic_cum.empty:
            dedup_non_cryptic = pd.DataFrame()
            exported_non_cryptic_rows = 0
            exported_non_cryptic_snv_count = 0
            exported_non_cryptic_sns_count = 0
            exported_non_cryptic_ac0_count = 0
            exported_positions = []
            exported_allele_counts = []
            exported_rows_replayed = []
        else:
            dedup_non_cryptic = (
                non_cryptic_cum.sort_values("mm", kind="mergesort")
                .drop_duplicates(subset=["scaffold", "position"], keep="last")
                .sort_index()
            )
            exported_non_cryptic_rows = int(len(dedup_non_cryptic))
            exported_non_cryptic_snv_count = int((dedup_non_cryptic["allele_count"] > 1).sum())
            exported_non_cryptic_sns_count = int((dedup_non_cryptic["allele_count"] == 1).sum())
            exported_non_cryptic_ac0_count = int((dedup_non_cryptic["allele_count"] == 0).sum())
            exported_positions = [int(x) for x in dedup_non_cryptic["position"].tolist()]
            exported_allele_counts = [int(x) for x in dedup_non_cryptic["allele_count"].tolist()]
            exported_rows_replayed = [
                {
                    "scaffold": str(r["scaffold"]),
                    "position": int(r["position"]),
                    "position_coverage": int(r["position_coverage"]),
                    "allele_count": int(r["allele_count"]),
                    "ref_base": str(r["ref_base"]).strip(),
                    "con_base": str(r["con_base"]).strip(),
                    "var_base": str(r["var_base"]).strip(),
                    "ref_freq": float(r["ref_freq"]),
                    "con_freq": float(r["con_freq"]),
                    "var_freq": float(r["var_freq"]),
                    "A": int(r["A"]),
                    "C": int(r["C"]),
                    "T": int(r["T"]),
                    "G": int(r["G"]),
                    "cryptic": False,
                    "class": str(r["class"]).strip(),
                }
                for _, r in dedup_non_cryptic.iterrows()
            ]

        cryptic_dedup = dedup_all_mm[dedup_all_mm["_cryptic_bool"]].sort_values("position", kind="mergesort")
        excluded_cryptic_positions = [
            {
                "position": int(crow["position"]),
                "mm": int(crow["mm"]),
                "cryptic": True,
                "allele_count": int(crow["allele_count"]),
                "class": str(crow["class"]).strip(),
                "ref_base": str(crow["ref_base"]).strip(),
                "con_base": str(crow["con_base"]).strip(),
                "var_base": str(crow["var_base"]).strip(),
                "position_coverage": int(crow["position_coverage"]),
                "var_freq": float(crow["var_freq"]),
            }
            for _, crow in cryptic_dedup.iterrows()
        ]

    if summary_divergent_unfiltered != summary_sns_unfiltered + summary_snv_unfiltered + summary_ac0_unfiltered:
        raise ValueError(
            f"raw cumulative summary partition mismatch for {representative}: "
            f"divergent_site_count={summary_divergent_unfiltered} != "
            f"SNS={summary_sns_unfiltered} + SNV={summary_snv_unfiltered} + ac0={summary_ac0_unfiltered}"
        )
    if summary_divergent_unfiltered != exported_non_cryptic_rows + len(excluded_cryptic_positions):
        raise ValueError(
            f"raw cumulative reconciliation arithmetic mismatch for {representative}: "
            f"summary_divergent_unfiltered={summary_divergent_unfiltered} != "
            f"exported_non_cryptic_rows={exported_non_cryptic_rows} + excluded_cryptic={len(excluded_cryptic_positions)}"
        )

    return {
        "cumulative_snv_table_path": str(raw_snv_path),
        "cumulative_snv_table_sha256": stream_sha256(raw_snv_path),
        "cumulative_scaffold_table_path": str(raw_scaff_path) if raw_scaff_path.exists() else None,
        "cumulative_scaffold_table_sha256": stream_sha256(raw_scaff_path) if raw_scaff_path.exists() else None,
        "raw_snv_rows_target": int(len(target_cum)),
        "raw_scaffold_mm_rows_target": len(summary_counts_all_mm),
        "summary_counts_all_mm": summary_counts_all_mm,
        "summary_at_max_mm": {
            "mm": max_mm,
            "SNV_count": summary_snv_unfiltered,
            "SNS_count": summary_sns_unfiltered,
            "zero_allele_divergent_site_count": summary_ac0_unfiltered,
            "divergent_site_count": summary_divergent_unfiltered,
        },
        "exported_noncryptic_rows_exactly_replayed": exported_non_cryptic_rows,
        "exported_noncryptic_snv_count": exported_non_cryptic_snv_count,
        "exported_noncryptic_sns_count": exported_non_cryptic_sns_count,
        "exported_noncryptic_zero_allele_divergent_site_count": exported_non_cryptic_ac0_count,
        "exported_noncryptic_positions_zero_based": exported_positions,
        "exported_noncryptic_allele_counts": exported_allele_counts,
        "exported_noncryptic_rows_replayed": exported_rows_replayed,
        "excluded_cryptic_positions_zero_based": excluded_cryptic_positions,
        "source_semantics": (
            "calc_snps counts last cumulative record per position at mm<=threshold, without cryptic filter. "
            "SNVprofile.get_nonredundant_snv_table filters cryptic==False BEFORE sorting mm and keeping last per scaffold+position."
        ),
    }


def qualify_instrain_profile_output(
    output_dir: Path,
    representative: str,
    min_cov: int = 5,
    min_cov_provenance: str = "explicit_cli_flag_min_cov_5",
) -> dict[str, Any]:
    """Posthoc qualification of inStrain 1.9.1 profile output per official inStrain documentation
    (https://instrain.readthedocs.io/en/latest/example_output.html) and installed inStrain 1.9.1 source:
    - Requires terminal completion banner `..:: inStrain profile finished ::..` in `log/log.log`.
    - Distinguishes `breadth` (>=1 mapped read) from `breadth_minCov` (positions with >= min_cov depth).
    - Records `round(length * breadth_minCov)` as `inferred_depth_qualified_positions_ge_mincov`
      (and never falls back to `breadth` if `breadth_minCov` is absent).
    - Reports `SNV_count` (polymorphic sites), `SNS_count` (fixed consensus substitutions), and
      `divergent_site_count` separately from `scaffold_info.tsv` (which `calc_snps` computes across
      `mm <= max_mm` WITHOUT filtering `cryptic == True` positions, and where `divergent_site_count`
      equals `SNV_count + SNS_count + zero_allele_divergent_site_count`).
    - Validates `<prefix>_SNVs.tsv` (which `SNVprofile.get_nonredundant_snv_table(cryptic=False)`
      computes by filtering `cryptic == False` BEFORE deduplicating by `(scaffold, position)`) and
      reconciles both distinct counts and all common export fields against `raw_data/cumulative_snv_table.csv.gz`
      when present. If `*SNVs.tsv` is missing on disk, reports `unknown_snv_table_missing` without inventing zero.
    """
    log_path = _require_instrain_completion_banner(
        output_dir,
        INSTRAIN_PROFILE_FINISHED_BANNER,
        f"inStrain profile ({representative})",
    )

    candidates = sorted(output_dir.glob("**/*scaffold_info.tsv"))
    if len(candidates) != 1:
        raise ValueError(f"expected one inStrain scaffold_info table in {output_dir}, found={candidates}")
    scaffold_info_path = candidates[0]
    frame = pd.read_csv(scaffold_info_path, sep="\t")
    require_columns(frame, ["scaffold", "length", "breadth"], "inStrain scaffold_info")
    selected = frame[frame["scaffold"].astype(str) == representative]
    if len(selected) != 1:
        raise ValueError(f"inStrain scaffold_info lacks exactly one row for {representative}")
    row = selected.iloc[0]

    ctx_scaff_info = f"scaffold_info.tsv ({scaffold_info_path}, {representative})"
    length = _parse_strict_int(row["length"], "length", ctx_scaff_info)
    breadth_ge1x = _parse_strict_float(row["breadth"], "breadth", ctx_scaff_info)
    if length <= 0 or not (0.0 <= breadth_ge1x <= 1.0):
        raise ValueError(f"invalid inStrain length/breadth for {representative}: length={length}, breadth={breadth_ge1x}")

    read_covered_positions_ge1x = int(round(length * breadth_ge1x))

    # Depth-qualified positions from actual breadth_minCov (never fallback to breadth)
    if "breadth_minCov" in frame.columns and pd.notna(row["breadth_minCov"]):
        breadth_mincov = _parse_strict_float(row["breadth_minCov"], "breadth_minCov", ctx_scaff_info)
        if not (0.0 <= breadth_mincov <= 1.0):
            raise ValueError(f"invalid inStrain breadth_minCov for {representative}: {breadth_mincov}")
        inferred_depth_qualified: int | str = int(round(length * breadth_mincov))
        breadth_mincov_out: float | str = breadth_mincov
        depth_method = (
            f"round(inStrain scaffold_info length * breadth_minCov) at min_cov={min_cov} ({min_cov_provenance}); "
            "inferred depth-qualified positions (>=min_cov) only — full variant callable status also "
            "requires base quality, read pair ANI, and allele frequency/FDR filters"
        )
    else:
        breadth_mincov_out = "unknown_breadth_minCov_column_absent"
        inferred_depth_qualified = "unknown_breadth_minCov_column_absent"
        depth_method = "unknown_breadth_minCov_column_absent; never falls back to breadth (>=1 read)"

    mean_coverage_instrain: float | str = (
        _parse_strict_float(row["coverage"], "coverage", ctx_scaff_info)
        if "coverage" in frame.columns and pd.notna(row["coverage"])
        else "unknown_coverage_column_absent"
    )
    nucl_diversity_instrain: float | str = (
        _parse_strict_float(row["nucl_diversity"], "nucl_diversity", ctx_scaff_info)
        if "nucl_diversity" in frame.columns and pd.notna(row["nucl_diversity"])
        else "NA_or_not_computed"
    )

    snv_count_scaffold_info: int | None = (
        _parse_strict_int(row["SNV_count"], "SNV_count", ctx_scaff_info)
        if "SNV_count" in frame.columns and pd.notna(row["SNV_count"])
        else None
    )
    sns_count_scaffold_info: int | None = (
        _parse_strict_int(row["SNS_count"], "SNS_count", ctx_scaff_info)
        if "SNS_count" in frame.columns and pd.notna(row["SNS_count"])
        else None
    )
    divergent_site_count_scaffold_info: int | None = (
        _parse_strict_int(row["divergent_site_count"], "divergent_site_count", ctx_scaff_info)
        if "divergent_site_count" in frame.columns and pd.notna(row["divergent_site_count"])
        else None
    )

    raw_reconciliation = reconcile_instrain_raw_cumulative_variants(output_dir, representative)
    if (
        raw_reconciliation is not None
        and snv_count_scaffold_info is not None
        and sns_count_scaffold_info is not None
        and divergent_site_count_scaffold_info is not None
    ):
        raw_max_summary = raw_reconciliation["summary_at_max_mm"]
        if (
            snv_count_scaffold_info != raw_max_summary["SNV_count"]
            or sns_count_scaffold_info != raw_max_summary["SNS_count"]
            or divergent_site_count_scaffold_info != raw_max_summary["divergent_site_count"]
        ):
            raise ValueError(
                f"scaffold_info.tsv variant counts for {representative} "
                f"(SNV={snv_count_scaffold_info}, SNS={sns_count_scaffold_info}, div={divergent_site_count_scaffold_info}) "
                f"contradict replayed raw cumulative_snv_table counts at max_mm={raw_max_summary['mm']} "
                f"({raw_max_summary})"
            )

    # Validate SNVs.tsv according to installed inStrain 1.9.1 SNVprofile.generate("SNVs") contract:
    snv_tables = sorted(output_dir.glob("**/*SNVs.tsv"))
    if len(snv_tables) > 1:
        raise ValueError(f"multiple SNVs.tsv files found in {output_dir}: {snv_tables}")
    exported_non_cryptic_snv_count_out: int | str = "unknown_snv_table_missing"
    exported_non_cryptic_sns_count_out: int | str = "unknown_snv_table_missing"
    exported_non_cryptic_ac0_count_out: int | str = "unknown_snv_table_missing"
    summary_ac0_unfiltered_out: int | str = (
        raw_reconciliation["summary_at_max_mm"]["zero_allele_divergent_site_count"]
        if raw_reconciliation is not None
        else "unknown"
    )
    cryptic_excluded_count_out: int | str = "unknown_snv_table_missing"
    excluded_cryptic_positions_out: list[dict[str, Any]] = (
        raw_reconciliation["excluded_cryptic_positions_zero_based"] if raw_reconciliation is not None else []
    )

    if snv_tables:
        snv_path = snv_tables[0]
        ctx_snv_tsv = f"SNVs.tsv ({snv_path}, {representative})"
        raw_text = snv_path.read_text(encoding="utf-8")
        if raw_text == "" or raw_text.strip() == "":
            # 1-byte newline written by pd.DataFrame().to_csv(index=False, sep='\t') when
            # get_nonredundant_snv_table(cryptic=False) returns an empty DataFrame.
            if (
                snv_count_scaffold_info is None
                or sns_count_scaffold_info is None
                or divergent_site_count_scaffold_info is None
            ):
                snv_table_status = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                snv_target_row_count: int | str = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                snv_target_non_cryptic_row_count: int | str = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                exported_non_cryptic_snv_count_out = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                exported_non_cryptic_sns_count_out = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                exported_non_cryptic_ac0_count_out = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                cryptic_excluded_count_out = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
            elif (
                snv_count_scaffold_info == 0
                and sns_count_scaffold_info == 0
                and divergent_site_count_scaffold_info == 0
            ):
                if raw_reconciliation is not None and raw_reconciliation["exported_noncryptic_rows_exactly_replayed"] != 0:
                    raise ValueError(
                        f"empty SNVs.tsv in {output_dir} contradicts non-empty raw cumulative non-cryptic rows"
                    )
                snv_table_status = "present_empty_dataframe_newline_verified_zero_variants_with_scaffold_info_and_log"
                snv_target_row_count = 0
                snv_target_non_cryptic_row_count = 0
                exported_non_cryptic_snv_count_out = 0
                exported_non_cryptic_sns_count_out = 0
                exported_non_cryptic_ac0_count_out = 0
                summary_ac0_unfiltered_out = 0
                cryptic_excluded_count_out = 0
            elif (
                raw_reconciliation is not None
                and raw_reconciliation["exported_noncryptic_rows_exactly_replayed"] == 0
                and len(raw_reconciliation["excluded_cryptic_positions_zero_based"]) == divergent_site_count_scaffold_info
                and snv_count_scaffold_info
                + sns_count_scaffold_info
                + raw_reconciliation["summary_at_max_mm"]["zero_allele_divergent_site_count"]
                == divergent_site_count_scaffold_info
            ):
                # Cryptic-only exported-empty case verified by raw cumulative table:
                # Zero exported non-cryptic variants in SNVs.tsv, but >0 unfiltered cumulative summary sites in scaffold_info.tsv.
                snv_table_status = "present_empty_dataframe_newline_cryptic_only_verified_with_raw_cumulative_table"
                snv_target_row_count = 0
                snv_target_non_cryptic_row_count = 0
                exported_non_cryptic_snv_count_out = 0
                exported_non_cryptic_sns_count_out = 0
                exported_non_cryptic_ac0_count_out = 0
                cryptic_excluded_count_out = len(raw_reconciliation["excluded_cryptic_positions_zero_based"])
            else:
                raise ValueError(
                    f"empty SNVs.tsv in {output_dir} contradicts scaffold_info counts "
                    f"(SNV_count={snv_count_scaffold_info}, SNS_count={sns_count_scaffold_info}, "
                    f"divergent_site_count={divergent_site_count_scaffold_info})"
                )
        else:
            snv_df = pd.read_csv(snv_path, sep="\t")
            missing_cols = sorted(REQUIRED_SNV_COLUMNS - set(snv_df.columns))
            if missing_cols:
                raise ValueError(f"SNVs.tsv in {output_dir} missing required columns {missing_cols}; actual={list(snv_df.columns)}")
            target_df = snv_df[snv_df["scaffold"].astype(str) == representative].copy()
            if not target_df.empty:
                target_df["position"] = _parse_strict_int_series(target_df["position"], "position", ctx_snv_tsv)
            if "cryptic" in target_df.columns and not target_df.empty:
                cryptic_flags = _parse_strict_bool_series(target_df["cryptic"], "cryptic", ctx_snv_tsv)
                if cryptic_flags.any():
                    raise ValueError(
                        f"malformed SNVs.tsv in {output_dir}: contains {int(cryptic_flags.sum())} cryptic=True rows "
                        f"for {representative}, violating get_nonredundant_snv_table(cryptic=False)"
                    )
                non_cryptic_df = target_df[~cryptic_flags].copy()
            else:
                non_cryptic_df = target_df.copy()
            if not target_df.empty and target_df["position"].duplicated().any():
                raise ValueError(f"malformed SNVs.tsv in {output_dir}: duplicate positions for {representative}")

            if "allele_count" in non_cryptic_df.columns:
                if not non_cryptic_df.empty:
                    ac_series = _parse_strict_int_series(non_cryptic_df["allele_count"], "allele_count", ctx_snv_tsv)
                    non_cryptic_df["allele_count"] = ac_series
                    for _, erow in non_cryptic_df.iterrows():
                        ac = int(erow["allele_count"])
                        if ac < 0:
                            raise ValueError(f"malformed SNVs.tsv row in {output_dir}: negative allele_count={ac}")
                        if "class" in non_cryptic_df.columns:
                            cls = str(erow["class"]).strip()
                            if ac == 0 and cls != "DivergentSite":
                                raise ValueError(f"malformed SNVs.tsv row in {output_dir}: allele_count=0 with class={cls!r}")
                            if ac == 1 and (cls != "SNS" or str(erow["con_base"]).strip() == str(erow["ref_base"]).strip()):
                                raise ValueError(f"malformed SNVs.tsv row in {output_dir}: allele_count=1 with class={cls!r}")
                            if ac > 1 and cls not in {"SNV", "con_SNV", "pop_SNV"}:
                                raise ValueError(f"malformed SNVs.tsv row in {output_dir}: allele_count={ac} with class={cls!r}")
                    exported_non_cryptic_snv_count_out = int((ac_series > 1).sum())
                    exported_non_cryptic_sns_count_out = int((ac_series == 1).sum())
                    exported_non_cryptic_ac0_count_out = int((ac_series == 0).sum())
                else:
                    exported_non_cryptic_snv_count_out = 0
                    exported_non_cryptic_sns_count_out = 0
                    exported_non_cryptic_ac0_count_out = 0
            else:
                exported_non_cryptic_snv_count_out = "unknown_allele_count_column_absent"
                exported_non_cryptic_sns_count_out = "unknown_allele_count_column_absent"
                exported_non_cryptic_ac0_count_out = "unknown_allele_count_column_absent"

            snv_target_row_count = int(len(target_df))
            snv_target_non_cryptic_row_count = int(len(non_cryptic_df))

            if snv_target_row_count == 0:
                if (
                    snv_count_scaffold_info is None
                    or sns_count_scaffold_info is None
                    or divergent_site_count_scaffold_info is None
                ):
                    snv_table_status = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                    snv_target_row_count = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                    snv_target_non_cryptic_row_count = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                    exported_non_cryptic_snv_count_out = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                    exported_non_cryptic_sns_count_out = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                    exported_non_cryptic_ac0_count_out = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                    cryptic_excluded_count_out = "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                elif (
                    snv_count_scaffold_info == 0
                    and sns_count_scaffold_info == 0
                    and divergent_site_count_scaffold_info == 0
                ):
                    if raw_reconciliation is not None and raw_reconciliation["exported_noncryptic_rows_exactly_replayed"] != 0:
                        raise ValueError(
                            f"zero-target-row SNVs.tsv for {representative} contradicts non-empty raw cumulative non-cryptic rows"
                        )
                    snv_table_status = "present_header_verified_zero_target_rows"
                    exported_non_cryptic_snv_count_out = 0
                    exported_non_cryptic_sns_count_out = 0
                    exported_non_cryptic_ac0_count_out = 0
                    summary_ac0_unfiltered_out = 0
                    cryptic_excluded_count_out = 0
                elif (
                    raw_reconciliation is not None
                    and raw_reconciliation["exported_noncryptic_rows_exactly_replayed"] == 0
                    and len(raw_reconciliation["excluded_cryptic_positions_zero_based"]) == divergent_site_count_scaffold_info
                    and snv_count_scaffold_info
                    + sns_count_scaffold_info
                    + raw_reconciliation["summary_at_max_mm"]["zero_allele_divergent_site_count"]
                    == divergent_site_count_scaffold_info
                ):
                    snv_table_status = "present_header_zero_target_rows_cryptic_only_verified_with_raw_cumulative_table"
                    exported_non_cryptic_snv_count_out = 0
                    exported_non_cryptic_sns_count_out = 0
                    exported_non_cryptic_ac0_count_out = 0
                    cryptic_excluded_count_out = len(raw_reconciliation["excluded_cryptic_positions_zero_based"])
                else:
                    raise ValueError(
                        f"header-only/zero-target-row SNVs.tsv for {representative} contradicts scaffold_info "
                        f"(SNV_count={snv_count_scaffold_info}, SNS_count={sns_count_scaffold_info}, "
                        f"divergent_site_count={divergent_site_count_scaffold_info})"
                    )
            else:
                if (
                    snv_count_scaffold_info is not None
                    and sns_count_scaffold_info is not None
                    and divergent_site_count_scaffold_info is not None
                ):
                    if raw_reconciliation is not None:
                        raw_ac0 = int(raw_reconciliation["summary_at_max_mm"]["zero_allele_divergent_site_count"])
                        summary_ac0_unfiltered_out = raw_ac0
                        if snv_count_scaffold_info + sns_count_scaffold_info + raw_ac0 != divergent_site_count_scaffold_info:
                            raise ValueError(
                                f"inconsistent scaffold_info variant counts for {representative}: "
                                f"SNV_count={snv_count_scaffold_info} + SNS_count={sns_count_scaffold_info} + "
                                f"zero_allele_divergent_site_count={raw_ac0} != "
                                f"divergent_site_count={divergent_site_count_scaffold_info}"
                            )
                        _verify_exported_snv_rows_match_raw_replayed(
                            non_cryptic_df,
                            raw_reconciliation["exported_noncryptic_rows_replayed"],
                            output_dir,
                            representative,
                        )
                        if (
                            isinstance(exported_non_cryptic_snv_count_out, int)
                            and isinstance(exported_non_cryptic_sns_count_out, int)
                            and isinstance(exported_non_cryptic_ac0_count_out, int)
                        ):
                            if (
                                exported_non_cryptic_snv_count_out != raw_reconciliation["exported_noncryptic_snv_count"]
                                or exported_non_cryptic_sns_count_out != raw_reconciliation["exported_noncryptic_sns_count"]
                                or exported_non_cryptic_ac0_count_out
                                != raw_reconciliation["exported_noncryptic_zero_allele_divergent_site_count"]
                            ):
                                raise ValueError(
                                    f"exported SNVs.tsv class counts for {representative} do not match replayed raw counts"
                                )
                        cryptic_excluded_count_out = len(raw_reconciliation["excluded_cryptic_positions_zero_based"])
                        if snv_target_non_cryptic_row_count + cryptic_excluded_count_out != divergent_site_count_scaffold_info:
                            raise ValueError(
                                f"SNVs.tsv non-cryptic rows ({snv_target_non_cryptic_row_count}) + "
                                f"raw cryptic sites ({cryptic_excluded_count_out}) != "
                                f"scaffold_info divergent_site_count ({divergent_site_count_scaffold_info})"
                            )
                    else:
                        exp_ac0 = (
                            exported_non_cryptic_ac0_count_out
                            if isinstance(exported_non_cryptic_ac0_count_out, int)
                            else 0
                        )
                        summary_ac0_unfiltered_out = (
                            exp_ac0 if isinstance(exported_non_cryptic_ac0_count_out, int) else "unknown_allele_count_column_absent"
                        )
                        if snv_count_scaffold_info + sns_count_scaffold_info + exp_ac0 != divergent_site_count_scaffold_info:
                            raise ValueError(
                                f"inconsistent scaffold_info variant counts for {representative}: "
                                f"SNV_count={snv_count_scaffold_info} + SNS_count={sns_count_scaffold_info} + "
                                f"zero_allele_divergent_site_count={exp_ac0} != "
                                f"divergent_site_count={divergent_site_count_scaffold_info}"
                            )
                        if snv_target_non_cryptic_row_count != divergent_site_count_scaffold_info:
                            raise ValueError(
                                f"SNVs.tsv non-cryptic row count ({snv_target_non_cryptic_row_count}) for {representative} "
                                f"does not match scaffold_info divergent_site_count ({divergent_site_count_scaffold_info}) "
                                "and no raw_data/cumulative_snv_table.csv.gz is present to verify cryptic sites"
                            )
                        if isinstance(exported_non_cryptic_snv_count_out, int) and isinstance(
                            exported_non_cryptic_sns_count_out, int
                        ):
                            if (
                                exported_non_cryptic_snv_count_out != snv_count_scaffold_info
                                or exported_non_cryptic_sns_count_out != sns_count_scaffold_info
                            ):
                                raise ValueError(
                                    f"exported SNVs.tsv SNV/SNS counts ({exported_non_cryptic_snv_count_out}/"
                                    f"{exported_non_cryptic_sns_count_out}) do not match scaffold_info "
                                    f"({snv_count_scaffold_info}/{sns_count_scaffold_info})"
                                )
                        cryptic_excluded_count_out = 0
                else:
                    cryptic_excluded_count_out = (
                        len(raw_reconciliation["excluded_cryptic_positions_zero_based"])
                        if raw_reconciliation is not None
                        else "unknown_scaffold_info_missing_sns_or_divergent_site_columns"
                    )
                if isinstance(cryptic_excluded_count_out, int) and cryptic_excluded_count_out > 0:
                    snv_table_status = "present_verified_target_rows_reconciled_with_raw_cumulative_cryptic_sites"
                else:
                    snv_table_status = "present_verified_target_rows"
        snv_table_path_out: str | None = str(snv_path)
        snv_table_sha256_out: str | None = stream_sha256(snv_path)
    else:
        # Missing SNVs.tsv file on disk: inStrain 1.9.1 always creates <prefix>_SNVs.tsv during generate("SNVs").
        # Never invent zero from a missing file.
        snv_table_path_out = None
        snv_table_sha256_out = None
        snv_table_status = "unknown_snv_table_missing"
        snv_target_row_count = "unknown_snv_table_missing"
        snv_target_non_cryptic_row_count = "unknown_snv_table_missing"
        exported_non_cryptic_snv_count_out = "unknown_snv_table_missing"
        exported_non_cryptic_sns_count_out = "unknown_snv_table_missing"
        exported_non_cryptic_ac0_count_out = "unknown_snv_table_missing"
        cryptic_excluded_count_out = "unknown_snv_table_missing"

    return {
        "scaffold_info_path": str(scaffold_info_path),
        "scaffold_info_sha256": stream_sha256(scaffold_info_path),
        "log_path": str(log_path),
        "log_sha256": stream_sha256(log_path),
        "log_terminal_banner_verified": INSTRAIN_PROFILE_FINISHED_BANNER,
        "instrain_length_bp": length,
        "breadth_ge1x_fraction": breadth_ge1x,
        "read_covered_positions_ge1x": read_covered_positions_ge1x,
        "read_covered_position_method": "round(inStrain scaffold_info length * breadth); positions with >=1 mapped read",
        "min_cov_threshold": int(min_cov),
        "min_cov_configuration_provenance": min_cov_provenance,
        "breadth_mincov_fraction": breadth_mincov_out,
        "inferred_depth_qualified_positions_ge_mincov": inferred_depth_qualified,
        "depth_qualified_positions_ge_mincov": inferred_depth_qualified,
        "depth_qualified_position_method": depth_method,
        "mean_coverage_instrain": mean_coverage_instrain,
        "nucl_diversity_instrain": nucl_diversity_instrain,
        "snv_count_scaffold_info": snv_count_scaffold_info if snv_count_scaffold_info is not None else "unknown",
        "sns_count_scaffold_info": sns_count_scaffold_info if sns_count_scaffold_info is not None else "unknown",
        "divergent_site_count_scaffold_info": (
            divergent_site_count_scaffold_info if divergent_site_count_scaffold_info is not None else "unknown"
        ),
        "snv_count_scaffold_info_unfiltered_cumulative": (
            snv_count_scaffold_info if snv_count_scaffold_info is not None else "unknown"
        ),
        "sns_count_scaffold_info_unfiltered_cumulative": (
            sns_count_scaffold_info if sns_count_scaffold_info is not None else "unknown"
        ),
        "zero_allele_divergent_site_count_scaffold_info_unfiltered_cumulative": summary_ac0_unfiltered_out,
        "divergent_site_count_scaffold_info_unfiltered_cumulative": (
            divergent_site_count_scaffold_info if divergent_site_count_scaffold_info is not None else "unknown"
        ),
        "snv_table_path": snv_table_path_out,
        "snv_table_sha256": snv_table_sha256_out,
        "snv_table_status": snv_table_status,
        "snv_row_count": snv_target_row_count,
        "snv_target_row_count": snv_target_row_count,
        "snv_target_non_cryptic_row_count": snv_target_non_cryptic_row_count,
        "exported_non_cryptic_snv_count": exported_non_cryptic_snv_count_out,
        "exported_non_cryptic_sns_count": exported_non_cryptic_sns_count_out,
        "exported_non_cryptic_zero_allele_divergent_site_count": exported_non_cryptic_ac0_count_out,
        "cryptic_variant_site_count_excluded_from_exported_snv_table": cryptic_excluded_count_out,
        "excluded_cryptic_positions_zero_based": excluded_cryptic_positions_out,
        "raw_cumulative_variant_reconciliation": raw_reconciliation,
        "variant_count_reconciliation_semantics": (
            "scaffold_info.tsv SNV_count/SNS_count/divergent_site_count (calc_snps in inStrain/profile/snv_utilities.py:249-272) "
            "count the last cumulative_snv_table record per position at mm<=max_mm WITHOUT filtering cryptic==True sites "
            "(where divergent_site_count = SNV_count + SNS_count + zero_allele_divergent_site_count); "
            "exported <prefix>_SNVs.tsv (SNVprofile.get_nonredundant_snv_table in inStrain/SNVprofile.py:495-510) filters "
            "cryptic==False BEFORE sorting mm and keeping last per (scaffold, position). Both counts are reproducible with "
            "distinct meanings: summary divergent_site_count = exported non-cryptic rows + excluded cryptic positions; "
            "zero exported variants does not imply zero summary sites when cryptic sites are present."
        ),
        "variant_interpretation": (
            "sequence polymorphism (SNV_count) and fixed consensus substitutions (SNS_count) only; "
            "inferred depth-qualified positions (round(length * breadth_minCov) at >=min_cov) are necessary "
            "but not sufficient for full variant callability (which also requires base quality, read pair ANI, "
            "and allele frequency/FDR filters); variants are never evidence of recombination"
        ),
    }


def verify_source_primary_profile_against_sample_receipt(
    res_item: dict[str, Any],
    sample_dir: Path,
    sample_receipt: dict[str, Any],
) -> tuple[Path, int, str]:
    """Verify that a source_primary inStrain profile listed in a sample mapping receipt has
    completed status, expected parameters, and intact output file hashes on disk.
    """
    if res_item.get("status") != "profiled_successfully":
        raise ValueError(f"source_primary profile status is not profiled_successfully: {res_item}")
    if res_item.get("analysis_mode") != "source_primary":
        raise ValueError(f"unexpected analysis_mode in sample receipt profile: {res_item.get('analysis_mode')}")
    if int(res_item.get("min_mapq_strict_greater_than", -1)) != 19:
        raise ValueError(f"unexpected min_mapq_strict_greater_than in sample receipt: {res_item}")
    if float(res_item.get("min_read_ani_strict_greater_than", -1.0)) != 0.95:
        raise ValueError(f"unexpected min_read_ani_strict_greater_than in sample receipt: {res_item}")

    prof_dir = Path(res_item["profile_dir"])
    if not prof_dir.exists():
        alt_dir = sample_dir / prof_dir.name
        if alt_dir.exists():
            prof_dir = alt_dir
        else:
            raise FileNotFoundError(f"receipted source_primary profile directory missing: {prof_dir}")

    output_hashes = res_item.get("output_sha256")
    if not isinstance(output_hashes, dict) or not output_hashes:
        raise ValueError(f"source_primary profile in {sample_dir} lacks output_sha256 dictionary")
    has_scaffold_info = False
    has_log = False
    for rel_path, expected_sha in output_hashes.items():
        fpath = prof_dir / rel_path
        if not fpath.exists():
            raise FileNotFoundError(f"receipted source_primary profile file missing: {fpath}")
        observed_sha = stream_sha256(fpath)
        if observed_sha != expected_sha:
            raise ValueError(f"receipted source_primary profile file hash mismatch for {fpath}: {observed_sha} != {expected_sha}")
        if rel_path.endswith("scaffold_info.tsv"):
            has_scaffold_info = True
        if rel_path.endswith("log.log"):
            has_log = True
    if not has_scaffold_info or not has_log:
        raise ValueError(
            f"receipted source_primary output_sha256 in {prof_dir} lacks scaffold_info.tsv or log.log: {sorted(output_hashes)}"
        )

    target_sha = str(res_item["target_raw_sha256"])
    min_cov, min_cov_prov = _extract_min_cov_provenance(
        sample_receipt.get("commands", []),
        target_sha,
        default_min_cov=5,
    )
    return prof_dir, min_cov, min_cov_prov


def verify_posthoc_sensitivity_profile_receipt(
    sens_dir: Path,
    sample: str,
    target_sha: str,
    representative: str,
    expected_bam_sha256: str,
    expected_fasta_sha256: str,
) -> dict[str, Any]:
    """Verify cryptographic receipt and output integrity for a posthoc crosssample_sensitivity
    (ANI > 0.92) inStrain profile directory. Never trusts directory existence alone.
    """
    done_path = sens_dir / "POSTHOC_INSTRAIN_PROFILE_DONE"
    receipt_path = sens_dir / "posthoc_instrain_profile_receipt.json"
    if not done_path.exists() or not receipt_path.exists():
        raise FileNotFoundError(f"sensitivity profile receipt/DONE missing in {sens_dir}")
    done_payload = json.loads(done_path.read_text(encoding="utf-8"))
    if done_payload.get("receipt_sha256") != stream_sha256(receipt_path):
        raise ValueError(f"POSTHOC_INSTRAIN_PROFILE_DONE hash mismatch in {sens_dir}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("status") != "completed_verified" or receipt.get("schema_version") != 1:
        raise ValueError(f"invalid sensitivity profile receipt status/schema in {sens_dir}")
    if (
        receipt.get("sample") != sample
        or receipt.get("target_raw_sha256") != target_sha
        or receipt.get("reference_contig") != representative
        or receipt.get("analysis_mode") != "crosssample_sensitivity"
    ):
        raise ValueError(f"sensitivity profile receipt identity mismatch in {sens_dir}")
    params = receipt.get("parameters", {})
    if (
        int(params.get("min_mapq", -1)) != 19
        or float(params.get("min_read_ani", -1.0)) != 0.92
        or int(params.get("min_cov", -1)) != 5
    ):
        raise ValueError(f"sensitivity profile parameters mismatch in {sens_dir}: {params}")
    inputs = receipt.get("inputs", {})
    if inputs.get("filtered_bam_sha256") != expected_bam_sha256:
        raise ValueError(f"sensitivity profile filtered_bam_sha256 mismatch in {sens_dir}")
    if inputs.get("focal_fasta_sha256") != expected_fasta_sha256:
        raise ValueError(f"sensitivity profile focal_fasta_sha256 mismatch in {sens_dir}")
    output_hashes = receipt.get("output_sha256", {})
    if not isinstance(output_hashes, dict) or not output_hashes:
        raise ValueError(f"sensitivity profile receipt has empty output_sha256 in {sens_dir}")
    for rel_path, expected_sha in output_hashes.items():
        fpath = sens_dir / rel_path
        if not fpath.exists():
            raise FileNotFoundError(f"sensitivity profile output missing: {fpath}")
        if stream_sha256(fpath) != expected_sha:
            raise ValueError(f"sensitivity profile output hash mismatch: {fpath}")
    return receipt


def _resolve_expected_instrain_bam_names(
    qualifying_profiles: list[dict[str, Any]] | list[str],
) -> tuple[dict[str, str], int | None]:
    """Resolve the exact BAM identifiers (`os.path.basename(ISP.get('bam_loc'))`) that
    inStrain 1.9.1 CompareController writes into `name1`/`name2` of `comparisonsTable.tsv`.
    """
    bam_to_sample: dict[str, str] = {}
    expected_lengths: set[int] = set()
    for item in qualifying_profiles:
        if isinstance(item, str):
            sample = item
            expected_bam_name = f"{sample}_mapq20.bam"
        elif isinstance(item, dict):
            sample = str(item["sample"])
            expected_bam_name = f"{sample}_mapq20.bam"
            if "instrain_length_bp" in item and item["instrain_length_bp"] is not None:
                expected_lengths.add(int(item["instrain_length_bp"]))
            prof_dir_str = item.get("profile_dir")
            if prof_dir_str:
                attr_path = Path(str(prof_dir_str)) / "raw_data/attributes.tsv"
                if attr_path.exists():
                    for ln in attr_path.read_text(encoding="utf-8").splitlines():
                        parts = ln.split("\t")
                        if len(parts) >= 2 and parts[0].strip() == "bam_loc":
                            attr_bam_name = os.path.basename(parts[1].strip())
                            if attr_bam_name != expected_bam_name:
                                raise ValueError(
                                    f"profile {prof_dir_str} attributes.tsv bam_loc basename {attr_bam_name!r} "
                                    f"does not match expected sample BAM {expected_bam_name!r}"
                                )
        else:
            raise TypeError(f"unexpected qualifying profile item: {item!r}")
        if expected_bam_name in bam_to_sample:
            raise ValueError(f"duplicate qualifying sample/BAM identifier for compare: {expected_bam_name}")
        bam_to_sample[expected_bam_name] = sample
    if len(expected_lengths) > 1:
        raise ValueError(f"inconsistent instrain_length_bp across qualifying profiles: {sorted(expected_lengths)}")
    expected_len = next(iter(expected_lengths)) if len(expected_lengths) == 1 else None
    return bam_to_sample, expected_len


def parse_instrain_comparisons_table(
    comp_dir: Path,
    representative: str,
    qualifying_samples: list[dict[str, Any]] | list[str],
) -> dict[str, Any]:
    """Verify inStrain 1.9.1 compare completion log and parse actual comparisonsTable.tsv rows
    for the exact representative scaffold:
    - Verifies exact compared pair identities (`name1`, `name2`) against the input profile BAM
      identifiers (`<sample>_mapq20.bam`), rejecting self-comparisons, unknown names, or duplicate/
      reciprocal pairs, and requiring the complete N choose 2 unordered pair set.
    - Verifies `0 <= compared_bases_count <= length` and `percent_genome_compared == compared_bases_count / length`.
    - Preserves explicit `NA_zero_coverage_overlap` for `consensus_SNPs`, `population_SNPs`, `popANI`,
      and `conANI` when `compared_bases_count == 0` (rejecting any non-null ANI at zero overlap).
    - Rejects missing `consensus_SNPs`, `population_SNPs`, `popANI`, or `conANI` when `compared_bases_count > 0`,
      and verifies `conANI == (compared_bases_count - consensus_SNPs) / compared_bases_count` and
      `popANI == (compared_bases_count - population_SNPs) / compared_bases_count`.
    - Qualifies comparisons with `<50%` genome compared as shared-region-only (not whole-genome similarity).
    """
    log_path = _require_instrain_completion_banner(
        comp_dir,
        INSTRAIN_COMPARE_FINISHED_BANNER,
        f"inStrain compare ({representative})",
    )
    comp_tables = sorted(comp_dir.glob("**/*comparisonsTable.tsv"))
    if len(comp_tables) != 1:
        raise ValueError(f"expected exactly one comparisonsTable.tsv in {comp_dir}, found={comp_tables}")
    comp_path = comp_tables[0]
    if comp_path.stat().st_size == 0:
        raise ValueError(f"empty comparisonsTable.tsv in {comp_dir}")

    bam_to_sample, expected_len = _resolve_expected_instrain_bam_names(qualifying_samples)
    expected_bam_names = sorted(bam_to_sample.keys())
    if len(expected_bam_names) < 2:
        raise ValueError(f"inStrain compare requires >=2 qualifying profiles, got {expected_bam_names}")
    expected_unordered_pairs = {
        frozenset({a, b})
        for idx, a in enumerate(expected_bam_names)
        for b in expected_bam_names[idx + 1 :]
    }

    df = pd.read_csv(comp_path, sep="\t")
    missing_cols = sorted(REQUIRED_COMPARE_COLUMNS - set(df.columns))
    if missing_cols:
        raise ValueError(f"comparisonsTable.tsv in {comp_dir} missing columns {missing_cols}; actual={list(df.columns)}")
    target_df = df[df["scaffold"].astype(str) == representative].copy()
    if len(target_df) != len(expected_unordered_pairs):
        raise ValueError(
            f"comparisonsTable.tsv in {comp_dir} has {len(target_df)} rows for {representative}, "
            f"expected {len(expected_unordered_pairs)} pairwise comparisons for {expected_bam_names}"
        )

    seen_unordered_pairs: set[frozenset[str]] = set()
    comparison_rows: list[dict[str, Any]] = []
    for _, row in target_df.iterrows():
        name1 = str(row["name1"]).strip()
        name2 = str(row["name2"]).strip()
        if name1 not in bam_to_sample or name2 not in bam_to_sample:
            raise ValueError(
                f"unexpected profile/BAM name in {comp_path}: ({name1!r}, {name2!r}); expected subset of {expected_bam_names}"
            )
        if name1 == name2:
            raise ValueError(f"self-comparison row ({name1!r}, {name2!r}) in {comp_path}")
        pair_key = frozenset({name1, name2})
        if pair_key in seen_unordered_pairs:
            raise ValueError(f"duplicate or reciprocal comparison row for pair ({name1!r}, {name2!r}) in {comp_path}")
        seen_unordered_pairs.add(pair_key)

        compared_bases = int(row["compared_bases_count"])
        pct_compared = float(row["percent_genome_compared"])
        cov_overlap = float(row["coverage_overlap"])
        length_bp = int(row["length"])
        if (
            length_bp <= 0
            or not (0 <= compared_bases <= length_bp)
            or not math.isfinite(pct_compared)
            or not math.isfinite(cov_overlap)
            or not (0.0 <= pct_compared <= 1.0)
            or not (0.0 <= cov_overlap <= 1.0)
        ):
            raise ValueError(f"invalid comparison length/overlap metrics in {comp_path}: {row.to_dict()}")
        if expected_len is not None and length_bp != expected_len:
            raise ValueError(
                f"comparisonsTable.tsv length ({length_bp}) != profile scaffold length ({expected_len}) for {representative}"
            )
        expected_fraction = compared_bases / length_bp
        if abs(pct_compared - expected_fraction) > 1e-6:
            raise ValueError(
                f"percent_genome_compared ({pct_compared}) != compared_bases_count/length ({expected_fraction}) in {comp_path}"
            )

        pop_ani_raw = row["popANI"]
        con_ani_raw = row["conANI"]
        con_snps_raw = row["consensus_SNPs"]
        pop_snps_raw = row["population_SNPs"]

        if compared_bases == 0:
            if pct_compared != 0.0 or cov_overlap != 0.0:
                raise ValueError(
                    f"non-zero percent_genome_compared ({pct_compared}) or coverage_overlap ({cov_overlap}) "
                    f"at compared_bases_count == 0 in {comp_path}"
                )
            if (pd.notna(pop_ani_raw) and str(pop_ani_raw).strip() != "") or (
                pd.notna(con_ani_raw) and str(con_ani_raw).strip() != ""
            ):
                raise ValueError(
                    f"non-null popANI ({pop_ani_raw!r}) or conANI ({con_ani_raw!r}) at compared_bases_count == 0 in {comp_path}"
                )
            con_snps_out: int | str = "NA_zero_coverage_overlap"
            pop_snps_out: int | str = "NA_zero_coverage_overlap"
            con_ani_out: float | str = "NA_zero_coverage_overlap"
            pop_ani_out: float | str = "NA_zero_coverage_overlap"
            estimability_status = "non_estimable_zero_coverage_overlap"
        else:
            if (
                pd.isna(con_snps_raw)
                or str(con_snps_raw).strip() == ""
                or pd.isna(pop_snps_raw)
                or str(pop_snps_raw).strip() == ""
                or pd.isna(pop_ani_raw)
                or str(pop_ani_raw).strip() == ""
                or pd.isna(con_ani_raw)
                or str(con_ani_raw).strip() == ""
            ):
                raise ValueError(
                    f"missing consensus_SNPs, population_SNPs, popANI, or conANI at positive compared_bases_count={compared_bases} "
                    f"in {comp_path}: {row.to_dict()}"
                )
            con_snps_int = int(float(con_snps_raw))
            pop_snps_int = int(float(pop_snps_raw))
            con_ani_flt = float(con_ani_raw)
            pop_ani_flt = float(pop_ani_raw)
            if (
                not (0 <= pop_snps_int <= con_snps_int <= compared_bases)
                or not math.isfinite(con_ani_flt)
                or not math.isfinite(pop_ani_flt)
                or not (0.0 <= con_ani_flt <= 1.0)
                or not (0.0 <= pop_ani_flt <= 1.0)
            ):
                raise ValueError(f"inconsistent SNP counts or ANI bounds at positive overlap in {comp_path}: {row.to_dict()}")
            expected_con_ani = (compared_bases - con_snps_int) / compared_bases
            expected_pop_ani = (compared_bases - pop_snps_int) / compared_bases
            if abs(con_ani_flt - expected_con_ani) > 1e-6 or abs(pop_ani_flt - expected_pop_ani) > 1e-6:
                raise ValueError(
                    f"conANI/popANI inconsistent with (compared_bases_count - SNPs)/compared_bases_count in {comp_path}: "
                    f"conANI={con_ani_flt} vs {expected_con_ani}, popANI={pop_ani_flt} vs {expected_pop_ani}"
                )
            con_snps_out = con_snps_int
            pop_snps_out = pop_snps_int
            con_ani_out = con_ani_flt
            pop_ani_out = pop_ani_flt
            estimability_status = (
                "estimable_positive_overlap_ge_50pct_genome_compared"
                if pct_compared >= 0.50
                else "positive_overlap_below_50pct_shared_region_only_not_whole_genome"
            )

        comparison_rows.append({
            "scaffold": str(row["scaffold"]),
            "name1": name1,
            "name2": name2,
            "sample1": bam_to_sample[name1],
            "sample2": bam_to_sample[name2],
            "length_bp": length_bp,
            "coverage_overlap_fraction": cov_overlap,
            "compared_bases_count": compared_bases,
            "percent_genome_compared_fraction": pct_compared,
            "pairwise_estimability_status": estimability_status,
            "consensus_SNPs": con_snps_out,
            "population_SNPs": pop_snps_out,
            "conANI": con_ani_out,
            "popANI": pop_ani_out,
            "raw_instrain_consensus_SNPs_column": int(float(con_snps_raw)) if pd.notna(con_snps_raw) and str(con_snps_raw).strip() != "" else None,
            "raw_instrain_population_SNPs_column": int(float(pop_snps_raw)) if pd.notna(pop_snps_raw) and str(pop_snps_raw).strip() != "" else None,
        })

    if seen_unordered_pairs != expected_unordered_pairs:
        raise ValueError(
            f"comparisonsTable.tsv pair set {seen_unordered_pairs} != expected {expected_unordered_pairs} in {comp_path}"
        )

    n_positive_overlap = sum(1 for r in comparison_rows if r["compared_bases_count"] > 0)
    n_zero_overlap = sum(1 for r in comparison_rows if r["compared_bases_count"] == 0)
    return {
        "comparisons_table_path": str(comp_path),
        "comparisons_table_sha256": stream_sha256(comp_path),
        "compare_log_path": str(log_path),
        "compare_log_sha256": stream_sha256(log_path),
        "compare_log_terminal_banner_verified": INSTRAIN_COMPARE_FINISHED_BANNER,
        "n_pairwise_comparisons": len(comparison_rows),
        "n_positive_overlap_pairs": n_positive_overlap,
        "n_zero_overlap_non_estimable_pairs": n_zero_overlap,
        "min_compared_bases_count": min(r["compared_bases_count"] for r in comparison_rows),
        "max_compared_bases_count": max(r["compared_bases_count"] for r in comparison_rows),
        "comparison_rows": comparison_rows,
    }


def verify_posthoc_compare_receipt(
    comp_dir: Path,
    target_sha: str,
    representative: str,
    analysis_mode: str,
    qualifying_profiles: list[dict[str, Any]],
) -> dict[str, Any]:
    """Verify cryptographic receipt and output integrity for an inStrain compare directory."""
    done_path = comp_dir / "POSTHOC_INSTRAIN_COMPARE_DONE"
    receipt_path = comp_dir / "posthoc_instrain_compare_receipt.json"
    if not done_path.exists() or not receipt_path.exists():
        raise FileNotFoundError(f"inStrain compare receipt/DONE missing in {comp_dir}")
    done_payload = json.loads(done_path.read_text(encoding="utf-8"))
    if done_payload.get("receipt_sha256") != stream_sha256(receipt_path):
        raise ValueError(f"POSTHOC_INSTRAIN_COMPARE_DONE hash mismatch in {comp_dir}")
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    if receipt.get("status") != "completed_verified" or receipt.get("schema_version") != 1:
        raise ValueError(f"invalid compare receipt status/schema in {comp_dir}")
    expected_samples = sorted(q["sample"] for q in qualifying_profiles)
    expected_scaff_hashes = {q["sample"]: q["scaffold_info_sha256"] for q in qualifying_profiles}
    if (
        receipt.get("target_seq_sha256") != target_sha
        or receipt.get("reference_contig") != representative
        or receipt.get("analysis_mode") != analysis_mode
        or receipt.get("qualifying_samples") != expected_samples
        or receipt.get("input_profile_scaffold_info_sha256") != expected_scaff_hashes
    ):
        raise ValueError(f"compare receipt identity/input binding mismatch in {comp_dir}")
    output_hashes = receipt.get("output_sha256", {})
    if not isinstance(output_hashes, dict) or not output_hashes:
        raise ValueError(f"compare receipt has empty output_sha256 in {comp_dir}")
    for rel_path, expected_sha in output_hashes.items():
        fpath = comp_dir / rel_path
        if not fpath.exists():
            raise FileNotFoundError(f"compare output file missing: {fpath}")
        if stream_sha256(fpath) != expected_sha:
            raise ValueError(f"compare output hash mismatch: {fpath}")
    return receipt


def _preserve_stale_or_partial_dir(target_dir: Path, reason: str, preserved_list: list[dict[str, Any]]) -> Path:
    """Rename a stale or partial profile/compare directory without deleting historical evidence."""
    stamp = f"{int(time.time())}_{os.getpid()}"
    preserved_path = target_dir.with_name(f"{target_dir.name}.stale_preserved_{stamp}")
    counter = 1
    while preserved_path.exists():
        preserved_path = target_dir.with_name(f"{target_dir.name}.stale_preserved_{stamp}_{counter}")
        counter += 1
    os.replace(target_dir, preserved_path)
    preserved_list.append({
        "original_path": str(target_dir),
        "preserved_path": str(preserved_path),
        "reason": reason,
    })
    return preserved_path


def _extract_faidx_sequence(ref_fna: Path, representative: str) -> str:
    if shutil.which("samtools") is None:
        raise RuntimeError(f"samtools unavailable to extract {representative} from {ref_fna}")
    res = subprocess.run(
        ["samtools", "faidx", str(ref_fna), representative],
        capture_output=True,
        text=True,
        check=True,
    )
    lines = res.stdout.splitlines()
    if not lines or not lines[0].startswith(">") or lines[0][1:].split()[0] != representative:
        raise ValueError(f"samtools faidx returned unexpected header for {representative}")
    return "".join(ln.strip() for ln in lines[1:] if not ln.startswith(">"))


def _verify_focal_fasta_against_target_sha(focal_fasta: Path, representative: str, target_sha: str) -> str:
    """Verify that focal_fasta contains exact header `representative` and sequence matching `target_sha`
    (either direct or reverse-complement, matching the reference orientation contract).
    """
    lines = [ln.strip() for ln in focal_fasta.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if len(lines) < 2 or lines[0].split()[0] != f">{representative}":
        raise ValueError(f"invalid focal FASTA header in {focal_fasta}; expected >{representative}")
    if any(ln.startswith(">") for ln in lines[1:]):
        raise ValueError(f"multiple FASTA records in {focal_fasta}")
    seq = "".join(lines[1:]).upper()
    direct_sha = sha256_sequence(seq)
    rc_sha = sha256_sequence(reverse_complement(seq))
    if target_sha not in (direct_sha, rc_sha):
        raise ValueError(
            f"focal FASTA {focal_fasta} sequence hash ({direct_sha} / rc={rc_sha}) does not match target_sha {target_sha}"
        )
    return stream_sha256(focal_fasta)


def run_crosssample_sensitivity_and_comparisons_if_eligible(
    root: Path,
    threads: int = 16,
    require_complete: bool = True,
    instrain_sif: Path = SIF_INSTRAIN,
    singularity_bin: str | None = None,
) -> dict[str, Any]:
    """Verify exact 14-sample competitive mapping completion, qualify all `source_primary`
    (ANI > 0.95) inStrain profiles with truthful `read_covered_positions_ge1x` vs
    `inferred_depth_qualified_positions_ge_mincov` semantics, execute/verify receipted
    `crosssample_sensitivity` (ANI > 0.92, min_cov=5) profiles from `<sample>_mapq20.bam`
    where eligible, and execute/verify receipted `inStrain compare` runs (or document
    measured 14-sample quantitative non-estimability where <2 samples qualify).
    """
    res_dir = root / "results"
    if not (res_dir / "stage3_tier1_pairs_individual_evidence_audit.tsv").exists():
        alt_dir = root / "compact_evidence"
        if (alt_dir / "stage3_tier1_pairs_individual_evidence_audit.tsv").exists():
            res_dir = alt_dir

    map_dir = res_dir / "stage4_competitive_mapping"
    ref_dir = res_dir / "stage4_competitive_reference"
    ref_fna = ref_dir / "unified_competitive_reference.fna"
    pair_ledger = res_dir / "stage3_tier1_pairs_individual_evidence_audit.tsv"
    boundary_cfg = root / "inputs/stage4_boundary_sensitive_intervals.json"

    pairs, reviewed_pairs, boundary_pairs = load_analysis_pairs(pair_ledger, boundary_cfg)
    planned_samples = load_planned_samples(root, pairs)

    # Check reference validation receipt if present
    val_receipt_path = ref_dir / "stage4_reference_validation_receipt.json"
    val_done_path = ref_dir / "REFERENCE_VALIDATION_DONE"
    expected_val_receipt_sha: str | None = None
    val_receipt_payload: dict[str, Any] | None = None
    if val_receipt_path.exists() or val_done_path.exists():
        if not val_receipt_path.exists() or not val_done_path.exists():
            raise FileNotFoundError(f"incomplete reference validation receipt/sentinel in {ref_dir}")
        val_receipt_payload = json.loads(val_receipt_path.read_text(encoding="utf-8"))
        if val_receipt_payload.get("status") != "completed_verified" or val_receipt_payload.get("schema_version") != 3:
            raise ValueError(f"invalid reference validation receipt in {val_receipt_path}")
        expected_val_receipt_sha = stream_sha256(val_receipt_path)

    # Scan and verify every present sample receipt first (corrupt receipts always fail closed immediately)
    receipt_paths = sorted(map_dir.glob("*/*_mapping_receipt.json")) if map_dir.exists() else []
    verified_receipts_by_sample: dict[str, dict[str, Any]] = {}
    observed_val_binding_shas: set[str] = set()

    for rpath in receipt_paths:
        receipt = verify_sample_receipt(rpath)
        sample = str(receipt["sample"])
        if sample in verified_receipts_by_sample:
            raise ValueError(f"duplicate sample mapping receipt for {sample}")
        if sample not in planned_samples:
            raise ValueError(f"unexpected sample {sample!r} not in planned 14-sample manifest {planned_samples}")

        ref_bind = receipt.get("reference_binding", {})
        bind_sha = str(ref_bind.get("validation_receipt_sha256", ""))
        if len(bind_sha) != 64 or not set(bind_sha) <= SHA256_HEX:
            raise ValueError(f"sample {sample} receipt lacks valid validation_receipt_sha256: {bind_sha!r}")
        if expected_val_receipt_sha is not None and bind_sha != expected_val_receipt_sha:
            raise ValueError(
                f"sample {sample} validation_receipt_sha256 ({bind_sha}) != current validation receipt ({expected_val_receipt_sha})"
            )
        if val_receipt_payload is not None:
            if ref_bind.get("reference_fna_sha256") != val_receipt_payload["reference_files"]["reference_fna"]["sha256"]:
                raise ValueError(f"sample {sample} reference_fna_sha256 mismatch with validation receipt")
            if ref_bind.get("crosswalk_sha256") != val_receipt_payload["reference_files"]["crosswalk_tsv"]["sha256"]:
                raise ValueError(f"sample {sample} crosswalk_sha256 mismatch with validation receipt")
            if ref_bind.get("build_receipt_sha256") != val_receipt_payload["build_receipt"]["sha256"]:
                raise ValueError(f"sample {sample} build_receipt_sha256 mismatch with validation receipt")
        observed_val_binding_shas.add(bind_sha)

        # Reconcile focal_evaluation_tsv rows and inStrain eligibility for this sample
        frame = receipt["_evaluation_frame"]
        for pair in pairs.to_dict(orient="records"):
            _ = one_row(frame, sample=sample, entity_type="host_mag", entity_id=pair["host_mag"])
            virus_row = one_row(frame, sample=sample, entity_type="virus_target", entity_sha256=pair["target_seq_sha256"])
            if str(virus_row["analysis_role"]) != str(pair["analysis_role"]):
                raise ValueError(f"analysis_role mismatch for {pair['target_seq_sha256']} in {sample}")
            depth = float(virus_row["mean_depth"])
            breadth = float(virus_row["horizontal_breadth_pct"])
            unique = str(virus_row["uniquely_attributable"]).strip().lower() in {"true", "1"}
            role = str(virus_row["analysis_role"])
            expected_qualifies = bool(depth >= 10.0 and breadth >= 50.0 and unique and role == "primary")
            reported_qualifies = bool(int(virus_row["qualifies_10x_50pct_instrain"]))
            if reported_qualifies != expected_qualifies:
                raise ValueError(
                    f"qualifies_10x_50pct_instrain flag mismatch for {pair['target_seq_sha256']} in {sample}: "
                    f"reported={reported_qualifies}, recomputed={expected_qualifies}"
                )
            prof_entry = next(
                (r for r in receipt.get("instrain", {}).get("results", []) if r.get("target_raw_sha256") == pair["target_seq_sha256"]),
                None,
            )
            if expected_qualifies and prof_entry is None:
                raise ValueError(
                    f"eligible primary target {pair['target_seq_sha256']} in {sample} lacks an inStrain profile in "
                    "receipt['instrain']['results']; an eligible missing profile is never a zero eligible population"
                )
            if not expected_qualifies and prof_entry is not None:
                raise ValueError(
                    f"ineligible target {pair['target_seq_sha256']} in {sample} unexpectedly has an inStrain profile"
                )

        verified_receipts_by_sample[sample] = receipt

    if len(observed_val_binding_shas) > 1:
        raise ValueError(f"inconsistent validation_receipt_sha256 across sample receipts: {sorted(observed_val_binding_shas)}")

    missing_samples = [s for s in planned_samples if s not in verified_receipts_by_sample]
    if missing_samples:
        if require_complete:
            raise RuntimeError(
                f"incomplete 14-sample competitive mapping receipt set: verified {len(verified_receipts_by_sample)}/14 "
                f"samples, missing={missing_samples}"
            )
        pending_status = (
            "pending_sample_mapping_completion"
            if len(verified_receipts_by_sample) == 0
            else "pending_incomplete_sample_mapping"
        )
        return {
            "schema_version": 2,
            "status": pending_status,
            "planned_samples": planned_samples,
            "n_planned_samples": len(planned_samples),
            "n_sample_receipts_verified": len(verified_receipts_by_sample),
            "verified_samples": [s for s in planned_samples if s in verified_receipts_by_sample],
            "missing_planned_samples": missing_samples,
            "source_primary_profiles": [],
            "crosssample_sensitivity_profiles": [],
            "population_comparisons": [],
        }

    # All 14 planned samples are verified!
    singularity_exec = singularity_bin if singularity_bin is not None else shutil.which("singularity")
    is_real_container_runtime = bool(
        instrain_sif == SIF_INSTRAIN and instrain_sif.exists() and singularity_exec and singularity_bin is None
    )
    overlay_info = verify_and_resolve_instrain_namespaced_header_overlay(
        root=root,
        instrain_sif=instrain_sif,
        singularity_exec=singularity_exec,
        verify_container_runtime=is_real_container_runtime,
    )
    instrain_version = "not_invoked_no_runtime_execution_needed"
    if instrain_sif.exists() and singularity_exec:
        ver_proc = subprocess.run(
            [singularity_exec, "exec", str(instrain_sif), "inStrain", "profile", "--version"],
            capture_output=True,
            text=True,
        )
        if ver_proc.returncode == 0:
            instrain_version = (ver_proc.stdout.strip() or ver_proc.stderr.strip()).splitlines()[0]

    posthoc_script_sha = stream_sha256(Path(__file__))
    preserved_stale_attempts: list[dict[str, Any]] = []
    profiles_qualified: list[dict[str, Any]] = []
    crosssample_profiles: list[dict[str, Any]] = []
    qualifying_primary_by_sha: dict[str, list[dict[str, Any]]] = {}
    qualifying_sensitivity_by_sha: dict[str, list[dict[str, Any]]] = {}

    for sample in planned_samples:
        receipt = verified_receipts_by_sample[sample]
        sample_dir = map_dir / sample
        filtered_bam = Path(receipt["outputs"]["filtered_bam"])
        filtered_bam_sha = str(receipt["outputs"]["filtered_bam_sha256"])

        for res_item in receipt.get("instrain", {}).get("results", []):
            target_sha = str(res_item["target_raw_sha256"])
            rep = str(res_item["reference_contig"])
            prof_dir, min_cov, min_cov_prov = verify_source_primary_profile_against_sample_receipt(
                res_item, sample_dir, receipt
            )
            qual = qualify_instrain_profile_output(
                prof_dir,
                rep,
                min_cov=min_cov,
                min_cov_provenance=min_cov_prov,
            )
            rec = {
                "sample": sample,
                "analysis_mode": "source_primary",
                "min_read_ani_strict_greater_than": 0.95,
                "min_mapq_strict_greater_than": 19,
                "effective_minimum_integer_mapq": 20,
                "target_raw_sha256": target_sha,
                "reference_contig": rep,
                "profile_dir": str(prof_dir),
                "tool_version": str(res_item.get("tool_version", instrain_version)),
                "legacy_mapper_callable_positions_ge1x": res_item.get("callable_positions"),
                **qual,
            }
            profiles_qualified.append(rec)
            qualifying_primary_by_sha.setdefault(target_sha, []).append(rec)

            # Now verify or execute crosssample_sensitivity (ANI > 0.92, min_cov=5) for this qualifying sample+target
            sens_dir = sample_dir / f"inStrain_crosssample_sensitivity_{target_sha[:16]}"
            focal_fasta = sample_dir / f"focal_{target_sha[:16]}.fna"
            if not focal_fasta.exists():
                if not ref_fna.exists():
                    raise RuntimeError(
                        f"eligible crosssample_sensitivity target {target_sha} in {sample} lacks {focal_fasta} and {ref_fna}"
                    )
                seq = _extract_faidx_sequence(ref_fna, rep)
                atomic_text(focal_fasta, f">{rep}\n{seq}\n")
            focal_fasta_sha = _verify_focal_fasta_against_target_sha(focal_fasta, rep, target_sha)

            sens_receipt: dict[str, Any] | None = None
            if sens_dir.exists():
                try:
                    sens_receipt = verify_posthoc_sensitivity_profile_receipt(
                        sens_dir=sens_dir,
                        sample=sample,
                        target_sha=target_sha,
                        representative=rep,
                        expected_bam_sha256=filtered_bam_sha,
                        expected_fasta_sha256=focal_fasta_sha,
                    )
                    _ = qualify_instrain_profile_output(
                        sens_dir,
                        rep,
                        min_cov=5,
                        min_cov_provenance="explicit_cli_flag_min_cov_5",
                    )
                except Exception as exc:
                    _preserve_stale_or_partial_dir(
                        sens_dir,
                        f"rejected stale/partial sensitivity profile: {exc}",
                        preserved_stale_attempts,
                    )
                    sens_receipt = None

            if sens_receipt is None:
                if not instrain_sif.exists() or not singularity_exec:
                    raise RuntimeError(
                        f"eligible crosssample_sensitivity profile required for {sample} / {target_sha}, "
                        f"but runtime container/singularity is unavailable (sif={instrain_sif}, singularity={singularity_exec})"
                    )
                cmd = [
                    singularity_exec,
                    "exec",
                    "-B",
                    "/project",
                    "-B",
                    "/scratch",
                    "-B",
                    overlay_info["singularity_bind_arg"],
                    str(instrain_sif),
                    "inStrain",
                    "profile",
                    str(filtered_bam),
                    str(focal_fasta),
                    "-o",
                    str(sens_dir),
                    "-p",
                    str(threads),
                    "--min_mapq",
                    "19",
                    "--min_read_ani",
                    "0.92",
                    "--min_cov",
                    "5",
                ]
                proc = subprocess.run(cmd, capture_output=True, text=True)
                if proc.returncode != 0:
                    if sens_dir.exists():
                        _preserve_stale_or_partial_dir(
                            sens_dir,
                            f"failed sensitivity profile execution (rc={proc.returncode})",
                            preserved_stale_attempts,
                        )
                    raise RuntimeError(
                        f"inStrain crosssample_sensitivity profile failed for {sample} / {target_sha} "
                        f"(rc={proc.returncode}): {proc.stderr}"
                    )
                # Verify output before writing receipt
                _ = qualify_instrain_profile_output(
                    sens_dir,
                    rep,
                    min_cov=5,
                    min_cov_provenance="explicit_cli_flag_min_cov_5",
                )
                out_hashes = {
                    str(p.relative_to(sens_dir)): stream_sha256(p)
                    for p in sorted(sens_dir.rglob("*"))
                    if p.is_file() and p.name not in {"posthoc_instrain_profile_receipt.json", "POSTHOC_INSTRAIN_PROFILE_DONE"}
                }
                sens_receipt_payload = {
                    "status": "completed_verified",
                    "schema_version": 1,
                    "sample": sample,
                    "target_raw_sha256": target_sha,
                    "reference_contig": rep,
                    "analysis_mode": "crosssample_sensitivity",
                    "tool_version": instrain_version,
                    "command": cmd,
                    "parameters": {
                        "min_mapq": 19,
                        "min_read_ani": 0.92,
                        "min_cov": 5,
                    },
                    "inputs": {
                        "filtered_bam": str(filtered_bam),
                        "filtered_bam_sha256": filtered_bam_sha,
                        "focal_fasta": str(focal_fasta),
                        "focal_fasta_sha256": focal_fasta_sha,
                        "container_path": str(instrain_sif),
                        "instrain_namespaced_header_overlay": overlay_info,
                        "posthoc_script_sha256": posthoc_script_sha,
                    },
                    "output_sha256": out_hashes,
                }
                sens_receipt_path = sens_dir / "posthoc_instrain_profile_receipt.json"
                atomic_json(sens_receipt_path, sens_receipt_payload)
                atomic_text(
                    sens_dir / "POSTHOC_INSTRAIN_PROFILE_DONE",
                    json.dumps({"receipt_sha256": stream_sha256(sens_receipt_path)}) + "\n",
                )
                sens_receipt = verify_posthoc_sensitivity_profile_receipt(
                    sens_dir=sens_dir,
                    sample=sample,
                    target_sha=target_sha,
                    representative=rep,
                    expected_bam_sha256=filtered_bam_sha,
                    expected_fasta_sha256=focal_fasta_sha,
                )

            sens_qual = qualify_instrain_profile_output(
                sens_dir,
                rep,
                min_cov=5,
                min_cov_provenance="explicit_cli_flag_min_cov_5",
            )
            sens_rec = {
                "sample": sample,
                "analysis_mode": "crosssample_sensitivity",
                "min_read_ani_strict_greater_than": 0.92,
                "min_mapq_strict_greater_than": 19,
                "effective_minimum_integer_mapq": 20,
                "target_raw_sha256": target_sha,
                "reference_contig": rep,
                "profile_dir": str(sens_dir),
                "profile_receipt_path": str(sens_dir / "posthoc_instrain_profile_receipt.json"),
                "profile_receipt_sha256": stream_sha256(sens_dir / "posthoc_instrain_profile_receipt.json"),
                "tool_version": str(sens_receipt.get("tool_version", instrain_version)),
                **sens_qual,
            }
            crosssample_profiles.append(sens_rec)
            qualifying_sensitivity_by_sha.setdefault(target_sha, []).append(sens_rec)

    # Evaluate population comparisons across all 14 samples for each primary viral target SHA
    # in BOTH modes: source_primary (ANI > 0.95) and crosssample_sensitivity (ANI > 0.92)
    population_comparisons: list[dict[str, Any]] = []
    primary_shas = sorted(set(reviewed_pairs.loc[reviewed_pairs["analysis_role"] == "primary", "target_seq_sha256"]))

    for target_sha in primary_shas:
        t_rows = reviewed_pairs[reviewed_pairs["target_seq_sha256"] == target_sha]
        target_id = str(t_rows.iloc[0]["target_id"])
        host_mags = sorted(set(t_rows["host_mag"].astype(str)))

        # Collect measured 14-sample evaluation rows for this target_sha from the verified sample receipts
        measured_sample_rows: list[dict[str, Any]] = []
        rep_contig: str | None = None
        for sample in planned_samples:
            frame = verified_receipts_by_sample[sample]["_evaluation_frame"]
            vrow = one_row(frame, sample=sample, entity_type="virus_target", entity_sha256=target_sha)
            rep_contig = str(vrow["reference_contig"])
            measured_sample_rows.append({
                "sample": sample,
                "reference_contig": rep_contig,
                "filtered_reads_mapped": int(vrow["filtered_reads_mapped"]),
                "horizontal_breadth_pct": float(vrow["horizontal_breadth_pct"]),
                "mean_depth": float(vrow["mean_depth"]),
                "uniquely_attributable": str(vrow["uniquely_attributable"]).strip().lower() in {"true", "1"},
                "ambiguity_status": str(vrow["ambiguity_status"]),
                "qualifies_10x_50pct_instrain": int(vrow["qualifies_10x_50pct_instrain"]),
            })
        assert rep_contig is not None

        for mode, min_ani, qual_map in [
            ("source_primary", 0.95, qualifying_primary_by_sha),
            ("crosssample_sensitivity", 0.92, qualifying_sensitivity_by_sha),
        ]:
            qual_samples = sorted(qual_map.get(target_sha, []), key=lambda x: str(x["sample"]))
            if len(qual_samples) >= 2:
                comp_dir = map_dir / f"inStrain_compare_{mode}_{target_sha[:16]}"
                comp_receipt: dict[str, Any] | None = None
                if comp_dir.exists():
                    try:
                        comp_receipt = verify_posthoc_compare_receipt(
                            comp_dir=comp_dir,
                            target_sha=target_sha,
                            representative=rep_contig,
                            analysis_mode=mode,
                            qualifying_profiles=qual_samples,
                        )
                        _ = parse_instrain_comparisons_table(
                            comp_dir=comp_dir,
                            representative=rep_contig,
                            qualifying_samples=qual_samples,
                        )
                    except Exception as exc:
                        _preserve_stale_or_partial_dir(
                            comp_dir,
                            f"rejected stale/partial compare dir: {exc}",
                            preserved_stale_attempts,
                        )
                        comp_receipt = None

                if comp_receipt is None:
                    if not instrain_sif.exists() or not singularity_exec:
                        raise RuntimeError(
                            f"eligible inStrain compare required for {target_id} ({target_sha}, mode={mode}) "
                            f"across {[q['sample'] for q in qual_samples]}, but container/singularity unavailable"
                        )
                    input_dirs = [str(q["profile_dir"]) for q in qual_samples]
                    cmd = [
                        singularity_exec,
                        "exec",
                        "-B",
                        "/project",
                        "-B",
                        "/scratch",
                        "-B",
                        overlay_info["singularity_bind_arg"],
                        str(instrain_sif),
                        "inStrain",
                        "compare",
                        "-i",
                        *input_dirs,
                        "-o",
                        str(comp_dir),
                        "-p",
                        str(threads),
                    ]
                    proc = subprocess.run(cmd, capture_output=True, text=True)
                    if proc.returncode != 0:
                        if comp_dir.exists():
                            _preserve_stale_or_partial_dir(
                                comp_dir,
                                f"failed compare execution (rc={proc.returncode})",
                                preserved_stale_attempts,
                            )
                        raise RuntimeError(
                            f"inStrain compare failed for {target_id} ({target_sha}, mode={mode}, rc={proc.returncode}): {proc.stderr}"
                        )
                    parsed_comp = parse_instrain_comparisons_table(
                        comp_dir=comp_dir,
                        representative=rep_contig,
                        qualifying_samples=qual_samples,
                    )
                    out_hashes = {
                        str(p.relative_to(comp_dir)): stream_sha256(p)
                        for p in sorted(comp_dir.rglob("*"))
                        if p.is_file() and p.name not in {"posthoc_instrain_compare_receipt.json", "POSTHOC_INSTRAIN_COMPARE_DONE"}
                    }
                    comp_receipt_payload = {
                        "status": "completed_verified",
                        "schema_version": 1,
                        "target_id": target_id,
                        "target_seq_sha256": target_sha,
                        "reference_contig": rep_contig,
                        "analysis_mode": mode,
                        "min_read_ani_strict_greater_than": min_ani,
                        "qualifying_samples": [q["sample"] for q in qual_samples],
                        "input_profile_dirs": input_dirs,
                        "input_profile_scaffold_info_sha256": {
                            q["sample"]: q["scaffold_info_sha256"] for q in qual_samples
                        },
                        "tool_version": instrain_version,
                        "command": cmd,
                        "inputs": {
                            "container_path": str(instrain_sif),
                            "instrain_namespaced_header_overlay": overlay_info,
                            "posthoc_script_sha256": posthoc_script_sha,
                        },
                        "output_sha256": out_hashes,
                    }
                    comp_receipt_path = comp_dir / "posthoc_instrain_compare_receipt.json"
                    atomic_json(comp_receipt_path, comp_receipt_payload)
                    atomic_text(
                        comp_dir / "POSTHOC_INSTRAIN_COMPARE_DONE",
                        json.dumps({"receipt_sha256": stream_sha256(comp_receipt_path)}) + "\n",
                    )
                    comp_receipt = verify_posthoc_compare_receipt(
                        comp_dir=comp_dir,
                        target_sha=target_sha,
                        representative=rep_contig,
                        analysis_mode=mode,
                        qualifying_profiles=qual_samples,
                    )

                parsed_comp = parse_instrain_comparisons_table(
                    comp_dir=comp_dir,
                    representative=rep_contig,
                    qualifying_samples=qual_samples,
                )
                population_comparisons.append({
                    "target_id": target_id,
                    "target_seq_sha256": target_sha,
                    "reference_contig": rep_contig,
                    "host_mags": host_mags,
                    "analysis_mode": mode,
                    "min_read_ani_strict_greater_than": min_ani,
                    "n_qualifying_samples_10x_50pct": len(qual_samples),
                    "qualifying_samples": [q["sample"] for q in qual_samples],
                    "comparison_status": "completed_verified_instrain_compare",
                    "compare_receipt_path": str(comp_dir / "posthoc_instrain_compare_receipt.json"),
                    "compare_receipt_sha256": stream_sha256(comp_dir / "posthoc_instrain_compare_receipt.json"),
                    **parsed_comp,
                    "measured_14_sample_evaluation_rows": measured_sample_rows,
                    "scientific_interpretation": (
                        "Pairwise inStrain population comparison across qualifying samples at actual overlapping "
                        "depth-qualified positions (compared_bases_count); popANI/conANI and SNP counts measure "
                        "sequence similarity/polymorphism only and are never evidence of recombination."
                    ),
                })
            else:
                population_comparisons.append({
                    "target_id": target_id,
                    "target_seq_sha256": target_sha,
                    "reference_contig": rep_contig,
                    "host_mags": host_mags,
                    "analysis_mode": mode,
                    "min_read_ani_strict_greater_than": min_ani,
                    "n_qualifying_samples_10x_50pct": len(qual_samples),
                    "qualifying_samples": [q["sample"] for q in qual_samples],
                    "comparison_status": "quantitatively_non_estimable_fewer_than_two_qualifying_sample_populations",
                    "reason": (
                        f"Across all 14 competitively mapped samples, {len(qual_samples)} sample(s) "
                        f"({[q['sample'] for q in qual_samples]}) met the >=10x mean depth, >=50% breadth, "
                        "and unique viral ownership gate for this primary viral target; cross-sample inStrain "
                        "population comparison requires >=2 qualifying sample populations."
                    ),
                    "measured_14_sample_evaluation_rows": measured_sample_rows,
                })

    # Final aggregate completion gate: every comparison must be either completed_verified_instrain_compare
    # or quantitatively_non_estimable_fewer_than_two_qualifying_sample_populations
    allowed_comparison_statuses = {
        "completed_verified_instrain_compare",
        "quantitatively_non_estimable_fewer_than_two_qualifying_sample_populations",
    }
    incomplete_comparisons = [
        c for c in population_comparisons if c.get("comparison_status") not in allowed_comparison_statuses
    ]
    if incomplete_comparisons:
        raise RuntimeError(f"incomplete population comparisons prevent completed_verified status: {incomplete_comparisons}")

    if len(crosssample_profiles) != len(profiles_qualified):
        raise RuntimeError(
            f"mismatch between qualified source_primary profiles ({len(profiles_qualified)}) "
            f"and crosssample_sensitivity profiles ({len(crosssample_profiles)})"
        )

    summary_out = {
        "schema_version": 2,
        "status": "completed_verified",
        "instrain_tool_version": instrain_version,
        "instrain_namespaced_header_overlay": overlay_info,
        "planned_samples": planned_samples,
        "n_planned_samples": len(planned_samples),
        "n_sample_receipts_verified": len(verified_receipts_by_sample),
        "verified_samples": planned_samples,
        "missing_planned_samples": [],
        "reference_validation_receipt_sha256": next(iter(observed_val_binding_shas)) if observed_val_binding_shas else None,
        "sample_mapping_receipt_sha256": {
            s: stream_sha256(map_dir / s / f"{s}_mapping_receipt.json") for s in planned_samples
        },
        "n_source_primary_profiles_qualified": len(profiles_qualified),
        "n_crosssample_sensitivity_profiles_qualified": len(crosssample_profiles),
        "source_primary_profiles": profiles_qualified,
        "crosssample_sensitivity_profiles": crosssample_profiles,
        "population_comparisons": population_comparisons,
        "preserved_stale_or_partial_attempts": preserved_stale_attempts,
        "callable_semantics_note": (
            "Official inStrain 1.9.1 scaffold_info output distinguishes breadth (>=1 mapped read) from "
            "breadth_minCov (positions with >= min_cov depth, verified default 5x in inStrain 1.9.1 argumentParser "
            "and explicitly set via --min_cov 5 in posthoc sensitivity runs). Legacy mapper field `callable_positions` "
            "computed round(length * breadth) and is preserved here as `read_covered_positions_ge1x`, while "
            "`inferred_depth_qualified_positions_ge_mincov` is derived strictly from round(length * breadth_minCov) "
            "(or reported as explicit unknown if absent). Full variant callability also depends on base quality, "
            "read ANI, and allele frequency/FDR filters. In installed inStrain 1.9.1, scaffold_info.tsv summary columns "
            "(`SNV_count`, `SNS_count`, `divergent_site_count`) are computed by `calc_snps` across `mm <= max_mm` WITHOUT "
            "filtering `cryptic == True` sites, whereas exported `<prefix>_SNVs.tsv` (`SNVprofile.get_nonredundant_snv_table`) "
            "filters `cryptic == False` BEFORE deduplicating by `(scaffold, position)`. Both counts are retained with "
            "explicit labels (`*_unfiltered_cumulative` vs `exported_non_cryptic_*` plus `cryptic_variant_site_count_excluded_from_exported_snv_table`) "
            "and validated against `raw_data/cumulative_snv_table.csv.gz`. Profile completion requires the exact "
            "terminal banner `..:: inStrain profile finished ::..` in log/log.log, and missing `*SNVs.tsv` files on "
            "disk are reported as `unknown_snv_table_missing` without inventing zero."
        ),
    }
    out_json = res_dir / "stage4_instrain_posthoc_callable_and_sensitivity_summary.json"
    atomic_json(out_json, summary_out)

    tsv_cols = [
        "sample",
        "analysis_mode",
        "min_read_ani_strict_greater_than",
        "effective_minimum_integer_mapq",
        "target_raw_sha256",
        "reference_contig",
        "instrain_length_bp",
        "breadth_ge1x_fraction",
        "read_covered_positions_ge1x",
        "min_cov_threshold",
        "min_cov_configuration_provenance",
        "breadth_mincov_fraction",
        "inferred_depth_qualified_positions_ge_mincov",
        "depth_qualified_positions_ge_mincov",
        "mean_coverage_instrain",
        "nucl_diversity_instrain",
        "snv_count_scaffold_info",
        "sns_count_scaffold_info",
        "divergent_site_count_scaffold_info",
        "snv_count_scaffold_info_unfiltered_cumulative",
        "sns_count_scaffold_info_unfiltered_cumulative",
        "zero_allele_divergent_site_count_scaffold_info_unfiltered_cumulative",
        "divergent_site_count_scaffold_info_unfiltered_cumulative",
        "snv_table_status",
        "snv_target_row_count",
        "snv_target_non_cryptic_row_count",
        "exported_non_cryptic_snv_count",
        "exported_non_cryptic_sns_count",
        "exported_non_cryptic_zero_allele_divergent_site_count",
        "cryptic_variant_site_count_excluded_from_exported_snv_table",
        "variant_count_reconciliation_semantics",
    ]
    out_tsv = res_dir / "stage4_instrain_posthoc_profiles_table.tsv"
    all_rows = profiles_qualified + crosssample_profiles
    with out_tsv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=tsv_cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    return summary_out


if __name__ == "__main__":
    out = run_crosssample_sensitivity_and_comparisons_if_eligible(BRANCH_ROOT, require_complete=True)
    run_host_array_provenance_audit(BRANCH_ROOT)
    print(f"Completed posthoc inStrain & host-array followup: status={out['status']}")
