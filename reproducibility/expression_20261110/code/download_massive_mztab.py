#!/usr/bin/env python3
"""
Download and strictly validate all 30 MS-GF+ mzTab files for MassIVE MSV000097826 / PXD063748
on the Shaheen login node (ROOT-DOWNLOAD-REVIEW.txt compliant).

Validation controls enforced on BOTH cached and newly downloaded files:
1. Filename -> sample_id provenance check against the 30 samples in essdive_3365729/samples.csv.
2. Trailing newline check (file must end with '\\n' and not be truncated mid-line).
3. HTML/error-page rejection (rejects files starting with '<' or lacking 'MTD\\tmzTab-version').
4. Exact PSH header column verification: requires 'sequence', 'PSM_ID', 'accession',
   'opt_global_cv_MS:1002217_decoy_peptide', 'opt_global_cv_MS:1002354_PSM-level_q-value'.
5. Exact row-count verification against upstream MassIVE result metadata (PSM_rows and Protein_rows
   from massive_MSV000097826_result_list.json), plus last-line column-count check against PSH.
6. Re-validation against previously recorded SHA-256 (if present in prior manifest) to detect
   same-length file corruption.
"""
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from typing import Dict, Optional, Set, Tuple

EXP_ROOT = "/scratch/project/k10018/eq-protein-external-20261110/expression"
BASE_DIR = os.path.join(EXP_ROOT, "inputs", "karaoz2026", "massive_MSV000097826")
MANIFEST_JSON = os.path.join(BASE_DIR, "massive_MSV000097826_result_list.json")
SAMPLES_CSV = os.path.join(EXP_ROOT, "inputs", "karaoz2026", "essdive_3365729", "samples.csv")
OUT_MANIFEST_TSV = os.path.join(BASE_DIR, "mztab_download_manifest.tsv")

REQUIRED_PSH_COLS = {
    "sequence",
    "PSM_ID",
    "accession",
    "opt_global_cv_MS:1002217_decoy_peptide",
    "opt_global_cv_MS:1002354_PSM-level_q-value",
}


def load_expected_samples(samples_csv_path: str) -> Set[str]:
    expected = set()
    if not os.path.exists(samples_csv_path):
        return expected
    with open(samples_csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sid = row.get("proteome-EMSL-sampleid", "").strip()
            if sid:
                expected.add(sid)
    return expected


def load_previous_verified_hashes(manifest_tsv_path: str) -> Dict[str, str]:
    prev = {}
    if not os.path.exists(manifest_tsv_path):
        return prev
    try:
        with open(manifest_tsv_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("status") in ("cached_verified_ok", "downloaded_verified_ok") and row.get("sha256"):
                    prev[row["sample_file"]] = row["sha256"].strip()
    except Exception:
        pass
    return prev


def validate_mztab_file(
    path: str,
    expected_psm_rows: int,
    expected_prt_rows: int,
    expected_sample_id: str,
    previous_sha256: Optional[str] = None,
) -> Tuple[bool, str, dict]:
    """
    Stream-validate an mzTab file against upstream MassIVE metadata and format integrity rules.
    Returns (is_valid, reason, stats_dict).
    """
    stats = {
        "bytes": 0,
        "sha256": "",
        "md5": "",
        "actual_psm_rows": 0,
        "actual_prt_rows": 0,
        "psh_col_count": 0,
    }
    if not os.path.exists(path):
        return False, "file_missing", stats

    size = os.path.getsize(path)
    stats["bytes"] = size
    if size == 0:
        return False, "empty_file", stats

    # Check trailing newline at end of file
    with open(path, "rb") as bf:
        bf.seek(max(0, size - 64))
        tail_bytes = bf.read()
        if not tail_bytes.endswith(b"\n"):
            return False, "missing_trailing_newline_truncated", stats

    h_sha = hashlib.sha256()
    h_md5 = hashlib.md5()
    has_mtd_version = False
    has_sample_in_mtd = False
    psh_cols = None
    actual_prt = 0
    actual_psm = 0
    last_nonempty_line = ""

    try:
        with open(path, "rb") as bf:
            for idx, raw_line in enumerate(bf):
                h_sha.update(raw_line)
                h_md5.update(raw_line)
                line = raw_line.decode("utf-8", errors="replace").rstrip("\r\n")
                if not line:
                    continue
                if idx == 0 and (line.lstrip().startswith("<") or "<!DOCTYPE" in line or "<html" in line.lower()):
                    return False, "html_error_payload_not_mztab", stats
                last_nonempty_line = line
                if line.startswith("MTD\t"):
                    if "mzTab-version" in line:
                        has_mtd_version = True
                    if expected_sample_id and expected_sample_id in line:
                        has_sample_in_mtd = True
                elif line.startswith("PRT\t"):
                    actual_prt += 1
                elif line.startswith("PSH\t"):
                    psh_cols = line.split("\t")
                    stats["psh_col_count"] = len(psh_cols)
                    missing_cols = REQUIRED_PSH_COLS - set(psh_cols)
                    if missing_cols:
                        return False, f"missing_required_psh_columns:{sorted(missing_cols)}", stats
                elif line.startswith("PSM\t"):
                    actual_psm += 1
    except Exception as exc:
        return False, f"read_exception:{exc}", stats

    sha_hex = h_sha.hexdigest()
    md5_hex = h_md5.hexdigest()
    stats["sha256"] = sha_hex
    stats["md5"] = md5_hex
    stats["actual_prt_rows"] = actual_prt
    stats["actual_psm_rows"] = actual_psm

    if not has_mtd_version:
        return False, "missing_mtd_mztab_version", stats
    if expected_sample_id and not has_sample_in_mtd:
        return False, f"sample_id_not_in_mtd:{expected_sample_id}", stats
    if psh_cols is None:
        return False, "missing_psh_header", stats
    if not last_nonempty_line.startswith("PSM\t"):
        return False, "last_line_not_psm_row", stats
    last_cols = last_nonempty_line.split("\t")
    if len(last_cols) != len(psh_cols):
        return False, f"last_psm_col_count_mismatch:{len(last_cols)}_vs_{len(psh_cols)}", stats
    if expected_psm_rows > 0 and actual_psm != expected_psm_rows:
        return False, f"psm_row_count_mismatch:actual_{actual_psm}_vs_expected_{expected_psm_rows}", stats
    if expected_prt_rows > 0 and actual_prt != expected_prt_rows:
        return False, f"prt_row_count_mismatch:actual_{actual_prt}_vs_expected_{expected_prt_rows}", stats
    if previous_sha256 and sha_hex != previous_sha256:
        return False, f"sha256_mismatch_with_previously_verified:{sha_hex}_vs_{previous_sha256}", stats

    return True, "verified_ok", stats


def download_and_validate_one(
    row: dict,
    expected_samples: Set[str],
    prev_hashes: Dict[str, str],
) -> dict:
    file_desc = row["File_descriptor"]
    fname = os.path.basename(file_desc)
    suffix = "_msgfplus.mzTab"
    sample_id = fname[:-len(suffix)] if fname.endswith(suffix) else fname
    in_expected_30 = (not expected_samples) or (sample_id in expected_samples)
    exp_psm = int(row.get("PSM_rows", "0") or 0)
    exp_prt = int(row.get("Protein_rows", "0") or 0)
    out_path = os.path.join(BASE_DIR, fname)
    https_url = f"https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file={file_desc}&forceDownload=true"
    prev_sha = prev_hashes.get(fname)

    if not in_expected_30:
        return {
            "sample_file": fname,
            "sample_id": sample_id,
            "mztab_id": row.get("MzTab_file", ""),
            "bytes": 0,
            "sha256": "",
            "md5": "",
            "expected_psm_rows": exp_psm,
            "actual_psm_rows": 0,
            "expected_prt_rows": exp_prt,
            "actual_prt_rows": 0,
            "status": f"error:unexpected_sample_id:{sample_id}",
            "url": https_url,
            "path": out_path,
        }

    # 1. Validate existing cached file if present
    if os.path.exists(out_path):
        ok, reason, stats = validate_mztab_file(out_path, exp_psm, exp_prt, sample_id, prev_sha)
        if ok:
            return {
                "sample_file": fname,
                "sample_id": sample_id,
                "mztab_id": row.get("MzTab_file", ""),
                "bytes": stats["bytes"],
                "sha256": stats["sha256"],
                "md5": stats["md5"],
                "expected_psm_rows": exp_psm,
                "actual_psm_rows": stats["actual_psm_rows"],
                "expected_prt_rows": exp_prt,
                "actual_prt_rows": stats["actual_prt_rows"],
                "status": "cached_verified_ok",
                "url": https_url,
                "path": out_path,
            }
        else:
            print(f"  [CACHE INVALID] {fname}: {reason} -> removing and re-downloading", flush=True)
            try:
                os.remove(out_path)
            except OSError:
                pass

    # 2. Sequential HTTPS download with bounded backoff on 429
    tmp_path = out_path + ".part"
    last_err = "unknown"
    for attempt in range(1, 9):
        try:
            subprocess.run(
                [
                    "/usr/bin/curl",
                    "-fsSL",
                    "--retry",
                    "3",
                    "--retry-delay",
                    "5",
                    "--max-time",
                    "300",
                    "-o",
                    tmp_path,
                    https_url,
                ],
                check=True,
            )
            ok, reason, stats = validate_mztab_file(tmp_path, exp_psm, exp_prt, sample_id, None)
            if not ok:
                last_err = f"downloaded_file_failed_validation:{reason}"
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass
                time.sleep(4 * attempt)
                continue

            os.replace(tmp_path, out_path)
            time.sleep(2.5)
            return {
                "sample_file": fname,
                "sample_id": sample_id,
                "mztab_id": row.get("MzTab_file", ""),
                "bytes": stats["bytes"],
                "sha256": stats["sha256"],
                "md5": stats["md5"],
                "expected_psm_rows": exp_psm,
                "actual_psm_rows": stats["actual_psm_rows"],
                "expected_prt_rows": exp_prt,
                "actual_prt_rows": stats["actual_prt_rows"],
                "status": "downloaded_verified_ok",
                "url": https_url,
                "path": out_path,
            }
        except Exception as exc:
            last_err = f"download_error:{exc}"
            time.sleep(5 * attempt)

    return {
        "sample_file": fname,
        "sample_id": sample_id,
        "mztab_id": row.get("MzTab_file", ""),
        "bytes": 0,
        "sha256": "",
        "md5": "",
        "expected_psm_rows": exp_psm,
        "actual_psm_rows": 0,
        "expected_prt_rows": exp_prt,
        "actual_prt_rows": 0,
        "status": f"error:{last_err}",
        "url": https_url,
        "path": out_path,
    }


def main() -> int:
    with open(MANIFEST_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)
    rows = data.get("blockData", [])
    expected_samples = load_expected_samples(SAMPLES_CSV)
    prev_hashes = load_previous_verified_hashes(OUT_MANIFEST_TSV)

    print(
        f"Starting sequential download & full-file PSM/PRT validation of {len(rows)} mzTab files "
        f"(expected samples={len(expected_samples)})...",
        flush=True,
    )
    results = []
    for idx, r in enumerate(rows, start=1):
        res = download_and_validate_one(r, expected_samples, prev_hashes)
        results.append(res)
        print(
            f"  [{idx}/{len(rows)}] {res['sample_file']} -> {res['status']} "
            f"(bytes={res['bytes']}, PSMs={res['actual_psm_rows']}/{res['expected_psm_rows']}, "
            f"PRTs={res['actual_prt_rows']}/{res['expected_prt_rows']})",
            flush=True,
        )

    results.sort(key=lambda x: x["sample_file"])
    cols = [
        "sample_file",
        "sample_id",
        "mztab_id",
        "bytes",
        "sha256",
        "md5",
        "expected_psm_rows",
        "actual_psm_rows",
        "expected_prt_rows",
        "actual_prt_rows",
        "status",
        "url",
        "path",
    ]
    with open(OUT_MANIFEST_TSV, "w", encoding="utf-8") as out_f:
        out_f.write("\t".join(cols) + "\n")
        for r in results:
            out_f.write("\t".join(str(r[c]) for c in cols) + "\n")

    verified_samples = {r["sample_id"] for r in results if r["status"].endswith("_verified_ok")}
    n_ok = len(verified_samples)
    missing_samples = sorted(expected_samples - verified_samples) if expected_samples else []
    print(f"Completed: {n_ok}/{len(rows)} verified. Missing expected samples: {len(missing_samples)}", flush=True)
    return 0 if (n_ok == len(rows) and not missing_samples) else 1


if __name__ == "__main__":
    sys.exit(main())
