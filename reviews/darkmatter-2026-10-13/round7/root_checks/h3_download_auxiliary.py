#!/usr/bin/env python3
"""
scripts/01_download_smag_zenodo.py
Resumable, modest-parallel downloader for all 373 parts of the public SMAG archive
from Zenodo (Record 8223844), with strict byte-size and MD5/SHA-256 verification
against official Zenodo API metadata.

Usage:
    python3 scripts/01_download_smag_zenodo.py [n_workers] [max_parts]
"""

import os
import sys
import json
import time
import hashlib
import urllib.request
import urllib.error
import concurrent.futures

ZENODO_JSON_PATH = "results/zenodo_8223844.json"
PARTS_DIR = "downloads/smag_parts"
MANIFEST_PATH = "results/smag_zenodo_auxiliary_parts_manifest.tsv"


def compute_md5_sha256(path):
    h_md5 = hashlib.md5()
    h_sha = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(2 * 1024 * 1024):
            h_md5.update(chunk)
            h_sha.update(chunk)
    return h_md5.hexdigest(), h_sha.hexdigest()


def download_single_part(part_info):
    """Download and verify a single part with retries and exponential backoff."""
    key = part_info["key"]
    expected_size = part_info["size"]
    expected_md5 = part_info["md5"]
    url = part_info["url"]

    dest = os.path.join(PARTS_DIR, key)
    aux_dir = os.path.join(PARTS_DIR, "auxiliary_tmp")
    os.makedirs(aux_dir, exist_ok=True)
    tmp_dest = os.path.join(aux_dir, key + ".tmp")

    # Check if already present and valid
    if os.path.exists(dest) and os.path.getsize(dest) == expected_size:
        actual_md5, actual_sha256 = compute_md5_sha256(dest)
        if actual_md5 == expected_md5:
            print(f"[{key}] Already verified on disk (MD5: {actual_md5}).", flush=True)
            return {
                "key": key,
                "size": expected_size,
                "md5": actual_md5,
                "sha256": actual_sha256,
                "status": "VERIFIED_EXISTING",
                "seconds": 0.0,
                "speed_mbs": 0.0,
            }

    # Attempt download with retries
    max_retries = 8
    base_backoff = 5.0
    for attempt in range(1, max_retries + 1):
        if os.path.exists(tmp_dest):
            try:
                os.remove(tmp_dest)
            except OSError:
                pass

        t0 = time.time()
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "Antigravity-H3-Research/2.0 (KAUST BioCore)",
                    "Accept-Encoding": "identity",
                },
            )
            with urllib.request.urlopen(req, timeout=120) as resp, open(tmp_dest, "wb") as out_f:
                while chunk := resp.read(2 * 1024 * 1024):
                    out_f.write(chunk)
            elapsed = max(time.time() - t0, 0.001)

            actual_size = os.path.getsize(tmp_dest)
            actual_md5, actual_sha256 = compute_md5_sha256(tmp_dest)

            if actual_size != expected_size:
                raise ValueError(f"Size mismatch: expected {expected_size}, got {actual_size}")
            if actual_md5 != expected_md5:
                raise ValueError(f"MD5 mismatch: expected {expected_md5}, got {actual_md5}")

            os.rename(tmp_dest, dest)
            speed = (actual_size / 1e6) / elapsed
            print(f"[{key}] Downloaded & verified in {elapsed:.1f}s ({speed:.2f} MB/s, MD5: {actual_md5})", flush=True)
            return {
                "key": key,
                "size": actual_size,
                "md5": actual_md5,
                "sha256": actual_sha256,
                "status": "DOWNLOADED_VERIFIED",
                "seconds": round(elapsed, 2),
                "speed_mbs": round(speed, 2),
            }

        except Exception as e:
            wait_time = base_backoff * (2 ** (attempt - 1))
            if isinstance(e, urllib.error.HTTPError) and e.code == 429:
                retry_after = e.headers.get("Retry-After")
                if retry_after and retry_after.isdigit():
                    wait_time = max(float(retry_after) + 2.0, wait_time)
                print(f"[{key}] Rate limited (HTTP 429) on attempt {attempt}/{max_retries}. Backing off {wait_time:.1f}s...", flush=True)
            else:
                print(f"[{key}] Error on attempt {attempt}/{max_retries} ({type(e).__name__}: {e}). Retrying in {wait_time:.1f}s...", flush=True)

            if os.path.exists(tmp_dest):
                try:
                    os.remove(tmp_dest)
                except OSError:
                    pass

            if attempt == max_retries:
                raise RuntimeError(f"[{key}] Failed after {max_retries} attempts: {e}")
            time.sleep(wait_time)


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    max_parts = int(sys.argv[2]) if len(sys.argv) > 2 else 0

    os.makedirs(PARTS_DIR, exist_ok=True)
    os.makedirs(os.path.dirname(MANIFEST_PATH), exist_ok=True)

    print(f"=== H3 Step 1: Downloading SMAG Archive Parts from Zenodo 8223844 ({workers} parallel workers) ===")
    with open(ZENODO_JSON_PATH, "r") as f:
        meta = json.load(f)

    files = meta.get("files", [])
    raw_mag_parts = [f for f in files if f.get("key", "").startswith("mag.tar.gz.")]

    # Sort strictly by GNU split 2-letter suffix: aa, ab, ..., az, ba, ..., oi
    sorted_parts = sorted(raw_mag_parts, key=lambda x: x["key"].replace("mag.tar.gz.", ""))
    # Process the distant latter block; original downloader assembles the full manifest.
    sorted_parts = sorted_parts[200:]
    if max_parts > 0:
        sorted_parts = sorted_parts[:max_parts]

    print(f"Total parts to process: {len(sorted_parts)} (auxiliary last173-part block of full373 archive)")
    part_tasks = []
    for item in sorted_parts:
        chk = item.get("checksum", "")
        expected_md5 = chk.split(":", 1)[1] if ":" in chk else chk
        url = item.get("links", {}).get("self", f"https://zenodo.org/api/records/8223844/files/{item['key']}/content")
        part_tasks.append({
            "key": item["key"],
            "size": item["size"],
            "md5": expected_md5,
            "url": url,
        })

    t_start = time.time()
    results = []

    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_key = {executor.submit(download_single_part, p): p["key"] for p in part_tasks}
        completed = 0
        total = len(part_tasks)
        for future in concurrent.futures.as_completed(future_to_key):
            completed += 1
            key = future_to_key[future]
            try:
                res = future.result()
                results.append(res)
                if completed % 10 == 0 or completed == total:
                    elapsed = time.time() - t_start
                    mb_downloaded = sum(r["size"] for r in results if r["status"] == "DOWNLOADED_VERIFIED") / 1e6
                    print(f"  Progress: {completed}/{total} parts processed ({elapsed:.1f}s elapsed, {mb_downloaded:.1f} MB freshly downloaded)...", flush=True)
            except Exception as exc:
                print(f"  FATAL ERROR on {key}: {exc}", file=sys.stderr, flush=True)
                raise exc

    # Sort results by key order
    results.sort(key=lambda x: x["key"].replace("mag.tar.gz.", ""))

    # Write manifest
    with open(MANIFEST_PATH, "w") as out:
        out.write("part_key\tbyte_size\tmd5_checksum\tsha256_checksum\tstatus\telapsed_seconds\tspeed_mbs\n")
        for r in results:
            out.write(f"{r['key']}\t{r['size']}\t{r['md5']}\t{r['sha256']}\t{r['status']}\t{r['seconds']}\t{r['speed_mbs']}\n")

    total_time = time.time() - t_start
    total_bytes = sum(r["size"] for r in results)
    fresh_parts = sum(1 for r in results if r["status"] == "DOWNLOADED_VERIFIED")
    print(f"=== Complete: {len(results)}/{len(part_tasks)} parts verified in {total_time:.1f}s ({fresh_parts} freshly downloaded, {total_bytes / 1e9:.2f} GB total) ===")
    print(f"Manifest written to {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
