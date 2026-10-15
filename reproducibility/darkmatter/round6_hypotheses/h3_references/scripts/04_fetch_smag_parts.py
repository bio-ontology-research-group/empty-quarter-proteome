#!/usr/bin/env python3
"""
scripts/04_fetch_smag_parts.py
Download the first ordered parts of the public SMAG archive from Zenodo (Record 8223844)
with strict MD5 verification against official Zenodo API metadata.

Usage:
  python3 scripts/04_fetch_smag_parts.py [n_parts]
"""

import os
import sys
import hashlib
import urllib.request
import time

# Official checksums from https://zenodo.org/api/records/8223844
# (retained in reviews/darkmatter-2026-10-13/round6/public_source_checks/zenodo_8223844.json)
ZENODO_PARTS = [
    {
        "key": "mag.tar.gz.aa",
        "size": 104857600,
        "md5": "5fb8a52b97940d1d1c14aecda57ac5bc",
        "url": "https://zenodo.org/api/records/8223844/files/mag.tar.gz.aa/content"
    },
    {
        "key": "mag.tar.gz.ab",
        "size": 104857600,
        "md5": "de5f1ecb60b8ceac42b616e2780b2d04",
        "url": "https://zenodo.org/api/records/8223844/files/mag.tar.gz.ab/content"
    }
]

def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()

def main():
    n_parts = int(sys.argv[1]) if len(sys.argv) > 1 else len(ZENODO_PARTS)
    parts_to_fetch = ZENODO_PARTS[:n_parts]

    h3_dir = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references"
    parts_dir = os.path.join(h3_dir, "downloads/smag_parts")
    out_dir = os.path.join(h3_dir, "results")
    os.makedirs(parts_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    print(f"=== Fetching {len(parts_to_fetch)} public SMAG archive parts from Zenodo ===")
    manifest_rows = []

    for item in parts_to_fetch:
        key = item["key"]
        expected_size = item["size"]
        expected_md5 = item["md5"]
        url = item["url"]
        dest = os.path.join(parts_dir, key)

        # Check if already present and valid
        if os.path.exists(dest) and os.path.getsize(dest) == expected_size:
            actual_md5 = md5sum(dest)
            if actual_md5 == expected_md5:
                print(f"[{key}] Already downloaded and verified (MD5: {actual_md5}).")
                manifest_rows.append((key, expected_size, actual_md5, "verified_existing", url))
                continue

        print(f"[{key}] Downloading from {url} ({expected_size / 1e6:.1f} MB)...")
        part_tmp = dest + ".tmp"
        if os.path.exists(part_tmp):
            os.remove(part_tmp)

        t0 = time.time()
        # Download with stream
        req = urllib.request.Request(url, headers={"User-Agent": "Antigravity-H3-Research/1.0"})
        with urllib.request.urlopen(req) as resp, open(part_tmp, "wb") as out_f:
            while chunk := resp.read(1024 * 1024):
                out_f.write(chunk)
        elapsed = time.time() - t0

        actual_size = os.path.getsize(part_tmp)
        actual_md5 = md5sum(part_tmp)

        if actual_size != expected_size:
            raise ValueError(f"[{key}] Size mismatch: expected {expected_size}, got {actual_size}")
        if actual_md5 != expected_md5:
            raise ValueError(f"[{key}] MD5 mismatch: expected {expected_md5}, got {actual_md5}")

        os.rename(part_tmp, dest)
        speed = (actual_size / 1e6) / elapsed if elapsed > 0 else 0
        print(f"[{key}] Downloaded and verified in {elapsed:.1f}s ({speed:.2f} MB/s). MD5: {actual_md5}")
        manifest_rows.append((key, actual_size, actual_md5, "downloaded_verified", url))

    # Write download manifest
    manifest_path = os.path.join(out_dir, "smag_zenodo_download_manifest.tsv")
    with open(manifest_path, "w") as out:
        out.write("part_key\tbyte_size\tmd5_checksum\tstatus\tsource_url\n")
        for row in manifest_rows:
            out.write("\t".join(str(x) for x in row) + "\n")
    print(f"Wrote download manifest to {manifest_path}")
    print("DONE_04_FETCH_SMAG_PARTS")

if __name__ == "__main__":
    main()
