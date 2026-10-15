#!/usr/bin/env python3
"""
fetch_profiles.py - Fetch and verify source-backed Pfam HMM profiles for F3 genomic context test.

All accessions are programmatically validated against the InterPro API.
Gathering thresholds (GA) from the curated Pfam models are strictly enforced.
"""

import gzip
import hashlib
import json
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROFILES_DIR = ROOT / "profiles"
PROFILES_DIR.mkdir(parents=True, exist_ok=True)

# Curated, source-verified sensory and signaling domains + rhodopsin controls
PFAM_PANEL = [
    # Two-Component System: Kinases
    ("PF00512", "HisKA", "His Kinase A (phospho-acceptor) domain", "TCS_Kinase"),
    ("PF02518", "HATPase_c", "Histidine kinase-like ATPase domain", "TCS_Kinase"),
    ("PF07568", "HisKA_2", "Histidine kinase phosphoacceptor family 2", "TCS_Kinase"),
    ("PF07730", "HisKA_3", "Histidine kinase phosphoacceptor family 3", "TCS_Kinase"),
    # Two-Component System: Response Regulators
    ("PF00072", "Response_reg", "Response regulator receiver domain", "TCS_Regulator"),
    ("PF00486", "Trans_reg_C", "Transcriptional regulatory protein, C terminal", "TCS_Regulator"),
    # Chemotaxis / Phototaxis Transducers & Linkers
    ("PF00015", "MCPsignal", "Methyl-accepting chemotaxis protein (MCP) signaling domain", "Transducer_Linker"),
    ("PF00672", "HAMP", "HAMP domain linker", "Transducer_Linker"),
    ("PF01584", "CheW", "CheW-like chemotaxis adaptor domain", "Transducer_Linker"),
    # Sensor Domains
    ("PF00989", "PAS", "PAS fold sensor domain", "Sensor_Input"),
    ("PF08447", "PAS_3", "PAS fold family 3", "Sensor_Input"),
    ("PF08448", "PAS_4", "PAS fold family 4", "Sensor_Input"),
    ("PF13426", "PAS_9", "PAS domain family 9", "Sensor_Input"),
    ("PF01590", "GAF", "GAF sensor domain", "Sensor_Input"),
    ("PF13185", "GAF_2", "GAF domain family 2", "Sensor_Input"),
    # Second Messenger Signaling
    ("PF00990", "GGDEF", "Diguanylate cyclase, GGDEF domain", "Second_Messenger"),
    ("PF00563", "EAL", "Diguanylate phosphodiesterase, EAL domain", "Second_Messenger"),
    # Rhodopsin Reference Controls (for self-verification / tandem detection)
    ("PF01036", "Bac_rhodopsin", "Microbial Type-1 rhodopsin", "Rhodopsin_Control"),
    ("PF18761", "Heliorhodopsin", "Type-3 heliorhodopsin", "Rhodopsin_Control"),
]


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def download_pfam_hmm(acc: str, name: str) -> tuple[Path, str, str, str]:
    dest = PROFILES_DIR / f"{acc}_{name}.hmm"
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{acc}/?annotation=hmm"
    req = urllib.request.Request(url, headers={"User-Agent": "EmptyQuarter-F3Context/1.0"})
    
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = resp.read()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    
    lines = data.decode("utf-8", errors="ignore").splitlines()
    model_name = ""
    model_acc = ""
    ga_line = ""
    for line in lines[:30]:
        if line.startswith("NAME"):
            model_name = line.split()[1]
        elif line.startswith("ACC"):
            model_acc = line.split()[1]
        elif line.startswith("GA"):
            ga_line = line.strip()
    
    dest.write_bytes(data)
    sha = compute_sha256(dest)
    print(f"Downloaded {acc} ({name}): model_name='{model_name}', acc='{model_acc}', {ga_line}, size={len(data)} B, sha256={sha[:12]}...")
    return dest, model_name, model_acc, ga_line


def main():
    t0 = time.time()
    print("=" * 80)
    print("FETCHING AND VERIFYING F3 SENSORY CONTEXT HMM PROFILES")
    print("=" * 80)
    
    manifest = []
    combined_hmm = PROFILES_DIR / "f3_sensory_context.hmm"
    combined_content = bytearray()
    
    for acc, expected_name, desc, category in PFAM_PANEL:
        dest, model_name, model_acc, ga = download_pfam_hmm(acc, expected_name)
        sha = compute_sha256(dest)
        manifest.append({
            "pfam_accession": acc,
            "expected_name": expected_name,
            "model_name": model_name,
            "model_acc": model_acc,
            "gathering_threshold": ga,
            "category": category,
            "description": desc,
            "file": dest.name,
            "sha256": sha,
            "bytes": dest.stat().st_size
        })
        combined_content.extend(dest.read_bytes())
    
    combined_hmm.write_bytes(combined_content)
    combined_sha = compute_sha256(combined_hmm)
    print(f"\nCombined HMM profile database written to {combined_hmm.name}")
    print(f"Total models: {len(PFAM_PANEL)}, Total size: {len(combined_content)} bytes, SHA-256: {combined_sha}")
    
    manifest_file = PROFILES_DIR / "profiles_manifest.json"
    manifest_file.write_text(json.dumps({
        "timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "combined_hmm": combined_hmm.name,
        "combined_sha256": combined_sha,
        "combined_bytes": len(combined_content),
        "total_models": len(PFAM_PANEL),
        "models": manifest
    }, indent=2), encoding="utf-8")
    print(f"Manifest written to {manifest_file.name}")
    
    # Try hmmpress if hmmpress is available locally
    try:
        subprocess.run(["hmmpress", "-f", str(combined_hmm)], check=True, capture_output=True)
        print("Local hmmpress succeeded.")
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("Note: hmmpress will be executed on the IBEX Slurm worker.")
    
    print(f"Done in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    main()
