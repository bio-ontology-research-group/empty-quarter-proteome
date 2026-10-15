#!/usr/bin/env python3
"""
Step 2: Track F1 Full Census Across All 5,438 ANI95 Representative Proteomes
- Reads all 5,438 species representatives from Wdb.csv.
- Concatenates the 5,438 proteomes from drep_rep_proteomes/{mag}.faa.gz.
- Runs profile-first search using the 3 core CPF profiles (PF00875, PF03441, PF04244).
- Applies audited gates: domain iE-value <= 1e-5, profile coverage >= 0.70.
- Emits eq_cpf_candidates.tsv and eq_cpf_candidates.faa.
"""

import csv
import gzip
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12")
if not ROOT.exists():
    ROOT = Path("/home/leechuck/Documents/papers/empty-quarter-mags/reproducibility/darkmatter/family_discovery_20261015/f12")

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
OUT_DIR = ROOT / "results/full_v1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PROFILES = ["PF00875", "PF03441", "PF04244"]

def prepare_core_hmm():
    """Concatenate and press the 3 core CPF profiles."""
    hmm_out = OUT_DIR / "cpf_core.hmm"
    with open(hmm_out, "wb") as out:
        for acc in PROFILES:
            src = ROOT / f"references_v2/{acc}.hmm"
            assert src.exists(), f"Missing profile {src}"
            b = src.read_bytes()
            assert f"ACC   {acc}".encode() in b
            out.write(b)
    # Run hmmpress
    subprocess.run(["hmmpress", "-f", str(hmm_out)], check=True)
    return hmm_out

def build_cohort_faa():
    """Concatenate all 5,438 representative proteomes."""
    wdb_path = W / "05_refine/_drep/out/data_tables/Wdb.csv"
    rows = list(csv.DictReader(open(wdb_path)))
    assert len(rows) == 5438 and len({r["cluster"] for r in rows}) == 5438, f"Invalid Wdb count {len(rows)}"
    
    cohort_faa = OUT_DIR / "all_5438_reps.faa"
    cohort_manifest = OUT_DIR / "all_5438_manifest.tsv"
    
    print(f"Concatenating {len(rows)} representative proteomes into {cohort_faa}...")
    manifest = []
    total_proteins = 0
    
    with open(cohort_faa, "w") as out:
        for r in rows:
            mag = r["genome"].removesuffix(".fa")
            p = W / f"06_darkmatter/drep_rep_proteomes/{mag}.faa.gz"
            assert p.exists(), f"Proteome {p} missing"
            ng = 0
            with gzip.open(p, "rt") as f:
                for line in f:
                    if line.startswith(">"):
                        ng += 1
                        assert line[1:].startswith(mag + "::"), f"Invalid header in {p}: {line[:40]}"
                    out.write(line)
            manifest.append({
                "mag": mag,
                "cluster": r["cluster"],
                "path": str(p),
                "proteins": ng
            })
            total_proteins += ng
            if len(manifest) % 1000 == 0:
                print(f"  Processed {len(manifest)} / {len(rows)} proteomes ({total_proteins} proteins)...")
                
    with open(cohort_manifest, "w") as f:
        w = csv.DictWriter(f, fieldnames=["mag", "cluster", "path", "proteins"], delimiter="\t")
        w.writeheader()
        w.writerows(manifest)
        
    print(f"[OK] Concatenated {len(rows)} proteomes: {total_proteins} total proteins.")
    return cohort_faa, total_proteins

def run_hmmsearch(cohort_faa, hmm_path):
    """Run hmmsearch against the full proteome catalog."""
    domtbl = OUT_DIR / "eq_domains.tsv"
    report = OUT_DIR / "eq_hmm_report.txt"
    cpus = os.environ.get("SLURM_CPUS_PER_TASK", "16")
    
    cmd = [
        "hmmsearch",
        "--cpu", cpus,
        "--domE", "1e-3",
        "-E", "1e-3",
        "--domtblout", str(domtbl),
        "-o", str(report),
        str(hmm_path),
        str(cohort_faa)
    ]
    print(f"Running hmmsearch: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"[OK] hmmsearch complete. Parsing {domtbl}...")
    return domtbl

def parse_hits(domtbl_path):
    """Parse domain hits requiring dom_iE <= 1e-5 and coverage >= 0.70."""
    hits = {}
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 22:
                continue
            acc = parts[4].split(".")[0]
            if acc not in PROFILES:
                continue
            hmm_len = int(parts[5])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            cov = (hmm_to - hmm_from + 1) / hmm_len
            dom_ie = float(parts[12])
            dom_score = float(parts[13])
            
            if dom_ie <= 1e-5 and cov >= 0.70:
                rec = {
                    "id": parts[0],
                    "profile": acc,
                    "profile_name": parts[3],
                    "domain_ievalue": dom_ie,
                    "domain_bitscore": dom_score,
                    "profile_coverage": cov,
                    "hmm_from": hmm_from,
                    "hmm_to": hmm_to,
                    "hmm_len": hmm_len
                }
                pid = parts[0]
                if pid not in hits or dom_score > hits[pid]["domain_bitscore"]:
                    hits[pid] = rec
    print(f"[OK] Identified {len(hits)} unique proteins passing profile gates.")
    return hits

def emit_candidates(hits, cohort_faa):
    """Extract FASTA sequences and metadata for all hits."""
    candidates_faa = OUT_DIR / "eq_cpf_candidates.faa"
    candidates_tsv = OUT_DIR / "eq_cpf_candidates.tsv"
    
    records = []
    current_header = None
    current_seq = []
    
    with open(candidates_faa, "w") as out:
        def process_entry(h, s_list):
            pid = h.split()[0]
            if pid in hits:
                seq_str = "".join(s_list)
                rec = hits[pid].copy()
                rec["mag"] = pid.split("::")[0]
                rec["length"] = len(seq_str)
                rec["sha256"] = hashlib.sha256(seq_str.encode()).hexdigest()
                rec["partial"] = "unknown"
                rec["header"] = h
                for item in h.split(";"):
                    if "partial=" in item:
                        rec["partial"] = item.split("partial=")[-1].strip()
                records.append(rec)
                out.write(f">{h}\n{seq_str}\n")
                
        with open(cohort_faa) as f:
            for line in f:
                if line.startswith(">"):
                    if current_header:
                        process_entry(current_header, current_seq)
                    current_header = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_header:
                process_entry(current_header, current_seq)
                
    assert len(records) == len(hits), f"Mismatch between parsed hits ({len(hits)}) and emitted records ({len(records)})"
    
    with open(candidates_tsv, "w") as f:
        fields = ["id", "profile", "profile_name", "domain_ievalue", "domain_bitscore", 
                  "profile_coverage", "hmm_from", "hmm_to", "hmm_len", "mag", "length", 
                  "sha256", "partial", "header"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(records)
        
    complete_cds = sum(1 for r in records if r["partial"] == "00")
    carrier_mags = len({r["mag"] for r in records})
    summary = {
        "catalogue_representatives": 5438,
        "total_proteins_scanned": None,
        "profile_supported_cpf_candidates": len(records),
        "complete_cds_partial00": complete_cds,
        "carrier_genomes": carrier_mags,
        "profiles_used": PROFILES,
        "gates": {"domain_ievalue": 1e-5, "profile_coverage": 0.70}
    }
    with open(OUT_DIR / "eq_cpf_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
        
    print(f"[SUCCESS] Emitted {len(records)} CPF candidates ({complete_cds} complete CDS) across {carrier_mags} representative MAGs.")
    return records, summary

def main():
    print("=== Step 2: Track F1 Full Census Across All 5,438 Proteomes ===")
    if (OUT_DIR / "eq_cpf_summary.json").exists() and (OUT_DIR / "eq_cpf_candidates.tsv").exists() and (OUT_DIR / "eq_cpf_candidates.faa").exists():
        print("[OK] Stage 2 full census candidates and summary already exist. Skipping recomputation.")
        return
    hmm_path = prepare_core_hmm()
    cohort_faa, total_proteins = build_cohort_faa()
    domtbl_path = run_hmmsearch(cohort_faa, hmm_path)
    hits = parse_hits(domtbl_path)
    records, summary = emit_candidates(hits, cohort_faa)
    summary["total_proteins_scanned"] = total_proteins
    with open(OUT_DIR / "eq_cpf_summary.json", "w") as f:
        json.dump(summary, f, indent=2)
    print("=== Step 2 Complete ===")

if __name__ == "__main__":
    main()
