#!/usr/bin/env python3
"""
Step 3: Track F1 Environmental Comparator Screening (SMAG & GEM Pools)
- Screens same-caller SMAG (marine) and GEM (global terrestrial/host) reference pools.
- Uses identical 3 core CPF profiles (PF00875, PF03441, PF04244) and audited gates:
  domain iE-value <= 1e-5, profile coverage >= 0.70.
- Emits candidate tables and FASTAs for both comparative environmental cohorts.
"""

import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path

ROOT = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12")
if not ROOT.exists():
    ROOT = Path("/home/leechuck/Documents/papers/empty-quarter-mags/reproducibility/darkmatter/family_discovery_20261015/f12")

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
H3_DIR = W / "06_darkmatter/p7_analyses_clean/h3_references"
OUT_DIR = ROOT / "results/full_v1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PROFILES = ["PF00875", "PF03441", "PF04244"]

COHORTS = {
    "SMAG": H3_DIR / "results/smag_pool_canonical.faa",
    "GEM": H3_DIR / "results/gem_harmonized_pool_canonical_seed20261014.faa"
}

def run_hmm_search(cohort_name, faa_path, hmm_path):
    """Run hmmsearch on a comparator pool."""
    domtbl = OUT_DIR / f"{cohort_name.lower()}_domains.tsv"
    report = OUT_DIR / f"{cohort_name.lower()}_hmm_report.txt"
    cpus = os.environ.get("SLURM_CPUS_PER_TASK", "16")
    
    cmd = [
        "hmmsearch",
        "--cpu", cpus,
        "--domE", "1e-3",
        "-E", "1e-3",
        "--domtblout", str(domtbl),
        "-o", str(report),
        str(hmm_path),
        str(faa_path)
    ]
    print(f"Running hmmsearch on {cohort_name}: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
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
    return hits

def emit_candidates(cohort_name, hits, faa_path):
    """Extract candidate FASTA and TSV for comparator pool."""
    out_faa = OUT_DIR / f"{cohort_name.lower()}_cpf_candidates.faa"
    out_tsv = OUT_DIR / f"{cohort_name.lower()}_cpf_candidates.tsv"
    
    records = []
    current_header = None
    current_seq = []
    total_scanned = 0
    
    with open(out_faa, "w") as out:
        def process_entry(h, s_list):
            pid = h.split()[0]
            if pid in hits:
                seq_str = "".join(s_list)
                rec = hits[pid].copy()
                rec["cohort"] = cohort_name
                rec["length"] = len(seq_str)
                rec["sha256"] = hashlib.sha256(seq_str.encode()).hexdigest()
                rec["header"] = h
                records.append(rec)
                out.write(f">{h}\n{seq_str}\n")
                
        with open(faa_path) as f:
            for line in f:
                if line.startswith(">"):
                    total_scanned += 1
                    if current_header:
                        process_entry(current_header, current_seq)
                    current_header = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_header:
                process_entry(current_header, current_seq)
                
    assert len(records) == len(hits), f"Mismatch in {cohort_name}: {len(hits)} hits vs {len(records)} records"
    
    with open(out_tsv, "w") as f:
        fields = ["id", "cohort", "profile", "profile_name", "domain_ievalue", "domain_bitscore",
                  "profile_coverage", "hmm_from", "hmm_to", "hmm_len", "length", "sha256", "header"]
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(records)
        
    print(f"[OK] {cohort_name}: Scanned {total_scanned} proteins, identified {len(records)} CPF candidates.")
    return {
        "cohort": cohort_name,
        "source_file": str(faa_path),
        "total_proteins_scanned": total_scanned,
        "cpf_candidates_identified": len(records),
        "candidates_faa": str(out_faa),
        "candidates_tsv": str(out_tsv)
    }

def main():
    print("=== Step 3: Track F1 Environmental Comparator Screening ===")
    summary_file = OUT_DIR / "environmental_comparators_summary.json"
    if summary_file.exists():
        print(f"[OK] Environmental comparators summary {summary_file} already exists. Skipping recomputation.")
        return
    hmm_path = OUT_DIR / "cpf_core.hmm"
    assert hmm_path.exists(), "cpf_core.hmm must be prepared by Step 2 first"
    
    summaries = {}
    for cname, cpath in COHORTS.items():
        assert cpath.exists(), f"Missing comparator pool {cpath}"
        domtbl = run_hmm_search(cname, cpath, hmm_path)
        hits = parse_hits(domtbl)
        summary = emit_candidates(cname, hits, cpath)
        summaries[cname] = summary
        
    summary_file = OUT_DIR / "environmental_comparators_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summaries, f, indent=2)
    print(f"[SUCCESS] Environmental comparators processed. Summary -> {summary_file}")

if __name__ == "__main__":
    main()
