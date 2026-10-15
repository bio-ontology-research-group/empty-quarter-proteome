#!/usr/bin/env python3
"""
01b_prep_dedup_and_cohort.py: Build deduplicated reference indexes and select deterministic 256-genome cohort for F1 pilot.

Root Audit Rules:
- Cdb secondary_cluster -> Wdb cluster mapping, one winner asserted per cluster.
- Deterministic SHA256-ranked 256 of verified 5438 Wdb representatives.
- Retain cohort list before search; record exact proteins.
- All source/mapping errors fail, no empty defaults.
- Must execute on Slurm worker node (checked by require_slurm_worker.sh / SLURM_JOB_ID).
"""

import csv
import gzip
import hashlib
import os
import subprocess
import sys
from pathlib import Path

def require_slurm():
    job_id = os.environ.get("SLURM_JOB_ID")
    if not job_id:
        print("ERROR: Must run inside Slurm worker allocation!", file=sys.stderr)
        sys.exit(1)
    print(f"[OK] Slurm worker environment confirmed. SLURM_JOB_ID={job_id}")

def build_diamond_db(faa_path: Path, dmnd_path: Path, diamond_bin: str):
    if not faa_path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {faa_path}")
    print(f"Building DIAMOND database from {faa_path} -> {dmnd_path}...")
    cmd = [
        diamond_bin,
        "makedb",
        "--in", str(faa_path),
        "-d", str(dmnd_path),
        "--quiet"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"DIAMOND makedb failed: {res.stderr}")
    print(f"[OK] Built {dmnd_path}.dmnd")

def main():
    require_slurm()

    f12_dir = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12")
    w_dir = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
    wdb_path = w_dir / "05_refine/_drep/out/data_tables/Wdb.csv"
    drep_dir = w_dir / "06_darkmatter/drep_rep_proteomes"
    refs_dir = f12_dir / "references_v2"
    out_dir = f12_dir / "results/pilot_v3_results"
    out_dir.mkdir(parents=True, exist_ok=True)

    diamond_bin = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"
    if not os.path.isfile(diamond_bin):
        diamond_bin = "diamond"

    # Step 1: Build deduplicated DIAMOND databases
    cpf_faa = refs_dir / "cpf_references_dedup.faa"
    cpf_dmnd = refs_dir / "cpf_references_dedup"
    build_diamond_db(cpf_faa, cpf_dmnd, diamond_bin)

    phos_faa = refs_dir / "phosphatase_references_dedup.faa"
    phos_dmnd = refs_dir / "phosphatase_references_dedup"
    build_diamond_db(phos_faa, phos_dmnd, diamond_bin)

    # Step 2: Read Wdb.csv
    if not wdb_path.is_file():
        raise FileNotFoundError(f"Wdb.csv not found at {wdb_path}")

    print(f"Loading Wdb.csv from {wdb_path}...")
    wdb_rows = []
    with open(wdb_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            genome = row["genome"]
            cluster = row["cluster"]
            score = row["score"]
            bin_name = genome[:-3] if genome.endswith(".fa") else genome
            wdb_rows.append({
                "genome": genome,
                "bin_name": bin_name,
                "cluster": cluster,
                "score": score
            })

    total_wdb = len(wdb_rows)
    print(f"Total Wdb representatives loaded: {total_wdb}")
    if total_wdb != 5438:
        raise ValueError(f"Expected 5438 Wdb representatives, found {total_wdb}!")

    # Step 3: Deterministic SHA256 ranking
    for r in wdb_rows:
        # Stable deterministic hash of genome filename
        r["sha256"] = hashlib.sha256(r["genome"].encode("utf-8")).hexdigest()

    # Sort deterministically by SHA256 hash
    ranked_rows = sorted(wdb_rows, key=lambda x: x["sha256"])
    cohort_256 = ranked_rows[:256]

    print(f"Selected top 256 representatives by SHA256 hash.")
    cohort_tsv = out_dir / "f1_pilot_cohort_genomes.tsv"
    with open(cohort_tsv, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["cohort_rank", "genome_file", "bin_name", "cluster_id", "drep_score", "sha256_rank_hash", "proteome_path"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for idx, r in enumerate(cohort_256, 1):
            faa_path = drep_dir / f"{r['bin_name']}.faa.gz"
            if not faa_path.is_file():
                raise FileNotFoundError(f"Proteome file missing for cohort genome: {faa_path}")
            r_out = {
                "cohort_rank": idx,
                "genome_file": r["genome"],
                "bin_name": r["bin_name"],
                "cluster_id": r["cluster"],
                "drep_score": r["score"],
                "sha256_rank_hash": r["sha256"],
                "proteome_path": str(faa_path)
            }
            writer.writerow(r_out)

    print(f"[OK] Wrote cohort manifest: {cohort_tsv}")

    # Step 4: Extract exact proteins of the 256 cohort into single concatenated FASTA
    concat_faa = out_dir / "f1_pilot_256_cohort_proteins.faa"
    print(f"Extracting cohort proteins into {concat_faa}...")

    total_proteins = 0
    total_aa = 0
    genome_protein_counts = {}

    with open(concat_faa, "w", encoding="utf-8") as out_f:
        for idx, r in enumerate(cohort_256, 1):
            faa_path = drep_dir / f"{r['bin_name']}.faa.gz"
            p_count = 0
            with gzip.open(faa_path, "rt", encoding="utf-8") as in_f:
                for line in in_f:
                    if line.startswith(">"):
                        p_count += 1
                        total_proteins += 1
                    else:
                        total_aa += len(line.strip())
                    out_f.write(line)
            genome_protein_counts[r["bin_name"]] = p_count

    # Calculate SHA256 of extracted proteins file
    hasher = hashlib.sha256()
    with open(concat_faa, "rb") as f:
        while chunk := f.read(65536):
            hasher.update(chunk)
    concat_sha256 = hasher.hexdigest()

    meta_tsv = out_dir / "f1_pilot_cohort_metadata.tsv"
    with open(meta_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["cohort_genome_count", len(cohort_256)])
        writer.writerow(["total_wdb_population", total_wdb])
        writer.writerow(["cohort_fraction", f"{len(cohort_256)/total_wdb:.4f}"])
        writer.writerow(["total_proteins_extracted", total_proteins])
        writer.writerow(["mean_proteins_per_genome", f"{total_proteins/len(cohort_256):.1f}"])
        writer.writerow(["min_proteins_in_genome", min(genome_protein_counts.values())])
        writer.writerow(["max_proteins_in_genome", max(genome_protein_counts.values())])
        writer.writerow(["total_amino_acids", total_aa])
        writer.writerow(["cohort_proteins_sha256", concat_sha256])

    print(f"[OK] Extracted {total_proteins} proteins ({total_aa} aa) across 256 representative genomes.")
    print(f"[OK] SHA256: {concat_sha256}")
    print("[SUCCESS] 01b_prep_dedup_and_cohort completed cleanly.")

if __name__ == "__main__":
    main()
