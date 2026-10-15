#!/usr/bin/env python3
"""
scripts/05_extract_and_genecall_smag.py
Extract complete genome members from the downloaded Zenodo SMAG parts,
guarding strictly against truncated trailing archive members.
Then call genes with Pyrodigal v3.7.1 (-p meta --no-stop-codon) in parallel,
matching the project's exact Empty Quarter and SMAG gene-calling convention.

Usage:
  python3 scripts/05_extract_and_genecall_smag.py [threads]
"""

import os
import sys
import glob
import subprocess
import time
import concurrent.futures

def process_single_genome(gpath, pyrodigal_bin):
    """Run Pyrodigal on a single genome and return manifest entry + fasta chunk."""
    gname = os.path.basename(gpath)
    mag_id = os.path.splitext(gname)[0]
    tmp_trans = gpath + ".trans.faa"
    
    cmd_call = [
        pyrodigal_bin,
        "-i", gpath,
        "-a", tmp_trans,
        "-p", "meta",
        "--no-stop-codon",
        "-o", "/dev/null"
    ]
    subprocess.run(cmd_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    records = []
    g_prots = 0
    g_residues = 0
    current_gene_idx = 0
    if os.path.exists(tmp_trans):
        with open(tmp_trans, "r") as tf:
            for line in tf:
                if line.startswith(">"):
                    current_gene_idx += 1
                    g_prots += 1
                    records.append(f">{mag_id}::{mag_id}_{current_gene_idx}\n")
                else:
                    clean_seq = line.strip()
                    g_residues += len(clean_seq)
                    records.append(clean_seq + "\n")
        os.remove(tmp_trans)

    manifest_entry = (mag_id, os.path.getsize(gpath), g_prots, g_residues)
    return manifest_entry, "".join(records)

def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    h3_dir = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references"
    parts_dir = os.path.join(h3_dir, "downloads/smag_parts")
    tmp_genomes_dir = os.path.join(h3_dir, "tmp/smag_pilot_genomes")
    out_dir = os.path.join(h3_dir, "results")
    os.makedirs(tmp_genomes_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    pyrodigal_bin = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"

    # Step 1: Find downloaded parts in order
    part_files = sorted(glob.glob(os.path.join(parts_dir, "mag.tar.gz.*")))
    print(f"=== Extracting from {len(part_files)} SMAG parts: {[os.path.basename(p) for p in part_files]} ===")
    assert len(part_files) > 0, "No parts found to extract!"

    # Step 2: Extract genomes using pipe, capturing any trailing truncated file
    # If already extracted, skip re-extracting
    extracted_files = []
    for root, _, files in os.walk(tmp_genomes_dir):
        for f in files:
            f_lower = f.lower()
            if any(f_lower.endswith(ext) for ext in [".fa", ".fasta", ".fna"]):
                extracted_files.append(os.path.join(root, f))

    if not extracted_files:
        print(f"Extracting tar stream into {tmp_genomes_dir}...")
        t0 = time.time()
        cmd_cat = f"cat {' '.join(part_files)} | tar -xz -C {tmp_genomes_dir} 2>/dev/null || true"
        subprocess.run(cmd_cat, shell=True, check=False)
        print(f"Extraction stream finished in {time.time() - t0:.1f}s.")

        for root, _, files in os.walk(tmp_genomes_dir):
            for f in files:
                f_lower = f.lower()
                if any(f_lower.endswith(ext) for ext in [".fa", ".fasta", ".fna"]):
                    extracted_files.append(os.path.join(root, f))

    extracted_files.sort()
    print(f"Found {len(extracted_files)} candidate genome files.")

    # Prune any invalid or truncated file (e.g. the last file if cut off mid-record)
    valid_genomes = []
    pruned_count = 0
    for fpath in extracted_files:
        try:
            size = os.path.getsize(fpath)
            if size < 1000:
                os.remove(fpath)
                pruned_count += 1
                continue
            with open(fpath, "r", errors="ignore") as f:
                first_char = f.read(1)
                if first_char != ">":
                    os.remove(fpath)
                    pruned_count += 1
                    continue
                f.seek(max(0, size - 200))
                tail = f.read()
                if not tail.endswith("\n"):
                    print(f"Pruning truncated file: {os.path.basename(fpath)}")
                    os.remove(fpath)
                    pruned_count += 1
                    continue
            valid_genomes.append(fpath)
        except Exception:
            if os.path.exists(fpath):
                os.remove(fpath)
            pruned_count += 1

    print(f"Validated complete genomes: {len(valid_genomes)} (pruned {pruned_count} incomplete/trailing).")
    assert len(valid_genomes) > 0, "No valid complete genomes found!"

    # Step 3: Run Pyrodigal gene calling in parallel with ProcessPoolExecutor
    print(f"=== Calling genes with Pyrodigal v3.7.1 on {len(valid_genomes)} genomes using {threads} parallel workers ===")
    t_call0 = time.time()
    recalled_faa_path = os.path.join(out_dir, "smag_pilot_recalled.faa")
    manifest_rows = []
    total_proteins = 0
    total_residues = 0

    with open(recalled_faa_path, "w") as out_faa:
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            future_to_genome = {executor.submit(process_single_genome, gpath, pyrodigal_bin): gpath for gpath in valid_genomes}
            completed_count = 0
            for future in concurrent.futures.as_completed(future_to_genome):
                completed_count += 1
                m_entry, fasta_str = future.result()
                out_faa.write(fasta_str)
                manifest_rows.append(m_entry)
                total_proteins += m_entry[2]
                total_residues += m_entry[3]

                if completed_count % 25 == 0 or completed_count == len(valid_genomes):
                    elapsed = time.time() - t_call0
                    rate = completed_count / elapsed if elapsed > 0 else 0
                    print(f"  Processed {completed_count}/{len(valid_genomes)} MAGs ({total_proteins} proteins, {rate:.1f} MAGs/sec)...", flush=True)

    elapsed_call = time.time() - t_call0
    print(f"Gene calling complete in {elapsed_call:.1f}s.")
    print(f"Total SMAG pilot proteins: {total_proteins} across {len(valid_genomes)} MAGs ({total_residues} residues).")

    # Step 4: Write genome manifest (sorted by mag_id)
    manifest_rows.sort(key=lambda x: x[0])
    manifest_path = os.path.join(out_dir, "smag_pilot_genome_manifest.tsv")
    with open(manifest_path, "w") as out_man:
        out_man.write("mag_id\tnucleotide_bytes\tprotein_count\tresidue_count\n")
        for row in manifest_rows:
            out_man.write("\t".join(str(x) for x in row) + "\n")
    print(f"Wrote SMAG genome manifest to {manifest_path}")

    # Clean up uncompressed nucleotide genomes to save disk
    print(f"Cleaning up {len(valid_genomes)} uncompressed nucleotide FASTA files...")
    for gpath in valid_genomes:
        if os.path.exists(gpath):
            os.remove(gpath)
    print("DONE_05_EXTRACT_AND_GENECALL_SMAG")

if __name__ == "__main__":
    main()
