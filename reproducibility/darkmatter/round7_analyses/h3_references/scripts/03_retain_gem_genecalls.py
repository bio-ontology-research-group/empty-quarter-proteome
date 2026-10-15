#!/usr/bin/env python3
"""
scripts/03_retain_gem_genecalls.py (Candidate Standalone Rerun Tool)

MINIMAL RETENTION-ONLY RERUN FOR EXECUTED GEM COHORT:
This script performs a targeted retention-only pass over the existing, already-downloaded
primary GEM assemblies (seed 20261014) to emit and preserve:
1. Raw protein translations: genecalls/gem/{gid}.raw.faa
2. Gene coordinates (GFF):   genecalls/gem/{gid}.gff
3. Coordinate ID mappings:   genecalls/gem/{gid}.id_map.tsv
4. Verified manifest with paths and hashes: results/gem_harmonized_manifest_seed20261014.tsv

Key Operational Constraints:
- NO downloads: reuses existing local genomes/gem_nucleotide/{gid}.fna.gz files.
- NO search or matching reruns: downstream DIAMOND databases and searches are not invalidated.
- NO published proteome re-parsing: skips seeds 20261015/20261016 published draws.
- Strict provenance verification:
    * Verifies that newly called raw FAA SHA256 matches existing raw_faa_sha256 in manifest.
    * Verifies that newly called canonical protein sequences and hashes match existing
      results/gem_harmonized_pool_canonical_seed20261014.faa byte-for-byte.
- Fails closed on any hash or count discrepancy.

STATUS: CANDIDATE REPAIR SCRIPT (DO NOT EXECUTE UNTIL ROOT REVIEW).
"""

import os
import sys
import gzip
import time
import hashlib
import tempfile
import subprocess
import concurrent.futures

GEM_FNA_DIR = os.environ.get("GEM_FNA_DIR", "genomes/gem_nucleotide")
GENECALLS_DIR = os.environ.get("GENECALLS_DIR", "genecalls/gem")
RESULTS_DIR = os.environ.get("RESULTS_DIR", "results")
PYRODIGAL_BIN = os.environ.get("PYRODIGAL_BIN", "/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal")

EXISTING_CANON_FAA = os.path.join(RESULTS_DIR, "gem_harmonized_pool_canonical_seed20261014.faa")
EXISTING_MANIFEST = os.path.join(RESULTS_DIR, "gem_harmonized_manifest_seed20261014.tsv")

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def parse_prodigal_header(header_line):
    hdr = header_line.lstrip(">").strip()
    tokens = hdr.split()
    orig_id = tokens[0] if tokens else ""
    start_coord = ""
    end_coord = ""
    strand = ""
    partial_val = ""

    hash_parts = [p.strip() for p in hdr.split("#")]
    if len(hash_parts) >= 4:
        start_coord = hash_parts[1]
        end_coord = hash_parts[2]
        strand_num = hash_parts[3]
        strand = "+" if strand_num == "1" else ("-" if strand_num == "-1" else strand_num)
        if len(hash_parts) >= 5:
            attrs = hash_parts[4]
            for attr in attrs.split(";"):
                if attr.strip().startswith("partial="):
                    partial_val = attr.strip().split("=", 1)[1]

    if "_" in orig_id:
        contig_id = orig_id.rsplit("_", 1)[0]
    else:
        contig_id = orig_id

    return orig_id, contig_id, start_coord, end_coord, strand, partial_val


def process_single_gem_retention(args):
    gid, fna_path, pyrodigal_bin, genecalls_dir, expected_raw_sha256 = args

    raw_faa_path = os.path.join(genecalls_dir, f"{gid}.raw.faa")
    gff_path = os.path.join(genecalls_dir, f"{gid}.gff")
    id_map_path = os.path.join(genecalls_dir, f"{gid}.id_map.tsv")

    tmp_raw_faa = raw_faa_path + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"
    tmp_gff = gff_path + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"
    tmp_id_map = id_map_path + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"

    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp_fa:
        tmp_fa_name = tmp_fa.name
        with gzip.open(fna_path, "rb") as gz_in:
            while chunk := gz_in.read(1024 * 1024):
                tmp_fa.write(chunk)
        tmp_fa.flush()

    try:
        cmd = [
            pyrodigal_bin,
            "-i", tmp_fa_name,
            "-a", tmp_raw_faa,
            "-p", "meta",
            "--no-stop-codon",
            "-f", "gff",
            "-o", tmp_gff,
        ]
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"Pyrodigal failed on {gid}: {res.stderr.strip()}")
    finally:
        if os.path.exists(tmp_fa_name):
            try: os.remove(tmp_fa_name)
            except OSError: pass

    with open(tmp_raw_faa, "rb") as bf:
        raw_faa_bytes = bf.read()
    raw_faa_sha256 = hashlib.sha256(raw_faa_bytes).hexdigest()

    # Verification against existing manifest raw_faa_sha256
    if expected_raw_sha256 and raw_faa_sha256 != expected_raw_sha256:
        raise ValueError(
            f"Provenance mismatch for GEM {gid}: newly called raw FAA SHA256 ({raw_faa_sha256}) "
            f"does not match existing manifest ({expected_raw_sha256})!"
        )

    with open(tmp_gff, "rb") as gf:
        gff_bytes = gf.read()
    gff_sha256 = hashlib.sha256(gff_bytes).hexdigest()

    os.replace(tmp_raw_faa, raw_faa_path)
    os.replace(tmp_gff, gff_path)

    text = raw_faa_bytes.decode("ascii")
    records = []
    mapping_records = []
    cur_idx = 0
    cur_header = ""
    cur_lines = []

    def flush_protein(gene_idx, header_line, lines):
        raw_seq = "".join(lines)
        orig_id, contig_id, start_c, end_c, strand, partial_v = parse_prodigal_header(header_line)
        canon_seq_id = f"{gid}::{gid}_{gene_idx}"
        stripped = raw_seq.rstrip("*")

        status = "retained"
        canon_seq = ""
        seq_sha256 = ""

        if "*" in stripped:
            status = "filtered_internal_star"
        else:
            canon_seq = stripped.upper()
            bad_chars = sum(1 for c in canon_seq if c not in VALID_AA)
            if not canon_seq:
                status = "filtered_empty"
            elif bad_chars > 0:
                status = "filtered_invalid_aa"
            else:
                seq_sha256 = hashlib.sha256(canon_seq.encode("ascii")).hexdigest()
                records.append((canon_seq_id, canon_seq, seq_sha256))

        mapping_records.append({
            "genome_id": gid,
            "gene_idx": gene_idx,
            "canonical_seq_id": canon_seq_id if status == "retained" else "",
            "original_protein_id": orig_id,
            "contig_id": contig_id,
            "start_coord": start_c,
            "end_coord": end_c,
            "strand": strand,
            "partial": partial_v,
            "length_aa": len(canon_seq) if status == "retained" else len(stripped),
            "canonical_seq_sha256": seq_sha256,
            "status": status,
            "raw_header": header_line,
        })

    for line in text.splitlines():
        if line.startswith(">"):
            if cur_idx > 0:
                flush_protein(cur_idx, cur_header, cur_lines)
            cur_idx += 1
            cur_header = line
            cur_lines = []
        else:
            cur_lines.append(line.strip())
    if cur_idx > 0:
        flush_protein(cur_idx, cur_header, cur_lines)

    # Write per-genome coordinate mapping table
    with open(tmp_id_map, "w", encoding="utf-8") as out_m:
        out_m.write(
            "genome_id\tgene_idx\tcanonical_seq_id\toriginal_protein_id\tcontig_id\tstart_coord\tend_coord\t"
            "strand\tpartial\tlength_aa\tcanonical_seq_sha256\tstatus\traw_header\n"
        )
        for mr in mapping_records:
            out_m.write(
                f"{mr['genome_id']}\t{mr['gene_idx']}\t{mr['canonical_seq_id']}\t{mr['original_protein_id']}\t"
                f"{mr['contig_id']}\t{mr['start_coord']}\t{mr['end_coord']}\t{mr['strand']}\t{mr['partial']}\t"
                f"{mr['length_aa']}\t{mr['canonical_seq_sha256']}\t{mr['status']}\t{mr['raw_header']}\n"
            )
        out_m.flush()
        os.fsync(out_m.fileno())
    os.replace(tmp_id_map, id_map_path)
    id_map_sha256 = sha256_file(id_map_path)

    stats = {
        "genome_id": gid,
        "raw_faa_sha256": raw_faa_sha256,
        "raw_faa_path": os.path.abspath(raw_faa_path),
        "gff_path": os.path.abspath(gff_path),
        "gff_sha256": gff_sha256,
        "id_map_path": os.path.abspath(id_map_path),
        "id_map_sha256": id_map_sha256,
        "canonical_protein_count": len(records),
    }
    return gid, stats, records, mapping_records


def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    print(f"=== Minimal Retention-Only Pass for Executed GEM Cohort ({threads} threads) ===")

    # Verification of input files
    for p in (EXISTING_MANIFEST, EXISTING_CANON_FAA, PYRODIGAL_BIN):
        if not os.path.exists(p):
            raise FileNotFoundError(f"Required existing artifact missing: '{p}'.")

    # Read existing manifest to get genome IDs and declared raw_faa_sha256
    manifest_info = {}
    with open(EXISTING_MANIFEST, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        gid_idx = header.index("genome_id")
        sha_idx = header.index("raw_faa_sha256")
        for line in f:
            parts = line.strip().split("\t")
            if parts and len(parts) > max(gid_idx, sha_idx):
                manifest_info[parts[gid_idx]] = parts[sha_idx]

    primary_gids = sorted(manifest_info.keys())
    print(f"Target genomes from manifest: {len(primary_gids):,}")

    os.makedirs(GENECALLS_DIR, exist_ok=True)

    tasks = []
    for gid in primary_gids:
        fna_path = os.path.join(GEM_FNA_DIR, f"{gid}.fna.gz")
        if not os.path.exists(fna_path):
            raise FileNotFoundError(f"Required downloaded assembly missing: '{fna_path}'.")
        tasks.append((gid, fna_path, PYRODIGAL_BIN, GENECALLS_DIR, manifest_info[gid]))

    print(f"Executing Pyrodigal coordinate extraction on {len(tasks):,} existing assemblies...")
    t0 = time.time()
    all_recs_by_gid = {}
    all_mapping_recs = []
    updated_stats = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(process_single_gem_retention, t) for t in tasks]
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            done += 1
            gid, st, recs, map_recs = fut.result()
            all_recs_by_gid[gid] = recs
            all_mapping_recs.extend(map_recs)
            updated_stats[gid] = st
            if done % 200 == 0 or done == len(tasks):
                print(f"  Processed & verified {done}/{len(tasks)} GEM genomes ({time.time() - t0:.1f}s)...", flush=True)

    # Verification against existing canonical pool
    print("Verifying canonical protein sequences against existing pool...")
    expected_hashes = {}
    for gid, recs in all_recs_by_gid.items():
        for sid, seq, seq_hash in recs:
            if sid in expected_hashes:
                raise ValueError(f"Duplicate retained canonical ID: {sid}")
            expected_hashes[sid] = seq_hash
    total_expected = len(expected_hashes)
    pool_seq_count = 0
    def verify_pool_record(sid, chunks):
        nonlocal pool_seq_count
        if sid not in expected_hashes:
            raise ValueError(f"Unknown or duplicate existing pool ID: {sid}")
        observed = hashlib.sha256("".join(chunks).encode("ascii")).hexdigest()
        expected = expected_hashes.pop(sid)
        if observed != expected:
            raise ValueError(f"Canonical sequence SHA256 mismatch for {sid}: {observed} != {expected}")
        pool_seq_count += 1
    with open(EXISTING_CANON_FAA, "r", encoding="utf-8") as f:
        cur_id, chunks = None, []
        for line in f:
            if line.startswith(">"):
                if cur_id is not None:
                    verify_pool_record(cur_id, chunks)
                cur_id, chunks = line[1:].strip().split()[0], []
            else:
                chunks.append(line.strip())
        if cur_id is not None:
            verify_pool_record(cur_id, chunks)
    if expected_hashes:
        raise ValueError(f"Existing pool missing {len(expected_hashes)} retained canonical IDs")
    assert pool_seq_count == total_expected
    print(f"  Verified exact IDs and sequence SHA256 for all {pool_seq_count:,} canonical proteins.", flush=True)

    # Write aggregated coordinate mapping
    agg_map_path = os.path.join(RESULTS_DIR, "gem_harmonized_protein_coordinate_mapping_seed20261014.tsv")
    tmp_agg_map = agg_map_path + f".tmp_{os.getpid()}"
    with open(tmp_agg_map, "w", encoding="utf-8") as out:
        out.write(
            "genome_id\tgene_idx\tcanonical_seq_id\toriginal_protein_id\tcontig_id\tstart_coord\tend_coord\t"
            "strand\tpartial\tlength_aa\tcanonical_seq_sha256\tstatus\traw_header\n"
        )
        for mr in all_mapping_recs:
            out.write(
                f"{mr['genome_id']}\t{mr['gene_idx']}\t{mr['canonical_seq_id']}\t{mr['original_protein_id']}\t"
                f"{mr['contig_id']}\t{mr['start_coord']}\t{mr['end_coord']}\t{mr['strand']}\t{mr['partial']}\t"
                f"{mr['length_aa']}\t{mr['canonical_seq_sha256']}\t{mr['status']}\t{mr['raw_header']}\n"
            )
    os.replace(tmp_agg_map, agg_map_path)

    # Rewrite manifest with added path and hash columns
    manifest_updated_path = os.path.join(RESULTS_DIR, "gem_harmonized_manifest_seed20261014.with_coordinates.tsv")
    tmp_man_updated = manifest_updated_path + f".tmp_{os.getpid()}"
    with open(EXISTING_MANIFEST, "r", encoding="utf-8") as in_m, open(tmp_man_updated, "w", encoding="utf-8") as out_m:
        orig_header = in_m.readline().strip()
        new_header = orig_header + "\traw_faa_path\tgff_path\tgff_sha256\tid_map_path\tid_map_sha256\n"
        out_m.write(new_header)
        for line in in_m:
            row = line.strip().split("\t")
            gid = row[0]
            st = updated_stats[gid]
            out_m.write(line.strip() + f"\t{st['raw_faa_path']}\t{st['gff_path']}\t{st['gff_sha256']}\t{st['id_map_path']}\t{st['id_map_sha256']}\n")
    os.replace(tmp_man_updated, manifest_updated_path)

    print(f"\nSUCCESS: Minimal retention pass complete!")
    print(f"  Retained files in: {GENECALLS_DIR}/")
    print(f"  Aggregated mapping: {agg_map_path}")
    print(f"  Updated manifest:   {manifest_updated_path}")


if __name__ == "__main__":
    main()
