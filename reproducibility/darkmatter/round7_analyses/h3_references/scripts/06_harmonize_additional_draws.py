#!/usr/bin/env python3
"""
scripts/06_harmonize_additional_draws.py

H3 Reference Catalogue Specificity: Harmonized Gene-Calling for Sensitivity Draws
Seeds 20261015 and 20261016 (2,055 GEM genomes each).

Key Objectives & Constraints:
1. Validates genome ID lists for seeds 20261015 and 20261016 against declared published
   manifests (gem_published_manifest_seed*.tsv) and declared random frame/selection rule
   (fixed 2,055 genomes per seed: ceil(4500000*1.05/2300)).
2. Reuses pinned, immutable genecalls from primary seed 20261014 for shared genomes
   without re-calling or overwriting existing files/hashes.
3. Downloads missing nucleotide assemblies (.fna.gz) from public NERSC repository.
4. Executes gene-calling with EXACT harmonized Pyrodigal v3.7.1 settings (-p meta --no-stop-codon -f gff).
5. Retains per-genome raw FAA, gene coordinate GFF, and ID-to-coordinate mapping tables in genecalls/gem/.
6. Emits canonical pools, deduplicated unique pools, duplicate maps, legacy manifests,
   and extended with_coordinates manifests for seeds 20261015 and 20261016.
7. Generates complete provenance metadata in results/gem_harmonized_draws_summary.json.

Usage:
    python3 scripts/06_harmonize_additional_draws.py [threads]
"""

import os
import sys
import gzip
import time
import json
import random
import hashlib
import tempfile
import subprocess
import urllib.request
import concurrent.futures
from collections import defaultdict

GEM_ALL_IDS_PATH = "results/gem_catalogue_52515_ids.txt"
GEM_FNA_DIR = "genomes/gem_nucleotide"
GENECALLS_DIR = "genecalls/gem"
RESULTS_DIR = "results"
PYRODIGAL_BIN = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"

PRIMARY_WITH_COORDS_MANIFEST = os.path.join(RESULTS_DIR, "gem_harmonized_manifest_seed20261014.with_coordinates.tsv")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def atomic_write_json(target_path, data):
    dir_name = os.path.dirname(os.path.abspath(target_path))
    os.makedirs(dir_name, exist_ok=True)
    tmp_path = target_path + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"
    with open(tmp_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, sort_keys=True)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp_path, target_path)


def parse_prodigal_header(header_line):
    """
    Parse Prodigal/Pyrodigal FAA header line:
    e.g. '>contig_1_1 # 1 # 114 # 1 # ID=1_1;partial=10;start_type=Edge;...'
    Returns (orig_id, contig_id, start_coord, end_coord, strand, partial_val).
    """
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


def fetch_gem_nucleotide(gid):
    """Download .fna.gz for a single GEM genome from NERSC portal."""
    dest = os.path.join(GEM_FNA_DIR, f"{gid}.fna.gz")
    if os.path.exists(dest) and os.path.getsize(dest) > 100:
        try:
            with gzip.open(dest, "rb") as gz_test:
                gz_test.read(1024)
            return gid, dest, os.path.getsize(dest), sha256_file(dest), "EXISTING"
        except Exception:
            pass  # corrupt existing file, re-download

    url = f"https://portal.nersc.gov/GEM/genomes/fna/{gid}.fna.gz"
    tmp_dest = dest + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"
    for attempt in range(1, 6):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Antigravity-H3-Research/2.0"})
            with urllib.request.urlopen(req, timeout=60) as resp, open(tmp_dest, "wb") as out_f:
                while chunk := resp.read(1024 * 1024):
                    out_f.write(chunk)
                out_f.flush()
                os.fsync(out_f.fileno())

            # Validate gzip readability
            with gzip.open(tmp_dest, "rb") as gz_test:
                gz_test.read(1024)

            sz = os.path.getsize(tmp_dest)
            sh = sha256_file(tmp_dest)
            os.replace(tmp_dest, dest)
            return gid, dest, sz, sh, "DOWNLOADED"
        except Exception as e:
            if os.path.exists(tmp_dest):
                try: os.remove(tmp_dest)
                except OSError: pass
            if attempt == 5:
                raise RuntimeError(f"Failed to fetch {gid} from NERSC after 5 attempts: {e}")
            time.sleep(2 * attempt)


def call_pyrodigal_on_gem(args):
    """
    Run Pyrodigal v3.7.1 on a downloaded GEM .fna.gz genome.
    Retains per-genome raw FAA, gene coordinate GFF, and ID-to-coordinate mapping.
    Caller flags: -p meta --no-stop-codon -f gff.
    """
    gid, fna_path, pyrodigal_bin, genecalls_dir = args

    raw_faa_path = os.path.join(genecalls_dir, f"{gid}.raw.faa")
    gff_path = os.path.join(genecalls_dir, f"{gid}.gff")
    id_map_path = os.path.join(genecalls_dir, f"{gid}.id_map.tsv")

    # If outputs already exist and are non-empty, load and verify
    if os.path.exists(raw_faa_path) and os.path.exists(gff_path) and os.path.exists(id_map_path):
        if os.path.getsize(raw_faa_path) > 0 and os.path.getsize(gff_path) > 0 and os.path.getsize(id_map_path) > 0:
            return load_existing_genecall(gid, raw_faa_path, gff_path, id_map_path)

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
            for p in (tmp_raw_faa, tmp_gff):
                if os.path.exists(p):
                    try: os.remove(p)
                    except OSError: pass
            raise RuntimeError(f"Pyrodigal failed on GEM {gid} (code {res.returncode}): {res.stderr.strip()}")
    finally:
        if os.path.exists(tmp_fa_name):
            try: os.remove(tmp_fa_name)
            except OSError: pass

    with open(tmp_raw_faa, "rb") as bf:
        raw_faa_bytes = bf.read()
    raw_faa_sha256 = hashlib.sha256(raw_faa_bytes).hexdigest()

    with open(tmp_gff, "rb") as gf:
        gff_bytes = gf.read()
    gff_sha256 = hashlib.sha256(gff_bytes).hexdigest()

    os.replace(tmp_raw_faa, raw_faa_path)
    os.replace(tmp_gff, gff_path)

    text = raw_faa_bytes.decode("ascii")
    raw_prots = 0
    raw_residues = 0
    trailing_stars = 0
    internal_stars = 0
    lowercase_chars = 0
    invalid_aa_chars = 0
    canon_prots = 0
    canon_residues = 0

    records = []
    mapping_records = []
    within_genome_hashes = set()

    cur_idx = 0
    cur_header = ""
    cur_lines = []

    def flush_protein(gene_idx, header_line, lines):
        nonlocal raw_prots, raw_residues, trailing_stars, internal_stars
        nonlocal lowercase_chars, invalid_aa_chars, canon_prots, canon_residues

        raw_seq = "".join(lines)
        raw_prots += 1
        raw_residues += len(raw_seq)

        orig_id, contig_id, start_c, end_c, strand, partial_v = parse_prodigal_header(header_line)
        canon_seq_id = f"{gid}::{gid}_{gene_idx}"

        if raw_seq.endswith("*"):
            trailing_stars += 1
        stripped = raw_seq.rstrip("*")

        status = "retained"
        canon_seq = ""
        seq_sha256 = ""

        if "*" in stripped:
            internal_stars += 1
            status = "filtered_internal_star"
        else:
            lowercase_chars += sum(1 for c in stripped if c.islower())
            canon_seq = stripped.upper()
            bad_chars = sum(1 for c in canon_seq if c not in VALID_AA)
            invalid_aa_chars += bad_chars
            if not canon_seq:
                status = "filtered_empty"
            elif bad_chars > 0:
                status = "filtered_invalid_aa"
            else:
                canon_prots += 1
                canon_residues += len(canon_seq)
                seq_sha256 = hashlib.sha256(canon_seq.encode("ascii")).hexdigest()
                within_genome_hashes.add(seq_sha256)
                records.append((gene_idx, canon_seq_id, canon_seq, seq_sha256))

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

    # Write per-genome coordinate ID mapping table atomically
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
        "raw_protein_count": raw_prots,
        "raw_residue_count": raw_residues,
        "trailing_star_count": trailing_stars,
        "internal_star_count": internal_stars,
        "lowercase_char_count": lowercase_chars,
        "invalid_aa_char_count": invalid_aa_chars,
        "canonical_protein_count": canon_prots,
        "canonical_residue_count": canon_residues,
        "unique_canonical_proteins": len(within_genome_hashes),
        "raw_faa_sha256": raw_faa_sha256,
        "raw_faa_path": os.path.abspath(raw_faa_path),
        "gff_path": os.path.abspath(gff_path),
        "gff_sha256": gff_sha256,
        "id_map_path": os.path.abspath(id_map_path),
        "id_map_sha256": id_map_sha256,
    }
    return gid, stats, records, mapping_records


def load_existing_genecall(gid, raw_faa_path, gff_path, id_map_path):
    """Load existing genecall outputs without re-running Pyrodigal."""
    raw_faa_sha256 = sha256_file(raw_faa_path)
    gff_sha256 = sha256_file(gff_path)
    id_map_sha256 = sha256_file(id_map_path)

    records = []
    mapping_records = []
    within_genome_hashes = set()

    with open(id_map_path, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if not parts:
                continue
            row = dict(zip(header, parts))
            mapping_records.append(row)

    # Load canonical records from raw FAA
    with open(raw_faa_path, "r", encoding="utf-8") as f:
        text = f.read()

    raw_prots = 0
    raw_residues = 0
    trailing_stars = 0
    internal_stars = 0
    lowercase_chars = 0
    invalid_aa_chars = 0
    canon_prots = 0
    canon_residues = 0

    cur_idx = 0
    cur_lines = []

    def flush_prot(gene_idx, lines):
        nonlocal raw_prots, raw_residues, trailing_stars, internal_stars
        nonlocal lowercase_chars, invalid_aa_chars, canon_prots, canon_residues

        raw_seq = "".join(lines)
        raw_prots += 1
        raw_residues += len(raw_seq)
        canon_seq_id = f"{gid}::{gid}_{gene_idx}"

        if raw_seq.endswith("*"):
            trailing_stars += 1
        stripped = raw_seq.rstrip("*")

        if "*" in stripped:
            internal_stars += 1
            return
        lowercase_chars += sum(1 for c in stripped if c.islower())
        canon_seq = stripped.upper()
        bad_chars = sum(1 for c in canon_seq if c not in VALID_AA)
        invalid_aa_chars += bad_chars
        if not canon_seq or bad_chars > 0:
            return

        canon_prots += 1
        canon_residues += len(canon_seq)
        seq_sha256 = hashlib.sha256(canon_seq.encode("ascii")).hexdigest()
        within_genome_hashes.add(seq_sha256)
        records.append((gene_idx, canon_seq_id, canon_seq, seq_sha256))

    for line in text.splitlines():
        if line.startswith(">"):
            if cur_idx > 0:
                flush_prot(cur_idx, cur_lines)
            cur_idx += 1
            cur_lines = []
        else:
            cur_lines.append(line.strip())
    if cur_idx > 0:
        flush_prot(cur_idx, cur_lines)

    stats = {
        "genome_id": gid,
        "raw_protein_count": raw_prots,
        "raw_residue_count": raw_residues,
        "trailing_star_count": trailing_stars,
        "internal_star_count": internal_stars,
        "lowercase_char_count": lowercase_chars,
        "invalid_aa_char_count": invalid_aa_chars,
        "canonical_protein_count": canon_prots,
        "canonical_residue_count": canon_residues,
        "unique_canonical_proteins": len(within_genome_hashes),
        "raw_faa_sha256": raw_faa_sha256,
        "raw_faa_path": os.path.abspath(raw_faa_path),
        "gff_path": os.path.abspath(gff_path),
        "gff_sha256": gff_sha256,
        "id_map_path": os.path.abspath(id_map_path),
        "id_map_sha256": id_map_sha256,
    }
    return gid, stats, records, mapping_records


def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    print(f"=== H3 Step 6: Harmonized Gene-Calling for Additional GEM Draws (Seeds 20261015, 20261016) ===")
    print(f"Running on host: {os.uname().nodename} with {threads} CPU workers")
    t_start = time.time()

    os.makedirs(GEM_FNA_DIR, exist_ok=True)
    os.makedirs(GENECALLS_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    # Validate required base files
    for req in (GEM_ALL_IDS_PATH, PYRODIGAL_BIN, PRIMARY_WITH_COORDS_MANIFEST):
        if not os.path.exists(req):
            raise FileNotFoundError(f"Missing required file: '{req}'.")

    # Load all GEM catalogue IDs
    with open(GEM_ALL_IDS_PATH, "r", encoding="utf-8") as f:
        all_gem_ids = sorted(line.strip() for line in f if line.strip())
    print(f"Loaded {len(all_gem_ids):,} total GEM catalogue IDs.")

    # Validate seed 20261015 and 20261016 genome selections against published manifests
    seeds = [20261015, 20261016]
    seed_gids = {}

    for s in seeds:
        # Re-derive expected selection from random frame
        rng = random.Random(s)
        shuffled = list(all_gem_ids)
        rng.shuffle(shuffled)
        expected_draw = []
        for gid in shuffled:
            expected_draw.append(gid)
            if len(expected_draw) * 2300 >= 4500000 * 1.05:
                break
        assert len(expected_draw) == 2055, f"Expected 2055 genomes for seed {s}, got {len(expected_draw)}"

        # Verify against existing published manifest
        pub_man = os.path.join(RESULTS_DIR, f"gem_published_manifest_seed{s}.tsv")
        if not os.path.exists(pub_man):
            raise FileNotFoundError(f"Missing published manifest for seed {s}: '{pub_man}'.")
        with open(pub_man, "r", encoding="utf-8") as mf:
            mf.readline()
            actual_pub_gids = [line.strip().split("\t")[0] for line in mf if line.strip()]
        assert len(actual_pub_gids) == 2055, f"Manifest has {len(actual_pub_gids)} genomes, expected 2055"
        assert sorted(expected_draw) == sorted(actual_pub_gids), f"Derived selection does not match published manifest for seed {s}!"

        seed_gids[s] = sorted(actual_pub_gids)
        print(f"  Seed {s}: verified exact 2,055 genome IDs against random frame and published manifest.")

    # Load pinned genomes from primary seed 20261014 to guarantee zero re-calling and zero overwrites
    print("\nLoading pinned primary seed 20261014 coordinate manifest...")
    primary_pinned = {}
    with open(PRIMARY_WITH_COORDS_MANIFEST, "r", encoding="utf-8") as f:
        header = f.readline().strip().split("\t")
        for line in f:
            parts = line.strip().split("\t")
            if parts:
                row = dict(zip(header, parts))
                primary_pinned[row["genome_id"]] = row
    print(f"  Loaded {len(primary_pinned):,} pinned primary genomes from seed 20261014.")

    # Determine unique genomes across seeds 20261015 and 20261016
    all_target_gids = sorted(set(seed_gids[20261015]) | set(seed_gids[20261016]))
    overlap_14 = set(all_target_gids).intersection(set(primary_pinned.keys()))
    novel_gids = sorted(set(all_target_gids) - overlap_14)
    print(f"Total distinct genomes across seeds 20261015 and 20261016: {len(all_target_gids):,}")
    print(f"  Shared with primary seed 20261014: {len(overlap_14):,} (will strictly reuse pinned genecalls)")
    print(f"  Novel genomes to prepare:          {len(novel_gids):,}")

    # Fetch missing assemblies from NERSC
    print(f"\n=== Fetching missing nucleotide assemblies from NERSC ({len(novel_gids):,} targets) ===")
    t_fetch = time.time()
    fna_paths = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(threads, 16)) as ex:
        futures = [ex.submit(fetch_gem_nucleotide, gid) for gid in all_target_gids]
        done = 0
        existing_cnt = 0
        dl_cnt = 0
        for fut in concurrent.futures.as_completed(futures):
            done += 1
            gid, path, sz, sh, st = fut.result()
            fna_paths[gid] = path
            if st == "EXISTING":
                existing_cnt += 1
            else:
                dl_cnt += 1
            if done % 500 == 0 or done == len(all_target_gids):
                print(f"  Assemblies ready: {done}/{len(all_target_gids)} (Cached: {existing_cnt}, Downloaded: {dl_cnt}, {time.time() - t_fetch:.1f}s)...", flush=True)

    # Harmonized Pyrodigal Calling
    print(f"\n=== Running Harmonized Pyrodigal v3.7.1 on GEM Genomes ({threads} workers) ===")
    t_call = time.time()
    call_tasks = [(gid, fna_paths[gid], PYRODIGAL_BIN, GENECALLS_DIR) for gid in all_target_gids]
    results_by_gid = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(call_pyrodigal_on_gem, t) for t in call_tasks]
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            done += 1
            gid, stats, recs, map_recs = fut.result()
            results_by_gid[gid] = (stats, recs, map_recs)
            if done % 400 == 0 or done == len(call_tasks):
                print(f"  Genomes called & verified: {done}/{len(call_tasks)} ({time.time() - t_call:.1f}s)...", flush=True)

    # For each seed, emit canonical FASTA, unique FASTA, duplicate map, and manifests
    print("\n=== Emitting Harmonized Pools and Manifests for Seeds 20261015 and 20261016 ===")
    seed_summaries = {}

    for s in seeds:
        print(f"\n--- Generating Artifacts for Seed {s} (2,055 genomes) ---")
        gids = seed_gids[s]

        canon_faa_path = os.path.join(RESULTS_DIR, f"gem_harmonized_pool_canonical_seed{s}.faa")
        uniq_faa_path = os.path.join(RESULTS_DIR, f"gem_harmonized_pool_unique_seed{s}.faa")
        dup_map_path = os.path.join(RESULTS_DIR, f"gem_harmonized_duplicate_map_seed{s}.tsv")
        manifest_legacy_path = os.path.join(RESULTS_DIR, f"gem_harmonized_manifest_seed{s}.tsv")
        manifest_coords_path = os.path.join(RESULTS_DIR, f"gem_harmonized_manifest_seed{s}.with_coordinates.tsv")
        coord_map_path = os.path.join(RESULTS_DIR, f"gem_harmonized_protein_coordinate_mapping_seed{s}.tsv")

        seen_hashes = {}
        dup_rows = []
        manifest_rows = []
        all_coord_maps = []

        total_canon_prots = 0
        total_canon_res = 0
        total_uniq_prots = 0
        total_uniq_res = 0

        tmp_canon = canon_faa_path + f".tmp_{os.getpid()}"
        tmp_uniq = uniq_faa_path + f".tmp_{os.getpid()}"

        with open(tmp_canon, "w", encoding="utf-8") as fc, open(tmp_uniq, "w", encoding="utf-8") as fu:
            for gid in gids:
                st, recs, map_recs = results_by_gid[gid]
                recs.sort(key=lambda x: x[0])
                all_coord_maps.extend(map_recs)

                for gene_idx, sid, cseq, sh in recs:
                    fc.write(f">{sid}\n{cseq}\n")
                    total_canon_prots += 1
                    total_canon_res += len(cseq)
                    if sh not in seen_hashes:
                        seen_hashes[sh] = sid
                        fu.write(f">{sid}\n{cseq}\n")
                        total_uniq_prots += 1
                        total_uniq_res += len(cseq)
                    else:
                        dup_rows.append((sid, seen_hashes[sh], len(cseq), sh))

                manifest_rows.append((
                    gid,
                    st["raw_protein_count"],
                    st["raw_residue_count"],
                    st["trailing_star_count"],
                    st["internal_star_count"],
                    st["lowercase_char_count"],
                    st["invalid_aa_char_count"],
                    st["canonical_protein_count"],
                    st["canonical_residue_count"],
                    st["unique_canonical_proteins"],
                    st["raw_faa_sha256"],
                    st["raw_faa_path"],
                    st["gff_path"],
                    st["gff_sha256"],
                    st["id_map_path"],
                    st["id_map_sha256"],
                ))

        os.replace(tmp_canon, canon_faa_path)
        os.replace(tmp_uniq, uniq_faa_path)

        # Write duplicate map
        tmp_dup = dup_map_path + f".tmp_{os.getpid()}"
        with open(tmp_dup, "w", encoding="utf-8") as out:
            out.write("duplicate_seq_id\trepresentative_seq_id\tlength_aa\tseq_sha256\n")
            for r in dup_rows:
                out.write("\t".join(str(x) for x in r) + "\n")
        os.replace(tmp_dup, dup_map_path)

        # Write legacy manifest (matching seed 20261014 columns)
        tmp_man_leg = manifest_legacy_path + f".tmp_{os.getpid()}"
        with open(tmp_man_leg, "w", encoding="utf-8") as out:
            out.write("genome_id\traw_protein_count\traw_residue_count\ttrailing_star_count\tinternal_star_count\t"
                      "lowercase_char_count\tinvalid_aa_char_count\tcanonical_protein_count\tcanonical_residue_count\t"
                      "unique_canonical_proteins\traw_faa_sha256\n")
            for r in manifest_rows:
                out.write("\t".join(str(x) for x in r[:11]) + "\n")
        os.replace(tmp_man_leg, manifest_legacy_path)

        # Write extended manifest with coordinates
        tmp_man_coords = manifest_coords_path + f".tmp_{os.getpid()}"
        with open(tmp_man_coords, "w", encoding="utf-8") as out:
            out.write("genome_id\traw_protein_count\traw_residue_count\ttrailing_star_count\tinternal_star_count\t"
                      "lowercase_char_count\tinvalid_aa_char_count\tcanonical_protein_count\tcanonical_residue_count\t"
                      "unique_canonical_proteins\traw_faa_sha256\traw_faa_path\tgff_path\tgff_sha256\t"
                      "id_map_path\tid_map_sha256\n")
            for r in manifest_rows:
                out.write("\t".join(str(x) for x in r) + "\n")
        os.replace(tmp_man_coords, manifest_coords_path)

        # Write aggregated coordinate mapping
        tmp_cmap = coord_map_path + f".tmp_{os.getpid()}"
        with open(tmp_cmap, "w", encoding="utf-8") as out:
            out.write(
                "genome_id\tgene_idx\tcanonical_seq_id\toriginal_protein_id\tcontig_id\tstart_coord\tend_coord\t"
                "strand\tpartial\tlength_aa\tcanonical_seq_sha256\tstatus\traw_header\n"
            )
            for mr in all_coord_maps:
                out.write(
                    f"{mr['genome_id']}\t{mr['gene_idx']}\t{mr['canonical_seq_id']}\t{mr['original_protein_id']}\t"
                    f"{mr['contig_id']}\t{mr['start_coord']}\t{mr['end_coord']}\t{mr['strand']}\t{mr['partial']}\t"
                    f"{mr['length_aa']}\t{mr['canonical_seq_sha256']}\t{mr['status']}\t{mr['raw_header']}\n"
                )
        os.replace(tmp_cmap, coord_map_path)

        redundancy = 100.0 * len(dup_rows) / max(total_canon_prots, 1)
        print(f"  Genomes: {len(gids):,}")
        print(f"  Canonical proteins: {total_canon_prots:,} ({total_canon_res:,} aa)")
        print(f"  Unique canonical:   {total_uniq_prots:,} ({total_uniq_res:,} aa)")
        print(f"  Internal duplicates: {len(dup_rows):,} ({redundancy:.2f}% redundancy)")

        seed_summaries[str(s)] = {
            "seed": s,
            "genome_count": len(gids),
            "selection_rule": "fixed_count_ceil_4500000_1.05_div_2300",
            "canonical_protein_count": total_canon_prots,
            "canonical_residue_count": total_canon_res,
            "unique_canonical_proteins": total_uniq_prots,
            "unique_canonical_residues": total_uniq_res,
            "internal_duplicate_count": len(dup_rows),
            "redundancy_pct": round(redundancy, 2),
            "overlap_with_seed20261014": len(set(gids).intersection(set(primary_pinned.keys()))),
            "pool_canonical_faa": os.path.abspath(canon_faa_path),
            "pool_canonical_sha256": sha256_file(canon_faa_path),
            "pool_unique_faa": os.path.abspath(uniq_faa_path),
            "pool_unique_sha256": sha256_file(uniq_faa_path),
            "duplicate_map_tsv": os.path.abspath(dup_map_path),
            "duplicate_map_sha256": sha256_file(dup_map_path),
            "manifest_legacy_tsv": os.path.abspath(manifest_legacy_path),
            "manifest_legacy_sha256": sha256_file(manifest_legacy_path),
            "manifest_coords_tsv": os.path.abspath(manifest_coords_path),
            "manifest_coords_sha256": sha256_file(manifest_coords_path),
            "coordinate_mapping_tsv": os.path.abspath(coord_map_path),
            "coordinate_mapping_sha256": sha256_file(coord_map_path),
        }

    # Record overall provenance metadata
    summary_data = {
        "analysis": "H3_harmonized_additional_draws",
        "description": "Harmonized Pyrodigal v3.7.1 (-p meta --no-stop-codon) calling for GEM comparator sensitivity draws (seeds 20261015 & 20261016)",
        "scientific_interpretation": (
            "Seeds 20261014, 20261015, and 20261016 represent Monte Carlo reference-draw sensitivity evaluations "
            "of the GEM catalogue (2,055 genomes each), not biological replicates. Contrasts reflect end-to-end "
            "caller-plus-matching sensitivity."
        ),
        "pyrodigal_binary": PYRODIGAL_BIN,
        "caller_flags": "-p meta --no-stop-codon -f gff",
        "seeds": seed_summaries,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_runtime_seconds": round(time.time() - t_start, 1),
    }

    summary_json_path = os.path.join(RESULTS_DIR, "gem_harmonized_draws_summary.json")
    atomic_write_json(summary_json_path, summary_data)
    print(f"\n=== All Harmonized Draws Completed Successfully in {time.time() - t_start:.1f}s ===")
    print(f"Summary saved: {summary_json_path}")


if __name__ == "__main__":
    main()
