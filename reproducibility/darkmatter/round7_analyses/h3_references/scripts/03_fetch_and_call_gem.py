#!/usr/bin/env python3
"""
scripts/03_fetch_and_call_gem.py (Candidate Repaired Version)

GEM Reference Cohort Processing:
1. Deterministically sample GEM genomes across 3 seeds (20261014, 20261015, 20261016).
2. For primary seed 20261014:
   - Download nucleotide contigs (.fna.gz) from public NERSC repository (https://portal.nersc.gov/GEM/genomes/fna/)
   - Call genes with Pyrodigal v3.7.1 (-p meta --no-stop-codon -f gff) to achieve TRUE CALLER HARMONIZATION with SMAG.
   - Retain per-genome raw FAA and gene coordinate GFF files with coordinate mappings in genecalls/gem/.
   - Also parse published GEM Prodigal v2.6.3 proteomes (.faa.gz) for caller-sensitivity comparison.
3. For sensitivity seeds 20261015 & 20261016:
   - Parse published proteomes to evaluate Monte Carlo reference-draw sensitivity.
4. Canonicalize (strip trailing '*', uppercase, audit internal '*' and invalid chars), deduplicate,
   and emit complete file-level manifests and FASTA pools.

Repair Specifications:
- Gene-call retention:
    * preserves per-genome raw protein translations: genecalls/gem/{gid}.raw.faa
    * preserves per-genome gene coordinate outputs:  genecalls/gem/{gid}.gff
    * records original-ID-to-canonical-ID mapping:   genecalls/gem/{gid}.id_map.tsv
      linking matched canonical protein IDs to original contigs, start/end coordinates, and strand.
    * caller flags remain biologically identical: -p meta --no-stop-codon.
    * preserves both original raw outputs and canonical deduplicated pools.
    * records paths and SHA256 hashes in results/gem_harmonized_manifest_seed20261014.tsv
      and emits results/gem_harmonized_protein_coordinate_mapping_seed20261014.tsv.
- Fail-closed checks on required input paths and executables.

Usage:
    python3 scripts/03_fetch_and_call_gem.py [target_unique_count] [threads]
"""

import os
import sys
import gzip
import time
import random
import hashlib
import tempfile
import subprocess
import urllib.request
import concurrent.futures
from collections import defaultdict

GEM_ALL_IDS_PATH = os.environ.get("GEM_ALL_IDS_PATH", "results/gem_catalogue_52515_ids.txt")
GEM_FAA_DIR = os.environ.get("GEM_FAA_DIR", "/ibex/scratch/projects/c2014/rob/metagenomics/databases/gem/faa/faa")
GEM_FNA_DIR = os.environ.get("GEM_FNA_DIR", "genomes/gem_nucleotide")
GENECALLS_DIR = os.environ.get("GENECALLS_DIR", "genecalls/gem")
RESULTS_DIR = os.environ.get("RESULTS_DIR", "results")
PYRODIGAL_BIN = os.environ.get("PYRODIGAL_BIN", "/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal")

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def parse_prodigal_header(header_line):
    """
    Parse Prodigal/Pyrodigal FAA header line:
    e.g. '>2004178001.a:gws2_d1_0448_1 # 1 # 81 # -1 # ID=1_1;partial=10;...'
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
        return gid, dest, os.path.getsize(dest), sha256_file(dest), "EXISTING"

    url = f"https://portal.nersc.gov/GEM/genomes/fna/{gid}.fna.gz"
    tmp_dest = dest + f".tmp_{os.getpid()}"
    for attempt in range(1, 6):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Antigravity-H3-Research/2.0"})
            with urllib.request.urlopen(req, timeout=60) as resp, open(tmp_dest, "wb") as out_f:
                while chunk := resp.read(1024 * 1024):
                    out_f.write(chunk)
                out_f.flush()
                os.fsync(out_f.fileno())
            sz = os.path.getsize(tmp_dest)
            sh = sha256_file(tmp_dest)
            os.replace(tmp_dest, dest)
            return gid, dest, sz, sh, "DOWNLOADED"
        except Exception as e:
            if os.path.exists(tmp_dest):
                try: os.remove(tmp_dest)
                except OSError: pass
            if attempt == 5:
                raise RuntimeError(f"Failed to fetch {gid} from NERSC: {e}")
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
            print(f"Pyrodigal failed on GEM {gid} (code {res.returncode}):\nSTDERR:\n{res.stderr}\nSTDOUT:\n{res.stdout}", file=sys.stderr)
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


def load_published_gem_proteome(gid):
    """Load published Prodigal v2.6.3 proteome for a single GEM genome."""
    faa_path = os.path.join(GEM_FAA_DIR, f"{gid}.faa.gz")
    if not os.path.exists(faa_path):
        raise FileNotFoundError(f"Published GEM proteome missing: '{faa_path}'.")

    raw_prots = 0
    raw_residues = 0
    trailing_stars = 0
    internal_stars = 0
    lowercase_chars = 0
    invalid_aa_chars = 0
    canon_prots = 0
    canon_residues = 0
    records = []
    within_genome_hashes = set()

    cur_id = None
    cur_header = ""
    cur_lines = []

    def flush_prot(sid, header_line, lines):
        nonlocal raw_prots, raw_residues, trailing_stars, internal_stars
        nonlocal lowercase_chars, invalid_aa_chars, canon_prots, canon_residues

        raw_seq = "".join(lines)
        raw_prots += 1
        raw_residues += len(raw_seq)

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
        records.append((sid, canon_seq, seq_sha256))

    with gzip.open(faa_path, "rt", encoding="utf-8") as gz:
        for line in gz:
            if line.startswith(">"):
                if cur_id is not None:
                    flush_prot(cur_id, cur_header, cur_lines)
                cur_header = line.strip()
                cur_id = line[1:].strip().split()[0]
                cur_lines = []
            else:
                cur_lines.append(line.strip())
        if cur_id is not None:
            flush_prot(cur_id, cur_header, cur_lines)

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
        "faa_gz_bytes": os.path.getsize(faa_path),
        "faa_gz_sha256": sha256_file(faa_path),
    }
    return gid, stats, records


def main():
    target_unique = int(sys.argv[1]) if len(sys.argv) > 1 else 4500000
    threads = int(sys.argv[2]) if len(sys.argv) > 2 else 16

    os.makedirs(GEM_FNA_DIR, exist_ok=True)
    os.makedirs(GENECALLS_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print(f"=== H3 Step 3: Processing GEM Comparator Cohorts (Target unique proteins: {target_unique:,}, Threads: {threads}) ===")

    # Fail closed on missing inputs
    if not os.path.exists(GEM_ALL_IDS_PATH):
        raise FileNotFoundError(f"Missing GEM all-IDs file: '{GEM_ALL_IDS_PATH}'.")
    if not os.path.exists(PYRODIGAL_BIN):
        raise FileNotFoundError(f"Pyrodigal binary not found at '{PYRODIGAL_BIN}'.")
    if not os.path.exists(GEM_FAA_DIR):
        raise FileNotFoundError(f"GEM published FAA directory not found at '{GEM_FAA_DIR}'.")

    with open(GEM_ALL_IDS_PATH, "r", encoding="utf-8") as f:
        all_gem_ids = sorted(line.strip() for line in f if line.strip())
    print(f"Total available GEM genomes: {len(all_gem_ids):,}")

    seeds = [20261014, 20261015, 20261016]
    sampled_cohorts = {}

    for seed in seeds:
        rng = random.Random(seed)
        shuffled = list(all_gem_ids)
        rng.shuffle(shuffled)

        # Select a fixed genome count using the planning estimate of 2,300 proteins per genome.
        print(f"\n--- Sampling GEM candidate draw for seed {seed} ---")
        selected_gids = []
        for gid in shuffled:
            selected_gids.append(gid)
            if len(selected_gids) * 2300 >= target_unique * 1.05:
                break
        print(f"Seed {seed}: Selected {len(selected_gids):,} GEM genomes.")
        sampled_cohorts[seed] = selected_gids

    # For primary seed 20261014: Fetch nucleotide assemblies and call Pyrodigal v3.7.1
    primary_gids = sampled_cohorts[20261014]
    print(f"\n=== Primary Cohort (Seed 20261014): Fetching {len(primary_gids)} GEM Nucleotide Assemblies from NERSC ===")
    t_fetch = time.time()
    fna_paths = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(threads, 12)) as ex:
        futures = [ex.submit(fetch_gem_nucleotide, gid) for gid in primary_gids]
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            done += 1
            gid, path, sz, sh, st = fut.result()
            fna_paths[gid] = path
            if done % 200 == 0 or done == len(primary_gids):
                print(f"  Fetched {done}/{len(primary_gids)} GEM assemblies ({time.time() - t_fetch:.1f}s)...", flush=True)

    print(f"=== Primary Cohort: Running Pyrodigal v3.7.1 on {len(primary_gids)} GEM Assemblies (Retaining raw FAA & GFF) ===")
    t_call = time.time()
    call_tasks = [(gid, fna_paths[gid], PYRODIGAL_BIN, GENECALLS_DIR) for gid in primary_gids]
    harmonized_results = {}
    all_harmonized_maps = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(call_pyrodigal_on_gem, t) for t in call_tasks]
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            done += 1
            gid, stats, recs, map_recs = fut.result()
            harmonized_results[gid] = (stats, recs)
            all_harmonized_maps.extend(map_recs)
            if done % 200 == 0 or done == len(primary_gids):
                print(f"  Called genes on {done}/{len(primary_gids)} GEM genomes ({time.time() - t_call:.1f}s)...", flush=True)

    # Write harmonized FASTA pools & manifest
    h_canon_faa = os.path.join(RESULTS_DIR, "gem_harmonized_pool_canonical_seed20261014.faa")
    h_uniq_faa = os.path.join(RESULTS_DIR, "gem_harmonized_pool_unique_seed20261014.faa")
    h_man_path = os.path.join(RESULTS_DIR, "gem_harmonized_manifest_seed20261014.tsv")
    h_map_path = os.path.join(RESULTS_DIR, "gem_harmonized_protein_coordinate_mapping_seed20261014.tsv")

    h_seen_hashes = {}
    h_canon_cnt = 0
    h_canon_res = 0
    h_uniq_cnt = 0
    h_uniq_res = 0
    h_man_rows = []

    tmp_h_canon = h_canon_faa + f".tmp_{os.getpid()}"
    tmp_h_uniq = h_uniq_faa + f".tmp_{os.getpid()}"

    with open(tmp_h_canon, "w", encoding="utf-8") as fc, open(tmp_h_uniq, "w", encoding="utf-8") as fu:
        for gid in sorted(primary_gids):
            st, recs = harmonized_results[gid]
            recs.sort(key=lambda x: x[0])
            for gene_idx, sid, cseq, sh in recs:
                fc.write(f">{sid}\n{cseq}\n")
                h_canon_cnt += 1
                h_canon_res += len(cseq)
                if sh not in h_seen_hashes:
                    h_seen_hashes[sh] = (sid, cseq)
                    fu.write(f">{sid}\n{cseq}\n")
                    h_uniq_cnt += 1
                    h_uniq_res += len(cseq)
            h_man_rows.append((
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

    os.replace(tmp_h_canon, h_canon_faa)
    os.replace(tmp_h_uniq, h_uniq_faa)

    # Write harmonized manifest
    tmp_h_man = h_man_path + f".tmp_{os.getpid()}"
    with open(tmp_h_man, "w", encoding="utf-8") as out:
        out.write("genome_id\traw_protein_count\traw_residue_count\ttrailing_star_count\tinternal_star_count\t"
                  "lowercase_char_count\tinvalid_aa_char_count\tcanonical_protein_count\tcanonical_residue_count\t"
                  "unique_canonical_proteins\traw_faa_sha256\traw_faa_path\tgff_path\tgff_sha256\t"
                  "id_map_path\tid_map_sha256\n")
        for r in h_man_rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp_h_man, h_man_path)

    # Write aggregated harmonized mapping
    tmp_h_map = h_map_path + f".tmp_{os.getpid()}"
    with open(tmp_h_map, "w", encoding="utf-8") as out_m:
        out_m.write(
            "genome_id\tgene_idx\tcanonical_seq_id\toriginal_protein_id\tcontig_id\tstart_coord\tend_coord\t"
            "strand\tpartial\tlength_aa\tcanonical_seq_sha256\tstatus\traw_header\n"
        )
        for mr in all_harmonized_maps:
            out_m.write(
                f"{mr['genome_id']}\t{mr['gene_idx']}\t{mr['canonical_seq_id']}\t{mr['original_protein_id']}\t"
                f"{mr['contig_id']}\t{mr['start_coord']}\t{mr['end_coord']}\t{mr['strand']}\t{mr['partial']}\t"
                f"{mr['length_aa']}\t{mr['canonical_seq_sha256']}\t{mr['status']}\t{mr['raw_header']}\n"
            )
    os.replace(tmp_h_map, h_map_path)

    print(f"=== Harmonized GEM Pool (Seed 20261014) Summary ===")
    print(f"  Genomes: {len(primary_gids):,}")
    print(f"  Canonical proteins: {h_canon_cnt:,} ({h_canon_res:,} aa)")
    print(f"  Unique canonical proteins: {h_uniq_cnt:,} ({h_uniq_res:,} aa)")
    print(f"  Internal duplicates: {h_canon_cnt - h_uniq_cnt:,} ({100.0 * (h_canon_cnt - h_uniq_cnt) / max(h_canon_cnt, 1):.2f}% redundancy)")
    print(f"  Retained genecalls in: {GENECALLS_DIR}/ ({len(primary_gids):,} raw FAA, GFF, and ID-map files)")

    # Parse published Prodigal proteomes for seeds 20261014, 20261015, 20261016
    for seed in seeds:
        gids = sampled_cohorts[seed]
        print(f"\n=== Processing Published Prodigal Proteomes for Seed {seed} ({len(gids)} genomes) ===")
        t_pub = time.time()
        pub_results = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
            futures = [ex.submit(load_published_gem_proteome, gid) for gid in gids]
            done = 0
            for fut in concurrent.futures.as_completed(futures):
                done += 1
                gid, stats, recs = fut.result()
                pub_results[gid] = (stats, recs)
                if done % 400 == 0 or done == len(gids):
                    print(f"  Parsed {done}/{len(gids)} published proteomes ({time.time() - t_pub:.1f}s)...", flush=True)

        p_canon_faa = os.path.join(RESULTS_DIR, f"gem_published_pool_canonical_seed{seed}.faa")
        p_uniq_faa = os.path.join(RESULTS_DIR, f"gem_published_pool_unique_seed{seed}.faa")
        p_man_path = os.path.join(RESULTS_DIR, f"gem_published_manifest_seed{seed}.tsv")

        p_seen_hashes = {}
        p_canon_cnt = 0
        p_canon_res = 0
        p_uniq_cnt = 0
        p_uniq_res = 0
        p_man_rows = []

        tmp_p_canon = p_canon_faa + f".tmp_{os.getpid()}"
        tmp_p_uniq = p_uniq_faa + f".tmp_{os.getpid()}"

        with open(tmp_p_canon, "w", encoding="utf-8") as fc, open(tmp_p_uniq, "w", encoding="utf-8") as fu:
            for gid in sorted(gids):
                st, recs = pub_results[gid]
                recs.sort(key=lambda x: x[0])
                for sid, cseq, sh in recs:
                    fc.write(f">{sid}\n{cseq}\n")
                    p_canon_cnt += 1
                    p_canon_res += len(cseq)
                    if sh not in p_seen_hashes:
                        p_seen_hashes[sh] = (sid, cseq)
                        fu.write(f">{sid}\n{cseq}\n")
                        p_uniq_cnt += 1
                        p_uniq_res += len(cseq)
                p_man_rows.append((
                    gid,
                    st["faa_gz_bytes"],
                    st["faa_gz_sha256"],
                    st["raw_protein_count"],
                    st["raw_residue_count"],
                    st["trailing_star_count"],
                    st["internal_star_count"],
                    st["lowercase_char_count"],
                    st["invalid_aa_char_count"],
                    st["canonical_protein_count"],
                    st["canonical_residue_count"],
                    st["unique_canonical_proteins"],
                ))

        os.replace(tmp_p_canon, p_canon_faa)
        os.replace(tmp_p_uniq, p_uniq_faa)

        tmp_p_man = p_man_path + f".tmp_{os.getpid()}"
        with open(tmp_p_man, "w", encoding="utf-8") as out:
            out.write("genome_id\tfaa_gz_bytes\tfaa_gz_sha256\traw_protein_count\traw_residue_count\t"
                      "trailing_star_count\tinternal_star_count\tlowercase_char_count\tinvalid_aa_char_count\t"
                      "canonical_protein_count\tcanonical_residue_count\tunique_canonical_proteins\n")
            for r in p_man_rows:
                out.write("\t".join(str(x) for x in r) + "\n")
        os.replace(tmp_p_man, p_man_path)

        print(f"=== Published GEM Pool (Seed {seed}) Summary ===")
        print(f"  Genomes: {len(gids):,}")
        print(f"  Canonical proteins: {p_canon_cnt:,} ({p_canon_res:,} aa)")
        print(f"  Unique canonical proteins: {p_uniq_cnt:,} ({p_uniq_res:,} aa)")
        print(f"  Internal duplicates: {p_canon_cnt - p_uniq_cnt:,} ({100.0 * (p_canon_cnt - p_uniq_cnt) / max(p_canon_cnt, 1):.2f}% redundancy)")


if __name__ == "__main__":
    main()
