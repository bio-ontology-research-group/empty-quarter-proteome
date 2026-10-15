#!/usr/bin/env python3
"""
scripts/02_extract_and_call_smag.py (Candidate Repaired Version)

Extract the pre-declared random sample of 1,500 SMAG MAGs from the streamed
multi-part Zenodo archive (Record 8223844), verify individual nucleotide FASTAs,
call genes with Pyrodigal v3.7.1 (-p meta --no-stop-codon -f gff), canonicalize,
deduplicate, retain per-genome raw FAA and gene coordinate GFF files with coordinate
mappings in genecalls/smag/, and emit complete manifests and FASTA pools.

Repair Specifications:
- Full sample guard: exactly 1,500 prespecified IDs, all 373 archive parts verified
  at expected byte sizes, and exact recovered target ID set (retaining root patch).
- Checked cache semantics for existing extracted genome assemblies:
    * requires verified extraction manifest (smag_extraction_manifest.json) recording
      declared bytes, observed bytes, nucleotide SHA256, contig count, and total bp.
    * rejects merely syntactically valid or partial FASTA files without verified manifest;
      forces re-extraction from archive stream if unverified.
    * fails closed immediately on missing required inputs (Zenodo JSON, sampled IDs, etc.).
- Gene-call retention:
    * preserves per-genome raw protein translations: genecalls/smag/{mag_id}.raw.faa
    * preserves per-genome gene coordinate outputs:  genecalls/smag/{mag_id}.gff
    * records original-ID-to-canonical-ID mapping:   genecalls/smag/{mag_id}.id_map.tsv
      linking matched canonical protein IDs to original contigs, start/end coordinates, and strand.
    * caller flags remain biologically identical: -p meta --no-stop-codon.
    * preserves both original raw outputs and canonical deduplicated pools.
    * records paths and SHA256 hashes in results/smag_genome_manifest.tsv and writes
      an aggregated coordinate mapping table in results/smag_protein_coordinate_mapping.tsv.

Usage:
    python3 scripts/02_extract_and_call_smag.py [threads]
"""

import os
import sys
import io
import tarfile
import hashlib
import subprocess
import time
import json
import concurrent.futures

ZENODO_JSON_PATH = os.environ.get("ZENODO_JSON_PATH", "results/zenodo_8223844.json")
SAMPLED_IDS_PATH = os.environ.get("SAMPLED_IDS_PATH", "results/sampled_smag_1500_ids_seed20261014.txt")
PARTS_DIR = os.environ.get("PARTS_DIR", "downloads/smag_parts")
GENOMES_DIR = os.environ.get("GENOMES_DIR", "genomes/smag_extracted")
GENECALLS_DIR = os.environ.get("GENECALLS_DIR", "genecalls/smag")
RESULTS_DIR = os.environ.get("RESULTS_DIR", "results")
PYRODIGAL_BIN = os.environ.get("PYRODIGAL_BIN", "/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal")

EXTRACTION_MANIFEST_PATH = os.path.join(RESULTS_DIR, "smag_extraction_manifest.json")

VALID_NUC = set(b"ACGTURYKMSWBDHVNacgturykmswbdhvn")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")


class MultiFileStream(io.RawIOBase):
    """Stream multiple ordered archive part files as a single contiguous byte stream."""

    def __init__(self, paths):
        self.paths = list(paths)
        self.current_idx = 0
        self.f = open(self.paths[self.current_idx], "rb") if self.paths else None

    def readinto(self, b):
        if self.f is None:
            return 0
        n = self.f.readinto(b)
        while n == 0 and self.current_idx + 1 < len(self.paths):
            self.f.close()
            self.current_idx += 1
            self.f = open(self.paths[self.current_idx], "rb")
            n = self.f.readinto(b)
        return n

    def readable(self):
        return True

    def close(self):
        if self.f and not self.f.closed:
            self.f.close()
        super().close()


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


def validate_fasta_bytes(buf):
    starts_gt = buf.startswith(b">")
    ends_nl = buf.endswith(b"\n")
    if not starts_gt or not ends_nl:
        return False, 0, 0, starts_gt, ends_nl
    contig_count = 0
    total_bp = 0
    for line in buf.splitlines():
        if not line:
            continue
        if line.startswith(b">"):
            contig_count += 1
        else:
            total_bp += len(line)
            if not set(line).issubset(VALID_NUC):
                return False, contig_count, total_bp, starts_gt, ends_nl
    return (contig_count > 0 and total_bp > 0), contig_count, total_bp, starts_gt, ends_nl


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

    # Prodigal standard format: ID # start # end # strand # attributes
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

    # Derive contig ID from original gene ID (typically contig_geneidx)
    if "_" in orig_id:
        contig_id = orig_id.rsplit("_", 1)[0]
    else:
        contig_id = orig_id

    return orig_id, contig_id, start_coord, end_coord, strand, partial_val


def call_genes_single_genome(args):
    """
    Run Pyrodigal v3.7.1 on a single genome assembly.
    Retains per-genome raw FAA, gene coordinate GFF, and ID-to-coordinate mapping.
    Caller flags: -p meta --no-stop-codon -f gff.
    """
    gpath, mag_id, pyrodigal_bin, genecalls_dir = args

    raw_faa_path = os.path.join(genecalls_dir, f"{mag_id}.raw.faa")
    gff_path = os.path.join(genecalls_dir, f"{mag_id}.gff")
    id_map_path = os.path.join(genecalls_dir, f"{mag_id}.id_map.tsv")

    tmp_raw_faa = raw_faa_path + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"
    tmp_gff = gff_path + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"
    tmp_id_map = id_map_path + f".tmp_{os.getpid()}_{int(time.time() * 1000)}"

    cmd = [
        pyrodigal_bin,
        "-i", gpath,
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
        print(f"Pyrodigal failed on {mag_id} (code {res.returncode}):\nSTDERR:\n{res.stderr}\nSTDOUT:\n{res.stdout}", file=sys.stderr)
        raise RuntimeError(f"Pyrodigal failed for {mag_id} (code {res.returncode}): {res.stderr.strip()}")

    # Read bytes and compute hashes before renaming
    with open(tmp_raw_faa, "rb") as bf:
        raw_faa_bytes = bf.read()
    raw_faa_sha256 = hashlib.sha256(raw_faa_bytes).hexdigest()

    with open(tmp_gff, "rb") as gf:
        gff_bytes = gf.read()
    gff_sha256 = hashlib.sha256(gff_bytes).hexdigest()

    # Atomically place raw FAA and GFF into final retained paths
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
        canon_seq_id = f"{mag_id}::{mag_id}_{gene_idx}"

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
            "mag_id": mag_id,
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
            "mag_id\tgene_idx\tcanonical_seq_id\toriginal_protein_id\tcontig_id\tstart_coord\tend_coord\t"
            "strand\tpartial\tlength_aa\tcanonical_seq_sha256\tstatus\traw_header\n"
        )
        for mr in mapping_records:
            out_m.write(
                f"{mr['mag_id']}\t{mr['gene_idx']}\t{mr['canonical_seq_id']}\t{mr['original_protein_id']}\t"
                f"{mr['contig_id']}\t{mr['start_coord']}\t{mr['end_coord']}\t{mr['strand']}\t{mr['partial']}\t"
                f"{mr['length_aa']}\t{mr['canonical_seq_sha256']}\t{mr['status']}\t{mr['raw_header']}\n"
            )
        out_m.flush()
        os.fsync(out_m.fileno())
    os.replace(tmp_id_map, id_map_path)
    id_map_sha256 = sha256_file(id_map_path)

    stats = {
        "mag_id": mag_id,
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
    return mag_id, stats, records, mapping_records


def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    os.makedirs(GENOMES_DIR, exist_ok=True)
    os.makedirs(GENECALLS_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print(f"=== H3 Step 2: Extracting & Gene Calling 1,500 SMAG Genomes ({threads} threads) ===")

    # Fail closed on missing required inputs
    for required_file in (SAMPLED_IDS_PATH, ZENODO_JSON_PATH):
        if not os.path.exists(required_file):
            raise FileNotFoundError(f"Missing required input file: '{required_file}'.")
    if not os.path.exists(PYRODIGAL_BIN):
        raise FileNotFoundError(f"Pyrodigal binary not found at '{PYRODIGAL_BIN}'.")

    with open(SAMPLED_IDS_PATH, "r", encoding="utf-8") as f:
        target_ids = sorted(line.strip() for line in f if line.strip())
    assert len(target_ids) == 1500, f"Expected exactly 1,500 prespecified IDs, got {len(target_ids)}"
    target_set = set(target_ids)
    print(f"Target sampled MAGs: {len(target_ids):,}")

    # Checked cache semantics for existing extracted assemblies
    # Requires verified extraction manifest; does not treat unverified or partial FASTA as proof
    existing_manifest = {}
    if os.path.exists(EXTRACTION_MANIFEST_PATH):
        try:
            with open(EXTRACTION_MANIFEST_PATH, "r", encoding="utf-8") as mf:
                existing_manifest = json.load(mf)
            print(f"Loaded existing extraction manifest with {len(existing_manifest):,} entries.")
        except Exception as e:
            print(f"Warning: Failed to load extraction manifest '{EXTRACTION_MANIFEST_PATH}': {e}. Re-verifying all assemblies.")
            existing_manifest = {}

    retained_genomes = {}
    for mid in target_ids:
        gpath = os.path.join(GENOMES_DIR, f"{mid}.fa")
        if os.path.exists(gpath) and mid in existing_manifest:
            m_info = existing_manifest[mid]
            with open(gpath, "rb") as bf:
                buf = bf.read()
            obs_bytes = len(buf)
            obs_sha = hashlib.sha256(buf).hexdigest()
            is_valid, n_contigs, total_bp, _, _ = validate_fasta_bytes(buf)

            # Strict verification against pinned extraction manifest
            if (
                is_valid
                and obs_bytes == m_info.get("declared_bytes")
                and obs_bytes == m_info.get("observed_bytes")
                and obs_sha == m_info.get("nucleotide_sha256")
                and n_contigs == m_info.get("contig_count")
                and total_bp == m_info.get("total_bp")
            ):
                retained_genomes[mid] = {
                    "mag_id": mid,
                    "tar_member": m_info.get("tar_member", f"{mid}.fa"),
                    "declared_bytes": m_info["declared_bytes"],
                    "observed_bytes": obs_bytes,
                    "nucleotide_sha256": obs_sha,
                    "contig_count": n_contigs,
                    "total_bp": total_bp,
                    "fasta_path": gpath,
                }
            else:
                print(f"  Existing FASTA for '{mid}' failed manifest checksum/size validation. Marking for re-extraction.")
        elif os.path.exists(gpath):
            print(f"  Existing FASTA for '{mid}' lacks verified extraction manifest record. Marking for re-extraction.")

    print(f"Verified cached extractions on disk: {len(retained_genomes)}/{len(target_ids)}")

    remaining_targets = target_set - set(retained_genomes.keys())
    if remaining_targets:
        print(f"Streaming archive to extract remaining {len(remaining_targets)} genomes...")
        with open(ZENODO_JSON_PATH, "r", encoding="utf-8") as f:
            meta = json.load(f)
        raw_parts = [p for p in meta.get("files", []) if p.get("key", "").startswith("mag.tar.gz.")]
        sorted_parts = sorted(raw_parts, key=lambda x: x["key"].replace("mag.tar.gz.", ""))
        part_paths = [os.path.join(PARTS_DIR, p["key"]) for p in sorted_parts]

        available_parts = [p for p in part_paths if os.path.exists(p)]
        print(f"Available archive parts for streaming: {len(available_parts)}/{len(part_paths)}")
        if len(part_paths) != 373 or len(available_parts) != len(part_paths):
            raise FileNotFoundError(f"Full analysis requires all 373 archive parts; found {len(available_parts)}/{len(part_paths)}")
        for info, part_path in zip(sorted_parts, part_paths):
            if os.path.getsize(part_path) != info["size"]:
                raise ValueError(f"Archive part size mismatch: '{part_path}' (expected {info['size']}, got {os.path.getsize(part_path)})")

        raw_stream = MultiFileStream(available_parts)
        buf_stream = io.BufferedReader(raw_stream, buffer_size=4 * 1024 * 1024)

        t0 = time.time()
        member_cnt = 0
        try:
            with tarfile.open(fileobj=buf_stream, mode="r|gz") as tf:
                while remaining_targets:
                    try:
                        m = tf.next()
                    except Exception as e:
                        print(f"  Tar stream ended: {e}")
                        break
                    if m is None:
                        break
                    if not m.isfile():
                        continue
                    member_cnt += 1
                    base_name = os.path.basename(m.name)
                    mid = os.path.splitext(base_name)[0]

                    if mid in remaining_targets:
                        fobj = tf.extractfile(m)
                        buf = fobj.read()
                        is_valid, n_contigs, total_bp, s_gt, e_nl = validate_fasta_bytes(buf)
                        if not (is_valid and len(buf) == m.size):
                            raise RuntimeError(f"Extracted member '{m.name}' failed byte size or FASTA validation!")

                        dest_fa = os.path.join(GENOMES_DIR, f"{mid}.fa")
                        tmp_dest_fa = dest_fa + f".tmp_{os.getpid()}"
                        with open(tmp_dest_fa, "wb") as out_fa:
                            out_fa.write(buf)
                            out_fa.flush()
                            os.fsync(out_fa.fileno())
                        os.replace(tmp_dest_fa, dest_fa)

                        entry = {
                            "mag_id": mid,
                            "tar_member": m.name,
                            "declared_bytes": m.size,
                            "observed_bytes": len(buf),
                            "nucleotide_sha256": hashlib.sha256(buf).hexdigest(),
                            "contig_count": n_contigs,
                            "total_bp": total_bp,
                            "fasta_path": dest_fa,
                        }
                        retained_genomes[mid] = entry
                        existing_manifest[mid] = entry
                        remaining_targets.remove(mid)
                        if len(retained_genomes) % 100 == 0 or not remaining_targets:
                            print(f"  Extracted {len(retained_genomes)}/{len(target_ids)} target genomes ({time.time() - t0:.1f}s)...", flush=True)
        finally:
            buf_stream.close()

    # Fail closed if target set is incomplete
    if set(retained_genomes.keys()) != target_set:
        missing = sorted(target_set - set(retained_genomes.keys()))
        raise RuntimeError(f"Incomplete prespecified sample: {len(retained_genomes)}/1500; missing IDs: {missing}")

    # Save updated extraction manifest atomically
    atomic_write_json(EXTRACTION_MANIFEST_PATH, existing_manifest)
    print(f"Extraction manifest verified and saved: '{EXTRACTION_MANIFEST_PATH}'.")

    # Call genes on all 1,500 genomes with per-genome retention
    sorted_mids = sorted(retained_genomes.keys())
    tasks = [(retained_genomes[m]["fasta_path"], m, PYRODIGAL_BIN, GENECALLS_DIR) for m in sorted_mids]
    print(f"=== Calling genes on {len(tasks)} SMAG genomes with {threads} workers (Retaining raw FAA & GFF) ===")
    t_call = time.time()
    results_by_mag = {}
    all_mapping_records = []

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(call_genes_single_genome, t) for t in tasks]
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            done += 1
            mag_id, st, recs, map_recs = fut.result()
            results_by_mag[mag_id] = (st, recs)
            all_mapping_records.extend(map_recs)
            if done % 100 == 0 or done == len(tasks):
                print(f"  Called {done}/{len(tasks)} MAGs ({time.time() - t_call:.1f}s)...", flush=True)

    # Write canonical FASTA, deduplicated unique FASTA, duplicate map, and genome manifest
    canon_faa_path = os.path.join(RESULTS_DIR, "smag_pool_canonical.faa")
    uniq_faa_path = os.path.join(RESULTS_DIR, "smag_pool_unique.faa")
    dup_map_path = os.path.join(RESULTS_DIR, "smag_duplicate_map.tsv")
    manifest_path = os.path.join(RESULTS_DIR, "smag_genome_manifest.tsv")
    agg_map_path = os.path.join(RESULTS_DIR, "smag_protein_coordinate_mapping.tsv")

    seen_seq_sha256 = {}
    dup_rows = []
    manifest_rows = []

    total_canon_prots = 0
    total_canon_res = 0
    total_uniq_prots = 0
    total_uniq_res = 0

    tmp_canon_path = canon_faa_path + f".tmp_{os.getpid()}"
    tmp_uniq_path = uniq_faa_path + f".tmp_{os.getpid()}"

    with open(tmp_canon_path, "w", encoding="utf-8") as f_canon, open(tmp_uniq_path, "w", encoding="utf-8") as f_uniq:
        for mid in sorted_mids:
            meta = retained_genomes[mid]
            st, recs = results_by_mag[mid]
            recs.sort(key=lambda x: x[0])
            for gene_idx, seq_id, canon_seq, seq_sha256 in recs:
                f_canon.write(f">{seq_id}\n{canon_seq}\n")
                total_canon_prots += 1
                total_canon_res += len(canon_seq)
                if seq_sha256 not in seen_seq_sha256:
                    seen_seq_sha256[seq_sha256] = seq_id
                    f_uniq.write(f">{seq_id}\n{canon_seq}\n")
                    total_uniq_prots += 1
                    total_uniq_res += len(canon_seq)
                else:
                    dup_rows.append((seq_id, seen_seq_sha256[seq_sha256], len(canon_seq), seq_sha256))

            manifest_rows.append((
                mid,
                meta["declared_bytes"],
                meta["observed_bytes"],
                meta["nucleotide_sha256"],
                meta["contig_count"],
                meta["total_bp"],
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

    os.replace(tmp_canon_path, canon_faa_path)
    os.replace(tmp_uniq_path, uniq_faa_path)

    # Write duplicate map
    tmp_dup_path = dup_map_path + f".tmp_{os.getpid()}"
    with open(tmp_dup_path, "w", encoding="utf-8") as out:
        out.write("duplicate_seq_id\trepresentative_seq_id\tlength_aa\tseq_sha256\n")
        for r in dup_rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp_dup_path, dup_map_path)

    # Write genome manifest
    tmp_manifest_path = manifest_path + f".tmp_{os.getpid()}"
    with open(tmp_manifest_path, "w", encoding="utf-8") as out:
        out.write("mag_id\tdeclared_bytes\tobserved_bytes\tnucleotide_sha256\tcontig_count\ttotal_bp\t"
                  "raw_protein_count\traw_residue_count\ttrailing_star_count\tinternal_star_count\t"
                  "lowercase_char_count\tinvalid_aa_char_count\tcanonical_protein_count\tcanonical_residue_count\t"
                  "unique_canonical_proteins\traw_faa_sha256\traw_faa_path\tgff_path\tgff_sha256\t"
                  "id_map_path\tid_map_sha256\n")
        for r in manifest_rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    os.replace(tmp_manifest_path, manifest_path)

    # Write aggregated coordinate mapping
    tmp_agg_path = agg_map_path + f".tmp_{os.getpid()}"
    with open(tmp_agg_path, "w", encoding="utf-8") as out_agg:
        out_agg.write(
            "mag_id\tgene_idx\tcanonical_seq_id\toriginal_protein_id\tcontig_id\tstart_coord\tend_coord\t"
            "strand\tpartial\tlength_aa\tcanonical_seq_sha256\tstatus\traw_header\n"
        )
        for mr in all_mapping_records:
            out_agg.write(
                f"{mr['mag_id']}\t{mr['gene_idx']}\t{mr['canonical_seq_id']}\t{mr['original_protein_id']}\t"
                f"{mr['contig_id']}\t{mr['start_coord']}\t{mr['end_coord']}\t{mr['strand']}\t{mr['partial']}\t"
                f"{mr['length_aa']}\t{mr['canonical_seq_sha256']}\t{mr['status']}\t{mr['raw_header']}\n"
            )
    os.replace(tmp_agg_path, agg_map_path)

    print(f"=== SMAG Processing Summary ===")
    print(f"  Genomes processed: {len(sorted_mids):,}")
    print(f"  Canonical proteins: {total_canon_prots:,} ({total_canon_res:,} aa)")
    print(f"  Unique canonical proteins: {total_uniq_prots:,} ({total_uniq_res:,} aa)")
    print(f"  Duplicate proteins: {len(dup_rows):,} ({100.0 * len(dup_rows) / max(total_canon_prots, 1):.2f}% redundancy)")
    print(f"  Retained genecalls in: {GENECALLS_DIR}/ ({len(sorted_mids):,} raw FAA, GFF, and ID-map files)")
    print(f"  Files written: {canon_faa_path}, {uniq_faa_path}, {dup_map_path}, {manifest_path}, {agg_map_path}")


if __name__ == "__main__":
    main()
