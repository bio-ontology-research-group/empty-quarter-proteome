#!/usr/bin/env python3
"""
scripts/02_extract_and_call_smag.py
Extract the pre-declared random sample of 1,500 SMAG MAGs from the streamed
multi-part Zenodo archive (Record 8223844), verify individual nucleotide FASTAs,
call genes with Pyrodigal v3.7.1 (-p meta --no-stop-codon), canonicalize,
deduplicate, and emit complete manifests.

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

ZENODO_JSON_PATH = "results/zenodo_8223844.json"
SAMPLED_IDS_PATH = "results/sampled_smag_1500_ids_seed20261014.txt"
PARTS_DIR = "downloads/smag_parts"
GENOMES_DIR = "genomes/smag_extracted"
RESULTS_DIR = "results"
PYRODIGAL_BIN = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"

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


def is_safe_tar_path(name):
    if name.startswith("/") or name.startswith("\\"):
        return False
    parts = name.replace("\\", "/").split("/")
    if ".." in parts:
        return False
    return True


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


def call_genes_single_genome(args):
    gpath, mag_id, pyrodigal_bin = args
    tmp_trans = gpath + ".call_tmp.faa"
    cmd = [
        pyrodigal_bin,
        "-i", gpath,
        "-a", tmp_trans,
        "-p", "meta",
        "--no-stop-codon",
        "-o", "/dev/null",
    ]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    with open(tmp_trans, "rb") as bf:
        raw_faa_bytes = bf.read()
    raw_faa_sha256 = hashlib.sha256(raw_faa_bytes).hexdigest()
    try:
        os.remove(tmp_trans)
    except OSError:
        pass

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
    within_genome_hashes = set()

    cur_idx = 0
    cur_lines = []

    def flush_protein(gene_idx, lines):
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
        seq_id = f"{mag_id}::{mag_id}_{gene_idx}"
        seq_sha256 = hashlib.sha256(canon_seq.encode("ascii")).hexdigest()
        within_genome_hashes.add(seq_sha256)
        records.append((gene_idx, seq_id, canon_seq, seq_sha256))

    for line in text.splitlines():
        if line.startswith(">"):
            if cur_idx > 0:
                flush_protein(cur_idx, cur_lines)
            cur_idx += 1
            cur_lines = []
        else:
            cur_lines.append(line.strip())
    if cur_idx > 0:
        flush_protein(cur_idx, cur_lines)

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
    }
    return mag_id, stats, records


def main():
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    os.makedirs(GENOMES_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)

    print(f"=== H3 Step 2: Extracting & Gene Calling 1,500 SMAG Genomes ({threads} threads) ===")
    with open(SAMPLED_IDS_PATH, "r") as f:
        target_ids = set(line.strip() for line in f if line.strip())
    print(f"Target sampled MAGs: {len(target_ids)}")

    # Check which ones are already extracted
    retained_genomes = {}
    for mid in sorted(target_ids):
        gpath = os.path.join(GENOMES_DIR, f"{mid}.fa")
        if os.path.exists(gpath):
            with open(gpath, "rb") as bf:
                buf = bf.read()
            is_valid, n_contigs, total_bp, _, _ = validate_fasta_bytes(buf)
            if is_valid:
                retained_genomes[mid] = {
                    "mag_id": mid,
                    "declared_bytes": len(buf),
                    "observed_bytes": len(buf),
                    "nucleotide_sha256": hashlib.sha256(buf).hexdigest(),
                    "contig_count": n_contigs,
                    "total_bp": total_bp,
                    "fasta_path": gpath,
                }
    print(f"Already extracted on disk: {len(retained_genomes)}/{len(target_ids)}")

    remaining_targets = target_ids - set(retained_genomes.keys())
    if remaining_targets:
        print(f"Streaming archive to extract remaining {len(remaining_targets)} genomes...")
        # Get list of all parts ordered by suffix
        with open(ZENODO_JSON_PATH, "r") as f:
            meta = json.load(f)
        raw_parts = [p for p in meta.get("files", []) if p.get("key", "").startswith("mag.tar.gz.")]
        sorted_parts = sorted(raw_parts, key=lambda x: x["key"].replace("mag.tar.gz.", ""))
        part_paths = [os.path.join(PARTS_DIR, p["key"]) for p in sorted_parts]

        # Check existing parts
        available_parts = [p for p in part_paths if os.path.exists(p)]
        print(f"Available archive parts for streaming: {len(available_parts)}/{len(part_paths)}")
        if not available_parts:
            raise FileNotFoundError(f"No archive parts found in {PARTS_DIR}!")

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
                        if is_valid and len(buf) == m.size:
                            dest_fa = os.path.join(GENOMES_DIR, f"{mid}.fa")
                            with open(dest_fa, "wb") as out_fa:
                                out_fa.write(buf)
                            retained_genomes[mid] = {
                                "mag_id": mid,
                                "declared_bytes": m.size,
                                "observed_bytes": len(buf),
                                "nucleotide_sha256": hashlib.sha256(buf).hexdigest(),
                                "contig_count": n_contigs,
                                "total_bp": total_bp,
                                "fasta_path": dest_fa,
                            }
                            remaining_targets.remove(mid)
                            if len(retained_genomes) % 100 == 0 or not remaining_targets:
                                print(f"  Extracted {len(retained_genomes)}/{len(target_ids)} target genomes ({time.time() - t0:.1f}s)...", flush=True)
        finally:
            buf_stream.close()

    if len(retained_genomes) < len(target_ids):
        print(f"WARNING: Extracted {len(retained_genomes)}/{len(target_ids)} target genomes (need remaining parts to complete all 1500).")

    # Call genes on all retained genomes
    sorted_mids = sorted(retained_genomes.keys())
    tasks = [(retained_genomes[m]["fasta_path"], m, PYRODIGAL_BIN) for m in sorted_mids]
    print(f"=== Calling genes on {len(tasks)} SMAG genomes with {threads} workers ===")
    t_call = time.time()
    results_by_mag = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(call_genes_single_genome, t) for t in tasks]
        done = 0
        for fut in concurrent.futures.as_completed(futures):
            done += 1
            mag_id, st, recs = fut.result()
            results_by_mag[mag_id] = (st, recs)
            if done % 100 == 0 or done == len(tasks):
                print(f"  Called {done}/{len(tasks)} MAGs ({time.time() - t_call:.1f}s)...", flush=True)

    # Write canonical FASTA, deduplicated unique FASTA, duplicate map, and genome manifest
    canon_faa_path = os.path.join(RESULTS_DIR, "smag_pool_canonical.faa")
    uniq_faa_path = os.path.join(RESULTS_DIR, "smag_pool_unique.faa")
    dup_map_path = os.path.join(RESULTS_DIR, "smag_duplicate_map.tsv")
    manifest_path = os.path.join(RESULTS_DIR, "smag_genome_manifest.tsv")

    seen_seq_sha256 = {}
    dup_rows = []
    manifest_rows = []

    total_canon_prots = 0
    total_canon_res = 0
    total_uniq_prots = 0
    total_uniq_res = 0

    with open(canon_faa_path, "w") as f_canon, open(uniq_faa_path, "w") as f_uniq:
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
            ))

    with open(dup_map_path, "w") as out:
        out.write("duplicate_seq_id\trepresentative_seq_id\tlength_aa\tseq_sha256\n")
        for r in dup_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    with open(manifest_path, "w") as out:
        out.write("mag_id\tdeclared_bytes\tobserved_bytes\tnucleotide_sha256\tcontig_count\ttotal_bp\t"
                  "raw_protein_count\traw_residue_count\ttrailing_star_count\tinternal_star_count\t"
                  "lowercase_char_count\tinvalid_aa_char_count\tcanonical_protein_count\tcanonical_residue_count\t"
                  "unique_canonical_proteins\traw_faa_sha256\n")
        for r in manifest_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    print(f"=== SMAG Processing Summary ===")
    print(f"  Genomes processed: {len(sorted_mids)}")
    print(f"  Canonical proteins: {total_canon_prots:,} ({total_canon_res:,} aa)")
    print(f"  Unique canonical proteins: {total_uniq_prots:,} ({total_uniq_res:,} aa)")
    print(f"  Duplicate proteins: {len(dup_rows):,} ({100.0 * len(dup_rows) / max(total_canon_prots, 1):.2f}% redundancy)")
    print(f"  Files written: {canon_faa_path}, {uniq_faa_path}, {dup_map_path}, {manifest_path}")


if __name__ == "__main__":
    main()
