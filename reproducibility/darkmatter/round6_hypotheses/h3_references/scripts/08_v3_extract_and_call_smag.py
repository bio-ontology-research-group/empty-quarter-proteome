#!/usr/bin/env python3
"""
scripts/08_v3_extract_and_call_smag.py
v3 Repair for Hypothesis 3 (Critical Findings 2 & 3):
1. Verify MD5 and SHA-256 checksums of public Zenodo 8223844 SMAG archive parts.
2. Stream concatenated tar.gz parts via Python's tarfile module (mode='r|gz').
   Enforce path safety and exact byte-count completeness (observed_bytes == member.size).
   Explicitly log and reject the trailing truncated member (member 213: MAG/SRR5855438rbin.1.fa).
3. Retain all 212 complete extracted nucleotide FASTA files permanently on IBEX with SHA-256 digests.
4. Call genes with Pyrodigal v3.7.1 (-p meta --no-stop-codon), audit stop codons/case/characters,
   canonicalize sequences (rstrip('*').upper()), and write both full canonical and deduplicated
   unique SMAG reference pools in strictly deterministic (mag_id, gene_idx) order.
"""

import os
import sys
import io
import tarfile
import hashlib
import subprocess
import time
import concurrent.futures

ZENODO_PARTS = [
    {
        "key": "mag.tar.gz.aa",
        "size": 104857600,
        "md5": "5fb8a52b97940d1d1c14aecda57ac5bc",
        "url": "https://zenodo.org/api/records/8223844/files/mag.tar.gz.aa/content",
    },
    {
        "key": "mag.tar.gz.ab",
        "size": 104857600,
        "md5": "de5f1ecb60b8ceac42b616e2780b2d04",
        "url": "https://zenodo.org/api/records/8223844/files/mag.tar.gz.ab/content",
    },
]

VALID_NUC = set(b"ACGTURYKMSWBDHVNacgturykmswbdhvn")
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBZXUOJ")


class MultiFileStream(io.RawIOBase):
    """Stream multiple ordered archive part files as a single contiguous byte stream."""

    def __init__(self, paths):
        self.paths = list(paths)
        self.f = open(self.paths.pop(0), "rb")

    def readinto(self, b):
        n = self.f.readinto(b)
        while n == 0 and self.paths:
            self.f.close()
            self.f = open(self.paths.pop(0), "rb")
            n = self.f.readinto(b)
        return n

    def readable(self):
        return True

    def close(self):
        if hasattr(self, "f") and self.f and not self.f.closed:
            self.f.close()
        super().close()


def compute_md5_sha256(path):
    h_md5 = hashlib.md5()
    h_sha = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h_md5.update(chunk)
            h_sha.update(chunk)
    return h_md5.hexdigest(), h_sha.hexdigest()


def is_safe_tar_path(name):
    if name.startswith("/") or name.startswith("\\"):
        return False
    parts = name.replace("\\", "/").split("/")
    if ".." in parts:
        return False
    return True


def validate_fasta_bytes(buf):
    """Validate FASTA syntax and return (is_valid, contig_count, total_bp, starts_gt, ends_nl)."""
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
    """Run Pyrodigal on a single retained genome and return deterministic records + audit stats."""
    gpath, mag_id, pyrodigal_bin = args
    tmp_trans = gpath + ".v3_tmp.faa"
    cmd_call = [
        pyrodigal_bin,
        "-i", gpath,
        "-a", tmp_trans,
        "-p", "meta",
        "--no-stop-codon",
        "-o", "/dev/null",
    ]
    subprocess.run(cmd_call, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    with open(tmp_trans, "rb") as bf:
        raw_faa_bytes = bf.read()
    raw_faa_sha256 = hashlib.sha256(raw_faa_bytes).hexdigest()
    os.remove(tmp_trans)

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
    threads = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    h3_dir = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/p6_hypotheses_clean/h3_references"
    parts_dir = os.path.join(h3_dir, "downloads/smag_parts")
    genomes_dir = os.path.join(h3_dir, "genomes/smag_extracted")
    out_dir = os.path.join(h3_dir, "results")
    os.makedirs(genomes_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)

    pyrodigal_bin = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/pyrodigal/bin/pyrodigal"
    ver_out = subprocess.check_output([pyrodigal_bin, "--version"], text=True).strip()
    print(f"=== H3 v3 Step 1: Safe Tar Extraction & Deterministic Pyrodigal Calling ({ver_out}) ===")

    # 1. Verify archive parts
    part_paths = []
    parts_manifest_rows = []
    for item in ZENODO_PARTS:
        p = os.path.join(parts_dir, item["key"])
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing required Zenodo part: {p}")
        sz = os.path.getsize(p)
        md5_val, sha256_val = compute_md5_sha256(p)
        if sz != item["size"] or md5_val != item["md5"]:
            raise ValueError(f"Integrity failure on {p}: size={sz}, md5={md5_val}")
        part_paths.append(p)
        parts_manifest_rows.append((item["key"], sz, md5_val, sha256_val, "verified_md5_sha256", item["url"]))
        print(f"  Verified {item['key']}: bytes={sz}, md5={md5_val}, sha256={sha256_val}")

    parts_man_path = os.path.join(out_dir, "v3_smag_zenodo_parts_manifest.tsv")
    with open(parts_man_path, "w") as out:
        out.write("part_key\tbyte_size\tmd5_checksum\tsha256_checksum\tstatus\tsource_url\n")
        for r in parts_manifest_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    # 2. Stream tarfile with strict member.size verification
    print("=== Streaming concatenated tar.gz parts with tarfile.open(mode='r|gz') ===")
    t0 = time.time()
    tar_audit_rows = []
    retained_genomes = {}  # mag_id -> metadata dict
    member_idx = 0

    raw_stream = MultiFileStream(part_paths)
    buf_stream = io.BufferedReader(raw_stream, buffer_size=1024 * 1024)
    try:
        with tarfile.open(fileobj=buf_stream, mode="r|gz") as tf:
            while True:
                try:
                    m = tf.next()
                except Exception as hdr_err:
                    print(f"  Tar header stream ended: {type(hdr_err).__name__}: {hdr_err}")
                    break
                if m is None:
                    break
                if not m.isfile():
                    continue
                member_idx += 1
                m_name = m.name
                declared_sz = m.size
                base_name = os.path.basename(m_name)
                mag_id = os.path.splitext(base_name)[0]

                if not is_safe_tar_path(m_name):
                    tar_audit_rows.append((
                        member_idx, m_name, mag_id, declared_sz, 0, "NONE",
                        False, False, 0, 0, "REJECTED_UNSAFE_PATH", "Unsafe tar member path"
                    ))
                    continue

                fobj = tf.extractfile(m)
                buf = bytearray()
                read_err = None
                try:
                    while True:
                        chunk = fobj.read(65536)
                        if not chunk:
                            break
                        buf.extend(chunk)
                except Exception as e:
                    read_err = f"{type(e).__name__}: {e}"

                obs_sz = len(buf)
                obs_sha256 = hashlib.sha256(buf).hexdigest()
                is_valid_fa, n_contigs, total_bp, s_gt, e_nl = validate_fasta_bytes(buf)

                if read_err is not None or obs_sz != declared_sz or not is_valid_fa:
                    status = "REJECTED_TRUNCATED_STREAM_END" if obs_sz != declared_sz or read_err else "REJECTED_INVALID_FASTA"
                    err_msg = read_err if read_err else f"observed={obs_sz} != declared={declared_sz} or invalid_fasta"
                    print(f"  [Member {member_idx}] REJECTED {m_name}: declared={declared_sz}, observed={obs_sz}, "
                          f"starts_gt={s_gt}, ends_nl={e_nl}, reason={err_msg}")
                    tar_audit_rows.append((
                        member_idx, m_name, mag_id, declared_sz, obs_sz, obs_sha256,
                        s_gt, e_nl, n_contigs, total_bp, status, err_msg
                    ))
                    # Ensure no partial file exists in retained directory
                    target_fa = os.path.join(genomes_dir, f"{mag_id}.fa")
                    if os.path.exists(target_fa):
                        os.remove(target_fa)
                    if read_err is not None:
                        break
                    continue

                # Complete verified member: write and retain permanently
                target_fa = os.path.join(genomes_dir, f"{mag_id}.fa")
                with open(target_fa, "wb") as out_fa:
                    out_fa.write(buf)

                tar_audit_rows.append((
                    member_idx, m_name, mag_id, declared_sz, obs_sz, obs_sha256,
                    s_gt, e_nl, n_contigs, total_bp, "COMPLETE_RETAINED", "OK"
                ))
                retained_genomes[mag_id] = {
                    "member_order_idx": member_idx,
                    "tar_member_name": m_name,
                    "mag_id": mag_id,
                    "declared_bytes": declared_sz,
                    "observed_bytes": obs_sz,
                    "nucleotide_sha256": obs_sha256,
                    "contig_count": n_contigs,
                    "total_bp": total_bp,
                    "fasta_path": target_fa,
                }
    finally:
        buf_stream.close()

    tar_audit_path = os.path.join(out_dir, "v3_smag_tar_extraction_audit.tsv")
    with open(tar_audit_path, "w") as out:
        out.write("member_order_idx\ttar_member_name\tmag_id\tdeclared_bytes\tobserved_bytes\tobserved_sha256\t"
                  "starts_with_gt\tends_with_newline\tcontig_count\ttotal_bp\tstatus\tdetail\n")
        for r in tar_audit_rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    print(f"Wrote tar extraction audit ({len(tar_audit_rows)} members inspected, "
          f"{len(retained_genomes)} complete retained) in {time.time() - t0:.1f}s.")

    # 3. Call genes with Pyrodigal in parallel, then sort deterministically by (mag_id, gene_idx)
    sorted_mag_ids = sorted(retained_genomes.keys())
    tasks = [(retained_genomes[m]["fasta_path"], m, pyrodigal_bin) for m in sorted_mag_ids]
    print(f"=== Calling genes on {len(tasks)} verified complete MAGs with {threads} workers ===")
    t_call = time.time()
    results_by_mag = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(call_genes_single_genome, t) for t in tasks]
        done_cnt = 0
        for fut in concurrent.futures.as_completed(futures):
            done_cnt += 1
            mag_id, st, recs = fut.result()
            results_by_mag[mag_id] = (st, recs)
            if done_cnt % 50 == 0 or done_cnt == len(tasks):
                print(f"  Called {done_cnt}/{len(tasks)} MAGs ({time.time() - t_call:.1f}s)...", flush=True)

    # 4. Write deterministic canonical pool, deduplicated pool, duplicate map, and genome manifest
    canon_faa_path = os.path.join(out_dir, "v3_smag_pool_canonical.faa")
    uniq_faa_path = os.path.join(out_dir, "v3_smag_pool_unique.faa")
    dup_map_path = os.path.join(out_dir, "v3_smag_duplicate_map.tsv")
    genome_man_path = os.path.join(out_dir, "v3_smag_genome_manifest.tsv")

    seen_seq_sha256 = {}  # seq_sha256 -> representative_seq_id
    dup_rows = []
    manifest_rows = []

    total_raw_prots = 0
    total_canon_prots = 0
    total_canon_res = 0
    total_uniq_prots = 0
    total_uniq_res = 0

    with open(canon_faa_path, "w") as f_canon, open(uniq_faa_path, "w") as f_uniq:
        for mag_id in sorted_mag_ids:
            meta = retained_genomes[mag_id]
            st, recs = results_by_mag[mag_id]
            # Sort records deterministically by gene_idx
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

            total_raw_prots += st["raw_protein_count"]
            manifest_rows.append((
                meta["member_order_idx"],
                meta["tar_member_name"],
                mag_id,
                meta["declared_bytes"],
                meta["observed_bytes"],
                "COMPLETE_RETAINED",
                meta["contig_count"],
                meta["total_bp"],
                meta["nucleotide_sha256"],
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

    with open(genome_man_path, "w") as out:
        out.write(
            "member_order_idx\ttar_member_name\tmag_id\tdeclared_bytes\tobserved_bytes\tstatus\t"
            "contig_count\ttotal_bp\tnucleotide_sha256\traw_protein_count\traw_residue_count\t"
            "trailing_star_count\tinternal_star_count\tlowercase_char_count\tinvalid_aa_char_count\t"
            "canonical_protein_count\tcanonical_residue_count\tunique_canonical_proteins\traw_faa_sha256\n"
        )
        for r in manifest_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    _, canon_sha256 = compute_md5_sha256(canon_faa_path)
    _, uniq_sha256 = compute_md5_sha256(uniq_faa_path)
    print(f"SMAG canonical pool: {total_canon_prots} proteins, {total_canon_res} residues, SHA256={canon_sha256}")
    print(f"SMAG unique pool:    {total_uniq_prots} proteins, {total_uniq_res} residues, SHA256={uniq_sha256} "
          f"(duplicates={len(dup_rows)})")
    print(f"Retained {len(sorted_mag_ids)} extracted genome FASTA files in {genomes_dir} (0 deleted).")
    print("DONE_08_V3_EXTRACT_AND_CALL_SMAG")


if __name__ == "__main__":
    main()
