#!/usr/bin/env python3
"""02_extract_neighbourhoods.py - Proven cross-caller neighbourhood extraction & locus eligibility (v3 protocol).

Resolves C1, C2, C3, C4, and C5 from round6/h2_review1:
1. Measures real nucleotide contig lengths (contig_len_bp) from verified bin FASTAs
   and evaluates caller-level flank completeness BEFORE outcome-independent primary-locus
   selection per (family, species).
2. Proves 1:1 raw-bin-to-Bakta contig identity via exact uppercase nucleotide SHA-256
   matching between output_bins/{parent}.fa.gz and bakta/{parent}/{parent}.fna PLUS
   id=contig_N, orig-id=<raw_contig> verification in bakta/{parent}/{parent}.log,
   failing closed on any hash mismatch, duplicate sequence hash without log disambiguation,
   or reverse-complement mismatch.
3. Explicitly records bins with missing Bakta/eggNOG files in bin_contig_verification.tsv
   and pilot_join_metrics.json.
4. Classifies every flanking position (m2, m1, p1, p2) into 7 mutually exclusive audit
   categories:
   - mapped_KO (exact (norm_contig, start, stop, strand) match with KEGG KO)
   - mapped_noKO (exact (norm_contig, start, stop, strand) match without KEGG KO)
   - stop3prime_matched_KO (no exact match; unique same-strand same-3'-stop match with KEGG KO)
   - stop3prime_matched_noKO (no exact match; unique same-strand same-3'-stop match without KEGG KO)
   - coordinate_mismatch_unmapped (proven bin/contig; no exact or unique 3'-stop match)
   - source_file_unavailable (unbinned, unclustered bin, or missing Bakta/eggNOG files)
   - contig_edge_truncated (flank index outside [0, L_c - 1])
5. Precomputes same-contig eligible pseudo-focal position sets I_c^prim and I_c^stop3
   applying the exact same 4-flank completeness and Bakta-mapping gates as observed loci.
6. Computes exact byte sizes and SHA-256 checksums for every fixed and dynamic input file
   read (assembly .faa.gz, bin .fa.gz, Bakta .fna/.log/.tsv, eggNOG .emapper.annotations)
   with zero DYNAMIC placeholders.
"""
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
T = Path("/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites")
META_PATH = W / "metadata/sample_metadata.tsv"
ASM_PROT_DIR = W / "04_protein_catalog/assembly_proteins"

FIXED_INPUTS = [
    ("per_rep_structures", W / "04_structure/00_structure_clean/05_validate/per_rep_structures.tsv"),
    ("fam_reps", W / "04_structure/00_structure_clean/05_validate/fam_reps.tsv"),
    ("clusters_90rep", W / "04_structure/00_structure_clean/04_cluster/clusters.tsv"),
    ("contig_to_bin_refined", W / "04_protein_catalog/contig_to_bin_refined.tsv.gz"),
    ("drep_cdb", W / "05_refine/_drep/out/data_tables/Cdb.csv"),
    ("drep_wdb", W / "05_refine/_drep/out/data_tables/Wdb.csv"),
    ("gtdbtk_r232_bac120", W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"),
    ("gtdbtk_r232_ar53", W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"),
    ("gtdbtk_refined_bac120", W / "05_refine/_gtdbtk_refined/gtdbtk.bac120.summary.tsv"),
    ("gtdbtk_refined_ar53", W / "05_refine/_gtdbtk_refined/gtdbtk.ar53.summary.tsv"),
    ("sample_metadata", META_PATH),
]

results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results")
n_workers = int(sys.argv[2]) if len(sys.argv) > 2 else 16
update_top_inputs = (sys.argv[3].lower() != "false") if len(sys.argv) > 3 else True

cand_file = results_dir / "pilot_members_candidates.tsv"
if not cand_file.exists():
    sys.exit(f"Error: {cand_file} does not exist. Run 01_select_cohort.py first.")

RC_TABLE = str.maketrans("ACGTURYKMSWBDHVNacgturykmswbdhvn", "TGCAAYRMKSWVHDBNTGCAAYRMKSWVHDBN")


def normalize_contig(name):
    return name.split("::")[-1] if name else ""


def file_sha256_and_size(path_obj):
    h = hashlib.sha256()
    size = 0
    with open(path_obj, "rb") as f:
        while True:
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            size += len(chunk)
            h.update(chunk)
    return size, h.hexdigest()


# Read sample metadata
meta = {r["sample"]: r for r in csv.DictReader(open(META_PATH), delimiter="\t")}

# Load candidate members from 01_select_cohort.py
candidates = []
samples_needed = set()
bins_needed = set()
with open(cand_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        candidates.append(row)
        samples_needed.add(row["sample"])
        if row["bin"]:
            bins_needed.add(row["bin"])

print(
    f"[flank:v3] Loaded {len(candidates)} candidate member rows across "
    f"{len(samples_needed)} samples and {len(bins_needed)} unique refined bins (workers={n_workers})",
    flush=True,
)

# Group candidate queries by sample and contig
sample_contig_queries = defaultdict(lambda: defaultdict(list))
for m in candidates:
    sample_contig_queries[m["sample"]][m["contig"]].append(m)


def process_single_sample(args):
    sample, contig_map = args
    asm_path = ASM_PROT_DIR / sample / f"{sample}.faa.gz"
    if not asm_path.exists():
        return [], 0, {}, None

    f_size, f_sha256 = file_sha256_and_size(asm_path)
    manifest_entry = {
        "role": f"assembly_faa:{sample}",
        "path": str(asm_path),
        "size_bytes": f_size,
        "sha256": f_sha256,
        "status": "verified_read",
    }

    needed_contigs = set(contig_map.keys())
    genes_by_contig = defaultdict(list)

    proc = subprocess.Popen(
        ["gzip", "-dc", str(asm_path)],
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )
    for line in proc.stdout:
        if not line.startswith(">"):
            continue
        parts = line[1:].strip().split(" # ")
        pid = parts[0].split()[0]
        if "::" not in pid:
            pid = sample + "::" + pid
        contig = pid.rsplit("_", 1)[0]
        if contig in needed_contigs:
            start = int(parts[1])
            stop = int(parts[2])
            strand = int(parts[3])
            partial = "00"
            if len(parts) > 4:
                for token in parts[4].split(";"):
                    if token.startswith("partial="):
                        partial = token.split("=")[1]
            genes_by_contig[contig].append({
                "pid": pid,
                "start": start,
                "stop": stop,
                "strand": strand,
                "partial": partial,
            })
    proc.wait()

    rows = []
    matched = 0
    contig_cds_records = {}

    for contig, q_list in contig_map.items():
        c_genes = sorted(genes_by_contig[contig], key=lambda g: (g["start"], g["stop"], g["pid"]))
        contig_cds_records[contig] = c_genes
        pid_to_idx = {g["pid"]: i for i, g in enumerate(c_genes)}
        max_cds_stop = max((g["stop"] for g in c_genes), default=0)

        for q in q_list:
            focal_pid = q["member"]
            if focal_pid not in pid_to_idx:
                continue
            matched += 1
            idx = pid_to_idx[focal_pid]
            focal_gene = c_genes[idx]

            row = dict(q)
            row["focal_start"] = focal_gene["start"]
            row["focal_stop"] = focal_gene["stop"]
            row["focal_strand"] = focal_gene["strand"]
            row["focal_partial"] = focal_gene["partial"]
            row["total_genes_on_contig"] = len(c_genes)
            row["focal_gene_index"] = idx
            row["max_cds_stop_bp"] = max_cds_stop

            for offset in [-2, -1, 1, 2]:
                flank_key = f"m{abs(offset)}" if offset < 0 else f"p{offset}"
                target_idx = idx + offset
                if 0 <= target_idx < len(c_genes):
                    tg = c_genes[target_idx]
                    dist = (
                        focal_gene["start"] - tg["stop"] - 1
                        if offset < 0
                        else tg["start"] - focal_gene["stop"] - 1
                    )
                    same_strand = 1 if tg["strand"] == focal_gene["strand"] else 0
                    row[f"{flank_key}_pid"] = tg["pid"]
                    row[f"{flank_key}_start"] = tg["start"]
                    row[f"{flank_key}_stop"] = tg["stop"]
                    row[f"{flank_key}_strand"] = tg["strand"]
                    row[f"{flank_key}_dist"] = dist
                    row[f"{flank_key}_same_strand"] = same_strand
                    row[f"{flank_key}_truncated"] = 0
                else:
                    row[f"{flank_key}_pid"] = "CONTIG_EDGE"
                    row[f"{flank_key}_start"] = -1
                    row[f"{flank_key}_stop"] = -1
                    row[f"{flank_key}_strand"] = 0
                    row[f"{flank_key}_dist"] = -1
                    row[f"{flank_key}_same_strand"] = -1
                    row[f"{flank_key}_truncated"] = 1
            rows.append(row)

    return rows, matched, contig_cds_records, manifest_entry


t0 = time.time()
sample_items = sorted(sample_contig_queries.items(), key=lambda x: x[0])
all_raw_rows = []
total_matched = 0
all_contig_cds = {}
dynamic_manifest = []

with ProcessPoolExecutor(max_workers=n_workers) as executor:
    for rows, matched, c_records, m_entry in executor.map(process_single_sample, sample_items):
        all_raw_rows.extend(rows)
        total_matched += matched
        all_contig_cds.update(c_records)
        if m_entry:
            dynamic_manifest.append(m_entry)

t1 = time.time()
print(
    f"[flank:v3] Extracted {len(all_raw_rows)} raw neighbourhood rows "
    f"({total_matched}/{len(candidates)} matched) in {t1 - t0:.1f}s",
    flush=True,
)
assert total_matched == len(candidates), f"Focal coordinate match failure: {total_matched} != {len(candidates)}"


def parse_fasta_contigs(path_obj, is_gzip=False, compute_rc=False):
    contigs = []
    op = (lambda p: gzip.open(p, "rt")) if is_gzip else (lambda p: open(p, "r"))
    with op(path_obj) as f:
        cur_hdr = None
        seq_parts = []
        for line in f:
            if line.startswith(">"):
                if cur_hdr is not None:
                    seq = "".join(seq_parts).upper()
                    sha = hashlib.sha256(seq.encode("ascii")).hexdigest()
                    rc_sha = (
                        hashlib.sha256(seq.translate(RC_TABLE)[::-1].encode("ascii")).hexdigest()
                        if compute_rc
                        else ""
                    )
                    contigs.append((cur_hdr, len(seq), sha, rc_sha))
                cur_hdr = line[1:].strip().split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if cur_hdr is not None:
            seq = "".join(seq_parts).upper()
            sha = hashlib.sha256(seq.encode("ascii")).hexdigest()
            rc_sha = (
                hashlib.sha256(seq.translate(RC_TABLE)[::-1].encode("ascii")).hexdigest()
                if compute_rc
                else ""
            )
            contigs.append((cur_hdr, len(seq), sha, rc_sha))
    return contigs


def process_single_bin(bin_name):
    s = bin_name.split("__")[0]
    m = meta.get(s)
    if not m:
        return bin_name, {
            "verification_status": "sample_metadata_missing",
            "contig_lengths": {},
            "exact_ann": {},
            "stop3_ann": {},
            "manifest": [],
            "audit_row": {
                "bin": bin_name,
                "sample": s,
                "parent_bin": "",
                "verification_status": "sample_metadata_missing",
                "n_raw_contigs": 0,
                "n_bakta_contigs": 0,
                "n_log_mapped": 0,
                "n_sha256_matched": 0,
                "n_duplicate_lengths_in_bin": 0,
                "n_duplicate_hashes_in_bin": 0,
                "n_rc_mismatches": 0,
                "n_hash_mismatches": 0,
                "bin_fa_path": "",
                "bakta_fna_path": "",
                "bakta_log_path": "",
                "bakta_tsv_path": "",
                "eggnog_path": "",
            },
        }

    rest = bin_name.split("__", 1)[1]
    if rest.startswith("single_"):
        modes = ["single_sample"]
        local = rest[len("single_") :]
    elif rest.startswith("self_"):
        modes = ["self_supervised"]
        local = rest[len("self_") :]
    else:
        modes = ["self_supervised", "single_sample"]
        local = rest

    parent = local[:-4] if local.endswith("_sub") else local
    trip = m["trip"]
    site = m["site"]
    base = T / trip / f"site{site}" / s / "wgs"

    bin_fa, bakta_fna, bakta_log, bakta_tsv, eggnog_ann = None, None, None, None, None
    for mode in modes:
        b_fa = base / "binning" / mode / "output_bins" / f"{parent}.fa.gz"
        b_fna = base / "annotation" / mode / "bakta" / parent / f"{parent}.fna"
        b_log = base / "annotation" / mode / "bakta" / parent / f"{parent}.log"
        b_tsv = base / "annotation" / mode / "bakta" / parent / f"{parent}.tsv"
        e_ann = base / "annotation" / mode / "eggnog" / parent / f"{parent}.emapper.annotations"
        if b_fa.exists():
            bin_fa = b_fa
            if b_fna.exists() and b_tsv.exists():
                bakta_fna = b_fna
                bakta_log = b_log if b_log.exists() else None
                bakta_tsv = b_tsv
                eggnog_ann = e_ann if e_ann.exists() else None
            break

    manifest_entries = []
    contig_lengths = {}
    raw_contigs = []
    if bin_fa and bin_fa.exists():
        sz, sha = file_sha256_and_size(bin_fa)
        manifest_entries.append({
            "role": f"bin_fasta:{bin_name}",
            "path": str(bin_fa),
            "size_bytes": sz,
            "sha256": sha,
            "status": "verified_read",
        })
        raw_contigs = parse_fasta_contigs(bin_fa, is_gzip=True, compute_rc=True)
        for r_hdr, r_len, _, _ in raw_contigs:
            contig_lengths[normalize_contig(r_hdr)] = r_len

    len_counts = Counter(r_len for _, r_len, _, _ in raw_contigs)
    sha_counts = Counter(r_sha for _, _, r_sha, _ in raw_contigs)
    n_dup_lengths = sum(1 for c in len_counts.values() if c > 1)
    n_dup_hashes = sum(1 for c in sha_counts.values() if c > 1)

    if not bin_fa or not bakta_fna or not bakta_tsv:
        status = "annotation_files_missing" if bin_fa else "bin_fasta_missing"
        return bin_name, {
            "verification_status": status,
            "contig_lengths": contig_lengths,
            "exact_ann": {},
            "stop3_ann": {},
            "manifest": manifest_entries,
            "audit_row": {
                "bin": bin_name,
                "sample": s,
                "parent_bin": parent,
                "verification_status": status,
                "n_raw_contigs": len(raw_contigs),
                "n_bakta_contigs": 0,
                "n_log_mapped": 0,
                "n_sha256_matched": 0,
                "n_duplicate_lengths_in_bin": n_dup_lengths,
                "n_duplicate_hashes_in_bin": n_dup_hashes,
                "n_rc_mismatches": 0,
                "n_hash_mismatches": 0,
                "bin_fa_path": str(bin_fa) if bin_fa else "",
                "bakta_fna_path": "",
                "bakta_log_path": "",
                "bakta_tsv_path": "",
                "eggnog_path": "",
            },
        }

    for role_prefix, p_obj in [
        ("bakta_fna", bakta_fna),
        ("bakta_log", bakta_log),
        ("bakta_tsv", bakta_tsv),
        ("eggnog_annotations", eggnog_ann),
    ]:
        if p_obj and p_obj.exists():
            sz, sha = file_sha256_and_size(p_obj)
            manifest_entries.append({
                "role": f"{role_prefix}:{bin_name}",
                "path": str(p_obj),
                "size_bytes": sz,
                "sha256": sha,
                "status": "verified_read",
            })

    bakta_contigs = parse_fasta_contigs(bakta_fna, is_gzip=False, compute_rc=False)

    # Parse Bakta log for explicit id=contig_N, orig-id=<raw_contig> mapping
    log_map = {}
    if bakta_log and bakta_log.exists():
        with open(bakta_log, errors="replace") as f:
            for line in f:
                if "qc: revised sequence: id=" in line and "orig-id=" in line:
                    after_id = line.split("qc: revised sequence: id=", 1)[1]
                    b_id = after_id.split(",", 1)[0].strip()
                    orig_part = after_id.split("orig-id=", 1)[1]
                    orig_id = orig_part.split(",", 1)[0].split(";", 1)[0].strip().split()[0].rstrip(",")
                    log_map[b_id] = normalize_contig(orig_id)

    raw_by_norm_id = {normalize_contig(r_hdr): (r_hdr, r_len, r_sha, r_rc_sha) for r_hdr, r_len, r_sha, r_rc_sha in raw_contigs}
    raw_by_sha = defaultdict(list)
    for r_hdr, r_len, r_sha, _ in raw_contigs:
        raw_by_sha[r_sha].append((normalize_contig(r_hdr), r_len))
    bakta_sha_counts = Counter(b_sha for _, _, b_sha, _ in bakta_contigs)

    bakta_to_raw = {}
    n_log_mapped = 0
    n_sha_matched = 0
    n_rc_mismatches = 0
    n_hash_mismatches = 0

    for b_hdr, b_len, b_sha, _ in bakta_contigs:
        if b_hdr in log_map:
            target_norm = log_map[b_hdr]
            n_log_mapped += 1
            if target_norm in raw_by_norm_id:
                _, r_len, r_sha, r_rc_sha = raw_by_norm_id[target_norm]
                if r_len == b_len and r_sha == b_sha:
                    bakta_to_raw[b_hdr] = target_norm
                    n_sha_matched += 1
                elif r_len == b_len and r_rc_sha == b_sha:
                    n_rc_mismatches += 1
                else:
                    n_hash_mismatches += 1
            else:
                n_hash_mismatches += 1
        else:
            # Fallback only if SHA-256 is unique in both raw and Bakta FASTAs
            if len(raw_by_sha.get(b_sha, [])) == 1 and bakta_sha_counts[b_sha] == 1:
                target_norm, r_len = raw_by_sha[b_sha][0]
                if r_len == b_len:
                    bakta_to_raw[b_hdr] = target_norm
                    n_sha_matched += 1
                else:
                    n_hash_mismatches += 1
            else:
                n_hash_mismatches += 1

    # Fail closed if any contig in the bin failed verification
    if (
        len(raw_contigs) != len(bakta_contigs)
        or len(bakta_to_raw) != len(raw_contigs)
        or n_rc_mismatches > 0
        or n_hash_mismatches > 0
    ):
        return bin_name, {
            "verification_status": "failed_closed_contig_identity",
            "contig_lengths": contig_lengths,
            "exact_ann": {},
            "stop3_ann": {},
            "manifest": manifest_entries,
            "audit_row": {
                "bin": bin_name,
                "sample": s,
                "parent_bin": parent,
                "verification_status": "failed_closed_contig_identity",
                "n_raw_contigs": len(raw_contigs),
                "n_bakta_contigs": len(bakta_contigs),
                "n_log_mapped": n_log_mapped,
                "n_sha256_matched": n_sha_matched,
                "n_duplicate_lengths_in_bin": n_dup_lengths,
                "n_duplicate_hashes_in_bin": n_dup_hashes,
                "n_rc_mismatches": n_rc_mismatches,
                "n_hash_mismatches": n_hash_mismatches,
                "bin_fa_path": str(bin_fa),
                "bakta_fna_path": str(bakta_fna),
                "bakta_log_path": str(bakta_log) if bakta_log else "",
                "bakta_tsv_path": str(bakta_tsv),
                "eggnog_path": str(eggnog_ann) if eggnog_ann else "",
            },
        }

    # Read eggNOG annotations
    eggnog_data = {}
    if eggnog_ann and eggnog_ann.exists():
        with open(eggnog_ann) as f:
            hdr = None
            for line in f:
                if line.startswith("#query"):
                    hdr = line.lstrip("#").rstrip("\n").split("\t")
                    continue
                if line.startswith("#") or not hdr:
                    continue
                row = dict(zip(hdr, line.rstrip("\n").split("\t")))
                locus = row.get("query")
                ko = row.get("KEGG_ko", "-")
                ko_clean = (
                    ",".join(sorted({k.replace("ko:", "").strip() for k in ko.split(",") if k.strip() and k.strip() != "-"}))
                    if ko and ko != "-"
                    else ""
                )
                eggnog_data[locus] = {
                    "ko": ko_clean,
                    "cog": row.get("COG_category", ""),
                    "desc": row.get("Description", ""),
                    "gene": row.get("Preferred_name", ""),
                }

    exact_ann = {}
    stop3_candidates = defaultdict(list)
    with open(bakta_tsv) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) >= 8 and c[1].lower() == "cds":
                b_contig = c[0]
                norm_r_contig = bakta_to_raw.get(b_contig)
                if not norm_r_contig:
                    continue
                start = int(c[2])
                stop = int(c[3])
                strand = 1 if c[4] == "+" else -1
                stop_3prime = stop if strand == 1 else start
                locus = c[5]
                prod = c[7] if len(c) > 7 else ""
                egg = eggnog_data.get(locus, {})
                rec = {
                    "locus": locus,
                    "start": start,
                    "stop": stop,
                    "strand": strand,
                    "product": prod,
                    "ko": egg.get("ko", ""),
                    "cog": egg.get("cog", ""),
                    "desc": egg.get("desc", ""),
                    "gene": egg.get("gene", ""),
                }
                exact_ann[(norm_r_contig, start, stop, strand)] = rec
                stop3_candidates[(norm_r_contig, stop_3prime, strand)].append(rec)

    # Keep only unique 3'-stop+strand matches
    stop3_ann = {k: v[0] for k, v in stop3_candidates.items() if len(v) == 1}

    ver_status = (
        "verified_sha256_and_log"
        if n_log_mapped == len(raw_contigs)
        else "verified_sha256_unique"
    )
    return bin_name, {
        "verification_status": ver_status,
        "contig_lengths": contig_lengths,
        "exact_ann": exact_ann,
        "stop3_ann": stop3_ann,
        "manifest": manifest_entries,
        "audit_row": {
            "bin": bin_name,
            "sample": s,
            "parent_bin": parent,
            "verification_status": ver_status,
            "n_raw_contigs": len(raw_contigs),
            "n_bakta_contigs": len(bakta_contigs),
            "n_log_mapped": n_log_mapped,
            "n_sha256_matched": n_sha_matched,
            "n_duplicate_lengths_in_bin": n_dup_lengths,
            "n_duplicate_hashes_in_bin": n_dup_hashes,
            "n_rc_mismatches": n_rc_mismatches,
            "n_hash_mismatches": n_hash_mismatches,
            "bin_fa_path": str(bin_fa),
            "bakta_fna_path": str(bakta_fna),
            "bakta_log_path": str(bakta_log) if bakta_log else "",
            "bakta_tsv_path": str(bakta_tsv),
            "eggnog_path": str(eggnog_ann) if eggnog_ann else "",
        },
    }


print(f"[flank:v3] Verifying contig SHA-256 + Bakta logs and reading annotations for {len(bins_needed)} bins...", flush=True)
t2 = time.time()
bin_results = {}
bin_audit_rows = []
with ProcessPoolExecutor(max_workers=n_workers) as executor:
    for b_name, res in executor.map(process_single_bin, sorted(bins_needed)):
        bin_results[b_name] = res
        bin_audit_rows.append(res["audit_row"])
        dynamic_manifest.extend(res["manifest"])

t3 = time.time()
print(
    f"[flank:v3] Completed bin verification in {t3 - t2:.1f}s: "
    f"{Counter(r['verification_status'] for r in bin_audit_rows)}",
    flush=True,
)

# Attach species info to bin_audit_rows and write bin_contig_verification.tsv
bin_to_sp_map = {r["bin"]: r["species"] for r in all_raw_rows if r["bin"]}
for r in bin_audit_rows:
    r["species"] = bin_to_sp_map.get(r["bin"], "")

bin_audit_rows.sort(key=lambda x: x["bin"])
with open(results_dir / "bin_contig_verification.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=list(bin_audit_rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(bin_audit_rows)

# Build per-contig gene audit and precompute eligible pseudo-focal position sets I_c^prim and I_c^stop3
contig_to_bin = {}
contig_to_species = {}
for r in all_raw_rows:
    contig_to_bin[r["contig"]] = r["bin"]
    contig_to_species[r["contig"]] = r["species"]

contig_eligible_prim_indices = {}
contig_eligible_stop3_indices = {}
contig_gene_rows = []

for contig, c_genes in sorted(all_contig_cds.items()):
    c_bin = contig_to_bin.get(contig, "")
    c_sp = contig_to_species.get(contig, "")
    b_res = bin_results.get(c_bin) if c_bin else None
    has_verified_ann = bool(
        b_res
        and b_res["verification_status"] in ("verified_sha256_and_log", "verified_sha256_unique")
        and c_sp != ""
    )
    exact_ann = b_res["exact_ann"] if has_verified_ann else {}
    stop3_ann = b_res["stop3_ann"] if has_verified_ann else {}
    norm_c = normalize_contig(contig)

    per_gene_records = []
    for idx, g in enumerate(c_genes):
        start, stop, strand = g["start"], g["stop"], g["strand"]
        stop3 = stop if strand == 1 else start
        exact_key = (norm_c, start, stop, strand)
        stop3_key = (norm_c, stop3, strand)

        if not has_verified_ann:
            audit = "source_file_unavailable"
            ko, cog, prod, desc = "", "", "", ""
        elif exact_key in exact_ann:
            ann = exact_ann[exact_key]
            ko = ann["ko"]
            cog = ann["cog"]
            prod = ann["product"]
            desc = ann["desc"]
            audit = "mapped_KO" if ko else "mapped_noKO"
        elif (
            stop3_key in stop3_ann
            and min(stop, stop3_ann[stop3_key]["stop"]) - max(start, stop3_ann[stop3_key]["start"]) + 1 > 0
            and (start % 3 == stop3_ann[stop3_key]["start"] % 3)
        ):
            ann = stop3_ann[stop3_key]
            ko = ann["ko"]
            cog = ann["cog"]
            prod = ann["product"]
            desc = ann["desc"]
            audit = "stop3prime_matched_KO" if ko else "stop3prime_matched_noKO"
        else:
            audit = "coordinate_mismatch_unmapped"
            ko, cog, prod, desc = "", "", "", ""

        per_gene_records.append({
            "contig": contig,
            "bin": c_bin,
            "species": c_sp,
            "gene_index": idx,
            "pid": g["pid"],
            "start": start,
            "stop": stop,
            "strand": strand,
            "ko": ko,
            "cog": cog,
            "product": prod,
            "desc": desc,
            "audit_category": audit,
        })

    n_g = len(per_gene_records)
    prim_indices = []
    stop3_indices = []
    for j in range(n_g):
        if 2 <= j <= n_g - 3 and has_verified_ann:
            flank_cats = [per_gene_records[j + w]["audit_category"] for w in (-2, -1, 1, 2)]
            if all(cat in ("mapped_KO", "mapped_noKO") for cat in flank_cats):
                prim_indices.append(j)
            if all(
                cat in ("mapped_KO", "mapped_noKO", "stop3prime_matched_KO", "stop3prime_matched_noKO")
                for cat in flank_cats
            ):
                stop3_indices.append(j)

    prim_set = set(prim_indices)
    stop3_set = set(stop3_indices)
    contig_eligible_prim_indices[contig] = prim_indices
    contig_eligible_stop3_indices[contig] = stop3_indices

    for rec in per_gene_records:
        j = rec["gene_index"]
        rec["is_eligible_focal_primary"] = 1 if j in prim_set else 0
        rec["is_eligible_focal_stop3"] = 1 if j in stop3_set else 0
        contig_gene_rows.append(rec)

# Annotate all candidate member rows and evaluate outcome-independent locus eligibility
for r in all_raw_rows:
    b = r["bin"]
    sp = r["species"]
    c = r["contig"]
    norm_c = normalize_contig(c)
    b_res = bin_results.get(b) if b else None

    if b_res and norm_c in b_res["contig_lengths"]:
        r["contig_len_bp"] = b_res["contig_lengths"][norm_c]
        r["contig_len_source"] = "verified_bin_fasta"
    else:
        r["contig_len_bp"] = r["max_cds_stop_bp"]
        r["contig_len_source"] = "unbinned_or_missing_max_cds_stop"
    del r["max_cds_stop_bp"]

    has_verified_ann = bool(
        b_res
        and b_res["verification_status"] in ("verified_sha256_and_log", "verified_sha256_unique")
        and sp != ""
    )
    exact_ann = b_res["exact_ann"] if has_verified_ann else {}
    stop3_ann = b_res["stop3_ann"] if has_verified_ann else {}

    # Focal gene annotation
    f_exact_key = (norm_c, r["focal_start"], r["focal_stop"], r["focal_strand"])
    f_stop3 = r["focal_stop"] if r["focal_strand"] == 1 else r["focal_start"]
    f_stop3_key = (norm_c, f_stop3, r["focal_strand"])
    if has_verified_ann and f_exact_key in exact_ann:
        f_ann = exact_ann[f_exact_key]
        r["focal_audit_category"] = "mapped_KO" if f_ann["ko"] else "mapped_noKO"
        r["focal_ko"] = f_ann["ko"]
        r["focal_product"] = f_ann["product"]
    elif has_verified_ann and f_stop3_key in stop3_ann:
        f_ann = stop3_ann[f_stop3_key]
        r["focal_audit_category"] = "stop3prime_matched_KO" if f_ann["ko"] else "stop3prime_matched_noKO"
        r["focal_ko"] = f_ann["ko"]
        r["focal_product"] = f_ann["product"]
    elif not has_verified_ann:
        r["focal_audit_category"] = "source_file_unavailable"
        r["focal_ko"] = ""
        r["focal_product"] = ""
    else:
        r["focal_audit_category"] = "coordinate_mismatch_unmapped"
        r["focal_ko"] = ""
        r["focal_product"] = ""

    n_present = 0
    n_exact_mapped = 0
    n_stop3_or_exact_mapped = 0

    for flank_key in ["m2", "m1", "p1", "p2"]:
        if r[f"{flank_key}_truncated"] == 1:
            audit_cat = "contig_edge_truncated"
            ko, cog, desc, prod = "", "", "", "CONTIG_EDGE"
        elif not has_verified_ann:
            n_present += 1
            audit_cat = "source_file_unavailable"
            ko, cog, desc, prod = "", "", "", "UNAVAILABLE"
        else:
            n_present += 1
            f_start = r[f"{flank_key}_start"]
            f_stop = r[f"{flank_key}_stop"]
            f_strand = r[f"{flank_key}_strand"]
            f_s3 = f_stop if f_strand == 1 else f_start
            ex_key = (norm_c, f_start, f_stop, f_strand)
            s3_key = (norm_c, f_s3, f_strand)
            if ex_key in exact_ann:
                t_ann = exact_ann[ex_key]
                ko = t_ann["ko"]
                cog = t_ann["cog"]
                desc = t_ann["desc"]
                prod = t_ann["product"]
                audit_cat = "mapped_KO" if ko else "mapped_noKO"
                n_exact_mapped += 1
                n_stop3_or_exact_mapped += 1
            elif (
                s3_key in stop3_ann
                and min(f_stop, stop3_ann[s3_key]["stop"]) - max(f_start, stop3_ann[s3_key]["start"]) + 1 > 0
                and (f_start % 3 == stop3_ann[s3_key]["start"] % 3)
            ):
                t_ann = stop3_ann[s3_key]
                ko = t_ann["ko"]
                cog = t_ann["cog"]
                desc = t_ann["desc"]
                prod = t_ann["product"]
                audit_cat = "stop3prime_matched_KO" if ko else "stop3prime_matched_noKO"
                n_stop3_or_exact_mapped += 1
            else:
                audit_cat = "coordinate_mismatch_unmapped"
                ko, cog, desc, prod = "", "", "", "COORD_MISMATCH"

        r[f"{flank_key}_audit_category"] = audit_cat
        r[f"{flank_key}_ko"] = ko
        r[f"{flank_key}_cog"] = cog
        r[f"{flank_key}_desc"] = desc
        r[f"{flank_key}_product"] = prod

        del r[f"{flank_key}_start"]
        del r[f"{flank_key}_stop"]
        del r[f"{flank_key}_strand"]

    r["n_present_flanks"] = n_present
    r["n_exact_mapped_flanks"] = n_exact_mapped
    r["n_stop3_or_exact_mapped_flanks"] = n_stop3_or_exact_mapped
    r["eligible_locus_primary"] = 1 if (has_verified_ann and n_exact_mapped == 4) else 0
    r["eligible_locus_stop3_sens"] = 1 if (has_verified_ann and n_stop3_or_exact_mapped == 4) else 0
    r["eligible_locus_partial_sens"] = 1 if (has_verified_ann and n_exact_mapped >= 1) else 0
    r["n_eligible_positions_primary"] = len(contig_eligible_prim_indices.get(c, []))
    r["n_eligible_positions_stop3"] = len(contig_eligible_stop3_indices.get(c, []))

# Outcome-independent deterministic primary-locus selection per (family, species)
by_fam_sp = defaultdict(list)
for r in all_raw_rows:
    if r["species"]:
        by_fam_sp[(r["family"], r["species"])].append(r)
    else:
        r["is_primary_locus"] = 0
        r["locus_copy_count"] = 0

for (fam, sp), group in by_fam_sp.items():
    group.sort(
        key=lambda x: (
            -x["eligible_locus_primary"],
            -x["eligible_locus_stop3_sens"],
            -x["n_exact_mapped_flanks"],
            -x["n_stop3_or_exact_mapped_flanks"],
            -x["n_present_flanks"],
            -x["contig_len_bp"],
            -x["total_genes_on_contig"],
            x["member"],
        )
    )
    for idx, r in enumerate(group):
        r["is_primary_locus"] = 1 if idx == 0 else 0
        r["locus_copy_count"] = len(group)

all_raw_rows.sort(key=lambda x: (x["arm"], x["family"], x["species"], -x["is_primary_locus"], x["member"]))

# Emit pilot_members.tsv (concise locus-level membership & eligibility table)
member_cols = [
    "family",
    "chosen_rep",
    "arm",
    "member",
    "sample",
    "contig",
    "contig_len_bp",
    "contig_len_source",
    "total_genes_on_contig",
    "bin",
    "bin_status",
    "species",
    "species_winner_genome",
    "r232_taxonomy",
    "r232_phylum",
    "is_primary_locus",
    "locus_copy_count",
    "n_present_flanks",
    "n_exact_mapped_flanks",
    "n_stop3_or_exact_mapped_flanks",
    "eligible_locus_primary",
    "eligible_locus_stop3_sens",
    "eligible_locus_partial_sens",
    "n_eligible_positions_primary",
    "n_eligible_positions_stop3",
]
with open(results_dir / "pilot_members.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=member_cols, delimiter="\t")
    w.writeheader()
    for r in all_raw_rows:
        w.writerow({k: r[k] for k in member_cols})

# Emit pilot_neighbourhoods.tsv
out_tsv = results_dir / "pilot_neighbourhoods.tsv"
with open(out_tsv, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(all_raw_rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(all_raw_rows)

# Emit pilot_contig_genes.tsv
out_genes_tsv = results_dir / "pilot_contig_genes.tsv"
with open(out_genes_tsv, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(contig_gene_rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(contig_gene_rows)

# Compute comprehensive join & eligibility metrics
flank_cats_all = Counter()
flank_cats_primary_loci = Counter()
flank_cats_non_primary = Counter()
flank_cats_eligible_primary = Counter()
flank_cats_ineligible_primary = Counter()
flank_cats_eligible_stop3 = Counter()
flank_cats_ineligible_stop3 = Counter()
n_flanks_product_all = 0

for r in all_raw_rows:
    for fk in ("m2", "m1", "p1", "p2"):
        cat = r[f"{fk}_audit_category"]
        flank_cats_all[cat] += 1
        if r["is_primary_locus"] == 1:
            flank_cats_primary_loci[cat] += 1
            if r["eligible_locus_primary"] == 1:
                flank_cats_eligible_primary[cat] += 1
            else:
                flank_cats_ineligible_primary[cat] += 1
            if r["eligible_locus_stop3_sens"] == 1:
                flank_cats_eligible_stop3[cat] += 1
            else:
                flank_cats_ineligible_stop3[cat] += 1
        else:
            flank_cats_non_primary[cat] += 1
        if r[f"{fk}_product"] not in ("", "CONTIG_EDGE", "UNAVAILABLE", "COORD_MISMATCH"):
            n_flanks_product_all += 1

bin_ver_counts = Counter(r["verification_status"] for r in bin_audit_rows)
missing_ann_bins = sorted(r["bin"] for r in bin_audit_rows if r["verification_status"] != "verified_sha256_and_log" and r["verification_status"] != "verified_sha256_unique")

n_prim_total = sum(1 for r in all_raw_rows if r["is_primary_locus"] == 1)
n_prim_t1 = sum(1 for r in all_raw_rows if r["is_primary_locus"] == 1 and r["eligible_locus_primary"] == 1)
n_prim_t2 = sum(1 for r in all_raw_rows if r["is_primary_locus"] == 1 and r["eligible_locus_stop3_sens"] == 1)
n_prim_t3 = sum(1 for r in all_raw_rows if r["is_primary_locus"] == 1 and r["eligible_locus_partial_sens"] == 1)

metrics = {
    "total_queries": len(all_raw_rows),
    "exact_focal_coordinate_matches": total_matched,
    "samples_processed": len(samples_needed),
    "unique_refined_bins": len(bins_needed),
    "bin_verification_counts": dict(bin_ver_counts),
    "bins_missing_annotations": missing_ann_bins,
    "n_bins_with_duplicate_contig_lengths": sum(1 for r in bin_audit_rows if r["n_duplicate_lengths_in_bin"] > 0),
    "n_bins_with_duplicate_contig_hashes": sum(1 for r in bin_audit_rows if r["n_duplicate_hashes_in_bin"] > 0),
    "n_contig_rc_mismatches": sum(r["n_rc_mismatches"] for r in bin_audit_rows),
    "n_contig_hash_mismatches": sum(r["n_hash_mismatches"] for r in bin_audit_rows),
    "member_bin_status_counts": dict(Counter(r["bin_status"] for r in all_raw_rows)),
    "primary_loci_total": n_prim_total,
    "non_primary_rows_total": len(all_raw_rows) - n_prim_total,
    "primary_loci_eligible_tier1_primary": n_prim_t1,
    "primary_loci_ineligible_tier1_primary": n_prim_total - n_prim_t1,
    "primary_loci_eligible_tier2_stop3_sens": n_prim_t2,
    "primary_loci_ineligible_tier2_stop3_sens": n_prim_total - n_prim_t2,
    "primary_loci_eligible_tier3_partial_sens": n_prim_t3,
    "primary_loci_ineligible_tier3_partial_sens": n_prim_total - n_prim_t3,
    "flanks_total_evaluated": 4 * len(all_raw_rows),
    "flanks_all_rows_by_category": {
        "mapped_KO": flank_cats_all["mapped_KO"],
        "mapped_noKO": flank_cats_all["mapped_noKO"],
        "stop3prime_matched_KO": flank_cats_all["stop3prime_matched_KO"],
        "stop3prime_matched_noKO": flank_cats_all["stop3prime_matched_noKO"],
        "coordinate_mismatch_unmapped": flank_cats_all["coordinate_mismatch_unmapped"],
        "source_file_unavailable": flank_cats_all["source_file_unavailable"],
        "contig_edge_truncated": flank_cats_all["contig_edge_truncated"],
    },
    "flanks_primary_loci_by_category": {
        "mapped_KO": flank_cats_primary_loci["mapped_KO"],
        "mapped_noKO": flank_cats_primary_loci["mapped_noKO"],
        "stop3prime_matched_KO": flank_cats_primary_loci["stop3prime_matched_KO"],
        "stop3prime_matched_noKO": flank_cats_primary_loci["stop3prime_matched_noKO"],
        "coordinate_mismatch_unmapped": flank_cats_primary_loci["coordinate_mismatch_unmapped"],
        "source_file_unavailable": flank_cats_primary_loci["source_file_unavailable"],
        "contig_edge_truncated": flank_cats_primary_loci["contig_edge_truncated"],
    },
    "flanks_non_primary_rows_by_category": dict(flank_cats_non_primary),
    "flanks_tier1_eligible_primary_loci_by_category": dict(flank_cats_eligible_primary),
    "flanks_tier1_ineligible_primary_loci_by_category": dict(flank_cats_ineligible_primary),
    "flanks_tier2_eligible_stop3_loci_by_category": dict(flank_cats_eligible_stop3),
    "flanks_tier2_ineligible_stop3_loci_by_category": dict(flank_cats_ineligible_stop3),
    "flanks_annotated_product": n_flanks_product_all,
}

out_json = results_dir / "pilot_join_metrics.json"
out_json.write_text(json.dumps(metrics, indent=2) + "\n")

# Build complete executed_inputs_manifest.tsv and INPUTS.tsv with zero DYNAMIC placeholders
print("[flank:v3] Computing SHA-256 checksums for fixed upstream tables...", flush=True)
fixed_manifest = []
for role, p_obj in FIXED_INPUTS:
    if p_obj.exists():
        sz, sha = file_sha256_and_size(p_obj)
        fixed_manifest.append({
            "role": role,
            "path": str(p_obj),
            "size_bytes": sz,
            "sha256": sha,
            "status": "verified_read",
        })

# Deduplicate dynamic_manifest by (role, path) and sort deterministically
seen_paths = set()
all_manifest = []
for entry in fixed_manifest + sorted(dynamic_manifest, key=lambda x: (x["role"], x["path"])):
    key = (entry["role"], entry["path"])
    if key not in seen_paths:
        seen_paths.add(key)
        all_manifest.append(entry)

manifest_path = results_dir / "executed_inputs_manifest.tsv"
with open(manifest_path, "w") as f:
    w = csv.DictWriter(f, fieldnames=["role", "path", "size_bytes", "sha256", "status"], delimiter="\t")
    w.writeheader()
    w.writerows(all_manifest)

if update_top_inputs:
    top_inputs_path = results_dir.parent / "INPUTS.tsv"
    with open(top_inputs_path, "w") as f:
        w = csv.DictWriter(f, fieldnames=["role", "path", "size_bytes", "sha256", "status"], delimiter="\t")
        w.writeheader()
        w.writerows(all_manifest)
    print(f"[flank:v3] Wrote {len(all_manifest)} verified input hashes to {top_inputs_path}", flush=True)

print(f"[flank:v3] Wrote {len(all_raw_rows)} rows to {out_tsv}", flush=True)
print(f"[flank:v3] Wrote {len(contig_gene_rows)} contig genes to {out_genes_tsv}", flush=True)
print(f"[flank:v3] Wrote {len(all_manifest)} verified input hashes to {manifest_path}", flush=True)
