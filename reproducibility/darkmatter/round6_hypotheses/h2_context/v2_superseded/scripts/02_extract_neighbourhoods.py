#!/usr/bin/env python3
"""02_extract_neighbourhoods.py - Proven cross-caller neighbourhood extraction (v2 protocol).

Implements audited coordinate joins between Pyrodigal assembly proteins and Bakta/eggNOG:
1. Fast C-piped header extraction of assembly proteins via gzip -dc.
2. Contig sequence matching between raw bin FASTA and Bakta FNA:
   - Verifies 1:1 contig length correspondence.
   - Maps Bakta renamed contig_N back to raw_contig.
3. Fully-qualified annotation indexing by (raw_contig, start, stop, strand).
4. Strict join auditing into 5 mutually exclusive categories:
   - mapped_KO
   - mapped_noKO
   - coordinate_mismatch
   - source_file_unavailable
   - contig_edge_truncated
5. Emits:
   - pilot_neighbourhoods.tsv
   - pilot_contig_genes.tsv (for contig-preserving permutation null)
   - pilot_join_metrics.json
"""
import csv
import gzip
import json
import os
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
T = Path("/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites")
META_PATH = W / "metadata/sample_metadata.tsv"
ASM_PROT_DIR = W / "04_protein_catalog/assembly_proteins"

results_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results")
pilot_members_file = results_dir / "pilot_members.tsv"

if not pilot_members_file.exists():
    sys.exit(f"Error: {pilot_members_file} does not exist. Run 01_select_cohort.py first.")

# Read sample metadata
meta = {r["sample"]: r for r in csv.DictReader(open(META_PATH), delimiter="\t")}

# Load pilot members
members = []
samples_needed = set()
with open(pilot_members_file) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        members.append(row)
        samples_needed.add(row["sample"])

print(f"[flank] Processing {len(members)} members across {len(samples_needed)} samples (8 worker pool)")

# Group members by sample and contig
sample_contig_queries = defaultdict(lambda: defaultdict(list))
for m in members:
    sample_contig_queries[m["sample"]][m["contig"]].append(m)

# Worker function to extract contig CDSs for a single sample
def process_single_sample(args):
    sample, contig_map = args
    asm_path = ASM_PROT_DIR / sample / f"{sample}.faa.gz"
    if not asm_path.exists():
        return [], 0, 0, {}
    
    needed_contigs = set(contig_map.keys())
    genes_by_contig = defaultdict(list)
    
    try:
        proc = subprocess.Popen(["gzip", "-dc", str(asm_path)], stdout=subprocess.PIPE, text=True, bufsize=1024*1024)
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
                    "partial": partial
                })
        proc.wait()
    except Exception as e:
        return [], 0, 0, {}
    
    rows = []
    matched = 0
    contig_cds_records = {}
    
    for contig, q_list in contig_map.items():
        c_genes = sorted(genes_by_contig[contig], key=lambda g: g["start"])
        contig_cds_records[contig] = c_genes
        pid_to_idx = {g["pid"]: i for i, g in enumerate(c_genes)}
        
        for q in q_list:
            focal_pid = q["member"]
            if focal_pid not in pid_to_idx:
                continue
            
            matched += 1
            idx = pid_to_idx[focal_pid]
            focal_gene = c_genes[idx]
            
            row = {
                "family": q["family"],
                "chosen_rep": q["chosen_rep"],
                "arm": q["arm"],
                "member": focal_pid,
                "sample": sample,
                "contig": contig,
                "contig_len": q["contig_len"],
                "bin": q["bin"],
                "species": q["species"],
                "taxonomy": q["taxonomy"],
                "is_primary_locus": q["is_primary_locus"],
                "locus_copy_count": q["locus_copy_count"],
                "focal_start": focal_gene["start"],
                "focal_stop": focal_gene["stop"],
                "focal_strand": focal_gene["strand"],
                "focal_partial": focal_gene["partial"],
                "total_genes_on_contig": len(c_genes),
                "focal_gene_index": idx
            }
            
            for offset in [-2, -1, 1, 2]:
                flank_key = f"m{abs(offset)}" if offset < 0 else f"p{offset}"
                target_idx = idx + offset
                if 0 <= target_idx < len(c_genes):
                    tg = c_genes[target_idx]
                    dist = focal_gene["start"] - tg["stop"] - 1 if offset < 0 else tg["start"] - focal_gene["stop"] - 1
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
            
    return rows, matched, 1, contig_cds_records

# Execute extraction across all samples in parallel
t0 = time.time()
sample_items = list(sample_contig_queries.items())
all_raw_rows = []
total_matched = 0
all_contig_cds = {}

with ProcessPoolExecutor(max_workers=8) as executor:
    for rows, matched, done_sample, c_records in executor.map(process_single_sample, sample_items):
        all_raw_rows.extend(rows)
        total_matched += matched
        all_contig_cds.update(c_records)

t1 = time.time()
print(f"[flank] Extracted {len(all_raw_rows)} raw neighbourhood rows ({total_matched}/{len(members)} matched) in {t1 - t0:.1f}s", flush=True)

def normalize_contig(name):
    return name.split("::")[-1] if name else ""

# Helper to parse Bakta TSV and eggNOG annotations with fail-closed contig mapping
bin_ann_cache = {}
bin_contig_map_cache = {}

def get_bin_annotations(bin_name):
    if not bin_name:
        return {}
    if bin_name in bin_ann_cache:
        return bin_ann_cache[bin_name]
    
    s = bin_name.split("__")[0]
    m = meta.get(s)
    if not m:
        bin_ann_cache[bin_name] = {}
        return {}
    
    rest = bin_name.split("__", 1)[1]
    if rest.startswith("single_"):
        modes = ["single_sample"]
        local = rest[len("single_"):]
    elif rest.startswith("self_"):
        modes = ["self_supervised"]
        local = rest[len("self_"):]
    else:
        modes = ["self_supervised", "single_sample"]
        local = rest
    
    parent = local[:-4] if local.endswith("_sub") else local
    trip = m["trip"]
    site = m["site"]
    base = T / trip / f"site{site}" / s / "wgs"
    
    bin_fa, bakta_fna, bakta_tsv, eggnog_ann = None, None, None, None
    for mode in modes:
        b_fa = base / "binning" / mode / "output_bins" / f"{parent}.fa.gz"
        b_fna = base / "annotation" / mode / "bakta" / parent / f"{parent}.fna"
        b_tsv = base / "annotation" / mode / "bakta" / parent / f"{parent}.tsv"
        e_ann = base / "annotation" / mode / "eggnog" / parent / f"{parent}.emapper.annotations"
        
        if b_fa.exists() and b_fna.exists() and b_tsv.exists():
            bin_fa = b_fa
            bakta_fna = b_fna
            bakta_tsv = b_tsv
            eggnog_ann = e_ann if e_ann.exists() else None
            break
            
    if not bakta_tsv:
        bin_ann_cache[bin_name] = {}
        return {}
        
    # Build 1:1 contig mapping from raw bin to Bakta
    raw_contigs = []
    try:
        proc = subprocess.Popen(["gzip", "-dc", str(bin_fa)], stdout=subprocess.PIPE, text=True)
        cur_hdr, cur_len = None, 0
        for line in proc.stdout:
            if line.startswith(">"):
                if cur_hdr:
                    raw_contigs.append((cur_hdr, cur_len))
                cur_hdr = line[1:].strip().split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_hdr:
            raw_contigs.append((cur_hdr, cur_len))
        proc.wait()
    except Exception:
        bin_ann_cache[bin_name] = {}
        return {}
        
    bakta_contigs = []
    with open(bakta_fna) as f:
        cur_hdr, cur_len = None, 0
        for line in f:
            if line.startswith(">"):
                if cur_hdr:
                    bakta_contigs.append((cur_hdr, cur_len))
                cur_hdr = line[1:].strip().split()[0]
                cur_len = 0
            else:
                cur_len += len(line.strip())
        if cur_hdr:
            bakta_contigs.append((cur_hdr, cur_len))
            
    if len(raw_contigs) != len(bakta_contigs):
        print(f"[flank:warn] Contig count mismatch for bin {bin_name}: raw={len(raw_contigs)} bakta={len(bakta_contigs)}")
        bin_ann_cache[bin_name] = {}
        return {}
        
    bakta_to_raw = {}
    for i, ((r_hdr, r_len), (b_hdr, b_len)) in enumerate(zip(raw_contigs, bakta_contigs)):
        if r_len != b_len:
            print(f"[flank:warn] Contig length mismatch for bin {bin_name} idx {i}: {r_len} != {b_len}")
            bin_ann_cache[bin_name] = {}
            return {}
        bakta_to_raw[b_hdr] = r_hdr
        
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
                ko_clean = ",".join(k.replace("ko:", "") for k in ko.split(",")) if ko and ko != "-" else ""
                eggnog_data[locus] = {
                    "ko": ko_clean,
                    "cog": row.get("COG_category", ""),
                    "desc": row.get("Description", ""),
                    "gene": row.get("Preferred_name", "")
                }
                
    # Read Bakta TSV and map to raw contig coordinates
    ann_by_contig_coord = {}
    with open(bakta_tsv) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) >= 8 and c[1].lower() == "cds":
                b_contig = c[0]
                r_contig = bakta_to_raw.get(b_contig)
                if not r_contig:
                    continue
                start = int(c[2])
                stop = int(c[3])
                strand = 1 if c[4] == "+" else -1
                locus = c[5]
                prod = c[7] if len(c) > 7 else ""
                egg = eggnog_data.get(locus, {})
                
                # Full key: (normalized raw_contig, start, stop, strand)
                ann_by_contig_coord[(normalize_contig(r_contig), start, stop, strand)] = {
                    "locus": locus,
                    "product": prod,
                    "ko": egg.get("ko", ""),
                    "cog": egg.get("cog", ""),
                    "desc": egg.get("desc", ""),
                    "gene": egg.get("gene", "")
                }
                
    bin_ann_cache[bin_name] = ann_by_contig_coord
    return ann_by_contig_coord

print("[flank] Joining functional annotations with verified contig keys...", flush=True)
unique_bins = {r["bin"] for r in all_raw_rows if r["bin"]}
print(f"[flank] Resolving annotations for {len(unique_bins)} unique bins...", flush=True)

for b in unique_bins:
    get_bin_annotations(b)

metrics = {
    "total_queries": len(members),
    "exact_focal_coordinate_matches": total_matched,
    "binned_queries": sum(1 for r in all_raw_rows if r["bin"]),
    "flanks_total_evaluated": 0,
    "flanks_edge_truncated": 0,
    "flanks_source_file_unavailable": 0,
    "flanks_coordinate_mismatch": 0,
    "flanks_mapped_noKO": 0,
    "flanks_mapped_KO": 0,
    "flanks_annotated_product": 0,
    "samples_processed": len(samples_needed)
}

final_rows = []
for r in all_raw_rows:
    b = r["bin"]
    c = r["contig"]
    norm_c = normalize_contig(c)
    bin_ann = bin_ann_cache.get(b, {}) if b else {}
    
    # Focal gene annotation
    f_key = (norm_c, r["focal_start"], r["focal_stop"], r["focal_strand"])
    f_ann = bin_ann.get(f_key, {})
    r["focal_ko"] = f_ann.get("ko", "")
    r["focal_product"] = f_ann.get("product", "")
    
    # Flank annotations and audit categories
    for flank_key in ["m2", "m1", "p1", "p2"]:
        metrics["flanks_total_evaluated"] += 1
        if r[f"{flank_key}_truncated"] == 1:
            metrics["flanks_edge_truncated"] += 1
            audit_cat = "contig_edge_truncated"
            ko = ""
            cog = ""
            desc = ""
            prod = "CONTIG_EDGE"
        elif not b or not bin_ann:
            metrics["flanks_source_file_unavailable"] += 1
            audit_cat = "source_file_unavailable"
            ko = ""
            cog = ""
            desc = ""
            prod = "UNAVAILABLE"
        else:
            flank_coord_key = (norm_c, r[f"{flank_key}_start"], r[f"{flank_key}_stop"], r[f"{flank_key}_strand"])
            if flank_coord_key in bin_ann:
                t_ann = bin_ann[flank_coord_key]
                ko = t_ann.get("ko", "")
                cog = t_ann.get("cog", "")
                desc = t_ann.get("desc", "")
                prod = t_ann.get("product", "")
                if ko:
                    metrics["flanks_mapped_KO"] += 1
                    audit_cat = "mapped_KO"
                else:
                    metrics["flanks_mapped_noKO"] += 1
                    audit_cat = "mapped_noKO"
                if prod:
                    metrics["flanks_annotated_product"] += 1
            else:
                metrics["flanks_coordinate_mismatch"] += 1
                audit_cat = "coordinate_mismatch"
                ko = ""
                cog = ""
                desc = ""
                prod = "COORD_MISMATCH"
                
        r[f"{flank_key}_audit_category"] = audit_cat
        r[f"{flank_key}_ko"] = ko
        r[f"{flank_key}_cog"] = cog
        r[f"{flank_key}_desc"] = desc
        r[f"{flank_key}_product"] = prod
        
        # Clean up temporary coord fields
        del r[f"{flank_key}_start"]
        del r[f"{flank_key}_stop"]
        del r[f"{flank_key}_strand"]
        
    final_rows.append(r)

out_tsv = results_dir / "pilot_neighbourhoods.tsv"
with open(out_tsv, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(final_rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(final_rows)

# Generate pilot_contig_genes.tsv for host/contig-preserving permutation null
contig_gene_rows = []
for contig, c_genes in sorted(all_contig_cds.items()):
    sample = contig.split("_")[0] if "_" in contig else ""
    # Find bin for this contig if any
    c_bin = ""
    for r in all_raw_rows:
        if r["contig"] == contig:
            c_bin = r["bin"]
            break
    bin_ann = bin_ann_cache.get(c_bin, {}) if c_bin else {}
    
    norm_contig = normalize_contig(contig)
    for idx, g in enumerate(c_genes):
        g_key = (norm_contig, g["start"], g["stop"], g["strand"])
        ann = bin_ann.get(g_key, {})
        ko = ann.get("ko", "")
        prod = ann.get("product", "")
        desc = ann.get("desc", "")
        
        if not c_bin or not bin_ann:
            audit = "source_file_unavailable"
        elif g_key in bin_ann:
            audit = "mapped_KO" if ko else "mapped_noKO"
        else:
            audit = "coordinate_mismatch"
            
        contig_gene_rows.append({
            "contig": contig,
            "bin": c_bin,
            "gene_index": idx,
            "pid": g["pid"],
            "start": g["start"],
            "stop": g["stop"],
            "strand": g["strand"],
            "ko": ko,
            "product": prod,
            "desc": desc,
            "audit_category": audit
        })

out_genes_tsv = results_dir / "pilot_contig_genes.tsv"
with open(out_genes_tsv, "w") as f:
    w = csv.DictWriter(f, fieldnames=list(contig_gene_rows[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(contig_gene_rows)

out_json = results_dir / "pilot_join_metrics.json"
out_json.write_text(json.dumps(metrics, indent=2) + "\n")

print(f"[flank] Successfully wrote {len(final_rows)} rows to {out_tsv}")
print(f"[flank] Successfully wrote {len(contig_gene_rows)} contig genes to {out_genes_tsv}")
print(f"[flank] Join metrics: {json.dumps(metrics, indent=2)}")
