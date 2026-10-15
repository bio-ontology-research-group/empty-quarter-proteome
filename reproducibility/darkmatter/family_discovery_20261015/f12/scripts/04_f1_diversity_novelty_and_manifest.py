#!/usr/bin/env python3
"""
Step 4: Track F1 Diversity, Novelty Quantification, Clustering, and ProstT5 Manifest
- Aligns all CPF candidates (EQ, SMAG, GEM) against curated references via DIAMOND --ultra-sensitive.
- Evaluates sequence identity and coverage distributions relative to curated subfamilies.
- Performs MMseqs2 clustering at 50% and 80% sequence identity to quantify added environmental clusters.
- Prepares combined deduplicated FASTA and manifest for Root's reusable GPU ProstT5 runner.
- Assesses internal exploratory compartment distribution without fabricated p-values.
"""

import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path
from collections import defaultdict

ROOT = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12")
if not ROOT.exists():
    ROOT = Path("/home/leechuck/Documents/papers/empty-quarter-mags/reproducibility/darkmatter/family_discovery_20261015/f12")

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
OUT_DIR = ROOT / "results/full_v1"
OUT_DIR.mkdir(parents=True, exist_ok=True)

REF_TSV = ROOT / "references_v2/verified_references_metadata_dedup.tsv"
REF_FAA = ROOT / "references_v2/cpf_references_dedup.faa"
REF_DMND = ROOT / "references_v2/cpf_references_dedup.dmnd"

def load_references():
    """Load curated CPF reference metadata."""
    ref_map = {}
    with open(REF_TSV) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r["track"] == "F1_CPF":
                ref_map[r["accession"]] = r
    print(f"[OK] Loaded {len(ref_map)} curated CPF reference metadata records.")
    return ref_map

def load_gtdb_taxonomy():
    """Load GTDB-Tk R232 taxonomy for representative MAGs."""
    tax_map = {}
    gtdb_files = [
        W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv",
        W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
    ]
    for p in gtdb_files:
        if p.exists():
            with open(p) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for r in reader:
                    tax_map[r["user_genome"]] = r["classification"]
    print(f"[OK] Loaded taxonomy for {len(tax_map)} representative genomes.")
    return tax_map

def load_sample_compartments():
    """Load sample metadata for surface/deep compartment mapping."""
    comp_map = {}
    meta_p = W / "metadata/sample_metadata.tsv"
    if meta_p.exists():
        with open(meta_p) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                # header: sample, compartment (Surface, Deep, etc.)
                s = r.get("sample", r.get("sample_id", ""))
                c = r.get("compartment", r.get("Compartment", "Unknown"))
                if s:
                    comp_map[s] = c
    return comp_map

def align_candidates_vs_references(all_candidates_faa):
    """Run DIAMOND blastp --ultra-sensitive against curated CPF references."""
    out_tsv = OUT_DIR / "all_candidates_vs_cpf_refs.tsv"
    diamond_bin = "diamond"
    eggnog_bin = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")
    if eggnog_bin.exists():
        diamond_bin = str(eggnog_bin)
        
    cmd = [
        diamond_bin, "blastp",
        "--threads", os.environ.get("SLURM_CPUS_PER_TASK", "16"),
        "--query", str(all_candidates_faa),
        "--db", str(REF_DMND),
        "--out", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
        "--evalue", "1e-3",
        "--max-target-seqs", "10",
        "--ultra-sensitive"
    ]
    print(f"Running DIAMOND vs references: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    
    top_hits = {}
    with open(out_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if not parts or len(parts) < 14:
                continue
            qid = parts[0]
            bitscore = float(parts[11])
            pident = float(parts[2])
            aln_len = int(parts[3])
            qlen = int(parts[12])
            cov = aln_len / qlen if qlen > 0 else 0.0
            
            if qid not in top_hits or bitscore > top_hits[qid]["bitscore"]:
                top_hits[qid] = {
                    "ref_acc": parts[1],
                    "pident": pident,
                    "aln_len": aln_len,
                    "qlen": qlen,
                    "coverage": cov,
                    "evalue": float(parts[10]),
                    "bitscore": bitscore
                }
    print(f"[OK] DIAMOND alignment completed. {len(top_hits)} candidates aligned to references.")
    return top_hits

def run_clustering(combined_faa):
    """Run MMseqs2 clustering at 50% and 80% sequence identity."""
    tmp_dir = OUT_DIR / "mmseqs_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    threads = os.environ.get("SLURM_CPUS_PER_TASK", "16")
    
    # 50% ID / 80% cov
    out_prefix_50 = OUT_DIR / "cpf_cluster_id50"
    cmd_50 = [
        "mmseqs", "easy-cluster", str(combined_faa), str(out_prefix_50), str(tmp_dir),
        "--min-seq-id", "0.5", "-c", "0.8", "--cov-mode", "0", "--threads", threads
    ]
    print(f"Running MMseqs2 50% clustering: {' '.join(cmd_50)}")
    subprocess.run(cmd_50, check=True)
    
    # 80% ID / 80% cov
    out_prefix_80 = OUT_DIR / "cpf_cluster_id80"
    cmd_80 = [
        "mmseqs", "easy-cluster", str(combined_faa), str(out_prefix_80), str(tmp_dir),
        "--min-seq-id", "0.8", "-c", "0.8", "--cov-mode", "0", "--threads", threads
    ]
    print(f"Running MMseqs2 80% clustering: {' '.join(cmd_80)}")
    subprocess.run(cmd_80, check=True)
    
    # Parse 50% clusters
    clusters_50 = defaultdict(list)
    tsv_50 = OUT_DIR / "cpf_cluster_id50_cluster.tsv"
    if tsv_50.exists():
        with open(tsv_50) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 2:
                    clusters_50[p[0]].append(p[1])
                    
    # Parse 80% clusters
    clusters_80 = defaultdict(list)
    tsv_80 = OUT_DIR / "cpf_cluster_id80_cluster.tsv"
    if tsv_80.exists():
        with open(tsv_80) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) >= 2:
                    clusters_80[p[0]].append(p[1])
                    
    print(f"[OK] MMseqs2 complete: {len(clusters_50)} clusters at 50% ID; {len(clusters_80)} clusters at 80% ID.")
    return clusters_50, clusters_80

def main():
    print("=== Step 4: Track F1 Diversity, Novelty Quantification, Clustering, and Manifest ===")
    ref_map = load_references()
    tax_map = load_gtdb_taxonomy()
    comp_map = load_sample_compartments()
    
    # Load candidate lists from Step 2 (EQ) and Step 3 (SMAG, GEM)
    eq_tsv = OUT_DIR / "eq_cpf_candidates.tsv"
    smag_tsv = OUT_DIR / "smag_cpf_candidates.tsv"
    gem_tsv = OUT_DIR / "gem_cpf_candidates.tsv"
    
    eq_records = list(csv.DictReader(open(eq_tsv), delimiter="\t")) if eq_tsv.exists() else []
    smag_records = list(csv.DictReader(open(smag_tsv), delimiter="\t")) if smag_tsv.exists() else []
    gem_records = list(csv.DictReader(open(gem_tsv), delimiter="\t")) if gem_tsv.exists() else []
    
    print(f"Candidates loaded: EQ={len(eq_records)}, SMAG={len(smag_records)}, GEM={len(gem_records)}")
    
    # Build combined FASTA of all candidates
    all_cand_faa = OUT_DIR / "all_environmental_cpf_candidates.faa"
    with open(all_cand_faa, "w") as out:
        for p in [OUT_DIR / "eq_cpf_candidates.faa", OUT_DIR / "smag_cpf_candidates.faa", OUT_DIR / "gem_cpf_candidates.faa"]:
            if p.exists():
                out.write(p.read_text())
                
    # Run alignment vs references
    top_ref_hits = align_candidates_vs_references(all_cand_faa)
    
    # Build combined deduplicated FASTA and manifest including references
    combined_dedup_faa = OUT_DIR / "f1_combined_dedup.faa"
    combined_manifest_tsv = OUT_DIR / "f1_combined_manifest.tsv"
    
    seen_seqs = {}
    unique_seq_dict = {}
    manifest_rows = []
    
    def add_sequence(seq_id, source, seq_str, tax_str, cls_str, extra_info=""):
        seq_clean = seq_str.strip().upper().replace("*", "")
        h = hashlib.sha256(seq_clean.encode()).hexdigest()
        if h in seen_seqs:
            canonical_id = seen_seqs[h]
            return canonical_id
        seen_seqs[h] = seq_id
        unique_seq_dict[seq_id] = seq_clean
        manifest_rows.append({
            "id": seq_id,
            "source": source,
            "taxonomy": tax_str,
            "subfamily_class": cls_str,
            "length": len(seq_clean),
            "sha256": h,
            "extra_info": extra_info
        })
        return seq_id

    # 1. Add curated references
    cur_header = None
    cur_seq = []
    with open(REF_FAA) as f:
        for line in f:
            if line.startswith(">"):
                if cur_header:
                    acc = cur_header.split()[0]
                    r_meta = ref_map.get(acc, {})
                    add_sequence(acc, "Curated_Reference", "".join(cur_seq), r_meta.get("organism", "Reference"), r_meta.get("broad_family", "CPF"))
                cur_header = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_header:
            acc = cur_header.split()[0]
            r_meta = ref_map.get(acc, {})
            add_sequence(acc, "Curated_Reference", "".join(cur_seq), r_meta.get("organism", "Reference"), r_meta.get("broad_family", "CPF"))

    # 2. Add EQ candidates
    cur_header = None
    cur_seq = []
    with open(OUT_DIR / "eq_cpf_candidates.faa") as f:
        for line in f:
            if line.startswith(">"):
                if cur_header:
                    pid = cur_header.split()[0]
                    mag = pid.split("::")[0]
                    tax = tax_map.get(mag, "Unclassified")
                    hit = top_ref_hits.get(pid, {})
                    r_acc = hit.get("ref_acc", "")
                    r_meta = ref_map.get(r_acc, {})
                    cls = r_meta.get("broad_family", "Unclassified_CPF") if hit else "Divergent_CPF_candidate"
                    add_sequence(pid, "EmptyQuarter_MAG", "".join(cur_seq), tax, cls, f"top_ref={r_acc};pident={hit.get('pident','')}")
                cur_header = line[1:].strip()
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_header:
            pid = cur_header.split()[0]
            mag = pid.split("::")[0]
            tax = tax_map.get(mag, "Unclassified")
            hit = top_ref_hits.get(pid, {})
            r_acc = hit.get("ref_acc", "")
            r_meta = ref_map.get(r_acc, {})
            cls = r_meta.get("broad_family", "Unclassified_CPF") if hit else "Divergent_CPF_candidate"
            add_sequence(pid, "EmptyQuarter_MAG", "".join(cur_seq), tax, cls, f"top_ref={r_acc};pident={hit.get('pident','')}")

    # 3. Add SMAG candidates
    if (OUT_DIR / "smag_cpf_candidates.faa").exists():
        cur_header = None; cur_seq = []
        with open(OUT_DIR / "smag_cpf_candidates.faa") as f:
            for line in f:
                if line.startswith(">"):
                    if cur_header:
                        pid = cur_header.split()[0]
                        hit = top_ref_hits.get(pid, {})
                        r_acc = hit.get("ref_acc", "")
                        r_meta = ref_map.get(r_acc, {})
                        cls = r_meta.get("broad_family", "Unclassified_CPF") if hit else "Divergent_CPF_candidate"
                        add_sequence(pid, "SMAG_Marine", "".join(cur_seq), "Marine_plankton", cls, f"top_ref={r_acc};pident={hit.get('pident','')}")
                    cur_header = line[1:].strip(); cur_seq = []
                else:
                    cur_seq.append(line.strip())
            if cur_header:
                pid = cur_header.split()[0]
                hit = top_ref_hits.get(pid, {})
                r_acc = hit.get("ref_acc", "")
                r_meta = ref_map.get(r_acc, {})
                cls = r_meta.get("broad_family", "Unclassified_CPF") if hit else "Divergent_CPF_candidate"
                add_sequence(pid, "SMAG_Marine", "".join(cur_seq), "Marine_plankton", cls, f"top_ref={r_acc};pident={hit.get('pident','')}")

    # 4. Add GEM candidates
    if (OUT_DIR / "gem_cpf_candidates.faa").exists():
        cur_header = None; cur_seq = []
        with open(OUT_DIR / "gem_cpf_candidates.faa") as f:
            for line in f:
                if line.startswith(">"):
                    if cur_header:
                        pid = cur_header.split()[0]
                        hit = top_ref_hits.get(pid, {})
                        r_acc = hit.get("ref_acc", "")
                        r_meta = ref_map.get(r_acc, {})
                        cls = r_meta.get("broad_family", "Unclassified_CPF") if hit else "Divergent_CPF_candidate"
                        add_sequence(pid, "GEM_Global", "".join(cur_seq), "Global_environmental", cls, f"top_ref={r_acc};pident={hit.get('pident','')}")
                    cur_header = line[1:].strip(); cur_seq = []
                else:
                    cur_seq.append(line.strip())
            if cur_header:
                pid = cur_header.split()[0]
                hit = top_ref_hits.get(pid, {})
                r_acc = hit.get("ref_acc", "")
                r_meta = ref_map.get(r_acc, {})
                cls = r_meta.get("broad_family", "Unclassified_CPF") if hit else "Divergent_CPF_candidate"
                add_sequence(pid, "GEM_Global", "".join(cur_seq), "Global_environmental", cls, f"top_ref={r_acc};pident={hit.get('pident','')}")

    # Emit f1_combined_dedup.faa and f1_combined_manifest.tsv
    print(f"Total unique sequences in combined set: {len(manifest_rows)}")
    with open(combined_manifest_tsv, "w") as out:
        fields = ["id", "source", "taxonomy", "subfamily_class", "length", "sha256", "extra_info"]
        w = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(manifest_rows)

    # Write deduplicated FASTA for MMseqs2 and Root ProstT5 runner
    with open(combined_dedup_faa, "w") as out:
        for r in manifest_rows:
            seq_val = unique_seq_dict.get(r["id"], "")
            out.write(f">{r['id']}\n{seq_val}\n")

    c50, c80 = run_clustering(combined_dedup_faa)
    
    # Analyze clusters: how many EQ clusters are novel relative to references?
    ref_ids = set(ref_map.keys())
    eq_ids = {r["id"] for r in eq_records}
    smag_ids = {r["id"] for r in smag_records}
    gem_ids = {r["id"] for r in gem_records}
    
    eq_clusters_50 = 0
    eq_novel_clusters_50 = 0
    for rep, members in c50.items():
        m_set = set(members)
        if m_set & eq_ids:
            eq_clusters_50 += 1
            if not (m_set & ref_ids) and not (m_set & smag_ids) and not (m_set & gem_ids):
                eq_novel_clusters_50 += 1
                
    eq_clusters_80 = 0
    eq_novel_clusters_80 = 0
    for rep, members in c80.items():
        m_set = set(members)
        if m_set & eq_ids:
            eq_clusters_80 += 1
            if not (m_set & ref_ids) and not (m_set & smag_ids) and not (m_set & gem_ids):
                eq_novel_clusters_80 += 1

    # Internal exploratory compartment distribution of EQ candidates
    surface_count = 0
    deep_count = 0
    other_comp = 0
    for r in eq_records:
        mag = r["mag"]
        sample = mag.split("__")[0]
        c = comp_map.get(sample, "Unknown")
        if "Surface" in c or c == "S":
            surface_count += 1
        elif "Deep" in c or c == "D":
            deep_count += 1
        else:
            other_comp += 1

    diversity_metrics = {
        "curated_references_count": len(ref_map),
        "eq_candidates_count": len(eq_records),
        "smag_candidates_count": len(smag_records),
        "gem_candidates_count": len(gem_records),
        "total_unique_combined_proteins": len(manifest_rows),
        "mmseqs_50_clusters_total": len(c50),
        "mmseqs_50_eq_clusters": eq_clusters_50,
        "mmseqs_50_eq_novel_clusters": eq_novel_clusters_50,
        "mmseqs_80_clusters_total": len(c80),
        "mmseqs_80_eq_clusters": eq_clusters_80,
        "mmseqs_80_eq_novel_clusters": eq_novel_clusters_80,
        "exploratory_eq_compartment": {
            "surface_candidates": surface_count,
            "deep_candidates": deep_count,
            "other_or_unknown": other_comp,
            "note": "Representative-isolate compartment carriage is an internal exploratory metric. Community abundance and site-paired niche partitioning require sample metagenomic read depth profiles and are marked NOT_TESTED_IN_PILOT."
        }
    }
    
    with open(OUT_DIR / "f1_diversity_and_novelty_metrics.json", "w") as f:
        json.dump(diversity_metrics, f, indent=2)
        
    print(f"[SUCCESS] Diversity & novelty analysis complete.")
    print(f"  EQ clusters: {eq_clusters_50} (50% ID, {eq_novel_clusters_50} EQ-specific)")
    print(f"  ProstT5 combined manifest ready: {combined_manifest_tsv}")

if __name__ == "__main__":
    main()
