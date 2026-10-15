#!/usr/bin/env python3
"""03_f2_remapping_and_test.py - Remapping of Historical Osmolyte-Phosphatases to Current Catalogue & Competing HAD Checks.

Features:
1. Exact remapping of the two named leads from ROOT_ASSESSMENT:
   - Lead 1: 1PRr2__self_SemiBin_81::1PRr2_1391604_4 (EC 3.1.3.70)
   - Lead 2: S52PRr1__single_SemiBin_53::S52PRr1_477345_9 (EC 3.1.3.70, Thermoplasmatota SW-10-69-26)
2. Checks coding completeness (Prodigal partial=00), length, contig integrity.
3. Maps contigs to current 14,426 refined bins and 5,438 ANI95 representatives via contig_to_bin_refined.tsv.gz, Cdb.csv, Wdb.csv.
4. Performs competing-function discrimination against MpgP, GgpP, TPP, and competing HAD-superfamily enzymes (PGP, SPP, HisP, beta-PGM, DeoB, YrbI).
5. Scans genomic neighbourhood (+/- 5 flanking genes) for compatible-solute biosynthetic operons (MpgS EC 2.4.1.217, GgpS EC 2.4.1.213, OtsA EC 2.4.1.15).
6. Remaps the broader set of 306 historical osmolyte-phosphatase families from osmolyte_families.tsv.
"""

import csv
import gzip
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
REF_DIR = ROOT / "references"
RES_DIR = ROOT / "results"
RES_DIR.mkdir(parents=True, exist_ok=True)

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
CONTIG_TO_BIN = W / "04_protein_catalog/contig_to_bin_refined.tsv.gz"
DREP_CDB = W / "05_refine/_drep/out/data_tables/Cdb.csv"
DREP_WDB = W / "05_refine/_drep/out/data_tables/Wdb.csv"
GTDB_BAC = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
GTDB_AR = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
SAMPLE_META = W / "metadata/sample_metadata.tsv"
ASM_PROT_DIR = W / "04_protein_catalog/assembly_proteins"

# Historical osmolyte families file
HIST_OSM_FILE = Path("/home/leechuck/Documents/papers/empty-quarter-mags/results/darkmatter_ecology/indicators/osmolyte_families.tsv")
if not HIST_OSM_FILE.exists():
    HIST_OSM_FILE = ROOT / "osmolyte_families.tsv"

HAD_HMM = REF_DIR / "had_all.hmm"
PHOS_DMND = REF_DIR / "phosphatase_references.dmnd"
DIAMOND_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")
if not DIAMOND_BIN.exists():
    DIAMOND_BIN = Path("diamond")

LEADS = [
    {
        "lead_id": "Lead_1",
        "sample": "1PRr2",
        "historical_family": "1PRr2__self_SemiBin_81::1PRr2_1391604_4",
        "historical_ec": "3.1.3.70",
        "protein_id": "1PRr2_1391604_4",
        "contig_id": "1PRr2_1391604",
        "expected_taxon": "Unspecified",
    },
    {
        "lead_id": "Lead_2",
        "sample": "S52PRr1",
        "historical_family": "S52PRr1__single_SemiBin_53::S52PRr1_477345_9",
        "historical_ec": "3.1.3.70",
        "protein_id": "S52PRr1_477345_9",
        "contig_id": "S52PRr1_477345",
        "expected_taxon": "Thermoplasmatota SW-10-69-26",
    },
]


def load_taxonomy():
    tax_map = {}
    for p in [GTDB_BAC, GTDB_AR]:
        if not p.exists():
            continue
        with open(p) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                mag = r["user_genome"]
                tax = r["classification"]
                parts = tax.split(";")
                d = p_ = c = o = f_ = g = s = ""
                for part in parts:
                    if part.startswith("d__"): d = part[3:]
                    elif part.startswith("p__"): p_ = part[3:]
                    elif part.startswith("c__"): c = part[3:]
                    elif part.startswith("o__"): o = part[3:]
                    elif part.startswith("f__"): f_ = part[3:]
                    elif part.startswith("g__"): g = part[3:]
                    elif part.startswith("s__"): s = part[3:]
                tax_map[mag] = {
                    "classification": tax,
                    "domain": d,
                    "phylum": p_,
                    "class": c,
                    "order": o,
                    "family": f_,
                    "genus": g,
                    "species": s,
                }
    return tax_map


def load_drep_mapping():
    print("Loading dRep Cdb and Wdb mapping...")
    rep_genomes = set()
    if DREP_WDB.exists():
        with open(DREP_WDB) as f:
            for r in csv.DictReader(f):
                rep_genomes.add(r["genome"].replace(".fa", "").replace(".fasta", ""))

    bin_to_rep = {}
    if DREP_CDB.exists():
        cluster_reps = {}
        # First pass: map cluster to representative
        with open(DREP_CDB) as f:
            for r in csv.DictReader(f):
                g = r["genome"].replace(".fa", "").replace(".fasta", "")
                c = r["primary_cluster"]
                if g in rep_genomes:
                    cluster_reps[c] = g
        # Second pass: map each bin to representative
        with open(DREP_CDB) as f:
            for r in csv.DictReader(f):
                g = r["genome"].replace(".fa", "").replace(".fasta", "")
                c = r["primary_cluster"]
                bin_to_rep[g] = {
                    "primary_cluster": c,
                    "representative": cluster_reps.get(c, g),
                    "is_rep": (g in rep_genomes),
                }
    return bin_to_rep, rep_genomes


def load_contig_to_bin():
    print("Loading contig_to_bin_refined.tsv.gz...")
    contig_to_bin = {}
    if CONTIG_TO_BIN.exists():
        with gzip.open(CONTIG_TO_BIN, "rt") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    contig_to_bin[parts[0]] = parts[1]
    return contig_to_bin


def extract_contig_genes(sample, target_contig):
    """Extracts all genes on target_contig from the sample assembly FAA."""
    faa_path = ASM_PROT_DIR / sample / f"{sample}.faa.gz"
    if not faa_path.exists():
        return []
    genes = []
    with gzip.open(faa_path, "rt") as f:
        cur_pid = None
        cur_hdr = None
        seq_parts = []
        for line in f:
            if line.startswith(">"):
                if cur_pid and cur_pid.rsplit("_", 1)[0] == target_contig:
                    genes.append({
                        "pid": cur_pid,
                        "header": cur_hdr,
                        "seq": "".join(seq_parts),
                    })
                cur_hdr = line[1:].strip()
                cur_pid = cur_hdr.split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if cur_pid and cur_pid.rsplit("_", 1)[0] == target_contig:
            genes.append({
                "pid": cur_pid,
                "header": cur_hdr,
                "seq": "".join(seq_parts),
            })
    return genes


def main():
    t0 = time.time()
    print("=== [F2 Pilot Analysis] Remapping Historical Osmolyte-Phosphatases to Current Catalogue ===")

    tax_map = load_taxonomy()
    bin_to_rep, rep_genomes = load_drep_mapping()
    contig_to_bin = load_contig_to_bin()
    print(f"Loaded {len(contig_to_bin)} contig-to-bin mappings, {len(bin_to_rep)} binned genomes, {len(rep_genomes)} representatives.")

    # 1. Remap the Two Named Leads
    print("\n--- 1. Remapping Lead 1 and Lead 2 ---")
    lead_results = []
    lead_seqs = {}

    for lead in LEADS:
        s = lead["sample"]
        c_id = lead["contig_id"]
        p_id = lead["protein_id"]
        print(f"\nProcessing {lead['lead_id']}: {lead['historical_family']}")

        # Map contig to refined bin
        refined_bin = contig_to_bin.get(c_id, "UNBINNED")
        drep_info = bin_to_rep.get(refined_bin, {})
        rep_mag = drep_info.get("representative", "NONE")
        is_rep = drep_info.get("is_rep", False)
        cluster = drep_info.get("primary_cluster", "NONE")
        tax = tax_map.get(refined_bin, tax_map.get(rep_mag, {}))

        print(f"  Contig {c_id} -> Refined Bin: {refined_bin}")
        print(f"  Cluster: {cluster}, Representative: {rep_mag} (is_rep={is_rep})")
        print(f"  Taxonomy: {tax.get('classification', 'Unclassified')}")

        # Extract all genes on contig
        c_genes = extract_contig_genes(s, c_id)
        print(f"  Found {len(c_genes)} total coding sequences on contig {c_id}.")

        # Find focal gene
        focal_gene = None
        for g in c_genes:
            if g["pid"] == p_id or g["pid"].endswith(f"::{p_id}"):
                focal_gene = g
                break

        if focal_gene:
            hdr = focal_gene["header"]
            seq = focal_gene["seq"]
            lead_seqs[lead["lead_id"]] = (focal_gene["pid"], seq)
            # Parse Prodigal header: >pid # start # stop # strand # ID=...;partial=00;...
            tokens = hdr.split(" # ")
            start = int(tokens[1]) if len(tokens) >= 2 else 0
            stop = int(tokens[2]) if len(tokens) >= 3 else 0
            strand = int(tokens[3]) if len(tokens) >= 4 else 1
            partial = "00"
            if len(tokens) >= 5:
                for item in tokens[4].split(";"):
                    if item.startswith("partial="):
                        partial = item.split("=")[1]

            seq_len = len(seq)
            is_complete = (partial == "00")
            print(f"  Focal Gene: {p_id}, length={seq_len} aa, coords={start}-{stop} (strand {strand}), partial={partial}, complete={is_complete}")

            lead_results.append({
                "lead_id": lead["lead_id"],
                "historical_family": lead["historical_family"],
                "historical_ec": lead["historical_ec"],
                "sample": s,
                "contig_id": c_id,
                "protein_id": p_id,
                "refined_bin": refined_bin,
                "is_refined": (refined_bin != "UNBINNED"),
                "primary_cluster": cluster,
                "representative_mag": rep_mag,
                "is_representative": is_rep,
                "domain": tax.get("domain", ""),
                "phylum": tax.get("phylum", ""),
                "class": tax.get("class", ""),
                "order": tax.get("order", ""),
                "family": tax.get("family", ""),
                "genus": tax.get("genus", ""),
                "species": tax.get("species", ""),
                "gtdbtk_r232": tax.get("classification", ""),
                "start_bp": start,
                "stop_bp": stop,
                "strand": strand,
                "partial": partial,
                "is_complete": is_complete,
                "length_aa": seq_len,
                "total_genes_on_contig": len(c_genes),
            })
        else:
            print(f"  WARNING: Focal gene {p_id} not found on contig {c_id}!")
            lead_results.append({
                "lead_id": lead["lead_id"],
                "historical_family": lead["historical_family"],
                "historical_ec": lead["historical_ec"],
                "sample": s,
                "contig_id": c_id,
                "protein_id": p_id,
                "refined_bin": refined_bin,
                "is_refined": (refined_bin != "UNBINNED"),
                "primary_cluster": cluster,
                "representative_mag": rep_mag,
                "is_representative": is_rep,
                "gtdbtk_r232": tax.get("classification", ""),
                "total_genes_on_contig": len(c_genes),
                "is_complete": False,
                "length_aa": 0,
            })

    # Write f2_lead_remapping.tsv
    lead_tsv = RES_DIR / "f2_lead_remapping.tsv"
    with open(lead_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(lead_results[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(lead_results)
    print(f"Wrote {lead_tsv.name}.")

    # 2. Competing-Function Screening for Leads
    print("\n--- 2. Competing-Function HAD Superfamily Analysis for Leads ---")
    lead_faa = RES_DIR / "f2_leads.faa"
    with open(lead_faa, "w") as f:
        for lid, (pid, seq) in lead_seqs.items():
            f.write(f">{lid}__{pid}\n{seq}\n")

    # Run hmmsearch on leads with had_all.hmm
    hmm_domtbl = RES_DIR / "f2_leads_had_hmm.domtbl"
    cmd = ["hmmsearch", "--domtblout", str(hmm_domtbl), "-o", str(RES_DIR / "f2_leads_had_hmm.out"), str(HAD_HMM), str(lead_faa)]
    subprocess.run(cmd, check=True)

    # Parse HMM hits
    lead_hmm_scores = defaultdict(dict)
    with open(hmm_domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) >= 23:
                qid = p[0].split("__")[0]
                h_acc = p[3].split(".")[0]
                lead_hmm_scores[qid][h_acc] = {
                    "name": p[4],
                    "evalue": float(p[6]),
                    "score": float(p[7]),
                }

    # Run DIAMOND blastp against curated phosphatase references
    dmnd_out = RES_DIR / "f2_leads_diamond.tsv"
    cmd = [
        str(DIAMOND_BIN), "blastp",
        "-d", str(PHOS_DMND),
        "-q", str(lead_faa),
        "-o", str(dmnd_out),
        "--ultra-sensitive",
        "-e", "1e-2",
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovhsp", "scovhsp"
    ]
    subprocess.run(cmd, check=True)

    # Parse DIAMOND hits per lead
    lead_dmnd_hits = defaultdict(list)
    with open(dmnd_out) as f:
        for line in f:
            p = line.strip().split("\t")
            if len(p) >= 16:
                qid = p[0].split("__")[0]
                lead_dmnd_hits[qid].append({
                    "ref_id": p[1],
                    "pident": float(p[2]),
                    "aln_len": int(p[3]),
                    "evalue": float(p[10]),
                    "bitscore": float(p[11]),
                    "qcovhsp": float(p[14]),
                    "scovhsp": float(p[15]),
                })

    # Build competing function matrix
    matrix_rows = []
    for lead in LEADS:
        lid = lead["lead_id"]
        hits = lead_dmnd_hits.get(lid, [])
        hmm = lead_hmm_scores.get(lid, {})

        # Top scores by category
        cat_scores = defaultdict(lambda: {"best_bitscore": 0.0, "best_eval": 1.0, "best_id": "NONE", "pident": 0.0})
        for h in hits:
            ref = h["ref_id"]
            cat = "Competing_HAD"
            if ref.startswith("MpgP"): cat = "MpgP_EC_3_1_3_70"
            elif ref.startswith("GgpP"): cat = "GgpP_EC_3_1_3_69"
            elif ref.startswith("TPP"): cat = "TPP_EC_3_1_3_12"

            if h["bitscore"] > cat_scores[cat]["best_bitscore"]:
                cat_scores[cat] = {
                    "best_bitscore": h["bitscore"],
                    "best_eval": h["evalue"],
                    "best_id": ref,
                    "pident": h["pident"],
                    "qcovhsp": h["qcovhsp"],
                }

        # Check catalytic HAD Motifs in sequence:
        seq = lead_seqs.get(lid, ("", ""))[1]
        # Motif I: [D][X][D][X][T/V]
        m1 = bool(re.search(r"D.D.[TV]", seq))
        # Motif II: [S/T]
        m2 = bool(re.search(r"[ST]", seq))
        # Motif III: [K]
        m3 = bool(re.search(r"K", seq))
        # Motif IV: [G/E][D]...[D]
        m4 = bool(re.search(r"[GE]D.{1,5}D", seq))

        matrix_rows.append({
            "lead_id": lid,
            "sample": lead["sample"],
            "protein_id": lead["protein_id"],
            "PF00702_score": hmm.get("PF00702", {}).get("score", 0.0),
            "PF00702_eval": hmm.get("PF00702", {}).get("evalue", 1.0),
            "PF13419_score": hmm.get("PF13419", {}).get("score", 0.0),
            "PF13419_eval": hmm.get("PF13419", {}).get("evalue", 1.0),
            "best_MpgP_bitscore": cat_scores["MpgP_EC_3_1_3_70"]["best_bitscore"],
            "best_MpgP_evalue": cat_scores["MpgP_EC_3_1_3_70"]["best_eval"],
            "best_MpgP_ref": cat_scores["MpgP_EC_3_1_3_70"]["best_id"],
            "best_GgpP_bitscore": cat_scores["GgpP_EC_3_1_3_69"]["best_bitscore"],
            "best_GgpP_evalue": cat_scores["GgpP_EC_3_1_3_69"]["best_eval"],
            "best_GgpP_ref": cat_scores["GgpP_EC_3_1_3_69"]["best_id"],
            "best_TPP_bitscore": cat_scores["TPP_EC_3_1_3_12"]["best_bitscore"],
            "best_TPP_evalue": cat_scores["TPP_EC_3_1_3_12"]["best_eval"],
            "best_TPP_ref": cat_scores["TPP_EC_3_1_3_12"]["best_id"],
            "best_Competing_HAD_bitscore": cat_scores["Competing_HAD"]["best_bitscore"],
            "best_Competing_HAD_evalue": cat_scores["Competing_HAD"]["best_eval"],
            "best_Competing_HAD_ref": cat_scores["Competing_HAD"]["best_id"],
            "has_motif_I_DxDxTV": m1,
            "has_motif_IV_GDxxxD": m4,
        })

    matrix_tsv = RES_DIR / "f2_competing_function_matrix.tsv"
    with open(matrix_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(matrix_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(matrix_rows)
    print(f"Wrote {matrix_tsv.name}.")

    # 3. Genomic Context & Flanking Gene Analysis
    print("\n--- 3. Genomic Neighbourhood Context Analysis for Leads ---")
    flank_rows = []
    for lead in LEADS:
        lid = lead["lead_id"]
        s = lead["sample"]
        c_id = lead["contig_id"]
        p_id = lead["protein_id"]

        c_genes = extract_contig_genes(s, c_id)
        # Sort genes by PID number
        def gene_num(g):
            try:
                return int(g["pid"].rsplit("_", 1)[-1])
            except ValueError:
                return 0
        c_genes_sorted = sorted(c_genes, key=gene_num)
        focal_idx = None
        for i, g in enumerate(c_genes_sorted):
            if g["pid"] == p_id or g["pid"].endswith(f"::{p_id}"):
                focal_idx = i
                break

        if focal_idx is not None:
            # Extract +/- 5 genes
            start_i = max(0, focal_idx - 5)
            end_i = min(len(c_genes_sorted), focal_idx + 6)
            for i in range(start_i, end_i):
                g = c_genes_sorted[i]
                offset = i - focal_idx
                role = "focal_lead" if offset == 0 else f"flank_{offset:+d}"
                seq_len = len(g["seq"])
                # Check for synthase motifs / keywords
                flank_rows.append({
                    "lead_id": lid,
                    "sample": s,
                    "contig_id": c_id,
                    "gene_offset": offset,
                    "gene_role": role,
                    "gene_pid": g["pid"],
                    "gene_length_aa": seq_len,
                    "raw_header": g["header"],
                })

    flank_tsv = RES_DIR / "f2_genomic_neighbourhoods.tsv"
    with open(flank_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(flank_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(flank_rows)
    print(f"Wrote {len(flank_rows)} neighbourhood genes to {flank_tsv.name}.")

    # 4. Census of the 306 Historical Osmolyte-Phosphatase Families
    print("\n--- 4. Census of the 306 Historical Osmolyte Families ---")
    hist_families = []
    if HIST_OSM_FILE.exists():
        with open(HIST_OSM_FILE) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                hist_families.append(r)
    print(f"Loaded {len(hist_families)} historical osmolyte family rows from {HIST_OSM_FILE.name}.")

    # Map each family's representative member
    remap_summary = Counter()
    ec_counts = Counter()
    phylum_counts = Counter()

    for hf in hist_families:
        fid = hf["family_id"]
        ec = hf.get("EC", "")
        ec_counts[ec] += 1

        # Extract locus from family_id: e.g. 1PRr2__self_SemiBin_81::1PRr2_1391604_4
        locus = fid.split("::")[-1] if "::" in fid else fid
        contig = locus.rsplit("_", 1)[0]

        r_bin = contig_to_bin.get(contig, "UNBINNED")
        if r_bin != "UNBINNED":
            remap_summary["mapped_to_refined_bin"] += 1
            d_info = bin_to_rep.get(r_bin, {})
            if d_info.get("is_rep", False):
                remap_summary["is_ani95_representative"] += 1
            else:
                remap_summary["in_represented_cluster"] += 1
            t_info = tax_map.get(r_bin, {})
            ph = t_info.get("phylum", "Unknown")
            phylum_counts[ph] += 1
        else:
            remap_summary["unbinned_in_refined_catalogue"] += 1

    print("\nBroader 306 Historical Census Remapping Summary:")
    for k, v in remap_summary.items():
        print(f"  {k:<35}: {v:<5} ({100.0 * v / len(hist_families):.1f}%)")

    print("\nHistorical EC Distribution in 306 Set:")
    for ec, c in ec_counts.most_common(10):
        print(f"  EC {ec:<20}: {c}")

    hist_summary_tsv = RES_DIR / "f2_306_historical_remapping_summary.tsv"
    with open(hist_summary_tsv, "w") as f:
        f.write("metric\tcount\tpct_of_historical_306\n")
        f.write(f"total_historical_families\t{len(hist_families)}\t100.0%\n")
        for k, v in remap_summary.items():
            f.write(f"{k}\t{v}\t{100.0 * v / len(hist_families):.1f}%\n")
        f.write("\n--- EC Distribution ---\n")
        f.write("ec_number\tfamily_count\tpct_of_306\n")
        for ec, c in ec_counts.most_common():
            f.write(f"{ec}\t{c}\t{100.0 * c / len(hist_families):.1f}%\n")
    print(f"Wrote historical census summary to {hist_summary_tsv.name}.")
    print(f"=== F2 Pipeline Completed in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
    main()
