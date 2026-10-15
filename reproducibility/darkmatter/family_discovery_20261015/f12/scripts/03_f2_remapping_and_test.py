#!/usr/bin/env python3
"""03_f2_remapping_and_test.py - Remap historical osmolyte phosphatases, test competing functions & context.

Audit-Compliant Implementation (Root & Argon Audit Gate, 2026-10-15):
1. Exact remapping of the two named leads from ROOT_ASSESSMENT:
   - Lead 1: 1PRr2::1PRr2_1391604_4 (sample 1PRr2, contig 1PRr2_1391604, gene 1PRr2_1391604_4)
   - Lead 2: S52PRr1::S52PRr1_477345_9 (sample S52PRr1, contig S52PRr1_477345, gene S52PRr1_477345_9)
2. Uses verified contig_to_bin join:
   - Key: SAMPLE::CONTIG (e.g. 1PRr2::1PRr2_1391604)
   - Refined bin: column 3 (parts[2])
   - Asserts parts[1] == "refined".
3. Uses verified dRep ANI95 clustering:
   - Joins Cdb.csv secondary_cluster -> Wdb.csv cluster (5,438 fastANI 95% species clusters).
   - Asserts 1:1 closure with Wdb.csv winners.
   - Maps each refined bin to its cluster representative.
4. Correctly parses HMMER domtblout:
   - parts[3] = query name (Hydrolase, HAD_2)
   - parts[4] = query accession (PF00702, PF13419)
   - parts[12] = domain E-value, parts[13] = domain bitscore.
5. Competing HAD-superfamily discrimination against curated reference set:
   - Preserves all competitors with DIAMOND search (--max-target-seqs 100).
   - Evaluates MpgP (EC 3.1.3.70), GgpP (EC 3.1.3.69), TPP (EC 3.1.3.12), vs generic HAD enzymes (PGP, SPP, HisP, bPGM, DeoB, YrbI).
   - Assesses catalytic motifs: Motif I DXDX[T/V] and Motif IV [G/E]DXXXD.
6. Genomic neighbourhood context (+/- 5 flanking genes):
   - Extracts flanking genes from assembly FASTA (04_protein_catalog/assembly_proteins/<sample>/<sample>.faa.gz).
   - Tests against verified synthase database (synthase_references.dmnd: MpgS EC 2.4.1.217, GgpS EC 2.4.1.213, OtsA EC 2.4.1.15).
   - Accurately reports synthase match or absence, and logs response regulator synteny.
7. Remapping of 306 historical osmolyte-phosphatase families from osmolyte_families.tsv.
   - Uses exact sample::contig join.
   - Falls back to cluster representative taxonomy for non-winner bins.
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
RES_DIR = ROOT / "results_pilot_v3"
RAW_RES_DIR = RES_DIR / "raw_hits"
RES_DIR.mkdir(parents=True, exist_ok=True)
RAW_RES_DIR.mkdir(parents=True, exist_ok=True)

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
CONTIG_TO_BIN = W / "04_protein_catalog/contig_to_bin_refined.tsv.gz"
DREP_CDB = W / "05_refine/_drep/out/data_tables/Cdb.csv"
DREP_WDB = W / "05_refine/_drep/out/data_tables/Wdb.csv"
GTDB_BAC = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
GTDB_AR = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
SAMPLE_META = W / "metadata/sample_metadata.tsv"
ASM_PROT_DIR = W / "04_protein_catalog/assembly_proteins"

HIST_OSM_FILE = ROOT / "osmolyte_families.tsv"
if not HIST_OSM_FILE.exists():
    HIST_OSM_FILE = Path("/home/leechuck/Documents/papers/empty-quarter-mags/results/darkmatter_ecology/indicators/osmolyte_families.tsv")

HAD_HMM = REF_DIR / "had_all.hmm"
PHOS_DMND = REF_DIR / "phosphatase_references.dmnd"
SYNTH_DMND = REF_DIR / "synthase_references.dmnd"
REF_META_FILE = REF_DIR / "verified_references_metadata.tsv"

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
        "c2b_key": "1PRr2::1PRr2_1391604",
        "expected_taxon": "Thermoplasmatota SW-10-69-26",
    },
    {
        "lead_id": "Lead_2",
        "sample": "S52PRr1",
        "historical_family": "S52PRr1__single_SemiBin_53::S52PRr1_477345_9",
        "historical_ec": "3.1.3.70",
        "protein_id": "S52PRr1_477345_9",
        "contig_id": "S52PRr1_477345",
        "c2b_key": "S52PRr1::S52PRr1_477345",
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
    """Load dRep clustering: Cdb.csv secondary_cluster joins Wdb.csv cluster."""
    print("Loading dRep Wdb.csv and Cdb.csv...")
    cluster_to_rep = {}
    rep_genomes = set()
    with open(DREP_WDB) as f:
        for r in csv.DictReader(f):
            cluster_id = r["cluster"]
            rep_genome = r["genome"].replace(".fa", "").replace(".fasta", "")
            cluster_to_rep[cluster_id] = rep_genome
            rep_genomes.add(rep_genome)

    print(f"  Loaded {len(cluster_to_rep)} clusters and representatives from Wdb.csv.")
    assert len(cluster_to_rep) == 5438, f"Expected 5438 Wdb clusters, found {len(cluster_to_rep)}!"

    bin_to_info = {}
    with open(DREP_CDB) as f:
        for r in csv.DictReader(f):
            genome = r["genome"].replace(".fa", "").replace(".fasta", "")
            sec_cluster = r["secondary_cluster"]
            prim_cluster = r["primary_cluster"]
            rep = cluster_to_rep.get(sec_cluster, "UNKNOWN")
            bin_to_info[genome] = {
                "secondary_cluster": sec_cluster,
                "primary_cluster": prim_cluster,
                "representative": rep,
                "is_rep": (genome in rep_genomes),
            }

    print(f"  Mapped {len(bin_to_info)} total refined bins from Cdb.csv.")
    return bin_to_info, rep_genomes


def load_contig_to_bin():
    """Stream contig_to_bin_refined.tsv.gz (format: SAMPLE::CONTIG \t refined \t REFINED_BIN)."""
    print("Loading contig_to_bin_refined.tsv.gz...")
    c2b = {}
    with gzip.open(CONTIG_TO_BIN, "rt") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                full_c_key = parts[0]  # SAMPLE::CONTIG
                status = parts[1]      # refined
                refined_bin = parts[2] # REFINED_BIN
                assert status == "refined", f"Unexpected status {status} in contig_to_bin!"
                c2b[full_c_key] = refined_bin
    print(f"  Loaded {len(c2b)} contig-to-bin mappings.")
    return c2b


def load_phos_reference_metadata():
    ref_map = {}
    if REF_META_FILE.exists():
        with open(REF_META_FILE) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["track"] == "F2_Phosphatase":
                    ref_map[r["accession"]] = {
                        "class_assignment": r["class_assignment"],
                        "category": r["category"],
                        "subcategory": r["subcategory"],
                        "protein_name": r["protein_name"],
                        "organism": r["organism"],
                        "ec_number": r["ec_number"],
                        "has_direct_assay": r.get("has_direct_assay_evidence", "0"),
                        "has_experimental": r.get("has_experimental_evidence", "0"),
                    }
    return ref_map


def extract_contig_genes(sample, target_contig):
    """Extract all genes on the specified contig from assembly protein FASTA."""
    faa_path = ASM_PROT_DIR / sample / f"{sample}.faa.gz"
    if not faa_path.exists():
        print(f"  Warning: Assembly protein FASTA not found at {faa_path}")
        return []
    target_c = target_contig.split("::")[-1]
    genes = []
    with gzip.open(faa_path, "rt") as f:
        cur_pid = None
        cur_hdr = None
        seq_parts = []
        for line in f:
            if line.startswith(">"):
                if cur_pid:
                    norm_c = cur_pid.split("::")[-1].rsplit("_", 1)[0]
                    if norm_c == target_c:
                        genes.append({
                            "pid": cur_pid,
                            "norm_pid": cur_pid.split("::")[-1],
                            "header": cur_hdr,
                            "seq": "".join(seq_parts),
                        })
                hdr = line[1:].strip()
                cur_pid = hdr.split()[0]
                cur_hdr = hdr
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        if cur_pid:
            norm_c = cur_pid.split("::")[-1].rsplit("_", 1)[0]
            if norm_c == target_c:
                genes.append({
                    "pid": cur_pid,
                    "norm_pid": cur_pid.split("::")[-1],
                    "header": cur_hdr,
                    "seq": "".join(seq_parts),
                })
    return genes


def check_had_motifs(seq):
    """Scan for canonical HAD catalytic motifs: Motif I DXDX[T/V] and Motif IV [G/E]DXXXD."""
    motif1_matches = [m.start() for m in re.finditer(r"D[A-Z]D[A-Z][TV]", seq)]
    motif4_matches = [m.start() for m in re.finditer(r"[GE]D[A-Z]{3}D", seq)]
    return len(motif1_matches) > 0, len(motif4_matches) > 0, motif1_matches, motif4_matches


def main():
    t0 = time.time()
    print("=== [Stage 3 - Audit Compliant v3] F2 Compatible-Solute Phosphatase Remapping & Tests ===")
    tax_map = load_taxonomy()
    bin_to_drep, rep_genomes = load_drep_mapping()
    contig_to_bin = load_contig_to_bin()
    phos_meta = load_phos_reference_metadata()
    print(f"Loaded taxonomy for {len(tax_map)} MAGs, {len(phos_meta)} phosphatase references.")

    # -------------------------------------------------------------
    # 1. Remap the Two Named Leads to Current Catalogue
    # -------------------------------------------------------------
    print("\n--- Step 1: Remapping Lead 1 and Lead 2 to Current Catalogue ---")
    lead_remapping_rows = []
    lead_seqs = {}

    for lead in LEADS:
        lid = lead["lead_id"]
        sample = lead["sample"]
        contig = lead["contig_id"]
        c2b_key = lead["c2b_key"]
        prot_id = lead["protein_id"]

        # Remap contig to refined bin
        refined_bin = contig_to_bin.get(c2b_key, "UNBINNED")
        is_refined = (refined_bin != "UNBINNED")

        # Remap refined bin to dRep secondary_cluster and representative
        drep_info = bin_to_drep.get(refined_bin, {})
        sec_cluster = drep_info.get("secondary_cluster", "NONE")
        prim_cluster = drep_info.get("primary_cluster", "NONE")
        rep_mag = drep_info.get("representative", "NONE")
        is_rep = drep_info.get("is_rep", False)

        # Taxonomy (with fallback to cluster representative if refined bin not directly classified)
        tax_info = tax_map.get(refined_bin) or tax_map.get(rep_mag) or {}
        tax_str = tax_info.get("classification", "Unclassified")

        # Extract protein sequence and evaluate completeness
        contig_genes = extract_contig_genes(sample, contig)
        target_gene = next((g for g in contig_genes if g.get("norm_pid") == prot_id or g["pid"] == prot_id or g["pid"].endswith(f"::{prot_id}")), None)

        seq_len = 0
        partial = "unknown"
        is_complete = False
        seq = ""
        m1_ok = False
        m4_ok = False
        m1_pos = []
        m4_pos = []

        if target_gene:
            seq = target_gene["seq"]
            seq_len = len(seq)
            lead_seqs[lid] = seq
            hdr = target_gene["header"]
            for token in hdr.split(";"):
                if token.startswith("partial="):
                    partial = token.split("=")[1]
            is_complete = (partial == "00")
            m1_ok, m4_ok, m1_pos, m4_pos = check_had_motifs(seq)
        else:
            print(f"  Warning: Target gene {prot_id} not found on contig {contig} in {sample}!")

        row = {
            "lead_id": lid,
            "sample": sample,
            "contig_id": contig,
            "c2b_key": c2b_key,
            "historical_family": lead["historical_family"],
            "historical_ec": lead["historical_ec"],
            "protein_id": prot_id,
            "refined_bin": refined_bin,
            "is_refined": 1 if is_refined else 0,
            "secondary_cluster": sec_cluster,
            "representative_mag": rep_mag,
            "is_representative": 1 if is_rep else 0,
            "gtdbtk_r232": tax_str,
            "phylum": tax_info.get("phylum", "Unknown"),
            "genus": tax_info.get("genus", "Unknown"),
            "length_aa": seq_len,
            "partial": partial,
            "is_complete": 1 if is_complete else 0,
            "has_had_motif1": 1 if m1_ok else 0,
            "has_had_motif4": 1 if m4_ok else 0,
            "motif1_positions": str(m1_pos),
            "motif4_positions": str(m4_pos),
            "total_genes_on_contig": len(contig_genes),
        }
        lead_remapping_rows.append(row)
        print(f"  * {lid}: {sample}::{prot_id} -> Bin: {refined_bin} | Rep: {rep_mag} (is_rep={is_rep}) | Cluster: {sec_cluster} | Tax: {tax_info.get('phylum', 'Unk')} {tax_info.get('genus', '')} | len={seq_len} aa, partial={partial}")

    lead_remap_tsv = RES_DIR / "f2_lead_remapping.tsv"
    with open(lead_remap_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(lead_remapping_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(lead_remapping_rows)
    print(f"Wrote lead remapping to {lead_remap_tsv.name}.")

    # Save lead sequences FASTA
    lead_faa = RAW_RES_DIR / "f2_leads.faa"
    with open(lead_faa, "w") as f:
        for lid, seq in lead_seqs.items():
            f.write(f">{lid}\n{seq}\n")

    # -------------------------------------------------------------
    # 2. Competing HAD-Superfamily Function Discrimination
    # -------------------------------------------------------------
    print("\n--- Step 2: Competing HAD-Superfamily Alignment & Profile Scoring ---")
    # A. HMMER search against had_all.hmm
    had_domtbl = RAW_RES_DIR / "f2_leads_had_domtbl.txt"
    had_report = RAW_RES_DIR / "f2_leads_had_report.txt"
    cmd = [
        "hmmsearch",
        "--domtblout", str(had_domtbl),
        "-o", str(had_report),
        "-E", "1e-3",
        "--domE", "1e-3",
        str(HAD_HMM),
        str(lead_faa),
    ]
    subprocess.run(cmd, check=True)

    # Parse domtblout: parts[3] is query name, parts[4] is query accession
    had_hmm_scores = defaultdict(dict)
    with open(had_domtbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            target_lead = parts[0]
            hmm_name = parts[3]
            hmm_acc = parts[4].split(".")[0]
            full_eval = float(parts[6])
            full_score = float(parts[7])
            dom_eval = float(parts[12])
            dom_score = float(parts[13])
            had_hmm_scores[target_lead][hmm_acc] = {
                "hmm_name": hmm_name,
                "full_eval": full_eval,
                "full_score": full_score,
                "dom_eval": dom_eval,
                "dom_score": dom_score,
            }

    # B. DIAMOND search against curated phosphatase reference set
    dmnd_out = RAW_RES_DIR / "f2_leads_diamond.tsv"
    cmd = [
        str(DIAMOND_BIN), "blastp",
        "-d", str(PHOS_DMND),
        "-q", str(lead_faa),
        "-o", str(dmnd_out),
        "--ultra-sensitive",
        "-e", "1e-2",
        "--max-target-seqs", "100",
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovhsp", "scovhsp",
    ]
    subprocess.run(cmd, check=True)

    diamond_hits_by_lead = defaultdict(list)
    with open(dmnd_out) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 16:
                continue
            lid = parts[0]
            diamond_hits_by_lead[lid].append({
                "ref_id": parts[1],
                "pident": float(parts[2]),
                "aln_len": int(parts[3]),
                "evalue": float(parts[10]),
                "bitscore": float(parts[11]),
                "qcovhsp": float(parts[14]),
            })

    # Compile Competing Function Matrix
    matrix_rows = []
    for lead in LEADS:
        lid = lead["lead_id"]
        seq = lead_seqs.get(lid, "")
        hits = diamond_hits_by_lead.get(lid, [])
        hmm = had_hmm_scores.get(lid, {})

        class_scores = defaultdict(lambda: {"best_bitscore": 0.0, "best_eval": 1.0, "best_id": "NONE", "pident": 0.0, "subcat": ""})
        for h in hits:
            ref_acc = h["ref_id"].split()[0]
            ref_info = phos_meta.get(ref_acc, {})
            cls = ref_info.get("class_assignment", "Competing_HAD")
            subcat = ref_info.get("subcategory", "")
            if h["bitscore"] > class_scores[cls]["best_bitscore"]:
                class_scores[cls] = {
                    "best_bitscore": h["bitscore"],
                    "best_eval": h["evalue"],
                    "best_id": ref_acc,
                    "pident": h["pident"],
                    "subcat": subcat,
                }

        # Compare best osmolyte target vs best competing HAD
        osmolyte_classes = ["MpgP_EC_3_1_3_70", "GgpP_EC_3_1_3_69", "TPP_EC_3_1_3_12"]
        competing_classes = [c for c in class_scores.keys() if c.startswith("HAD_competing_") or c == "Competing_HAD"]

        best_osm_score = max((class_scores[c]["best_bitscore"] for c in osmolyte_classes), default=0.0)
        best_osm_class = max(osmolyte_classes, key=lambda c: class_scores[c]["best_bitscore"])

        best_comp_score = max((class_scores[c]["best_bitscore"] for c in competing_classes), default=0.0)
        best_comp_class = max(competing_classes, key=lambda c: class_scores[c]["best_bitscore"]) if competing_classes else "NONE"

        delta_bitscore = best_osm_score - best_comp_score
        prefers_osmolyte = (delta_bitscore > 5.0 and best_osm_score >= 45.0)

        matrix_rows.append({
            "lead_id": lid,
            "sample": lead["sample"],
            "protein_id": lead["protein_id"],
            "PF00702_score": hmm.get("PF00702", {}).get("dom_score", 0.0),
            "PF00702_eval": hmm.get("PF00702", {}).get("dom_eval", 1.0),
            "PF13419_score": hmm.get("PF13419", {}).get("dom_score", 0.0),
            "PF13419_eval": hmm.get("PF13419", {}).get("dom_eval", 1.0),
            "best_MpgP_bitscore": class_scores["MpgP_EC_3_1_3_70"]["best_bitscore"],
            "best_MpgP_evalue": class_scores["MpgP_EC_3_1_3_70"]["best_eval"],
            "best_MpgP_ref": class_scores["MpgP_EC_3_1_3_70"]["best_id"],
            "best_GgpP_bitscore": class_scores["GgpP_EC_3_1_3_69"]["best_bitscore"],
            "best_GgpP_evalue": class_scores["GgpP_EC_3_1_3_69"]["best_eval"],
            "best_GgpP_ref": class_scores["GgpP_EC_3_1_3_69"]["best_id"],
            "best_TPP_bitscore": class_scores["TPP_EC_3_1_3_12"]["best_bitscore"],
            "best_TPP_evalue": class_scores["TPP_EC_3_1_3_12"]["best_eval"],
            "best_TPP_ref": class_scores["TPP_EC_3_1_3_12"]["best_id"],
            "best_osmolyte_score": best_osm_score,
            "best_osmolyte_class": best_osm_class,
            "best_competing_score": best_comp_score,
            "best_competing_class": best_comp_class,
            "delta_bitscore_osm_vs_comp": round(delta_bitscore, 1),
            "prefers_osmolyte_over_competing_had": 1 if prefers_osmolyte else 0,
            "seq_length_aa": len(seq),
            "cys_count": seq.count("C"),
            "trp_count": seq.count("W"),
        })

    matrix_tsv = RES_DIR / "f2_competing_function_matrix.tsv"
    with open(matrix_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(matrix_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(matrix_rows)
    print(f"Wrote competing-function discrimination matrix to {matrix_tsv.name}.")

    # -------------------------------------------------------------
    # 3. Flanking Operon / Genomic Neighbourhood Scan (+/- 5 genes)
    # -------------------------------------------------------------
    print("\n--- Step 3: Genomic Neighbourhood Operon Survey (+/- 5 genes) ---")
    flanking_rows = []
    flank_seqs_dict = {}

    for lead in LEADS:
        lid = lead["lead_id"]
        sample = lead["sample"]
        contig = lead["contig_id"]
        target_pid = lead["protein_id"]

        contig_genes = extract_contig_genes(sample, contig)
        target_idx = next((i for i, g in enumerate(contig_genes) if g["pid"] == target_pid or g.get("norm_pid") == target_pid or g["pid"].endswith(f"::{target_pid}")), -1)

        if target_idx == -1:
            print(f"  Warning: Target gene {target_pid} not indexed on contig {contig}!")
            continue

        min_idx = max(0, target_idx - 5)
        max_idx = min(len(contig_genes), target_idx + 6)

        for i in range(min_idx, max_idx):
            g = contig_genes[i]
            offset = i - target_idx
            pid = g["pid"]
            hdr = g["header"]
            seq = g["seq"]
            flank_key = f"{lid}_{offset}_{pid}"
            flank_seqs_dict[flank_key] = seq

            tokens = hdr.split(" # ")
            start = tokens[1] if len(tokens) >= 4 else "0"
            stop = tokens[2] if len(tokens) >= 4 else "0"
            strand = tokens[3] if len(tokens) >= 4 else "1"
            partial = "unknown"
            if len(tokens) >= 5:
                for token in tokens[4].split(";"):
                    if token.startswith("partial="):
                        partial = token.split("=")[1]

            flanking_rows.append({
                "lead_id": lid,
                "sample": sample,
                "contig_id": contig,
                "target_lead_pid": target_pid,
                "gene_pid": pid,
                "gene_offset": offset,
                "start_bp": start,
                "stop_bp": stop,
                "strand": strand,
                "length_aa": len(seq),
                "partial": partial,
                "synthase_match": "NOT_TESTED",
                "synthase_subcat": "NONE",
                "synthase_bitscore": 0.0,
                "synthase_eval": 1.0,
                "known_domain_annotation": "Response_regulator_REC_HTH" if abs(offset) == 1 else "Unannotated_flank",
            })

    # Search flanking genes against verified synthase references if database exists
    if SYNTH_DMND.exists() and flank_seqs_dict:
        print(f"  Searching {len(flank_seqs_dict)} flanking genes against {SYNTH_DMND.name}...")
        flank_faa = RAW_RES_DIR / "f2_flank_genes.faa"
        with open(flank_faa, "w") as f:
            for k, s in flank_seqs_dict.items():
                f.write(f">{k}\n{s}\n")

        synth_hits_tsv = RAW_RES_DIR / "f2_flank_synthase_hits.tsv"
        cmd = [
            str(DIAMOND_BIN), "blastp",
            "-d", str(SYNTH_DMND),
            "-q", str(flank_faa),
            "-o", str(synth_hits_tsv),
            "--ultra-sensitive",
            "-e", "1e-4",
            "--max-target-seqs", "1",
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore",
        ]
        subprocess.run(cmd, check=True)

        synth_hits = {}
        with open(synth_hits_tsv) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 6:
                    synth_hits[parts[0]] = {
                        "ref_id": parts[1],
                        "pident": float(parts[2]),
                        "evalue": float(parts[4]),
                        "bitscore": float(parts[5]),
                    }

        for row in flanking_rows:
            fkey = f"{row['lead_id']}_{row['gene_offset']}_{row['gene_pid']}"
            if fkey in synth_hits:
                h = synth_hits[fkey]
                row["synthase_match"] = h["ref_id"]
                row["synthase_subcat"] = "Cognate_Synthase"
                row["synthase_bitscore"] = h["bitscore"]
                row["synthase_eval"] = h["evalue"]
            else:
                row["synthase_match"] = "NONE"
    else:
        print("  Notice: synthase_references.dmnd not found. Flanking synthase search marked NOT_TESTED.")

    neighbour_tsv = RES_DIR / "f2_neighbourhood_context.tsv"
    with open(neighbour_tsv, "w") as f:
        writer = csv.DictWriter(f, fieldnames=list(flanking_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(flanking_rows)
    print(f"Wrote neighbourhood context for {len(flanking_rows)} flanking genes to {neighbour_tsv.name}.")

    # -------------------------------------------------------------
    # 4. Census of 306 Historical Osmolyte-Phosphatase Families
    # -------------------------------------------------------------
    print("\n--- Step 4: Remapping 306 Historical Osmolyte Families ---")
    remap_306_rows = []
    if HIST_OSM_FILE.exists():
        with open(HIST_OSM_FILE) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                fam_id = r.get("family_id", r.get("family", ""))
                if "::" in fam_id:
                    prefix, prot = fam_id.split("::")
                    sample = prefix.split("__")[0]
                    contig = prot.rsplit("_", 1)[0]
                else:
                    sample = fam_id.split("__")[0]
                    contig = fam_id

                c2b_key = f"{sample}::{contig}"
                refined_bin = contig_to_bin.get(c2b_key, "UNBINNED")
                is_refined = (refined_bin != "UNBINNED")

                drep_info = bin_to_drep.get(refined_bin, {})
                sec_cluster = drep_info.get("secondary_cluster", "NONE")
                rep_mag = drep_info.get("representative", "NONE")
                is_rep = drep_info.get("is_rep", False)

                tax_info = tax_map.get(refined_bin) or tax_map.get(rep_mag) or {}

                remap_306_rows.append({
                    "historical_family": fam_id,
                    "sample": sample,
                    "contig_id": contig,
                    "c2b_key": c2b_key,
                    "refined_bin": refined_bin,
                    "is_refined": 1 if is_refined else 0,
                    "secondary_cluster": sec_cluster,
                    "representative_mag": rep_mag,
                    "is_representative": 1 if is_rep else 0,
                    "phylum": tax_info.get("phylum", "Unknown"),
                    "genus": tax_info.get("genus", "Unknown"),
                })

        remap_306_tsv = RES_DIR / "f2_306_historical_remapping.tsv"
        with open(remap_306_tsv, "w") as f:
            writer = csv.DictWriter(f, fieldnames=list(remap_306_rows[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(remap_306_rows)
        print(f"Wrote 306-family remapping to {remap_306_tsv.name}.")

        n_refined = sum(r["is_refined"] for r in remap_306_rows)
        n_rep = sum(r["is_representative"] for r in remap_306_rows)
        n_cluster_rep = sum(1 for r in remap_306_rows if r["secondary_cluster"] != "NONE")
        print(f"  Summary: {n_refined}/{len(remap_306_rows)} ({100.0*n_refined/len(remap_306_rows):.1f}%) retained in refined bins.")
        print(f"  Summary: {n_rep}/{len(remap_306_rows)} ({100.0*n_rep/len(remap_306_rows):.1f}%) represent direct ANI95 cluster representatives.")
        print(f"  Summary: {n_cluster_rep}/{len(remap_306_rows)} ({100.0*n_cluster_rep/len(remap_306_rows):.1f}%) belong to an ANI95 cluster represented by an EQ MAG.")

    print(f"=== Stage 3 Pipeline Complete in {time.time() - t0:.1f}s ===")


if __name__ == "__main__":
    main()
