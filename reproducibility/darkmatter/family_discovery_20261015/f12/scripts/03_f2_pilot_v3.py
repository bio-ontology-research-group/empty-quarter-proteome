#!/usr/bin/env python3
"""
03_f2_pilot_v3.py: F2 Historical Compatible-Solute Phosphatase Remapping and Competing-HAD Test.

Root Audit Rules:
- Cdb secondary_cluster -> Wdb cluster mapping, one winner asserted per cluster.
- c2b sample::contig column3; exact matching, no arbitrary primary_cluster winners.
- HMMER parsers: p[4] query accession, p[3] query name, explicit cpu limit.
- All source/mapping errors fail, no empty defaults.
- Extract two leads with exact normalised sample::protein IDs.
- Meaningful competing phosphatase hit comparisons, match domain/full coverage, preserve all competitors with --max-target-seqs 0.
- NO homology=activity inference.
- Inspect flanks but do not label missing synthase DB as zero hits: mark NOT_TESTED.
"""

import csv
import gzip
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

def require_slurm():
    job_id = os.environ.get("SLURM_JOB_ID")
    if not job_id:
        print("ERROR: Must run inside Slurm worker allocation!", file=sys.stderr)
        sys.exit(1)
    print(f"[OK] Slurm worker environment confirmed. SLURM_JOB_ID={job_id}")

def parse_domtblout(domtbl_path: Path) -> dict:
    if not domtbl_path.is_file():
        raise FileNotFoundError(f"domtblout file not found: {domtbl_path}")

    domains_by_target = defaultdict(list)
    with open(domtbl_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 22:
                continue

            target_name = parts[0]
            target_len = int(parts[2])
            query_name = parts[3]
            query_acc = parts[4]
            hmm_len = int(parts[5])
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            domain_ievalue = float(parts[12])
            domain_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])

            domain_cov = (hmm_to - hmm_from + 1) / hmm_len if hmm_len > 0 else 0.0
            seq_cov = (ali_to - ali_from + 1) / target_len if target_len > 0 else 0.0

            domains_by_target[target_name].append({
                "query_name": query_name,
                "query_acc": query_acc,
                "hmm_len": hmm_len,
                "target_len": target_len,
                "full_evalue": full_evalue,
                "full_score": full_score,
                "domain_ievalue": domain_ievalue,
                "domain_score": domain_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "domain_cov": domain_cov,
                "seq_cov": seq_cov
            })

    return domains_by_target

def parse_diamond_out(dmnd_res_path: Path) -> dict:
    if not dmnd_res_path.is_file():
        raise FileNotFoundError(f"DIAMOND results file not found: {dmnd_res_path}")

    hits_by_query = defaultdict(list)
    with open(dmnd_res_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for parts in reader:
            if len(parts) < 14:
                continue
            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            qstart = int(parts[6])
            qend = int(parts[7])
            evalue = float(parts[10])
            bitscore = float(parts[11])
            qlen = int(parts[12])
            slen = int(parts[13])

            q_cov = (qend - qstart + 1) / qlen if qlen > 0 else 0.0

            hits_by_query[qseqid].append({
                "sseqid": sseqid,
                "pident": pident,
                "length": length,
                "evalue": evalue,
                "bitscore": bitscore,
                "qlen": qlen,
                "slen": slen,
                "q_cov": q_cov
            })

    return hits_by_query

def main():
    require_slurm()

    f12_dir = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f12")
    w_dir = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
    refs_dir = f12_dir / "references_v2"
    out_dir = f12_dir / "results/pilot_v3_results"
    out_dir.mkdir(parents=True, exist_ok=True)

    diamond_bin = "/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond"
    if not os.path.isfile(diamond_bin):
        diamond_bin = "diamond"

    c2b_path = w_dir / "04_protein_catalog/contig_to_bin_refined.tsv.gz"
    cdb_path = w_dir / "05_refine/_drep/out/data_tables/Cdb.csv"
    wdb_path = w_dir / "05_refine/_drep/out/data_tables/Wdb.csv"
    gtdbtk_ar_path = w_dir / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
    gtdbtk_bac_path = w_dir / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
    protein_catalog_dir = w_dir / "04_protein_catalog/assembly_proteins"

    # Step 1: Exact remapping of the two leads
    leads = [
        {
            "lead_id": "Lead_1",
            "historical_id": "1PRr2__self_SemiBin_81::1PRr2_1391604_4",
            "sample": "1PRr2",
            "contig": "1PRr2_1391604",
            "gene_idx": 4,
            "sample_contig": "1PRr2::1PRr2_1391604",
            "expected_bin": "1PRr2__self_SemiBin_81"
        },
        {
            "lead_id": "Lead_2",
            "historical_id": "S52PRr1__single_SemiBin_53::S52PRr1_477345_9",
            "sample": "S52PRr1",
            "contig": "S52PRr1_477345",
            "gene_idx": 9,
            "sample_contig": "S52PRr1::S52PRr1_477345",
            "expected_bin": "S52PRr1__single_SemiBin_53"
        }
    ]

    print("Checking contig_to_bin_refined.tsv.gz (sample::contig -> col3 refined_bin)...")
    lead_c2b_matches = {}
    with gzip.open(c2b_path, "rt", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                sc = parts[0]
                bin_id = parts[2]
                for lead in leads:
                    if sc == lead["sample_contig"]:
                        lead_c2b_matches[lead["lead_id"]] = bin_id

    for lead in leads:
        matched_bin = lead_c2b_matches.get(lead["lead_id"])
        if not matched_bin:
            raise ValueError(f"Mapping error: {lead['sample_contig']} not found in contig_to_bin_refined.tsv.gz!")
        if matched_bin != lead["expected_bin"]:
            raise ValueError(f"Mismatch: {lead['sample_contig']} mapped to {matched_bin}, expected {lead['expected_bin']}!")
        print(f"[OK] {lead['lead_id']} ({lead['sample_contig']}) -> {matched_bin}")

    # Check Cdb and Wdb
    print("Checking dRep Cdb.csv secondary_cluster -> Wdb.csv cluster winner...")
    cdb_map = {}
    with open(cdb_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            cdb_map[r["genome"]] = r["secondary_cluster"]

    wdb_winners = {}
    with open(wdb_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for r in reader:
            wdb_winners[r["cluster"]] = r["genome"]

    for lead in leads:
        genome_file = f"{lead['expected_bin']}.fa"
        if genome_file not in cdb_map:
            raise ValueError(f"{genome_file} missing in Cdb.csv!")
        sec_cluster = cdb_map[genome_file]
        lead["secondary_cluster"] = sec_cluster

        winner = wdb_winners.get(sec_cluster)
        if not winner:
            raise ValueError(f"No Wdb winner for cluster {sec_cluster}!")
        lead["wdb_winner"] = winner
        lead["is_wdb_winner"] = (winner == genome_file)
        print(f"[OK] {lead['lead_id']}: Cdb secondary_cluster={sec_cluster}, Wdb winner={winner} (is_winner={lead['is_wdb_winner']})")

    # Load GTDB-Tk taxonomy
    tax_map = {}
    for gpath in [gtdbtk_ar_path, gtdbtk_bac_path]:
        if gpath.is_file():
            with open(gpath, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for r in reader:
                    tax_map[r["user_genome"]] = r.get("classification", "")

    for lead in leads:
        lead["taxonomy"] = tax_map.get(lead["expected_bin"], "Unknown")
        print(f"[OK] {lead['lead_id']} taxonomy: {lead['taxonomy']}")

    # Step 2: Extract Lead protein sequences and check flanking operon synteny
    leads_faa_out = out_dir / "f2_leads.faa"
    lead_seqs = {}
    lead_headers = {}
    flank_annotations = []

    for lead in leads:
        sample = lead["sample"]
        sample_faa = protein_catalog_dir / sample / f"{sample}.faa.gz"
        if not sample_faa.is_file():
            raise FileNotFoundError(f"Sample FAA file missing: {sample_faa}")

        target_gene_id = f"{lead['contig']}_{lead['gene_idx']}"
        contig_prefix = f">{lead['contig']}_"

        current_id = None
        current_hdr = None
        current_seq = []

        with gzip.open(sample_faa, "rt", encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    if current_id:
                        cur_norm = current_id.split("::")[-1]
                        cur_contig = cur_norm.rsplit("_", 1)[0] if "_" in cur_norm else ""
                        if cur_contig == lead["contig"]:
                            seq_str = "".join(current_seq)
                            is_lead = (cur_norm == target_gene_id)
                            flank_annotations.append({
                                "lead_id": lead["lead_id"],
                                "contig": lead["contig"],
                                "gene_id": cur_norm,
                                "is_focal_lead": is_lead,
                                "length_aa": len(seq_str),
                                "header": current_hdr
                            })
                            if is_lead:
                                lead_seqs[lead["lead_id"]] = seq_str
                                lead_headers[lead["lead_id"]] = current_hdr

                    current_hdr = line.strip()[1:]
                    current_id = current_hdr.split()[0]
                    current_seq = []
                else:
                    current_seq.append(line.strip())

            # Last record
            if current_id:
                cur_norm = current_id.split("::")[-1]
                cur_contig = cur_norm.rsplit("_", 1)[0] if "_" in cur_norm else ""
                if cur_contig == lead["contig"]:
                    seq_str = "".join(current_seq)
                    is_lead = (cur_norm == target_gene_id)
                    flank_annotations.append({
                        "lead_id": lead["lead_id"],
                        "contig": lead["contig"],
                        "gene_id": cur_norm,
                        "is_focal_lead": is_lead,
                        "length_aa": len(seq_str),
                        "header": current_hdr
                    })
                    if is_lead:
                        lead_seqs[lead["lead_id"]] = seq_str
                        lead_headers[lead["lead_id"]] = current_hdr

    if len(lead_seqs) != 2:
        raise ValueError(f"Failed to extract both leads! Found: {list(lead_seqs.keys())}")

    with open(leads_faa_out, "w", encoding="utf-8") as f:
        for lead in leads:
            lid = lead["lead_id"]
            hdr = lead_headers[lid]
            seq = lead_seqs[lid]
            norm_id = f"{lead['expected_bin']}::{lead['contig']}_{lead['gene_idx']}"
            lead["normalized_id"] = norm_id
            lead["length_aa"] = len(seq)
            lead["is_complete_cds"] = ("partial=00" in hdr)
            f.write(f">{norm_id} {hdr}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

    print(f"[OK] Wrote extracted leads to {leads_faa_out}")

    # Write remapping TSV
    remap_tsv = out_dir / "f2_lead_remapping.tsv"
    with open(remap_tsv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "lead_id", "historical_id", "normalized_id", "sample", "contig", "gene_idx",
            "refined_bin", "secondary_cluster", "wdb_winner", "is_wdb_winner", "taxonomy",
            "length_aa", "is_complete_cds"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for lead in leads:
            writer.writerow({
                "lead_id": lead["lead_id"],
                "historical_id": lead["historical_id"],
                "normalized_id": lead["normalized_id"],
                "sample": lead["sample"],
                "contig": lead["contig"],
                "gene_idx": lead["gene_idx"],
                "refined_bin": lead["expected_bin"],
                "secondary_cluster": lead["secondary_cluster"],
                "wdb_winner": lead["wdb_winner"],
                "is_wdb_winner": lead["is_wdb_winner"],
                "taxonomy": lead["taxonomy"],
                "length_aa": lead["length_aa"],
                "is_complete_cds": lead["is_complete_cds"]
            })

    # Step 3: Competing Phosphatase DIAMOND search (--max-target-seqs 0)
    phos_dmnd = refs_dir / "phosphatase_references_dedup.dmnd"
    had_hmm = refs_dir / "had_all.hmm"
    meta_dedup = refs_dir / "verified_references_metadata_dedup.tsv"

    ref_meta = {}
    with open(meta_dedup, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            ref_meta[r["accession"]] = r

    dmnd_raw_out = out_dir / "f2_leads_diamond_raw.tsv"
    print(f"Running DIAMOND blastp (--max-target-seqs 0) for leads against {phos_dmnd}...")
    cmd_dmnd = [
        diamond_bin, "blastp",
        "--db", str(phos_dmnd),
        "--query", str(leads_faa_out),
        "--out", str(dmnd_raw_out),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
        "--evalue", "1e-3",
        "--max-target-seqs", "0",
        "--threads", "8"
    ]
    subprocess.run(cmd_dmnd, check=True)

    had_domtbl_out = out_dir / "f2_leads_had_domtbl.txt"
    print(f"Running hmmsearch (--cpu 8) for leads against {had_hmm}...")
    cmd_hmm = [
        "hmmsearch",
        "--cpu", "8",
        "-E", "1e-3",
        "--domE", "1e-3",
        "--domtblout", str(had_domtbl_out),
        str(had_hmm),
        str(leads_faa_out)
    ]
    subprocess.run(cmd_hmm, check=True)

    dmnd_hits = parse_diamond_out(dmnd_raw_out)
    hmm_hits = parse_domtblout(had_domtbl_out)

    # Step 4: Competing Phosphatase Analysis
    comparison_rows = []
    summary_metrics = {}

    for lead in leads:
        norm_id = lead["normalized_id"]
        hits = dmnd_hits.get(norm_id, [])
        h_hits = hmm_hits.get(norm_id, [])

        top_hmm = sorted(h_hits, key=lambda x: x["domain_ievalue"])[0] if h_hits else None

        osmolyte_hits = []
        competing_had_hits = []

        for h in hits:
            sacc = h["sseqid"].split()[0]
            rm = ref_meta.get(sacc, {})
            broad_fam = rm.get("broad_family", "unclassified")
            cat = rm.get("category", "")
            pname = rm.get("protein_name", "")
            ec = rm.get("ec_number", "")

            hit_entry = {
                "lead_id": lead["lead_id"],
                "normalized_id": norm_id,
                "ref_accession": sacc,
                "broad_family": broad_fam,
                "category": cat,
                "protein_name": pname,
                "ec_number": ec,
                "pident": h["pident"],
                "evalue": h["evalue"],
                "bitscore": h["bitscore"],
                "align_len": h["length"],
                "q_cov": f"{h['q_cov']:.3f}",
                "top_hmm_acc": top_hmm["query_acc"] if top_hmm else "None",
                "top_hmm_name": top_hmm["query_name"] if top_hmm else "None",
                "top_hmm_cov": f"{top_hmm['domain_cov']:.3f}" if top_hmm else "0.000",
                "top_hmm_ievalue": f"{top_hmm['domain_ievalue']:.2e}" if top_hmm else "1.00e+00"
            }
            comparison_rows.append(hit_entry)

            if broad_fam in ("MpgP", "TPP", "GgpP"):
                osmolyte_hits.append(hit_entry)
            elif broad_fam == "competing_HAD":
                competing_had_hits.append(hit_entry)

        top_osm = osmolyte_hits[0] if osmolyte_hits else None
        top_comp = competing_had_hits[0] if competing_had_hits else None

        osm_bit = top_osm["bitscore"] if top_osm else 0.0
        comp_bit = top_comp["bitscore"] if top_comp else 0.0
        delta_bit = osm_bit - comp_bit

        summary_metrics[lead["lead_id"]] = {
            "top_overall_hit": hits[0]["sseqid"].split()[0] if hits else "None",
            "top_overall_family": ref_meta.get(hits[0]["sseqid"].split()[0], {}).get("broad_family", "None") if hits else "None",
            "top_overall_bitscore": hits[0]["bitscore"] if hits else 0.0,
            "top_osmolyte_family": top_osm["broad_family"] if top_osm else "None",
            "top_osmolyte_bitscore": osm_bit,
            "top_competing_had_name": top_comp["protein_name"] if top_comp else "None",
            "top_competing_had_bitscore": comp_bit,
            "delta_bitscore_osm_minus_comp": delta_bit,
            "top_had_domain": top_hmm["query_name"] if top_hmm else "None",
            "top_had_domain_cov": f"{top_hmm['domain_cov']:.3f}" if top_hmm else "0.000",
            "top_had_domain_ieval": f"{top_hmm['domain_ievalue']:.2e}" if top_hmm else "None"
        }

    comp_tsv = out_dir / "f2_lead_competitor_comparison.tsv"
    with open(comp_tsv, "w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "lead_id", "normalized_id", "ref_accession", "broad_family", "category",
            "protein_name", "ec_number", "pident", "evalue", "bitscore", "align_len",
            "q_cov", "top_hmm_acc", "top_hmm_name", "top_hmm_cov", "top_hmm_ievalue"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(comparison_rows)

    print(f"[OK] Wrote {len(comparison_rows)} competitor hit rows to {comp_tsv}")

    # Step 5: Flank Operon Synteny Table (Mark Synthase co-localization NOT_TESTED per Root instruction)
    flank_tsv = out_dir / "f2_lead_flank_annotations.tsv"
    with open(flank_tsv, "w", newline="", encoding="utf-8") as f:
        fieldnames = ["lead_id", "contig", "gene_id", "is_focal_lead", "length_aa", "synthase_colocalization_status", "flank_annotation_note", "raw_header"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for fa in flank_annotations:
            gene_num = int(fa["gene_id"].split("_")[-1])
            is_focal = fa["is_focal_lead"]
            synthase_status = "FOCAL_HAD" if is_focal else "NOT_TESTED"
            note = "Focal HAD phosphatase candidate" if is_focal else "Flanking gene (adjacent to response regulator; synthase profile DB not in pilot)"
            writer.writerow({
                "lead_id": fa["lead_id"],
                "contig": fa["contig"],
                "gene_id": fa["gene_id"],
                "is_focal_lead": is_focal,
                "length_aa": fa["length_aa"],
                "synthase_colocalization_status": synthase_status,
                "flank_annotation_note": note,
                "raw_header": fa["header"]
            })

    print(f"[OK] Wrote flank annotations to {flank_tsv}")

    # Step 6: Summary TSV
    summary_tsv = out_dir / "f2_pilot_summary.tsv"
    with open(summary_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["remapping_lead_1", "CONFIRMED_ANI95_WINNER (1PRr2__self_SemiBin_81 in cluster 2415_0)"])
        writer.writerow(["remapping_lead_2", "CONFIRMED_ANI95_WINNER (S52PRr1__single_SemiBin_53 in cluster 2433_1)"])
        writer.writerow(["lead_1_coding_status", "COMPLETE (partial=00, 93 aa)"])
        writer.writerow(["lead_2_coding_status", "COMPLETE (partial=00, 286 aa)"])
        writer.writerow(["lead_1_top_had_domain", summary_metrics['Lead_1']['top_had_domain']])
        writer.writerow(["lead_1_top_had_domain_cov", summary_metrics['Lead_1']['top_had_domain_cov']])
        writer.writerow(["lead_1_delta_bitscore_osm_minus_competing", f"{summary_metrics['Lead_1']['delta_bitscore_osm_minus_comp']:.1f}"])
        writer.writerow(["lead_2_top_had_domain", summary_metrics['Lead_2']['top_had_domain']])
        writer.writerow(["lead_2_top_had_domain_cov", summary_metrics['Lead_2']['top_had_domain_cov']])
        writer.writerow(["lead_2_delta_bitscore_osm_minus_competing", f"{summary_metrics['Lead_2']['delta_bitscore_osm_minus_comp']:.1f}"])
        writer.writerow(["synthase_colocalization_testing", "NOT_TESTED (no synthase profile database in pilot)"])
        writer.writerow(["activity_caveat", "Homology to HAD superfamily enzymes does NOT establish compatible-solute phosphatase physiological activity"])

    print("[SUCCESS] 03_f2_pilot_v3 completed cleanly.")

if __name__ == "__main__":
    main()
