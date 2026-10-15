#!/usr/bin/env python3
"""
scripts/recompute_pilot_summaries.py
Recomputes f{3,4,5}_candidates_summary.tsv from raw pilot DIAMOND hits and candidate data.
Fixes:
1. Identifier semantics: FASTA header was >{pid}|{mag}|{class}|{cluster}, so DIAMOND qseqid
   retained the whole pipe ID, while summary joined by bare protein_id.
   Normalizes qseqid to bare pid and asserts every qseqid matches a candidate row.
2. Renames reciprocal blast labels to one-way reference blastp (ref_blastp_*).
3. F4 synteny: Explicitly separates nearby roles EctA, EctB, EctD with distance and strand.
   Differentiates complete operon (ABC, ABCD) vs partial clusters (BCD, BC, CD) vs unlinked
   pathway carriers, and incorporates CheckM MAG completeness to respect incompleteness.
4. Preserves original broken versions as *_v0_brokenjoin.tsv.
"""

import csv
import os
import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
PILOT_DIR = BASE_DIR / "results" / "pilot_v1"
META_DIR = BASE_DIR / "metadata"

def load_checkm_metadata():
    info_path = META_DIR / "genomeInformation.csv"
    res = {}
    if info_path.exists():
        with open(info_path, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                mag = r["genome"].replace(".fa", "")
                res[mag] = {
                    "completeness": float(r["completeness"]) if r["completeness"] else None,
                    "contamination": float(r["contamination"]) if r["contamination"] else None,
                }
    return res

def parse_raw_diamond(dmnd_path, valid_pids):
    """
    Parses raw DIAMOND output (6 columns: qseqid, sseqid, pident, length, evalue, bitscore).
    Normalizes qseqid to bare protein_id and asserts membership in valid_pids.
    Returns dict: pid -> hit dict.
    """
    hits = {}
    if not dmnd_path.exists():
        return hits
    with open(dmnd_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            assert len(parts) >= 6, f"Malformed line in {dmnd_path}: {line}"
            qseqid = parts[0]
            pid = qseqid.split("|")[0]
            assert pid in valid_pids, f"Assertion failed: {pid} (from {qseqid}) not in candidate set!"
            if pid not in hits:
                hits[pid] = {
                    "ref_blastp_hit_id": parts[1],
                    "ref_blastp_pident": float(parts[2]),
                    "ref_blastp_align_len": int(parts[3]),
                    "ref_blastp_evalue": float(parts[4]),
                    "ref_blastp_bitscore": float(parts[5]),
                }
    return hits

def recompute_f3(checkm_info):
    print(">>> Recomputing F3 Rhodopsin Summary...")
    track_dir = PILOT_DIR / "f3_rhodopsin"
    v0_file = track_dir / "f3_candidates_summary_v0_brokenjoin.tsv"
    out_file = track_dir / "f3_candidates_summary.tsv"
    dmnd_file = track_dir / "f3_diamond_reciprocal.tsv"
    
    with open(v0_file, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    valid_pids = {r["protein_id"] for r in rows}
    
    dmnd_hits = parse_raw_diamond(dmnd_file, valid_pids)
    print(f"  Loaded {len(rows)} candidates, {len(dmnd_hits)} DIAMOND hits.")
    
    out_rows = []
    for r in rows:
        pid = r["protein_id"]
        hit = dmnd_hits.get(pid, {
            "ref_blastp_hit_id": "NO_HIT",
            "ref_blastp_pident": 0.0,
            "ref_blastp_align_len": 0,
            "ref_blastp_evalue": 1.0,
            "ref_blastp_bitscore": 0.0,
        })
        new_row = dict(r)
        # Update both new ref_blastp_* and legacy top_ref_* columns
        new_row["ref_blastp_hit_id"] = hit["ref_blastp_hit_id"]
        new_row["ref_blastp_pident"] = hit["ref_blastp_pident"]
        new_row["ref_blastp_align_len"] = hit["ref_blastp_align_len"]
        new_row["ref_blastp_evalue"] = hit["ref_blastp_evalue"]
        new_row["ref_blastp_bitscore"] = hit["ref_blastp_bitscore"]
        new_row["top_ref_id"] = hit["ref_blastp_hit_id"]
        new_row["top_ref_pident"] = hit["ref_blastp_pident"]
        new_row["top_ref_evalue"] = hit["ref_blastp_evalue"]
        new_row["top_ref_bitscore"] = hit["ref_blastp_bitscore"]
        out_rows.append(new_row)
        
    fieldnames = list(out_rows[0].keys())
    with open(out_file, "w", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"  Wrote {len(out_rows)} rows to {out_file.name}")

def recompute_f4(checkm_info):
    print(">>> Recomputing F4 Ectoine Synthase Summary with Granular Synteny...")
    track_dir = PILOT_DIR / "f4_ectc"
    v0_file = track_dir / "f4_candidates_summary_v0_brokenjoin.tsv"
    out_file = track_dir / "f4_candidates_summary.tsv"
    dmnd_file = track_dir / "f4_diamond_reciprocal.tsv"
    
    with open(v0_file, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    valid_pids = {r["protein_id"] for r in rows}
    
    dmnd_hits = parse_raw_diamond(dmnd_file, valid_pids)
    print(f"  Loaded {len(rows)} candidates, {len(dmnd_hits)} DIAMOND hits.")
    
    out_rows = []
    for r in rows:
        pid = r["protein_id"]
        mag = r["mag"]
        hit = dmnd_hits.get(pid, {
            "ref_blastp_hit_id": "NO_HIT",
            "ref_blastp_pident": 0.0,
            "ref_blastp_align_len": 0,
            "ref_blastp_evalue": 1.0,
            "ref_blastp_bitscore": 0.0,
        })
        
        # Granular synteny analysis
        neighbors_str = r.get("neighbor_gene_ids", "")
        # Parse nearby roles: EctA, EctB, EctD
        nearby_a = "NONE"
        nearby_b = "NONE"
        nearby_d = "NONE"
        has_nearby_a = False
        has_nearby_b = False
        has_nearby_d = False
        
        # Check genome-wide pathway roles from unlinked or neighbor info
        # If neighbors_str contains "unlinked_contig", it's an unlinked carrier
        is_unlinked = "unlinked_contig" in neighbors_str
        
        has_mag_a = "EctA:" in neighbors_str
        has_mag_b = "EctB:" in neighbors_str
        has_mag_d = "EctD:" in neighbors_str
        
        if not is_unlinked and neighbors_str != "none_detected" and neighbors_str:
            items = neighbors_str.split(";")
            for item in items:
                if item.startswith("EctA:"):
                    nearby_a = item.replace("EctA:", "")
                    has_nearby_a = True
                elif item.startswith("EctB:"):
                    nearby_b = item.replace("EctB:", "")
                    has_nearby_b = True
                elif item.startswith("EctD:"):
                    nearby_d = item.replace("EctD:", "")
                    has_nearby_d = True
                    
        # Classify exact operon architecture
        if has_nearby_a and has_nearby_b and has_nearby_d:
            operon_arch = "complete_operon_ABCD"
            synteny_detail = "Complete canonical operon (EctA+B+C+D co-localized on same contig)"
        elif has_nearby_a and has_nearby_b and not has_nearby_d:
            operon_arch = "complete_operon_ABC"
            synteny_detail = "Complete canonical operon (EctA+B+C co-localized on same contig; EctD absent/optional)"
        elif has_nearby_b and has_nearby_d and not has_nearby_a:
            operon_arch = "partial_cluster_BCD"
            synteny_detail = "Partial cluster (EctB+C+D co-localized; EctA absent on contig)"
        elif has_nearby_b and not has_nearby_a and not has_nearby_d:
            operon_arch = "partial_cluster_BC"
            synteny_detail = "Partial cluster (EctB+C co-localized; EctA and EctD absent on contig)"
        elif has_nearby_d and not has_nearby_a and not has_nearby_b:
            operon_arch = "partial_cluster_CD"
            synteny_detail = "Partial cluster (EctC+D co-localized; EctA and EctB absent on contig)"
        elif has_nearby_a and not has_nearby_b:
            operon_arch = "partial_cluster_AC"
            synteny_detail = "Partial cluster (EctA+C co-localized; EctB absent on contig)"
        elif is_unlinked:
            if has_mag_a and has_mag_b and has_mag_d:
                operon_arch = "unlinked_carrier_mag_has_ABCD"
                synteny_detail = "Solitary EctC on contig; complete pathway (A+B+D) present elsewhere in MAG"
            elif has_mag_a and has_mag_b:
                operon_arch = "unlinked_carrier_mag_has_ABC"
                synteny_detail = "Solitary EctC on contig; complete pathway (A+B) present elsewhere in MAG"
            elif has_mag_b:
                operon_arch = "unlinked_carrier_mag_missing_ectA"
                synteny_detail = "Solitary EctC on contig; EctB present elsewhere in MAG, but EctA missing in assembly"
            else:
                operon_arch = "unlinked_carrier_incomplete"
                synteny_detail = "Solitary EctC on contig; incomplete markers elsewhere in MAG"
        else:
            operon_arch = "solitary_ectC_no_pathway_in_assembly"
            synteny_detail = "Solitary EctC candidate; no pathway markers detected in MAG assembly"
            
        # MAG CheckM quality
        cinfo = checkm_info.get(mag, {})
        comp = cinfo.get("completeness")
        cont = cinfo.get("contamination")
        comp_str = f"{comp:.1f}%" if comp is not None else "Unknown"
        cont_str = f"{cont:.1f}%" if cont is not None else "Unknown"
        
        # MAG pathway completeness
        if has_mag_a and has_mag_b and has_mag_d:
            mag_pw_comp = "complete_pathway_ABCD"
        elif has_mag_a and has_mag_b:
            mag_pw_comp = "complete_pathway_ABC"
        elif has_mag_b:
            mag_pw_comp = "incomplete_pathway_missing_EctA"
        elif has_mag_a:
            mag_pw_comp = "incomplete_pathway_missing_EctB"
        else:
            mag_pw_comp = "isolated_EctC_only"
            
        new_row = dict(r)
        new_row["operon_architecture"] = operon_arch
        new_row["synteny_detail"] = synteny_detail
        new_row["nearby_ectA"] = nearby_a
        new_row["nearby_ectB"] = nearby_b
        new_row["nearby_ectD"] = nearby_d
        new_row["mag_has_ectA"] = "1" if has_mag_a else "0"
        new_row["mag_has_ectB"] = "1" if has_mag_b else "0"
        new_row["mag_has_ectD"] = "1" if has_mag_d else "0"
        new_row["mag_pathway_completeness"] = mag_pw_comp
        new_row["mag_checkm_completeness"] = comp_str
        new_row["mag_checkm_contamination"] = cont_str
        
        # Update both new ref_blastp_* and legacy top_ref_* columns
        new_row["ref_blastp_hit_id"] = hit["ref_blastp_hit_id"]
        new_row["ref_blastp_pident"] = hit["ref_blastp_pident"]
        new_row["ref_blastp_align_len"] = hit["ref_blastp_align_len"]
        new_row["ref_blastp_evalue"] = hit["ref_blastp_evalue"]
        new_row["ref_blastp_bitscore"] = hit["ref_blastp_bitscore"]
        new_row["top_ref_id"] = hit["ref_blastp_hit_id"]
        new_row["top_ref_pident"] = hit["ref_blastp_pident"]
        new_row["top_ref_evalue"] = hit["ref_blastp_evalue"]
        new_row["top_ref_bitscore"] = hit["ref_blastp_bitscore"]
        out_rows.append(new_row)
        
    fieldnames = list(out_rows[0].keys())
    with open(out_file, "w", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"  Wrote {len(out_rows)} rows to {out_file.name}")

def recompute_f5(checkm_info):
    print(">>> Recomputing F5 Carbonic Anhydrase Summary...")
    track_dir = PILOT_DIR / "f5_carbonic_anhydrase"
    v0_file = track_dir / "f5_candidates_summary_v0_brokenjoin.tsv"
    out_file = track_dir / "f5_candidates_summary.tsv"
    dmnd_file = track_dir / "f5_diamond_reciprocal.tsv"
    
    with open(v0_file, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    valid_pids = {r["protein_id"] for r in rows}
    
    dmnd_hits = parse_raw_diamond(dmnd_file, valid_pids)
    print(f"  Loaded {len(rows)} candidates, {len(dmnd_hits)} DIAMOND hits.")
    
    out_rows = []
    for r in rows:
        pid = r["protein_id"]
        hit = dmnd_hits.get(pid, {
            "ref_blastp_hit_id": "NO_HIT",
            "ref_blastp_pident": 0.0,
            "ref_blastp_align_len": 0,
            "ref_blastp_evalue": 1.0,
            "ref_blastp_bitscore": 0.0,
        })
        new_row = dict(r)
        new_row["ref_blastp_hit_id"] = hit["ref_blastp_hit_id"]
        new_row["ref_blastp_pident"] = hit["ref_blastp_pident"]
        new_row["ref_blastp_align_len"] = hit["ref_blastp_align_len"]
        new_row["ref_blastp_evalue"] = hit["ref_blastp_evalue"]
        new_row["ref_blastp_bitscore"] = hit["ref_blastp_bitscore"]
        new_row["top_ref_id"] = hit["ref_blastp_hit_id"]
        new_row["top_ref_pident"] = hit["ref_blastp_pident"]
        new_row["top_ref_evalue"] = hit["ref_blastp_evalue"]
        new_row["top_ref_bitscore"] = hit["ref_blastp_bitscore"]
        out_rows.append(new_row)
        
    fieldnames = list(out_rows[0].keys())
    with open(out_file, "w", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"  Wrote {len(out_rows)} rows to {out_file.name}")

def main():
    checkm_info = load_checkm_metadata()
    print(f"Loaded CheckM metadata for {len(checkm_info)} genomes.")
    recompute_f3(checkm_info)
    recompute_f4(checkm_info)
    recompute_f5(checkm_info)
    print("\nAll pilot candidate summaries successfully recomputed from raw tables!")

if __name__ == "__main__":
    main()
