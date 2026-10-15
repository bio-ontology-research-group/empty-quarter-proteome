#!/usr/bin/env python3
"""
02_extract_neighborhoods.py - Extract same-contig +/- 5kb neighborhoods and matched controls for EQ rhodopsins.

Features:
1. Loads EQ complete rhodopsin candidates (140 Type-1, 65 heliorhodopsin).
2. Parses individual MAG proteomes (drep_rep_proteomes/<mag>.faa.gz).
3. Computes exact Prodigal coordinates, contig boundaries, and flank spans.
4. Identifies contig-edge limitations (upstream/downstream truncation).
5. Extracts non-rhodopsin neighbors within +/- 5kb (excluding self).
6. Identifies matched non-overlapping control loci (same-contig preferred; same-MAG surrogate if contig is too short).
7. Emits FASTA files for:
   - rhodopsin candidates
   - rhodopsin neighborhood proteins
   - matched control proteins
   - complete contig background proteins
8. Outputs comprehensive loci_manifest.tsv and metadata audit JSON.
Supports --pilot (subset of 20 loci: 10 Type-1, 10 helio) and --full (all 205 complete loci).
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS_DIR = ROOT / "results"
METADATA_DIR = ROOT / "metadata"
for d in [RESULTS_DIR, METADATA_DIR]:
    d.mkdir(parents=True, exist_ok=True)

DEFAULT_CANDIDATES_TSV = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/family_discovery_20261015_clean/f345/results/full_v1/f3_rhodopsin/f3_candidates_summary.tsv")
DEFAULT_PROTEOMES_DIR = Path("/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/drep_rep_proteomes")

# Fallback for local testing if running in repository
LOCAL_CANDIDATES_TSV = ROOT.parent / "f345/results/full_v1/f3_rhodopsin/f3_candidates_summary.tsv"


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class ProteinRecord:
    protein_id: str
    mag: str
    contig: str
    gene_idx: int
    start: int
    stop: int
    strand: int
    partial: str
    start_type: str
    rbs_motif: str
    gc_cont: float
    seq: str


def parse_proteome_faa_gz(faa_gz_path: Path) -> List[ProteinRecord]:
    """Parse a gzipped Prodigal FASTA proteome file."""
    records = []
    current_header = ""
    current_seq = []
    
    opener = gzip.open(faa_gz_path, "rt", encoding="utf-8") if str(faa_gz_path).endswith(".gz") else open(faa_gz_path, "r", encoding="utf-8")
    with opener as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    rec = _parse_single_header(current_header, "".join(current_seq))
                    if rec:
                        records.append(rec)
                    current_seq = []
                current_header = line
            else:
                current_seq.append(line)
        if current_header:
            rec = _parse_single_header(current_header, "".join(current_seq))
            if rec:
                records.append(rec)
    return records


def _parse_single_header(header: str, seq: str) -> Optional[ProteinRecord]:
    # Format: >MAG::contig_geneidx # start # stop # strand # ID=...;partial=00;...
    try:
        parts = header[1:].split(" # ")
        protein_id = parts[0].strip()
        start = int(parts[1].strip())
        stop = int(parts[2].strip())
        strand = int(parts[3].strip())
        attr_str = parts[4].strip() if len(parts) > 4 else ""
        
        attrs = {}
        for item in attr_str.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                attrs[k.strip()] = v.strip()
        
        mag, gene_id = protein_id.split("::")
        contig, gene_idx_str = gene_id.rsplit("_", 1)
        gene_idx = int(gene_idx_str)
        
        return ProteinRecord(
            protein_id=protein_id,
            mag=mag,
            contig=contig,
            gene_idx=gene_idx,
            start=start,
            stop=stop,
            strand=strand,
            partial=attrs.get("partial", ""),
            start_type=attrs.get("start_type", ""),
            rbs_motif=attrs.get("rbs_motif", ""),
            gc_cont=float(attrs.get("gc_cont", 0.0)),
            seq=seq
        )
    except Exception as e:
        print(f"Warning: could not parse header '{header}': {e}", file=sys.stderr)
        return None


def select_pilot_subset(candidates: List[dict]) -> List[dict]:
    """Select 20 diverse representative rhodopsins: 10 Type-1, 10 Heliorhodopsin."""
    type1 = [c for c in candidates if c["assigned_class"] == "Type-1_Microbial_Rhodopsin"]
    helio = [c for c in candidates if c["assigned_class"] == "Heliorhodopsin_Type-3"]
    
    # Sort deterministically by protein_id
    type1.sort(key=lambda x: x["protein_id"])
    helio.sort(key=lambda x: x["protein_id"])
    
    # Pick 10 spread across the list
    step_t1 = max(1, len(type1) // 10)
    step_h = max(1, len(helio) // 10)
    
    pilot_t1 = [type1[i * step_t1] for i in range(10)]
    pilot_h = [helio[i * step_h] for i in range(10)]
    return pilot_t1 + pilot_h


def main():
    parser = argparse.ArgumentParser(description="Extract F3 genomic neighborhoods and matched controls")
    parser.add_argument("--mode", choices=["pilot", "full"], default="full", help="Execution mode: pilot (20 loci) or full (all 205 loci)")
    parser.add_argument("--candidates-tsv", type=Path, default=None, help="Path to f3_candidates_summary.tsv")
    parser.add_argument("--proteomes-dir", type=Path, default=DEFAULT_PROTEOMES_DIR, help="Path to drep_rep_proteomes directory")
    parser.add_argument("--flank-window", type=int, default=5000, help="Flanking distance in bp (default: 5000)")
    args = parser.parse_args()

    t0 = time.time()
    print("=" * 80)
    print(f"F3 RHODOPSIN GENOMIC CONTEXT EXTRACTION (Mode: {args.mode.upper()})")
    print(f"Flank Window: +/- {args.flank_window} bp")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print("=" * 80)

    # 1. Resolve candidates file
    cand_path = args.candidates_tsv
    if cand_path is None:
        if DEFAULT_CANDIDATES_TSV.exists():
            cand_path = DEFAULT_CANDIDATES_TSV
        elif LOCAL_CANDIDATES_TSV.exists():
            cand_path = LOCAL_CANDIDATES_TSV
        else:
            sys.exit(f"Error: candidates file not found at {DEFAULT_CANDIDATES_TSV} or {LOCAL_CANDIDATES_TSV}")

    assert cand_path.exists(), f"Missing candidates file: {cand_path}"
    cand_sha = compute_sha256(cand_path)
    print(f"Loaded candidates file: {cand_path} (SHA-256: {cand_sha[:16]}...)")

    # 2. Filter complete EQ candidates
    with open(cand_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        all_candidates = list(reader)

    eq_complete = [
        c for c in all_candidates
        if c.get("cohort") == "EQ_MAG" and str(c.get("is_complete_gene", "")).strip() in ("1", "True")
    ]
    print(f"Total rows in summary: {len(all_candidates)}")
    print(f"Filtered EQ complete rhodopsins: {len(eq_complete)}")
    
    t1_count = sum(1 for c in eq_complete if c["assigned_class"] == "Type-1_Microbial_Rhodopsin")
    h_count = sum(1 for c in eq_complete if c["assigned_class"] == "Heliorhodopsin_Type-3")
    print(f"  Type-1 Microbial Rhodopsins: {t1_count}")
    print(f"  Heliorhodopsin Type-3: {h_count}")
    assert len(eq_complete) == 205, f"Expected exactly 205 complete EQ rhodopsins, found {len(eq_complete)}"
    assert t1_count == 140 and h_count == 65, f"Expected 140 Type1 and 65 helio, got {t1_count} and {h_count}"

    # Select candidates according to mode
    if args.mode == "pilot":
        target_candidates = select_pilot_subset(eq_complete)
        print(f"\n[PILOT MODE] Selected {len(target_candidates)} representative loci (10 Type-1, 10 Helio)")
        prefix = "pilot_"
    else:
        target_candidates = eq_complete
        print(f"\n[FULL MODE] Processing all {len(target_candidates)} loci")
        prefix = "full_"

    target_mags = sorted(list({c["mag_or_genome"] for c in target_candidates}))
    all_rho_ids = {c["protein_id"] for c in all_candidates}
    print(f"Unique MAGs to load: {len(target_mags)}")

    # 3. Load proteomes for target MAGs
    print("\nLoading MAG proteomes...")
    mag_proteomes: Dict[str, List[ProteinRecord]] = {}
    mag_contigs: Dict[str, Dict[str, List[ProteinRecord]]] = {}

    for mag in target_mags:
        faa_gz = args.proteomes_dir / f"{mag}.faa.gz"
        if not faa_gz.exists():
            # Try uncompressed
            faa_gz = args.proteomes_dir / f"{mag}.faa"
        if not faa_gz.exists():
            sys.exit(f"CRITICAL: Proteome file not found for MAG '{mag}' at {faa_gz}")
        
        recs = parse_proteome_faa_gz(faa_gz)
        mag_proteomes[mag] = recs
        contig_dict: Dict[str, List[ProteinRecord]] = {}
        for r in recs:
            contig_dict.setdefault(r.contig, []).append(r)
        # Sort each contig by start position
        for ctg in contig_dict:
            contig_dict[ctg].sort(key=lambda x: x.start)
        mag_contigs[mag] = contig_dict

    total_loaded_proteins = sum(len(p) for p in mag_proteomes.values())
    print(f"Loaded {len(mag_proteomes)} MAG proteomes ({total_loaded_proteins} total proteins).")

    # 4. Extract neighborhoods and matched controls
    print(f"\nExtracting +/- {args.flank_window} bp neighborhoods and matched controls...")
    
    loci_manifest = []
    rho_proteins_out: Dict[str, ProteinRecord] = {}
    neighbor_proteins_out: Dict[str, ProteinRecord] = {}
    control_proteins_out: Dict[str, ProteinRecord] = {}
    all_contig_proteins_out: Dict[str, ProteinRecord] = {}

    for cand in target_candidates:
        pid = cand["protein_id"]
        mag = cand["mag_or_genome"]
        contig = cand["contig"]
        c_class = cand["assigned_class"]
        c_tax = cand.get("gtdbtk_taxonomy", "")
        c_start = int(cand["start"])
        c_stop = int(cand["stop"])
        c_strand = int(cand["strand"])
        c_cluster = cand.get("secondary_cluster_ani95", "")

        contig_genes = mag_contigs[mag].get(contig, [])
        if not contig_genes:
            print(f"Warning: no genes found for contig '{contig}' in MAG '{mag}'", file=sys.stderr)
            continue

        for cg in contig_genes:
            all_contig_proteins_out[cg.protein_id] = cg

        # Contig span metrics based on predicted gene coordinates
        contig_min_start = contig_genes[0].start
        contig_max_stop = contig_genes[-1].stop
        contig_gene_count = len(contig_genes)
        contig_span_bp = contig_max_stop - contig_min_start + 1

        # Flank distances
        flank_up_bp = c_start - contig_min_start
        flank_down_bp = contig_max_stop - c_stop
        is_full_flank = (flank_up_bp >= args.flank_window) and (flank_down_bp >= args.flank_window)

        if flank_up_bp >= args.flank_window and flank_down_bp >= args.flank_window:
            edge_cat = "full_flank"
        elif flank_up_bp < args.flank_window and flank_down_bp >= args.flank_window:
            edge_cat = "upstream_truncated"
        elif flank_up_bp >= args.flank_window and flank_down_bp < args.flank_window:
            edge_cat = "downstream_truncated"
        else:
            edge_cat = "both_truncated"

        # Neighborhood coordinates
        win_start = max(1, c_start - args.flank_window)
        win_stop = c_stop + args.flank_window

        # Identify neighbors (exclude self)
        neighbors = []
        tandem_rhodopsins = []
        for g in contig_genes:
            if g.protein_id == pid:
                rho_proteins_out[g.protein_id] = g
                continue
            # Check overlap with window
            if g.start <= win_stop and g.stop >= win_start:
                neighbors.append(g)
                neighbor_proteins_out[g.protein_id] = g
                if g.protein_id in all_rho_ids:
                    tandem_rhodopsins.append(g.protein_id)

        # ------------------------------------------------------------------
        # Matched control selection:
        # Candidate control genes on SAME CONTIG that do not overlap the
        # rhodopsin +/- flank window (distance >= 2 * flank_window).
        # ------------------------------------------------------------------
        eligible_same_contig_controls = []
        for g in contig_genes:
            if g.protein_id in all_rho_ids:
                continue
            # Ensure window of control gene does not overlap rhodopsin window
            ctrl_win_start = max(1, g.start - args.flank_window)
            ctrl_win_stop = g.stop + args.flank_window
            if not (ctrl_win_start <= win_stop and ctrl_win_stop >= win_start):
                eligible_same_contig_controls.append(g)

        control_gene: Optional[ProteinRecord] = None
        control_scope = "NONE"
        control_blocker = ""
        control_neighbors: List[ProteinRecord] = []

        if eligible_same_contig_controls:
            # We have valid non-overlapping control genes on the same contig!
            control_scope = "same_contig"
            # Prioritize matching full_flank status if candidate has full_flank
            if is_full_flank:
                ff_controls = [
                    g for g in eligible_same_contig_controls
                    if (g.start - contig_min_start >= args.flank_window) and (contig_max_stop - g.stop >= args.flank_window)
                ]
                pool = ff_controls if ff_controls else eligible_same_contig_controls
            else:
                pool = eligible_same_contig_controls
            
            # Pick the median control gene deterministically
            pool.sort(key=lambda x: x.start)
            control_gene = pool[len(pool) // 2]
        else:
            control_blocker = f"CONTIG_SPAN_TOO_SHORT({contig_span_bp}bp)"
            # Find a surrogate control from another contig in the SAME MAG
            all_other_contigs = [
                ctg for ctg, genes in mag_contigs[mag].items()
                if ctg != contig and not any(cg.protein_id in all_rho_ids for cg in genes)
            ]
            if all_other_contigs:
                # Pick the longest other contig in this MAG
                all_other_contigs.sort(key=lambda ctg: len(mag_contigs[mag][ctg]), reverse=True)
                surrogate_ctg = all_other_contigs[0]
                s_genes = mag_contigs[mag][surrogate_ctg]
                control_gene = s_genes[len(s_genes) // 2]
                control_scope = "same_mag_surrogate"
            else:
                control_scope = "NO_CONTROL_AVAILABLE"
                control_blocker = "NO_OTHER_CONTIGS_IN_MAG"

        if control_gene:
            control_proteins_out[control_gene.protein_id] = control_gene
            c_ctg = control_gene.contig
            c_ctg_genes = mag_contigs[mag].get(c_ctg, [])
            c_win_start = max(1, control_gene.start - args.flank_window)
            c_win_stop = control_gene.stop + args.flank_window
            for g in c_ctg_genes:
                if g.protein_id == control_gene.protein_id:
                    continue
                if g.start <= c_win_stop and g.stop >= c_win_start:
                    control_neighbors.append(g)
                    control_proteins_out[g.protein_id] = g

        loci_manifest.append({
            "rhodopsin_id": pid,
            "assigned_class": c_class,
            "mag": mag,
            "secondary_cluster_ani95": c_cluster,
            "gtdbtk_taxonomy": c_tax,
            "contig": contig,
            "rho_start": c_start,
            "rho_stop": c_stop,
            "rho_strand": c_strand,
            "contig_min_start": contig_min_start,
            "contig_max_stop": contig_max_stop,
            "contig_span_bp": contig_span_bp,
            "contig_gene_count": contig_gene_count,
            "flank_upstream_bp": flank_up_bp,
            "flank_downstream_bp": flank_down_bp,
            "is_full_flank": int(is_full_flank),
            "edge_category": edge_cat,
            "window_start": win_start,
            "window_stop": win_stop,
            "window_span_bp": win_stop - win_start + 1,
            "n_neighbors": len(neighbors),
            "neighbor_ids": ";".join(n.protein_id for n in neighbors),
            "has_tandem_rhodopsin": int(len(tandem_rhodopsins) > 0),
            "tandem_rhodopsin_ids": ";".join(tandem_rhodopsins),
            "control_scope": control_scope,
            "control_gene_id": control_gene.protein_id if control_gene else "",
            "control_contig": control_gene.contig if control_gene else "",
            "control_start": control_gene.start if control_gene else "",
            "control_stop": control_gene.stop if control_gene else "",
            "control_n_neighbors": len(control_neighbors),
            "control_neighbor_ids": ";".join(cn.protein_id for cn in control_neighbors),
            "control_blocker": control_blocker,
        })

    print(f"\nExtraction completed for {len(loci_manifest)} loci.")
    print(f"  Rhodopsin proteins: {len(rho_proteins_out)}")
    print(f"  Neighbor proteins: {len(neighbor_proteins_out)}")
    print(f"  Control proteins: {len(control_proteins_out)}")
    print(f"  All candidate contig proteins: {len(all_contig_proteins_out)}")

    # Contig edge summary
    edge_counts = {}
    for r in loci_manifest:
        edge_counts[r["edge_category"]] = edge_counts.get(r["edge_category"], 0) + 1
    print("\nContig-Edge Limitations Distribution:")
    for cat, cnt in sorted(edge_counts.items()):
        print(f"  {cat}: {cnt} ({cnt / len(loci_manifest) * 100:.1f}%)")

    ctrl_counts = {}
    for r in loci_manifest:
        ctrl_counts[r["control_scope"]] = ctrl_counts.get(r["control_scope"], 0) + 1
    print("\nControl Locus Matching Distribution:")
    for c_scope, cnt in sorted(ctrl_counts.items()):
        print(f"  {c_scope}: {cnt} ({cnt / len(loci_manifest) * 100:.1f}%)")

    # 5. Write outputs
    manifest_tsv = RESULTS_DIR / f"{prefix}loci_manifest.tsv"
    with open(manifest_tsv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(loci_manifest[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(loci_manifest)
    print(f"\nWrote loci manifest: {manifest_tsv} ({manifest_tsv.stat().st_size} bytes)")

    # Write FASTA files
    def write_fasta(out_path: Path, recs: Dict[str, ProteinRecord]):
        with open(out_path, "w", encoding="utf-8") as f:
            for pid, rec in recs.items():
                f.write(f">{rec.protein_id} # {rec.start} # {rec.stop} # {rec.strand} # partial={rec.partial};start_type={rec.start_type}\n{rec.seq}\n")
        print(f"Wrote FASTA: {out_path.name} ({len(recs)} seqs, {out_path.stat().st_size} bytes)")

    write_fasta(RESULTS_DIR / f"{prefix}rhodopsin_candidates.faa", rho_proteins_out)
    write_fasta(RESULTS_DIR / f"{prefix}neighborhood_proteins.faa", neighbor_proteins_out)
    write_fasta(RESULTS_DIR / f"{prefix}matched_control_proteins.faa", control_proteins_out)
    write_fasta(RESULTS_DIR / f"{prefix}all_contig_proteins.faa", all_contig_proteins_out)

    # Combined target FASTA for single hmmsearch run
    combined_targets: Dict[str, ProteinRecord] = {}
    combined_targets.update(rho_proteins_out)
    combined_targets.update(neighbor_proteins_out)
    combined_targets.update(control_proteins_out)
    combined_targets_path = RESULTS_DIR / f"{prefix}combined_search_proteins.faa"
    write_fasta(combined_targets_path, combined_targets)

    # Write audit JSON
    audit_json = RESULTS_DIR / f"{prefix}extraction_summary.json"
    audit_json.write_text(json.dumps({
        "mode": args.mode,
        "timestamp_utc": time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
        "flank_window_bp": args.flank_window,
        "n_loci_processed": len(loci_manifest),
        "n_type1": sum(1 for r in loci_manifest if r["assigned_class"] == "Type-1_Microbial_Rhodopsin"),
        "n_helio": sum(1 for r in loci_manifest if r["assigned_class"] == "Heliorhodopsin_Type-3"),
        "n_unique_mags": len(target_mags),
        "n_unique_contigs": len(set(r["contig"] for r in loci_manifest)),
        "edge_distribution": edge_counts,
        "control_distribution": ctrl_counts,
        "n_rho_proteins": len(rho_proteins_out),
        "n_neighbor_proteins": len(neighbor_proteins_out),
        "n_control_proteins": len(control_proteins_out),
        "n_combined_search_proteins": len(combined_targets),
        "n_all_contig_proteins": len(all_contig_proteins_out),
        "manifest_sha256": compute_sha256(manifest_tsv),
        "combined_faa_sha256": compute_sha256(combined_targets_path)
    }, indent=2), encoding="utf-8")
    print(f"Wrote extraction audit: {audit_json.name}")
    print(f"Extraction pipeline completed in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    main()
