#!/usr/bin/env python3
"""
scripts/09_f4_ectc_full.py
Track F4: Ectoine Synthase (EctC) & Pathway Synteny Full Discovery across:
1. All 5,438 Empty Quarter representative proteomes (cohort.faa, 20,451,986 proteins)
2. Harmonized Environmental SMAG cohort (smag_pool_canonical.faa)
3. Harmonized Environmental GEM cohort (gem_harmonized_pool_canonical_seed20261014.faa)

Strict Provenance & Audit Resolutions (Argon Section 4 & Root Requirements):
- Reuses verified pilot profile collection: references/profiles/f4_ectc.hmm
  - Target EctC: PF06339 (Ectoine_synth), NF009806 (ectoine_synthase)
  - Pathway markers: TIGR02406 (EctA), TIGR02407 (EctB), TIGR02408 (EctD)
  - Generic Cupin controls: PF00190 (Cupin_1), PF07883 (Cupin_2), PF02678 (Pirin)
- Never infer complete CDS from length>=180 or heuristic. Keep UNKNOWN distinct.
  For EQ MAGs: parse partial=00 from Prodigal header.
  For SMAG/GEM: look up partial from source coordinate mapping TSVs.
  Only partial == "00" is flagged is_complete_gene = 1.
- Detailed Synteny & Pathway Context:
  - Measures exact same-contig intergenic distance (bp) and strand orientation.
  - Distinguishes Complete_ABCD_Operon, Complete_ABC_Operon, Partial_Operon_BC,
    Partial_Operon_AC, Partial_Operon_CD, and Isolated_EctC.
  - Incorporates CheckM completeness to distinguish true orphans from incomplete assemblies:
    If EctA/B absent in MAG: <90% completeness -> Incomplete_MAG_cannot_exclude_unlinked_pathway.
    If EctA/B present on other contigs -> Unlinked_Pathway_Carrier.
- One-way DIAMOND blastp (--ultra-sensitive) against references/all_f345_references.dmnd.
  Query and target coverage calculated from coordinate spans.
  Evaluates match against Czech 2019 PDB 5ONN anchor A0A385TEI2 vs close homolog YP_003245677.1.
- MMseqs2 clustering (prefix_cluster.tsv) via envs/mmseqs2/bin/mmseqs. Fail-closed.
- Sequence-level deduplication across cohorts preserving cross-source memberships.
- ProstT5 input package for EctC including curated references.
- Retains all rejected and ambiguous hits with explicit reasons.
- Emits results to results/full_v1/f4_ectc/ and freezes executed code.
"""

import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
REFS_DIR = BASE_DIR / "references"
PROFILES_DIR = REFS_DIR / "profiles"
COHORTS_DIR = BASE_DIR / "cohorts"
META_DIR = BASE_DIR / "metadata"
RESULTS_DIR = BASE_DIR / "results" / "full_v1"
F4_DIR = RESULTS_DIR / "f4_ectc"
PROSTT5_DIR = RESULTS_DIR / "prostt5_inputs"
FROZEN_DIR = RESULTS_DIR / "frozen_code"

for d in [F4_DIR, PROSTT5_DIR, FROZEN_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# Profile constants (matching pilot exactly!)
F4_HMM = PROFILES_DIR / "f4_ectc.hmm"
ECTC_ACCS = {"PF06339", "NF009806"}
ECTA_ACCS = {"TIGR02406"}
ECTB_ACCS = {"TIGR02407"}
ECTD_ACCS = {"TIGR02408"}
CUPIN_ACCS = {"PF00190", "PF07883", "PF02678"}

# References & Tools
ALL_REFS_DMND = REFS_DIR / "all_f345_references.dmnd"
ALL_REFS_FAA = REFS_DIR / "all_f345_references.faa"
REFS_META_TSV = REFS_DIR / "verified_references_metadata.tsv"

MMSEQS_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/mmseqs2/bin/mmseqs")
DIAMOND_BIN = Path("/ibex/scratch/projects/c2014/rob/miniconda3/envs/eggnog/bin/diamond")

W = Path("/ibex/scratch/projects/c2014/rob/metagenomics")
WDB_CSV = W / "05_refine/_drep/out/data_tables/Wdb.csv"
EQ_COHORT_FAA = W / "06_darkmatter/family_discovery_20261015_clean/root_audit/f1_full/cohort.faa"
EQ_COHORT_TSV = W / "06_darkmatter/family_discovery_20261015_clean/root_audit/f1_full/cohort.tsv"

H3_DIR = W / "06_darkmatter/p7_analyses_clean/h3_references/results"
SMAG_CANONICAL = H3_DIR / "smag_pool_canonical.faa"
SMAG_COORD_MAP = H3_DIR / "smag_protein_coordinate_mapping.tsv"
GEM_CANONICAL = H3_DIR / "gem_harmonized_pool_canonical_seed20261014.faa"
GEM_COORD_MAP = H3_DIR / "gem_harmonized_protein_coordinate_mapping_seed20261014.tsv"

GTDB_BAC = W / "05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv"
GTDB_AR = W / "05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"
CHECKM_CSV = META_DIR / "genomeInformation.csv"


def load_metadata():
    print("Loading Wdb.csv, GTDB-Tk R232 taxonomies, and CheckM metrics...")
    assert WDB_CSV.exists(), f"CRITICAL: Missing Wdb.csv: {WDB_CSV}"
    mag_to_cluster = {}
    with open(WDB_CSV, "r", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            mag = r["genome"].removesuffix(".fa").removesuffix(".fasta")
            mag_to_cluster[mag] = r.get("cluster", "Unknown")

    tax_map = {}
    for p in [GTDB_BAC, GTDB_AR]:
        assert p.exists(), f"CRITICAL: Missing GTDB-Tk summary: {p}"
        with open(p, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                ug = r.get("user_genome", "").removesuffix(".fa").removesuffix(".fasta")
                cls = r.get("classification", "").strip()
                if ug and cls:
                    tax_map[ug] = cls

    checkm_map = {}
    if CHECKM_CSV.exists():
        with open(CHECKM_CSV, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                g = r.get("genome", "").removesuffix(".fa").removesuffix(".fasta")
                try:
                    comp = float(r.get("completeness", 0.0))
                    cont = float(r.get("contamination", 0.0))
                except (ValueError, TypeError):
                    comp, cont = 0.0, 0.0
                checkm_map[g] = {"completeness": comp, "contamination": cont}

    print(f"Loaded {len(mag_to_cluster)} clusters, {len(tax_map)} taxonomies, and {len(checkm_map)} CheckM records.")
    return mag_to_cluster, tax_map, checkm_map


def parse_domtblout(domtbl_path):
    hits = []
    with open(domtbl_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 23:
                continue
            target_name = parts[0]
            target_len = int(parts[2])
            profile_name = parts[3]
            raw_acc = parts[4]
            profile_acc = raw_acc.split(".")[0]
            qlen = int(parts[5])
            full_evalue = float(parts[6])
            full_score = float(parts[7])
            dom_ievalue = float(parts[12])
            dom_score = float(parts[13])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            hmm_cov = (hmm_to - hmm_from + 1) / float(qlen) if qlen > 0 else 0.0

            hits.append({
                "target_name": target_name,
                "target_len": target_len,
                "profile_name": profile_name,
                "profile_acc": profile_acc,
                "profile_acc_raw": raw_acc,
                "qlen": qlen,
                "full_evalue": full_evalue,
                "full_score": full_score,
                "dom_ievalue": dom_ievalue,
                "dom_score": dom_score,
                "hmm_from": hmm_from,
                "hmm_to": hmm_to,
                "ali_from": ali_from,
                "ali_to": ali_to,
                "hmm_coverage": hmm_cov,
            })
    return hits


def evaluate_f4_hits(hits_by_protein):
    ectc_candidates = {}
    pathway_markers = {"EctA": {}, "EctB": {}, "EctD": {}}
    rejections = []

    for prot_id, prot_hits in hits_by_protein.items():
        ectc_hits = [h for h in prot_hits if h["profile_acc"] in ECTC_ACCS]
        ecta_hits = [h for h in prot_hits if h["profile_acc"] in ECTA_ACCS]
        ectb_hits = [h for h in prot_hits if h["profile_acc"] in ECTB_ACCS]
        ectd_hits = [h for h in prot_hits if h["profile_acc"] in ECTD_ACCS]
        cupin_hits = [h for h in prot_hits if h["profile_acc"] in CUPIN_ACCS]

        best_ectc = max(ectc_hits, key=lambda x: x["dom_score"]) if ectc_hits else None
        best_ecta = max(ecta_hits, key=lambda x: x["dom_score"]) if ecta_hits else None
        best_ectb = max(ectb_hits, key=lambda x: x["dom_score"]) if ectb_hits else None
        best_ectd = max(ectd_hits, key=lambda x: x["dom_score"]) if ectd_hits else None
        best_cupin = max(cupin_hits, key=lambda x: x["dom_score"]) if cupin_hits else None

        score_ectc = best_ectc["dom_score"] if best_ectc else 0.0
        score_cupin = best_cupin["dom_score"] if best_cupin else 0.0

        # Record pathway markers (EctA, EctB, EctD)
        if best_ecta and best_ecta["dom_ievalue"] <= 1e-5 and best_ecta["hmm_coverage"] >= 0.50:
            pathway_markers["EctA"][prot_id] = best_ecta
        if best_ectb and best_ectb["dom_ievalue"] <= 1e-5 and best_ectb["hmm_coverage"] >= 0.50:
            pathway_markers["EctB"][prot_id] = best_ectb
        if best_ectd and best_ectd["dom_ievalue"] <= 1e-5 and best_ectd["hmm_coverage"] >= 0.50:
            pathway_markers["EctD"][prot_id] = best_ectd

        # Evaluate EctC candidate gates
        if not best_ectc:
            if best_cupin and not (best_ecta or best_ectb or best_ectd):
                rejections.append({
                    "protein_id": prot_id,
                    "reason": "Only_generic_cupin_control_hit",
                    "target_score": 0.0,
                    "control_score": score_cupin,
                    "hmm_coverage": 0.0,
                    "dom_ievalue": 1.0,
                    "profile_acc": best_cupin["profile_acc"],
                })
            continue

        cov = best_ectc["hmm_coverage"]
        ie = best_ectc["dom_ievalue"]

        # Gate 1: dom_ievalue <= 1e-10
        if ie > 1e-10:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Failed_EctC_dom_ievalue_gate_{ie:.2e}>1e-10",
                "target_score": score_ectc,
                "control_score": score_cupin,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_ectc["profile_acc"],
            })
            continue

        # Gate 2: domain coverage >= 0.70
        if cov < 0.70:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Failed_EctC_domain_coverage_gate_{cov:.2f}<0.70",
                "target_score": score_ectc,
                "control_score": score_cupin,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_ectc["profile_acc"],
            })
            continue

        # Gate 3: Score discrimination vs generic Cupins (EctC >= Cupin + 25.0 bits)
        if (score_ectc - score_cupin) < 25.0:
            rejections.append({
                "protein_id": prot_id,
                "reason": f"Excluded_generic_cupin_cross_reactivity_diff_{score_ectc - score_cupin:.1f}<25b",
                "target_score": score_ectc,
                "control_score": score_cupin,
                "hmm_coverage": cov,
                "dom_ievalue": ie,
                "profile_acc": best_ectc["profile_acc"],
            })
            continue

        # Check for profile-supported same-ORF EctC+EctD fusion with non-overlapping domain intervals
        is_fusion = 0
        fusion_c_coords = "NA"
        fusion_d_coords = "NA"
        fusion_d_score = "NA"
        fusion_d_ievalue = "NA"

        if best_ectd and best_ectd["dom_ievalue"] <= 1e-5 and best_ectd["hmm_coverage"] >= 0.50:
            c_from, c_to = best_ectc["ali_from"], best_ectc["ali_to"]
            d_from, d_to = best_ectd["ali_from"], best_ectd["ali_to"]
            # Non-overlap criterion: min of ends < max of starts
            if min(c_to, d_to) < max(c_from, d_from):
                is_fusion = 1
                fusion_c_coords = f"{c_from}..{c_to}"
                fusion_d_coords = f"{d_from}..{d_to}"
                fusion_d_score = best_ectd["dom_score"]
                fusion_d_ievalue = best_ectd["dom_ievalue"]

        ectc_candidates[prot_id] = {
            "assigned_class": "Ectoine_Synthase_EctC",
            "best_hit": best_ectc,
            "best_ectc_score": score_ectc,
            "cupin_control_score": score_cupin,
            "is_same_orf_fusion": is_fusion,
            "fusion_ectc_coords": fusion_c_coords,
            "fusion_ectd_coords": fusion_d_coords,
            "fusion_ectd_score": fusion_d_score,
            "fusion_ectd_ievalue": fusion_d_ievalue,
        }

    return ectc_candidates, pathway_markers, rejections


def stream_sequences_from_faa(faa_path, target_ids):
    seqs = {}
    with open(faa_path, "r", encoding="utf-8") as f:
        curr_id = None
        curr_hdr = ""
        curr_seq = []
        for line in f:
            if line.startswith(">"):
                if curr_id and curr_id in target_ids:
                    seqs[curr_id] = (curr_hdr, "".join(curr_seq))
                hdr = line.strip()[1:]
                curr_id = hdr.split()[0]
                curr_hdr = hdr
                curr_seq = []
            else:
                curr_seq.append(line.strip())
        if curr_id and curr_id in target_ids:
            seqs[curr_id] = (curr_hdr, "".join(curr_seq))
    return seqs


def stream_h3_coordinates(mapping_path, target_ids):
    coords = {}
    with open(mapping_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            cid = r.get("canonical_seq_id", "")
            if cid in target_ids:
                coords[cid] = {
                    "mag": r.get("mag_id") or r.get("genome_id", "Unknown"),
                    "contig": r.get("contig_id", "Unknown"),
                    "start": r.get("start_coord", "0"),
                    "stop": r.get("end_coord", "0"),
                    "strand": r.get("strand", "+"),
                    "partial": r.get("partial", "unknown"),
                }
    return coords


def calc_intergenic_dist(start1, stop1, start2, stop2):
    s1, e1 = min(start1, stop1), max(start1, stop1)
    s2, e2 = min(start2, stop2), max(start2, stop2)
    if e1 < s2:
        return s2 - e1 - 1
    elif e2 < s1:
        return s1 - e2 - 1
    else:
        return 0  # Overlapping


def run_diamond_search(candidates_faa, out_tsv):
    print("Running one-way DIAMOND blastp (--ultra-sensitive) against all_f345_references.dmnd...")
    cmd = [
        str(DIAMOND_BIN), "blastp",
        "--query", str(candidates_faa),
        "--db", str(ALL_REFS_DMND),
        "--out", str(out_tsv),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qstart", "qend", "sstart", "send", "qlen", "slen",
        "--max-target-seqs", "5",
        "--evalue", "1e-5",
        "--ultra-sensitive",
        "--threads", "4",
        "--quiet",
    ]
    subprocess.run(cmd, check=True)

    hits_by_q = {}
    if out_tsv.exists():
        with open(out_tsv, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 12:
                    qid = parts[0]
                    qstart, qend = int(parts[6]), int(parts[7])
                    sstart, send = int(parts[8]), int(parts[9])
                    qlen, slen = int(parts[10]), int(parts[11])
                    qcov = (abs(qend - qstart) + 1) / float(qlen) if qlen > 0 else 0.0
                    scov = (abs(send - sstart) + 1) / float(slen) if slen > 0 else 0.0

                    rec = {
                        "sseqid": parts[1],
                        "pident": float(parts[2]),
                        "length": int(parts[3]),
                        "evalue": float(parts[4]),
                        "bitscore": float(parts[5]),
                        "qcov": round(qcov, 4),
                        "scov": round(scov, 4),
                    }
                    hits_by_q.setdefault(qid, []).append(rec)
    return hits_by_q


def run_mmseqs_clustering(candidates_faa, prefix_path, tmp_dir):
    print("Running MMseqs2 easy-linclust at 95% sequence identity...")
    cmd = [
        str(MMSEQS_BIN), "easy-linclust",
        str(candidates_faa),
        str(prefix_path),
        str(tmp_dir),
        "--min-seq-id", "0.95",
        "-c", "0.80",
        "--cov-mode", "0",
        "--threads", "4",
    ]
    subprocess.run(cmd, check=True)

    clu_tsv = Path(f"{prefix_path}_cluster.tsv")
    assert clu_tsv.exists(), f"CRITICAL: Expected MMseqs2 cluster file {clu_tsv} not found!"
    clu_map = {}
    with open(clu_tsv, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                rep, member = parts[0], parts[1]
                clu_map[member] = rep
    assert len(clu_map) > 0, "MMseqs2 cluster mapping is empty!"
    print(f"MMseqs2 clustering complete: {len(set(clu_map.values()))} clusters across {len(clu_map)} members.")
    return clu_map


def main():
    start_time = time.time()
    print("=" * 80)
    print("TRACK F4: ECTOINE SYNTHASE (EctC) & SYNTENY FULL DISCOVERY CENSUS")
    print(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")
    print("Cohorts: 5,438 EQ MAGs, SMAG environmental, GEM environmental")
    print("=" * 80)

    # Freeze executed code BEFORE execution
    FROZEN_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(__file__, FROZEN_DIR / "executed_09_f4_full.py")

    # Assert required inputs
    assert F4_HMM.exists(), f"CRITICAL: Missing {F4_HMM}"
    assert EQ_COHORT_FAA.exists(), f"CRITICAL: Missing EQ cohort FASTA: {EQ_COHORT_FAA}"
    assert SMAG_CANONICAL.exists(), f"CRITICAL: Missing SMAG FASTA: {SMAG_CANONICAL}"
    assert GEM_CANONICAL.exists(), f"CRITICAL: Missing GEM FASTA: {GEM_CANONICAL}"

    mag_to_cluster, tax_map, checkm_map = load_metadata()

    cohort_configs = [
        {"name": "EQ_MAG", "faa": EQ_COHORT_FAA, "type": "eq"},
        {"name": "Environmental_SMAG", "faa": SMAG_CANONICAL, "type": "smag", "coord_map": SMAG_COORD_MAP},
        {"name": "Environmental_GEM", "faa": GEM_CANONICAL, "type": "gem", "coord_map": GEM_COORD_MAP},
    ]

    all_candidates = []
    all_rejections = []
    cpus = os.environ.get("SLURM_CPUS_PER_TASK", "16")

    for cfg in cohort_configs:
        cname = cfg["name"]
        faa = cfg["faa"]
        print(f"\n--- Processing Cohort: {cname} ({faa.name}) ---")

        domtbl_out = F4_DIR / f"f4_domtbl_{cname}.tsv"
        hmm_log = F4_DIR / f"f4_hmmsearch_{cname}.log"

        if domtbl_out.exists() and domtbl_out.stat().st_size > 0:
            print(f"  Reusing existing hmmsearch domtblout: {domtbl_out}")
        else:
            hmm_cmd = [
                "hmmsearch",
                "--cpu", cpus,
                "-E", "1e-3",
                "--domtblout", str(domtbl_out),
                "-o", str(hmm_log),
                str(F4_HMM),
                str(faa),
            ]
            t0 = time.time()
            subprocess.run(hmm_cmd, check=True)
            print(f"  hmmsearch complete in {time.time() - t0:.1f}s.")

        raw_hits = parse_domtblout(domtbl_out)
        print(f"  Raw domain hits: {len(raw_hits)}")

        hits_by_prot = {}
        for h in raw_hits:
            hits_by_prot.setdefault(h["target_name"], []).append(h)

        cands, pmarkers, rejs = evaluate_f4_hits(hits_by_prot)
        for r in rejs:
            r["cohort"] = cname
        all_rejections.extend(rejs)
        print(f"  EctC candidates passing gates: {len(cands)}")
        print(f"  Pathway markers detected: EctA={len(pmarkers['EctA'])}, EctB={len(pmarkers['EctB'])}, EctD={len(pmarkers['EctD'])}")

        if not cands:
            continue

        # We also need coordinates for pathway markers to measure synteny distances!
        needed_all_ids = set(cands.keys()) | set(pmarkers["EctA"].keys()) | set(pmarkers["EctB"].keys()) | set(pmarkers["EctD"].keys())

        print(f"  Extracting sequences & headers for {len(needed_all_ids)} proteins from {faa.name}...")
        seq_data = stream_sequences_from_faa(faa, needed_all_ids)

        coord_data = {}
        if cfg["type"] in ("smag", "gem"):
            print(f"  Streaming coordinates from {cfg['coord_map'].name}...")
            coord_data = stream_h3_coordinates(cfg["coord_map"], needed_all_ids)

        # Parse gene coordinates for all detected proteins (EctC and markers)
        gene_annotations = {}
        for pid in needed_all_ids:
            if pid not in seq_data:
                continue
            hdr, seq = seq_data[pid]
            if cfg["type"] == "eq":
                pmatch = re.search(r"partial=(\d{2})", hdr)
                pflag = pmatch.group(1) if pmatch else "unknown"
                fields = hdr.split(" # ")
                s_coord = int(fields[1]) if len(fields) > 1 and fields[1].isdigit() else 0
                e_coord = int(fields[2]) if len(fields) > 2 and fields[2].isdigit() else 0
                strnd = fields[3] if len(fields) > 3 else "+"
                mag = pid.split("::")[0] if "::" in pid else pid.split("_")[0]
                contig = pid.split("::")[1].rsplit("_", 1)[0] if "::" in pid else pid.rsplit("_", 1)[0]
            else:
                cd = coord_data.get(pid, {})
                pflag = cd.get("partial", "unknown")
                s_coord = int(cd.get("start", 0))
                e_coord = int(cd.get("stop", 0))
                strnd = cd.get("strand", "+")
                mag = cd.get("mag", pid.split("::")[0] if "::" in pid else "Unknown")
                contig = cd.get("contig", "Unknown")

            gene_annotations[pid] = {
                "mag": mag,
                "contig": contig,
                "start": s_coord,
                "stop": e_coord,
                "strand": strnd,
                "partial": pflag,
                "seq": seq.rstrip("*"),
                "header": hdr,
            }

        # Index pathway markers by (mag, contig) and by mag
        markers_by_contig = {"EctA": {}, "EctB": {}, "EctD": {}}
        markers_by_mag = {"EctA": set(), "EctB": set(), "EctD": set()}

        for role in ["EctA", "EctB", "EctD"]:
            for mpid in pmarkers[role]:
                if mpid in gene_annotations:
                    ginfo = gene_annotations[mpid]
                    key = (ginfo["mag"], ginfo["contig"])
                    markers_by_contig[role].setdefault(key, []).append((mpid, ginfo))
                    markers_by_mag[role].add(ginfo["mag"])

        # Process EctC candidates with granular synteny
        for pid, cinfo in cands.items():
            if pid not in gene_annotations:
                all_rejections.append({
                    "protein_id": pid,
                    "cohort": cname,
                    "reason": "Missing_in_FASTA_sequence_stream",
                    "target_score": cinfo["best_ectc_score"],
                    "control_score": cinfo["cupin_control_score"],
                    "hmm_coverage": cinfo["best_hit"]["hmm_coverage"],
                    "dom_ievalue": cinfo["best_hit"]["dom_ievalue"],
                    "profile_acc": cinfo["best_hit"]["profile_acc"],
                })
                continue

            ginfo = gene_annotations[pid]
            clean_seq = ginfo["seq"]
            seq_sha = hashlib.sha256(clean_seq.encode("utf-8")).hexdigest()
            bh = cinfo["best_hit"]
            mag = ginfo["mag"]
            contig = ginfo["contig"]
            is_comp = 1 if ginfo["partial"] == "00" else 0

            # Find nearest same-contig markers
            nearby_info = {}
            for role in ["EctA", "EctB", "EctD"]:
                c_markers = markers_by_contig[role].get((mag, contig), [])
                best_dist = 999999999
                best_mpid = "None"
                best_rel_strand = "None"

                for mpid, minfo in c_markers:
                    if mpid == pid:
                        continue
                    d = calc_intergenic_dist(ginfo["start"], ginfo["stop"], minfo["start"], minfo["stop"])
                    if d < best_dist:
                        best_dist = d
                        best_mpid = mpid
                        best_rel_strand = "Same_strand" if ginfo["strand"] == minfo["strand"] else "Opposite_strand"

                nearby_info[role] = {
                    "nearest_id": best_mpid,
                    "dist_bp": best_dist if best_dist < 999999999 else "Unlinked",
                    "rel_strand": best_rel_strand,
                    "is_colocalized": (best_dist <= 5000),
                }

            has_nearby_a = nearby_info["EctA"]["is_colocalized"]
            has_nearby_b = nearby_info["EctB"]["is_colocalized"]
            has_nearby_d = nearby_info["EctD"]["is_colocalized"]

            # Fused EctC-EctD same-ORF detection
            is_fusion = cinfo.get("is_same_orf_fusion", 0)
            fusion_c_coords = cinfo.get("fusion_ectc_coords", "NA")
            fusion_d_coords = cinfo.get("fusion_ectd_coords", "NA")
            fusion_d_score = cinfo.get("fusion_ectd_score", "NA")
            fusion_d_ievalue = cinfo.get("fusion_ectd_ievalue", "NA")

            # Granular operon / fusion classification
            if is_fusion:
                if has_nearby_a and has_nearby_b:
                    pathway_org = "Complete_ABCD_Operon_with_sameORF_EctCD_fusion"
                    is_colocalized_operon = 1
                elif has_nearby_b:
                    pathway_org = "Partial_Operon_BCD_with_sameORF_EctCD_fusion"
                    is_colocalized_operon = 1
                elif has_nearby_a:
                    pathway_org = "Partial_Operon_ACD_with_sameORF_EctCD_fusion"
                    is_colocalized_operon = 1
                else:
                    pathway_org = "Solitary_sameORF_EctCD_fusion"
                    is_colocalized_operon = 0
            else:
                if has_nearby_a and has_nearby_b and has_nearby_d:
                    pathway_org = "Complete_ABCD_Operon"
                    is_colocalized_operon = 1
                elif has_nearby_a and has_nearby_b:
                    pathway_org = "Complete_ABC_Operon"
                    is_colocalized_operon = 1
                elif has_nearby_b and has_nearby_d:
                    pathway_org = "Partial_Operon_BCD"
                    is_colocalized_operon = 1
                elif has_nearby_b:
                    pathway_org = "Partial_Operon_BC"
                    is_colocalized_operon = 1
                elif has_nearby_a:
                    pathway_org = "Partial_Operon_AC"
                    is_colocalized_operon = 1
                elif has_nearby_d:
                    pathway_org = "Partial_Operon_CD"
                    is_colocalized_operon = 1
                else:
                    pathway_org = "Isolated_EctC_Candidate"
                    is_colocalized_operon = 0

            # Whole-MAG pathway presence
            mag_has_a = mag in markers_by_mag["EctA"]
            mag_has_b = mag in markers_by_mag["EctB"]
            mag_has_d = (mag in markers_by_mag["EctD"]) or (is_fusion == 1)

            # MAG completeness & assembly limitation context
            if cfg["type"] != "eq":
                mag_comp_str = "NA"
                mag_cont_str = "NA"
                if is_colocalized_operon and ("Complete_ABCD" in pathway_org or "Complete_ABC" in pathway_org):
                    mag_pathway_status = "Operon_Co_localized"
                elif mag_has_a and mag_has_b:
                    mag_pathway_status = "Unlinked_Pathway_Carrier"
                else:
                    mag_pathway_status = "External_Cohort_CheckM_Not_Loaded_cannot_exclude_unlinked_pathway"
            else:
                chk = checkm_map.get(mag, {})
                mag_comp = chk.get("completeness", 0.0)
                mag_cont = chk.get("contamination", 0.0)
                mag_comp_str = f"{mag_comp:.1f}%"
                mag_cont_str = f"{mag_cont:.1f}%"

                if is_colocalized_operon and ("Complete_ABCD" in pathway_org or "Complete_ABC" in pathway_org):
                    mag_pathway_status = "Operon_Co_localized"
                elif mag_has_a and mag_has_b:
                    mag_pathway_status = "Unlinked_Pathway_Carrier"
                elif mag_comp < 90.0:
                    mag_pathway_status = "Incomplete_MAG_cannot_exclude_unlinked_pathway"
                else:
                    mag_pathway_status = "Unlinked_EctC_in_high_completeness_MAG_absence_conditional"

            sec_clust = mag_to_cluster.get(mag, "Unknown") if cfg["type"] == "eq" else "Environmental_Cohort"
            tax = tax_map.get(mag, "Unclassified") if cfg["type"] == "eq" else f"{cname}_Environmental"

            all_candidates.append({
                "protein_id": pid,
                "cohort": cname,
                "mag_or_genome": mag,
                "contig": contig,
                "start": ginfo["start"],
                "stop": ginfo["stop"],
                "strand": ginfo["strand"],
                "secondary_cluster_ani95": sec_clust,
                "gtdbtk_taxonomy": tax,
                "assigned_class": "Ectoine_Synthase_EctC",
                "best_profile_name": bh["profile_name"],
                "profile_accession": bh["profile_acc"],
                "domain_score": bh["dom_score"],
                "domain_ievalue": bh["dom_ievalue"],
                "hmm_coverage": round(bh["hmm_coverage"], 4),
                "best_ectc_score": cinfo["best_ectc_score"],
                "cupin_control_score": cinfo["cupin_control_score"],
                "prodigal_partial": ginfo["partial"],
                "is_complete_gene": is_comp,
                "length_aa": len(clean_seq),
                "seq_sha256": seq_sha,
                # Fusion details
                "is_same_orf_ectcd_fusion": is_fusion,
                "same_orf_ectc_domain_coords": fusion_c_coords,
                "same_orf_ectd_domain_coords": fusion_d_coords,
                "same_orf_ectd_score": fusion_d_score,
                "same_orf_ectd_ievalue": fusion_d_ievalue,
                # Nearby synteny details
                "has_nearby_ectA": int(has_nearby_a),
                "nearby_ectA_id": nearby_info["EctA"]["nearest_id"],
                "nearby_ectA_dist_bp": nearby_info["EctA"]["dist_bp"],
                "nearby_ectA_strand": nearby_info["EctA"]["rel_strand"],
                "has_nearby_ectB": int(has_nearby_b),
                "nearby_ectB_id": nearby_info["EctB"]["nearest_id"],
                "nearby_ectB_dist_bp": nearby_info["EctB"]["dist_bp"],
                "nearby_ectB_strand": nearby_info["EctB"]["rel_strand"],
                "has_nearby_ectD": int(has_nearby_d),
                "nearby_ectD_id": nearby_info["EctD"]["nearest_id"],
                "nearby_ectD_dist_bp": nearby_info["EctD"]["dist_bp"],
                "nearby_ectD_strand": nearby_info["EctD"]["rel_strand"],
                "is_colocalized_operon": is_colocalized_operon,
                "pathway_organization": pathway_org,
                "mag_has_ectA": int(mag_has_a),
                "mag_has_ectB": int(mag_has_b),
                "mag_has_ectD": int(mag_has_d),
                "mag_completeness": mag_comp_str,
                "mag_contamination": mag_cont_str,
                "mag_pathway_status": mag_pathway_status,
                "sequence": clean_seq,
            })

    print(f"\nTotal EctC candidates across all cohorts: {len(all_candidates)}")
    total_complete = sum(c["is_complete_gene"] == 1 for c in all_candidates)
    print(f"Total complete coding genes (partial=00): {total_complete}")

    # Emit combined candidate FASTA with bare protein_id headers
    candidates_faa = F4_DIR / "f4_all_candidates.faa"
    with open(candidates_faa, "w", encoding="utf-8") as f_out:
        for c in all_candidates:
            f_out.write(f">{c['protein_id']}\n{c['sequence']}\n")

    # One-way DIAMOND blastp
    diamond_tsv = F4_DIR / "f4_diamond_hits.tsv"
    if diamond_tsv.exists() and diamond_tsv.stat().st_size > 0:
        print(f"Reusing existing DIAMOND hits from {diamond_tsv}...")
        diamond_hits = {}
        with open(diamond_tsv, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 12:
                    qid = parts[0]
                    qstart, qend = int(parts[6]), int(parts[7])
                    sstart, send = int(parts[8]), int(parts[9])
                    qlen, slen = int(parts[10]), int(parts[11])
                    qcov = (abs(qend - qstart) + 1) / float(qlen) if qlen > 0 else 0.0
                    scov = (abs(send - sstart) + 1) / float(slen) if slen > 0 else 0.0
                    rec = {
                        "sseqid": parts[1],
                        "pident": float(parts[2]),
                        "length": int(parts[3]),
                        "evalue": float(parts[4]),
                        "bitscore": float(parts[5]),
                        "qcov": round(qcov, 4),
                        "scov": round(scov, 4),
                    }
                    diamond_hits.setdefault(qid, []).append(rec)
    else:
        diamond_hits = run_diamond_search(candidates_faa, diamond_tsv) if all_candidates else {}

    for c in all_candidates:
        dh = diamond_hits.get(c["protein_id"], [])
        if dh:
            top = dh[0]
            c["top_ref_id"] = top["sseqid"]
            c["top_ref_pident"] = top["pident"]
            c["top_ref_evalue"] = top["evalue"]
            c["top_ref_bitscore"] = top["bitscore"]
            c["top_ref_qcov"] = top["qcov"]
            c["top_ref_scov"] = top["scov"]
            c["top_competing_refs"] = ";".join([f"{h['sseqid']}(id={h['pident']}%,bs={h['bitscore']})" for h in dh[1:4]])
        else:
            c["top_ref_id"] = "NO_HIT"
            c["top_ref_pident"] = 0.0
            c["top_ref_evalue"] = 1.0
            c["top_ref_bitscore"] = 0.0
            c["top_ref_qcov"] = 0.0
            c["top_ref_scov"] = 0.0
            c["top_competing_refs"] = "None"

    # MMseqs2 95% clustering
    if all_candidates:
        clu_prefix = F4_DIR / "f4_clu95"
        clu_tsv = Path(f"{clu_prefix}_cluster.tsv")
        if clu_tsv.exists() and clu_tsv.stat().st_size > 0:
            print(f"Reusing existing MMseqs2 cluster TSV: {clu_tsv}...")
            clu_map = {}
            with open(clu_tsv, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        clu_map[parts[1]] = parts[0]
        else:
            tmp_mmseqs = F4_DIR / "tmp_mmseqs"
            tmp_mmseqs.mkdir(exist_ok=True)
            clu_map = run_mmseqs_clustering(candidates_faa, clu_prefix, tmp_mmseqs)
        for c in all_candidates:
            c["mmseqs_cluster_rep_95"] = clu_map.get(c["protein_id"], c["protein_id"])
    else:
        for c in all_candidates:
            c["mmseqs_cluster_rep_95"] = c["protein_id"]

    # Sequence deduplication & cross-source membership tracking
    seq_to_members = {}
    for c in all_candidates:
        sha = c["seq_sha256"]
        if sha not in seq_to_members:
            seq_to_members[sha] = {
                "sequence": c["sequence"],
                "length_aa": c["length_aa"],
                "assigned_class": c["assigned_class"],
                "cohorts": set(),
                "protein_ids": [],
                "mags": set(),
                "taxonomies": set(),
                "has_complete": False,
            }
        seq_to_members[sha]["cohorts"].add(c["cohort"])
        seq_to_members[sha]["protein_ids"].append(c["protein_id"])
        seq_to_members[sha]["mags"].add(c["mag_or_genome"])
        seq_to_members[sha]["taxonomies"].add(c["gtdbtk_taxonomy"])
        if c["is_complete_gene"] == 1:
            seq_to_members[sha]["has_complete"] = True

    for c in all_candidates:
        sm = seq_to_members[c["seq_sha256"]]
        c["shared_across_cohorts"] = ";".join(sorted(sm["cohorts"]))

    # Write Candidate Summary TSV
    summary_tsv = F4_DIR / "f4_candidates_summary.tsv"
    fieldnames = [
        "protein_id", "cohort", "mag_or_genome", "contig", "start", "stop", "strand",
        "secondary_cluster_ani95", "gtdbtk_taxonomy", "assigned_class",
        "best_profile_name", "profile_accession", "domain_score", "domain_ievalue", "hmm_coverage",
        "best_ectc_score", "cupin_control_score",
        "top_ref_id", "top_ref_pident", "top_ref_evalue", "top_ref_bitscore", "top_ref_qcov", "top_ref_scov",
        "top_competing_refs", "mmseqs_cluster_rep_95", "shared_across_cohorts",
        "prodigal_partial", "is_complete_gene", "length_aa", "seq_sha256",
        "is_same_orf_ectcd_fusion", "same_orf_ectc_domain_coords", "same_orf_ectd_domain_coords",
        "same_orf_ectd_score", "same_orf_ectd_ievalue",
        "has_nearby_ectA", "nearby_ectA_id", "nearby_ectA_dist_bp", "nearby_ectA_strand",
        "has_nearby_ectB", "nearby_ectB_id", "nearby_ectB_dist_bp", "nearby_ectB_strand",
        "has_nearby_ectD", "nearby_ectD_id", "nearby_ectD_dist_bp", "nearby_ectD_strand",
        "is_colocalized_operon", "pathway_organization", "mag_has_ectA", "mag_has_ectB", "mag_has_ectD",
        "mag_completeness", "mag_contamination", "mag_pathway_status"
    ]
    with open(summary_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_candidates)
    print(f"Saved candidates summary to {summary_tsv}")

    # Write Rejections TSV
    rejections_tsv = F4_DIR / "f4_rejections_and_ambiguous.tsv"
    rej_fields = ["protein_id", "cohort", "reason", "target_score", "control_score", "hmm_coverage", "dom_ievalue", "profile_acc"]
    with open(rejections_tsv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=rej_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rejections)
    print(f"Saved rejections and ambiguous cases to {rejections_tsv}")

    # Prepare ProstT5 Inputs: Ectoine Synthase EctC including curated references
    print("\nPreparing ProstT5 input package for EctC with curated references...")
    ref_records = []
    if ALL_REFS_FAA.exists() and REFS_META_TSV.exists():
        ref_seqs = {}
        curr_acc = None
        curr_seq = []
        with open(ALL_REFS_FAA, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if curr_acc:
                        ref_seqs[curr_acc] = "".join(curr_seq)
                    hdr = line[1:].split()[0]
                    curr_acc = hdr.split("|")[0]
                    curr_seq = []
                else:
                    curr_seq.append(line)
            if curr_acc:
                ref_seqs[curr_acc] = "".join(curr_seq)

        with open(REFS_META_TSV, "r", encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r.get("track") == "F4" and r.get("category") == "Ectoine_synthase":
                    acc = r["accession"]
                    if acc in ref_seqs:
                        r["sequence"] = ref_seqs[acc]
                        ref_records.append(r)

    assert len(ref_records) > 0, f"CRITICAL: 0 curated references matched for EctC"

    complete_fam_seqs = {sha: data for sha, data in seq_to_members.items() if data["has_complete"]}
    partial_fam_seqs = {sha: data for sha, data in seq_to_members.items() if not data["has_complete"]}

    # 1. PRIMARY ProstT5 Package: Curated references + Complete (partial=00) candidates
    fam_faa = PROSTT5_DIR / "prostt5_ectc.faa"
    fam_manifest_tsv = PROSTT5_DIR / "prostt5_ectc_manifest.tsv"

    manifest_rows = []
    with open(fam_faa, "w", encoding="utf-8") as out_f:
        # 1a. Curated references
        for ref in ref_records:
            rid = ref["accession"]
            seq = ref["sequence"]
            r_sha = hashlib.sha256(seq.encode("utf-8")).hexdigest()
            out_f.write(f">{rid}\n{seq}\n")
            manifest_rows.append({
                "seq_id": rid,
                "seq_sha256": r_sha,
                "length_aa": len(seq),
                "source_type": "Curated_Reference",
                "member_count": 1,
                "member_ids": rid,
                "carrier_mag_or_genomes": ref.get("organism", "Reference"),
                "taxonomies": ref.get("lineage", "Reference"),
                "assigned_family": "Ectoine_Synthase_EctC",
                "is_complete_gene": 1,
            })

        # 1b. Complete candidate sequences (deduplicated)
        for sha, sdata in complete_fam_seqs.items():
            rep_id = sdata["protein_ids"][0]
            out_f.write(f">{rep_id}\n{sdata['sequence']}\n")
            src_type = ";".join(sorted(sdata["cohorts"]))
            manifest_rows.append({
                "seq_id": rep_id,
                "seq_sha256": sha,
                "length_aa": sdata["length_aa"],
                "source_type": src_type,
                "member_count": len(sdata["protein_ids"]),
                "member_ids": ";".join(sdata["protein_ids"]),
                "carrier_mag_or_genomes": ";".join(sorted(sdata["mags"])),
                "taxonomies": ";".join(sorted(sdata["taxonomies"])),
                "assigned_family": "Ectoine_Synthase_EctC",
                "is_complete_gene": 1,
            })

    m_fields = [
        "seq_id", "seq_sha256", "length_aa", "source_type", "member_count",
        "member_ids", "carrier_mag_or_genomes", "taxonomies", "assigned_family", "is_complete_gene"
    ]
    with open(fam_manifest_tsv, "w", newline="", encoding="utf-8") as mf:
        w = csv.DictWriter(mf, fieldnames=m_fields, delimiter="\t")
        w.writeheader()
        w.writerows(manifest_rows)

    # Validate PRIMARY FASTA ID closure
    fasta_ids = []
    with open(fam_faa, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                fasta_ids.append(line[1:].strip().split()[0])
    manifest_ids = [m["seq_id"] for m in manifest_rows]
    assert fasta_ids == manifest_ids, f"ID closure failure for primary EctC: {len(fasta_ids)} in FASTA vs {len(manifest_ids)} in manifest"
    assert len(set(manifest_ids)) == len(manifest_ids), f"Duplicate seq_ids in primary manifest for EctC"

    # 2. SEPARATE Partial Candidates Package (retained, NOT silently dropped)
    part_faa = PROSTT5_DIR / "prostt5_ectc_partial_excluded.faa"
    part_manifest_tsv = PROSTT5_DIR / "prostt5_ectc_partial_excluded_manifest.tsv"
    part_rows = []
    with open(part_faa, "w", encoding="utf-8") as out_pf:
        for sha, sdata in partial_fam_seqs.items():
            rep_id = sdata["protein_ids"][0]
            out_pf.write(f">{rep_id}\n{sdata['sequence']}\n")
            src_type = ";".join(sorted(sdata["cohorts"]))
            part_rows.append({
                "seq_id": rep_id,
                "seq_sha256": sha,
                "length_aa": sdata["length_aa"],
                "source_type": src_type,
                "member_count": len(sdata["protein_ids"]),
                "member_ids": ";".join(sdata["protein_ids"]),
                "carrier_mag_or_genomes": ";".join(sorted(sdata["mags"])),
                "taxonomies": ";".join(sorted(sdata["taxonomies"])),
                "assigned_family": "Ectoine_Synthase_EctC",
                "is_complete_gene": 0,
                "exclusion_reason": "Incomplete_gene_fragment_partial_not_00",
            })
    p_fields = m_fields + ["exclusion_reason"]
    with open(part_manifest_tsv, "w", newline="", encoding="utf-8") as pf:
        w = csv.DictWriter(pf, fieldnames=p_fields, delimiter="\t")
        w.writeheader()
        w.writerows(part_rows)

    part_fasta_ids = []
    with open(part_faa, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                part_fasta_ids.append(line[1:].strip().split()[0])
    part_manifest_ids = [m["seq_id"] for m in part_rows]
    assert part_fasta_ids == part_manifest_ids, f"ID closure failure for partial EctC"

    print(f"Validated PRIMARY ProstT5 package for EctC: {len(manifest_rows)} sequences ({len(ref_records)} refs + {len(complete_fam_seqs)} complete candidates) to {fam_faa.name}")
    print(f"Retained {len(part_rows)} partial candidates in separate table {part_manifest_tsv.name}")

    # Emit Metrics JSON
    operon_counts = {}
    for c in all_candidates:
        org = c["pathway_organization"]
        operon_counts[org] = operon_counts.get(org, 0) + 1

    metrics = {
        "track": "F4",
        "family": "Ectoine_Synthase_EctC",
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_candidates": len(all_candidates),
        "total_complete_candidates": total_complete,
        "by_cohort": {
            "EQ_MAG": sum(1 for c in all_candidates if c["cohort"] == "EQ_MAG"),
            "Environmental_SMAG": sum(1 for c in all_candidates if c["cohort"] == "Environmental_SMAG"),
            "Environmental_GEM": sum(1 for c in all_candidates if c["cohort"] == "Environmental_GEM"),
        },
        "complete_by_cohort": {
            "EQ_MAG": sum(1 for c in all_candidates if c["cohort"] == "EQ_MAG" and c["is_complete_gene"] == 1),
            "Environmental_SMAG": sum(1 for c in all_candidates if c["cohort"] == "Environmental_SMAG" and c["is_complete_gene"] == 1),
            "Environmental_GEM": sum(1 for c in all_candidates if c["cohort"] == "Environmental_GEM" and c["is_complete_gene"] == 1),
        },
        "pathway_organization_counts": operon_counts,
        "same_orf_ectcd_fusions": sum(1 for c in all_candidates if c.get("is_same_orf_ectcd_fusion") == 1),
        "solitary_same_orf_fusions": sum(1 for c in all_candidates if c.get("pathway_organization") == "Solitary_sameORF_EctCD_fusion"),
        "unique_sequences": len(seq_to_members),
        "primary_prostt5_candidates": len(complete_fam_seqs),
        "partial_excluded_candidates": len(partial_fam_seqs),
        "total_rejections_and_ambiguous": len(all_rejections),
        "execution_time_seconds": round(time.time() - start_time, 2),
    }
    metrics_json = F4_DIR / "f4_full_metrics.json"
    with open(metrics_json, "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nSaved metrics to {metrics_json}")

    # Freeze code
    shutil.copy2(Path(__file__), FROZEN_DIR / "executed_09_f4_full.py")
    print(f"Code frozen to {FROZEN_DIR / 'executed_09_f4_full.py'}")
    print("=" * 80)
    print("TRACK F4 FULL DISCOVERY COMPLETE.")
    print("=" * 80)


if __name__ == "__main__":
    main()
