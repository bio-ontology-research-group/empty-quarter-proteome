#!/usr/bin/env python3
"""Bounded Slurm extraction of existing contig/bin/marker/CRISPR provenance for Experiment 3 host MAGs.

Runs on a Shaheen compute node via Slurm (`sbatch_extract_host_array_provenance.sh`) to extract:
1. Exact FASTA contig lengths, GC content, length ranks, raw & canonical SHA-256s, and literal
   array/spacer substring verification for the 6 Stage 4 host MAGs (4 primary credible hosts,
   1 boundary-sensitive provisional host 43PRr3__self_SemiBin_469, and 1 unresolved sensitivity
   host 41Sr3__self_SemiBin_108).
2. Cross-MAG contig membership (`cellular_contig_to_mag.tsv.gz`) and competitive reference
   ownership (`competitive_reference_crosswalk.tsv`: unique vs shared/ambiguous contigs and bp).
3. Authoritative Stage 1 CRISPR arrays, Cas operons (same-contig vs other-contig), canonical
   spacers, and Stage 2 spacer-to-target alignments.
4. MAG-level GTDB r232 taxonomy, CheckM2 completeness/contamination, and GUNC chimerism/clade
   separation metrics (`c3_gunc_per_genome.tsv` and `batch_XX.detail.tar.gz` `all_levels.tsv`).
5. Existing contig-level geNomad marker taxonomy (`batch_taxonomy.tsv`), gene-level marker calls
   and taxonomic breakdown (`batch_genes.tsv`), chromosome/plasmid/virus aggregated scores
   (`batch_aggregated_classification.tsv`), provirus calls (`batch_provirus.tsv`), and
   DefenseFinder/PADLOC systems from `/project/k10018/work/eq-annot-20261030/results/c8_defense/`.
6. Existing gene catalogue darkness/membership stats on array-bearing contigs vs remaining MAG
   from `/scratch/project/k10018/eq-archaeal-function-20261106/inputs/archaeal_genes.tsv`.
7. Any completed per-sample competitive mapping contig coverage (`per_sample/<sample>_contig_coverage.tsv`)
   if available, or explicit pending status with pre-verified reference ownership crosswalk.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import os
import tarfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

EXP3_ROOT = Path(os.environ.get("EXP3_ROOT", "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses"))
ANNOT_ROOT = Path("/project/k10018/work/eq-annot-20261030")
ARCH_FN_ROOT = Path("/scratch/project/k10018/eq-archaeal-function-20261106")

FOCAL_HOSTS: dict[str, dict[str, Any]] = {
    "41PRr2__single_SemiBin_236": {
        "host_role": "primary_credible_host",
        "source_sample": "41PRr2",
        "paired_targets": ["V60PRr2::V60PRr2_715976"],
    },
    "41PRr2__single_SemiBin_97": {
        "host_role": "primary_credible_host",
        "source_sample": "41PRr2",
        "paired_targets": ["55Dr2::55Dr2_75333", "9Sr2::9Sr2_284237"],
    },
    "S49Dr1__self_SemiBin_1106": {
        "host_role": "primary_credible_host",
        "source_sample": "S49Dr1",
        "paired_targets": ["48Sr2::48Sr2_190170", "coasm_48_S::k127_1210458"],
    },
    "V48Dr1__self_SemiBin_3215": {
        "host_role": "primary_credible_host",
        "source_sample": "V48Dr1",
        "paired_targets": ["coasm_48_PR::k127_4172870"],
    },
    "43PRr3__self_SemiBin_469": {
        "host_role": "provisional_boundary_sensitive_host",
        "source_sample": "43PRr3",
        "paired_targets": ["44Dr2::44Dr2_4002037"],
    },
    "41Sr3__self_SemiBin_108": {
        "host_role": "sensitivity_unresolved_virality_host",
        "source_sample": "41Sr3",
        "paired_targets": ["41Sr3::41Sr3_56779"],
    },
}

ALL_PAIRED_TARGETS = {
    t for meta in FOCAL_HOSTS.values() for t in meta["paired_targets"]
}

DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1]


def sha256_str(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def canonical_seq_sha256(seq: str) -> str:
    up = seq.upper()
    rc = revcomp(up)
    return sha256_str(up if up <= rc else rc)


def stream_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    curr: str | None = None
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                curr = line[1:].split()[0]
                seqs[curr] = []
            elif curr is not None:
                seqs[curr].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def compute_n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2.0
    running = 0
    for length in s:
        running += length
        if running >= half:
            return length
    return s[-1]


def read_tsv_filtered(path: Path, filter_col: str, allowed: set[str]) -> list[dict[str, str]]:
    out: list[dict[str, str]] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get(filter_col, "") in allowed:
                out.append(dict(row))
    return out


def main() -> int:
    source_hashes: dict[str, dict[str, Any]] = {}

    def record_source(label: str, path: Path, compute_hash: bool = True) -> None:
        if path.exists():
            source_hashes[label] = {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": stream_sha256(path) if compute_hash else "large_archive_skipped_full_stream",
            }

    host_set = set(FOCAL_HOSTS.keys())

    # 1. Stage 1 CRISPR arrays, Cas operons, verified spacers, and Stage 2 alignments
    arrays_path = EXP3_ROOT / "results/stage1_raw_crispr_arrays.tsv"
    cas_path = EXP3_ROOT / "results/stage1_cas_operons.tsv"
    spacers_path = EXP3_ROOT / "results/stage1_verified_archaeal_spacers.tsv"
    align_path = EXP3_ROOT / "results/stage2_all_spacer_alignments_scored.tsv"
    s0_ver_path = EXP3_ROOT / "results/stage0_arch591_fasta_verification.tsv"
    s0_dup_path = EXP3_ROOT / "results/stage0_arch591_sample_duplicate_contig_audit.tsv"
    gm_path = EXP3_ROOT / "inputs/genome_manifest.tsv"
    ref_rep_path = EXP3_ROOT / "inputs/refined_representatives.tsv"
    tax_all_path = EXP3_ROOT / "inputs/taxonomy_reconciliation_all5438.tsv.gz"
    tax_th_path = EXP3_ROOT / "inputs/taxonomy_reconciliation_thermo237.tsv"

    for lbl, p in [
        ("stage1_raw_crispr_arrays", arrays_path),
        ("stage1_cas_operons", cas_path),
        ("stage1_verified_archaeal_spacers", spacers_path),
        ("stage2_all_spacer_alignments_scored", align_path),
        ("stage0_arch591_fasta_verification", s0_ver_path),
        ("stage0_arch591_sample_duplicate_contig_audit", s0_dup_path),
        ("genome_manifest", gm_path),
        ("refined_representatives", ref_rep_path),
        ("taxonomy_reconciliation_all5438", tax_all_path),
        ("taxonomy_reconciliation_thermo237", tax_th_path),
    ]:
        record_source(lbl, p, compute_hash=True)

    arrays_rows = read_tsv_filtered(arrays_path, "mag_id", host_set)
    cas_rows = read_tsv_filtered(cas_path, "mag_id", host_set)
    spacer_rows = read_tsv_filtered(spacers_path, "mag_id", host_set)

    # Filter alignments where host_mag is in host_set and target_id is in ALL_PAIRED_TARGETS
    pair_alignments: list[dict[str, str]] = []
    with align_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("host_mag") in host_set and row.get("target_id") in ALL_PAIRED_TARGETS:
                pair_alignments.append(dict(row))

    # Map spacer_id -> spacer row
    spacer_by_id = {r["unique_spacer_id"]: r for r in spacer_rows}

    # Identify all array-bearing contigs and viral-hit array-bearing contigs per host
    all_array_contigs_by_host: dict[str, set[str]] = defaultdict(set)
    for r in arrays_rows:
        all_array_contigs_by_host[r["mag_id"]].add(r["raw_contig_id"])

    viral_hit_contigs_by_host: dict[str, set[str]] = defaultdict(set)
    for aln in pair_alignments:
        h = aln["host_mag"]
        ns_ctg = aln.get("host_namespaced_contig_id", "")
        raw_ctg = ns_ctg.split("::", 1)[-1] if "::" in ns_ctg else ns_ctg
        if raw_ctg:
            viral_hit_contigs_by_host[h].add(raw_ctg)

    # 2. Parse FASTAs for the 6 host MAGs and compute exact contig metrics & substring checks
    fasta_data_by_host: dict[str, dict[str, Any]] = {}
    focal_contig_raw_ids: set[str] = set()
    focal_contig_namespaced_ids: set[str] = set()
    focal_contig_canonical_shas: dict[str, str] = {}

    for host_mag, hmeta in FOCAL_HOSTS.items():
        fa_path = EXP3_ROOT / "inputs/arch591_fastas" / f"{host_mag}.fa"
        record_source(f"fasta_{host_mag}", fa_path, compute_hash=True)
        seqs = parse_fasta(fa_path)
        lengths_sorted = sorted((len(s) for s in seqs.values()), reverse=True)
        total_bp = sum(lengths_sorted)
        n50 = compute_n50(lengths_sorted)

        # Rank map (1-based rank by descending length)
        contig_ranks: dict[str, int] = {}
        for idx, (cid, cseq) in enumerate(sorted(seqs.items(), key=lambda kv: (-len(kv[1]), kv[0])), start=1):
            contig_ranks[cid] = idx

        target_contigs = sorted(all_array_contigs_by_host[host_mag] | viral_hit_contigs_by_host[host_mag])
        contig_details: dict[str, Any] = {}
        for cid in target_contigs:
            cseq = seqs[cid]
            c_up = cseq.upper()
            gc_count = c_up.count("G") + c_up.count("C")
            raw_sha = sha256_str(c_up)
            can_sha = canonical_seq_sha256(c_up)
            ns_id = f"{hmeta['source_sample']}::{cid}"
            focal_contig_raw_ids.add(cid)
            focal_contig_namespaced_ids.add(ns_id)
            focal_contig_canonical_shas[ns_id] = can_sha

            # Verify CRISPR array intervals on this contig
            ctg_arrays = [a for a in arrays_rows if a["mag_id"] == host_mag and a["raw_contig_id"] == cid]
            verified_intervals = []
            for a in ctg_arrays:
                astart = int(a["array_start"])
                aend = int(a["array_end"])
                in_bounds = (1 <= astart < aend <= len(c_up))
                subseq = c_up[astart - 1 : aend] if in_bounds else ""
                rep = a["consensus_repeat"].upper()
                rep_rc = revcomp(rep)
                rep_occurrences = subseq.count(rep) + (subseq.count(rep_rc) if rep_rc != rep else 0)
                verified_intervals.append({
                    "array_id": a["array_id"],
                    "caller": a["caller"],
                    "subtype_or_label": a["subtype_or_label"],
                    "array_trust_tier": a["array_trust_tier"],
                    "array_start_1based": astart,
                    "array_end_1based": aend,
                    "array_span_bp": aend - astart + 1,
                    "in_contig_bounds": in_bounds,
                    "distance_from_contig_5p_bp": astart - 1,
                    "distance_from_contig_3p_bp": len(c_up) - aend,
                    "interval_seq_sha256": sha256_str(subseq) if in_bounds else None,
                    "consensus_repeat": rep,
                    "exact_consensus_repeat_occurrences_in_interval": rep_occurrences,
                    "n_repeats_reported": int(a["n_repeats"]),
                    "n_raw_spacers_reported": int(a["n_raw_spacers"]),
                    "n_valid_spacers_reported": int(a["n_valid_spacers"]),
                })

            # Verify viral-matching spacers on this contig
            ctg_spacers_verified = {}
            for aln in pair_alignments:
                if aln["host_mag"] != host_mag:
                    continue
                sp_id = aln["unique_spacer_id"]
                sp_row = spacer_by_id.get(sp_id)
                if not sp_row or sp_row["primary_namespaced_contig_id"] != ns_id:
                    continue
                if sp_id in ctg_spacers_verified:
                    continue
                sp_seq = sp_row["representative_spacer_seq"].upper()
                sp_rc = revcomp(sp_seq)
                astart = int(sp_row["primary_array_start"])
                aend = int(sp_row["primary_array_end"])
                arr_sub = c_up[astart - 1 : aend]
                in_arr_fwd = sp_seq in arr_sub
                in_arr_rc = sp_rc in arr_sub
                pos_in_contig_1based = (c_up.find(sp_seq) + 1) if sp_seq in c_up else ((c_up.find(sp_rc) + 1) if sp_rc in c_up else None)
                ctg_spacers_verified[sp_id] = {
                    "unique_spacer_id": sp_id,
                    "representative_spacer_seq": sp_seq,
                    "canonical_spacer_seq_rc": sp_row["canonical_spacer_seq_rc"],
                    "spacer_len": int(sp_row["spacer_len"]),
                    "callers": sp_row["callers"],
                    "array_trust_tiers": sp_row["array_trust_tiers"],
                    "primary_array_start": astart,
                    "primary_array_end": aend,
                    "all_array_loci": sp_row["all_array_loci"],
                    "repeat_seq": sp_row["repeat_seq"],
                    "verified_literal_in_array_interval": bool(in_arr_fwd or in_arr_rc),
                    "strand_in_contig": "+" if in_arr_fwd else ("-" if in_arr_rc else "not_found"),
                    "contig_1based_start": pos_in_contig_1based,
                }

            contig_details[cid] = {
                "raw_contig_id": cid,
                "namespaced_contig_id": ns_id,
                "length_bp": len(c_up),
                "gc_fraction": round(gc_count / len(c_up), 6) if c_up else 0.0,
                "length_rank_in_mag": contig_ranks[cid],
                "mag_total_contigs": len(seqs),
                "fraction_of_mag_bases": round(len(c_up) / total_bp, 6) if total_bp else 0.0,
                "raw_seq_sha256": raw_sha,
                "canonical_seq_sha256": can_sha,
                "carries_viral_matching_spacer": cid in viral_hit_contigs_by_host[host_mag],
                "verified_crispr_arrays": verified_intervals,
                "verified_viral_matching_spacers": list(ctg_spacers_verified.values()),
            }

        fasta_data_by_host[host_mag] = {
            "fasta_sha256": source_hashes[f"fasta_{host_mag}"]["sha256"],
            "n_contigs": len(seqs),
            "total_bases": total_bp,
            "n50_bp": n50,
            "max_contig_bp": lengths_sorted[0] if lengths_sorted else 0,
            "min_contig_bp": lengths_sorted[-1] if lengths_sorted else 0,
            "array_bearing_contigs": contig_details,
            "all_contig_canonical_shas": {cid: canonical_seq_sha256(s) for cid, s in seqs.items()},
            "all_contig_lengths": {cid: len(s) for cid, s in seqs.items()},
        }

    # 3. Check cellular_contig_to_mag.tsv.gz for cross-MAG sharing of array-bearing contigs
    cell_map_path = EXP3_ROOT / "inputs/cellular_contig_to_mag.tsv.gz"
    record_source("cellular_contig_to_mag", cell_map_path, compute_hash=True)
    contig_to_mags: dict[str, list[str]] = defaultdict(list)
    with gzip.open(cell_map_path, "rt", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # Check columns
            cid = row.get("contig_id") or row.get("namespaced_contig_id") or row.get("contig") or ""
            raw_cid = cid.split("::", 1)[-1] if "::" in cid else cid
            if raw_cid in focal_contig_raw_ids or cid in focal_contig_namespaced_ids:
                mag = row.get("mag_id") or row.get("mag") or row.get("genome") or ""
                contig_to_mags[raw_cid].append(mag)

    # 4. Stream stage4_reference_sequence_crosswalk.tsv once to get ownership of all contigs in the 6 host MAGs
    crosswalk_path = EXP3_ROOT / "results/stage4_competitive_reference/stage4_reference_sequence_crosswalk.tsv"
    record_source("stage4_reference_sequence_crosswalk", crosswalk_path, compute_hash=True)

    all_host_canonical_shas: set[str] = set()
    for hdata in fasta_data_by_host.values():
        all_host_canonical_shas.update(hdata["all_contig_canonical_shas"].values())

    sha_to_crosswalk_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    with crosswalk_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            csha = row["canonical_sha256"]
            if csha in all_host_canonical_shas:
                sha_to_crosswalk_rows[csha].append({
                    "representative_header": row["representative_header"],
                    "source_category": row["source_category"],
                    "entity_id": row["entity_id"],
                    "contig_id": row["contig_id"],
                    "original_header": row["original_header"],
                    "seq_length": row["seq_length"],
                    "is_reference_representative": row["is_reference_representative"],
                    "dedup_reason": row["dedup_reason"],
                    "ambiguity_status": row["ambiguity_status"],
                })

    reference_ownership_by_host: dict[str, Any] = {}
    for host_mag, hdata in fasta_data_by_host.items():
        n_unique_contigs = 0
        n_ambig_contigs = 0
        unique_bp = 0
        ambig_bp = 0
        array_contig_ownership: dict[str, Any] = {}
        for cid, csha in hdata["all_contig_canonical_shas"].items():
            clen = hdata["all_contig_lengths"][cid]
            crows = sha_to_crosswalk_rows.get(csha, [])
            entities = sorted({(r["source_category"], r["entity_id"]) for r in crows})
            is_unique = (len(entities) == 1 and entities[0] == ("host_mag", host_mag))
            if is_unique:
                n_unique_contigs += 1
                unique_bp += clen
            else:
                n_ambig_contigs += 1
                ambig_bp += clen
            if cid in hdata["array_bearing_contigs"]:
                array_contig_ownership[cid] = {
                    "canonical_sha256": csha,
                    "representative_header": crows[0]["representative_header"] if crows else None,
                    "ownership_status": "unique_to_host_mag" if is_unique else "shared_ambiguous",
                    "crosswalk_ambiguity_status": crows[0]["ambiguity_status"] if crows else None,
                    "crosswalk_entries": crows,
                }
        reference_ownership_by_host[host_mag] = {
            "mag_total_contigs": hdata["n_contigs"],
            "mag_unique_contigs": n_unique_contigs,
            "mag_ambiguous_contigs": n_ambig_contigs,
            "mag_total_bp": hdata["total_bases"],
            "mag_unique_bp": unique_bp,
            "mag_ambiguous_bp": ambig_bp,
            "array_contig_ownership": array_contig_ownership,
        }

    # 5. Extract MAG-level metadata (GTDB r232, CheckM2, GUNC per genome, GUNC all_levels)
    s0_ver_rows = {r["genome"]: r for r in read_tsv_filtered(s0_ver_path, "genome", host_set)}
    s0_dup_rows = {r["genome"]: r for r in read_tsv_filtered(s0_dup_path, "genome", host_set)}
    ref_rep_rows = {r["mag"]: r for r in read_tsv_filtered(ref_rep_path, "mag", host_set)}
    tax_all_rows = {r["genome"]: r for r in read_tsv_filtered(tax_all_path, "genome", host_set)}

    gunc_per_genome_path = ARCH_FN_ROOT / "inputs/c3_gunc_per_genome.tsv"
    genome_qual_path = ARCH_FN_ROOT / "inputs/genome_quality.tsv"
    record_source("c3_gunc_per_genome", gunc_per_genome_path, compute_hash=True)
    record_source("genome_quality", genome_qual_path, compute_hash=True)
    gunc_rows = {r["genome"]: r for r in read_tsv_filtered(gunc_per_genome_path, "genome", host_set)}
    qual_rows = {r["genome"]: r for r in read_tsv_filtered(genome_qual_path, "genome", host_set)}

    # Find exact batch index (00..15) for each host MAG from /project/k10018/work/eq-annot-20261030/inputs/genomes.txt
    genomes_txt_path = ANNOT_ROOT / "inputs/genomes.txt"
    record_source("eq_annot_genomes_txt", genomes_txt_path, compute_hash=True)
    host_to_batch: dict[str, str] = {}
    with genomes_txt_path.open("r", encoding="utf-8") as fh:
        for idx, line in enumerate(fh, start=1):
            g = line.strip()
            if g in host_set:
                host_to_batch[g] = f"{idx % 16:02d}"

    batches_needed = sorted(set(host_to_batch.values()))

    # Extract GUNC all_levels.tsv and any contig_assignments.tsv from batch_XX.detail.tar.gz for the needed batches
    gunc_all_levels_by_host: dict[str, list[dict[str, str]]] = {}
    gunc_contig_assignments_by_host: dict[str, list[dict[str, str]]] = defaultdict(list)
    gunc_tar_members_by_host: dict[str, list[str]] = defaultdict(list)
    for b in batches_needed:
        gunc_tar_path = ANNOT_ROOT / f"results/c3_gunc/batch_{b}.detail.tar.gz"
        record_source(f"c3_gunc_batch_{b}_detail_tar", gunc_tar_path, compute_hash=True)
        batch_hosts = {h for h, hb in host_to_batch.items() if hb == b}
        with tarfile.open(gunc_tar_path, "r:gz") as tf:
            for m in tf.getmembers():
                for host_mag in batch_hosts:
                    if host_mag in m.name:
                        gunc_tar_members_by_host[host_mag].append(m.name)
                        fobj = tf.extractfile(m)
                        if fobj is None:
                            continue
                        text = fobj.read().decode("utf-8")
                        rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
                        if m.name.endswith(".all_levels.tsv"):
                            gunc_all_levels_by_host[host_mag] = rows
                        elif "contig" in m.name:
                            # Filter to array-bearing contigs or keep all if compact
                            for r in rows:
                                cid = r.get("contig", "")
                                if cid in focal_contig_raw_ids:
                                    gunc_contig_assignments_by_host[host_mag].append(r)

    # 6. Extract contig-level geNomad taxonomy, gene markers, aggregated scores, proviruses, and DefenseFinder/PADLOC
    genomad_contig_taxonomy: dict[str, dict[str, str]] = {}
    genomad_contig_scores: dict[str, dict[str, str]] = {}
    genomad_proviruses_by_host: dict[str, list[dict[str, str]]] = defaultdict(list)
    genomad_array_contig_genes: dict[str, list[dict[str, str]]] = defaultdict(list)
    genomad_mag_gene_tax_counts: dict[str, Counter] = defaultdict(Counter)
    genomad_array_contig_gene_tax_counts: dict[str, Counter] = defaultdict(Counter)
    defensefinder_by_host: dict[str, list[dict[str, str]]] = defaultdict(list)
    padloc_by_host: dict[str, list[dict[str, str]]] = defaultdict(list)

    for b in batches_needed:
        batch_hosts = {h for h, hb in host_to_batch.items() if hb == b}
        # DefenseFinder systems & genes
        df_sys_path = ANNOT_ROOT / f"results/c8_defense/batch_{b}.defensefinder_systems.tsv.gz"
        df_gen_path = ANNOT_ROOT / f"results/c8_defense/batch_{b}.defensefinder_genes.tsv.gz"
        pad_path = ANNOT_ROOT / f"results/c8_defense/batch_{b}.padloc.csv.gz"
        for lbl, p in [
            (f"c8_batch_{b}_defensefinder_systems", df_sys_path),
            (f"c8_batch_{b}_defensefinder_genes", df_gen_path),
            (f"c8_batch_{b}_padloc", pad_path),
        ]:
            record_source(lbl, p, compute_hash=True)

        if df_sys_path.exists() and df_sys_path.stat().st_size > 100:
            with gzip.open(df_sys_path, "rt", encoding="utf-8") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if parts and parts[0] in batch_hosts:
                        defensefinder_by_host[parts[0]].append({"raw_line": line.rstrip("\n")})

        if pad_path.exists() and pad_path.stat().st_size > 100:
            with gzip.open(pad_path, "rt", encoding="utf-8") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split(",")
                    if parts and parts[0] in batch_hosts:
                        padloc_by_host[parts[0]].append({"raw_line": line.rstrip("\n")})

        # geNomad tar archive (stream sequentially via r|gz for speed)
        gn_tar_path = ANNOT_ROOT / f"results/c8_defense/batch_{b}.genomad.tar.gz"
        record_source(f"c8_batch_{b}_genomad_tar", gn_tar_path, compute_hash=False)
        wanted_members = {
            "./batch_annotate/batch_taxonomy.tsv",
            "./batch_annotate/batch_genes.tsv",
            "./batch_aggregated_classification/batch_aggregated_classification.tsv",
            "./batch_find_proviruses/batch_provirus.tsv",
        }
        with tarfile.open(gn_tar_path, "r|gz") as tf:
            for member in tf:
                if member.name not in wanted_members:
                    continue
                fobj = tf.extractfile(member)
                if fobj is None:
                    continue
                reader = csv.DictReader((line.decode("utf-8") for line in fobj), delimiter="\t")
                if member.name == "./batch_annotate/batch_taxonomy.tsv":
                    for row in reader:
                        seq_name = row.get("seq_name", "")
                        if "|" not in seq_name:
                            continue
                        g, raw_cid = seq_name.split("|", 1)
                        if g in batch_hosts:
                            genomad_contig_taxonomy[seq_name] = dict(row)
                elif member.name == "./batch_aggregated_classification/batch_aggregated_classification.tsv":
                    for row in reader:
                        seq_name = row.get("seq_name", "")
                        if "|" not in seq_name:
                            continue
                        g, raw_cid = seq_name.split("|", 1)
                        if g in batch_hosts:
                            genomad_contig_scores[seq_name] = dict(row)
                elif member.name == "./batch_find_proviruses/batch_provirus.tsv":
                    for row in reader:
                        seq_name = row.get("source_seq", "") or row.get("seq_name", "")
                        if "|" not in seq_name:
                            continue
                        g = seq_name.split("|", 1)[0]
                        if g in batch_hosts:
                            genomad_proviruses_by_host[g].append(dict(row))
                elif member.name == "./batch_annotate/batch_genes.tsv":
                    for row in reader:
                        gene_id = row.get("gene", "")
                        if "|" not in gene_id:
                            continue
                        g, rest = gene_id.split("|", 1)
                        if g not in batch_hosts:
                            continue
                        raw_cid = rest.rsplit("_", 1)[0]
                        tname = row.get("taxname", "") or "unannotated"
                        genomad_mag_gene_tax_counts[g][tname] += 1
                        if raw_cid in focal_contig_raw_ids:
                            genomad_array_contig_gene_tax_counts[f"{g}|{raw_cid}"][tname] += 1
                            genomad_array_contig_genes[f"{g}|{raw_cid}"].append(dict(row))

    # 7. Scan archaeal_genes.tsv for gene counts and dark fraction on array contigs vs whole MAG
    arch_genes_path = ARCH_FN_ROOT / "inputs/archaeal_genes.tsv"
    record_source("archaeal_genes_tsv", arch_genes_path, compute_hash=False)
    mag_gene_catalogue_stats: dict[str, Counter] = defaultdict(Counter)
    contig_gene_catalogue_stats: dict[str, Counter] = defaultdict(Counter)
    with arch_genes_path.open("r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            g = row.get("genome", "")
            if g not in host_set:
                continue
            cid = row.get("contig", "")
            is_dark = float(row.get("dark", "0") or 0.0) > 0.0
            is_dark_strict = float(row.get("dark_strict", "0") or 0.0) > 0.0
            in_cat = int(float(row.get("contig_in_catalogue", "0") or 0.0))
            mag_gene_catalogue_stats[g]["n_genes"] += 1
            mag_gene_catalogue_stats[g]["n_dark"] += int(is_dark)
            mag_gene_catalogue_stats[g]["n_dark_strict"] += int(is_dark_strict)
            mag_gene_catalogue_stats[g]["n_in_catalogue"] += in_cat
            if cid in focal_contig_raw_ids:
                key = f"{g}|{cid}"
                contig_gene_catalogue_stats[key]["n_genes"] += 1
                contig_gene_catalogue_stats[key]["n_dark"] += int(is_dark)
                contig_gene_catalogue_stats[key]["n_dark_strict"] += int(is_dark_strict)
                contig_gene_catalogue_stats[key]["n_in_catalogue"] += in_cat

    # 8. Check if any per-sample competitive mapping outputs exist yet
    mapping_dir = EXP3_ROOT / "results/stage4_competitive_mapping/per_sample"
    completed_sample_coverage: dict[str, Any] = {}
    for host_mag, hmeta in FOCAL_HOSTS.items():
        samp = hmeta["source_sample"]
        cov_path = mapping_dir / f"{samp}_contig_coverage.tsv"
        done_path = mapping_dir / f"{samp}.done"
        if cov_path.exists() and done_path.exists():
            record_source(f"sample_coverage_{samp}", cov_path, compute_hash=True)
            rows = []
            with cov_path.open("r", encoding="utf-8") as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("entity_id") == host_mag or r.get("contig_id") in focal_contig_namespaced_ids:
                        rows.append(dict(r))
            completed_sample_coverage[samp] = rows

    # Assemble output JSON
    hosts_payload: dict[str, Any] = {}
    for host_mag, hmeta in FOCAL_HOSTS.items():
        finfo = fasta_data_by_host[host_mag]
        b = host_to_batch.get(host_mag, "")
        # Summarize geNomad contig taxonomy across all contigs in the MAG
        mag_contig_tax_counter = Counter()
        for cid in finfo["all_contig_lengths"]:
            trow = genomad_contig_taxonomy.get(f"{host_mag}|{cid}")
            lineage = trow["lineage"] if trow else "no_genomad_marker_taxonomy"
            mag_contig_tax_counter[lineage] += 1

        # Enrich array-bearing contigs with geNomad contig taxonomy, scores, proviruses, and gene markers
        enriched_contigs: dict[str, Any] = {}
        for cid, cdetail in finfo["array_bearing_contigs"].items():
            key = f"{host_mag}|{cid}"
            trow = genomad_contig_taxonomy.get(key)
            srow = genomad_contig_scores.get(key)
            genes = genomad_array_contig_genes.get(key, [])
            marker_genes = [
                {
                    "gene": g["gene"],
                    "start": int(g["start"]),
                    "end": int(g["end"]),
                    "strand": int(g["strand"]),
                    "marker": g["marker"],
                    "evalue": g["evalue"],
                    "bitscore": g["bitscore"],
                    "uscg": int(g["uscg"]) if g.get("uscg", "").isdigit() else 0,
                    "plasmid_hallmark": int(g["plasmid_hallmark"]) if g.get("plasmid_hallmark", "").isdigit() else 0,
                    "virus_hallmark": int(g["virus_hallmark"]) if g.get("virus_hallmark", "").isdigit() else 0,
                    "taxid": g["taxid"],
                    "taxname": g["taxname"],
                    "annotation_accessions": g["annotation_accessions"],
                    "annotation_description": g["annotation_description"],
                }
                for g in genes
                if g.get("marker") and g.get("marker") != "NA"
            ]
            # Cas operons on this contig
            ctg_cas = [c for c in cas_rows if c["mag_id"] == host_mag and c["raw_contig_id"] == cid]
            other_cas = [c for c in cas_rows if c["mag_id"] == host_mag and c["raw_contig_id"] != cid]
            enriched_contigs[cid] = {
                **cdetail,
                "cellular_contig_to_mag_memberships": sorted(set(contig_to_mags.get(cid, []))),
                "competitive_reference_ownership": reference_ownership_by_host[host_mag]["array_contig_ownership"].get(cid),
                "cas_operons_on_same_contig": ctg_cas,
                "cas_operons_on_other_mag_contigs_count": len(other_cas),
                "genomad_contig_taxonomy": trow,
                "genomad_aggregated_classification": srow,
                "genomad_contig_total_orfs": len(genes),
                "genomad_contig_marker_orfs_count": len(marker_genes),
                "genomad_contig_uscg_count": sum(1 for m in marker_genes if m["uscg"] == 1),
                "genomad_contig_gene_taxname_counts": dict(genomad_array_contig_gene_tax_counts.get(key, {})),
                "genomad_contig_marker_genes": marker_genes,
                "catalogue_gene_stats": dict(contig_gene_catalogue_stats.get(key, {})),
            }

        hosts_payload[host_mag] = {
            "host_mag": host_mag,
            "host_role": hmeta["host_role"],
            "source_sample": hmeta["source_sample"],
            "paired_targets": hmeta["paired_targets"],
            "eq_annot_batch_id": b,
            "mag_fasta_verification": {
                "fasta_sha256": finfo["fasta_sha256"],
                "n_contigs": finfo["n_contigs"],
                "total_bases": finfo["total_bases"],
                "n50_bp": finfo["n50_bp"],
                "max_contig_bp": finfo["max_contig_bp"],
                "min_contig_bp": finfo["min_contig_bp"],
                "stage0_record": s0_ver_rows.get(host_mag),
                "stage0_sample_duplicate_audit": s0_dup_rows.get(host_mag),
            },
            "mag_taxonomy_and_quality": {
                "refined_representatives_record": ref_rep_rows.get(host_mag),
                "taxonomy_reconciliation_record": tax_all_rows.get(host_mag),
                "gunc_per_genome_record": gunc_rows.get(host_mag),
                "gunc_tar_members": gunc_tar_members_by_host.get(host_mag, []),
                "gunc_all_levels": gunc_all_levels_by_host.get(host_mag, []),
                "gunc_contig_assignments_on_array_contigs": gunc_contig_assignments_by_host.get(host_mag, []),
                "genome_quality_record": qual_rows.get(host_mag),
            },
            "mag_competitive_reference_ownership_summary": {
                k: v for k, v in reference_ownership_by_host[host_mag].items() if k != "array_contig_ownership"
            },
            "mag_genomad_contig_taxonomy_summary": dict(mag_contig_tax_counter),
            "mag_genomad_gene_taxname_counts": dict(genomad_mag_gene_tax_counts.get(host_mag, {})),
            "mag_genomad_proviruses": genomad_proviruses_by_host.get(host_mag, []),
            "mag_catalogue_gene_stats": dict(mag_gene_catalogue_stats.get(host_mag, {})),
            "mag_defensefinder_systems_count": len(defensefinder_by_host.get(host_mag, [])),
            "mag_padloc_rows_count": len(padloc_by_host.get(host_mag, [])),
            "mag_cas_operons": [c for c in cas_rows if c["mag_id"] == host_mag],
            "array_bearing_contigs": enriched_contigs,
        }

    out_payload = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID", "unset"),
        "slurm_nodelist": os.environ.get("SLURM_NODELIST", "unset"),
        "extraction_provenance": {
            "script_path": "code/extract_stage4_host_array_contig_provenance.py",
            "script_sha256": stream_sha256(Path(__file__).resolve()),
            "source_files": source_hashes,
        },
        "hosts": hosts_payload,
        "pair_spacer_alignments": pair_alignments,
        "completed_sample_coverage": completed_sample_coverage,
    }

    out_path = EXP3_ROOT / "results/stage4_host_array_contig_provenance_extract.json"
    out_path.write_text(json.dumps(out_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"Wrote {out_path} ({out_path.stat().st_size} bytes, sha256={stream_sha256(out_path)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
