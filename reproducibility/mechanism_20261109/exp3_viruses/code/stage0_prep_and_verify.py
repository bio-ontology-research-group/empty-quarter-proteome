#!/usr/bin/env python3
"""Stage 0: One-Pass Refined Archive Extraction, Canonical 591-Archaeal-MAG Verification,
Sample/Duplicate Contig Overlap Audit, and Input Provenance Hashing.

Runs on a Shaheen compute node (`workq`, never on login node).
Performs ONE sequential pass over `/project/k10018/empty-quarter/data/mobilome/tars/refined_bins_all.tar`
to:
1. Extract all 591 canonical archaeal MAG FASTAs (`237` Thermoplasmatota + `354` other archaea)
   and verify each against `genome_manifest.tsv` (exact SHA-256, contig count, base count).
2. Compare against `/project/k10018/empty-quarter/data/catalog/mags/dereplicated/<mag>.fa.gz`
   under registered header transforms vs genuine contig-membership differences, retaining
   any unrefined/mismatched product as explicitly rejected (`rejected_wrong_product_contig_membership_differs`)
   while staging only the verified refined FASTA for downstream CRISPR/viral analysis.
3. Extract non-representative refined archaeal bins (mapped via `strict_mag_to_species.tsv`,
   `host_crosswalk.tsv`, or `refined_mag_manifest.tsv` `ar53` / `d__Archaea`) into a separate
   directory (`inputs/nonrep_arch_fastas/`) so their spacers can be tracked without conflating
   them with the 591 canonical representatives.
4. Build the complete cellular contig index (`inputs/cellular_contig_to_mag.tsv.gz`) across all
   refined bins in `refined_bins_all.tar` for Stage 2 cellular-chromosome decoy and provirus-boundary
   verification.
5. Audit sample-specific or duplicate contig/sequence overlap across all 591 archaeal MAGs
   (e.g. `self_SemiBin_*` vs `single_SemiBin_*` from the same sample sharing identical contigs).
6. Generate deterministic Stage 1 pilot (`inputs/pilot_30_mags.txt`) and 12 array chunk lists
   (`inputs/chunks/chunk_01.txt` .. `chunk_12.txt`) covering all 591 canonical archaeal MAGs.
"""
from __future__ import annotations

import collections
import csv
import gzip
import json
import os
import sys
import tarfile
import time
from pathlib import Path
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exp3_common import (
    parse_fasta_bytes,
    qualify_namespaced_contig_id,
    sha256_bytes,
    sha256_file,
    verify_and_reconcile_fasta,
)

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
IN_DIR = BRANCH_ROOT / "inputs"
RES_DIR = BRANCH_ROOT / "results"
PROV_DIR = BRANCH_ROOT / "provenance"

P_SCI_IN = Path("/scratch/project/k10018/eq-science-20261108/inputs")
P_ARCH_IN = Path("/scratch/project/k10018/eq-archaeal-function-20261106/inputs")
P_MOB_DATA = Path("/project/k10018/empty-quarter/data/mobilome")
P_DEREP_GZ = Path("/project/k10018/empty-quarter/data/catalog/mags/dereplicated")


def main() -> None:
    t0 = time.time()
    for d in (
        IN_DIR / "arch591_fastas",
        IN_DIR / "nonrep_arch_fastas",
        IN_DIR / "chunks",
        RES_DIR,
        PROV_DIR,
    ):
        d.mkdir(parents=True, exist_ok=True)

    # 1. Load canonical taxonomy (5,438 representatives -> 591 Archaea, 237 Thermoplasmatota)
    tax_path = IN_DIR / "taxonomy_reconciliation_all5438.tsv.gz"
    if not tax_path.exists():
        tax_path = Path(
            "/scratch/project/k10018/eq-science-20261109-canonical-taxonomy/results/taxonomy_reconciliation_all5438.tsv.gz"
        )
    tax_df = pd.read_csv(tax_path, sep="\t")
    if len(tax_df) != 5438:
        raise ValueError(f"Expected 5438 rows in taxonomy_reconciliation_all5438, got {len(tax_df)}")

    rep_5438_set = set(tax_df["genome"])
    arch_df = tax_df[tax_df["canonical_domain"] == "d__Archaea"].copy().sort_values("genome")
    if len(arch_df) != 591:
        raise ValueError(f"Expected 591 archaeal MAGs, got {len(arch_df)}")
    thermo_df = arch_df[arch_df["canonical_phylum"] == "p__Thermoplasmatota"]
    if len(thermo_df) != 237:
        raise ValueError(f"Expected 237 Thermoplasmatota MAGs, got {len(thermo_df)}")

    arch591_list = list(arch_df["genome"])
    arch591_set = set(arch591_list)
    tax_lookup = tax_df.set_index("genome").to_dict("index")

    # 2. Load authoritative genome_manifest.tsv
    gman_path = IN_DIR / "genome_manifest.tsv"
    if not gman_path.exists():
        gman_path = P_SCI_IN / "genome_manifest.tsv"
    gman_df = pd.read_csv(gman_path, sep="\t")
    gman_lookup = gman_df.set_index("mag").to_dict("index")

    # 3. Load GUNC, CheckM2 quality, and refined_mag_manifest (14,426 bins) + species crosswalks
    gunc_path = P_ARCH_IN / "c3_gunc_per_genome.tsv"
    gunc_lookup = {}
    if gunc_path.exists():
        gdf = pd.read_csv(gunc_path, sep="\t")
        gunc_lookup = gdf.set_index("genome").to_dict("index")

    nonrep_to_species: dict[str, str] = {}
    hcw_path = IN_DIR / "host_crosswalk.tsv"
    if hcw_path.exists():
        for r in csv.DictReader(open(hcw_path), delimiter="\t"):
            if r["host_mag"] not in rep_5438_set and r["species_rep"] in arch591_set:
                nonrep_to_species[r["host_mag"]] = r["species_rep"]

    skani_sp_path = P_MOB_DATA / "results_20261104/curated/e2/mag2sp/strict_mag_to_species.tsv"
    if skani_sp_path.exists():
        with open(skani_sp_path) as fh:
            first_line = fh.readline().strip().split("\t")
            fh.seek(0)
            if "species_rep" in first_line or "rep" in first_line:
                for r in csv.DictReader(fh, delimiter="\t"):
                    m_id = r.get("mag_id") or r.get("host_mag") or r.get("mag") or ""
                    s_id = r.get("species_rep") or r.get("rep") or ""
                    if m_id and s_id and m_id not in rep_5438_set and s_id in arch591_set:
                        nonrep_to_species.setdefault(m_id, s_id)
            else:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        m_id, s_id = parts[0], parts[1]
                        if m_id not in rep_5438_set and s_id in arch591_set:
                            nonrep_to_species.setdefault(m_id, s_id)

    ref_mag_man_path = P_MOB_DATA / "canonical_manifest_v1/refined_mag_manifest.tsv"
    refined_bin_meta: dict[str, dict[str, str]] = {}
    nonrep_arch_bins: set[str] = set(nonrep_to_species.keys())
    if ref_mag_man_path.exists():
        for r in csv.DictReader(open(ref_mag_man_man := ref_mag_man_path), delimiter="\t"):
            mid = r["mag_id"]
            refined_bin_meta[mid] = r
            if mid not in arch591_set:
                if r.get("gtdb_marker_set") == "ar53" or "d__Archaea" in r.get("gtdb_classification", ""):
                    nonrep_arch_bins.add(mid)

    print(
        f"[Stage 0] Canonical archaeal reps: {len(arch591_set)} (Thermoplasmatota: {len(thermo_df)}); "
        f"non-rep refined archaeal bins: {len(nonrep_arch_bins)}",
         flush=True,
    )

    # 4. Single-pass stream over refined_bins_all.tar
    tar_path = P_MOB_DATA / "tars/refined_bins_all.tar"
    if not tar_path.exists():
        raise FileNotFoundError(f"Missing refined_bins_all.tar at {tar_path}")

    verification_rows: list[dict[str, Any]] = []
    found_arch591: set[str] = set()
    found_nonrep_arch: set[str] = set()

    # Track contig -> MAGs and seq_sha256 -> MAGs across all 591 archaeal MAGs for duplicate/sample audit
    arch_contig_id_to_mags: dict[str, list[tuple[str, int, str]]] = collections.defaultdict(list)
    arch_seq_sha_to_mags: dict[str, list[tuple[str, str, int]]] = collections.defaultdict(list)

    contig_index_path = IN_DIR / "cellular_contig_to_mag.tsv.gz"
    n_tar_bins = 0
    n_tar_contigs = 0

    with gzip.open(contig_index_path, "wt") as gz_idx, tarfile.open(tar_path, "r:") as tf:
        gz_idx.write(
            "contig_id\tnamespaced_contig_id\tsource_sample\tmag_id\tis_5438_rep\tis_arch591_rep\tis_nonrep_arch\tspecies_rep\tdomain\tphylum\tcontig_len\tcontig_sha256\n"
        )
        for member in tf:
            if not member.isfile() or not member.name.endswith(".fa"):
                continue
            bname = Path(member.name).name
            mag_id = bname[:-3]
            source_sample = mag_id.split("__")[0] if "__" in mag_id else ""
            fobj = tf.extractfile(member)
            if fobj is None:
                continue
            raw_bytes = fobj.read()
            n_tar_bins += 1

            is_5438 = mag_id in rep_5438_set
            is_arch591 = mag_id in arch591_set
            is_nonrep_a = mag_id in nonrep_arch_bins
            sp_rep = mag_id if is_5438 else nonrep_to_species.get(mag_id, "")

            if is_5438:
                dom = str(tax_lookup[mag_id]["canonical_domain"])
                phy = str(tax_lookup[mag_id]["canonical_phylum"])
            else:
                rb = refined_bin_meta.get(mag_id, {})
                gtdb_str = rb.get("gtdb_classification", "")
                dom = "d__Archaea" if (is_nonrep_a or "d__Archaea" in gtdb_str) else "d__Bacteria_or_unplaced"
                phy = gtdb_str.split(";")[1] if ";" in gtdb_str else "unclassified"

            if is_arch591:
                found_arch591.add(mag_id)
                m_exp = gman_lookup[mag_id]

                # Strictly verify primary canonical archive bytes against genome_manifest FIRST
                # with verified_canonical_reference_bytes=None (never pass self as reference!)
                ver = verify_and_reconcile_fasta(
                    mag_id=mag_id,
                    raw_bytes=raw_bytes,
                    expected_sha256=str(m_exp["sha256"]),
                    expected_contigs=int(m_exp["contigs"]),
                    expected_bases=int(m_exp["bases"]),
                    verified_canonical_reference_bytes=None,
                )
                if (
                    not ver["is_verified_canonical"]
                    or ver["status"] != "verified_exact_match_genome_manifest"
                ):
                    raise RuntimeError(
                        f"Refined archive FASTA failed strict canonical manifest SHA-256 verification for {mag_id}: {ver}"
                    )

                out_fa = IN_DIR / "arch591_fastas" / f"{mag_id}.fa"
                out_fa.write_bytes(raw_bytes)
                ref_records, _ = parse_fasta_bytes(raw_bytes)

                # Only AFTER canonical raw_bytes passed strict manifest SHA-256 check may it
                # be passed as `verified_canonical_reference_bytes` to audit alternative products.
                alt_gz = P_DEREP_GZ / f"{mag_id}.fa.gz"
                alt_status = "missing"
                alt_sha256 = ""
                alt_n_contigs = -1
                alt_n_bases = -1
                alt_transform = "none"
                if alt_gz.exists():
                    alt_bytes = gzip.open(alt_gz, "rb").read()
                    alt_ver = verify_and_reconcile_fasta(
                        mag_id=mag_id,
                        raw_bytes=alt_bytes,
                        expected_sha256=str(m_exp["sha256"]),
                        expected_contigs=int(m_exp["contigs"]),
                        expected_bases=int(m_exp["bases"]),
                        verified_canonical_reference_bytes=raw_bytes,
                    )
                    alt_status = alt_ver["status"]
                    alt_sha256 = alt_ver["computed_sha256"]
                    alt_n_contigs = alt_ver["n_contigs"]
                    alt_n_bases = alt_ver["n_bases"]
                    alt_transform = alt_ver["header_transform_applied"]

                gunc_row = gunc_lookup.get(mag_id, {})
                t_row = tax_lookup[mag_id]
                verification_rows.append(
                    {
                        "genome": mag_id,
                        "source_sample": source_sample,
                        "binning_strategy": "self" if "__self_" in mag_id else ("single" if "__single_" in mag_id else "other"),
                        "has_sub_suffix": mag_id.endswith("_sub"),
                        "canonical_domain": t_row["canonical_domain"],
                        "canonical_phylum": t_row["canonical_phylum"],
                        "canonical_class": t_row["canonical_class"],
                        "canonical_order": t_row["canonical_order"],
                        "canonical_family": t_row["canonical_family"],
                        "canonical_genus": t_row["canonical_genus"],
                        "canonical_species": t_row["canonical_species"],
                        "checkm2_completeness": t_row["canonical_completeness"],
                        "checkm2_contamination": t_row["canonical_contamination"],
                        "pass_gunc": gunc_row.get("pass.GUNC", "unscored"),
                        "refined_tar_sha256": ver["computed_sha256"],
                        "manifest_expected_sha256": ver["expected_sha256"],
                        "refined_n_contigs": ver["n_contigs"],
                        "manifest_expected_contigs": ver["expected_contigs"],
                        "refined_n_bases": ver["n_bases"],
                        "manifest_expected_bases": ver["expected_bases"],
                        "refined_verification_status": ver["status"],
                        "catalog_derep_gz_sha256": alt_sha256,
                        "catalog_derep_gz_n_contigs": alt_n_contigs,
                        "catalog_derep_gz_n_bases": alt_n_bases,
                        "catalog_derep_gz_header_transform": alt_transform,
                        "catalog_derep_gz_disposition": alt_status,
                        "selected_authoritative_fasta": str(out_fa),
                    }
                )

                for h, seq in ref_records:
                    cid = h.split()[0]
                    _, _, q_cid, _, _ = qualify_namespaced_contig_id(cid, fallback_sample=source_sample)
                    c_sha = sha256_bytes(seq.encode("ascii"))
                    c_len = len(seq)
                    arch_contig_id_to_mags[q_cid].append((mag_id, c_len, c_sha))
                    arch_seq_sha_to_mags[c_sha].append((mag_id, q_cid, c_len))
                    gz_idx.write(
                        f"{cid}\t{q_cid}\t{source_sample}\t{mag_id}\t1\t1\t0\t{mag_id}\t{dom}\t{phy}\t{c_len}\t{c_sha}\n"
                    )
                    n_tar_contigs += 1

            elif is_nonrep_a:
                found_nonrep_arch.add(mag_id)
                out_fa = IN_DIR / "nonrep_arch_fastas" / f"{mag_id}.fa"
                out_fa.write_bytes(raw_bytes)
                records, _ = parse_fasta_bytes(raw_bytes)
                for h, seq in records:
                    cid = h.split()[0]
                    _, _, q_cid, _, _ = qualify_namespaced_contig_id(cid, fallback_sample=source_sample)
                    c_sha = sha256_bytes(seq.encode("ascii"))
                    c_len = len(seq)
                    gz_idx.write(
                        f"{cid}\t{q_cid}\t{source_sample}\t{mag_id}\t0\t0\t1\t{sp_rep}\t{dom}\t{phy}\t{c_len}\t{c_sha}\n"
                    )
                    n_tar_contigs += 1
            else:
                # Bacterial or other refined bin: fast header & length scan without hashing every sequence
                cur_cid = None
                cur_len = 0
                for line in raw_bytes.splitlines():
                    if line.startswith(b">"):
                        if cur_cid is not None:
                            _, _, q_cid, _, _ = qualify_namespaced_contig_id(cur_cid, fallback_sample=source_sample)
                            gz_idx.write(
                                f"{cur_cid}\t{q_cid}\t{source_sample}\t{mag_id}\t{1 if is_5438 else 0}\t0\t0\t{sp_rep}\t{dom}\t{phy}\t{cur_len}\tNA\n"
                            )
                            n_tar_contigs += 1
                        cur_cid = line[1:].split()[0].decode("utf-8", errors="replace")
                        cur_len = 0
                    else:
                        cur_len += len(line.strip())
                if cur_cid is not None:
                    _, _, q_cid, _, _ = qualify_namespaced_contig_id(cur_cid, fallback_sample=source_sample)
                    gz_idx.write(
                        f"{cur_cid}\t{q_cid}\t{source_sample}\t{mag_id}\t{1 if is_5438 else 0}\t0\t0\t{sp_rep}\t{dom}\t{phy}\t{cur_len}\tNA\n"
                    )
                    n_tar_contigs += 1

    if len(found_arch591) != 591:
        missing = sorted(arch591_set - found_arch591)
        raise RuntimeError(f"Missing {len(missing)} of 591 archaeal MAGs from refined_bins_all.tar: {missing[:10]}")

    # Write verification table
    ver_df = pd.DataFrame(verification_rows).sort_values("genome")
    ver_out = RES_DIR / "stage0_arch591_fasta_verification.tsv"
    ver_df.to_csv(ver_out, sep="\t", index=False)

    # 5. Audit sample-specific or duplicate MAG issues across the 591 archaeal MAGs
    # Check: (a) multiple archaeal MAGs from the same sample, (b) any shared contig IDs or identical contig sequences
    sample_to_mags: dict[str, list[str]] = collections.defaultdict(list)
    for g in arch591_list:
        sample_to_mags[g.split("__")[0]].append(g)

    pair_shared_contigs: dict[tuple[str, str], int] = collections.Counter()
    pair_shared_bp: dict[tuple[str, str], int] = collections.Counter()
    for c_sha, entries in arch_seq_sha_to_mags.items():
        distinct_mags = sorted({m for m, _, _ in entries})
        if len(distinct_mags) > 1:
            c_len = entries[0][2]
            for i in range(len(distinct_mags)):
                for j in range(i + 1, len(distinct_mags)):
                    pair_shared_contigs[(distinct_mags[i], distinct_mags[j])] += 1
                    pair_shared_bp[(distinct_mags[i], distinct_mags[j])] += c_len

    dup_audit_rows: list[dict[str, Any]] = []
    for g in arch591_list:
        samp = g.split("__")[0]
        co_mags = [m for m in sample_to_mags[samp] if m != g]
        same_genus_co = [
            m for m in co_mags
            if tax_lookup[m]["canonical_genus"] == tax_lookup[g]["canonical_genus"]
            and str(tax_lookup[g]["canonical_genus"]) not in ("g__", "")
        ]
        shared_partners = []
        shared_contig_count = 0
        shared_bp_count = 0
        for (m1, m2), cnt in pair_shared_contigs.items():
            if m1 == g or m2 == g:
                other = m2 if m1 == g else m1
                shared_partners.append(f"{other}({cnt})")
                shared_contig_count += cnt
                shared_bp_count += pair_shared_bp[(m1, m2)]
        dup_audit_rows.append(
            {
                "genome": g,
                "source_sample": samp,
                "binning_strategy": "self" if "__self_" in g else "single",
                "canonical_phylum": tax_lookup[g]["canonical_phylum"],
                "canonical_order": tax_lookup[g]["canonical_order"],
                "canonical_genus": tax_lookup[g]["canonical_genus"],
                "canonical_species": tax_lookup[g]["canonical_species"],
                "n_arch_mags_in_same_sample": len(sample_to_mags[samp]),
                "other_arch_mags_in_same_sample": ";".join(co_mags) if co_mags else "none",
                "same_sample_same_named_genus_mags": ";".join(same_genus_co) if same_genus_co else "none",
                "n_exact_shared_contig_seqs_with_other_arch591": shared_contig_count,
                "shared_contig_bp_with_other_arch591": shared_bp_count,
                "shared_contig_partner_mags": ";".join(shared_partners) if shared_partners else "none",
            }
        )
    dup_df = pd.DataFrame(dup_audit_rows)
    dup_out = RES_DIR / "stage0_arch591_sample_duplicate_contig_audit.tsv"
    dup_df.to_csv(dup_out, sep="\t", index=False)

    # 6. Build deterministic pilot (30 MAGs) and 12 array chunks (all 591 MAGs) + chunk 13 (nonrep archaeal bins)
    # Include the known archaeal MAGs with CRISPR hits from edges.tsv in the pilot so the pilot validates both positive and negative cases
    known_pos_candidates = [
        "41PRr2__single_SemiBin_236",
        "43PRr3__self_SemiBin_469",
        "S54PRr1__self_SemiBin_338",
        "T4PRr3__self_SemiBin_1415",
        "V45Dr1__single_SemiBin_136",
        "41PRr2__single_SemiBin_97",
        "48PRr2__self_SemiBin_737",
        "V55Dr1__self_SemiBin_1639",
        "41PRr2__self_SemiBin_3539",
        "S49Dr1__self_SemiBin_174",
    ]
    pilot_mags = [g for g in known_pos_candidates if g in arch591_set]
    for g in arch591_list:
        if len(pilot_mags) >= 30:
            break
        if g not in pilot_mags:
            pilot_mags.append(g)

    (IN_DIR / "pilot_30_mags.txt").write_text("\n".join(pilot_mags) + "\n")
    (IN_DIR / "arch591_mags.txt").write_text("\n".join(arch591_list) + "\n")

    n_chunks = 12
    for idx in range(n_chunks):
        chunk_mags = arch591_list[idx::n_chunks]
        (IN_DIR / "chunks" / f"chunk_{idx + 1:02d}.txt").write_text("\n".join(chunk_mags) + "\n")

    nonrep_sorted = sorted(found_nonrep_arch)
    (IN_DIR / "chunks" / "chunk_13_nonrep_arch.txt").write_text("\n".join(nonrep_sorted) + "\n")

    # 7. Record SHA-256 hashes of key input files in provenance/stage0_input_hashes.tsv
    files_to_hash = [
        tax_path,
        gman_path,
        gunc_path,
        ref_mag_man_path,
        skani_sp_path,
        P_MOB_DATA / "catalogues/votu_vclust_1.3.1_ani95_af85/votu_reps.fna",
        P_MOB_DATA / "catalogues/votu_vclust_1.3.1_ani95_af85/votu_membership.tsv",
        P_MOB_DATA / "results_20261104/curated/checkv_votu/quality_summary.tsv",
        P_MOB_DATA / "results_20261104/curated/checkv_votu/contamination.tsv",
        P_MOB_DATA / "results_20261104/curated/genomad_votu/votu_reps_annotate/votu_reps_genes.tsv",
        P_MOB_DATA / "results_20261104/curated/genomad_votu/votu_reps_annotate/votu_reps_taxonomy.tsv",
        P_MOB_DATA / "results_20261104/curated/decoy/decoys_derep.fna",
        P_MOB_DATA / "results_20261104/curated/e2/provirus_host_links.tsv",
        P_MOB_DATA / "results_20261104/curated/e5/repeats.tsv",
    ]
    hash_rows = []
    for p in files_to_hash:
        if p.exists():
            hash_rows.append(
                {
                    "path": str(p),
                    "size_bytes": p.stat().st_size,
                    "sha256": sha256_file(p),
                }
            )
    pd.DataFrame(hash_rows).to_csv(PROV_DIR / "stage0_input_hashes.tsv", sep="\t", index=False)

    summary = {
        "stage": "stage0_prep_and_verify",
        "status": "completed",
        "elapsed_seconds": round(time.time() - t0, 2),
        "n_canonical_5438_reps": len(rep_5438_set),
        "n_canonical_591_arch_mags": len(found_arch591),
        "n_canonical_237_thermo_mags": int(len(thermo_df)),
        "n_verified_exact_match_genome_manifest": int(
            (ver_df["refined_verification_status"] == "verified_exact_match_genome_manifest").sum()
        ),
        "n_catalog_derep_gz_rejected_wrong_product": int(
            (ver_df["catalog_derep_gz_disposition"] == "rejected_wrong_product_contig_membership_differs").sum()
        ),
        "n_catalog_derep_gz_verified_transform": int(
            ver_df["catalog_derep_gz_disposition"].isin(
                ["verified_exact_match_genome_manifest", "verified_identical_contigs_registered_header_transform"]
            ).sum()
        ),
        "n_nonrep_refined_arch_bins_extracted": len(found_nonrep_arch),
        "n_total_refined_bins_indexed": n_tar_bins,
        "n_total_cellular_contigs_indexed": n_tar_contigs,
        "n_arch591_with_shared_contigs": int((dup_df["n_exact_shared_contig_seqs_with_other_arch591"] > 0).sum()),
    }
    (RES_DIR / "stage0_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("[Stage 0] Summary:", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
