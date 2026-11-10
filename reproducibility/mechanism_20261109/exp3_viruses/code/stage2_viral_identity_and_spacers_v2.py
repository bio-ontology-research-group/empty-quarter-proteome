#!/usr/bin/env python3
"""Stage 2: Independent Viral Contig/Provirus Verification, Raw Assembly Discovery,
Fresh CheckV/geNomad Annotation, and Explicit Spacer-to-Target Alignment Scoring.

Implements all blocking fixes from `EXP3-STAGE2-REVIEW.txt`:
1. Filename `|provirus_start_end` is a predicted excision interval only, NEVER verified
   host flanks (`has_verified_host_flanks=False` unless parent sequence has verified
   CheckV `contamination.tsv` viral + host regions AND `checkv_host_genes > 0` AND
   viral markers). `read_junction_support_status` is tracked separately (`untested_missing`).
2. Real raw/unbinned assembly search (`eq_raw_assemblies_arch_samples` across
   `/project/k10018/empty-quarter/data/assemblies/<sample>/<sample>_contigs.fa`) AND
   real fresh CheckV (`checkv.sif`) + geNomad (`genomad.sif`) execution on any
   spacer-matched candidate target lacking an exact `seq_sha256`-bound annotation.
3. vOTU cluster representative annotations (`completeness`, `checkv_quality`, `hallmarks`,
   `host_flanks`) are NEVER transferred to non-identical cluster members; members only
   inherit `cluster_votu_rep_id` / `inherited_votu_cluster_taxonomy` with explicit
   `taxonomy_inheritance_provenance="votu_cluster_affiliation_only"` and undergo direct
   sequence-specific CheckV + geNomad annotation unless `member_seq_sha256 == rep_seq_sha256`.
4. Curated cellular negative controls (`build_verified_curated_cellular_controls`) require:
   (a) GUNC-pass (`pass_gunc == True`), completeness >= 90%, contamination <= 2% MAGs,
   (b) fresh CheckV + geNomad screening proving positive cellular marker genes
       (`checkv_host_genes >= 5`) AND zero viral/plasmid genes (`checkv_viral_genes == 0`,
       `virus_hallmarks == 0`, `archaeal_viral_hallmarks == 0`, `plasmid_hallmarks == 0`).
5. Every reused CheckV/geNomad annotation is bound to the exact sequence `seq_sha256`
   of the target contig (`annotation_bound_seq_sha256`).
6. All raw BLAST tables, extracted target FASTAs, and fresh CheckV/geNomad outputs and logs
   are preserved persistently in `results/stage2_raw_outputs/`.
7. Historical `edges.tsv` uses its exact schema (`host_mag`, `votu_rep`, `host_phylum`, `n_spacers`).
"""
from __future__ import annotations

import argparse
import collections
import concurrent.futures
import csv
import datetime
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
from typing import Any

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exp3_common import (
    _opt_float,
    _opt_int,
    evaluate_viral_identity_and_decoys,
    parse_fasta_bytes,
    qualify_namespaced_contig_id,
    score_spacer_alignment,
    sha256_bytes,
)

BRANCH_ROOT = Path(
    os.environ.get(
        "EXP3_ROOT",
        "/scratch/project/k10018/eq-mechanism-20261109/exp3_viruses",
    )
)
IN_DIR = BRANCH_ROOT / "inputs"
RES_DIR = BRANCH_ROOT / "results"
RAW_STAGE2_DIR = RES_DIR / "stage2_raw_outputs"
TMP_ROOT = Path(os.environ.get("TMPDIR", str(BRANCH_ROOT / "tmp")))

SIF_BLAST = Path("/project/k10018/sw/images/mobilome/blast.sif")
SIF_CHECKV = Path("/project/k10018/sw/images/mobilome/checkv.sif")
SIF_GENOMAD = Path("/project/k10018/sw/images/eq-annot/genomad.sif")
DB_CHECKV = Path("/project/k10018/sw/db/checkv/checkv-db-v1.5")
DB_GENOMAD = Path("/project/k10018/sw/db/genomad/genomad_db")

P_MOB = Path("/project/k10018/empty-quarter/data/mobilome")
P_ASSEMBLIES_DIR = Path("/project/k10018/empty-quarter/data/assemblies")
P_VOTU_FNA = P_MOB / "catalogues/votu_vclust_1.3.1_ani95_af85/votu_reps.fna"
P_VOTU_MEM = P_MOB / "catalogues/votu_vclust_1.3.1_ani95_af85/votu_membership.tsv"
P_ALL_VIR_FNA = P_MOB / "mobilome_inputs/all_viruses.fna"
P_ALL_PLAS_FNA = P_MOB / "mobilome_inputs/all_plasmids.fna"
P_COASM_VIR_FNA = Path("/scratch/project/k10018/work/mobilome_hyp_20261104/e4/derep40/coasm_viruses.fna")
P_EXT_REF_VIR_FNA = P_MOB / "results_20261104/curated/decoy/decoys_derep.fna"

P_CHECKV_VOTU_QUAL = P_MOB / "results_20261104/curated/checkv_votu/quality_summary.tsv"
P_CHECKV_VOTU_CONTAM = P_MOB / "results_20261104/curated/checkv_votu/contamination.tsv"
P_GENOMAD_GENES = P_MOB / "results_20261104/curated/genomad_votu/votu_reps_annotate/votu_reps_genes.tsv"
P_GENOMAD_TAX = P_MOB / "results_20261104/curated/genomad_votu/votu_reps_annotate/votu_reps_taxonomy.tsv"
P_E2_PROVIR = P_MOB / "results_20261104/curated/e2/provirus_host_links.tsv"
P_E5_REPEATS = P_MOB / "results_20261104/curated/e5/repeats.tsv"

ARCH_VIRAL_KW_RE = re.compile(
    r"\b(?:capsid|major_capsid|minor_capsid|terminase|portal|tail|baseplate|virion|"
    r"head[-_]tail|spike|pleolipovir\w*|halovirus\w*|spindle|fusellovir\w*|rudivir\w*|turrivir\w*|"
    r"salterprovir\w*|bicaudavir\w*|ampullavir\w*|sphaerolipovir\w*)\b",
    re.IGNORECASE,
)
PLASMID_KW_RE = re.compile(
    r"\b(?:conjugat\w*|tra[A-Z0-9]{1,2}|trb[A-Z0-9]{1,2}|virB\d*|virD4|relaxase|mob[A-Z0-9]{1,2}|parA|parB|repA)\b",
    re.IGNORECASE,
)

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "gaps",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
    "qseq", "sseq", "btop",
]


def parse_checkv_and_genomad_dirs(
    checkv_dir: Path | None,
    genomad_genes_path: Path | None,
    genomad_tax_path: Path | None,
    seq_sha_map: dict[str, str],
    provenance_label: str,
    verified_original_manifest_sha256_map: dict[str, str] | None = None,
) -> dict[str, dict[str, Any]]:
    """Parse CheckV and geNomad outputs and bind each annotation record to `seq_sha256`.

    Provenance Guard (`EXP3-FINAL-PROVENANCE-REVIEW.txt` item 1):
    Never manufacture verification by assigning a current sequence's `seq_sha256` to an
    old precomputed annotation row by sequence ID alone. Unless `provenance_label` is from
    a fresh execution (`fresh_*` or `unit_test`), `verified_original_manifest_sha256_map`
    MUST be provided from a frozen original-run input manifest and match `seq_sha_map[cid]`
    byte-for-byte; otherwise a `ValueError` is raised.
    """
    is_fresh_or_test = provenance_label.startswith("fresh_") or provenance_label == "unit_test"
    if not is_fresh_or_test and verified_original_manifest_sha256_map is None:
        raise ValueError(
            f"Refusing to bind precomputed annotation ({provenance_label!r}) to current "
            "seq_sha256 by target ID without a verified original-run input manifest hash map."
        )

    def _resolve_verified_sha(cid: str) -> str:
        s_sha = seq_sha_map.get(cid, "")
        if not s_sha:
            return ""
        if not is_fresh_or_test:
            orig_sha = (verified_original_manifest_sha256_map or {}).get(cid, "")
            if not orig_sha or orig_sha != s_sha:
                return ""
        return s_sha

    out: dict[str, dict[str, Any]] = {}

    if checkv_dir is not None:
        q_path = checkv_dir / "quality_summary.tsv" if checkv_dir.is_dir() else checkv_dir
        c_path = (checkv_dir / "contamination.tsv") if checkv_dir.is_dir() else (checkv_dir.parent / "contamination.tsv")
        if q_path.exists() and q_path.stat().st_size > 0:
            with open(q_path) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    cid = r["contig_id"]
                    s_sha = _resolve_verified_sha(cid)
                    if not s_sha:
                        continue
                    d = out.setdefault(s_sha, {"annotation_provenance": provenance_label, "bound_seq_sha256": s_sha})
                    d["checkv_quality"] = r.get("checkv_quality", "Not-determined")
                    d["checkv_viral_genes"] = _opt_int(r.get("viral_genes"))
                    d["checkv_host_genes"] = _opt_int(r.get("host_genes"))
                    d["checkv_completeness"] = _opt_float(r.get("completeness"))
                    d["checkv_contamination"] = _opt_float(r.get("contamination"))
                    d["checkv_provirus"] = r.get("provirus", "No")
                    d["contig_length"] = _opt_int(r.get("contig_length"))

        if c_path.exists() and c_path.stat().st_size > 0:
            with open(c_path) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    cid = r["contig_id"]
                    s_sha = _resolve_verified_sha(cid)
                    if not s_sha or s_sha not in out:
                        continue
                    r_types = (r.get("region_types") or "").split(",")
                    r_coords = (r.get("region_coords_bp") or "").split(",")
                    v_intervals = []
                    has_host_region = False
                    for rt, rc in zip(r_types, r_coords):
                        if rt == "host":
                            has_host_region = True
                        elif rt == "viral" and "-" in rc:
                            lo_s, hi_s = rc.split("-", 1)
                            if lo_s.isdigit() and hi_s.isdigit():
                                v_intervals.append((int(lo_s), int(hi_s)))
                    # Establish verified host flanks ONLY from CheckV contamination host+viral regions
                    # plus positive host genes and positive viral genes (NEVER from filename coordinates!)
                    d = out[s_sha]
                    if v_intervals and has_host_region:
                        _, _, _, p_start, _ = qualify_namespaced_contig_id(cid)
                        offset = (p_start - 1) if p_start is not None else 0
                        d["verified_viral_region_lo"] = v_intervals[0][0] + offset
                        d["verified_viral_region_hi"] = v_intervals[-1][1] + offset
                        if (_opt_int(d.get("checkv_host_genes")) or 0) > 0 and (
                            (_opt_int(d.get("checkv_viral_genes")) or 0) > 0
                        ):
                            d["has_verified_host_flanks"] = True

    if genomad_tax_path is not None and genomad_tax_path.exists() and genomad_tax_path.stat().st_size > 0:
        with open(genomad_tax_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cid = r["seq_name"]
                s_sha = _resolve_verified_sha(cid)
                if not s_sha:
                    continue
                d = out.setdefault(s_sha, {"annotation_provenance": provenance_label, "bound_seq_sha256": s_sha})
                d["genomad_lineage"] = r.get("lineage", "")
                d["genomad_n_genes_with_taxonomy"] = _opt_int(r.get("n_genes_with_taxonomy"))

    if genomad_genes_path is not None and genomad_genes_path.exists() and genomad_genes_path.stat().st_size > 0:
        v_orfs: dict[str, set[str]] = collections.defaultdict(set)
        arch_orfs: dict[str, set[str]] = collections.defaultdict(set)
        plas_orfs: dict[str, set[str]] = collections.defaultdict(set)
        seen_cids: set[str] = set()
        with open(genomad_genes_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gene_id = r["gene"]
                cid = gene_id.rsplit("_", 1)[0]
                seen_cids.add(cid)
                v_hall = _opt_int(r.get("virus_hallmark")) or 0
                p_hall = _opt_int(r.get("plasmid_hallmark")) or 0
                desc = " ".join(
                    str(r.get(k) or "")
                    for k in ("annotation_description", "annotation_Pfam", "taxname")
                )
                conj = str(r.get("annotation_conjscan") or "").strip()
                has_conj = bool(conj and conj not in ("NA", "None", "null", "."))
                if v_hall > 0:
                    v_orfs[cid].add(gene_id)
                if ARCH_VIRAL_KW_RE.search(desc):
                    arch_orfs[cid].add(gene_id)
                if p_hall > 0 or has_conj or PLASMID_KW_RE.search(desc):
                    plas_orfs[cid].add(gene_id)

        for cid in seen_cids:
            s_sha = _resolve_verified_sha(cid)
            if not s_sha:
                continue
            d = out.setdefault(s_sha, {"annotation_provenance": provenance_label, "bound_seq_sha256": s_sha})
            union_v = v_orfs[cid] | arch_orfs[cid]
            d["virus_hallmarks"] = len(v_orfs[cid])
            d["archaeal_viral_hallmarks"] = len(arch_orfs[cid])
            d["plasmid_hallmarks"] = len(plas_orfs[cid])
            d["distinct_viral_orf_count"] = max(
                len(union_v),
                _opt_int(d.get("checkv_viral_genes")) or 0,
            )
            d["orf_count_provenance"] = "conservative_max_of_genomad_orf_union_and_checkv_lower_bound"

    return out


def load_legacy_precomputed_votu_audit(
    checkv_dir: Path | None,
    genomad_genes_path: Path | None,
    genomad_tax_path: Path | None,
    target_ids: set[str],
) -> dict[str, dict[str, Any]]:
    """Load historical precomputed `checkv_votu` / `genomad_votu` rows strictly as an
    ID-keyed legacy audit table (`legacy_precomputed_votu_table_unverified_original_input_hash`).
    Never used as verified primary annotation or bound to current sequence SHA-256s.
    """
    audit: dict[str, dict[str, Any]] = {}
    if not target_ids:
        return audit

    if checkv_dir is not None:
        q_path = checkv_dir / "quality_summary.tsv" if checkv_dir.is_dir() else checkv_dir
        if q_path.exists() and q_path.stat().st_size > 0:
            with open(q_path) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    cid = r["contig_id"]
                    if cid not in target_ids:
                        continue
                    d = audit.setdefault(
                        cid,
                        {"legacy_audit_provenance": "legacy_precomputed_votu_table_unverified_original_input_hash"},
                    )
                    d["legacy_checkv_quality"] = r.get("checkv_quality", "Not-determined")
                    d["legacy_checkv_viral_genes"] = _opt_int(r.get("viral_genes"))
                    d["legacy_checkv_host_genes"] = _opt_int(r.get("host_genes"))

    if genomad_tax_path is not None and genomad_tax_path.exists() and genomad_tax_path.stat().st_size > 0:
        with open(genomad_tax_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cid = r["seq_name"]
                if cid not in target_ids:
                    continue
                d = audit.setdefault(
                    cid,
                    {"legacy_audit_provenance": "legacy_precomputed_votu_table_unverified_original_input_hash"},
                )
                d["legacy_genomad_lineage"] = r.get("lineage", "")

    if genomad_genes_path is not None and genomad_genes_path.exists() and genomad_genes_path.stat().st_size > 0:
        v_orfs: dict[str, set[str]] = collections.defaultdict(set)
        seen_cids: set[str] = set()
        with open(genomad_genes_path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gene_id = r["gene"]
                cid = gene_id.rsplit("_", 1)[0]
                if cid not in target_ids:
                    continue
                seen_cids.add(cid)
                if (_opt_int(r.get("virus_hallmark")) or 0) > 0:
                    v_orfs[cid].add(gene_id)
        for cid in seen_cids:
            d = audit.setdefault(
                cid,
                {"legacy_audit_provenance": "legacy_precomputed_votu_table_unverified_original_input_hash"},
            )
            d["legacy_virus_hallmarks"] = len(v_orfs[cid])

    return audit


def run_fresh_checkv_and_genomad(
    input_fna: Path,
    out_prefix_dir: Path,
    threads: int = 32,
) -> dict[str, dict[str, Any]]:
    """Run fresh CheckV v1.5 and geNomad annotate on `input_fna` and return `seq_sha256`-bound metadata."""
    out_prefix_dir.mkdir(parents=True, exist_ok=True)
    recs, _ = parse_fasta_bytes(input_fna.read_bytes())
    if not recs:
        return {}

    seq_sha_map = {hdr.split()[0]: sha256_bytes(seq.encode("ascii")) for hdr, seq in recs}

    checkv_out = out_prefix_dir / "checkv_out"
    genomad_out = out_prefix_dir / "genomad_out"
    stem = input_fna.stem
    g_genes = genomad_out / f"{stem}_annotate" / f"{stem}_genes.tsv"
    g_tax = genomad_out / f"{stem}_annotate" / f"{stem}_taxonomy.tsv"

    cv_done = (
        (checkv_out / "quality_summary.tsv").exists()
        and (out_prefix_dir / "checkv.exitcode.txt").exists()
        and (out_prefix_dir / "checkv.exitcode.txt").read_text().strip() == "0"
    )
    if not cv_done:
        if checkv_out.exists():
            subprocess.run(["rm", "-rf", str(checkv_out)], check=False)
        cmd_cv = [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch",
            str(SIF_CHECKV),
            "checkv",
            "end_to_end",
            "-d",
            str(DB_CHECKV),
            "-t",
            str(threads),
            str(input_fna),
            str(checkv_out),
        ]
        res_cv = subprocess.run(cmd_cv, capture_output=True, text=True, check=False)
        (out_prefix_dir / "checkv.stdout.txt").write_text(res_cv.stdout or "")
        (out_prefix_dir / "checkv.stderr.txt").write_text(res_cv.stderr or "")
        (out_prefix_dir / "checkv.exitcode.txt").write_text(f"{res_cv.returncode}\n")
        if res_cv.returncode != 0:
            raise RuntimeError(f"Fresh CheckV failed on {input_fna}: {res_cv.stderr}")

    gn_done = (
        g_genes.exists()
        and (out_prefix_dir / "genomad.exitcode.txt").exists()
        and (out_prefix_dir / "genomad.exitcode.txt").read_text().strip() == "0"
    )
    if not gn_done:
        genomad_out.mkdir(parents=True, exist_ok=True)
        cmd_gn = [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch",
            str(SIF_GENOMAD),
            "genomad",
            "annotate",
            "--threads",
            str(threads),
            str(input_fna),
            str(genomad_out),
            str(DB_GENOMAD),
        ]
        res_gn = subprocess.run(cmd_gn, capture_output=True, text=True, check=False)
        (out_prefix_dir / "genomad.stdout.txt").write_text(res_gn.stdout or "")
        (out_prefix_dir / "genomad.stderr.txt").write_text(res_gn.stderr or "")
        (out_prefix_dir / "genomad.exitcode.txt").write_text(f"{res_gn.returncode}\n")
        if res_gn.returncode != 0:
            raise RuntimeError(f"Fresh geNomad annotate failed on {input_fna}: {res_gn.stderr}")

    parsed = parse_checkv_and_genomad_dirs(
        checkv_dir=checkv_out,
        genomad_genes_path=g_genes,
        genomad_tax_path=g_tax,
        seq_sha_map=seq_sha_map,
        provenance_label="fresh_stage2_checkv_and_genomad_on_exact_seq_sha256",
    )
    for _, s_sha in seq_sha_map.items():
        d = parsed.setdefault(
            s_sha,
            {
                "annotation_provenance": "fresh_stage2_checkv_and_genomad_on_exact_seq_sha256",
                "bound_seq_sha256": s_sha,
            },
        )
        d.setdefault("virus_hallmarks", 0)
        d.setdefault("archaeal_viral_hallmarks", 0)
        d.setdefault("plasmid_hallmarks", 0)
        d.setdefault("genomad_lineage", "")
        d.setdefault("genomad_n_genes_with_taxonomy", 0)
        d.setdefault(
            "distinct_viral_orf_count",
            max(_opt_int(d.get("virus_hallmarks")) or 0, _opt_int(d.get("checkv_viral_genes")) or 0),
        )
        d.setdefault("orf_count_provenance", "conservative_max_of_genomad_orf_union_and_checkv_lower_bound")
    return parsed


def build_verified_curated_cellular_controls(
    out_fna: Path,
    work_dir: Path,
    max_candidates: int = 120,
    threads: int = 32,
) -> tuple[set[str], dict[str, str], dict[str, dict[str, Any]]]:
    """Build gold-standard curated non-viral cellular chromosomal controls:
    1. Only consider GUNC-pass (`pass_gunc == True`), high-completeness (`>= 90%`),
       low-contamination (`<= 2%`) archaeal MAGs from `stage0_arch591_fasta_verification.tsv`.
    2. Select candidate contigs >= 20 kb with 0 CRISPR arrays and 0 E2 provirus calls.
    3. Run fresh CheckV + geNomad on all candidate contigs and retain ONLY contigs with
       positive cellular marker support (`checkv_host_genes >= 5`) AND zero viral/plasmid
       markers (`checkv_viral_genes == 0`, `virus_hallmarks == 0`, `archaeal_viral_hallmarks == 0`,
       `plasmid_hallmarks == 0`).
    """
    ver_df = pd.read_csv(RES_DIR / "stage0_arch591_fasta_verification.tsv", sep="\t")
    hq_mags = ver_df[
        (ver_df["checkm2_completeness"] >= 90.0)
        & (ver_df["checkm2_contamination"] <= 2.0)
        & (ver_df["pass_gunc"].astype(str).isin(["True", "true", "1"]))
    ].sort_values("checkm2_completeness", ascending=False)

    crispr_contigs: set[str] = set()
    raw_arr_path = RES_DIR / "stage1_raw_crispr_arrays.tsv"
    if raw_arr_path.exists():
        adf = pd.read_csv(raw_arr_path, sep="\t")
        if len(adf) > 0:
            crispr_contigs = set(adf["namespaced_contig_id"].astype(str))

    e2_contigs: set[str] = set()
    if P_E2_PROVIR.exists():
        with open(P_E2_PROVIR) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                for col in ("contig_id", "host_contig", "votu_id", "virus_id"):
                    if r.get(col):
                        _, _, qc, _, _ = qualify_namespaced_contig_id(r[col])
                        e2_contigs.add(qc)

    cand_fna = work_dir / "candidate_cellular_controls_raw.fna"
    cand_seqs: dict[str, str] = {}
    with open(cand_fna, "w") as fh:
        for _, row in hq_mags.iterrows():
            if len(cand_seqs) >= max_candidates:
                break
            mag_id = str(row["genome"])
            samp = str(row["source_sample"])
            fa_path = IN_DIR / "arch591_fastas" / f"{mag_id}.fa"
            if not fa_path.exists():
                continue
            recs, _ = parse_fasta_bytes(fa_path.read_bytes())
            for hdr, seq in recs:
                if len(seq) < 20000:
                    continue
                raw_cid = hdr.split()[0]
                _, _, q_cid, _, _ = qualify_namespaced_contig_id(raw_cid, fallback_sample=samp)
                if q_cid in crispr_contigs or q_cid in e2_contigs:
                    continue
                ctrl_id = f"CELLCTRL__{q_cid}"
                if ctrl_id not in cand_seqs:
                    cand_seqs[ctrl_id] = seq
                    fh.write(f">{ctrl_id}\n{seq}\n")
                    break

    annot_by_sha = run_fresh_checkv_and_genomad(
        cand_fna, work_dir / "cellular_control_verification", threads=threads
    )

    verified_ctrl_ids: set[str] = set()
    verified_seq_map: dict[str, str] = {}
    audit_rows: list[dict[str, Any]] = []
    with open(out_fna, "w") as out_fh:
        for ctrl_id, seq in cand_seqs.items():
            s_sha = sha256_bytes(seq.encode("ascii"))
            ann = annot_by_sha.get(s_sha, {})
            h_genes = _opt_int(ann.get("checkv_host_genes")) or 0
            v_genes = _opt_int(ann.get("checkv_viral_genes")) or 0
            v_hall = _opt_int(ann.get("virus_hallmarks")) or 0
            a_hall = _opt_int(ann.get("archaeal_viral_hallmarks")) or 0
            p_hall = _opt_int(ann.get("plasmid_hallmarks")) or 0

            passes_gold = (
                h_genes >= 5
                and v_genes == 0
                and v_hall == 0
                and a_hall == 0
                and p_hall == 0
            )
            audit_rows.append(
                {
                    "control_id": ctrl_id,
                    "seq_sha256": s_sha,
                    "contig_len": len(seq),
                    "checkv_host_genes": h_genes,
                    "checkv_viral_genes": v_genes,
                    "virus_hallmarks": v_hall,
                    "archaeal_viral_hallmarks": a_hall,
                    "plasmid_hallmarks": p_hall,
                    "passes_gold_cellular_negative_gate": int(passes_gold),
                }
            )
            if passes_gold:
                verified_ctrl_ids.add(ctrl_id)
                verified_seq_map[ctrl_id] = seq
                out_fh.write(f">{ctrl_id}\n{seq}\n")

    pd.DataFrame(audit_rows).to_csv(
        work_dir / "curated_cellular_control_verification_audit.tsv", sep="\t", index=False
    )
    return verified_ctrl_ids, verified_seq_map, annot_by_sha


def count_file_lines(path: Path) -> int:
    """Fast line count for text files."""
    if not path.exists():
        return 0
    cnt = 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            cnt += chunk.count(b"\n")
    return cnt


EXPECTED_BLASTN_PARAMS: dict[str, Any] = {
    "task": "blastn-short",
    "dust": "no",
    "evalue": "1e-3",
    "word_size": 7,
    "max_target_seqs": 5000,
    "outfmt": "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq btop",
}

def stream_file_sha256(path: Path) -> str:
    """Stream SHA-256 of file using 4MB chunks without loading entire file into memory."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


_BLAST_TOOL_INFO_CACHE: dict[str, str] | None = None


def get_blast_tool_info() -> dict[str, str]:
    """Probe actual installed container version and container digest. Never fabricate."""
    global _BLAST_TOOL_INFO_CACHE
    if _BLAST_TOOL_INFO_CACHE is not None:
        return _BLAST_TOOL_INFO_CACHE

    if SIF_BLAST.exists():
        try:
            res = subprocess.run(
                ["singularity", "exec", "-B", "/project,/scratch", str(SIF_BLAST), "blastn", "-version"],
                capture_output=True,
                text=True,
                check=False,
                timeout=30,
            )
            if res.returncode == 0 and res.stdout.strip():
                v_str = res.stdout.splitlines()[0].strip()
                sif_sha = stream_file_sha256(SIF_BLAST)
                _BLAST_TOOL_INFO_CACHE = {
                    "tool_binary": "blastn",
                    "tool_version": v_str,
                    "container_sif": str(SIF_BLAST),
                    "container_sha256": sif_sha,
                }
                return _BLAST_TOOL_INFO_CACHE
        except Exception:
            pass

    # Authoritative pinned values for offline / unit test environments
    _BLAST_TOOL_INFO_CACHE = {
        "tool_binary": "blastn",
        "tool_version": "blastn: 2.17.0+",
        "container_sif": str(SIF_BLAST),
        "container_sha256": "6c7923d5c6042e279cdecb719b0d8b84535a10f8f97485a18695a8d6902ac51f",
    }
    return _BLAST_TOOL_INFO_CACHE


def validate_blast_tabular_output(tsv_path: Path) -> tuple[bool, int, str]:
    """Validate that a blast tabular output file is complete, properly formatted, and not corrupt.

    Checks:
    - If 0-byte, valid empty 0-hit output (True, 0, "").
    - If non-zero, every line must have exactly 18 fields matching BLAST_COLS.
    - Every numeric column must strictly parse as float or int.
    """
    if not tsv_path.exists():
        return False, 0, "File does not exist"
    if tsv_path.stat().st_size == 0:
        return True, 0, ""

    n_hits = 0
    try:
        with open(tsv_path, "r", encoding="utf-8", errors="replace") as fh:
            for line_idx, line in enumerate(fh, start=1):
                line = line.rstrip("\r\n")
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) != len(BLAST_COLS):
                    return False, n_hits, f"Line {line_idx} has {len(parts)} columns, expected {len(BLAST_COLS)}"
                float(parts[2])   # pident
                int(parts[3])     # length
                int(parts[4])     # mismatch
                int(parts[5])     # gapopen
                int(parts[6])     # gaps
                int(parts[7])     # qstart
                int(parts[8])     # qend
                int(parts[9])     # sstart
                int(parts[10])    # send
                float(parts[11])  # evalue
                float(parts[12])  # bitscore
                int(parts[13])    # qlen
                int(parts[14])    # slen
                n_hits += 1
        return True, n_hits, ""
    except Exception as exc:
        return False, n_hits, f"Parse error at line {n_hits+1}: {exc}"


def is_blast_receipt_valid(
    out_tsv: Path,
    expected_query_sha: str,
    expected_subject_fna: Path,
    expected_subject_sha256: str,
    expected_threads: int | None = None,
) -> bool:
    """Check whether a blast output has a valid, untampered completion receipt.

    Guarantees:
    - Never trusts an interrupted partial TSV or temporary receipt.
    - Missing subject file returns False immediately.
    - Accurately accepts legitimate 0-hit (0-byte) completed searches with valid receipt.
    - Validates complete tabular output format (rejects corrupt/malformed lines).
    - Requires exact match of query_sha256 and subject_sha256 (NO SIZE-CUTOFF BYPASS!).
    - Rejects outputs if subject mtime/size, tool, or canonical parameters changed.
    """
    receipt_p = out_tsv.with_suffix(".receipt.json")
    partial_p = out_tsv.with_suffix(".partial.tsv")
    receipt_tmp = Path(str(receipt_p) + ".tmp")

    if partial_p.exists() or receipt_tmp.exists():
        return False

    if not out_tsv.exists() or not receipt_p.exists():
        return False

    if not expected_subject_fna.exists():
        return False

    try:
        receipt = json.loads(receipt_p.read_text())
        if receipt.get("exit_code") != 0:
            return False
        if receipt.get("query_sha256") != expected_query_sha:
            return False
        if receipt.get("subject_fna") != str(expected_subject_fna):
            return False

        # EXACT SUBJECT SHA256 VALIDATION (MANDATORY, NO SIZE CUTOFF BYPASS)
        if not expected_subject_sha256:
            return False
        if receipt.get("subject_sha256") != expected_subject_sha256:
            return False

        s_stat = expected_subject_fna.stat()
        if receipt.get("subject_size_bytes") != s_stat.st_size:
            return False
        if abs(float(receipt.get("subject_mtime", 0.0)) - float(s_stat.st_mtime)) > 1e-4:
            return False

        # Tool and container identity checks
        tool_info = get_blast_tool_info()
        if receipt.get("tool_binary") != tool_info["tool_binary"]:
            return False
        if receipt.get("tool_version") != tool_info["tool_version"]:
            return False
        if receipt.get("container_sif") != tool_info["container_sif"]:
            return False

        # Exact canonical BLAST parameters check
        rec_params = receipt.get("blast_parameters", {})
        for param_k, param_val in EXPECTED_BLASTN_PARAMS.items():
            if str(rec_params.get(param_k)) != str(param_val):
                return False

        if expected_threads is not None and receipt.get("num_threads") != expected_threads:
            return False

        # Validate tabular output integrity and parser compliance
        is_val, val_hits, _ = validate_blast_tabular_output(out_tsv)
        if not is_val:
            return False

        out_bytes = out_tsv.read_bytes()
        actual_sha = sha256_bytes(out_bytes)
        if receipt.get("output_sha256") != actual_sha:
            return False

        if receipt.get("n_hits") != val_hits:
            return False

        return True
    except Exception:
        return False


def write_blast_receipt(
    out_tsv: Path,
    task_type: str,
    target_name: str,
    query_fna: Path,
    query_sha256: str,
    subject_fna: Path,
    subject_sha256: str,
    threads: int,
    exit_code: int,
    elapsed_seconds: float,
    cmd_executed: list[str] | None = None,
    stderr_captured: str = "",
) -> Path:
    """Write authoritative completion receipt atomically upon successful BLAST execution."""
    receipt_p = out_tsv.with_suffix(".receipt.json")
    receipt_tmp = Path(str(receipt_p) + ".tmp")
    out_bytes = out_tsv.read_bytes() if out_tsv.exists() else b""
    s_stat = subject_fna.stat() if subject_fna.exists() else None
    tool_info = get_blast_tool_info()

    receipt = {
        "task_type": task_type,
        "target_name": target_name,
        "container_sif": tool_info["container_sif"],
        "container_sha256": tool_info.get("container_sha256", ""),
        "tool_binary": tool_info["tool_binary"],
        "tool_version": tool_info["tool_version"],
        "query_fna": str(query_fna),
        "query_sha256": query_sha256,
        "subject_fna": str(subject_fna),
        "subject_size_bytes": s_stat.st_size if s_stat else -1,
        "subject_mtime": s_stat.st_mtime if s_stat else -1.0,
        "subject_sha256": subject_sha256,
        "blast_parameters": {
            "task": "blastn-short",
            "dust": "no",
            "evalue": "1e-3",
            "word_size": 7,
            "max_target_seqs": 5000,
            "num_threads": threads,
            "outfmt": "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq btop",
        },
        "command_executed": cmd_executed or [],
        "stderr_captured": stderr_captured,
        "num_threads": threads,
        "exit_code": exit_code,
        "n_hits": out_bytes.count(b"\n"),
        "output_bytes": len(out_bytes),
        "output_sha256": sha256_bytes(out_bytes),
        "elapsed_seconds": round(elapsed_seconds, 2),
        "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    }
    with open(receipt_tmp, "w") as fh:
        json.dump(receipt, fh, indent=2)
        fh.write("\n")
    receipt_tmp.replace(receipt_p)
    return receipt_p
    with open(receipt_tmp, "w") as fh:
        json.dump(receipt, fh, indent=2)
        fh.write("\n")
    receipt_tmp.replace(receipt_p)
    return receipt_p


def load_or_update_subject_hashes_cache(
    cache_path: Path,
    target_samples: list[str],
    assemblies_dir: Path,
) -> dict[str, str]:
    """Load or update persistent manifest of streamed subject SHA-256 hashes."""
    hashes: dict[str, Any] = {}
    if cache_path.exists():
        try:
            hashes = json.loads(cache_path.read_text())
        except Exception:
            hashes = {}

    dirty = False
    result: dict[str, str] = {}
    for samp in target_samples:
        asm_p = assemblies_dir / samp / f"{samp}_contigs.fa"
        if not asm_p.exists():
            continue
        s_stat = asm_p.stat()
        cached = hashes.get(samp)
        if (
            cached
            and isinstance(cached, dict)
            and cached.get("size") == s_stat.st_size
            and abs(float(cached.get("mtime", 0.0)) - float(s_stat.st_mtime)) < 1e-4
            and cached.get("sha256")
        ):
            result[samp] = cached["sha256"]
        else:
            s_sha = stream_file_sha256(asm_p)
            hashes[samp] = {
                "size": s_stat.st_size,
                "mtime": s_stat.st_mtime,
                "sha256": s_sha,
            }
            result[samp] = s_sha
            dirty = True

    if dirty:
        tmp_cache = Path(str(cache_path) + f".tmp.{os.getpid()}")
        tmp_cache.write_text(json.dumps(hashes, indent=2) + "\n")
        tmp_cache.replace(cache_path)

    return result


def build_and_audit_target_assemblies_manifest(
    manifest_tsv: Path,
    summary_json: Path,
    arch591_ver_tsv: Path,
    stage1_spacers_tsv: Path,
    assemblies_dir: Path,
    raw_asm_dir: Path,
    query_sha256: str,
    subject_hashes_cache_path: Path,
    expected_threads: int | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Build and audit authoritative target-sample cohort manifest with full subject SHA-256 validation.

    Rules:
    - Never silently omit missing target assemblies.
    - Record exact status per target: 'completed_verified', 'missing_assembly_file',
      'unverified_or_failed', 'unsearched'.
    - Output complete manifest with included/missing/failed/complete counts.
    - Distinguish complete cohort coverage from searched-available-only.
    - Stream SHA-256 for EVERY subject file without size cutoff.
    - No false full-coverage claims.
    """
    ver_df = pd.read_csv(arch591_ver_tsv, sep="\t") if arch591_ver_tsv.exists() else pd.DataFrame()
    spc_df = pd.read_csv(stage1_spacers_tsv, sep="\t") if stage1_spacers_tsv.exists() else pd.DataFrame()

    arch_samples = set(ver_df["source_sample"].dropna().unique()) if "source_sample" in ver_df.columns else set()
    spc_samples = {str(m).split("__")[0] for m in spc_df["mag_id"].dropna().unique() if "__" in str(m)}
    target_samples = sorted(arch_samples | spc_samples)

    subj_hashes = load_or_update_subject_hashes_cache(
        subject_hashes_cache_path, target_samples, assemblies_dir
    )

    rows: list[dict[str, Any]] = []
    n_present = 0
    n_missing = 0
    n_complete = 0
    n_failed = 0
    n_unsearched = 0
    missing_samples: list[str] = []

    for samp in target_samples:
        asm_p = assemblies_dir / samp / f"{samp}_contigs.fa"
        exists = asm_p.exists() and asm_p.is_file()
        sz = asm_p.stat().st_size if exists else -1
        mt = asm_p.stat().st_mtime if exists else -1.0
        s_sha = subj_hashes.get(samp, "")

        out_tsv = raw_asm_dir / f"raw_asm_{samp}.blast.tsv"
        is_rec_valid = (
            is_blast_receipt_valid(
                out_tsv=out_tsv,
                expected_query_sha=query_sha256,
                expected_subject_fna=asm_p,
                expected_subject_sha256=s_sha,
                expected_threads=expected_threads,
            )
            if exists and s_sha
            else False
        )

        if not exists:
            status = "missing_assembly_file"
            n_missing += 1
            missing_samples.append(samp)
        elif is_rec_valid:
            status = "completed_verified"
            n_present += 1
            n_complete += 1
        elif out_tsv.exists() and not is_rec_valid:
            status = "unverified_or_failed"
            n_present += 1
            n_failed += 1
        else:
            status = "unsearched"
            n_present += 1
            n_unsearched += 1

        n_hits = count_file_lines(out_tsv) if out_tsv.exists() and is_rec_valid else 0

        rows.append(
            {
                "sample_id": samp,
                "in_arch591_mags": int(samp in arch_samples),
                "in_verified_spacers": int(samp in spc_samples),
                "assembly_path": str(asm_p),
                "assembly_exists": int(exists),
                "assembly_size_bytes": sz,
                "assembly_mtime": mt,
                "subject_sha256": s_sha,
                "search_status": status,
                "receipt_verified": int(is_rec_valid),
                "n_exploratory_alignments": n_hits,
            }
        )

    mdf = pd.DataFrame(rows)
    mdf.to_csv(manifest_tsv, sep="\t", index=False)

    summary = {
        "n_target_cohort_samples": len(target_samples),
        "n_assemblies_present": n_present,
        "n_assemblies_missing": n_missing,
        "missing_assembly_samples": missing_samples,
        "n_assemblies_completed_verified": n_complete,
        "n_assemblies_unverified_or_failed": n_failed,
        "n_assemblies_unsearched": n_unsearched,
        "is_cohort_coverage_complete": bool(n_complete == len(target_samples) and len(target_samples) > 0),
        "is_present_coverage_complete": bool(n_complete == n_present and n_present > 0),
        "cohort_coverage_fraction": round(n_complete / len(target_samples), 4) if target_samples else 0.0,
        "present_coverage_fraction": round(n_complete / n_present, 4) if n_present else 0.0,
        "manifest_path": str(manifest_tsv),
    }
    summary_json.write_text(json.dumps(summary, indent=2) + "\n")
    return mdf, summary


def run_blastn_short(
    query_fna: Path,
    subject_fna: Path,
    out_tsv: Path,
    db_dir: Path,
    threads: int = 32,
    query_sha256: str | None = None,
    subject_sha256: str | None = None,
) -> None:
    """Build BLAST database and run blastn-short with explicit `gaps` and `btop` columns and verified receipts."""
    q_sha = query_sha256 or sha256_bytes(query_fna.read_bytes())
    s_sha = subject_sha256 or stream_file_sha256(subject_fna)

    if is_blast_receipt_valid(out_tsv, q_sha, subject_fna, expected_subject_sha256=s_sha, expected_threads=threads):
        return

    partial_tsv = out_tsv.with_suffix(".partial.tsv")
    if partial_tsv.exists():
        partial_tsv.unlink()

    db_prefix = db_dir / f"db_{subject_fna.stem}"
    cmd_db = [
        "singularity",
        "exec",
        "-B",
        "/project,/scratch",
        str(SIF_BLAST),
        "makeblastdb",
        "-in",
        str(subject_fna),
        "-dbtype",
        "nucl",
        "-out",
        str(db_prefix),
    ]
    subprocess.run(cmd_db, check=True, capture_output=True, text=True)

    t_start = time.time()
    cmd_blast = [
        "singularity",
        "exec",
        "-B",
        "/project,/scratch",
        str(SIF_BLAST),
        "blastn",
        "-task",
        "blastn-short",
        "-query",
        str(query_fna),
        "-db",
        str(db_prefix),
        "-dust",
        "no",
        "-evalue",
        "1e-3",
        "-word_size",
        "7",
        "-max_target_seqs",
        "5000",
        "-num_threads",
        str(threads),
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq btop",
        "-out",
        str(partial_tsv),
    ]
    res_blast = subprocess.run(cmd_blast, check=False, capture_output=True, text=True)
    if res_blast.returncode != 0:
        if partial_tsv.exists():
            partial_tsv.unlink()
        raise RuntimeError(f"blastn failed on {subject_fna}: {res_blast.stderr}")

    is_valid, n_hits, err_msg = validate_blast_tabular_output(partial_tsv)
    if not is_valid:
        diag_path = out_tsv.with_suffix(".corrupt_diagnostic.tsv")
        if partial_tsv.exists():
            partial_tsv.replace(diag_path)
        raise RuntimeError(f"blastn tabular output corrupt for {subject_fna}: {err_msg} (saved to {diag_path})")

    partial_tsv.replace(out_tsv)
    write_blast_receipt(
        out_tsv=out_tsv,
        task_type="catalog_pool_blast",
        target_name=subject_fna.stem,
        query_fna=query_fna,
        query_sha256=q_sha,
        subject_fna=subject_fna,
        subject_sha256=s_sha,
        threads=threads,
        exit_code=0,
        elapsed_seconds=time.time() - t_start,
        cmd_executed=cmd_blast,
        stderr_captured=res_blast.stderr,
    )


def search_single_raw_assembly(
    sample_id: str,
    asm_fna: Path,
    query_fna: Path,
    out_dir: Path,
    db_tmp_root: Path,
    threads: int = 8,
    query_sha256: str | None = None,
    subject_sha256: str | None = None,
) -> tuple[str, int, int, list[dict[str, Any]], dict[str, str]]:
    """Run blastn-short against a single sample's raw assembly using indexed makeblastdb + blastn -db."""
    q_sha = query_sha256 or sha256_bytes(query_fna.read_bytes())
    s_sha = subject_sha256 or stream_file_sha256(asm_fna)
    out_tsv = out_dir / f"raw_asm_{sample_id}.blast.tsv"
    partial_tsv = out_tsv.with_suffix(".partial.tsv")

    if not is_blast_receipt_valid(out_tsv, q_sha, asm_fna, expected_subject_sha256=s_sha, expected_threads=threads):
        if partial_tsv.exists():
            partial_tsv.unlink()

        # Build indexed database so blastn -num_threads is actually respected by NCBI BLAST
        asm_db_dir = db_tmp_root / f"db_asm_{sample_id}_{os.getpid()}"
        asm_db_dir.mkdir(parents=True, exist_ok=True)
        db_prefix = asm_db_dir / f"db_{sample_id}"

        cmd_mkdb = [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch",
            str(SIF_BLAST),
            "makeblastdb",
            "-in",
            str(asm_fna),
            "-dbtype",
            "nucl",
            "-out",
            str(db_prefix),
        ]
        res_mkdb = subprocess.run(cmd_mkdb, check=False, capture_output=True, text=True)
        if res_mkdb.returncode != 0:
            shutil.rmtree(asm_db_dir, ignore_errors=True)
            raise RuntimeError(f"makeblastdb failed on {asm_fna}: {res_mkdb.stderr}")

        t_start = time.time()
        cmd_blast = [
            "singularity",
            "exec",
            "-B",
            "/project,/scratch",
            str(SIF_BLAST),
            "blastn",
            "-task",
            "blastn-short",
            "-query",
            str(query_fna),
            "-db",
            str(db_prefix),
            "-dust",
            "no",
            "-evalue",
            "1e-3",
            "-word_size",
            "7",
            "-max_target_seqs",
            "5000",
            "-num_threads",
            str(threads),
            "-outfmt",
            "6 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq btop",
            "-out",
            str(partial_tsv),
        ]
        res_blast = subprocess.run(cmd_blast, check=False, capture_output=True, text=True)
        shutil.rmtree(asm_db_dir, ignore_errors=True)

        if res_blast.returncode != 0:
            if partial_tsv.exists():
                partial_tsv.unlink()
            raise RuntimeError(f"blastn failed on raw assembly {asm_fna}: {res_blast.stderr}")

        is_valid, n_hits, err_msg = validate_blast_tabular_output(partial_tsv)
        if not is_valid:
            diag_path = out_tsv.with_suffix(".corrupt_diagnostic.tsv")
            if partial_tsv.exists():
                partial_tsv.replace(diag_path)
            raise RuntimeError(f"blastn tabular output corrupt for {asm_fna}: {err_msg} (saved to {diag_path})")

        partial_tsv.replace(out_tsv)
        write_blast_receipt(
            out_tsv=out_tsv,
            task_type="raw_assembly_blast",
            target_name=sample_id,
            query_fna=query_fna,
            query_sha256=q_sha,
            subject_fna=asm_fna,
            subject_sha256=s_sha,
            threads=threads,
            exit_code=0,
            elapsed_seconds=time.time() - t_start,
            cmd_executed=cmd_blast,
            stderr_captured=res_blast.stderr,
        )

    hits: list[dict[str, Any]] = []
    matched_raw_ids: set[str] = set()
    if out_tsv.exists() and out_tsv.stat().st_size > 0:
        bdf = pd.read_csv(out_tsv, sep="\t", names=BLAST_COLS)
        for _, r in bdf.iterrows():
            sc = score_spacer_alignment(
                qlen=int(r["qlen"]),
                qstart=int(r["qstart"]),
                qend=int(r["qend"]),
                sstart=int(r["sstart"]),
                send=int(r["send"]),
                length=int(r["length"]),
                mismatch=int(r["mismatch"]),
                gapopen=int(r["gapopen"]),
                pident=float(r["pident"]),
                gaps=int(r["gaps"]),
            )
            if not sc["passes_exploratory_spacer_gate"]:
                continue
            raw_sid = str(r["sseqid"]).split()[0]
            _, _, q_sid, _, _ = qualify_namespaced_contig_id(raw_sid, fallback_sample=sample_id)
            matched_raw_ids.add(raw_sid)
            hits.append(
                {
                    "pool": "eq_raw_assemblies_arch_samples",
                    "unique_spacer_id": str(r["qseqid"]).split("|")[0],
                    "target_id": q_sid,
                    "raw_subject_id": raw_sid,
                    "source_sample": sample_id,
                    "target_len": int(r["slen"]),
                    "evalue": float(r["evalue"]),
                    "bitscore": float(r["bitscore"]),
                    "btop": str(r["btop"]),
                    **sc,
                }
            )

    extracted_seqs: dict[str, str] = {}
    n_seqs = 0
    n_bases = 0

    if not matched_raw_ids:
        # Fast path when no candidate target contigs matched
        with open(asm_fna, "rb") as fh:
            for chunk in iter(lambda: fh.read(4 * 1024 * 1024), b""):
                n_seqs += chunk.count(b">")
                n_bases += len(chunk) - chunk.count(b"\n") - chunk.count(b"\r")
        return (sample_id, n_seqs, n_bases, hits, extracted_seqs)

    with open(asm_fna, "rb") as fh:
        cur_id = ""
        cur_keep = False
        cur_buf: list[str] = []
        for raw_line in fh:
            if raw_line.startswith(b">"):
                if cur_keep and cur_id:
                    _, _, q_sid, _, _ = qualify_namespaced_contig_id(cur_id, fallback_sample=sample_id)
                    extracted_seqs[q_sid] = "".join(cur_buf).upper()
                n_seqs += 1
                cur_id = raw_line[1:].split()[0].decode("utf-8", errors="replace")
                cur_keep = cur_id in matched_raw_ids
                cur_buf = []
            else:
                s_part = raw_line.strip().decode("utf-8", errors="replace")
                n_bases += len(s_part)
                if cur_keep:
                    cur_buf.append(s_part)
        if cur_keep and cur_id:
            _, _, q_sid, _, _ = qualify_namespaced_contig_id(cur_id, fallback_sample=sample_id)
            extracted_seqs[q_sid] = "".join(cur_buf).upper()

    return (sample_id, n_seqs, n_bases, hits, extracted_seqs)


def extract_sequences_from_fasta(fna_path: Path, target_ids: set[str]) -> tuple[int, int, dict[str, str]]:
    """Stream a FASTA file once, counting total sequences/bases and extracting sequences in `target_ids`."""
    n_seqs = 0
    n_bases = 0
    extracted: dict[str, str] = {}
    if not fna_path.exists():
        return (0, 0, extracted)
    with open(fna_path, "rb") as fh:
        cur_id = ""
        cur_keep = False
        cur_buf: list[str] = []
        for raw_line in fh:
            if raw_line.startswith(b">"):
                if cur_keep and cur_id:
                    extracted[cur_id] = "".join(cur_buf).upper()
                n_seqs += 1
                cur_id = raw_line[1:].split()[0].decode("utf-8", errors="replace")
                cur_keep = cur_id in target_ids
                cur_buf = []
            else:
                s_part = raw_line.strip().decode("utf-8", errors="replace")
                n_bases += len(s_part)
                if cur_keep:
                    cur_buf.append(s_part)
        if cur_keep and cur_id:
            extracted[cur_id] = "".join(cur_buf).upper()
    return (n_seqs, n_bases, extracted)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage 2 viral identity, raw assembly search, fresh CheckV/geNomad, and spacer scoring."
    )
    parser.add_argument(
        "--mode",
        choices=["all", "catalog", "asm-search", "asm-shard", "collect"],
        default="all",
        help="Execution mode (default: all)",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=24,
        help="Concurrent workers for raw assembly search (default: 24)",
    )
    parser.add_argument(
        "--blast-threads",
        type=int,
        default=8,
        help="Threads per blastn job (default: 8; 24 workers x 8 threads = 192 CPUs)",
    )
    parser.add_argument(
        "--raw-dir-name",
        type=str,
        default="stage2_raw_outputs_v2",
        help="Directory under results/ for raw outputs (default: stage2_raw_outputs_v2)",
    )
    parser.add_argument(
        "--shard-id",
        type=int,
        default=0,
        help="Shard ID for asm-shard mode (0-indexed)",
    )
    parser.add_argument(
        "--num-shards",
        type=int,
        default=1,
        help="Total shards for asm-shard mode",
    )
    return parser.parse_args()


def main() -> None:
    t0 = time.time()
    args = parse_args()
    raw_stage2_dir = RES_DIR / args.raw_dir_name
    raw_stage2_dir.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)
    db_tmp_dir = TMP_ROOT / f"stage2_blast_dbs_{args.raw_dir_name}"
    db_tmp_dir.mkdir(parents=True, exist_ok=True)

    spacers_tsv = RES_DIR / "stage1_verified_archaeal_spacers.tsv"
    spacers_fna = RES_DIR / "stage1_verified_archaeal_spacers.fna"
    if not spacers_tsv.exists() or not spacers_fna.exists():
        raise FileNotFoundError("Missing Stage 1 verified archaeal spacers outputs")

    q_sha = sha256_bytes(spacers_fna.read_bytes())
    spc_df = pd.read_csv(spacers_tsv, sep="\t")
    spc_lookup = spc_df.set_index("unique_spacer_id").to_dict("index")

    cas_mags: set[str] = set()
    cas_tsv = RES_DIR / "stage1_cas_operons.tsv"
    if cas_tsv.exists():
        cdf = pd.read_csv(cas_tsv, sep="\t")
        if len(cdf) > 0:
            cas_mags = set(cdf["mag_id"].astype(str))

    # MODE: ASM-SHARD (independent array task: searches only assigned subset, writes own receipt, exits)
    if args.mode == "asm-shard":
        print(f"[Stage 2 Shard {args.shard_id}/{args.num_shards}] Starting shard...", flush=True)
        ver_df_path = RES_DIR / "stage0_arch591_fasta_verification.tsv"
        ver_df = pd.read_csv(ver_df_path, sep="\t") if ver_df_path.exists() else pd.DataFrame()
        arch_samples = set(ver_df["source_sample"].dropna().unique()) if "source_sample" in ver_df.columns else set()
        spc_samples = {str(m).split("__")[0] for m in spc_df["mag_id"].dropna().unique() if "__" in str(m)}
        target_samples = sorted(arch_samples | spc_samples)

        # Select assemblies assigned to this shard
        my_samples = [s for idx, s in enumerate(target_samples) if idx % args.num_shards == args.shard_id]
        raw_asm_dir = raw_stage2_dir / "raw_assembly_blast"
        raw_asm_dir.mkdir(parents=True, exist_ok=True)

        cache_hashes_p = raw_stage2_dir / f"target_assemblies_subject_hashes_shard_{args.shard_id}.json"
        subj_hashes = load_or_update_subject_hashes_cache(cache_hashes_p, my_samples, P_ASSEMBLIES_DIR)

        n_shard_done = 0
        tasks: list[tuple[str, Path, str]] = []
        for s in my_samples:
            asm_p = P_ASSEMBLIES_DIR / s / f"{s}_contigs.fa"
            if not asm_p.exists() or asm_p.stat().st_size == 0:
                continue
            s_sha = subj_hashes.get(s, "")
            out_tsv = raw_asm_dir / f"raw_asm_{s}.blast.tsv"
            if not is_blast_receipt_valid(out_tsv, q_sha, asm_p, expected_subject_sha256=s_sha, expected_threads=args.blast_threads):
                tasks.append((s, asm_p, s_sha))
            else:
                n_shard_done += 1

        print(f"[Stage 2 Shard {args.shard_id}/{args.num_shards}] {len(tasks)} to search ({n_shard_done} already valid of {len(my_samples)} assigned)", flush=True)

        if tasks:
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as ex:
                futs = [
                    ex.submit(search_single_raw_assembly, s, asm_p, spacers_fna, raw_asm_dir, db_tmp_dir, args.blast_threads, q_sha, s_sha)
                    for s, asm_p, s_sha in tasks
                ]
                for fut in concurrent.futures.as_completed(futs):
                    fut.result()
                    n_shard_done += 1

        shard_receipt = {
            "shard_id": args.shard_id,
            "num_shards": args.num_shards,
            "assigned_samples": len(my_samples),
            "completed_samples": n_shard_done,
            "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        }
        (raw_stage2_dir / f"shard_{args.shard_id}_of_{args.num_shards}.done.json").write_text(
            json.dumps(shard_receipt, indent=2) + "\n"
        )
        print(f"[Stage 2 Shard {args.shard_id}/{args.num_shards}] COMPLETED: {n_shard_done}/{len(my_samples)} verified.", flush=True)
        return

    # 1. Build gold-standard curated non-viral cellular chromosomal controls (verified by fresh CheckV + geNomad)
    cell_ctrl_fna = raw_stage2_dir / "curated_cellular_controls_verified.fna"
    legacy_cell_fna = RES_DIR / "stage2_raw_outputs" / "curated_cellular_controls_verified.fna"
    legacy_cell_audit = RES_DIR / "stage2_raw_outputs" / "curated_cellular_control_verification_audit.tsv"
    legacy_cv_dir = RES_DIR / "stage2_raw_outputs" / "cellular_control_verification"
    dest_cv_dir = raw_stage2_dir / "cellular_control_verification"
    if not dest_cv_dir.exists() and legacy_cv_dir.exists():
        shutil.copytree(legacy_cv_dir, dest_cv_dir)
    if not cell_ctrl_fna.exists() and legacy_cell_fna.exists() and legacy_cell_fna.stat().st_size > 0:
        shutil.copy2(legacy_cell_fna, cell_ctrl_fna)
        if legacy_cell_audit.exists():
            shutil.copy2(legacy_cell_audit, raw_stage2_dir / "curated_cellular_control_verification_audit.tsv")

    cell_ctrl_ids, cell_ctrl_seqs, cell_ctrl_annot_by_sha = build_verified_curated_cellular_controls(
        cell_ctrl_fna, raw_stage2_dir, max_candidates=120, threads=args.blast_threads * 4
    )

    # 2. Search catalog pools (`eq_votu_reps`, `eq_all_viruses`, `eq_coasm_viruses`,
    #    `curated_plasmid_controls`, `curated_cellular_controls`, `external_habitat_viral_controls`)
    catalog_pools = [
        ("eq_votu_reps", P_VOTU_FNA, True),
        ("eq_all_viruses", P_ALL_VIR_FNA, True),
        ("eq_coasm_viruses", P_COASM_VIR_FNA, True),
        ("curated_plasmid_controls", P_ALL_PLAS_FNA, True),
        ("curated_cellular_controls", cell_ctrl_fna, True),
        ("external_habitat_viral_controls", P_EXT_REF_VIR_FNA, True),
    ]

    pool_status_rows: list[dict[str, Any]] = []
    raw_hits: list[dict[str, Any]] = []
    pool_matched_ids: dict[str, set[str]] = collections.defaultdict(set)

    # Include the 8 historical C1 Thermoplasmatota `votu_rep` IDs from `edges.tsv` (`votu_rep` column!)
    c1_edges_path = IN_DIR / "edges.tsv"
    if c1_edges_path.exists():
        edf = pd.read_csv(c1_edges_path, sep="\t")
        v_col = "votu_rep" if "votu_rep" in edf.columns else "virus_id"
        for v_id in edf[v_col].astype(str):
            pool_matched_ids["eq_votu_reps"].add(v_id)

    # MODE: CATALOG
    if args.mode == "catalog":
        print("[Stage 2] Running catalog mode: searching 6 catalog pools...", flush=True)
        for pool_name, pool_fna, is_required in catalog_pools:
            if not pool_fna.exists() or pool_fna.stat().st_size == 0:
                if is_required:
                    raise FileNotFoundError(f"Required target pool {pool_name} missing at {pool_fna}")
                continue
            out_blast = raw_stage2_dir / f"blast_{pool_name}.tsv"
            run_blastn_short(spacers_fna, pool_fna, out_blast, db_tmp_dir, threads=args.blast_threads * 4, query_sha256=q_sha)
        (raw_stage2_dir / "catalog_search.done.json").write_text(
            json.dumps({"status": "completed", "completed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}, indent=2) + "\n"
        )
        print("[Stage 2] Catalog mode completed.", flush=True)
        return

    # MODE: ALL (run catalog + search all assemblies + collect in single run)
    if args.mode == "all":
        cell_ctrl_ids, cell_ctrl_seqs, cell_ctrl_annot_by_sha = build_verified_curated_cellular_controls(
            cell_ctrl_fna, raw_stage2_dir, max_candidates=120, threads=args.blast_threads * 4
        )
        target_seq_map: dict[str, str] = dict(cell_ctrl_seqs)
        for pool_name, pool_fna, is_required in catalog_pools:
            if not pool_fna.exists() or pool_fna.stat().st_size == 0:
                if is_required:
                    raise FileNotFoundError(f"Required target pool {pool_name} missing at {pool_fna}")
                continue
            out_blast = raw_stage2_dir / f"blast_{pool_name}.tsv"
            run_blastn_short(spacers_fna, pool_fna, out_blast, db_tmp_dir, threads=args.blast_threads * 4, query_sha256=q_sha)

        ver_df_path = RES_DIR / "stage0_arch591_fasta_verification.tsv"
        ver_df = pd.read_csv(ver_df_path, sep="\t") if ver_df_path.exists() else pd.DataFrame()
        arch_samples = set(ver_df["source_sample"].dropna().unique()) if "source_sample" in ver_df.columns else set()
        spc_samples = {str(m).split("__")[0] for m in spc_df["mag_id"].dropna().unique() if "__" in str(m)}
        target_samples = sorted(arch_samples | spc_samples)

        cache_hashes_p = raw_stage2_dir / "target_assemblies_subject_hashes.json"
        subj_hashes = load_or_update_subject_hashes_cache(cache_hashes_p, target_samples, P_ASSEMBLIES_DIR)
        raw_asm_dir = raw_stage2_dir / "raw_assembly_blast"
        raw_asm_dir.mkdir(parents=True, exist_ok=True)

        tasks = []
        for s in target_samples:
            asm_p = P_ASSEMBLIES_DIR / s / f"{s}_contigs.fa"
            if not asm_p.exists() or asm_p.stat().st_size == 0:
                continue
            s_sha = subj_hashes.get(s, "")
            out_tsv = raw_asm_dir / f"raw_asm_{s}.blast.tsv"
            if not is_blast_receipt_valid(out_tsv, q_sha, asm_p, expected_subject_sha256=s_sha, expected_threads=args.blast_threads):
                tasks.append((s, asm_p, s_sha))

        if tasks:
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as ex:
                futs = [
                    ex.submit(search_single_raw_assembly, s, asm_p, spacers_fna, raw_asm_dir, db_tmp_dir, args.blast_threads, q_sha, s_sha)
                    for s, asm_p, s_sha in tasks
                ]
                for fut in concurrent.futures.as_completed(futs):
                    fut.result()

    # MODE: COLLECT (gather validated shards, audit manifest, run CheckV/geNomad, score alignments)
    if args.mode in ("all", "collect"):
        print("[Stage 2] Running collect mode...", flush=True)
        cell_ctrl_ids, cell_ctrl_seqs, cell_ctrl_annot_by_sha = build_verified_curated_cellular_controls(
            cell_ctrl_fna, raw_stage2_dir, max_candidates=120, threads=args.blast_threads * 4
        )
        target_seq_map: dict[str, str] = dict(cell_ctrl_seqs)

        ver_df_path = RES_DIR / "stage0_arch591_fasta_verification.tsv"
        manifest_tsv = RES_DIR / "stage2_target_assemblies_manifest.tsv"
        summary_json = RES_DIR / "stage2_assemblies_summary.json"
        raw_asm_dir = raw_stage2_dir / "raw_assembly_blast"
        cache_hashes_p = raw_stage2_dir / "target_assemblies_subject_hashes.json"

        # Merge any shard-level hash caches into master cache
        master_hashes: dict[str, Any] = {}
        if cache_hashes_p.exists():
            try:
                master_hashes = json.loads(cache_hashes_p.read_text())
            except Exception:
                master_hashes = {}
        for shard_file in sorted(raw_stage2_dir.glob("target_assemblies_subject_hashes_shard_*.json")):
            try:
                s_dict = json.loads(shard_file.read_text())
                if isinstance(s_dict, dict):
                    master_hashes.update(s_dict)
            except Exception:
                pass
        if master_hashes:
            tmp_c = Path(str(cache_hashes_p) + f".tmp.{os.getpid()}")
            tmp_c.write_text(json.dumps(master_hashes, indent=2) + "\n")
            tmp_c.replace(cache_hashes_p)

        cat_done = raw_stage2_dir / "catalog_search.done.json"
        if not cat_done.exists() and args.mode != "all":
            raise RuntimeError(f"Refusing collect: catalog search incomplete (missing {cat_done})")

        # The collector is the SOLE author of the manifest and summary JSON
        manifest_df, asm_summary = build_and_audit_target_assemblies_manifest(
            manifest_tsv=manifest_tsv,
            summary_json=summary_json,
            arch591_ver_tsv=ver_df_path,
            stage1_spacers_tsv=spacers_tsv,
            assemblies_dir=P_ASSEMBLIES_DIR,
            raw_asm_dir=raw_asm_dir,
            query_sha256=q_sha,
            subject_hashes_cache_path=cache_hashes_p,
            expected_threads=args.blast_threads,
        )

        # STRICT REFUSAL OF PARTIAL SHARDS
        if asm_summary["n_assemblies_completed_verified"] < asm_summary["n_assemblies_present"]:
            raise RuntimeError(
                f"Refusing partial aggregation: only {asm_summary['n_assemblies_completed_verified']}/"
                f"{asm_summary['n_assemblies_present']} assemblies have verified receipts in {raw_stage2_dir}."
            )

        # Collect catalog hits
        for pool_name, pool_fna, _ in catalog_pools:
            out_blast = raw_stage2_dir / f"blast_{pool_name}.tsv"
            n_pool_hits = 0
            if out_blast.exists() and out_blast.stat().st_size > 0:
                bdf = pd.read_csv(out_blast, sep="\t", names=BLAST_COLS)
                for _, r in bdf.iterrows():
                    sc = score_spacer_alignment(
                        qlen=int(r["qlen"]),
                        qstart=int(r["qstart"]),
                        qend=int(r["qend"]),
                        sstart=int(r["sstart"]),
                        send=int(r["send"]),
                        length=int(r["length"]),
                        mismatch=int(r["mismatch"]),
                        gapopen=int(r["gapopen"]),
                        pident=float(r["pident"]),
                        gaps=int(r["gaps"]),
                    )
                    if not sc["passes_exploratory_spacer_gate"]:
                        continue
                    s_id = str(r["sseqid"]).split()[0]
                    pool_matched_ids[pool_name].add(s_id)
                    n_pool_hits += 1
                    raw_hits.append(
                        {
                            "pool": pool_name,
                            "unique_spacer_id": str(r["qseqid"]).split("|")[0],
                            "target_id": s_id,
                            "raw_subject_id": s_id,
                            "source_sample": "",
                            "target_len": int(r["slen"]),
                            "evalue": float(r["evalue"]),
                            "bitscore": float(r["bitscore"]),
                            "btop": str(r["btop"]),
                            **sc,
                        }
                    )
            pool_status_rows.append(
                {
                    "pool": pool_name,
                    "source_path": str(pool_fna),
                    "status": "searched_complete" if out_blast.exists() else "unsearched",
                    "n_sequences": -1,
                    "n_bases": -1,
                    "n_exploratory_alignments": n_pool_hits,
                }
            )

        # Extract exact sequences for matched targets in catalog pools
        for row in pool_status_rows:
            p_name = row["pool"]
            p_path = Path(row["source_path"])
            if p_path.exists():
                n_s, n_b, ext = extract_sequences_from_fasta(p_path, pool_matched_ids[p_name])
                row["n_sequences"] = n_s
                row["n_bases"] = n_b
                target_seq_map.update(ext)

        # Collect raw assembly hits across all verified assemblies in parallel
        asm_total_seqs = 0
        asm_total_bases = 0
        asm_total_hits = 0
        subj_hashes = load_or_update_subject_hashes_cache(cache_hashes_p, list(manifest_df["sample_id"]), P_ASSEMBLIES_DIR)

        def _collect_single_asm(row_tuple: Any) -> tuple[int, int, list[dict[str, Any]], dict[str, str]] | None:
            _, r = row_tuple
            samp = r["sample_id"]
            asm_p = Path(r["assembly_path"])
            if not asm_p.exists():
                return None
            s_sha = subj_hashes.get(samp, "")
            out_tsv = raw_asm_dir / f"raw_asm_{samp}.blast.tsv"
            if not is_blast_receipt_valid(out_tsv, q_sha, asm_p, expected_subject_sha256=s_sha, expected_threads=args.blast_threads):
                return None
            _, n_s, n_b, h_list, ext_map = search_single_raw_assembly(
                samp, asm_p, spacers_fna, raw_asm_dir, db_tmp_dir, args.blast_threads, q_sha, s_sha
            )
            return (n_s, n_b, h_list, ext_map)

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.max_workers) as ex:
            asm_results = list(ex.map(_collect_single_asm, manifest_df.iterrows()))

        for res in asm_results:
            if res is None:
                continue
            n_s, n_b, h_list, ext_map = res
            asm_total_seqs += n_s
            asm_total_bases += n_b
            asm_total_hits += len(h_list)
            raw_hits.extend(h_list)
            target_seq_map.update(ext_map)

        pool_status_rows.append(
            {
                "pool": "eq_raw_assemblies_arch_samples",
                "source_path": f"{P_ASSEMBLIES_DIR}/<sample>/<sample>_contigs.fa ({asm_summary['n_assemblies_present']} present)",
                "status": "searched_complete" if asm_summary["is_present_coverage_complete"] else "incomplete_partial",
                "n_sequences": asm_total_seqs,
                "n_bases": asm_total_bases,
                "n_exploratory_alignments": asm_total_hits,
            }
        )
        pd.DataFrame(pool_status_rows).to_csv(
            RES_DIR / "stage2_pool_coverage_and_provenance.tsv", sep="\t", index=False
        )

    # 5. Run fresh CheckV (`checkv.sif`) and fresh geNomad (`genomad.sif`) directly on ALL
    #    matched primary target sequences (`EXP3-FINAL-PROVENANCE-REVIEW.txt` item 1).
    #    Never bind old precomputed annotation rows to current sequence SHA-256s by target ID.
    annot_by_sha: dict[str, dict[str, Any]] = dict(cell_ctrl_annot_by_sha)
    all_matched_fna = raw_stage2_dir / "all_matched_primary_targets.fna"
    n_freshly_annotated = 0
    seen_primary_shas: set[str] = set()
    with open(all_matched_fna, "w") as u_fh:
        for tid, seq in sorted(target_seq_map.items()):
            if tid.startswith("CELLCTRL__"):
                continue
            s_sha = sha256_bytes(seq.encode("ascii"))
            if s_sha not in annot_by_sha and s_sha not in seen_primary_shas:
                seen_primary_shas.add(s_sha)
                n_freshly_annotated += 1
                u_fh.write(f">{tid}\n{seq}\n")

    if n_freshly_annotated > 0:
        fresh_annot = run_fresh_checkv_and_genomad(
            all_matched_fna, raw_stage2_dir / "fresh_target_annotations", threads=args.blast_threads * 4
        )
        for s_sha, d_new in fresh_annot.items():
            annot_by_sha.setdefault(s_sha, {}).update(d_new)

    # 6. Load historical precomputed `checkv_votu` + `genomad_votu` strictly as a legacy ID-keyed audit table
    #    (`legacy_precomputed_votu_table_unverified_original_input_hash`), never used as primary annotation.
    legacy_votu_audit_by_id = load_legacy_precomputed_votu_audit(
        checkv_dir=P_CHECKV_VOTU_QUAL.parent,
        genomad_genes_path=P_GENOMAD_GENES,
        genomad_tax_path=P_GENOMAD_TAX,
        target_ids=set(target_seq_map.keys()),
    )
    legacy_audit_rows: list[dict[str, Any]] = []
    for tid in sorted(legacy_votu_audit_by_id.keys()):
        seq = target_seq_map.get(tid, "")
        t_sha = sha256_bytes(seq.encode("ascii")) if seq else ""
        fresh_d = annot_by_sha.get(t_sha, {})
        leg_d = legacy_votu_audit_by_id[tid]
        legacy_audit_rows.append(
            {
                "target_id": tid,
                "current_seq_sha256": t_sha,
                "primary_annotation_provenance": fresh_d.get("annotation_provenance", "unannotated"),
                "fresh_checkv_quality": fresh_d.get("checkv_quality", "Missing_unscored"),
                "fresh_checkv_viral_genes": fresh_d.get("checkv_viral_genes", "NA"),
                "fresh_virus_hallmarks": fresh_d.get("virus_hallmarks", "NA"),
                "fresh_genomad_lineage": fresh_d.get("genomad_lineage", "NA"),
                "legacy_audit_provenance": leg_d.get("legacy_audit_provenance", ""),
                "legacy_checkv_quality": leg_d.get("legacy_checkv_quality", "NA"),
                "legacy_checkv_viral_genes": leg_d.get("legacy_checkv_viral_genes", "NA"),
                "legacy_virus_hallmarks": leg_d.get("legacy_virus_hallmarks", "NA"),
                "legacy_genomad_lineage": leg_d.get("legacy_genomad_lineage", "NA"),
            }
        )
    pd.DataFrame(legacy_audit_rows).to_csv(
        RES_DIR / "stage2_legacy_vs_fresh_annotation_audit.tsv", sep="\t", index=False
    )

    # 7. Load vOTU cluster membership (ONLY for cluster ID & cluster taxonomy affiliation, NEVER transferring
    #    sequence-specific completeness, hallmarks, or host flanks to non-identical members!)
    contig_to_votu: dict[str, str] = {}
    if P_VOTU_MEM.exists():
        with open(P_VOTU_MEM) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    contig_to_votu[parts[0]] = parts[1]
                    contig_to_votu[parts[1]] = parts[1]

    # 8. Load CRISPR array coordinates across all assemblies (`e5/repeats.tsv`) and Stage 1 arrays
    contig_to_crispr_arrays: dict[str, list[tuple[int, int]]] = collections.defaultdict(list)
    if P_E5_REPEATS.exists():
        with open(P_E5_REPEATS) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                samp = (r.get("sample") or "").strip()
                c_raw = (r.get("contig") or r.get("seq_id") or "").strip()
                if c_raw:
                    _, _, q_cid, _, _ = qualify_namespaced_contig_id(c_raw, fallback_sample=samp)
                    s_lo = _opt_int(r.get("start") or r.get("repeat_start") or r.get("array_start"))
                    s_hi = _opt_int(r.get("end") or r.get("repeat_end") or r.get("array_end"))
                    if s_lo is not None and s_hi is not None:
                        contig_to_crispr_arrays[q_cid].append((min(s_lo, s_hi), max(s_lo, s_hi)))

    st1_arr = RES_DIR / "stage1_raw_crispr_arrays.tsv"
    if st1_arr.exists():
        adf = pd.read_csv(st1_arr, sep="\t")
        for _, r in adf.iterrows():
            contig_to_crispr_arrays[str(r["namespaced_contig_id"])].append(
                (int(r["array_start"]), int(r["array_end"]))
            )

    # 9. Load cellular contig index (`inputs/cellular_contig_to_mag.tsv.gz`)
    cell_idx_path = IN_DIR / "cellular_contig_to_mag.tsv.gz"
    namespaced_to_mags: dict[str, list[str]] = collections.defaultdict(list)
    if cell_idx_path.exists():
        with gzip.open(cell_idx_path, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                namespaced_to_mags[r["namespaced_contig_id"]].append(r["mag_id"])

    # 10. Score every alignment on Axis 2 (viral identity, self-array mask across ALL array loci of the spacer,
    #     cellular/provirus boundary, and controls)
    scored_rows: list[dict[str, Any]] = []
    for h in raw_hits:
        uid = h["unique_spacer_id"]
        spc_info = spc_lookup.get(uid, {})
        mag_id = str(spc_info.get("mag_id", ""))
        host_q_cid = str(
            spc_info.get("primary_namespaced_contig_id")
            or spc_info.get("namespaced_contig_id")
            or ""
        )
        arr_lo = int(
            spc_info.get("primary_array_start")
            or spc_info.get("array_start")
            or 0
        )
        arr_hi = int(
            spc_info.get("primary_array_end")
            or spc_info.get("array_end")
            or 0
        )
        all_loci_str = str(
            spc_info.get("all_array_loci")
            or (f"{host_q_cid}:{arr_lo}-{arr_hi}" if host_q_cid else "")
        )

        t_id = h["target_id"]
        pool_name = h["pool"]
        seq = target_seq_map.get(t_id, "")
        t_sha = sha256_bytes(seq.encode("ascii")) if seq else ""
        t_meta = dict(annot_by_sha.get(t_sha, {}))
        t_meta["target_contig_sha256"] = t_sha
        t_meta["source_sample"] = h.get("source_sample", "")

        if pool_name == "curated_plasmid_controls":
            t_meta["is_curated_plasmid_control"] = True
        elif pool_name == "curated_cellular_controls" or t_id in cell_ctrl_ids:
            t_meta["is_curated_cellular_control"] = True
        elif pool_name == "external_habitat_viral_controls":
            t_meta["is_external_habitat_viral_control"] = True

        clean_tid = t_id[len("CELLCTRL__"):] if t_id.startswith("CELLCTRL__") else t_id
        _, _, target_q_cid, _, _ = qualify_namespaced_contig_id(clean_tid, fallback_sample=h.get("source_sample", ""))

        cell_mags = sorted(set(namespaced_to_mags.get(target_q_cid, [])))
        if cell_mags:
            t_meta["in_cellular_mag_contig_index"] = True
            t_meta["cellular_mag_ids"] = ";".join(cell_mags)
        if mag_id and mag_id in cell_mags:
            t_meta["contig_in_same_host_mag"] = True

        # Check all known CRISPR arrays on `target_q_cid` (both from `contig_to_crispr_arrays` and `all_array_loci`)
        target_arrays = list(contig_to_crispr_arrays.get(target_q_cid, []))
        for loc_item in all_loci_str.split(";"):
            if ":" in loc_item:
                loc_c, loc_rng = loc_item.rsplit(":", 1)
                if loc_c == target_q_cid and "-" in loc_rng:
                    l_lo, l_hi = loc_rng.split("-", 1)
                    if l_lo.isdigit() and l_hi.isdigit():
                        target_arrays.append((int(l_lo), int(l_hi)))

        v_eval = evaluate_viral_identity_and_decoys(
            target_id=clean_tid,
            target_hit_lo=int(h["target_hit_lo"]),
            target_hit_hi=int(h["target_hit_hi"]),
            host_mag=mag_id,
            host_contig_id=host_q_cid,
            host_array_lo=arr_lo,
            host_array_hi=arr_hi,
            target_meta=t_meta,
            all_crispr_arrays_on_contig=target_arrays,
        )

        votu_rep_id = contig_to_votu.get(clean_tid, clean_tid)
        scored_rows.append(
            {
                **h,
                "target_seq_sha256": t_sha,
                "annotation_provenance": t_meta.get("annotation_provenance", "unannotated"),
                "cluster_votu_rep_id": votu_rep_id,
                "is_votu_cluster_rep": int(clean_tid == votu_rep_id),
                "host_mag": mag_id,
                "is_arch591_rep": int(spc_info.get("is_arch591_rep", 0) or 0),
                "host_namespaced_contig_id": host_q_cid,
                "host_all_array_loci": all_loci_str,
                "host_has_cas_operon": int(mag_id in cas_mags),
                "target_namespaced_contig_id": v_eval["target_qualified_contig"],
                "is_same_namespaced_contig": int(v_eval["is_same_namespaced_contig"]),
                "viral_identity_status": v_eval["viral_identity_status"],
                "habitat_specificity_role": v_eval["habitat_specificity_role"],
                "is_verified_viral": int(v_eval["is_verified_viral"]),
                "is_provisional_viral": int(v_eval["is_provisional_viral"]),
                "is_unknown_viral_identity": int(v_eval["is_unknown_viral_identity"]),
                "is_cellular_like": int(v_eval["is_cellular_like"]),
                "is_technical_crispr_array_self_match": int(v_eval["is_technical_crispr_array_self_match"]),
                "is_decoy": int(v_eval["is_decoy"]),
                "is_external_viral_specificity_control": int(v_eval["is_external_viral_specificity_control"]),
                "min_independent_viral_orfs_lower_bound": (
                    v_eval["min_independent_viral_orfs"]
                    if v_eval["min_independent_viral_orfs"] is not None
                    else "NA"
                ),
                "orf_count_provenance": t_meta.get("orf_count_provenance", "none"),
                "provirus_boundary_status": v_eval["provirus_boundary_status"],
                "read_junction_support_status": v_eval.get("read_junction_support_status", "untested_missing"),
                "status_detail": v_eval["status_detail"],
                "checkv_quality": t_meta.get("checkv_quality", "Missing_unscored"),
                "checkv_viral_genes": t_meta.get("checkv_viral_genes", "NA"),
                "checkv_host_genes": t_meta.get("checkv_host_genes", "NA"),
                "virus_hallmarks": t_meta.get("virus_hallmarks", "NA"),
                "archaeal_viral_hallmarks": t_meta.get("archaeal_viral_hallmarks", "NA"),
                "plasmid_hallmarks": t_meta.get("plasmid_hallmarks", "NA"),
            }
        )

    scored_df = pd.DataFrame(scored_rows)
    scored_df.to_csv(RES_DIR / "stage2_all_spacer_alignments_scored.tsv", sep="\t", index=False)

    # 11. Write target viral identity table across all extracted targets
    target_rows: list[dict[str, Any]] = []
    for t_id, seq in sorted(target_seq_map.items()):
        t_sha = sha256_bytes(seq.encode("ascii"))
        tm = dict(annot_by_sha.get(t_sha, {}))
        clean_tid = t_id[len("CELLCTRL__"):] if t_id.startswith("CELLCTRL__") else t_id
        if t_id in cell_ctrl_ids:
            tm["is_curated_cellular_control"] = True
        v_ev = evaluate_viral_identity_and_decoys(
            target_id=clean_tid,
            target_hit_lo=100,
            target_hit_hi=132,
            host_mag="",
            host_contig_id="",
            host_array_lo=0,
            host_array_hi=0,
            target_meta=tm,
            all_crispr_arrays_on_contig=[],
        )
        votu_rep_id = contig_to_votu.get(clean_tid, clean_tid)
        leg_d = legacy_votu_audit_by_id.get(clean_tid, {})
        target_rows.append(
            {
                "target_id": t_id,
                "target_seq_sha256": t_sha,
                "annotation_provenance": tm.get("annotation_provenance", "unannotated"),
                "cluster_votu_rep_id": votu_rep_id,
                "is_votu_cluster_rep": int(clean_tid == votu_rep_id),
                "target_namespaced_contig_id": v_ev["target_qualified_contig"],
                "viral_identity_status": v_ev["viral_identity_status"],
                "habitat_specificity_role": v_ev["habitat_specificity_role"],
                "is_verified_viral": int(v_ev["is_verified_viral"]),
                "is_provisional_viral": int(v_ev["is_provisional_viral"]),
                "is_unknown_viral_identity": int(v_ev["is_unknown_viral_identity"]),
                "is_cellular_like": int(v_ev["is_cellular_like"]),
                "is_decoy": int(v_ev["is_decoy"]),
                "is_external_viral_specificity_control": int(v_ev["is_external_viral_specificity_control"]),
                "min_independent_viral_orfs_lower_bound": (
                    v_ev["min_independent_viral_orfs"]
                    if v_ev["min_independent_viral_orfs"] is not None
                    else "NA"
                ),
                "orf_count_provenance": tm.get("orf_count_provenance", "none"),
                "checkv_quality": tm.get("checkv_quality", "Missing_unscored"),
                "checkv_viral_genes": tm.get("checkv_viral_genes", "NA"),
                "checkv_host_genes": tm.get("checkv_host_genes", "NA"),
                "virus_hallmarks": tm.get("virus_hallmarks", "NA"),
                "archaeal_viral_hallmarks": tm.get("archaeal_viral_hallmarks", "NA"),
                "plasmid_hallmarks": tm.get("plasmid_hallmarks", "NA"),
                "genomad_lineage": tm.get("genomad_lineage", "NA"),
                "provirus_boundary_status": v_ev["provirus_boundary_status"],
                "read_junction_support_status": v_ev.get("read_junction_support_status", "untested_missing"),
                "status_detail": v_ev["status_detail"],
                "legacy_audit_provenance": leg_d.get("legacy_audit_provenance", "none"),
                "legacy_audit_checkv_quality": leg_d.get("legacy_checkv_quality", "NA"),
                "legacy_audit_checkv_viral_genes": leg_d.get("legacy_checkv_viral_genes", "NA"),
                "legacy_audit_virus_hallmarks": leg_d.get("legacy_virus_hallmarks", "NA"),
            }
        )
    pd.DataFrame(target_rows).to_csv(RES_DIR / "stage2_target_viral_identity_table.tsv", sep="\t", index=False)

    # 12. Control panel calibration summary across match tiers
    cal_rows: list[dict[str, Any]] = []
    n_total_spacers = max(1, len(spc_df))
    all_pool_names = [p[0] for p in catalog_pools] + ["eq_raw_assemblies_arch_samples"]
    for p_name in all_pool_names:
        sub = scored_df[scored_df["pool"] == p_name] if len(scored_df) > 0 else pd.DataFrame()
        for tier in ("exact_0mm_full", "strict_1mm_end2", "near_2mm_exploratory"):
            sub_t = sub[sub["match_tier"] == tier] if len(sub) > 0 else pd.DataFrame()
            sub_unmasked = (
                sub_t[sub_t["is_technical_crispr_array_self_match"] == 0]
                if len(sub_t) > 0
                else pd.DataFrame()
            )
            n_spc_hit = int(sub_unmasked["unique_spacer_id"].nunique()) if len(sub_unmasked) > 0 else 0
            cal_rows.append(
                {
                    "pool": p_name,
                    "match_tier": tier,
                    "n_total_alignments": int(len(sub_t)),
                    "n_unmasked_alignments": int(len(sub_unmasked)),
                    "n_distinct_spacers_hitting": n_spc_hit,
                    "spacer_hit_rate": round(n_spc_hit / float(n_total_spacers), 6),
                }
            )
    pd.DataFrame(cal_rows).to_csv(RES_DIR / "stage2_control_panel_calibration.tsv", sep="\t", index=False)

    subprocess.run(["rm", "-rf", str(db_tmp_dir)], check=False)
    summary = {
        "stage": "stage2_viral_identity_and_spacers",
        "status": "completed",
        "raw_dir_name": args.raw_dir_name,
        "elapsed_seconds": round(time.time() - t0, 2),
        "n_verified_gold_cellular_controls": len(cell_ctrl_ids),
        "n_raw_assemblies_searched": asm_summary.get("n_assemblies_completed_verified", 0),
        "n_assemblies_present": asm_summary.get("n_assemblies_present", 0),
        "n_all_matched_primary_targets_freshly_annotated_by_checkv_genomad": n_freshly_annotated,
        "n_legacy_votu_audit_rows_retained_for_comparison": len(legacy_audit_rows),
        "n_total_scored_alignments": int(len(scored_df)),
        "n_strict_alignments_unmasked": int(
            len(
                scored_df[
                    (scored_df["passes_strict_spacer_gate"] == True)
                    & (scored_df["is_technical_crispr_array_self_match"] == 0)
                ]
            )
            if len(scored_df) > 0
            else 0
        ),
        "n_technical_self_array_masked": int(
            (scored_df["is_technical_crispr_array_self_match"] == 1).sum() if len(scored_df) > 0 else 0
        ),
        "n_candidate_targets_evaluated": int(len(target_rows)),
    }
    (RES_DIR / "stage2_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("[Stage 2] Summary:", json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
