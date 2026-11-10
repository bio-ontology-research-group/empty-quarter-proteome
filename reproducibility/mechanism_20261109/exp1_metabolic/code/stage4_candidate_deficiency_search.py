#!/usr/bin/env python3
"""Stage 4 of `exp1_metabolic` (v5 Repaired): Candidate-Conditioned Missing-Enzyme & Transporter Verification
(`EXP1-CONTRACT.txt` Item 4, `EXP1-FOLLOWUP-AND-GATE-COUNTEREXAMPLES.txt` Item B,
`EXP1-DEFICIENCY-REFERENCE-NAMESPACE-REPAIR.txt`, and `EXP1-DEFICIENCY-ROOT-SEMANTIC-REVIEW.txt`).

Repairs over the unperformed v4 job (`15894447`):
1. Authoritative ModelSEED-to-Gapseq Reference Namespace Crosswalk:
   - Maps candidate ModelSEED reaction IDs (`rxnXXXXX`) to actual gapseq reference sequence keys
     (`rev/<ec>.fasta`, `unrev/<ec>.fasta`, `rxn/<MetaCyc_RXN>.fasta`, `rxn/<UniProt>.fasta`,
     `rev/<md5(name)>.fasta`, `unrev/<md5(name)>.fasta`, and `user/*.fasta`) using installed gapseq
     crosswalk tables and per-genome `<g>-all-Reactions.tbl` tables.
   - Deduplicates exact reference protein sequences for DIAMOND search while preserving full
     many-to-many sequence-to-reaction ownership (`stage4_v5_reference_sequence_ownership.tsv`)
     and flagging shared-reference ambiguity (`is_ambiguous_many_to_many_reference`).
   - Explicitly qualifies reactions with non-shared references among the 63 selected candidate
     reactions as `selected_candidate_reference_unique` (unique within the 63 selected candidate
     reaction mappings only; broad EC crosswalks may span multiple substrate specificities and do
     not prove reaction-specific catalysis).
2. Single-HSP Reference Coverage & Qualified Evidence Labels (`EXP1-DEFICIENCY-ROOT-SEMANTIC-REVIEW.txt`):
   - Computes single-HSP reference span coverage as `(abs(send - sstart) + 1) / slen` (strictly
     bounded in `(0, 1]`), never `alignment_length / slen` (which counts gaps and can exceed 1.0).
   - Retains `alignment_length`, `qstart/qend/qlen` (DNA contig coordinates), and `sstart/send/slen`
     (reference protein coordinates) in both MAG and unbinned sample assembly searches.
   - Uses qualified dispositions:
     * `broad_reference_alignment_near_contig_edge` (near contig edge <=300 bp with single-HSP
       reference coverage >=0.70; candidate homolog, not claimed as a complete enzyme),
     * `candidate_truncated_or_split_gene_at_contig_edge` (near contig edge <=300 bp with
       single-HSP reference coverage <0.70),
     * `homolog_candidate_detected_in_MAG` (interior hit with bitscore >=50, pident >=40, and
       single-HSP reference coverage >=0.70; uncalled vs annotation miss is unknown without
       gene-call coordinate overlap verification),
     * `domain_level_or_partial_homolog_candidate_in_MAG` (interior hit with bitscore >=35 that
       has single-HSP reference coverage <0.70 or bitscore <50 or pident <40),
     * `unresolved_sample_community_homolog_host_unassigned` (sample assembly hit with bitscore >=65
       and single-HSP reference coverage >=0.70 without host linkage),
     * `unresolved_sample_community_domain_level_hit_host_unassigned` (sample assembly hit with
       bitscore >=65 but single-HSP reference coverage <0.70).
3. Fail-Closed Resource & Schema Validation:
   - Validates existence and schema of all required reference directories, crosswalk tables,
     MAG archives, sample assemblies, and pathway-hole synteny tables (`gap_context_candidates.tsv.gz`
     and `c18_pathway_holes/candidates.tsv` with `hole_rxn`).
   - Fails closed (`RuntimeError`) if candidate instances are non-empty and zero reactions resolve
     to reference sequences.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import hashlib
import io
import json
import math
import re
import subprocess
import tarfile
import time
from pathlib import Path
from typing import Any

import pandas as pd

EC4_RE = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$")
SEED_OR_RXN_PREFIX_RE = re.compile(r"^(R_)?rxn[0-9]+$|^RXN")


def compute_file_sha256(path: Path) -> str:
    """Compute SHA-256 hex digest of a file."""
    if not path.exists() or not path.is_file():
        return ""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(65536):
            h.update(chunk)
    return h.hexdigest()


def compute_text_sha256(text: str) -> str:
    """Compute SHA-256 hex digest of a UTF-8 string."""
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def compute_validated_single_hsp_reference_coverage(sstart: int, send: int, slen: int) -> float:
    """Compute single-HSP reference protein span coverage `(abs(send - sstart) + 1) / slen`.
    Never uses gapped alignment `length / slen` (which can exceed 1.0 when gaps are present).
    Validates coordinate bounds strictly (`1 <= sstart <= slen`, `1 <= send <= slen`, `slen >= 1`).
    """
    if slen < 1 or sstart < 1 or send < 1 or sstart > slen or send > slen:
        raise ValueError(
            f"Invalid DIAMOND subject coordinates for single-HSP reference coverage: "
            f"sstart={sstart}, send={send}, slen={slen}"
        )
    span = abs(send - sstart) + 1
    cov = float(span) / float(slen)
    if not (0.0 < cov <= 1.0):
        raise ValueError(
            f"Single-HSP reference coverage out of (0, 1] bounds: cov={cov} "
            f"(sstart={sstart}, send={send}, slen={slen})"
        )
    return cov


def parse_fasta_records(fasta_path: Path) -> list[tuple[str, str, str]]:
    """Parse a FASTA file into a list of `(seq_id_token, full_header, cleaned_sequence)`.
    Only retains records with non-empty amino acid sequence.
    """
    if not fasta_path.exists() or not fasta_path.is_file():
        return []
    records: list[tuple[str, str, str]] = []
    curr_id = ""
    curr_hdr = ""
    curr_seq_parts: list[str] = []

    for raw_line in fasta_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if curr_id:
                seq_str = "".join(curr_seq_parts)
                if seq_str:
                    records.append((curr_id, curr_hdr, seq_str))
            curr_hdr = line[1:].strip()
            curr_id = curr_hdr.split()[0] if curr_hdr else "unnamed"
            curr_seq_parts = []
        else:
            cleaned = re.sub(r"[^A-Za-z*]", "", line).upper()
            if cleaned:
                curr_seq_parts.append(cleaned)
    if curr_id:
        seq_str = "".join(curr_seq_parts)
        if seq_str:
            records.append((curr_id, curr_hdr, seq_str))
    return records


def _extract_ec4_tokens(raw_val: str) -> list[str]:
    """Extract valid 4-part EC numbers from a delimited string."""
    if not raw_val:
        return []
    out: list[str] = []
    for tok in re.split(r"[,;|/\s]+", str(raw_val).strip()):
        tok_clean = tok.strip()
        if EC4_RE.match(tok_clean):
            out.append(tok_clean)
    return sorted(set(out))


def _extract_pathway_tokens(raw_val: str) -> list[str]:
    """Extract clean MetaCyc pathway IDs (stripped of `|`) from a string."""
    if not raw_val:
        return []
    out: list[str] = []
    for tok in re.split(r"[,;\s|]+", str(raw_val).strip()):
        tok_clean = tok.strip().strip("|")
        if tok_clean and tok_clean != "NA":
            out.append(tok_clean)
    return sorted(set(out))


def resolve_gapseq_reference_crosswalk(
    target_clean_rxns: list[str],
    rxn_meta_by_clean_id: dict[str, dict[str, Any]],
    candidate_genomes: list[str],
    gapseq_dat_dir: Path,
    gapseq_seq_archaea_dir: Path,
    gapseq_user_seq_dir: Path | None,
    gapseq_out_dir: Path | None,
    out_dedup_fasta: Path,
) -> tuple[
    int,
    dict[str, dict[str, Any]],
    dict[str, dict[str, Any]],
    pd.DataFrame,
    pd.DataFrame,
    dict[str, str],
]:
    """Resolve ModelSEED `rxnXXXXX` IDs to gapseq Archaea `.fasta` files using authoritative
    gapseq crosswalk tables and per-genome `<g>-all-Reactions.tbl` tables.
    """
    if not gapseq_seq_archaea_dir.exists() or not gapseq_seq_archaea_dir.is_dir():
        raise FileNotFoundError(
            f"Gapseq Archaea sequence directory does not exist: {gapseq_seq_archaea_dir}"
        )

    rev_dir = gapseq_seq_archaea_dir / "rev"
    unrev_dir = gapseq_seq_archaea_dir / "unrev"
    rxn_dir = gapseq_seq_archaea_dir / "rxn"
    if not (rev_dir.is_dir() and unrev_dir.is_dir() and rxn_dir.is_dir()):
        raise FileNotFoundError(
            f"Gapseq Archaea sequence directory {gapseq_seq_archaea_dir} is missing required "
            f"'rev/', 'unrev/', or 'rxn/' subdirectories."
        )

    if not gapseq_dat_dir.exists() or not gapseq_dat_dir.is_dir():
        raise FileNotFoundError(
            f"Gapseq data directory does not exist: {gapseq_dat_dir}"
        )

    required_dat_files = [
        "meta_rea_pwy-gapseq.tbl",
        "seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv",
        "seed_reactions_corrected.tsv",
        "mnxref_seed-other.tsv",
        "meta_rea.tbl",
        "meta_genes.csv",
        "vmh_reactions.tsv",
        "SEED2VMH_translation_edited.csv",
        "seed_pwy.tbl",
        "corrections_seed_reactionDB.tsv",
    ]
    crosswalk_source_hashes: dict[str, str] = {}
    for fname in required_dat_files:
        fpath = gapseq_dat_dir / fname
        if not fpath.exists():
            raise FileNotFoundError(
                f"Required gapseq crosswalk table missing from {gapseq_dat_dir}: {fname}"
            )
        crosswalk_source_hashes[fname] = compute_file_sha256(fpath)

    target_set = sorted(set(target_clean_rxns))
    info: dict[str, dict[str, Any]] = {
        r: {
            "direct_rel_files": set(),
            "metacyc_rxns": set(),
            "ec_numbers": set(),
            "uniprot_ids": set(),
            "name_md5s": set(),
            "pathways": set(),
            "evidence_sources": set(),
            "seed_name": rxn_meta_by_clean_id.get(r, {}).get("rxn_name", ""),
            "is_transport": bool(rxn_meta_by_clean_id.get(r, {}).get("is_transport", False)),
            "gapseq_notes": [],
            "in_all_reactions_tbl": False,
        }
        for r in target_set
    }

    # 1. Per-genome <g>-all-Reactions.tbl across candidate genomes (if gapseq_out_dir provided)
    if gapseq_out_dir is not None:
        if not gapseq_out_dir.exists() or not gapseq_out_dir.is_dir():
            raise FileNotFoundError(f"Gapseq output directory does not exist: {gapseq_out_dir}")
        for g in sorted(set(candidate_genomes)):
            tgz = gapseq_out_dir / f"{g}.tar.gz"
            if not tgz.exists():
                raise FileNotFoundError(f"Missing candidate genome gapseq archive: {tgz}")
            member_name = f"{g}/{g}-all-Reactions.tbl"
            with tarfile.open(tgz, "r:gz") as tf:
                try:
                    m = tf.getmember(member_name)
                except KeyError as exc:
                    raise FileNotFoundError(
                        f"Archive {tgz} missing expected member {member_name}"
                    ) from exc
                fobj = tf.extractfile(m)
                if fobj is None:
                    raise RuntimeError(f"Unable to read {member_name} from {tgz}")
                text_stream = io.TextIOWrapper(fobj, encoding="utf-8", errors="replace")
                reader = csv.DictReader(
                    (ln for ln in text_stream if not ln.startswith("#")),
                    delimiter="\t",
                )
                for row in reader:
                    dbhit = (row.get("dbhit") or "").strip()
                    if not dbhit:
                        continue
                    hits = [tok for tok in dbhit.split() if tok in info]
                    if not hits:
                        continue
                    f_col = (row.get("file") or "").strip()
                    rxn_col = (row.get("rxn") or "").strip()
                    name_col = (row.get("name") or "").strip()
                    ec_col = (row.get("ec") or "").strip()
                    pwy_col = (row.get("pathway") or "").strip()

                    for h in hits:
                        d = info[h]
                        d["in_all_reactions_tbl"] = True
                        if f_col:
                            for ftok in re.split(r"[,;\s]+", f_col):
                                ftok_cl = ftok.strip()
                                if ftok_cl:
                                    d["direct_rel_files"].add(ftok_cl)
                                    d["evidence_sources"].add("all_reactions_tbl_file")
                        if rxn_col and not rxn_col.startswith("rxn"):
                            d["metacyc_rxns"].add(rxn_col)
                            d["evidence_sources"].add("all_reactions_tbl_metacyc_rxn")
                        for ec_tok in _extract_ec4_tokens(ec_col):
                            d["ec_numbers"].add(ec_tok)
                            d["evidence_sources"].add("all_reactions_tbl_ec")
                        for pwy_tok in _extract_pathway_tokens(pwy_col):
                            d["pathways"].add(pwy_tok)
                        if name_col and not SEED_OR_RXN_PREFIX_RE.match(rxn_col):
                            md5_hex = hashlib.md5(name_col.encode("utf-8")).hexdigest()
                            d["name_md5s"].add(md5_hex)
                            d["evidence_sources"].add("all_reactions_tbl_name_md5")

    # 2. meta_rea_pwy-gapseq.tbl
    with open(gapseq_dat_dir / "meta_rea_pwy-gapseq.tbl", "r", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gseq_ids = [tok for tok in (row.get("gapseq") or "").split() if tok in info]
            if not gseq_ids:
                continue
            rxn_col = (row.get("rxn") or "").strip()
            name_col = (row.get("name") or "").strip()
            ec_col = (row.get("ec") or "").strip()
            pwy_col = (row.get("pathway") or "").strip()
            for h in gseq_ids:
                d = info[h]
                for rtok in re.split(r"[,;\s]+", rxn_col):
                    rtok_cl = rtok.strip()
                    if rtok_cl and not rtok_cl.startswith("rxn"):
                        d["metacyc_rxns"].add(rtok_cl)
                        d["evidence_sources"].add("meta_rea_pwy_gapseq_rxn")
                for ec_tok in _extract_ec4_tokens(ec_col):
                    d["ec_numbers"].add(ec_tok)
                    d["evidence_sources"].add("meta_rea_pwy_gapseq_ec")
                for pwy_tok in _extract_pathway_tokens(pwy_col):
                    d["pathways"].add(pwy_tok)
                if name_col and not SEED_OR_RXN_PREFIX_RE.match(rxn_col):
                    md5_hex = hashlib.md5(name_col.encode("utf-8")).hexdigest()
                    d["name_md5s"].add(md5_hex)
                    d["evidence_sources"].add("meta_rea_pwy_gapseq_name_md5")

    # 3. seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv
    with open(
        gapseq_dat_dir / "seed_Enzyme_Class_Reactions_Aliases_unique_edited.tsv",
        "r",
        encoding="utf-8",
        errors="replace",
    ) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ms_ids = set((row.get("MS ID") or "").split("|")) | set(
                (row.get("Old MS ID") or "").split("|")
            )
            matched = [r for r in ms_ids if r in info]
            if not matched:
                continue
            ecs = _extract_ec4_tokens(row.get("External ID") or "")
            for h in matched:
                for ec_tok in ecs:
                    info[h]["ec_numbers"].add(ec_tok)
                    info[h]["evidence_sources"].add("seed_ec_aliases_edited")

    # 4. seed_reactions_corrected.tsv & corrections_seed_reactionDB.tsv
    for tbl_name, ev_tag in [
        ("seed_reactions_corrected.tsv", "seed_reactions_corrected"),
        ("corrections_seed_reactionDB.tsv", "corrections_seed_reactionDB"),
    ]:
        with open(gapseq_dat_dir / tbl_name, "r", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rid = (row.get("id") or "").strip()
                if rid not in info:
                    continue
                d = info[rid]
                if not d["seed_name"] and (row.get("name") or "").strip():
                    d["seed_name"] = (row.get("name") or "").strip()
                for ec_tok in _extract_ec4_tokens(row.get("ec_numbers") or ""):
                    d["ec_numbers"].add(ec_tok)
                    d["evidence_sources"].add(f"{ev_tag}_ec")
                note_val = (row.get("notes") or "").strip()
                if note_val and note_val not in d["gapseq_notes"]:
                    d["gapseq_notes"].append(note_val)

    # 5. mnxref_seed-other.tsv + meta_rea.tbl
    meta_rea_by_id: dict[str, dict[str, str]] = {}
    meta_rea_by_kegg: dict[str, list[dict[str, str]]] = {}
    with open(gapseq_dat_dir / "meta_rea.tbl", "r", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            mid = (row.get("id") or "").replace("META:", "").strip()
            if mid:
                meta_rea_by_id[mid] = row
            for ktok in re.split(r"[,;\s]+", (row.get("kegg") or "").strip()):
                if ktok:
                    meta_rea_by_kegg.setdefault(ktok, []).append(row)

    with open(gapseq_dat_dir / "mnxref_seed-other.tsv", "r", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s_id = (row.get("seed") or "").strip()
            if s_id not in info:
                continue
            other = (row.get("other") or "").strip()
            if not other:
                continue
            d = info[s_id]
            if other in meta_rea_by_id:
                mr = meta_rea_by_id[other]
                d["metacyc_rxns"].add(other)
                d["evidence_sources"].add("mnxref_meta_rea_rxn")
                for ec_tok in _extract_ec4_tokens(mr.get("ec") or ""):
                    d["ec_numbers"].add(ec_tok)
                    d["evidence_sources"].add("mnxref_meta_rea_ec")
                for pwy_tok in _extract_pathway_tokens(mr.get("pwy") or ""):
                    d["pathways"].add(pwy_tok)
            elif other in meta_rea_by_kegg:
                for mr in meta_rea_by_kegg[other]:
                    mid = (mr.get("id") or "").replace("META:", "").strip()
                    if mid:
                        d["metacyc_rxns"].add(mid)
                        d["evidence_sources"].add("mnxref_kegg_meta_rea_rxn")
                    for ec_tok in _extract_ec4_tokens(mr.get("ec") or ""):
                        d["ec_numbers"].add(ec_tok)
                        d["evidence_sources"].add("mnxref_kegg_meta_rea_ec")
                    for pwy_tok in _extract_pathway_tokens(mr.get("pwy") or ""):
                        d["pathways"].add(pwy_tok)
            elif not other.startswith("R") and not other.startswith("rxn"):
                d["metacyc_rxns"].add(other)
                d["evidence_sources"].add("mnxref_seed_other_rxn")

    # 6. vmh_reactions.tsv + SEED2VMH_translation_edited.csv + seed_pwy.tbl
    vmh_to_seed: dict[str, set[str]] = {}
    with open(
        gapseq_dat_dir / "SEED2VMH_translation_edited.csv",
        "r",
        encoding="utf-8",
        errors="replace",
    ) as fh:
        for line in fh:
            parts = line.strip().split(",")
            if len(parts) >= 2 and parts[0].strip() in info:
                vmh_to_seed.setdefault(parts[1].strip(), set()).add(parts[0].strip())

    with open(gapseq_dat_dir / "vmh_reactions.tsv", "r", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s_id = (row.get("seed") or "").strip()
            abbr = (row.get("abbreviation") or "").strip()
            matched_seeds: set[str] = set()
            if s_id in info:
                matched_seeds.add(s_id)
            if abbr in vmh_to_seed:
                matched_seeds |= vmh_to_seed[abbr]
            if not matched_seeds:
                continue
            ecs = _extract_ec4_tokens(row.get("ecnumber") or "")
            for h in matched_seeds:
                for ec_tok in ecs:
                    info[h]["ec_numbers"].add(ec_tok)
                    info[h]["evidence_sources"].add("vmh_reactions_ec")

    with open(gapseq_dat_dir / "seed_pwy.tbl", "r", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gseq_ids = [tok for tok in (row.get("gapseq") or "").split() if tok in info]
            if not gseq_ids:
                continue
            ecs = _extract_ec4_tokens(row.get("ec") or "")
            for h in gseq_ids:
                for ec_tok in ecs:
                    info[h]["ec_numbers"].add(ec_tok)
                    info[h]["evidence_sources"].add("seed_pwy_ec")

    # 7. meta_genes.csv: map MetaCyc rxn -> UniProt IDs (excluding Q91194 per prepare_batch_alignments.R)
    metacyc_to_uniprot: dict[str, set[str]] = {}
    with open(gapseq_dat_dir / "meta_genes.csv", "r", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter=","):
            mr = (row.get("rxn") or "").replace("META:", "").strip()
            up = (row.get("uniprot") or "").strip()
            if mr and up and up != "Q91194":
                metacyc_to_uniprot.setdefault(mr, set()).add(up)

    for r in target_set:
        d = info[r]
        for mcyc in sorted(d["metacyc_rxns"]):
            for up in sorted(metacyc_to_uniprot.get(mcyc, set())):
                d["uniprot_ids"].add(up)
                d["evidence_sources"].add("meta_genes_uniprot")

    # 8. Resolve candidate .fasta paths per reaction and parse sequences
    rxn_to_valid_files: dict[str, list[tuple[str, Path, str, list[tuple[str, str, str]]]]] = {}
    fasta_cache: dict[Path, tuple[str, list[tuple[str, str, str]]]] = {}

    def _inspect_fasta(fpath: Path) -> tuple[str, list[tuple[str, str, str]]]:
        if fpath not in fasta_cache:
            if fpath.exists() and fpath.is_file():
                recs = parse_fasta_records(fpath)
                fsha = compute_file_sha256(fpath) if recs else ""
                fasta_cache[fpath] = (fsha, recs)
            else:
                fasta_cache[fpath] = ("", [])
        return fasta_cache[fpath]

    for r in target_set:
        d = info[r]
        candidate_paths: list[tuple[str, Path]] = []
        seen_paths: set[Path] = set()

        def _add_cand(rel_label: str, p: Path) -> None:
            if p not in seen_paths:
                seen_paths.add(p)
                candidate_paths.append((rel_label, p))

        for relf in sorted(d["direct_rel_files"]):
            _add_cand(relf, gapseq_seq_archaea_dir / relf)

        for ec in sorted(d["ec_numbers"]):
            _add_cand(f"rev/{ec}.fasta", rev_dir / f"{ec}.fasta")
            _add_cand(f"unrev/{ec}.fasta", unrev_dir / f"{ec}.fasta")
            if gapseq_user_seq_dir is not None:
                _add_cand(f"user/{ec}.fasta", gapseq_user_seq_dir / f"{ec}.fasta")

        for mcyc in sorted(d["metacyc_rxns"]):
            _add_cand(f"rxn/{mcyc}.fasta", rxn_dir / f"{mcyc}.fasta")
            if gapseq_user_seq_dir is not None:
                _add_cand(f"user/{mcyc}.fasta", gapseq_user_seq_dir / f"{mcyc}.fasta")

        for up in sorted(d["uniprot_ids"]):
            _add_cand(f"rxn/{up}.fasta", rxn_dir / f"{up}.fasta")

        for md5_hex in sorted(d["name_md5s"]):
            _add_cand(f"rev/{md5_hex}.fasta", rev_dir / f"{md5_hex}.fasta")
            _add_cand(f"unrev/{md5_hex}.fasta", unrev_dir / f"{md5_hex}.fasta")
            if gapseq_user_seq_dir is not None:
                _add_cand(f"user/{md5_hex}.fasta", gapseq_user_seq_dir / f"{md5_hex}.fasta")

        valid_entries: list[tuple[str, Path, str, list[tuple[str, str, str]]]] = []
        for rel_label, fpath in candidate_paths:
            fsha, recs = _inspect_fasta(fpath)
            if recs:
                valid_entries.append((rel_label, fpath, fsha, recs))
        rxn_to_valid_files[r] = valid_entries

    # 9. Build deduplicated sequence catalog and many-to-many ownership tables
    seq_by_sha: dict[str, dict[str, Any]] = {}
    ownership_rows: list[dict[str, Any]] = []

    for r in target_set:
        for rel_label, fpath, fsha, recs in rxn_to_valid_files[r]:
            for raw_id, raw_hdr, seq_str in recs:
                seq_sha = compute_text_sha256(seq_str)
                if seq_sha not in seq_by_sha:
                    idx = len(seq_by_sha) + 1
                    dedup_id = f"refseq_{idx:06d}|{seq_sha[:16]}"
                    seq_by_sha[seq_sha] = {
                        "dedup_seq_id": dedup_id,
                        "seq_sha256": seq_sha,
                        "seq_len": len(seq_str),
                        "sequence": seq_str,
                        "clean_rxn_ids": set(),
                        "source_relpaths_by_rxn": {},
                        "source_headers": set(),
                        "primary_raw_id": raw_id,
                        "primary_raw_header": raw_hdr,
                    }
                entry = seq_by_sha[seq_sha]
                entry["clean_rxn_ids"].add(r)
                entry["source_relpaths_by_rxn"].setdefault(r, set()).add(rel_label)
                entry["source_headers"].add(raw_hdr)

                ownership_rows.append({
                    "clean_rxn_id": r,
                    "dedup_seq_id": entry["dedup_seq_id"],
                    "seq_sha256": seq_sha,
                    "seq_len": len(seq_str),
                    "fasta_relpath": rel_label,
                    "fasta_sha256": fsha,
                    "raw_seq_id": raw_id,
                    "raw_header": raw_hdr,
                })

    # Write deduplicated FASTA file
    out_dedup_fasta.parent.mkdir(parents=True, exist_ok=True)
    seq_ownership: dict[str, dict[str, Any]] = {}
    with open(out_dedup_fasta, "w", encoding="utf-8") as out_fh:
        for seq_sha, entry in sorted(seq_by_sha.items(), key=lambda kv: kv[1]["dedup_seq_id"]):
            dedup_id = entry["dedup_seq_id"]
            out_fh.write(f">{dedup_id}\n{entry['sequence']}\n")
            owners_sorted = sorted(entry["clean_rxn_ids"])
            seq_ownership[dedup_id] = {
                "dedup_seq_id": dedup_id,
                "seq_sha256": seq_sha,
                "seq_len": entry["seq_len"],
                "clean_rxn_ids": owners_sorted,
                "is_shared_across_candidate_rxns": len(owners_sorted) > 1,
                "source_relpaths_by_rxn": {
                    k: sorted(v) for k, v in sorted(entry["source_relpaths_by_rxn"].items())
                },
                "primary_raw_id": entry["primary_raw_id"],
                "primary_raw_header": entry["primary_raw_header"],
            }

    for row in ownership_rows:
        owners = seq_ownership[row["dedup_seq_id"]]["clean_rxn_ids"]
        row["n_owning_candidate_rxns"] = len(owners)
        row["owning_candidate_rxn_ids_csv"] = ",".join(owners)
        row["is_shared_sequence_across_candidate_rxns"] = len(owners) > 1

    # 10. Build per-reaction ref_inventory and crosswalk_df
    ref_inventory: dict[str, dict[str, Any]] = {}
    crosswalk_rows: list[dict[str, Any]] = []

    for r in target_set:
        d = info[r]
        valid_entries = rxn_to_valid_files[r]
        relpaths = [e[0] for e in valid_entries]
        fshas = [e[2] for e in valid_entries]
        n_total_records = sum(len(e[3]) for e in valid_entries)

        rxn_seq_shas = sorted({
            compute_text_sha256(rec[2])
            for e in valid_entries
            for rec in e[3]
        })
        shared_rxns: set[str] = set()
        for s_sha in rxn_seq_shas:
            for other_r in seq_by_sha[s_sha]["clean_rxn_ids"]:
                if other_r != r:
                    shared_rxns.add(other_r)

        has_ref = len(rxn_seq_shas) > 0
        is_shared = len(shared_rxns) > 0
        has_any_keys = bool(d["metacyc_rxns"] or d["ec_numbers"] or d["uniprot_ids"] or d["direct_rel_files"])

        if has_ref:
            if is_shared:
                res_status = "resolved_shared_many_to_many_reference_set"
                unavail_reason = ""
            else:
                res_status = "resolved_selected_candidate_reference_unique"
                unavail_reason = ""
        else:
            if has_any_keys:
                res_status = "unresolved_mapped_keys_absent_from_archaea_seq_db"
                unavail_reason = "mapped_metacyc_or_ec_keys_have_no_fasta_in_archaea_seq_db"
            else:
                res_status = "unresolved_putative_or_unindexed_in_gapseq_seq_db"
                unavail_reason = "no_ec_or_metacyc_key_with_sequences_in_gapseq_archaea_db"

        notes_str = "; ".join(d["gapseq_notes"])
        entry_meta = {
            "clean_rxn_id": r,
            "rxn_name": d["seed_name"],
            "is_transport": d["is_transport"],
            "reference_available": has_ref,
            "reference_resolution_status": res_status,
            "reference_uniqueness_scope_note": (
                "selected_candidate_reference_unique denotes non-sharing among the 63 selected "
                "candidate reactions only; broad EC/MetaCyc reference sets may span multiple "
                "substrate specificities and do not prove reaction-specific catalysis."
            ),
            "unavailability_reason": unavail_reason,
            "is_ambiguous_many_to_many_reference": is_shared,
            "shared_with_candidate_rxn_ids": sorted(shared_rxns),
            "shared_with_candidate_rxn_ids_csv": ",".join(sorted(shared_rxns)),
            "in_per_genome_all_reactions_tbl": d["in_all_reactions_tbl"],
            "mapped_metacyc_rxn_ids": sorted(d["metacyc_rxns"]),
            "mapped_ec_numbers": sorted(d["ec_numbers"]),
            "mapped_uniprot_ids": sorted(d["uniprot_ids"]),
            "mapped_metacyc_pathways": sorted(d["pathways"]),
            "mapped_fasta_relpaths": relpaths,
            "mapped_fasta_sha256s": fshas,
            "mapping_evidence_sources": sorted(d["evidence_sources"]),
            "n_mapped_fasta_files": len(valid_entries),
            "n_reference_sequences": len(rxn_seq_shas),
            "n_raw_fasta_records": n_total_records,
            "gapseq_notes": notes_str,
        }
        ref_inventory[r] = entry_meta
        crosswalk_rows.append({
            "clean_rxn_id": r,
            "rxn_name": d["seed_name"],
            "is_transport": d["is_transport"],
            "reference_available": has_ref,
            "reference_resolution_status": res_status,
            "unavailability_reason": unavail_reason,
            "is_ambiguous_many_to_many_reference": is_shared,
            "shared_with_candidate_rxn_ids_csv": ",".join(sorted(shared_rxns)),
            "in_per_genome_all_reactions_tbl": d["in_all_reactions_tbl"],
            "mapped_metacyc_rxn_ids_csv": ",".join(sorted(d["metacyc_rxns"])),
            "mapped_ec_numbers_csv": ",".join(sorted(d["ec_numbers"])),
            "mapped_uniprot_ids_csv": ",".join(sorted(d["uniprot_ids"])),
            "mapped_metacyc_pathways_csv": ",".join(sorted(d["pathways"])),
            "mapped_fasta_relpaths_csv": ",".join(relpaths),
            "mapped_fasta_sha256s_csv": ",".join(fshas),
            "mapping_evidence_sources_csv": ",".join(sorted(d["evidence_sources"])),
            "n_mapped_fasta_files": len(valid_entries),
            "n_dedup_reference_sequences": len(rxn_seq_shas),
            "n_raw_fasta_records": n_total_records,
            "gapseq_notes": notes_str,
        })

    crosswalk_df = pd.DataFrame(crosswalk_rows)
    ownership_df = pd.DataFrame(
        ownership_rows,
        columns=[
            "clean_rxn_id",
            "dedup_seq_id",
            "seq_sha256",
            "seq_len",
            "fasta_relpath",
            "fasta_sha256",
            "raw_seq_id",
            "raw_header",
            "n_owning_candidate_rxns",
            "owning_candidate_rxn_ids_csv",
            "is_shared_sequence_across_candidate_rxns",
        ],
    )
    return (
        len(seq_by_sha),
        ref_inventory,
        seq_ownership,
        crosswalk_df,
        ownership_df,
        crosswalk_source_hashes,
    )


def load_pathway_hole_synteny_candidates(
    gap_context_tsv_gz: Path,
    c18_candidates_tsv: Path,
    c18_windows_tsv: Path | None = None,
) -> dict[str, Any]:
    """Validate schemas fail-closed and load pathway-hole neighbourhood candidates."""
    if not gap_context_tsv_gz.exists() or not gap_context_tsv_gz.is_file():
        raise FileNotFoundError(
            f"Required pathway-hole synteny context file missing: {gap_context_tsv_gz}"
        )
    if not c18_candidates_tsv.exists() or not c18_candidates_tsv.is_file():
        raise FileNotFoundError(
            f"Required c18 pathway-hole candidates file missing: {c18_candidates_tsv}"
        )

    gap_by_genome_mcyc_rxn_dark: dict[tuple[str, str], list[str]] = {}
    gap_by_genome_mcyc_rxn_all: dict[tuple[str, str], list[str]] = {}
    gap_by_genome_ec_dark: dict[tuple[str, str], list[str]] = {}
    gap_by_genome_pwy_dark: dict[tuple[str, str], list[str]] = {}
    gap_by_genome_pwy_all: dict[tuple[str, str], list[str]] = {}

    with gzip.open(gap_context_tsv_gz, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        req_gap_cols = {"genome", "pathway", "missing_reaction", "ec", "candidate_gene", "dark"}
        missing_cols = req_gap_cols - fieldnames
        if missing_cols:
            raise ValueError(
                f"Invalid schema in {gap_context_tsv_gz}: missing columns {sorted(missing_cols)}; "
                f"found {sorted(fieldnames)}"
            )
        for row in reader:
            g = (row.get("genome") or "").strip()
            pwy = (row.get("pathway") or "").strip().strip("|")
            m_rxn = (row.get("missing_reaction") or "").strip()
            ec = (row.get("ec") or "").strip()
            gene = (row.get("candidate_gene") or "").strip()
            is_dark = str(row.get("dark") or "").strip().lower() == "true"
            if not g or not gene:
                continue
            if m_rxn:
                gap_by_genome_mcyc_rxn_all.setdefault((g, m_rxn), []).append(gene)
                if is_dark:
                    gap_by_genome_mcyc_rxn_dark.setdefault((g, m_rxn), []).append(gene)
            if ec and is_dark:
                for ec_tok in _extract_ec4_tokens(ec):
                    gap_by_genome_ec_dark.setdefault((g, ec_tok), []).append(gene)
            if pwy:
                gap_by_genome_pwy_all.setdefault((g, pwy), []).append(gene)
                if is_dark:
                    gap_by_genome_pwy_dark.setdefault((g, pwy), []).append(gene)

    c18_by_genome_rxn_dark: dict[tuple[str, str], list[str]] = {}
    c18_by_genome_rxn_all: dict[tuple[str, str], list[str]] = {}

    with open(c18_candidates_tsv, "r", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        req_c18_cols = {"genome", "hole_rxn", "candidate_gene", "dark"}
        missing_c18 = req_c18_cols - fieldnames
        if missing_c18:
            raise ValueError(
                f"Invalid schema in {c18_candidates_tsv}: missing columns {sorted(missing_c18)}; "
                f"found {sorted(fieldnames)}"
            )
        for row in reader:
            g = (row.get("genome") or "").strip()
            hrxn = (row.get("hole_rxn") or "").strip()
            gene = (row.get("candidate_gene") or "").strip()
            is_dark = str(row.get("dark") or "").strip().lower() == "true"
            if not g or not hrxn or not gene:
                continue
            c18_by_genome_rxn_all.setdefault((g, hrxn), []).append(gene)
            if is_dark:
                c18_by_genome_rxn_dark.setdefault((g, hrxn), []).append(gene)

    win_by_genome_pwy_dark_fams: dict[tuple[str, str], list[str]] = {}
    if c18_windows_tsv is not None and c18_windows_tsv.exists():
        with open(c18_windows_tsv, "r", encoding="utf-8", errors="replace") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            fieldnames = set(reader.fieldnames or [])
            req_win_cols = {"genome", "pathway", "has_hole", "n_window_genes", "n_dark", "dark_families"}
            missing_win = req_win_cols - fieldnames
            if missing_win:
                raise ValueError(
                    f"Invalid schema in {c18_windows_tsv}: missing columns {sorted(missing_win)}; "
                    f"found {sorted(fieldnames)}"
                )
            for row in reader:
                g = (row.get("genome") or "").strip()
                pwy = (row.get("pathway") or "").strip().strip("|")
                has_hole = str(row.get("has_hole") or "").strip() == "1"
                dfams = (row.get("dark_families") or "").strip()
                if g and pwy and has_hole and dfams:
                    for fam in [t.strip() for t in dfams.split(",") if t.strip()]:
                        win_by_genome_pwy_dark_fams.setdefault((g, pwy), []).append(fam)

    return {
        "gap_by_genome_mcyc_rxn_dark": gap_by_genome_mcyc_rxn_dark,
        "gap_by_genome_mcyc_rxn_all": gap_by_genome_mcyc_rxn_all,
        "gap_by_genome_ec_dark": gap_by_genome_ec_dark,
        "gap_by_genome_pwy_dark": gap_by_genome_pwy_dark,
        "gap_by_genome_pwy_all": gap_by_genome_pwy_all,
        "c18_by_genome_rxn_dark": c18_by_genome_rxn_dark,
        "c18_by_genome_rxn_all": c18_by_genome_rxn_all,
        "win_by_genome_pwy_dark_fams": win_by_genome_pwy_dark_fams,
    }


def lookup_synteny_context_for_instance(
    genome: str,
    clean_rxn: str,
    instance_pathway_id: str,
    ref_meta: dict[str, Any],
    synteny_index: dict[str, Any],
) -> dict[str, Any]:
    """Join `(genome, clean_rxn)` to pathway-hole context across `c18_pathway_holes/candidates.tsv`
    (direct ModelSEED `hole_rxn`) and `gap_context_candidates.tsv.gz` / `windows.tsv`.
    """
    c18_dark = list(dict.fromkeys(synteny_index["c18_by_genome_rxn_dark"].get((genome, clean_rxn), [])))
    c18_all = list(dict.fromkeys(synteny_index["c18_by_genome_rxn_all"].get((genome, clean_rxn), [])))

    mcyc_rxns = ref_meta.get("mapped_metacyc_rxn_ids", [])
    ecs = ref_meta.get("mapped_ec_numbers", [])
    pwys = set(ref_meta.get("mapped_metacyc_pathways", []))
    if instance_pathway_id:
        for p in _extract_pathway_tokens(instance_pathway_id):
            pwys.add(p)
    pwys_sorted = sorted(pwys)

    gap_rxn_dark: list[str] = []
    gap_rxn_all: list[str] = []
    for mr in mcyc_rxns:
        gap_rxn_dark.extend(synteny_index["gap_by_genome_mcyc_rxn_dark"].get((genome, mr), []))
        gap_rxn_all.extend(synteny_index["gap_by_genome_mcyc_rxn_all"].get((genome, mr), []))
    for ec in ecs:
        gap_rxn_dark.extend(synteny_index["gap_by_genome_ec_dark"].get((genome, ec), []))

    gap_pwy_dark: list[str] = []
    gap_pwy_all: list[str] = []
    win_dark_fams: list[str] = []
    for pwy in pwys_sorted:
        gap_pwy_dark.extend(synteny_index["gap_by_genome_pwy_dark"].get((genome, pwy), []))
        gap_pwy_all.extend(synteny_index["gap_by_genome_pwy_all"].get((genome, pwy), []))
        win_dark_fams.extend(synteny_index["win_by_genome_pwy_dark_fams"].get((genome, pwy), []))

    combined_dark = list(dict.fromkeys(c18_dark + gap_rxn_dark + gap_pwy_dark))
    rxn_specific_dark = list(dict.fromkeys(c18_dark + gap_rxn_dark))
    combined_all = list(dict.fromkeys(c18_all + gap_rxn_all + gap_pwy_all))
    win_dark_fams_dedup = list(dict.fromkeys(win_dark_fams))

    evidence_types: list[str] = []
    if c18_dark:
        evidence_types.append("c18_hole_rxn_dark_candidate")
    elif c18_all:
        evidence_types.append("c18_hole_rxn_non_dark_candidate")
    if gap_rxn_dark:
        evidence_types.append("gap_context_missing_rxn_dark_candidate")
    if gap_pwy_dark:
        evidence_types.append("gap_context_pathway_dark_candidate")
    if win_dark_fams_dedup:
        evidence_types.append("c18_window_dark_family")

    return {
        "synteny_dark_genes": combined_dark,
        "rxn_specific_dark_genes": rxn_specific_dark,
        "c18_hole_dark_genes": c18_dark,
        "c18_hole_all_genes": c18_all,
        "gap_rxn_dark_genes": list(dict.fromkeys(gap_rxn_dark)),
        "gap_pwy_dark_genes": list(dict.fromkeys(gap_pwy_dark)),
        "all_context_candidate_genes": combined_all,
        "window_dark_families": win_dark_fams_dedup,
        "synteny_evidence_types_csv": ",".join(evidence_types),
    }


def extract_refined_mags_and_run_diamond_blastx(
    refined_mags_tar: Path,
    candidate_genomes: list[str],
    db_dmnd: Path,
    seq_ownership: dict[str, dict[str, Any]],
    raw_searches_dir: Path,
    tmp_dir: Path,
    threads: int,
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[tuple[str, str], int], str, dict[str, Any]]:
    """Extract refined MAG contigs, verify all candidate genomes are present in the archive,
    and run `diamond blastx` (`--max-target-seqs 0` to prevent pooled-DB hit masking).
    Computes single-HSP reference coverage as `(abs(send - sstart) + 1) / slen`.
    """
    if not refined_mags_tar.exists() or not refined_mags_tar.is_file():
        raise FileNotFoundError(f"Refined MAGs archive does not exist: {refined_mags_tar}")

    mags_fasta = tmp_dir / "all_refined_mags.fna"
    contig_lengths: dict[tuple[str, str], int] = {}
    archive_genomes: set[str] = set()

    with tarfile.open(refined_mags_tar, "r:gz") as tf, open(mags_fasta, "w", encoding="utf-8") as out_fh:
        for member in tf.getmembers():
            if not member.isfile() or not (
                member.name.endswith(".fa")
                or member.name.endswith(".fna")
                or member.name.endswith(".fasta")
            ):
                continue
            g_name = Path(member.name).stem
            archive_genomes.add(g_name)
            fobj = tf.extractfile(member)
            if fobj is None:
                continue
            text = fobj.read().decode("utf-8", errors="replace")
            curr_cid = ""
            curr_len = 0
            for ln in text.splitlines():
                if ln.startswith(">"):
                    if curr_cid:
                        contig_lengths[(g_name, curr_cid)] = curr_len
                    curr_cid = ln[1:].strip().split()[0]
                    curr_len = 0
                    out_fh.write(f">{g_name}||{curr_cid}\n")
                else:
                    s = ln.strip()
                    curr_len += len(s)
                    out_fh.write(s + "\n")
            if curr_cid:
                contig_lengths[(g_name, curr_cid)] = curr_len

    missing_cand_genomes = sorted(set(candidate_genomes) - archive_genomes)
    if missing_cand_genomes:
        raise RuntimeError(
            f"Refined MAGs archive {refined_mags_tar} is missing {len(missing_cand_genomes)} "
            f"candidate genomes: {missing_cand_genomes[:10]}"
        )

    out_tsv = raw_searches_dir / "mag_contigs_vs_target_rxns.tsv"
    cmd = [
        "diamond",
        "blastx",
        "--query",
        str(mags_fasta),
        "--db",
        str(db_dmnd),
        "--out",
        str(out_tsv),
        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "qstart",
        "qend",
        "qlen",
        "sstart",
        "send",
        "slen",
        "evalue",
        "bitscore",
        "--evalue",
        "1e-3",
        "--max-target-seqs",
        "0",
        "--threads",
        str(threads),
        "--tmpdir",
        str(tmp_dir),
    ]
    subprocess.run(cmd, check=True)

    best_hits: dict[tuple[str, str], dict[str, Any]] = {}
    n_raw_alignments = 0
    if out_tsv.exists():
        with open(out_tsv, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 12:
                    continue
                n_raw_alignments += 1
                qseqid, sseqid, pident, length, qstart, qend, qlen, sstart, send, slen, evalue, bitscore = parts
                if "||" not in qseqid:
                    continue
                g_name, cid = qseqid.split("||", 1)
                own_meta = seq_ownership.get(sseqid)
                if own_meta is None:
                    if "|" in sseqid and not sseqid.startswith("refseq_"):
                        rxn_id_fallback = sseqid.split("|", 1)[0]
                        own_meta = {
                            "dedup_seq_id": sseqid,
                            "seq_sha256": "",
                            "clean_rxn_ids": [rxn_id_fallback],
                            "is_shared_across_candidate_rxns": False,
                            "source_relpaths_by_rxn": {rxn_id_fallback: []},
                            "primary_raw_id": sseqid,
                            "primary_raw_header": sseqid,
                        }
                    else:
                        continue

                bs = float(bitscore)
                qs, qe, ql = int(qstart), int(qend), int(qlen)
                ss, se, sl = int(sstart), int(send), int(slen)
                aln_len = int(length)

                lo, hi = min(qs, qe), max(qs, qe)
                dist_edge = min(lo - 1, max(0, ql - hi))
                is_near_edge = dist_edge <= 300
                single_hsp_ref_cov = compute_validated_single_hsp_reference_coverage(ss, se, sl)
                is_truncated_near_edge = is_near_edge and (single_hsp_ref_cov < 0.70)
                is_broad_alignment_near_edge = is_near_edge and (single_hsp_ref_cov >= 0.70)

                for rxn_id in own_meta["clean_rxn_ids"]:
                    key = (g_name, rxn_id)
                    prev = best_hits.get(key)
                    if prev is None or bs > prev["bitscore"]:
                        src_files = own_meta["source_relpaths_by_rxn"].get(rxn_id, [])
                        best_hits[key] = {
                            "contig_id": cid,
                            "bitscore": bs,
                            "pident": float(pident),
                            "evalue": float(evalue),
                            "alignment_length": aln_len,
                            "subject_length": sl,
                            "single_hsp_reference_coverage": single_hsp_ref_cov,
                            "subject_coverage": single_hsp_ref_cov,
                            "qstart": qs,
                            "qend": qe,
                            "qlen": ql,
                            "sstart": ss,
                            "send": se,
                            "dist_to_contig_edge_bp": dist_edge,
                            "is_near_contig_edge": is_near_edge,
                            "is_truncated_near_edge": is_truncated_near_edge,
                            "is_broad_alignment_near_edge": is_broad_alignment_near_edge,
                            "best_hit_dedup_seq_id": own_meta["dedup_seq_id"],
                            "best_hit_seq_sha256": own_meta["seq_sha256"],
                            "best_hit_ref_files_csv": ",".join(src_files),
                            "best_hit_ref_header": own_meta["primary_raw_header"],
                            "best_hit_ref_is_shared_across_candidate_rxns": own_meta[
                                "is_shared_across_candidate_rxns"
                            ],
                            "best_hit_ref_shared_rxn_ids_csv": ",".join(
                                own_meta["clean_rxn_ids"]
                            ),
                        }

    tsv_hash = compute_file_sha256(out_tsv)
    mag_audit = {
        "n_archive_mags_total": len(archive_genomes),
        "n_candidate_mags_verified_in_archive": len(set(candidate_genomes)),
        "n_raw_diamond_alignments": n_raw_alignments,
        "diamond_command": cmd,
    }
    return best_hits, contig_lengths, tsv_hash, mag_audit


def search_sample_assemblies_for_community_homologs(
    samples: list[str],
    assemblies_dir: Path,
    db_dmnd: Path,
    seq_ownership: dict[str, dict[str, Any]],
    raw_searches_dir: Path,
    tmp_dir: Path,
    threads: int,
) -> tuple[dict[tuple[str, str], dict[str, Any]], dict[str, bool], dict[str, str]]:
    """Search unbinned sample assemblies against deduplicated target reactions database.
    Uses 12-column outfmt (including `qstart/qend/qlen` and `sstart/send/slen`) and computes
    validated `single_hsp_reference_coverage = (abs(send - sstart) + 1) / slen`.
    """
    if not assemblies_dir.exists() or not assemblies_dir.is_dir():
        raise FileNotFoundError(f"Sample assemblies directory does not exist: {assemblies_dir}")

    sample_hits: dict[tuple[str, str], dict[str, Any]] = {}
    sample_availability: dict[str, bool] = {}
    sample_tsv_hashes: dict[str, str] = {}

    unique_samples = sorted(set(samples))
    existing_samples: list[tuple[str, Path]] = []
    for sample in unique_samples:
        contigs_fa = assemblies_dir / sample / f"{sample}_contigs.fa"
        if not contigs_fa.exists():
            sample_availability[sample] = False
        else:
            sample_availability[sample] = True
            existing_samples.append((sample, contigs_fa))

    if not existing_samples:
        return sample_hits, sample_availability, sample_tsv_hashes

    max_workers = min(4, len(existing_samples), max(1, threads))
    threads_per_worker = max(1, threads // max_workers)

    def _run_one_sample(sample: str, contigs_fa: Path) -> tuple[str, Path, str]:
        out_tsv = raw_searches_dir / f"sample_{sample}_vs_target_rxns.tsv"
        sample_tmp = tmp_dir / f"tmp_sample_{sample}"
        sample_tmp.mkdir(parents=True, exist_ok=True)
        cmd = [
            "diamond",
            "blastx",
            "--query",
            str(contigs_fa),
            "--db",
            str(db_dmnd),
            "--out",
            str(out_tsv),
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qstart",
            "qend",
            "qlen",
            "sstart",
            "send",
            "slen",
            "evalue",
            "bitscore",
            "--evalue",
            "1e-5",
            "--max-target-seqs",
            "0",
            "--threads",
            str(threads_per_worker),
            "--tmpdir",
            str(sample_tmp),
        ]
        subprocess.run(cmd, check=True)
        return sample, out_tsv, compute_file_sha256(out_tsv)

    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {
            ex.submit(_run_one_sample, sample, contigs_fa): sample
            for sample, contigs_fa in existing_samples
        }
        for fut in as_completed(futs):
            sample, out_tsv, sha_hex = fut.result()
            sample_tsv_hashes[sample] = sha_hex
            if out_tsv.exists():
                with open(out_tsv, "r", encoding="utf-8") as fh:
                    for line in fh:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 12:
                            continue
                        (
                            qseqid,
                            sseqid,
                            pident,
                            aln_len,
                            qstart,
                            qend,
                            qlen,
                            sstart,
                            send,
                            slen,
                            evalue,
                            bitscore,
                        ) = parts
                        own_meta = seq_ownership.get(sseqid)
                        if own_meta is None:
                            if "|" in sseqid and not sseqid.startswith("refseq_"):
                                rxn_id_fallback = sseqid.split("|", 1)[0]
                                own_meta = {
                                    "dedup_seq_id": sseqid,
                                    "seq_sha256": "",
                                    "clean_rxn_ids": [rxn_id_fallback],
                                    "is_shared_across_candidate_rxns": False,
                                    "source_relpaths_by_rxn": {rxn_id_fallback: []},
                                    "primary_raw_id": sseqid,
                                    "primary_raw_header": sseqid,
                                }
                            else:
                                continue
                        bs = float(bitscore)
                        qs, qe, ql = int(qstart), int(qend), int(qlen)
                        ss, se, sl = int(sstart), int(send), int(slen)
                        single_hsp_ref_cov = compute_validated_single_hsp_reference_coverage(
                            ss, se, sl
                        )
                        for rxn_id in own_meta["clean_rxn_ids"]:
                            key = (sample, rxn_id)
                            prev = sample_hits.get(key)
                            if prev is None or bs > prev["sample_assembly_bitscore"]:
                                sample_hits[key] = {
                                    "sample_contig_id": qseqid,
                                    "sample_assembly_bitscore": bs,
                                    "sample_assembly_pident": float(pident),
                                    "sample_assembly_evalue": float(evalue),
                                    "sample_alignment_length": int(aln_len),
                                    "sample_subject_length": sl,
                                    "sample_single_hsp_reference_coverage": single_hsp_ref_cov,
                                    "sample_qstart": qs,
                                    "sample_qend": qe,
                                    "sample_qlen": ql,
                                    "sample_sstart": ss,
                                    "sample_send": se,
                                    "sample_best_hit_dedup_seq_id": own_meta["dedup_seq_id"],
                                    "sample_best_hit_seq_sha256": own_meta["seq_sha256"],
                                    "sample_best_hit_ref_header": own_meta["primary_raw_header"],
                                    "sample_best_hit_ref_is_shared": own_meta[
                                        "is_shared_across_candidate_rxns"
                                    ],
                                }

    return sample_hits, sample_availability, sample_tsv_hashes


def _clean_str_or_empty(val: Any) -> str:
    if val is None or (isinstance(val, float) and math.isnan(val)) or pd.isna(val):
        return ""
    s = str(val).strip()
    return "" if s.lower() == "nan" else s


def verify_candidate_genomes_against_inventory(
    inst_df: pd.DataFrame,
    inventory_tsv: Path | None,
) -> dict[str, Any]:
    """Verify candidate genome IDs, cohort membership (`EQ_desert`), and canonical taxonomy
    against `exp1_model_inventory_and_qc_484.tsv`.
    """
    cand_genomes = sorted(inst_df["genome"].dropna().astype(str).unique().tolist())
    if inventory_tsv is None:
        return {
            "inventory_verified": False,
            "n_candidate_genomes": len(cand_genomes),
        }
    if not inventory_tsv.exists() or not inventory_tsv.is_file():
        raise FileNotFoundError(f"Genome inventory TSV does not exist: {inventory_tsv}")

    inv_df = pd.read_csv(inventory_tsv, sep="\t")
    inv_by_g = {str(r["genome"]): r for _, r in inv_df.iterrows()}
    missing_g = [g for g in cand_genomes if g not in inv_by_g]
    if missing_g:
        raise ValueError(
            f"{len(missing_g)} candidate genomes missing from inventory {inventory_tsv}: {missing_g[:5]}"
        )

    non_eq = [g for g in cand_genomes if _clean_str_or_empty(inv_by_g[g].get("cohort")) != "EQ_desert"]
    if non_eq:
        raise ValueError(
            f"Candidate deficiency genomes must belong to EQ_desert cohort; found non-EQ: {non_eq[:5]}"
        )

    taxonomy_mismatches: list[str] = []
    for _, row in inst_df[["genome", "canonical_order", "canonical_family"]].drop_duplicates().iterrows():
        g = str(row["genome"])
        inv_r = inv_by_g[g]
        if _clean_str_or_empty(row.get("canonical_order")) != _clean_str_or_empty(
            inv_r.get("canonical_order")
        ):
            taxonomy_mismatches.append(f"{g}:order")
        if _clean_str_or_empty(row.get("canonical_family")) != _clean_str_or_empty(
            inv_r.get("canonical_family")
        ):
            taxonomy_mismatches.append(f"{g}:family")
    if taxonomy_mismatches:
        raise ValueError(
            f"Candidate genome taxonomy mismatch against inventory: {taxonomy_mismatches[:5]}"
        )

    return {
        "inventory_verified": True,
        "inventory_tsv": str(inventory_tsv),
        "inventory_sha256": compute_file_sha256(inventory_tsv),
        "n_candidate_genomes": len(cand_genomes),
        "all_candidate_genomes_in_eq_desert_484": True,
        "n_unique_orders": int(inst_df["canonical_order"].nunique()),
        "n_unique_families": int(inst_df["canonical_family"].nunique()),
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Stage 4 (v5 Repaired): Candidate-Conditioned Deficiency Verification."
    )
    ap.add_argument("--instances-tsv", type=Path, required=True)
    ap.add_argument("--summary-tsv", type=Path, required=True)
    ap.add_argument(
        "--inventory-tsv",
        type=Path,
        default=None,
        help="Optional path to exp1_model_inventory_and_qc_484.tsv for taxonomy/cohort verification.",
    )
    ap.add_argument(
        "--gapseq-seq-archaea-dir",
        type=Path,
        default=Path("/project/k10018/sw/db/gapseq_seq/Archaea"),
        help="Root of gapseq Archaea sequence DB containing rev/, unrev/, rxn/.",
    )
    ap.add_argument(
        "--gapseq-rxn-fasta-dir",
        type=Path,
        default=None,
        help="Legacy flag for backward compatibility; if pointing to .../Archaea/rxn, parent .../Archaea is used.",
    )
    ap.add_argument(
        "--gapseq-dat-dir",
        type=Path,
        default=Path("/project/k10018/sw/envs/gapseq/share/gapseq/dat"),
        help="Installed gapseq dat/ directory containing authoritative reaction crosswalk tables.",
    )
    ap.add_argument(
        "--gapseq-user-seq-dir",
        type=Path,
        default=Path("/project/k10018/sw/envs/gapseq/share/gapseq/dat/seq/Archaea/user"),
    )
    ap.add_argument(
        "--gapseq-out-dir",
        type=Path,
        default=Path("/project/k10018/work/eq-annot-20261030/results/c11_gapseq/out"),
        help="Directory of per-genome <g>.tar.gz archives containing <g>/<g>-all-Reactions.tbl.",
    )
    ap.add_argument(
        "--refined-mags-tar",
        type=Path,
        default=Path("/scratch/project/k10018/eq-science-20261108/inputs/thermo237_refined_mags.tar.gz"),
    )
    ap.add_argument(
        "--assemblies-dir",
        type=Path,
        default=Path("/project/k10018/empty-quarter/data/assemblies"),
    )
    ap.add_argument(
        "--gap-context-tsv-gz",
        type=Path,
        default=Path("/scratch/project/k10018/eq-archaeal-function-20261106/results/gap_context_candidates.tsv.gz"),
    )
    ap.add_argument(
        "--c18-candidates-tsv",
        type=Path,
        default=Path("/project/k10018/work/eq-annot-20261030/results/c18_pathway_holes/candidates.tsv"),
    )
    ap.add_argument(
        "--c18-windows-tsv",
        type=Path,
        default=Path("/project/k10018/work/eq-annot-20261109-nonexistent"),
        help="Optional path to c18_pathway_holes/windows.tsv.",
    )
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--tmp-dir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument(
        "--output-version",
        type=str,
        default="",
        help="Optional version tag (e.g. 'v5') to emit versioned filenames without overwriting 15894447 baseline.",
    )
    ap.add_argument(
        "--allow-all-unmapped-in-unit-test",
        action="store_true",
        help="Only for unit tests specifically testing single-reaction unmapped disposition.",
    )
    args = ap.parse_args()

    t0 = time.time()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.tmp_dir.mkdir(parents=True, exist_ok=True)
    prov_dir = args.out_dir / "provenance"
    prov_dir.mkdir(parents=True, exist_ok=True)

    ver_prefix = f"_{args.output_version}" if args.output_version else ""
    raw_searches_dir = prov_dir / f"raw_searches{ver_prefix}"
    raw_searches_dir.mkdir(parents=True, exist_ok=True)

    inst_out_path = args.out_dir / f"exp1_stage4{ver_prefix}_deficiency_verification_instances.tsv"
    sum_out_path = args.out_dir / f"exp1_stage4{ver_prefix}_deficiency_verification_summary.tsv"
    receipt_out_path = prov_dir / f"stage4{ver_prefix}_execution_receipts.json"
    crosswalk_out_path = prov_dir / f"stage4{ver_prefix}_reference_namespace_crosswalk.tsv"
    ownership_out_path = prov_dir / f"stage4{ver_prefix}_reference_sequence_ownership.tsv"

    if not args.instances_tsv.exists():
        raise FileNotFoundError(f"Instances TSV does not exist: {args.instances_tsv}")
    if not args.summary_tsv.exists():
        raise FileNotFoundError(f"Summary TSV does not exist: {args.summary_tsv}")

    inst_df = (
        pd.read_csv(args.instances_tsv, sep="\t")
        if args.instances_tsv.stat().st_size > 10
        else pd.DataFrame()
    )
    collector_sum_df = (
        pd.read_csv(args.summary_tsv, sep="\t")
        if args.summary_tsv.stat().st_size > 10
        else pd.DataFrame()
    )

    # Empty candidate handling per Root Item B: emit schema-valid completed nonestimable outputs
    if inst_df.empty or len(inst_df) == 0:
        empty_inst_cols = [
            "scenario_id", "rxn_id", "base_id", "rxn_name", "is_transport",
            "biomass_sensitivity_classification", "followup_selection_tier", "epistemic_status",
            "genome", "cohort", "canonical_order", "canonical_family", "canonical_genus",
            "completeness", "contamination", "sample_id", "clean_rxn_id",
            "search_execution_status", "reference_sequence_status", "reference_resolution_status",
            "is_ambiguous_many_to_many_reference", "shared_with_candidate_rxn_ids_csv",
            "n_reference_sequences", "mag_search_status", "mag_tblastn_bitscore",
            "mag_tblastn_contig_id", "mag_single_hsp_reference_coverage",
            "mag_tblastn_subject_coverage", "mag_alignment_length", "mag_dist_to_contig_edge_bp",
            "mag_hit_near_contig_edge", "mag_hit_truncated_near_edge",
            "mag_gene_call_overlap_status", "sample_assembly_status",
            "sample_assembly_search_status", "sample_assembly_bitscore",
            "sample_assembly_contig_id", "sample_single_hsp_reference_coverage",
            "sample_alignment_length", "verified_host_linkage", "n_synteny_dark_candidates",
            "synteny_dark_candidates_csv", "reference_ambiguity_qualification",
            "deficiency_verification_disposition",
        ]
        pd.DataFrame(columns=empty_inst_cols).to_csv(inst_out_path, sep="\t", index=False)
        pd.DataFrame(
            columns=["clean_rxn_id", "rxn_name", "deficiency_verification_disposition", "n_instances"]
        ).to_csv(sum_out_path, sep="\t", index=False)
        with open(receipt_out_path, "w", encoding="utf-8") as fh:
            json.dump(
                {
                    "status": "completed_zero_candidates",
                    "n_candidates": 0,
                    "elapsed_seconds": round(time.time() - t0, 2),
                },
                fh,
                indent=2,
            )
        print("[Stage 4 Complete] Zero candidates to verify; emitted schema-valid empty outputs.")
        return

    seq_archaea_dir = args.gapseq_seq_archaea_dir
    if args.gapseq_rxn_fasta_dir is not None:
        if (
            args.gapseq_rxn_fasta_dir.name == "rxn"
            and (args.gapseq_rxn_fasta_dir.parent / "rev").is_dir()
        ):
            seq_archaea_dir = args.gapseq_rxn_fasta_dir.parent
        elif (args.gapseq_rxn_fasta_dir / "rev").is_dir():
            seq_archaea_dir = args.gapseq_rxn_fasta_dir
        else:
            seq_archaea_dir = args.gapseq_rxn_fasta_dir

    inv_audit = verify_candidate_genomes_against_inventory(inst_df, args.inventory_tsv)

    c18_win = args.c18_windows_tsv if (args.c18_windows_tsv and args.c18_windows_tsv.exists()) else None
    synteny_index = load_pathway_hole_synteny_candidates(
        gap_context_tsv_gz=args.gap_context_tsv_gz,
        c18_candidates_tsv=args.c18_candidates_tsv,
        c18_windows_tsv=c18_win,
    )

    candidate_genomes = sorted(inst_df["genome"].dropna().astype(str).unique().tolist())
    rxn_meta_by_clean_id: dict[str, dict[str, Any]] = {}
    for _, row in inst_df.iterrows():
        bid = str(row["base_id"])
        clean_rxn = re.sub(r"_[ce]\d*$", "", bid)
        if clean_rxn not in rxn_meta_by_clean_id:
            rxn_meta_by_clean_id[clean_rxn] = {
                "rxn_name": str(row.get("rxn_name", "")),
                "is_transport": str(row.get("is_transport", "False")).lower() == "true",
            }
    target_clean_rxns = sorted(rxn_meta_by_clean_id.keys())

    gapseq_out = args.gapseq_out_dir if (args.gapseq_out_dir and args.gapseq_out_dir.exists()) else None
    gapseq_user = args.gapseq_user_seq_dir if (args.gapseq_user_seq_dir and args.gapseq_user_seq_dir.exists()) else None

    ref_fasta = args.tmp_dir / f"stage4{ver_prefix}_dedup_target_rxns.faa"
    (
        n_dedup_seqs,
        ref_inventory,
        seq_ownership,
        crosswalk_df,
        ownership_df,
        crosswalk_source_hashes,
    ) = resolve_gapseq_reference_crosswalk(
        target_clean_rxns=target_clean_rxns,
        rxn_meta_by_clean_id=rxn_meta_by_clean_id,
        candidate_genomes=candidate_genomes,
        gapseq_dat_dir=args.gapseq_dat_dir,
        gapseq_seq_archaea_dir=seq_archaea_dir,
        gapseq_user_seq_dir=gapseq_user,
        gapseq_out_dir=gapseq_out,
        out_dedup_fasta=ref_fasta,
    )

    crosswalk_df.to_csv(crosswalk_out_path, sep="\t", index=False)
    ownership_df.to_csv(ownership_out_path, sep="\t", index=False)

    n_mapped_rxns = sum(1 for r in target_clean_rxns if ref_inventory[r]["reference_available"])
    n_unmapped_rxns = len(target_clean_rxns) - n_mapped_rxns

    if n_mapped_rxns == 0 and not args.allow_all_unmapped_in_unit_test:
        raise RuntimeError(
            f"FAIL-CLOSED: 0 of {len(target_clean_rxns)} candidate reactions resolved to reference "
            f"sequences under {seq_archaea_dir} and {args.gapseq_dat_dir}. Refusing to emit "
            f"unperformed verification outputs."
        )

    mag_hits: dict[tuple[str, str], dict[str, Any]] = {}
    sample_hits: dict[tuple[str, str], dict[str, Any]] = {}
    sample_availability: dict[str, bool] = {}
    mag_tsv_hash = ""
    sample_tsv_hashes: dict[str, str] = {}
    mag_audit: dict[str, Any] = {}

    samples = sorted({
        str(
            row.get("sample_id")
            or (str(row["genome"]).split("__")[0] if "__" in str(row["genome"]) else str(row["genome"]))
        )
        for _, row in inst_df.iterrows()
    })

    if n_dedup_seqs > 0:
        db_dmnd = args.tmp_dir / f"stage4{ver_prefix}_target_rxns.dmnd"
        subprocess.run(
            [
                "diamond",
                "makedb",
                "--in",
                str(ref_fasta),
                "--db",
                str(db_dmnd),
                "--threads",
                str(args.threads),
            ],
            check=True,
        )
        mag_hits, _, mag_tsv_hash, mag_audit = extract_refined_mags_and_run_diamond_blastx(
            refined_mags_tar=args.refined_mags_tar,
            candidate_genomes=candidate_genomes,
            db_dmnd=db_dmnd,
            seq_ownership=seq_ownership,
            raw_searches_dir=raw_searches_dir,
            tmp_dir=args.tmp_dir,
            threads=args.threads,
        )
        sample_hits, sample_availability, sample_tsv_hashes = (
            search_sample_assemblies_for_community_homologs(
                samples=samples,
                assemblies_dir=args.assemblies_dir,
                db_dmnd=db_dmnd,
                seq_ownership=seq_ownership,
                raw_searches_dir=raw_searches_dir,
                tmp_dir=args.tmp_dir,
                threads=args.threads,
            )
        )

    verified_rows: list[dict[str, Any]] = []
    n_instances_searched = 0
    n_instances_unsearched = 0

    for _, row in inst_df.iterrows():
        g = str(row["genome"])
        sample = str(row.get("sample_id") or (g.split("__")[0] if "__" in g else g))
        bid = str(row["base_id"])
        clean_rxn = re.sub(r"_[ce]\d*$", "", str(bid))
        pwy = _clean_str_or_empty(row.get("pathway_id")).strip("|")
        is_trans = str(row.get("is_transport", "False")).lower() == "true"

        # Axis 1: Reference Availability & Many-to-Many Ambiguity (within selected 63 candidates)
        ref_meta = ref_inventory.get(clean_rxn, {})
        has_ref = bool(ref_meta.get("reference_available", False))
        is_shared_ref = bool(ref_meta.get("is_ambiguous_many_to_many_reference", False))
        shared_rxns_csv = str(ref_meta.get("shared_with_candidate_rxn_ids_csv", ""))
        res_status = str(ref_meta.get("reference_resolution_status", "unresolved"))

        if has_ref:
            n_instances_searched += 1
            search_exec_status = "searched_with_gapseq_archaea_reference_sequences"
            ref_seq_status = (
                "selected_candidate_reference_shared_many_to_many"
                if is_shared_ref
                else "selected_candidate_reference_unique"
            )
        else:
            n_instances_unsearched += 1
            search_exec_status = (
                "not_searched_missing_transporter_reference_sequences"
                if is_trans
                else "not_searched_no_gapseq_archaea_reference_sequences"
            )
            ref_seq_status = "reference_unavailable"

        # Axis 2: Refined MAG Contigs Search (single-HSP reference coverage)
        mh = mag_hits.get((g, clean_rxn), {})
        mag_evaluated = bool(has_ref and args.refined_mags_tar.exists())
        mag_bs = float(mh.get("bitscore", 0.0)) if mag_evaluated else math.nan
        mag_pident = float(mh.get("pident", 0.0)) if mag_evaluated and mh else math.nan
        mag_eval = float(mh.get("evalue", math.nan)) if mag_evaluated and mh else math.nan
        mag_aln_len = int(mh.get("alignment_length", 0)) if mag_evaluated and mh else -1
        mag_slen = int(mh.get("subject_length", 0)) if mag_evaluated and mh else -1
        near_edge = bool(mh.get("is_near_contig_edge", False)) if mag_evaluated else False
        dist_edge = int(mh.get("dist_to_contig_edge_bp", -1)) if mag_evaluated else -1

        raw_mag_cov = mh.get(
            "single_hsp_reference_coverage",
            mh.get("subject_coverage", math.nan),
        )
        mag_ref_cov = float(raw_mag_cov) if (mag_evaluated and mh) else math.nan
        trunc_near_edge = bool(near_edge and not math.isnan(mag_ref_cov) and mag_ref_cov < 0.70)
        broad_near_edge = bool(near_edge and not math.isnan(mag_ref_cov) and mag_ref_cov >= 0.70)
        mag_hit_is_shared = bool(mh.get("best_hit_ref_is_shared_across_candidate_rxns", False))

        # Axis 3: Unbinned Sample Metagenome Assembly Search (with single-HSP reference coverage)
        sample_fa_exists = bool(sample_availability.get(sample, False))
        sample_evaluated = bool(has_ref and sample_fa_exists)
        sh = sample_hits.get((sample, clean_rxn), {}) if sample_evaluated else {}
        samp_bs = float(sh.get("sample_assembly_bitscore", 0.0)) if sample_evaluated else math.nan
        samp_aln_len = int(sh.get("sample_alignment_length", 0)) if sample_evaluated and sh else -1
        samp_slen = int(sh.get("sample_subject_length", 0)) if sample_evaluated and sh else -1
        raw_samp_cov = sh.get("sample_single_hsp_reference_coverage", math.nan)
        samp_ref_cov = float(raw_samp_cov) if (sample_evaluated and sh) else math.nan
        samp_hit_is_shared = bool(sh.get("sample_best_hit_ref_is_shared", False))
        verified_host_linkage = False  # Contig provenance unverified from raw sample assembly

        # Axis 4: Synteny & Pathway-Hole Genomic Context
        synt_ctx = lookup_synteny_context_for_instance(
            genome=g,
            clean_rxn=clean_rxn,
            instance_pathway_id=pwy,
            ref_meta=ref_meta,
            synteny_index=synteny_index,
        )
        synt_genes = synt_ctx["synteny_dark_genes"]
        rxn_spec_dark = synt_ctx["rxn_specific_dark_genes"]
        c18_hole_dark = synt_ctx["c18_hole_dark_genes"]

        # Axis 5: Qualified Competing-Hypothesis Disposition (`EXP1-DEFICIENCY-ROOT-SEMANTIC-REVIEW.txt`)
        mag_gene_call_overlap_status = "not_applicable_no_mag_hit"
        if not has_ref:
            if is_trans:
                disposition = "unknown_missing_transporter_reference_sequences"
            else:
                disposition = "not_evaluable_missing_reference_sequences"
            ref_ambig_qual = "not_searched_no_reference"
        elif not args.refined_mags_tar.exists():
            disposition = "not_evaluable_missing_mag_assembly"
            ref_ambig_qual = "not_searched_missing_mag_archive"
        elif not math.isnan(mag_bs) and mag_bs >= 35.0:
            mag_gene_call_overlap_status = "unknown_not_checked_against_gene_coordinates"
            if trunc_near_edge:
                disposition = "candidate_truncated_or_split_gene_at_contig_edge"
            elif broad_near_edge:
                disposition = "broad_reference_alignment_near_contig_edge"
            elif (
                mag_bs >= 50.0
                and mh.get("pident", 0.0) >= 40.0
                and not math.isnan(mag_ref_cov)
                and mag_ref_cov >= 0.70
            ):
                disposition = "homolog_candidate_detected_in_MAG"
            else:
                disposition = "domain_level_or_partial_homolog_candidate_in_MAG"
            ref_ambig_qual = (
                "mag_hit_on_shared_many_to_many_reference_not_reaction_unique"
                if (mag_hit_is_shared or is_shared_ref)
                else "mag_hit_on_selected_candidate_reference_unique_no_catalytic_specificity_proof"
            )
        elif not math.isnan(samp_bs) and samp_bs >= 65.0:
            if not math.isnan(samp_ref_cov) and samp_ref_cov < 0.70:
                disposition = "unresolved_sample_community_domain_level_hit_host_unassigned"
            elif verified_host_linkage:
                disposition = "demonstrated_assembly_or_binning_dropout_with_host_linkage"
            else:
                disposition = "unresolved_sample_community_homolog_host_unassigned"
            ref_ambig_qual = (
                "sample_hit_on_shared_many_to_many_reference_not_reaction_unique"
                if (samp_hit_is_shared or is_shared_ref)
                else "sample_hit_on_selected_candidate_reference_unique_no_catalytic_specificity_proof"
            )
        elif len(synt_genes) > 0:
            disposition = "candidate_novel_or_non_orthologous_replacement_in_synteny"
            ref_ambig_qual = (
                "tested_no_homolog_hit_shared_reference_set"
                if is_shared_ref
                else "tested_no_homolog_hit_selected_candidate_reference_unique"
            )
        elif not sample_fa_exists:
            disposition = "tested_no_homolog_detected_in_mag_sample_assembly_unevaluated"
            ref_ambig_qual = (
                "tested_no_homolog_hit_shared_reference_set"
                if is_shared_ref
                else "tested_no_homolog_hit_selected_candidate_reference_unique"
            )
        else:
            disposition = "tested_no_homolog_detected_conditional_negative_evidence"
            ref_ambig_qual = (
                "tested_no_homolog_hit_shared_reference_set"
                if is_shared_ref
                else "tested_no_homolog_hit_selected_candidate_reference_unique"
            )

        out_r = dict(row)
        out_r["sample_id"] = sample
        out_r["clean_rxn_id"] = clean_rxn
        out_r["search_execution_status"] = search_exec_status
        out_r["reference_sequence_status"] = ref_seq_status
        out_r["reference_resolution_status"] = res_status
        out_r["is_ambiguous_many_to_many_reference"] = is_shared_ref
        out_r["shared_with_candidate_rxn_ids_csv"] = shared_rxns_csv
        out_r["mapped_metacyc_rxn_ids_csv"] = ",".join(ref_meta.get("mapped_metacyc_rxn_ids", []))
        out_r["mapped_ec_numbers_csv"] = ",".join(ref_meta.get("mapped_ec_numbers", []))
        out_r["mapped_fasta_relpaths_csv"] = ",".join(ref_meta.get("mapped_fasta_relpaths", []))
        out_r["n_reference_sequences"] = int(ref_meta.get("n_reference_sequences", 0))
        out_r["mag_search_status"] = "evaluated" if mag_evaluated else "not_searched_missing_reference"
        out_r["mag_tblastn_bitscore"] = mag_bs
        out_r["mag_tblastn_pident"] = mag_pident
        out_r["mag_tblastn_evalue"] = mag_eval
        out_r["mag_tblastn_contig_id"] = mh.get("contig_id", "")
        out_r["mag_alignment_length"] = mag_aln_len
        out_r["mag_subject_length"] = mag_slen
        out_r["mag_single_hsp_reference_coverage"] = mag_ref_cov
        out_r["mag_tblastn_subject_coverage"] = mag_ref_cov
        out_r["mag_qstart"] = int(mh.get("qstart", -1)) if mag_evaluated and mh else -1
        out_r["mag_qend"] = int(mh.get("qend", -1)) if mag_evaluated and mh else -1
        out_r["mag_qlen_contig_bp"] = int(mh.get("qlen", -1)) if mag_evaluated and mh else -1
        out_r["mag_sstart"] = int(mh.get("sstart", -1)) if mag_evaluated and mh else -1
        out_r["mag_send"] = int(mh.get("send", -1)) if mag_evaluated and mh else -1
        out_r["mag_dist_to_contig_edge_bp"] = dist_edge
        out_r["mag_hit_near_contig_edge"] = near_edge
        out_r["mag_hit_truncated_near_edge"] = trunc_near_edge
        out_r["mag_hit_broad_alignment_near_edge"] = broad_near_edge
        out_r["mag_gene_call_overlap_status"] = mag_gene_call_overlap_status
        out_r["mag_best_hit_dedup_seq_id"] = mh.get("best_hit_dedup_seq_id", "")
        out_r["mag_best_hit_seq_sha256"] = mh.get("best_hit_seq_sha256", "")
        out_r["mag_best_hit_ref_files_csv"] = mh.get("best_hit_ref_files_csv", "")
        out_r["mag_best_hit_ref_header"] = mh.get("best_hit_ref_header", "")
        out_r["mag_best_hit_ref_is_shared"] = mag_hit_is_shared
        out_r["sample_assembly_status"] = (
            "sample_assembly_available"
            if sample_fa_exists
            else "unknown_sample_assembly_unavailable"
        )
        out_r["sample_assembly_search_status"] = (
            "evaluated"
            if sample_evaluated
            else (
                "not_searched_missing_sample_assembly"
                if has_ref
                else "not_searched_missing_reference"
            )
        )
        out_r["sample_assembly_bitscore"] = samp_bs
        out_r["sample_assembly_contig_id"] = sh.get("sample_contig_id", "")
        out_r["sample_alignment_length"] = samp_aln_len
        out_r["sample_subject_length"] = samp_slen
        out_r["sample_single_hsp_reference_coverage"] = samp_ref_cov
        out_r["sample_qstart"] = int(sh.get("sample_qstart", -1)) if sample_evaluated and sh else -1
        out_r["sample_qend"] = int(sh.get("sample_qend", -1)) if sample_evaluated and sh else -1
        out_r["sample_qlen_contig_bp"] = int(sh.get("sample_qlen", -1)) if sample_evaluated and sh else -1
        out_r["sample_sstart"] = int(sh.get("sample_sstart", -1)) if sample_evaluated and sh else -1
        out_r["sample_send"] = int(sh.get("sample_send", -1)) if sample_evaluated and sh else -1
        out_r["sample_best_hit_dedup_seq_id"] = sh.get("sample_best_hit_dedup_seq_id", "")
        out_r["sample_best_hit_seq_sha256"] = sh.get("sample_best_hit_seq_sha256", "")
        out_r["sample_best_hit_ref_header"] = sh.get("sample_best_hit_ref_header", "")
        out_r["verified_host_linkage"] = verified_host_linkage
        out_r["n_c18_hole_dark_candidates"] = len(c18_hole_dark)
        out_r["c18_hole_dark_candidates_csv"] = ",".join(c18_hole_dark[:5])
        out_r["n_rxn_specific_synteny_dark_candidates"] = len(rxn_spec_dark)
        out_r["rxn_specific_synteny_dark_candidates_csv"] = ",".join(rxn_spec_dark[:5])
        out_r["n_synteny_dark_candidates"] = len(synt_genes)
        out_r["synteny_dark_candidates_csv"] = ",".join(synt_genes[:5])
        out_r["n_window_dark_families"] = len(synt_ctx["window_dark_families"])
        out_r["window_dark_families_csv"] = ",".join(synt_ctx["window_dark_families"][:5])
        out_r["synteny_evidence_types_csv"] = synt_ctx["synteny_evidence_types_csv"]
        out_r["reference_ambiguity_qualification"] = ref_ambig_qual
        out_r["deficiency_verification_disposition"] = disposition
        verified_rows.append(out_r)

    ver_df = pd.DataFrame(verified_rows)
    ver_df.to_csv(inst_out_path, sep="\t", index=False)

    if not ver_df.empty:
        ver_sum = (
            ver_df.groupby(
                [
                    "clean_rxn_id",
                    "rxn_name",
                    "reference_resolution_status",
                    "is_ambiguous_many_to_many_reference",
                    "shared_with_candidate_rxn_ids_csv",
                    "search_execution_status",
                    "biomass_sensitivity_classification",
                    "deficiency_verification_disposition",
                ],
                dropna=False,
            )
            .agg(
                n_instances=("genome", "count"),
                n_unique_genomes=("genome", "nunique"),
                n_unique_families=("canonical_family", "nunique")
                if "canonical_family" in ver_df.columns
                else ("genome", "count"),
                n_reference_sequences=("n_reference_sequences", "first"),
                median_mag_bitscore=("mag_tblastn_bitscore", "median"),
                median_mag_single_hsp_ref_coverage=("mag_single_hsp_reference_coverage", "median"),
                median_sample_assembly_bitscore=("sample_assembly_bitscore", "median"),
                median_sample_single_hsp_ref_coverage=(
                    "sample_single_hsp_reference_coverage",
                    "median",
                ),
                n_instances_with_c18_hole_dark=(
                    "n_c18_hole_dark_candidates",
                    lambda s: int((s > 0).sum()),
                ),
                n_instances_with_any_synteny_dark=(
                    "n_synteny_dark_candidates",
                    lambda s: int((s > 0).sum()),
                ),
            )
            .reset_index()
            .sort_values(["n_unique_genomes", "n_instances"], ascending=[False, False])
        )
    else:
        ver_sum = pd.DataFrame()
    ver_sum.to_csv(sum_out_path, sep="\t", index=False)

    disposition_counts = (
        ver_df["deficiency_verification_disposition"].value_counts().to_dict()
        if not ver_df.empty
        else {}
    )
    rxn_res_counts = (
        crosswalk_df["reference_resolution_status"].value_counts().to_dict()
        if not crosswalk_df.empty
        else {}
    )

    receipts = {
        "status": "completed_verified_search" if n_mapped_rxns > 0 else "completed_unit_test_unmapped",
        "output_version": args.output_version or "unversioned",
        "script_path": str(Path(__file__).resolve()),
        "script_sha256": compute_file_sha256(Path(__file__).resolve()),
        "evidence_semantics_notes": {
            "single_hsp_reference_coverage": (
                "Computed as (abs(send - sstart) + 1) / slen for the single best-bitscore HSP, "
                "strictly bounded in (0, 1]. Gapped alignment_length is recorded separately and "
                "never divided by slen. Query is a DNA contig (qlen in bp), so coverage is "
                "reference protein span coverage, never query protein coverage."
            ),
            "selected_candidate_reference_unique": (
                "Denotes that mapped reference sequences are not shared with another reaction "
                "among the 63 selected candidate reactions. Broad EC/MetaCyc crosswalks may "
                "include proteins with multiple substrate activities; this does not prove "
                "reaction-specific catalysis."
            ),
            "homolog_candidate_detected_in_MAG": (
                "Interior contig hit with bitscore >= 50, pident >= 40, and single-HSP reference "
                "coverage >= 0.70. Uncalled vs annotation-miss status is unknown because hits "
                "are not intersected against called ORF coordinates."
            ),
        },
        "input_files": {
            "instances_tsv": str(args.instances_tsv),
            "instances_tsv_sha256": compute_file_sha256(args.instances_tsv),
            "collector_summary_tsv": str(args.summary_tsv),
            "collector_summary_tsv_sha256": compute_file_sha256(args.summary_tsv),
            "refined_mags_tar": str(args.refined_mags_tar),
            "refined_mags_tar_sha256": compute_file_sha256(args.refined_mags_tar),
            "gap_context_tsv_gz": str(args.gap_context_tsv_gz),
            "gap_context_tsv_gz_sha256": compute_file_sha256(args.gap_context_tsv_gz),
            "c18_candidates_tsv": str(args.c18_candidates_tsv),
            "c18_candidates_tsv_sha256": compute_file_sha256(args.c18_candidates_tsv),
            "c18_windows_tsv": str(c18_win) if c18_win else None,
            "c18_windows_tsv_sha256": compute_file_sha256(c18_win) if c18_win else "",
        },
        "gapseq_crosswalk_source_hashes": crosswalk_source_hashes,
        "target_and_query_universe_accounting": {
            "n_candidate_instance_rows": len(ver_df),
            "n_collector_summary_strata_rows": len(collector_sum_df),
            "collector_summary_rows_explanation": (
                "244 rows in exp1_stage4_supported_candidate_deficiencies_summary.tsv represent "
                "(scenario_id, rxn_id, biomass_sensitivity_classification, followup_selection_tier) "
                "strata across 63 distinct ModelSEED reactions, NOT 244 distinct reactions."
            ),
            "n_distinct_candidate_reactions": len(target_clean_rxns),
            "n_distinct_candidate_genomes": len(candidate_genomes),
            "n_distinct_candidate_samples": len(samples),
            "n_sample_assemblies_available": sum(1 for v in sample_availability.values() if v),
            "inventory_verification": inv_audit,
            "mag_archive_verification": mag_audit,
        },
        "reference_namespace_resolution_accounting": {
            "n_reactions_with_reference_sequences_searched": n_mapped_rxns,
            "n_reactions_without_reference_sequences_unsearched": n_unmapped_rxns,
            "reaction_resolution_status_counts": {str(k): int(v) for k, v in rxn_res_counts.items()},
            "n_deduplicated_reference_sequences_total": n_dedup_seqs,
            "deduplicated_reference_fasta_sha256": compute_file_sha256(ref_fasta),
            "n_instances_searched": n_instances_searched,
            "n_instances_unsearched_missing_reference": n_instances_unsearched,
        },
        "disposition_instance_counts": {str(k): int(v) for k, v in disposition_counts.items()},
        "output_files": {
            "verification_instances_tsv": str(inst_out_path),
            "verification_instances_tsv_sha256": compute_file_sha256(inst_out_path),
            "verification_summary_tsv": str(sum_out_path),
            "verification_summary_tsv_sha256": compute_file_sha256(sum_out_path),
            "reference_namespace_crosswalk_tsv": str(crosswalk_out_path),
            "reference_namespace_crosswalk_tsv_sha256": compute_file_sha256(crosswalk_out_path),
            "reference_sequence_ownership_tsv": str(ownership_out_path),
            "reference_sequence_ownership_tsv_sha256": compute_file_sha256(ownership_out_path),
            "mag_search_raw_tsv_sha256": mag_tsv_hash,
            "sample_search_tsv_hashes": sample_tsv_hashes,
        },
        "ref_inventory": ref_inventory,
        "elapsed_seconds": round(time.time() - t0, 2),
    }
    with open(receipt_out_path, "w", encoding="utf-8") as fh:
        json.dump(receipts, fh, indent=2)

    dt = time.time() - t0
    print(
        f"[Stage 4 Complete in {dt:.1f}s] Candidate universe: {len(ver_df)} instances across "
        f"{len(target_clean_rxns)} distinct ModelSEED reactions ({len(collector_sum_df)} collector "
        f"scenario/tier summary strata) in {len(candidate_genomes)} MAGs ({len(samples)} samples). "
        f"Reference resolution: {n_mapped_rxns}/{len(target_clean_rxns)} reactions mapped to gapseq "
        f"Archaea sequences ({n_dedup_seqs} unique sequences); "
        f"{n_unmapped_rxns}/{len(target_clean_rxns)} reactions unsearched (no sequences in gapseq Archaea DB). "
        f"Instances searched: {n_instances_searched}/{len(ver_df)}; "
        f"instances unsearched (missing reference): {n_instances_unsearched}/{len(ver_df)}."
    )


if __name__ == "__main__":
    main()
