#!/usr/bin/env python3
"""Three-paper cross-manuscript consistency and submission-readiness test suite.

Reads the active source trees and included fragments of all three Empty Quarter
manuscripts:
  1. Data resource:    /home/leechuck/Documents/papers/empty-quarter-data-current/sn-article.tex
  2. Ecology paper:    /home/leechuck/Documents/papers/empty-quarter-ecology-current/alt_stats.tex
                       + supp_altstat.tex + included MAG/Life-on-air fragments
  3. Proteome paper:   /home/leechuck/Documents/papers/empty-quarter-proteome/main.tex
                       + supplementary_notes.tex + authors.tex

Modes:
  - Default pytest (`uv run pytest test_three_paper_consistency.py`):
    Runs all cross-paper consistency, source-provenance, cohort-membership,
    taxonomy-reconciliation, citation, and adversarial negative-fixture tests,
    and asserts that the submission-readiness check FAILS (`SubmissionNotReadyError`)
    while pending blockers remain.
  - Submission-readiness CLI gate (`uv run python test_three_paper_consistency.py --require-submission-ready`
    or `EQ_REQUIRE_SUBMISSION_READY=1 uv run pytest test_three_paper_consistency.py`):
    Exits non-zero if any `[PENDING:...]`, `DATA-BULK-ARCHIVE-DOI`, or `\todo{...}`
    placeholder remains in any of the three manuscripts.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

import pandas as pd
import pytest

DATA_MS_DIR = Path(
    os.environ.get("EQ_DATA_MS_DIR", "/home/leechuck/Documents/papers/empty-quarter-data-current")
)
DATA_REPO_DIR = Path(
    os.environ.get("EQ_DATA_REPO_DIR", "/home/leechuck/Public/software/empty-quarter-paper-repos/data-paper")
)
ECO_MS_DIR = Path(
    os.environ.get("EQ_ECO_MS_DIR", "/home/leechuck/Documents/papers/empty-quarter-ecology-current")
)
ECO_REPO_DIR = Path(
    os.environ.get("EQ_ECO_REPO_DIR", "/home/leechuck/Public/software/empty-quarter-paper-repos/ecology")
)
PROT_MS_DIR = Path(
    os.environ.get("EQ_PROT_MS_DIR", "/home/leechuck/Documents/papers/empty-quarter-proteome")
)
PROT_REPO_DIR = Path(
    os.environ.get("EQ_PROT_REPO_DIR", "/home/leechuck/Public/software/empty-quarter-proteome")
)


class SubmissionNotReadyError(AssertionError):
    """Raised when submission-readiness gate encounters unresolved pending markers."""


def strip_tex_comments(text: str) -> str:
    lines = []
    for line in text.splitlines():
        out = []
        i = 0
        while i < len(line):
            if line[i] == "%" and (i == 0 or line[i - 1] != "\\"):
                break
            out.append(line[i])
            i += 1
        lines.append("".join(out))
    return "\n".join(lines)


def expand_tex_tree(entry_tex: Path) -> Tuple[str, List[Path], List[Path]]:
    """Recursively resolve all \\input and \\include files and collect \\includegraphics paths."""
    if not entry_tex.is_file():
        raise AssertionError(f"Missing manuscript entrypoint: {entry_tex}")
    root_dir = entry_tex.parent.resolve()
    visited: List[Path] = []
    figures: List[Path] = []
    chunks: List[str] = []

    def _visit(tex_path: Path):
        resolved = tex_path.resolve()
        if not resolved.is_relative_to(root_dir):
            raise AssertionError(f"External include outside manuscript root: {tex_path}")
        if not resolved.is_file():
            raise AssertionError(f"Missing included TeX file: {tex_path}")
        if resolved in visited:
            return
        visited.append(resolved)
        raw = resolved.read_text(encoding="utf-8")
        uncommented = strip_tex_comments(raw)
        chunks.append(uncommented)

        for fig_rel in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", uncommented):
            fig_path = root_dir / fig_rel
            if not fig_path.is_file():
                # Try common extensions if omitted
                found = None
                for ext in (".pdf", ".png", ".jpg", ".eps"):
                    cand = root_dir / (fig_rel + ext)
                    if cand.is_file():
                        found = cand
                        break
                if found is None:
                    raise AssertionError(f"Missing figure {fig_rel!r} referenced in {resolved}")
                figures.append(found)
            else:
                figures.append(fig_path)

        for inc in re.findall(r"\\(?:input|include)\{([^}]+)\}", uncommented):
            inc_name = inc if inc.endswith(".tex") else f"{inc}.tex"
            _visit(root_dir / inc_name)

    _visit(entry_tex)
    return "\n".join(chunks), visited, figures


def normalize_tex_whitespace(s: str) -> str:
    s = re.sub(r"\\allowbreak\s*", "", s)
    s = re.sub(r"\\(?:textbf|emph|texttt|ttfamily|raggedright)\b", "", s)
    s = s.replace("\\_", "_")
    s = s.replace("{,}", ",")
    s = s.replace("~", " ")
    s = s.replace("\\ ", " ")
    s = s.replace("{", "").replace("}", "")
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def validate_shotgun_cohort_exclusions_in_text(label: str, text: str, manifest: dict) -> None:
    """Validate that shotgun sample exclusions match the foundation manifest and reject fabricated IDs."""
    norm = normalize_tex_whitespace(text)
    fabricated_ids = ["62Sr1", "64Dr1", "64Rr1", "64Sr1", "T345Dr2", "35Rr3", "T160Sr1"]
    for bad_id in fabricated_ids:
        if bad_id in norm:
            raise AssertionError(f"{label} contains fabricated sample exclusion identifier {bad_id!r}")

    cohorts = manifest["cohorts"]
    assert cohorts["EQ-WGS-NONZERO-550"]["excluded_zero_mapped"] == ["F5Dr1"]
    assert cohorts["EQ-WGS-CORE-546"]["excluded_samples"]["auxiliary_location"] == [
        "62PRr2",
        "64Dr3",
        "64PRr2",
        "64Sr2",
    ]
    assert cohorts["EQ-WGS-CORE-MODEL-545"]["excluded_missing_post_qc_reads"] == ["F45PRr3"]
    assert cohorts["EQ-MATCHED-MODEL-441"]["excluded_missing_post_qc_reads"] == ["F45PRr3"]

    if "F5Dr1" not in norm:
        raise AssertionError(f"{label} does not name F5Dr1 as the zero-mapped shotgun profile exclusion")
    for req_id in ["62PRr2", "64Dr3", "64PRr2", "64Sr2", "F45PRr3"]:
        if req_id not in norm:
            raise AssertionError(f"{label} does not name required core-transect/depth exclusion {req_id!r}")


def validate_mimag_boundary_in_text(label: str, text: str) -> None:
    """Validate 126 HQ + 5,311 strict MQ + 1 boundary (48Dr2__self_SemiBin_25 at 64.61% / 10.0%)."""
    norm = normalize_tex_whitespace(text)
    if "48Dr2__self_SemiBin_25" not in norm and "48Dr2\\_\\_self\\_SemiBin\\_25" not in text:
        raise AssertionError(f"{label} does not identify boundary genome 48Dr2__self_SemiBin_25")
    if "90.2\\%" in text or "90.2%" in norm:
        if "48Dr2" in norm:
            raise AssertionError(f"{label} reports fabricated 90.2% completeness for 48Dr2__self_SemiBin_25 (true value: 64.61%)")
    if "64.61" not in norm:
        raise AssertionError(f"{label} does not report true 64.61% completeness for boundary genome 48Dr2__self_SemiBin_25")
    if "5,311" not in norm:
        raise AssertionError(f"{label} does not report 5,311 strict medium-quality MAGs (<10% contamination)")
    if "126" not in norm:
        raise AssertionError(f"{label} does not report 126 high-quality MAGs")


def validate_accessions_and_deposit_honesty(label: str, text: str) -> None:
    """Validate ENA BioProjects, KG Zenodo DOI, and absence of false bulk-archive completion claims."""
    norm = normalize_tex_whitespace(text)
    for acc in ["PRJEB104209", "PRJEB106069", "PRJEB125153"]:
        if acc not in norm:
            raise AssertionError(f"{label} is missing required ENA BioProject accession {acc}")
    if "10.5281/zenodo.23134168" not in norm:
        raise AssertionError(f"{label} is missing released Knowledge Graph v3.0.1 Zenodo DOI 10.5281/zenodo.23134168")
    if "23132097" in norm:
        raise AssertionError(f"{label} contains superseded/wrong Zenodo DOI 10.5281/zenodo.23132097")

    # Reject assigning 1,186 amplicon profiles to PRJEB125153
    if re.search(r"1,186[^.]{0,120}PRJEB125153", norm) and "PRJEB106069" not in re.search(r"1,186[^.]{0,120}PRJEB125153", norm).group(0):
        raise AssertionError(f"{label} wrongly associates 1,186 amplicon profiles with shotgun project PRJEB125153")

    # Reject false assertions that the bulk MAG archive or Bakta protein release is already archived at Zenodo
    forbidden_patterns = [
        r"archived 19,186,792-protein Bakta",
        r"5,438[^.]{0,80}archived at Zenodo",
        r"bulk MAG[^.]{0,60}is archived at Zenodo",
    ]
    for pat in forbidden_patterns:
        if re.search(pat, norm, flags=re.I):
            raise AssertionError(f"{label} makes false completed-archive assertion matching {pat!r}")


def reconcile_canonical_vs_legacy_taxonomy(
    canonical_tsv: Path, legacy_tsv: Path
) -> Dict[str, int]:
    """Compare canonical refined_representatives.tsv against legacy analytical taxonomy 1:1."""
    can = pd.read_csv(canonical_tsv, sep="\t")
    leg = pd.read_csv(legacy_tsv, sep="\t")

    if len(can) != 5438 or len(leg) != 5438:
        raise AssertionError(f"Expected 5,438 rows in both tables, got {len(can)} and {len(leg)}")
    if not can["mag"].is_unique or not leg["genome"].is_unique:
        raise AssertionError("Duplicate genome IDs detected in taxonomy tables")
    if set(can["mag"]) != set(leg["genome"]):
        raise AssertionError("Genome ID sets do not match 1:1 between canonical and legacy tables")

    merged = can.merge(leg, left_on="mag", right_on="genome", how="inner", validate="one_to_one")
    assert len(merged) == 5438

    def split_ranks(s: str) -> Dict[str, str]:
        out = {}
        for part in str(s).split(";"):
            part = part.strip()
            if "__" in part:
                k, v = part.split("__", 1)
                out[k] = v
        return out

    can_ranks = merged["gtdb_r232_classification"].map(split_ranks)
    leg_col = "gtdb_classification" if "gtdb_classification" in merged.columns else "classification"
    leg_ranks = merged[leg_col].map(split_ranks)

    whole_diff = int((merged["gtdb_r232_classification"] != merged[leg_col]).sum())
    rank_diffs = {}
    for r_key, r_name in [("d", "domain"), ("p", "phylum"), ("c", "class"), ("o", "order"), ("f", "family"), ("g", "genus"), ("s", "species")]:
        c_val = can_ranks.map(lambda d: d.get(r_key, ""))
        l_val = leg_ranks.map(lambda d: d.get(r_key, ""))
        rank_diffs[r_name] = int((c_val != l_val).sum())

    can_thermo = set(merged.loc[can_ranks.map(lambda d: d.get("p") == "Thermoplasmatota"), "mag"])
    leg_thermo = set(merged.loc[leg_ranks.map(lambda d: d.get("p") == "Thermoplasmatota"), "mag"])
    if can_thermo != leg_thermo or len(can_thermo) != 237:
        raise AssertionError("Thermoplasmatota 237 genome membership differs between canonical and legacy tables")

    th_mask = merged["mag"].isin(can_thermo)
    th_whole_diff = int((merged.loc[th_mask, "gtdb_r232_classification"] != merged.loc[th_mask, leg_col]).sum())

    th_can_orders = can_ranks[th_mask].map(lambda d: d.get("o", ""))
    th_can_genera = can_ranks[th_mask].map(lambda d: d.get("g", ""))
    th_leg_genera = leg_ranks[th_mask].map(lambda d: d.get("g", ""))

    return {
        "n_genomes": len(merged),
        "whole_diff": whole_diff,
        "domain_diff": rank_diffs["domain"],
        "phylum_diff": rank_diffs["phylum"],
        "class_diff": rank_diffs["class"],
        "order_diff": rank_diffs["order"],
        "genus_diff": rank_diffs["genus"],
        "thermo_n": len(can_thermo),
        "thermo_whole_diff": th_whole_diff,
        "thermo_can_named_orders": int(th_can_orders[th_can_orders != ""].nunique()),
        "thermo_can_unassigned_orders": int((th_can_orders == "").sum()),
        "thermo_can_named_genera": int(th_can_genera[th_can_genera != ""].nunique()),
        "thermo_can_unassigned_genera": int((th_can_genera == "").sum()),
        "thermo_leg_named_genera": int(th_leg_genera[th_leg_genera != ""].nunique()),
        "thermo_leg_unassigned_genera": int((th_leg_genera == "").sum()),
    }


def parse_bib_entry(bib_text: str, key: str) -> Dict[str, str]:
    start_m = re.search(r"@\w+\{\s*" + re.escape(key) + r"\s*,", bib_text)
    if not start_m:
        raise AssertionError(f"Missing BibTeX entry {key!r}")
    pos = start_m.end()
    depth = 1
    while pos < len(bib_text) and depth > 0:
        if bib_text[pos] == "{":
            depth += 1
        elif bib_text[pos] == "}":
            depth -= 1
        pos += 1
    body = bib_text[start_m.end() : pos - 1]
    fields: Dict[str, str] = {}
    i = 0
    while i < len(body):
        km = re.search(r"(\w+)\s*=\s*\{", body[i:])
        if not km:
            break
        k = km.group(1).lower().strip()
        v_start = i + km.end()
        j = v_start
        d = 1
        while j < len(body) and d > 0:
            if body[j] == "{":
                d += 1
            elif body[j] == "}":
                d -= 1
            j += 1
        raw_val = body[v_start : j - 1]
        v = re.sub(r"\s+", " ", raw_val.replace("{", "").replace("}", "")).strip()
        fields[k] = v
        i = j
    return fields


def collect_pending_markers() -> List[Tuple[str, str]]:
    """Return list of (location, marker) across all active manuscript and bibliography files."""
    files = [
        DATA_MS_DIR / "sn-article.tex",
        DATA_MS_DIR / "01_introduction.tex",
        DATA_MS_DIR / "02_methods.tex",
        DATA_MS_DIR / "03_knowledge_representation.tex",
        DATA_MS_DIR / "04_data_records.tex",
        DATA_MS_DIR / "05_validation.tex",
        DATA_MS_DIR / "06_usage.tex",
        DATA_MS_DIR / "sn-bibliography.bib",
        ECO_MS_DIR / "alt_stats.tex",
        ECO_MS_DIR / "supp_altstat.tex",
        ECO_MS_DIR / "mag_lifeonair_results.tex",
        ECO_MS_DIR / "mag_lifeonair_methods.tex",
        ECO_MS_DIR / "mag_lifeonair_supplement.tex",
        ECO_MS_DIR / "sample.bib",
        PROT_MS_DIR / "main.tex",
        PROT_MS_DIR / "supplementary_notes.tex",
        PROT_MS_DIR / "references.bib",
    ]
    markers: List[Tuple[str, str]] = []
    for p in files:
        if not p.is_file():
            markers.append((str(p), "MISSING_ACTIVE_FILE"))
            continue
        txt = strip_tex_comments(p.read_text(encoding="utf-8"))
        for m in re.finditer(r"\[PENDING:[^\]]+\]|DATA-BULK-ARCHIVE-DOI|\\todo\{[^}]+\}", txt):
            clean_m = re.sub(r"\\allowbreak\\?\s*", "", m.group(0))
            clean_m = re.sub(r"\s+", " ", clean_m)
            markers.append((f"{p.parent.name}/{p.name}", clean_m))
    return markers


def assert_submission_ready() -> None:
    """Fail if any pending marker or unresolved archive placeholder remains across the 3 papers."""
    markers = collect_pending_markers()
    if markers:
        formatted = "\n  - ".join(f"{loc}: {m}" for loc, m in markers)
        raise SubmissionNotReadyError(
            f"Submission-readiness gate FAILED with {len(markers)} unresolved pending blocker(s):\n  - {formatted}"
        )


# =============================================================================
# Pytest Test Cases
# =============================================================================


def test_active_entrypoints_and_includes_resolved():
    data_text, data_files, data_figs = expand_tex_tree(DATA_MS_DIR / "sn-article.tex")
    eco_main_text, eco_main_files, eco_main_figs = expand_tex_tree(ECO_MS_DIR / "alt_stats.tex")
    eco_supp_text, eco_supp_files, eco_supp_figs = expand_tex_tree(ECO_MS_DIR / "supp_altstat.tex")
    prot_text, prot_files, prot_figs = expand_tex_tree(PROT_MS_DIR / "main.tex")

    assert len(data_files) >= 7, f"Expected sn-article.tex + >=6 sections, got {len(data_files)}"
    assert (PROT_MS_DIR / "supplementary_notes.tex").resolve() in prot_files
    assert (PROT_MS_DIR / "authors.tex").resolve() in prot_files
    assert len(prot_figs) >= 9, f"Expected 4 main + >=5 supplementary figures in proteome, got {len(prot_figs)}"

    # Ecology must actually include the MAG / Life-on-Air files rather than leaving them orphaned
    eco_all_files = {p.name for p in eco_main_files} | {p.name for p in eco_supp_files}
    assert "mag_lifeonair_methods.tex" in eco_all_files
    assert "mag_lifeonair_supplement.tex" in eco_all_files

    # No local /home filesystem paths in manuscript TeX files
    for label, txt in [
        ("data", data_text),
        ("ecology_main", eco_main_text),
        ("ecology_supp", eco_supp_text),
        ("proteome", prot_text),
    ]:
        assert "/home/" not in txt and "leechuck" not in txt, f"Local filesystem path leaked in {label}"


def test_core_shotgun_exclusions_and_cohort_hierarchy():
    manifest = json.loads(
        (DATA_REPO_DIR / "metadata/cohorts/cohort_product_manifest.json").read_text()
    )
    data_text, _, _ = expand_tex_tree(DATA_MS_DIR / "sn-article.tex")
    eco_main_text, _, _ = expand_tex_tree(ECO_MS_DIR / "alt_stats.tex")
    eco_supp_text, _, _ = expand_tex_tree(ECO_MS_DIR / "supp_altstat.tex")
    prot_text, _, _ = expand_tex_tree(PROT_MS_DIR / "main.tex")

    validate_shotgun_cohort_exclusions_in_text("data", data_text, manifest)
    validate_shotgun_cohort_exclusions_in_text("proteome", prot_text, manifest)
    validate_shotgun_cohort_exclusions_in_text("ecology", eco_main_text + "\n" + eco_supp_text, manifest)


def test_mimag_boundary_and_quality_tiers():
    manifest = json.loads(
        (DATA_REPO_DIR / "metadata/cohorts/cohort_product_manifest.json").read_text()
    )
    b_info = manifest["mimag_quality_accounting"]["strict_mimag_tiers_lt10_contamination"]["boundary_genome_at_10pct"]
    assert b_info["mag"] == "48Dr2__self_SemiBin_25"
    assert round(float(b_info["completeness"]), 2) == 64.61
    assert round(float(b_info["contamination"]), 2) == 10.00
    assert b_info["binned_assembly_accession"] == "ERZ29258845"

    can_reps = pd.read_csv(
        DATA_REPO_DIR / "metadata/metagenome/canonical/refined_representatives.tsv", sep="\t"
    ).set_index("mag")
    assert len(can_reps) == 5438
    b_row = can_reps.loc["48Dr2__self_SemiBin_25"]
    assert float(b_row["Completeness"]) >= 50.0
    assert round(float(b_row["Contamination"]), 2) == 10.00
    assert int((can_reps["Contamination"] == 10.0).sum()) == 1
    assert int((can_reps["Contamination"] < 10.0).sum()) == 5437

    mimag = pd.read_csv(
        PROT_REPO_DIR / "reproducibility/darkmatter/inputs/mimag_rrna_trna_clean/results/mimag_tiers.tsv.gz",
        sep="\t",
    )
    assert len(mimag) == 5438
    assert int((mimag["mimag_tier"] == "high").sum()) == 126
    assert int((mimag["mimag_tier"] == "medium").sum()) == 5312
    assert int(((mimag["mimag_tier"] == "medium") & (mimag["contamination"] < 10.0)).sum()) == 5311
    assert int(((mimag["mimag_tier"] == "medium") & (mimag["contamination"] == 10.0)).sum()) == 1

    data_text, _, _ = expand_tex_tree(DATA_MS_DIR / "sn-article.tex")
    eco_main_text, _, _ = expand_tex_tree(ECO_MS_DIR / "alt_stats.tex")
    eco_supp_text, _, _ = expand_tex_tree(ECO_MS_DIR / "supp_altstat.tex")
    prot_text, _, _ = expand_tex_tree(PROT_MS_DIR / "main.tex")

    validate_mimag_boundary_in_text("data", data_text)
    validate_mimag_boundary_in_text("proteome", prot_text)
    validate_mimag_boundary_in_text("ecology", eco_main_text + "\n" + eco_supp_text)


def test_distinct_protein_and_ko_products():
    manifest = json.loads(
        (DATA_REPO_DIR / "metadata/cohorts/cohort_product_manifest.json").read_text()
    )
    products = manifest["products"]
    assert products["EQ-PROT-BAKTA-v1.0.0"]["protein_count"] == 19186792
    assert products["EQ-PROT-PYRODIGAL-REP-v1.0.0"]["protein_count"] == 20451986

    rec = pd.read_csv(
        PROT_REPO_DIR / "reproducibility/reframed_20261108/inputs/compute_results/genecall_reconcile_20261102/per_genome_reconciliation.tsv",
        sep="\t",
    )
    assert int(rec["n_genes"].sum()) == 20451986
    assert int(rec["n_labelled"].sum()) == 20216302

    prov = json.loads(
        (PROT_REPO_DIR / "reproducibility/reframed_20261108/rep_ko_v2_provenance.json").read_text()
    )
    assert prov["upstream_ko_table"].endswith("04_functional/repko_clean/results/rep_ko_v2.tsv")
    assert prov["relationship_to_named_protein_products"]["EQ-PROT-BAKTA-v1.0.0"]["proteins_total"] == 19186792

    prot_text, _, _ = expand_tex_tree(PROT_MS_DIR / "main.tex")
    norm_prot = normalize_tex_whitespace(prot_text)
    assert "19,186,792" in norm_prot
    assert "20,451,986" in norm_prot
    assert "20,216,302" in norm_prot
    assert "rep\\_ko\\_v2.tsv" in prot_text or "rep_ko_v2.tsv" in norm_prot
    assert "byte- or sequence-level equivalence between the two products has not been established" in norm_prot


def test_canonical_vs_legacy_taxonomy_reconciliation():
    stats = reconcile_canonical_vs_legacy_taxonomy(
        DATA_REPO_DIR / "metadata/metagenome/canonical/refined_representatives.tsv",
        PROT_REPO_DIR / "reproducibility/reframed_20261108/inputs/compute_results/genecall_reconcile_20261102/per_genome_reconciliation.tsv",
    )
    assert stats == {
        "n_genomes": 5438,
        "whole_diff": 1898,
        "domain_diff": 0,
        "phylum_diff": 247,
        "class_diff": 34,
        "order_diff": 220,
        "genus_diff": 1124,
        "thermo_n": 237,
        "thermo_whole_diff": 42,
        "thermo_can_named_orders": 4,
        "thermo_can_unassigned_orders": 2,
        "thermo_can_named_genera": 19,
        "thermo_can_unassigned_genera": 135,
        "thermo_leg_named_genera": 14,
        "thermo_leg_unassigned_genera": 139,
    }

    prot_text, _, _ = expand_tex_tree(PROT_MS_DIR / "main.tex")
    norm_prot = normalize_tex_whitespace(prot_text)
    assert "1,898" in norm_prot
    assert "247 phylum" in norm_prot
    assert "1,124 genus" in norm_prot
    assert "4 named orders plus 2" in norm_prot
    assert "19 named genera plus 135" in norm_prot
    assert "TAXONOMY-SNAPSHOT-SENSITIVITY" in norm_prot
    assert "13 orders" not in norm_prot and "22 genera" not in norm_prot

    manifest_raw = (DATA_REPO_DIR / "metadata/cohorts/cohort_product_manifest.json").read_text(encoding="utf-8")
    assert "gtdb_r232_classification" in manifest_raw and "1898" in manifest_raw, (
        "Data cohort_product_manifest.json has not yet added the canonical vs legacy taxonomy reconciliation entry (1898 whole-string diffs / gtdb_r232_classification)"
    )


def test_ph_measurements_vs_samples_distinction():
    data_text, _, _ = expand_tex_tree(DATA_MS_DIR / "sn-article.tex")
    eco_main_text, _, _ = expand_tex_tree(ECO_MS_DIR / "alt_stats.tex")
    eco_supp_text, _, _ = expand_tex_tree(ECO_MS_DIR / "supp_altstat.tex")

    norm_data = normalize_tex_whitespace(data_text)
    norm_eco = normalize_tex_whitespace(eco_main_text + "\n" + eco_supp_text)

    assert "712" in norm_data and "767" in norm_data
    assert "712" in norm_eco and "767" in norm_eco
    # Must not claim 712 distinct biological samples in active text
    assert "712 distinct biological samples" not in norm_data
    assert "712 distinct biological samples" not in norm_eco


def test_accessions_citations_and_package_snapshots():
    data_text, _, _ = expand_tex_tree(DATA_MS_DIR / "sn-article.tex")
    eco_main_text, _, _ = expand_tex_tree(ECO_MS_DIR / "alt_stats.tex")
    eco_supp_text, _, _ = expand_tex_tree(ECO_MS_DIR / "supp_altstat.tex")
    prot_text, _, _ = expand_tex_tree(PROT_MS_DIR / "main.tex")

    validate_accessions_and_deposit_honesty("data", data_text)
    validate_accessions_and_deposit_honesty("ecology", eco_main_text + "\n" + eco_supp_text)
    validate_accessions_and_deposit_honesty("proteome", prot_text)

    # Check companion citations across the 3 .bib files
    data_bib = (DATA_MS_DIR / "sn-bibliography.bib").read_text(encoding="utf-8")
    eco_bib = (ECO_MS_DIR / "sample.bib").read_text(encoding="utf-8")
    prot_bib = (PROT_MS_DIR / "references.bib").read_text(encoding="utf-8")

    expected_data_title = "A linked resource of microbial genomes, community profiles and environmental data from the Rub' al-Khali"
    expected_eco_title = "Microbial biogeography and functional differentiation across the Rub' al-Khali"
    expected_prot_title = "Widespread desert archaea harbour an export-enriched unknown protein repertoire"

    expected_prot_author_prefix = "Lopez Velazquez, Alejandra and Abdelhakim, Marwa and Alajel, Sulaiman M."

    eco_dp = parse_bib_entry(eco_bib, "rubalkhali_datapaper")
    eco_pr = parse_bib_entry(eco_bib, "proteome_companion")
    assert eco_dp["title"] == expected_data_title
    assert eco_pr["title"] == expected_prot_title
    assert eco_pr["author"].startswith(expected_prot_author_prefix)

    prot_dp = parse_bib_entry(prot_bib, "rubalkhali_datapaper")
    prot_ec = parse_bib_entry(prot_bib, "tawfiq2026companion")
    assert prot_dp["title"] == expected_data_title
    assert prot_ec["title"] == expected_eco_title

    data_ec = parse_bib_entry(data_bib, "ecology_companion")
    data_pr = parse_bib_entry(data_bib, "proteome_companion")
    assert data_ec["title"] == expected_eco_title
    assert data_pr["title"] == expected_prot_title

    def _norm_authors(author_str: str) -> List[str]:
        return [
            a.strip()
            .replace("Lopez-Velazquez", "Lopez Velazquez")
            .replace("Zubair, Jood Kamal", "Zubair, Jood")
            for a in author_str.split(" and ")
            if a.strip()
        ]

    eco_pr_authors = _norm_authors(eco_pr["author"])
    data_pr_authors = _norm_authors(data_pr["author"])
    assert len(eco_pr_authors) == 19
    assert data_pr_authors == eco_pr_authors, (
        f"Ordered 19-author list mismatch in proteome_companion between sn-bibliography.bib and sample.bib/authors.tex: "
        f"{data_pr_authors} != {eco_pr_authors}"
    )

    # Verify protein manuscript snapshot byte-for-byte identity with reproducibility package
    for rel in ["main.tex", "supplementary_notes.tex", "authors.tex", "references.bib", "main.pdf"]:
        ms_bytes = (PROT_MS_DIR / rel).read_bytes()
        pkg_bytes = (PROT_REPO_DIR / "manuscript" / rel).read_bytes()
        assert hashlib.sha256(ms_bytes).hexdigest() == hashlib.sha256(pkg_bytes).hexdigest(), f"Protein package snapshot out of sync for {rel}"


def test_adversarial_mutations_rejected(tmp_path: Path):
    manifest = json.loads(
        (DATA_REPO_DIR / "metadata/cohorts/cohort_product_manifest.json").read_text()
    )
    prot_text, _, _ = expand_tex_tree(PROT_MS_DIR / "main.tex")

    # 1. Fabricated sample exclusion ID must fail
    with pytest.raises(AssertionError, match="fabricated sample exclusion"):
        validate_shotgun_cohort_exclusions_in_text("mutated_exclusions", prot_text + "\nExcluded 62Sr1 and 35Rr3.", manifest)

    # 2. Wrong completeness for boundary genome 48Dr2__self_SemiBin_25 must fail
    with pytest.raises(AssertionError, match="90.2% completeness"):
        validate_mimag_boundary_in_text(
            "mutated_mimag",
            prot_text.replace("64.61\\% completeness", "90.2\\% completeness"),
        )

    # 3. Wrong Zenodo DOI or false archived assertion must fail
    with pytest.raises(AssertionError, match="23134168|23132097"):
        validate_accessions_and_deposit_honesty(
            "mutated_doi",
            prot_text.replace("10.5281/zenodo.23134168", "10.5281/zenodo.23132097"),
        )
    with pytest.raises(AssertionError, match="false completed-archive assertion"):
        validate_accessions_and_deposit_honesty(
            "mutated_archive_claim",
            prot_text + "\nThe archived 19,186,792-protein Bakta release product is archived at Zenodo.",
        )

    # 4. Same-count wrong-genome-ID substitution in taxonomy table must fail
    can_path = DATA_REPO_DIR / "metadata/metagenome/canonical/refined_representatives.tsv"
    leg_path = PROT_REPO_DIR / "reproducibility/reframed_20261108/inputs/compute_results/genecall_reconcile_20261102/per_genome_reconciliation.tsv"
    bad_leg = pd.read_csv(leg_path, sep="\t")
    bad_leg.loc[0, "genome"] = "FABRICATED_MAG_ID_0001"
    bad_leg_path = tmp_path / "bad_legacy.tsv"
    bad_leg.to_csv(bad_leg_path, sep="\t", index=False)
    with pytest.raises(AssertionError, match="Genome ID sets do not match 1:1"):
        reconcile_canonical_vs_legacy_taxonomy(can_path, bad_leg_path)

    # 5. Missing included supplement file must fail
    bad_tex = tmp_path / "main.tex"
    bad_tex.write_text("\\documentclass{article}\n\\input{nonexistent_supplement.tex}\n\\begin{document}\\end{document}\n")
    with pytest.raises(AssertionError, match="Missing included TeX file"):
        expand_tex_tree(bad_tex)


def test_completed_protein_results_and_methods_source_verification():
    """Verify integrated adjusted-export, external-UniRef50, and candidate-audit results and methods against source TSVs and code helpers."""
    comp_dir = PROT_REPO_DIR / "reproducibility/reframed_20261108/inputs/compute_results/science_20261108_completed"
    c_exp = pd.read_csv(comp_dir / "canonical_task_b_export_stage1b_within_genome_and_genus_sensitivity.tsv", sep="\t")
    l_exp = pd.read_csv(comp_dir / "legacy_task_b_export_stage1b_within_genome_and_genus_sensitivity.tsv", sep="\t")
    u50 = pd.read_csv(comp_dir / "task_b_stage2b_uniref50_eq_vs_external_canonical_and_legacy_taxonomy.tsv", sep="\t")
    cand = pd.read_csv(comp_dir / "canonical_task_c_targeted_family_functional_evaluation.tsv", sep="\t")
    s3_sum = json.loads((comp_dir / "canonical_stage3_summary.json").read_text(encoding="utf-8"))

    # 1. Verify source TSV values
    glm_all = c_exp[
        (c_exp["cohort"] == "all_labelled_orfs")
        & (c_exp["endpoint"] == "sp_any")
        & (c_exp["estimator"] == "order_FE_categorical_len7_Binomial_GLM_genus_clustered_SE")
    ].iloc[0]
    glm_comp = c_exp[
        (c_exp["cohort"] == "complete_interior_orfs_only")
        & (c_exp["endpoint"] == "sp_any")
        & (c_exp["estimator"] == "order_FE_categorical_len7_Binomial_GLM_genus_clustered_SE")
    ].iloc[0]
    mh_all = c_exp[
        (c_exp["cohort"] == "all_labelled_orfs")
        & (c_exp["endpoint"] == "sp_any")
        & (c_exp["estimator"] == "within_genome_x_len7_Mantel_Haenszel")
    ].iloc[0]
    leg_glm_all = l_exp[
        (l_exp["cohort"] == "all_labelled_orfs")
        & (l_exp["endpoint"] == "sp_any")
        & (l_exp["estimator"] == "order_FE_categorical_len7_Binomial_GLM_genus_clustered_SE")
    ].iloc[0]

    u50_shared = u50[
        (u50["taxonomy_scheme"] == "canonical_gtdb_r232_data_foundation")
        & (u50["stratum"] == "shared_named_taxonomic_orders_only")
        & (u50["orf_subset"] == "all_orfs")
    ].iloc[0]

    assert len(cand) == 19
    assert int(s3_sum["n_evaluated_total_families"]) == 19
    assert int(s3_sum["n_evaluated_new_dark_candidates"]) == 9
    assert int(s3_sum["n_evaluated_prior_heat_stress_examples"]) == 3
    assert int(s3_sum["n_evaluated_annotated_positive_controls"]) == 0
    assert int(s3_sum["n_evaluated_sequence_annotated_comparators"]) == 2
    assert int(s3_sum["n_evaluated_matched_negative_exported_controls"]) == 5
    assert int(cand["n_total_member_loci"].sum()) == 295
    assert int(cand["n_verified_complete_interior_loci"].sum()) == 217

    prot_text, _, _ = expand_tex_tree(PROT_MS_DIR / "main.tex")
    norm_prot = normalize_tex_whitespace(prot_text)

    # Check source-derived numbers and table labels in manuscript + supplement
    assert f"{glm_all['or_estimate']:.2f} (95\\% genus-clustered CI {glm_all['or_ci95_low']:.2f}--{glm_all['or_ci95_high']:.2f})" in prot_text
    assert f"{glm_comp['or_estimate']:.2f} ({glm_comp['or_ci95_low']:.2f}--{glm_comp['or_ci95_high']:.2f})" in prot_text
    assert f"{mh_all['or_estimate']:.2f} ({mh_all['or_ci95_low']:.2f}--{mh_all['or_ci95_high']:.2f})" in prot_text
    assert f"{leg_glm_all['or_estimate']:.2f} ({leg_glm_all['or_ci95_low']:.2f}--{leg_glm_all['or_ci95_high']:.2f})" in prot_text
    assert f"{u50_shared['diff_mean_genome_pct_dark_eq_minus_ext']:.2f} percentage points" in norm_prot
    assert f"$-0.32$ to {u50_shared['diff_ci95_high']:.2f}" in prot_text
    for lbl in ["\\label{tab:export_adj_completed}", "\\label{tab:ext_uniref50_completed}", "\\label{tab:cand_audit_completed}"]:
        assert lbl in prot_text, f"Missing required supplementary table label {lbl}"

    # 2. Verify code helpers (len_bin_7, audit_orf_integrity, identify_mgnify_gtdb_duplicates) match manuscript Methods
    code_dir = PROT_REPO_DIR / "reproducibility/reframed_20261108/science_20261108_provisional/code"
    if str(code_dir) not in sys.path:
        sys.path.insert(0, str(code_dir))
    import common_science as cs
    import stage1b_refined_fasta_and_corrections as s1b

    probe_lengths = pd.Series([50, 100, 150, 200, 250, 350, 500])
    assert list(s1b.len_bin_7(probe_lengths)) == [
        "<100",
        "100-149",
        "150-199",
        "200-249",
        "250-349",
        "350-499",
        ">=500",
    ]
    assert "$<100$, 100--149, 150--199, 200--249, 250--349, 350--499 and $\\geq 500$" in prot_text

    probe_orfs = pd.DataFrame(
        {
            "partial": ["00", "00", "00", "01"],
            "internal_stop": [0, 0, 0, 0],
            "length": [80, 80, 80, 150],
            "gene_rank": [2, 3, 3, 3],
            "n_genes": [6, 6, 4, 6],
        }
    )
    audited = cs.audit_orf_integrity(
        probe_orfs,
        contig_rank_col="gene_rank",
        contig_n_genes_col="n_genes",
        edge_gene_buffer=2,
    )
    # Row 1 (rank 3 of 6, length 80aa < 100aa, partial=00, internal_stop=0) is complete AND interior
    assert list(audited["orf_status"]) == ["complete", "complete", "complete", "incomplete_or_internal_stop"]
    assert list(audited["contig_edge_status"]) == ["edge", "interior", "edge", "interior"]
    assert "excluding the first and last two genes on each contig" in norm_prot
    assert "99.9\\% skani ANI" in prot_text

    # 3. Ensure false/unsupported method phrases or geographic-restriction claims are absent
    forbidden_phrases = [
        "<50, 50--99",
        "472,882",
        "411,648",
        "1,188,825",
        ">50~bp from both contig ends",
        "consistent with their restriction to desert lineages",
        "assay pending",
        "no domain or assay proof",
        "no catalytic residue or assay",
    ]
    for bad in forbidden_phrases:
        assert bad not in prot_text and bad not in norm_prot, f"Forbidden/unsupported phrase {bad!r} found in protein manuscript"

    # 4. Verify longtable layout and zero Overfull \hbox warnings in compiled main.log
    assert "\\begin{longtable}" in prot_text, "Expected longtable layout for candidate audit table"
    for log_path in [PROT_MS_DIR / "main.log", PROT_REPO_DIR / "manuscript/main.log"]:
        log_text = log_path.read_text(encoding="utf-8", errors="replace")
        assert "Overfull \\hbox" not in log_text, f"Found Overfull \\hbox warning in {log_path}"

    # 5. Verify mechanistic_tests.tex inclusion and 47-family source evidence (38 insuff, 9 fold, 0 biochem)
    prot_all_text, prot_files, prot_figs = expand_tex_tree(PROT_MS_DIR / "main.tex")
    norm_prot_all = normalize_tex_whitespace(prot_all_text)
    assert any(f.name == "mechanistic_tests.tex" for f in prot_files), (
        "mechanistic_tests.tex must be included in protein manuscript tree"
    )

    exp4_json_candidates = [
        PROT_REPO_DIR / "reproducibility/mechanism_20261109/exp4_surface/compact_evidence/exp4_surface_final_summary_v2_audited.json",
        PROT_REPO_DIR / "reproducibility/integration_20261110/surface_figure/exp4_surface_final_summary_v2_audited.json",
    ]
    exp4_json_path = next((p for p in exp4_json_candidates if p.exists()), None)
    assert exp4_json_path is not None, "Missing exp4_surface_final_summary_v2_audited.json in reproducibility package"
    with open(exp4_json_path) as f:
        exp4_summary = json.load(f)

    deep_dive = exp4_summary["deep_dive_dossiers_summary"]
    assert deep_dive["n_dossiers_evaluated"] == 47
    assert deep_dive["n_members_audited_total"] == 428
    assert deep_dive["n_members_verified_interior"] == 326
    disp_counts = deep_dive["disposition_counts"]
    assert disp_counts["uncertain_insufficient_evidence"] == 38
    assert disp_counts["exploratory_fold_or_motif_candidate_no_operon_convergence"] == 9
    assert disp_counts.get("biochemical_function_validated", 0) == 0

    assert "47-family follow-up audited 428 loci" in norm_prot_all
    assert "326 were complete interior coding sequences" in norm_prot_all
    assert "38 families with insufficient evidence and nine with exploratory fold or motif hypotheses" in norm_prot_all



def test_working_draft_pending_registry_and_readiness_gate():
    """Verify required pending placeholders are present and that submission readiness fails until resolved."""
    require_ready = os.environ.get("EQ_REQUIRE_SUBMISSION_READY", "0") == "1"
    if require_ready:
        assert_submission_ready()
        return

    markers = collect_pending_markers()
    joined = "\n".join(f"{loc}:{m}" for loc, m in markers)
    required_tags = [
        "EXT-COMP-HARMON",
        "TAXONOMY-SNAPSHOT-SENSITIVITY",
        "BULK-MAG-ZENODO-DOI",
        "PROTEIN-STRUCTURE-BULK-RELEASE",
        "CODE-SNAPSHOT-RELEASE",
        "DATA-BULK-ARCHIVE-DOI",
    ]
    for tag in required_tags:
        assert tag in joined, f"Required pending blocker {tag!r} missing from active manuscripts"

    completed_tags = [
        "PHYL-EXPORT-ADJ",
        "CAND-AUDIT-A1",
        "EXPORT-RANK-A2",
        "MECH-FAMILIES-C",
    ]
    for ctag in completed_tags:
        assert ctag not in joined, f"Completed analysis tag {ctag!r} must not remain as a [PENDING:...] marker"

    with pytest.raises(SubmissionNotReadyError, match="Submission-readiness gate FAILED"):
        assert_submission_ready()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Three-paper consistency and submission-readiness verifier.")
    ap.add_argument(
        "--require-submission-ready",
        action="store_true",
        help="Fail with non-zero exit code if any [PENDING:...] or DATA-BULK-ARCHIVE-DOI marker remains.",
    )
    args, remaining = ap.parse_known_args()
    if args.require_submission_ready:
        try:
            assert_submission_ready()
        except SubmissionNotReadyError as exc:
            print(str(exc), file=sys.stderr)
            sys.exit(1)
        print("SUBMISSION READINESS: OK")
        sys.exit(0)
    sys.exit(pytest.main([__file__, *remaining]))
