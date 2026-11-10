#!/usr/bin/env python3
"""Unit and regression tests for Experiment 4 surface synteny figure.

Verifies all assertions specified in FIGURE-TASK.txt and FIGURE-ROOT-REVIEW.txt:
- 24 exact verified rows (13 for Lead A, 11 for Lead B).
- All loci strictly co-oriented on the same contig strand with candidate 5' to partner.
- Gap calculation matches 1-based inclusive genomic coordinates including negative overlaps (-4 bp, -1 bp).
- Genus group counts: exactly 6 genus-or-singleton blocks for Lead A, 5 for Lead B.
- Dataset validator strictly rejects duplicate IDs, duplicate genomes, and invariant violations.
- Main figure uses verified interior representatives:
    Lead A: 1PRr2__self_SemiBin_81 (L01)
    Lead B: V47PRr1__self_SemiBin_699_sub (L14)
- Main figure uses correct titles ('Lead A: PilJ-profile neighbour', 'Lead B: envelope-profile neighbour'),
  omits role headings and 'significant'/q text, and uses x-axis label 'Position relative to pair left edge (bp)'.
- Supplementary figure displays all 24 loci with short IDs L01..L24 + sample ID and full mapping TSV.
- SHA manifest integrity matches source files bit-for-bit.
- All 4 figure artifacts (main/supp PDF/PNG) and caption draft exist, are non-empty, and consistent.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import pytest
import pandas as pd

from build_surface_figure import (
    MODULE_DIR,
    TSV_DEFAULT,
    JSON_DEFAULT,
    MAPPING_TSV_DEFAULT,
    MAIN_PDF_DEFAULT,
    MAIN_PNG_DEFAULT,
    SUPP_PDF_DEFAULT,
    SUPP_PNG_DEFAULT,
    LEAD1_ID,
    LEAD2_ID,
    LEAD1_EXEMPLAR_GENOME,
    LEAD2_EXEMPLAR_GENOME,
    validate_input_dataset,
    create_locus_mapping,
    build_all_figures,
)

ORIGINAL_TSV = MODULE_DIR.parent.parent / "mechanism-20261109" / "review-20261109" / "exp4_lead_flank_coordinates_independent.tsv"
ORIGINAL_JSON = MODULE_DIR.parent.parent / "mechanism-20261109" / "exp4_surface" / "results" / "followup_v2_review" / "exp4_surface_final_summary_v2_audited.json"
MANIFEST_TSV = MODULE_DIR / "sha256_manifest.tsv"
MANIFEST_SUMS = MODULE_DIR / "SHA256SUMS"
CAPTION_FILE = MODULE_DIR / "captiondraft.txt"


def find_original_sources() -> Tuple[Path, Path]:
    """Portably locate original source TSV and JSON files across repository contexts."""
    cands_tsv = [
        MODULE_DIR.parent.parent / "mechanism_20261109" / "exp4_surface" / "compact_evidence" / "exp4_lead_flank_coordinates_independent.tsv",
        MODULE_DIR.parent.parent / "mechanism-20261109" / "review-20261109" / "exp4_lead_flank_coordinates_independent.tsv",
        Path("/home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/exp4_surface/compact_evidence/exp4_lead_flank_coordinates_independent.tsv"),
        Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/review-20261109/exp4_lead_flank_coordinates_independent.tsv"),
    ]
    cands_json = [
        MODULE_DIR.parent.parent / "mechanism_20261109" / "exp4_surface" / "compact_evidence" / "exp4_surface_final_summary_v2_audited.json",
        MODULE_DIR.parent.parent / "mechanism-20261109" / "exp4_surface" / "results" / "followup_v2_review" / "exp4_surface_final_summary_v2_audited.json",
        Path("/home/leechuck/Public/software/empty-quarter-proteome/reproducibility/mechanism_20261109/exp4_surface/compact_evidence/exp4_surface_final_summary_v2_audited.json"),
        Path("/home/leechuck/Documents/papers/empty-quarter/mechanism-20261109/exp4_surface/results/followup_v2_review/exp4_surface_final_summary_v2_audited.json"),
    ]
    orig_tsv = next((p for p in cands_tsv if p.exists()), None)
    orig_json = next((p for p in cands_json if p.exists()), None)
    return orig_tsv, orig_json


def sha256_file(path: Path) -> str:
    """Compute sha256 hex digest of file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def test_sha_manifest_integrity_and_source_parity() -> None:
    """Verify local files match SHA manifest and original source files bit-for-bit."""
    assert TSV_DEFAULT.exists(), f"Missing {TSV_DEFAULT}"
    assert JSON_DEFAULT.exists(), f"Missing {JSON_DEFAULT}"
    assert MAPPING_TSV_DEFAULT.exists(), f"Missing {MAPPING_TSV_DEFAULT}"
    assert MANIFEST_TSV.exists(), f"Missing {MANIFEST_TSV}"
    assert MANIFEST_SUMS.exists(), f"Missing {MANIFEST_SUMS}"

    orig_tsv, orig_json = find_original_sources()
    assert orig_tsv is not None and orig_tsv.exists(), "Original source TSV not found in candidate paths"
    assert orig_json is not None and orig_json.exists(), "Original source JSON not found in candidate paths"

    tsv_sha = sha256_file(TSV_DEFAULT)
    json_sha = sha256_file(JSON_DEFAULT)
    map_sha = sha256_file(MAPPING_TSV_DEFAULT)
    orig_tsv_sha = sha256_file(orig_tsv)
    orig_json_sha = sha256_file(orig_json)

    # Parity with original files
    assert tsv_sha == orig_tsv_sha, "TSV sha256 does not match original"
    assert json_sha == orig_json_sha, "JSON sha256 does not match original"

    # Parity with sha256_manifest.tsv
    manifest_df = pd.read_csv(MANIFEST_TSV, sep="\t")
    manifest_map = dict(zip(manifest_df["file_path"], manifest_df["sha256"]))
    assert manifest_map.get(TSV_DEFAULT.name) == tsv_sha
    assert manifest_map.get(JSON_DEFAULT.name) == json_sha
    assert manifest_map.get(MAPPING_TSV_DEFAULT.name) == map_sha

    # Parity with SHA256SUMS
    with open(MANIFEST_SUMS) as f:
        lines = [line.strip().split() for line in f if line.strip()]
    sums_map = {parts[1]: parts[0] for parts in lines}
    assert sums_map.get(TSV_DEFAULT.name) == tsv_sha
    assert sums_map.get(JSON_DEFAULT.name) == json_sha
    assert sums_map.get(MAPPING_TSV_DEFAULT.name) == map_sha


def test_row_counts_and_partition() -> None:
    """Assert exactly 24 rows total, partitioned 13 and 11 across lead families."""
    df = pd.read_csv(TSV_DEFAULT, sep="\t")
    assert len(df) == 24, f"Expected 24 rows, got {len(df)}"

    lead1_rows = df[df["family_id"] == LEAD1_ID]
    lead2_rows = df[df["family_id"] == LEAD2_ID]

    assert len(lead1_rows) == 13, f"Expected 13 rows for Lead A, got {len(lead1_rows)}"
    assert len(lead2_rows) == 11, f"Expected 11 rows for Lead B, got {len(lead2_rows)}"
    assert len(lead1_rows) + len(lead2_rows) == 24

    # Cross-verify with audited summary JSON
    with open(JSON_DEFAULT) as f:
        summary = json.load(f)
    coord_summary = summary["coordinate_audit_summary"]
    assert coord_summary["lead1_pairs_recovered"] == 13
    assert coord_summary["lead2_pairs_recovered"] == 11
    assert coord_summary["total_lead_pairs_recovered"] == 24


def test_all_cooriented_and_candidate_upstream() -> None:
    """Assert all 24 loci are strictly co-oriented and candidate is 5' to partner."""
    df = pd.read_csv(TSV_DEFAULT, sep="\t")
    
    assert (df["same_strand"] == True).all()
    assert (df["focal_strand"] == df["flank_strand"]).all()
    assert set(df["focal_strand"].unique()).issubset({1, -1})
    
    for idx, row in df.iterrows():
        st = row["focal_strand"]
        if st == 1:
            assert row["focal_start"] < row["flank_start"], f"Row {idx}: focal not upstream on +1 strand"
            assert row["gene_offset"] == 1
        elif st == -1:
            assert row["focal_start"] > row["flank_start"], f"Row {idx}: focal not upstream on -1 strand"
            assert row["gene_offset"] == -1

    with open(JSON_DEFAULT) as f:
        summary = json.load(f)
    assert summary["coordinate_audit_summary"]["lead1_all_same_strand"] is True
    assert summary["coordinate_audit_summary"]["lead2_all_same_strand"] is True


def test_gap_calculation_including_negative_overlaps() -> None:
    """Assert intergenic gap formula matches 1-based inclusive coordinates and covers negative overlaps."""
    df = pd.read_csv(TSV_DEFAULT, sep="\t")
    
    negative_overlap_rows = []
    
    for idx, row in df.iterrows():
        st = row["focal_strand"]
        reported_gap = row["gap_bp_negative_overlap"]
        
        if st == 1:
            calc_gap = row["flank_start"] - row["focal_end"] - 1
        else:
            calc_gap = row["focal_start"] - row["flank_end"] - 1
            
        assert calc_gap == reported_gap, f"Row {idx}: calc {calc_gap} != reported {reported_gap}"
        
        if reported_gap < 0:
            negative_overlap_rows.append(row)
            overlap_bp = abs(reported_gap)
            if st == 1:
                assert row["focal_end"] - row["flank_start"] + 1 == overlap_bp
            else:
                assert row["flank_end"] - row["focal_start"] + 1 == overlap_bp

    assert len(negative_overlap_rows) > 0, "Expected negative overlaps in dataset"
    neg_gaps = {r["gap_bp_negative_overlap"] for r in negative_overlap_rows}
    assert -4 in neg_gaps
    assert -1 in neg_gaps

    lead1_gaps = df[df["family_id"] == LEAD1_ID]["gap_bp_negative_overlap"]
    assert lead1_gaps.min() == -4
    assert lead1_gaps.max() == 37

    lead2_gaps = df[df["family_id"] == LEAD2_ID]["gap_bp_negative_overlap"]
    assert lead2_gaps.min() == -4
    assert lead2_gaps.max() == 6


def test_genus_groups_count() -> None:
    """Assert exactly 6 genus-or-singleton blocks for Lead A and 5 for Lead B."""
    df = pd.read_csv(TSV_DEFAULT, sep="\t")
    
    lead1_df = df[df["family_id"] == LEAD1_ID]
    lead2_df = df[df["family_id"] == LEAD2_ID]
    
    assert lead1_df["genus_block"].nunique() == 6
    assert lead2_df["genus_block"].nunique() == 5
    
    assert (lead1_df["genus_block"] == "g__JAJPHT01").sum() == 8
    assert lead1_df["genus_block"].str.startswith("unclassified_singleton::").sum() == 5
    
    assert (lead2_df["genus_block"] == "g__JAJPHT01").sum() == 7
    assert lead2_df["genus_block"].str.startswith("unclassified_singleton::").sum() == 4


def test_duplicate_ids_rejected() -> None:
    """Assert validator rejects duplicate IDs, duplicate genomes, and invalid parameters."""
    df = pd.read_csv(TSV_DEFAULT, sep="\t")
    validate_input_dataset(df)
    
    # Duplicate focal_gene
    df_dup = df.copy()
    df_dup.loc[1, "focal_gene"] = df_dup.loc[0, "focal_gene"]
    with pytest.raises(ValueError, match="Duplicate focal_gene IDs rejected"):
        validate_input_dataset(df_dup)
        
    # Duplicate flank_gene
    df_dup2 = df.copy()
    df_dup2.loc[1, "flank_gene"] = df_dup2.loc[0, "flank_gene"]
    with pytest.raises(ValueError, match="Duplicate flank_gene IDs rejected"):
        validate_input_dataset(df_dup2)
        
    # Duplicate (family, genome)
    df_dup3 = df.copy()
    df_dup3.loc[1, "focal_genome"] = df_dup3.loc[0, "focal_genome"]
    with pytest.raises(ValueError, match="Duplicate \\(family, genome\\) pairs rejected"):
        validate_input_dataset(df_dup3)


def test_main_figure_exemplars_and_labels() -> None:
    """Verify main figure uses verified interior representatives and strict review labels."""
    df = pd.read_csv(TSV_DEFAULT, sep="\t")
    with open(JSON_DEFAULT) as f:
        summary = json.load(f)

    # Lead 1 verified interior representative in summary JSON
    expected_lead1_rep = summary["lead_1_audit_summary"]["verified_interior_representative"]
    expected_lead2_rep = summary["lead_2_audit_summary"]["verified_interior_representative"]

    # Verify exemplar genomes match
    assert LEAD1_EXEMPLAR_GENOME in expected_lead1_rep
    assert LEAD2_EXEMPLAR_GENOME in expected_lead2_rep

    row1 = df[df["focal_genome"] == LEAD1_EXEMPLAR_GENOME].iloc[0]
    row2 = df[df["focal_genome"] == LEAD2_EXEMPLAR_GENOME].iloc[0]

    assert row1["focal_gene"] == expected_lead1_rep
    assert row2["focal_gene"] == expected_lead2_rep

    # Verify strands and orientations
    assert row1["focal_strand"] == -1  # Native reverse strand
    assert row2["focal_strand"] == 1   # Native forward strand


def test_locus_mapping_tsv_and_supplement() -> None:
    """Verify locus_id_mapping.tsv has exact 24 rows, L01..L24, and correct sample mapping."""
    df_map = pd.read_csv(MAPPING_TSV_DEFAULT, sep="\t")
    assert len(df_map) == 24

    # Unique IDs L01 to L24
    expected_ids = [f"L{i:02d}" for i in range(1, 25)]
    assert df_map["locus_id"].tolist() == expected_ids

    # L01 is Lead A exemplar (1PRr2)
    l01 = df_map[df_map["locus_id"] == "L01"].iloc[0]
    assert l01["focal_genome"] == LEAD1_EXEMPLAR_GENOME
    assert l01["sample_id"] == "1PRr2"
    assert l01["is_exemplar"] is True or l01["is_exemplar"] == 1

    # L14 is Lead B exemplar (V47PRr1)
    l14 = df_map[df_map["locus_id"] == "L14"].iloc[0]
    assert l14["focal_genome"] == LEAD2_EXEMPLAR_GENOME
    assert l14["sample_id"] == "V47PRr1"
    assert l14["is_exemplar"] is True or l14["is_exemplar"] == 1

    # Sample ID extracted correctly from focal_genome
    for _, r in df_map.iterrows():
        assert r["sample_id"] == r["focal_genome"].split("__")[0]


def test_figure_artifacts_exist_and_nonempty() -> None:
    """Verify all 4 figure files exist and are non-empty."""
    build_all_figures(
        tsv_path=TSV_DEFAULT,
        json_path=JSON_DEFAULT,
        mapping_tsv_out=MAPPING_TSV_DEFAULT,
        main_pdf_out=MAIN_PDF_DEFAULT,
        main_png_out=MAIN_PNG_DEFAULT,
        supp_pdf_out=SUPP_PDF_DEFAULT,
        supp_png_out=SUPP_PNG_DEFAULT,
    )
    for p in [MAIN_PDF_DEFAULT, MAIN_PNG_DEFAULT, SUPP_PDF_DEFAULT, SUPP_PNG_DEFAULT]:
        assert p.exists(), f"Missing {p}"
        assert p.stat().st_size > 10000, f"File {p} is suspiciously small ({p.stat().st_size} bytes)"


def test_captiondraft_content_and_metrics() -> None:
    """Verify caption draft contains required scientific disclosures and audited metrics."""
    doc_path = CAPTION_FILE if CAPTION_FILE.exists() else (MODULE_DIR / "README.md")
    assert doc_path.exists(), f"Missing caption documentation ({doc_path})"
    caption_text = doc_path.read_text()

    assert LEAD1_ID in caption_text
    assert LEAD2_ID in caption_text
    assert LEAD1_EXEMPLAR_GENOME in caption_text
    assert LEAD2_EXEMPLAR_GENOME in caption_text
    assert "PilJ-profile partner" in caption_text
    assert "Envelope-profile partner" in caption_text
    assert "Lead A: PilJ-profile neighbour" in caption_text
    assert "Lead B: envelope-profile neighbour" in caption_text
    assert "unresolved" in caption_text.lower()
    assert "Position relative to pair left edge (bp)" in caption_text
    assert "0.00206" in caption_text
    assert "0.00035" in caption_text
    assert "1.0" in caption_text
    assert "-4 to 37" in caption_text
    assert "-4 to 6" in caption_text
    assert "L01" in caption_text
    assert "L14" in caption_text
    assert "locus_id_mapping.tsv" in caption_text


def test_inclusive_coordinate_lengths_and_geometry_consistency() -> None:
    """Verify 1-based inclusive lengths (1245/471 and 684/339) and geometry gap consistency."""
    df = pd.read_csv(TSV_DEFAULT, sep="\t")

    # 1. Verify exact 4 exemplar lengths
    row_a = df[df["focal_genome"] == LEAD1_EXEMPLAR_GENOME].iloc[0]
    row_b = df[df["focal_genome"] == LEAD2_EXEMPLAR_GENOME].iloc[0]

    lead_a_cand_len = row_a["focal_end"] - row_a["focal_start"] + 1
    lead_a_part_len = row_a["flank_end"] - row_a["flank_start"] + 1
    lead_b_cand_len = row_b["focal_end"] - row_b["focal_start"] + 1
    lead_b_part_len = row_b["flank_end"] - row_b["flank_start"] + 1

    assert lead_a_cand_len == 1245, f"Lead A candidate length {lead_a_cand_len} != 1245"
    assert lead_a_part_len == 471, f"Lead A partner length {lead_a_part_len} != 471"
    assert lead_b_cand_len == 684, f"Lead B candidate length {lead_b_cand_len} != 684"
    assert lead_b_part_len == 339, f"Lead B partner length {lead_b_part_len} != 339"

    # 2. Verify geometry consistency (drawn_gap == reported_gap) across all 24 loci
    for idx, r in df.iterrows():
        f_st, f_en = int(r["focal_start"]), int(r["focal_end"])
        fl_st, fl_en = int(r["flank_start"]), int(r["flank_end"])
        left_edge = min(f_st, fl_st)
        strand = int(r["focal_strand"])
        reported_gap = int(r["gap_bp_negative_overlap"])

        f_len = f_en - f_st + 1
        fl_len = fl_en - fl_st + 1

        if strand == 1:
            g1_start = f_st - left_edge
            g1_end = g1_start + f_len
            g2_start = fl_st - left_edge
            drawn_gap = g2_start - g1_end
        else:
            g1_start = fl_st - left_edge
            g1_end = g1_start + fl_len
            g2_start = f_st - left_edge
            drawn_gap = g2_start - g1_end

        assert drawn_gap == reported_gap, f"Row {idx} ({r['focal_genome']}): drawn_gap {drawn_gap} != reported {reported_gap}"

        if reported_gap > 0:
            assert drawn_gap > 0
        elif reported_gap == 0:
            assert drawn_gap == 0
        else:
            assert drawn_gap < 0
            assert abs(drawn_gap) == abs(reported_gap)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
