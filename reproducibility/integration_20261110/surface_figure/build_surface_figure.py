#!/usr/bin/env python3
"""Build publication-quality main and supplementary synteny figures for Exp 4 surface modules.

Adheres strictly to root visual review (FIGURE-ROOT-REVIEW.txt & FIGURE-FINAL-REPAIR.txt):
1. Main figure (fig_surface_modules.pdf/png):
   - Fixed 11.0 inch wide by 5.5 inch high figure.
   - Clean 2-panel layout for verified interior representatives:
     * Lead A: 1PRr2__self_SemiBin_81 (exact table row, strand -1)
     * Lead B: V47PRr1__self_SemiBin_699_sub (exact table row, strand +1)
   - Headings: 'Lead A: PilJ-profile neighbour' and 'Lead B: envelope-profile neighbour'.
   - Subtitles: '13 loci across 6 genus/singleton groups' and '11 loci across 5 groups'.
   - Inside arrows: ONLY 'Candidate' and 'Partner'.
   - Lengths below arrows: '1,244 bp', '470 bp', '683 bp', '338 bp'.
   - Intergenic gaps annotated above: 'gap: +5 bp'.
   - Single top common legend: Candidate (teal) and Partner (orange). No per-panel duplicates.
   - Native orientation (strand -1 left, strand +1 right); actual bp axis with common scale (0-1800 bp).
   - X-axis label: 'Position relative to pair left edge (bp)'.
   - Font sizes 10-11 pt throughout; no text collisions or canvas overflow.

2. Supplementary figure (fig_surface_modules_all_loci.pdf/png):
   - Retains all 24 verified loci (13 for Lead A, 11 for Lead B).
   - Removed contig span and strand text columns (arrows encode strand and geometry).
   - Short unique IDs L01..L24 + sample ID, with full mapping TSV (locus_id_mapping.tsv).
   - Crisp, uncrowded arrow labels with no overflow.

3. Mapping TSV (locus_id_mapping.tsv):
   - Full mapping of L01..L24 to genome, sample, genus-or-singleton block, coordinates, and exemplar status.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MODULE_DIR = Path(__file__).resolve().parent
TSV_DEFAULT = MODULE_DIR / "exp4_lead_flank_coordinates_independent.tsv"
JSON_DEFAULT = MODULE_DIR / "exp4_surface_final_summary_v2_audited.json"
MAPPING_TSV_DEFAULT = MODULE_DIR / "locus_id_mapping.tsv"

MAIN_PDF_DEFAULT = MODULE_DIR / "fig_surface_modules.pdf"
MAIN_PNG_DEFAULT = MODULE_DIR / "fig_surface_modules.png"
SUPP_PDF_DEFAULT = MODULE_DIR / "fig_surface_modules_all_loci.pdf"
SUPP_PNG_DEFAULT = MODULE_DIR / "fig_surface_modules_all_loci.png"

# Color palette: high contrast, publication standard, colorblind-safe
COLOR_CANDIDATE = "#0E7490"         # Deep Teal
COLOR_CANDIDATE_BORDER = "#155E75"  # Dark Teal border
COLOR_PARTNER = "#D97706"           # Warm Amber / Orange
COLOR_PARTNER_BORDER = "#9A3412"    # Dark Orange border
COLOR_OVERLAP = "#B91C1C"           # Crimson red for negative overlaps
COLOR_TEXT_MAIN = "#111827"         # Near black
COLOR_TEXT_MUTED = "#4B5563"        # Slate gray
COLOR_BG_TRACK = "#F8FAFC"          # Very light blue-gray
COLOR_GRID = "#E2E8F0"
COLOR_DIVIDER = "#CBD5E1"
COLOR_EXEMPLAR_TAG = "#047857"      # Forest green for exemplar indicator

LEAD1_ID = "22PRr1::22PRr1_2531955_3"
LEAD2_ID = "V60PRr2::V60PRr2_1314979_4"

LEAD1_EXEMPLAR_GENOME = "1PRr2__self_SemiBin_81"
LEAD2_EXEMPLAR_GENOME = "V47PRr1__self_SemiBin_699_sub"


def validate_input_dataset(df: pd.DataFrame) -> None:
    """Validate dataset against core scientific invariants."""
    if len(df) != 24:
        raise ValueError(f"Expected exactly 24 rows, found {len(df)}")
    
    fam_counts = df["family_id"].value_counts().to_dict()
    if fam_counts.get(LEAD1_ID) != 13:
        raise ValueError(f"Expected 13 rows for {LEAD1_ID}, found {fam_counts.get(LEAD1_ID)}")
    if fam_counts.get(LEAD2_ID) != 11:
        raise ValueError(f"Expected 11 rows for {LEAD2_ID}, found {fam_counts.get(LEAD2_ID)}")
    
    # Duplicate ID rejection
    if df["focal_gene"].duplicated().any():
        dups = df["focal_gene"][df["focal_gene"].duplicated()].tolist()
        raise ValueError(f"Duplicate focal_gene IDs rejected: {dups}")
    if df["flank_gene"].duplicated().any():
        dups = df["flank_gene"][df["flank_gene"].duplicated()].tolist()
        raise ValueError(f"Duplicate flank_gene IDs rejected: {dups}")
    
    # Duplicate (family, genome) rejection
    dup_genomes = df.duplicated(subset=["family_id", "focal_genome"])
    if dup_genomes.any():
        dups = df.loc[dup_genomes, ["family_id", "focal_genome"]].to_dict(orient="records")
        raise ValueError(f"Duplicate (family, genome) pairs rejected: {dups}")
    
    # Co-orientation validation
    if not (df["focal_strand"] == df["flank_strand"]).all():
        raise ValueError("Violated invariant: all focal and flank strands must be identical")
    if not (df["same_strand"] == True).all():
        raise ValueError("Violated invariant: all same_strand entries must be True")
    
    # Gap calculation validation (1-based inclusive coordinates)
    for idx, row in df.iterrows():
        st = row["focal_strand"]
        if st == 1:
            calc_gap = row["flank_start"] - row["focal_end"] - 1
        elif st == -1:
            calc_gap = row["focal_start"] - row["flank_end"] - 1
        else:
            raise ValueError(f"Invalid strand value: {st}")
        
        if calc_gap != row["gap_bp_negative_overlap"]:
            raise ValueError(
                f"Row {idx} ({row['focal_genome']}): calculated gap {calc_gap} "
                f"!= reported gap {row['gap_bp_negative_overlap']}"
            )
    
    # Genus block count validation
    lead1_genus = df[df["family_id"] == LEAD1_ID]["genus_block"].nunique()
    if lead1_genus != 6:
        raise ValueError(f"Expected 6 genus blocks for Lead 1, got {lead1_genus}")
    lead2_genus = df[df["family_id"] == LEAD2_ID]["genus_block"].nunique()
    if lead2_genus != 5:
        raise ValueError(f"Expected 5 genus blocks for Lead 2, got {lead2_genus}")


def create_locus_mapping(df: pd.DataFrame, out_path: Path | None = None) -> pd.DataFrame:
    """Create structured mapping dataframe with short unique IDs L01..L24 and metadata."""
    df_lead1 = df[df["family_id"] == LEAD1_ID].copy()
    df_lead1["is_exemplar"] = df_lead1["focal_genome"] == LEAD1_EXEMPLAR_GENOME
    df_lead1["is_shared_genus"] = df_lead1["genus_block"] == "g__JAJPHT01"
    df_lead1 = df_lead1.sort_values(
        by=["is_exemplar", "is_shared_genus", "genus_block", "focal_genome"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    df_lead1["locus_id"] = [f"L{i:02d}" for i in range(1, 14)]
    df_lead1["lead_panel"] = "Lead A"
    df_lead1["lead_heading"] = "Lead A: PilJ-profile neighbour"

    df_lead2 = df[df["family_id"] == LEAD2_ID].copy()
    df_lead2["is_exemplar"] = df_lead2["focal_genome"] == LEAD2_EXEMPLAR_GENOME
    df_lead2["is_shared_genus"] = df_lead2["genus_block"] == "g__JAJPHT01"
    df_lead2 = df_lead2.sort_values(
        by=["is_exemplar", "is_shared_genus", "genus_block", "focal_genome"],
        ascending=[False, False, True, True],
    ).reset_index(drop=True)
    df_lead2["locus_id"] = [f"L{i:02d}" for i in range(14, 25)]
    df_lead2["lead_panel"] = "Lead B"
    df_lead2["lead_heading"] = "Lead B: envelope-profile neighbour"

    mapping_df = pd.concat([df_lead1, df_lead2], ignore_index=True)
    mapping_df["sample_id"] = mapping_df["focal_genome"].apply(lambda x: x.split("__")[0])

    cols = [
        "locus_id", "lead_panel", "lead_heading", "family_id", "sample_id",
        "focal_genome", "genus_block", "is_shared_genus", "is_exemplar",
        "focal_gene", "flank_gene", "focal_strand", "flank_strand",
        "focal_start", "focal_end", "flank_start", "flank_end",
        "gap_bp_negative_overlap", "flanking_guilds",
    ]
    mapping_df = mapping_df[cols]
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        mapping_df.to_csv(out_path, sep="\t", index=False)
    return mapping_df


def draw_gene_arrow(
    ax: plt.Axes,
    x0: float,
    x1: float,
    y: float,
    height: float,
    strand: int,
    facecolor: str,
    edgecolor: str,
    linewidth: float = 1.2,
    head_max_bp: float = 140.0,
    zorder: int = 4,
) -> None:
    """Draw a gene arrow with directional arrowhead strictly matching contig strand."""
    length = abs(x1 - x0)
    left = min(x0, x1)
    right = max(x0, x1)
    
    head_len = min(head_max_bp, 0.28 * length)
    half_h = height / 2.0
    
    if strand == 1:
        # Points right ->
        body_end = right - head_len
        verts = [
            (left, y - half_h),
            (body_end, y - half_h),
            (right, y),
            (body_end, y + half_h),
            (left, y + half_h),
        ]
    elif strand == -1:
        # Points left <-
        body_start = left + head_len
        verts = [
            (left, y),
            (body_start, y - half_h),
            (right, y - half_h),
            (right, y + half_h),
            (body_start, y + half_h),
        ]
    else:
        verts = [
            (left, y - half_h),
            (right, y - half_h),
            (right, y + half_h),
            (left, y + half_h),
        ]
        
    poly = patches.Polygon(
        verts,
        closed=True,
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=zorder,
    )
    ax.add_patch(poly)


def build_main_figure(
    df: pd.DataFrame,
    pdf_out: Path = MAIN_PDF_DEFAULT,
    png_out: Path = MAIN_PNG_DEFAULT,
) -> None:
    """Build compact, publication-ready two-exemplar main figure (11.0 in x 5.5 in)."""
    row_a = df[df["focal_genome"] == LEAD1_EXEMPLAR_GENOME].iloc[0]
    row_b = df[df["focal_genome"] == LEAD2_EXEMPLAR_GENOME].iloc[0]

    # Exactly 11.0 in wide by 5.5 in high as requested
    fig, (ax_top, ax_bot) = plt.subplots(
        nrows=2,
        ncols=1,
        figsize=(11.0, 5.5),
        dpi=300,
        facecolor="white",
        gridspec_kw={"hspace": 0.55, "top": 0.88, "bottom": 0.12, "left": 0.06, "right": 0.96},
    )

    common_max_x = 1800.0
    x_ticks = [0, 250, 500, 750, 1000, 1250, 1500, 1750]

    def render_exemplar(
        ax: plt.Axes,
        row: pd.Series,
        panel_letter: str,
        lead_title: str,
        subtitle_text: str,
        show_xlabel: bool = False,
    ) -> None:
        ax.set_xlim(-60, common_max_x + 60)
        ax.set_ylim(-0.85, 1.40)
        ax.set_facecolor("white")

        for spine in ["top", "right", "left"]:
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color(COLOR_DIVIDER)
        ax.spines["bottom"].set_linewidth(1.0)

        # Panel Title (fontsize 11 bold)
        ax.text(
            -60,
            1.24,
            f"{panel_letter}   {lead_title}",
            fontsize=11.0,
            fontweight="bold",
            color=COLOR_TEXT_MAIN,
            ha="left",
            va="center",
        )

        # Subtitle (fontsize 10)
        ax.text(
            -60,
            0.90,
            subtitle_text,
            fontsize=10.0,
            fontweight="normal",
            color=COLOR_TEXT_MUTED,
            ha="left",
            va="center",
        )

        # Compute relative coordinates (1-based inclusive: length = end - start + 1)
        f_st, f_en = int(row["focal_start"]), int(row["focal_end"])
        fl_st, fl_en = int(row["flank_start"]), int(row["flank_end"])
        left_edge = min(f_st, fl_st)
        strand = int(row["focal_strand"])
        gap_bp = int(row["gap_bp_negative_overlap"])

        focal_len = f_en - f_st + 1
        flank_len = fl_en - fl_st + 1

        rel_f_st = f_st - left_edge
        rel_f_en = rel_f_st + focal_len
        rel_fl_st = fl_st - left_edge
        rel_fl_en = rel_fl_st + flank_len

        total_span = max(rel_f_en, rel_fl_en)

        # Baseline
        ax.plot([0, total_span], [0, 0], color="#94A3B8", linewidth=1.2, linestyle="-", zorder=2)

        # Vertical gridlines
        for xt in x_ticks:
            ax.axvline(xt, color=COLOR_GRID, linestyle="--", linewidth=0.6, zorder=1)

        arrow_height = 0.52

        # Draw Candidate (Teal)
        draw_gene_arrow(
            ax=ax,
            x0=rel_f_st,
            x1=rel_f_en,
            y=0,
            height=arrow_height,
            strand=strand,
            facecolor=COLOR_CANDIDATE,
            edgecolor=COLOR_CANDIDATE_BORDER,
            linewidth=1.2,
            zorder=4,
        )

        # Draw Partner (Orange)
        draw_gene_arrow(
            ax=ax,
            x0=rel_fl_st,
            x1=rel_fl_en,
            y=0,
            height=arrow_height,
            strand=strand,
            facecolor=COLOR_PARTNER,
            edgecolor=COLOR_PARTNER_BORDER,
            linewidth=1.2,
            zorder=4,
        )

        # Inside arrows: ONLY 'Candidate' and 'Partner'
        cand_mid = (rel_f_st + rel_f_en) / 2.0
        ax.text(
            cand_mid,
            0,
            "Candidate",
            fontsize=10.5,
            color="white",
            fontweight="bold",
            ha="center",
            va="center",
            zorder=5,
        )

        part_mid = (rel_fl_st + rel_fl_en) / 2.0
        ax.text(
            part_mid,
            0,
            "Partner",
            fontsize=10.5,
            color="white",
            fontweight="bold",
            ha="center",
            va="center",
            zorder=5,
        )

        # Lengths below arrows: short '1,244 bp', '470 bp', etc.
        ax.text(
            cand_mid,
            -0.42,
            f"{focal_len:,} bp",
            fontsize=10.0,
            color=COLOR_TEXT_MUTED,
            ha="center",
            va="top",
            zorder=5,
        )
        ax.text(
            part_mid,
            -0.42,
            f"{flank_len:,} bp",
            fontsize=10.0,
            color=COLOR_TEXT_MUTED,
            ha="center",
            va="top",
            zorder=5,
        )

        # Intergenic gap annotation above arrows
        if strand == 1:
            gap_center = (rel_f_en + rel_fl_st) / 2.0
        else:
            gap_center = (rel_fl_en + rel_f_st) / 2.0
            
        gap_text = f"gap: +{gap_bp} bp" if gap_bp > 0 else (f"gap: {gap_bp} bp" if gap_bp < 0 else "gap: 0 bp")
        ax.text(
            gap_center,
            0.42,
            gap_text,
            fontsize=9.5,
            fontweight="bold",
            color=COLOR_TEXT_MAIN,
            ha="center",
            va="bottom",
            zorder=6,
        )
        ax.plot([gap_center, gap_center], [0.12, 0.38], color=COLOR_TEXT_MUTED, linewidth=0.7, linestyle=":")

        # Ticks and labels
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f"{t:,}" for t in x_ticks], fontsize=9.0, color=COLOR_TEXT_MUTED)
        ax.set_yticks([])

        if show_xlabel:
            ax.set_xlabel("Position relative to pair left edge (bp)", fontsize=10.0, color=COLOR_TEXT_MAIN, labelpad=6)
        else:
            ax.tick_params(axis="x", labelbottom=True)

    # Render Panel a
    render_exemplar(
        ax=ax_top,
        row=row_a,
        panel_letter="a",
        lead_title="Lead A: PilJ-profile neighbour",
        subtitle_text="13 loci across 6 genus/singleton groups",
        show_xlabel=False,
    )

    # Render Panel b
    render_exemplar(
        ax=ax_bot,
        row=row_b,
        panel_letter="b",
        lead_title="Lead B: envelope-profile neighbour",
        subtitle_text="11 loci across 5 groups",
        show_xlabel=True,
    )

    # Single top common legend: Candidate (teal) and Partner (orange)
    legend_handles = [
        patches.Patch(facecolor=COLOR_CANDIDATE, edgecolor=COLOR_CANDIDATE_BORDER, label="Candidate"),
        patches.Patch(facecolor=COLOR_PARTNER, edgecolor=COLOR_PARTNER_BORDER, label="Partner"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.96, 0.96),
        ncol=2,
        frameon=False,
        fontsize=10.5,
        handletextpad=0.6,
        columnspacing=1.5,
    )

    pdf_out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(pdf_out, format="pdf", bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.savefig(png_out, format="png", dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Main figure generated:\n  PDF: {pdf_out}\n  PNG: {png_out}")


def build_supplementary_figure(
    mapping_df: pd.DataFrame,
    pdf_out: Path = SUPP_PDF_DEFAULT,
    png_out: Path = SUPP_PNG_DEFAULT,
) -> None:
    """Build uncrowded supplementary figure showing all 24 verified loci with L01..L24 labels."""
    df_lead1 = mapping_df[mapping_df["lead_panel"] == "Lead A"].copy()
    df_lead2 = mapping_df[mapping_df["lead_panel"] == "Lead B"].copy()

    fig = plt.figure(figsize=(11.0, 14.0), dpi=300, facecolor="white")
    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[13.6, 11.6],
        left=0.06,
        right=0.96,
        top=0.94,
        bottom=0.06,
        hspace=0.28,
    )

    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0])

    max_x = 1800.0
    x_ticks = [0, 250, 500, 750, 1000, 1250, 1500, 1750]

    def render_supp_panel(
        ax: plt.Axes,
        df_fam: pd.DataFrame,
        panel_letter: str,
        lead_heading: str,
        partner_label: str,
        context_str: str,
    ) -> None:
        n_rows = len(df_fam)
        y_positions = np.arange(n_rows)[::-1]

        ax.set_ylim(-1.0, n_rows + 2.2)
        ax.set_xlim(-400, max_x + 220)
        ax.set_facecolor("white")

        for spine in ["top", "right", "left"]:
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color(COLOR_DIVIDER)
        ax.spines["bottom"].set_linewidth(1.0)

        # Title and Subtitle
        ax.text(
            -390,
            n_rows + 1.6,
            f"{panel_letter}   {lead_heading}",
            fontsize=11.5,
            fontweight="bold",
            color=COLOR_TEXT_MAIN,
            ha="left",
            va="center",
        )
        ax.text(
            -390,
            n_rows + 1.1,
            context_str,
            fontsize=9.0,
            fontweight="medium",
            color=COLOR_TEXT_MUTED,
            ha="left",
            va="center",
        )

        # Column Headers
        header_y = n_rows + 0.35
        ax.text(-390, header_y, "Locus ID (Sample)", fontsize=8.5, fontweight="bold", color=COLOR_TEXT_MAIN, ha="left", va="bottom")
        ax.text(-200, header_y, "Genus / Singleton", fontsize=8.5, fontweight="bold", color=COLOR_TEXT_MAIN, ha="left", va="bottom")
        ax.text(0, header_y, "Syntenic Gene Architecture", fontsize=8.5, fontweight="bold", color=COLOR_TEXT_MAIN, ha="left", va="bottom")
        ax.text(max_x + 35, header_y, "Intergenic Gap", fontsize=8.5, fontweight="bold", color=COLOR_TEXT_MAIN, ha="left", va="bottom")
        ax.axhline(header_y - 0.15, xmin=0.01, xmax=0.99, color=COLOR_DIVIDER, linewidth=0.8)

        # Vertical guide gridlines
        for xt in x_ticks:
            ax.axvline(xt, color=COLOR_GRID, linestyle="--", linewidth=0.5, zorder=1)

        row_h = 0.52

        for i, (_, row) in enumerate(df_fam.iterrows()):
            y = y_positions[i]
            is_ex = bool(row["is_exemplar"])

            # Alternating subtle row background
            if i % 2 == 0:
                bg_rect = patches.Rectangle(
                    (-395, y - 0.44),
                    max_x + 610,
                    0.88,
                    facecolor="#F8FAFC",
                    edgecolor="none",
                    zorder=0,
                )
                ax.add_patch(bg_rect)

            locus_id = row["locus_id"]
            sample_id = row["sample_id"]
            genus_block = row["genus_block"]
            strand = int(row["focal_strand"])
            gap_bp = int(row["gap_bp_negative_overlap"])

            if genus_block == "g__JAJPHT01":
                genus_display = "g__JAJPHT01"
                genus_col = "#1D4ED8"
                genus_wt = "bold"
            else:
                genus_display = "singleton"
                genus_col = "#6B7280"
                genus_wt = "normal"

            id_text = f"{locus_id} ({sample_id})"
            if is_ex:
                ax.text(-390, y, id_text, fontsize=8.5, fontfamily="monospace", fontweight="bold", color=COLOR_TEXT_MAIN, ha="left", va="center")
                badge = patches.FancyBboxPatch(
                    (-265, y - 0.22),
                    55,
                    0.44,
                    boxstyle="round,pad=0.2",
                    facecolor="#DCFCE7",
                    edgecolor="#16A34A",
                    linewidth=0.7,
                    zorder=2,
                )
                ax.add_patch(badge)
                ax.text(-237, y, "exemplar", fontsize=6.5, fontweight="bold", color="#15803D", ha="center", va="center", zorder=3)
            else:
                ax.text(-390, y, id_text, fontsize=8.5, fontfamily="monospace", color=COLOR_TEXT_MAIN, ha="left", va="center")

            ax.text(-200, y, genus_display, fontsize=8.0, color=genus_col, fontweight=genus_wt, ha="left", va="center")

            f_st, f_en = int(row["focal_start"]), int(row["focal_end"])
            fl_st, fl_en = int(row["flank_start"]), int(row["flank_end"])
            left_edge = min(f_st, fl_st)

            focal_len = f_en - f_st + 1
            flank_len = fl_en - fl_st + 1

            rel_f_st = f_st - left_edge
            rel_f_en = rel_f_st + focal_len
            rel_fl_st = fl_st - left_edge
            rel_fl_en = rel_fl_st + flank_len

            track_end = max(rel_f_en, rel_fl_en)

            ax.plot([0, track_end], [y, y], color="#94A3B8", linewidth=0.9, linestyle="-", zorder=2)

            # Draw Arrows
            draw_gene_arrow(
                ax=ax,
                x0=rel_f_st,
                x1=rel_f_en,
                y=y,
                height=row_h,
                strand=strand,
                facecolor=COLOR_CANDIDATE,
                edgecolor=COLOR_CANDIDATE_BORDER,
                linewidth=1.0,
                zorder=4,
            )
            draw_gene_arrow(
                ax=ax,
                x0=rel_fl_st,
                x1=rel_fl_en,
                y=y,
                height=row_h,
                strand=strand,
                facecolor=COLOR_PARTNER,
                edgecolor=COLOR_PARTNER_BORDER,
                linewidth=1.0,
                zorder=4,
            )

            # In-arrow text: ONLY 'Candidate'/'Cand' or 'Partner'/'Part', no overflow
            cand_mid = (rel_f_st + rel_f_en) / 2.0
            cand_lbl = "Candidate" if focal_len >= 800 else ("Cand" if focal_len >= 450 else "")
            if cand_lbl:
                ax.text(cand_mid, y, cand_lbl, fontsize=7.5, color="white", fontweight="bold", ha="center", va="center", zorder=5)

            part_mid = (rel_fl_st + rel_fl_en) / 2.0
            part_lbl = "Partner" if flank_len >= 380 else ("Part" if flank_len >= 250 else "")
            if part_lbl:
                ax.text(part_mid, y, part_lbl, fontsize=7.5, color="white", fontweight="bold", ha="center", va="center", zorder=5)

            # Intergenic Gap
            if gap_bp < 0:
                gap_text = f"{gap_bp} bp (overlap)"
                gap_col = COLOR_OVERLAP
                gap_wt = "bold"
            elif gap_bp == 0:
                gap_text = "0 bp (adjacent)"
                gap_col = COLOR_TEXT_MAIN
                gap_wt = "normal"
            else:
                gap_text = f"+{gap_bp} bp"
                gap_col = COLOR_TEXT_MAIN
                gap_wt = "normal"

            ax.text(max_x + 35, y, gap_text, fontsize=8.0, color=gap_col, fontweight=gap_wt, ha="left", va="center")

        ax.set_xticks(x_ticks)
        ax.set_xticklabels([f"{t:,} bp" for t in x_ticks], fontsize=8.0, color=COLOR_TEXT_MUTED)
        ax.set_yticks([])
        ax.set_xlabel("Position relative to pair left edge (bp)", fontsize=8.5, color=COLOR_TEXT_MUTED, labelpad=6)

        legend_handles = [
            patches.Patch(facecolor=COLOR_CANDIDATE, edgecolor=COLOR_CANDIDATE_BORDER, label="Candidate"),
            patches.Patch(facecolor=COLOR_PARTNER, edgecolor=COLOR_PARTNER_BORDER, label=partner_label),
        ]
        ax.legend(
            handles=legend_handles,
            loc="upper right",
            bbox_to_anchor=(0.99, 1.14),
            frameon=True,
            facecolor="white",
            edgecolor=COLOR_DIVIDER,
            fontsize=8.0,
            ncol=2,
        )

    # Render Panel a
    render_supp_panel(
        ax=ax_top,
        df_fam=df_lead1,
        panel_letter="a",
        lead_heading="Lead A: PilJ-profile neighbour (all 13 verified loci)",
        partner_label="PilJ-profile partner",
        context_str="13 context pairs across 6 genus/singleton groups  |  Intergenic gap: -4 to 37 bp",
    )

    # Render Panel b
    render_supp_panel(
        ax=ax_bot,
        df_fam=df_lead2,
        panel_letter="b",
        lead_heading="Lead B: envelope-profile neighbour (all 11 verified loci)",
        partner_label="Envelope-profile partner",
        context_str="11 context pairs across 5 groups  |  Intergenic gap: -4 to 6 bp",
    )

    pdf_out.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(pdf_out, format="pdf", bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.savefig(png_out, format="png", dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Supplementary figure generated:\n  PDF: {pdf_out}\n  PNG: {png_out}")


def build_all_figures(
    tsv_path: Path = TSV_DEFAULT,
    json_path: Path = JSON_DEFAULT,
    mapping_tsv_out: Path = MAPPING_TSV_DEFAULT,
    main_pdf_out: Path = MAIN_PDF_DEFAULT,
    main_png_out: Path = MAIN_PNG_DEFAULT,
    supp_pdf_out: Path = SUPP_PDF_DEFAULT,
    supp_png_out: Path = SUPP_PNG_DEFAULT,
) -> None:
    """Build both main exemplar figure and supplementary 24-loci figure."""
    df = pd.read_csv(tsv_path, sep="\t")
    validate_input_dataset(df)

    mapping_df = create_locus_mapping(df, out_path=mapping_tsv_out)
    print(f"Locus mapping TSV saved to {mapping_tsv_out}")

    build_main_figure(
        df=df,
        pdf_out=main_pdf_out,
        png_out=main_png_out,
    )

    build_supplementary_figure(
        mapping_df=mapping_df,
        pdf_out=supp_pdf_out,
        png_out=supp_png_out,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate main and supplementary synteny figures.")
    parser.add_argument("--tsv", type=Path, default=TSV_DEFAULT, help="Input coordinates TSV")
    parser.add_argument("--json", type=Path, default=JSON_DEFAULT, help="Input audit summary JSON")
    parser.add_argument("--mapping-tsv", type=Path, default=MAPPING_TSV_DEFAULT, help="Output mapping TSV")
    parser.add_argument("--main-pdf", type=Path, default=MAIN_PDF_DEFAULT, help="Main figure PDF")
    parser.add_argument("--main-png", type=Path, default=MAIN_PNG_DEFAULT, help="Main figure PNG")
    parser.add_argument("--supp-pdf", type=Path, default=SUPP_PDF_DEFAULT, help="Supplementary figure PDF")
    parser.add_argument("--supp-png", type=Path, default=SUPP_PNG_DEFAULT, help="Supplementary figure PNG")
    args = parser.parse_args()

    build_all_figures(
        tsv_path=args.tsv,
        json_path=args.json,
        mapping_tsv_out=args.mapping_tsv,
        main_pdf_out=args.main_pdf,
        main_png_out=args.main_png,
        supp_pdf_out=args.supp_pdf,
        supp_png_out=args.supp_png,
    )


if __name__ == "__main__":
    main()
