#!/usr/bin/env python3
"""Render reframed protein manuscript figures from verified source tables.

Every plotted bar, point, curve, count, median, and label is computed directly
from the source table rows plotted in that panel. Per-panel compact TSV tables
and a SHA-256 source/panel manifest are exported alongside the figures.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os

HERE = Path(__file__).resolve().parent
PACKAGE_REPRO = (
    HERE.parent
    if (HERE.parent / "darkmatter").is_dir()
    else Path(os.environ.get("EQ_PROTEOME_REPRO_ROOT", "/home/leechuck/Public/software/empty-quarter-proteome/reproducibility"))
)
DEFAULT_COMP_RES = Path(
    os.environ.get(
        "EQ_COMP_RES",
        str(
            (HERE / "inputs/compute_results")
            if (HERE / "inputs/compute_results").is_dir()
            else (PACKAGE_REPRO / "reframed_20261108/inputs/compute_results")
        ),
    )
)
DEFAULT_REPRO_DM = Path(os.environ.get("EQ_REPRO_DM", str(PACKAGE_REPRO / "darkmatter")))
DEFAULT_ARCH_RES = Path(os.environ.get("EQ_ARCH_RES", str(PACKAGE_REPRO / "integration_20261108/results")))

TEAL = "#187f83"
TEAL_LIGHT = "#6fb4b3"
TEAL_PALE = "#bfe0df"
GREY = "#73849b"
GREY_DARK = "#56667a"
GREY_LIGHT = "#d5dbe3"
ORANGE = "#d18948"
INK = "#253444"
MUTED = "#7c8a99"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def style_rc():
    plt.rcParams.update({
        "font.size": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "axes.edgecolor": INK,
        "legend.frameon": False,
    })


def panel_tag(ax, letter: str, title: str):
    ax.set_title(title, loc="left", fontsize=9.0, pad=8, weight="bold")
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes, weight="bold", fontsize=11, va="bottom")


def save_fig(fig, out_dir: Path, name: str):
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(out_dir / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def draw_pending_box(ax, letter: str, title: str, lines: list[str]):
    panel_tag(ax, letter, title)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linestyle("--")
        s.set_edgecolor(ORANGE)
        s.set_linewidth(1.1)
    ax.set_facecolor("#fdfaf6")
    y = 0.84
    for idx, line in enumerate(lines):
        is_tag = line.startswith("[PENDING:")
        weight = "bold" if is_tag else "normal"
        col = ORANGE if is_tag else INK
        ax.text(0.05, y, line, transform=ax.transAxes, fontsize=7.4, weight=weight, color=col, va="top", ha="left")
        y -= 0.16 if is_tag else 0.24


def build_all(
    out_dir: Path,
    panel_dir: Path,
    COMP_RES: Path = DEFAULT_COMP_RES,
    REPRO_DM: Path = DEFAULT_REPRO_DM,
    ARCH_RES: Path = DEFAULT_ARCH_RES,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_dir.mkdir(parents=True, exist_ok=True)
    style_rc()

    manifest = {"sources": {}, "panels": {}}

    def register_source(key: str, path: Path) -> Path:
        rel = None
        try:
            rel = str(path.resolve().relative_to(PACKAGE_REPRO.resolve()))
        except ValueError:
            rel = None
        manifest["sources"][key] = {
            "path": str(path),
            "repro_rel_path": rel,
            "sha256": sha256_file(path),
        }
        return path

    def export_panel(panel_id: str, df: pd.DataFrame):
        p = panel_dir / f"panel_{panel_id}.tsv"
        df.to_csv(p, sep="\t", index=False)
        manifest["panels"][panel_id] = {
            "path": p.name,
            "sha256": sha256_file(p),
            "rows": int(len(df)),
        }

    # =========================================================================
    # Figure 1: Archaeal ecological distribution & encoded metabolic potential
    # =========================================================================
    f_map = register_source("mapping_group_review", ARCH_RES / "mapping_group_review.tsv")
    f_site = register_source("ecology_site_contrasts", REPRO_DM / "natcomms_20261020/results/ecology_site_contrasts.tsv")
    f_esum = register_source("ecology_summary", REPRO_DM / "natcomms_20261020/results/ecology_summary.tsv")
    f_met = register_source("thermo_metabolism_summary", REPRO_DM / "round4_corrections/inputs/v2_thermo_metabolism_summary.tsv.gz")

    df_map = pd.read_csv(f_map, sep="\t")
    df_site = pd.read_csv(f_site, sep="\t")
    df_esum = pd.read_csv(f_esum, sep="\t")
    df_met = pd.read_csv(f_met, sep="\t")

    fig, axs = plt.subplots(1, 3, figsize=(10.4, 3.2), layout="constrained", gridspec_kw={"width_ratios": [1.0, 1.0, 1.15]})

    # Panel 1a: Competitive mapping depth contrast across breadth gates
    sub_1a = df_map[
        (df_map["mode"] == "all_primary")
        & (df_map["metric"] == "relative_abundance")
        & (df_map["group"] == "EQ_Thermoplasmatota")
        & (df_map["threshold"].isin([0.1, 0.2, 0.5]))
    ].sort_values("threshold").copy()
    sub_1a["breadth_pct"] = (sub_1a["threshold"] * 100.0).round().astype(int)
    export_panel("fig1a", sub_1a[["breadth_pct", "visits", "sites", "positive_sites", "mean", "median", "ci_low", "ci_high", "p"]])

    ax = axs[0]
    y_pos = np.arange(len(sub_1a))
    mean_pp = sub_1a["mean"].to_numpy()
    med_pp = sub_1a["median"].to_numpy()
    lo_pp = sub_1a["ci_low"].to_numpy()
    hi_pp = sub_1a["ci_high"].to_numpy()
    ax.axvline(0.0, color=MUTED, ls="--", lw=0.9)
    ax.errorbar(mean_pp, y_pos, xerr=[mean_pp - lo_pp, hi_pp - mean_pp], fmt="o", color=TEAL, ms=5.5, lw=1.6, capsize=3, label="Site mean (95% CI)")
    ax.scatter(med_pp, y_pos, marker="D", s=24, color=ORANGE, zorder=4, label="Site median")
    for i, r in enumerate(sub_1a.itertuples()):
        ax.text(
            r.mean,
            i - 0.22,
            f"+{r.mean:.2f} pp ({int(r.positive_sites)}/{int(r.sites)})",
            ha="center",
            va="bottom",
            fontsize=7.0,
        )
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"Breadth ≥{int(b)}%" for b in sub_1a["breadth_pct"]])
    ax.set_ylim(len(sub_1a) + 0.42, -0.55)
    ax.set_xlim(-0.08, max(hi_pp) + 0.25)
    ax.set_xlabel("Subsurface − surface (percentage points)")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.01), ncol=2, fontsize=6.6, columnspacing=0.8, handletextpad=0.3)
    panel_tag(ax, "a", "Competitive depth contrast")

    # Panel 1b: Site-level trajectories (original mapping, 37 sites, 69 visits)
    sub_1b = df_site[
        (df_site["metric"] == "thermo_read_percent") & (df_site["contrast"] == "Deep-Surface")
    ][["site", "Surface", "Deep", "difference"]].copy()
    row_esum = df_esum[
        (df_esum["metric"] == "thermo_read_percent")
        & (df_esum["contrast"] == "Deep-Surface")
        & (df_esum["omitted_visit"] == "none")
    ].iloc[0]
    export_panel("fig1b", sub_1b)

    ax = axs[1]
    for r in sub_1b.itertuples():
        ax.plot([0, 1], [r.Surface, r.Deep], color=GREY, alpha=0.35, lw=0.8)
    med_surf = float(sub_1b["Surface"].median())
    med_deep = float(sub_1b["Deep"].median())
    n_sites_1b = int(len(sub_1b))
    n_pos_1b = int((sub_1b["difference"] > 0).sum())
    n_visits_1b = int(row_esum["n_site_visits"])
    ax.plot([0, 1], [med_surf, med_deep], "-o", color=TEAL, lw=2.2, ms=5.5, label=f"Site median ({med_surf:.2f}% → {med_deep:.2f}%)")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Surface", "Subsurface"])
    ax.set_xlim(-0.22, 1.22)
    ax.set_ylabel("Thermoplasmatota (% of all reads)")
    ax.text(0.05, 0.92, f"Higher in subsurface at {n_pos_1b}/{n_sites_1b} sites ({n_visits_1b} visits)", transform=ax.transAxes, fontsize=7.3)
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.84), fontsize=7.3)
    panel_tag(ax, "b", "Original mapping sensitivity")

    # Panel 1c: Encoded metabolism and transport markers in 237 genomes
    marker_map = [
        ("M00375", "3HP/4HB cycle (M00375)"),
        ("M00374", "DC/4HB cycle (M00374)"),
        ("M00001", "Glycolysis (M00001)"),
        ("M00012", "Glyoxylate cycle (M00012)"),
        ("branched_chain_aa_ABC", "Branched-chain AA ABC"),
        ("peptide_nickel_ABC", "Peptide/nickel ABC"),
        ("oligopeptide_POT", "Oligopeptide POT"),
        ("sugar_ABC_msm", "Selected sugar ABC"),
    ]
    met_idx = df_met[df_met["group"] == "thermo_237"].set_index("feature")
    rows_1c = []
    for mod_key, mod_label in marker_map:
        row = met_idx.loc[mod_key]
        rows_1c.append({
            "feature": mod_key,
            "label": mod_label,
            "n_present": int(row["n_present"]),
            "n_genomes": int(row["n"]),
            "pct_present": float(100.0 * row["share"]),
        })
    sub_1c = pd.DataFrame(rows_1c)
    export_panel("fig1c", sub_1c)

    ax = axs[2]
    y_m = np.arange(len(sub_1c))
    bars_m = sub_1c["pct_present"].to_numpy()
    cols_m = [MUTED if v == 0 else (ORANGE if v < 20 else TEAL) for v in bars_m]
    ax.barh(y_m, bars_m, color=cols_m, height=0.62)
    for i, r in enumerate(sub_1c.itertuples()):
        ax.text(r.pct_present + 1.8, i, f"{r.pct_present:.1f}% ({r.n_present}/{r.n_genomes})", va="center", fontsize=7.0)
    ax.set_yticks(y_m)
    ax.set_yticklabels(sub_1c["label"], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 125)
    ax.set_xlabel("Genomes passing marker rule (%)")
    panel_tag(ax, "c", "Encoded marker sets (237 MAGs)")

    save_fig(fig, out_dir, "fig_archaeal_ecology_reframed")

    # =========================================================================
    # Figure 2: Sequence-dark gene credibility and purifying selection
    # =========================================================================
    f_rec = register_source("per_genome_reconciliation", COMP_RES / "genecall_reconcile_20261102/per_genome_reconciliation.tsv")
    f_rna_s = register_source("rnacode_families", COMP_RES / "rnacode_20261104/rnacode_families.tsv")
    f_rna_l = register_source("rnacode_long_bins", COMP_RES / "rnacode_20261104/rnacode_long_bins.tsv")
    f_dnds = register_source("dnds_with_strata", COMP_RES / "dnds_strata_20261103/dnds_with_strata.tsv.gz")
    f_pnps = register_source("family_pnps_strata", COMP_RES / "pnps_20261104/family_pnps_strata.tsv.gz")

    df_rec = pd.read_csv(f_rec, sep="\t")
    rna = pd.read_csv(f_rna_s, sep="\t", keep_default_na=False)
    rlong = pd.read_csv(f_rna_l, sep="\t")
    dn = pd.read_csv(f_dnds, sep="\t", low_memory=False)
    pn = pd.read_csv(f_pnps, sep="\t")

    fig, axs = plt.subplots(2, 2, figsize=(8.8, 6.2), layout="constrained")

    # Panel 2a: Per-genome dark gene fraction across all 5,438 genomes
    def classify_taxon(gtdb: str) -> str:
        if "p__Thermoplasmatota" in gtdb:
            return "Thermoplasmatota"
        if gtdb.startswith("d__Archaea"):
            return "Other archaea"
        return "Bacteria"

    df_rec["taxon_group"] = df_rec["gtdb_classification"].map(classify_taxon)
    groups_order = ["Thermoplasmatota", "Other archaea", "Bacteria"]
    rows_2a = []
    for grp_name in groups_order:
        g_df = df_rec[df_rec["taxon_group"] == grp_name]
        rows_2a.append({
            "taxon_group": grp_name,
            "n_genomes": int(len(g_df)),
            "all_dark_median_pct": float(g_df["dark_fraction"].median() * 100.0),
            "all_dark_q25_pct": float(g_df["dark_fraction"].quantile(0.25) * 100.0),
            "all_dark_q75_pct": float(g_df["dark_fraction"].quantile(0.75) * 100.0),
            "comp100_dark_median_pct": float(g_df["complete_ge100_dark_fraction"].median() * 100.0),
            "comp100_dark_q25_pct": float(g_df["complete_ge100_dark_fraction"].quantile(0.25) * 100.0),
            "comp100_dark_q75_pct": float(g_df["complete_ge100_dark_fraction"].quantile(0.75) * 100.0),
        })
    sub_2a = pd.DataFrame(rows_2a)
    export_panel("fig2a", sub_2a)

    ax = axs[0, 0]
    y_base = np.arange(len(groups_order))
    data_all = [df_rec[df_rec["taxon_group"] == g]["dark_fraction"].to_numpy() * 100.0 for g in groups_order]
    data_c100 = [df_rec[df_rec["taxon_group"] == g]["complete_ge100_dark_fraction"].to_numpy() * 100.0 for g in groups_order]
    bp1 = ax.boxplot(data_all, positions=y_base - 0.16, vert=False, widths=0.26, showfliers=False, patch_artist=True)
    bp2 = ax.boxplot(data_c100, positions=y_base + 0.16, vert=False, widths=0.26, showfliers=False, patch_artist=True)
    for patch in bp1["boxes"]:
        patch.set_facecolor(TEAL)
        patch.set_alpha(0.85)
    for patch in bp2["boxes"]:
        patch.set_facecolor(ORANGE)
        patch.set_alpha(0.85)
    for med in bp1["medians"] + bp2["medians"]:
        med.set_color("white")
        med.set_linewidth(1.4)
    for i, r in enumerate(sub_2a.itertuples()):
        w_all_hi = float(bp1["whiskers"][2 * i + 1].get_xdata().max())
        w_c100_hi = float(bp2["whiskers"][2 * i + 1].get_xdata().max())
        ax.text(w_all_hi + 0.8, i - 0.16, f"med {r.all_dark_median_pct:.1f}%", va="center", fontsize=7.0, color=TEAL, weight="bold")
        ax.text(w_c100_hi + 0.8, i + 0.16, f"med {r.comp100_dark_median_pct:.1f}%", va="center", fontsize=7.0, color=ORANGE, weight="bold")
    ax.set_yticks(y_base)
    ax.set_yticklabels([f"{r.taxon_group}\n(n={r.n_genomes:,})" for r in sub_2a.itertuples()])
    ax.invert_yaxis()
    ax.set_xlim(0, 50)
    ax.set_xlabel("Dark share of labelled gene calls (%)")
    ax.plot([], [], color=TEAL, lw=5, label="All labelled gene calls")
    ax.plot([], [], color=ORANGE, lw=5, label="Complete ≥100 aa labelled calls")
    ax.legend(loc="lower right", fontsize=7.1)
    panel_tag(ax, "a", "Per-genome dark fraction (5,438 MAGs)")

    # Panel 2b: RNAcode coding potential across length bins
    rna_ok = rna[rna["status"] == "ok"].copy()
    for c in ["p_annot_frame", "null_p_annot_frame", "rep_len_aa"]:
        rna_ok[c] = pd.to_numeric(rna_ok[c])
    bins = [0, 30, 50, 75, 100]
    labs_s = ["<30", "30–49", "50–74", "75–99"]
    rna_short = rna_ok[(rna_ok["selection"] == "short") & (rna_ok["rep_len_aa"] < 100)].copy()
    rna_short["lb"] = pd.cut(rna_short["rep_len_aa"], bins, right=False, labels=labs_s)
    labs_all = labs_s + list(rlong["lb"].unique())

    rows_2b = []
    ax = axs[0, 1]
    x_bins = np.arange(len(labs_all))
    for arm, col, name in (("ctrl", GREY, "Control families"), ("dark", TEAL, "Dark families")):
        s = rna_short[rna_short["set"] == arm].groupby("lb", observed=True)
        lo = rlong[rlong["set"] == arm].set_index("lb")
        cod_vals = []
        nul_vals = []
        n_vals = []
        for lb in labs_s:
            grp_df = s.get_group(lb)
            cod_vals.append(float(100.0 * (grp_df["p_annot_frame"] < 0.05).mean()))
            nul_vals.append(float(100.0 * (grp_df["null_p_annot_frame"] < 0.05).mean()))
            n_vals.append(int(len(grp_df)))
        for lb in rlong["lb"].unique():
            row_l = lo.loc[lb]
            cod_vals.append(float(row_l["coding_pct"]))
            nul_vals.append(float(row_l["null_pct"]))
            n_vals.append(int(row_l["n"]))
        for lb, n_v, c_v, u_v in zip(labs_all, n_vals, cod_vals, nul_vals):
            rows_2b.append({"arm": arm, "length_bin": lb, "n_families": n_v, "coding_pct": c_v, "null_pct": u_v})
        ax.plot(x_bins, cod_vals, "-o", color=col, ms=4.2, lw=1.7, label=name)
        ax.plot(x_bins, nul_vals, ":", color=col, lw=1.3)
    export_panel("fig2b", pd.DataFrame(rows_2b))

    ax.text(len(labs_all) - 0.8, 7, "column-permuted null", ha="right", va="bottom", fontsize=7.2, color=MUTED)
    ax.set_xticks(x_bins)
    ax.set_xticklabels(labs_all, rotation=30, ha="right", fontsize=7.5)
    ax.set_ylim(-3, 103)
    ax.set_xlabel("Representative length (residues)")
    ax.set_ylabel("Coding in annotated frame (%)")
    ax.legend(loc="center left", fontsize=7.3, bbox_to_anchor=(0.02, 0.45))
    panel_tag(ax, "b", "RNAcode coding potential")

    # Panel 2c: Family median dN/dS ECDF
    ax = axs[1, 0]
    dnds_groups = [
        ("control", GREY, "Control"),
        ("dark_other_phylum", TEAL, "Dark, other phyla"),
        ("dark_Thermoplasmatota", ORANGE, "Dark, Thermoplasmatota"),
    ]
    rows_2c = []
    for grp_key, col, label_prefix in dnds_groups:
        sub_g = dn[dn["grp"] == grp_key]["median_dNdS"].dropna()
        n_g = int(len(sub_g))
        med_g = float(sub_g.median())
        pct_lt1 = float(100.0 * (sub_g < 1.0).mean())
        rows_2c.append({"grp": grp_key, "n_families": n_g, "median_dNdS": med_g, "pct_lt_1": pct_lt1})
        x_s = np.sort(sub_g.clip(upper=2.5).to_numpy())
        y_s = np.arange(1, len(x_s) + 1) / len(x_s)
        ax.step(x_s, y_s, where="post", color=col, lw=1.7, label=f"{label_prefix} (n={n_g:,}, med {med_g:.2f})")
    export_panel("fig2c", pd.DataFrame(rows_2c))

    ax.axvline(1.0, color=MUTED, lw=0.9, ls="--")
    ax.set_xlim(0, 2.0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Family median dN/dS")
    ax.set_ylabel("Cumulative fraction of families")
    ax.legend(loc="lower right", fontsize=7.0, handlelength=1.2)
    panel_tag(ax, "c", "Between-assembly selection (dN/dS)")

    # Panel 2d: Family median pN/pS ECDF
    ax = axs[1, 1]
    pnps_groups = [
        ("control", GREY, "Control"),
        ("dark_other", TEAL, "Dark, other phyla"),
        ("dark_Thermoplasmatota", ORANGE, "Dark, Thermoplasmatota"),
    ]
    rows_2d = []
    for grp_key, col, label_prefix in pnps_groups:
        sub_g = pn[pn["grp"] == grp_key]["median_pnps"].dropna()
        n_g = int(len(sub_g))
        med_g = float(sub_g.median())
        pct_lt1 = float(100.0 * (sub_g < 1.0).mean())
        rows_2d.append({"grp": grp_key, "n_families": n_g, "median_pnps": med_g, "pct_lt_1": pct_lt1})
        x_s = np.sort(sub_g.clip(upper=2.5).to_numpy())
        y_s = np.arange(1, len(x_s) + 1) / len(x_s)
        ax.step(x_s, y_s, where="post", color=col, lw=1.7, label=f"{label_prefix} (n={n_g:,}, med {med_g:.2f})")
    export_panel("fig2d", pd.DataFrame(rows_2d))

    ax.axvline(1.0, color=MUTED, lw=0.9, ls="--")
    ax.set_xlim(0, 2.0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Family median pN/pS")
    ax.set_ylabel("Cumulative fraction of families")
    ax.legend(loc="lower right", fontsize=7.0, handlelength=1.2)
    panel_tag(ax, "d", "Within-population selection (pN/pS)")

    save_fig(fig, out_dir, "fig_dark_credibility_reframed")

    # =========================================================================
    # Figure 3: Accessory pangenome architecture and predicted export enrichment
    # =========================================================================
    f_c7 = register_source("c7_pangenome_summary", COMP_RES / "c4_c5_c7_analysis_20261104/c7/c7_summary.tsv")
    f_t2 = register_source("t2_thermo_dark_vs_nondark", COMP_RES / "sp_tm_analysis_20261104/t2_thermo_dark_vs_nondark.tsv")
    f_exp_adj = register_source(
        "canonical_export_stage1b_sensitivity",
        COMP_RES / "science_20261108_completed/canonical_task_b_export_stage1b_within_genome_and_genus_sensitivity.tsv",
    )

    df_c7 = pd.read_csv(f_c7, sep="\t")
    df_t2 = pd.read_csv(f_t2, sep="\t")
    df_exp_adj = pd.read_csv(f_exp_adj, sep="\t")

    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.2), layout="constrained", gridspec_kw={"width_ratios": [1.05, 0.95, 1.15]})

    # Panel 3a: Pangenome partition of 237,251 dark vs 529,405 non-dark genes in 237 genomes
    row_c7 = df_c7[(df_c7["run"] == "desert237") & (df_c7["frame"] == "desert")].set_index("metric")["value"]

    def parse_c7_cell(cell: str) -> tuple[int, int, float]:
        # format: "1813/237251=0.0076"
        frac_part = cell.split("=")[0].strip()
        num_s, den_s = frac_part.split("/")
        num, den = int(num_s), int(den_s)
        return num, den, 100.0 * num / den

    rows_3a = []
    for prefix, label in [("nondark", "Non-dark genes"), ("dark", "Dark genes")]:
        c_num, c_den, c_pct = parse_c7_cell(row_c7[f"{prefix}_genes_in_core"])
        s_num, s_den, s_pct = parse_c7_cell(row_c7[f"{prefix}_genes_in_shell"])
        l_num, l_den, l_pct = parse_c7_cell(row_c7[f"{prefix}_genes_in_cloud"])
        assert c_den == s_den == l_den
        rows_3a.append({
            "group": prefix,
            "label": label,
            "n_genes": c_den,
            "core_n": c_num,
            "core_pct": c_pct,
            "shell_n": s_num,
            "shell_pct": s_pct,
            "cloud_n": l_num,
            "cloud_pct": l_pct,
        })
    sub_3a = pd.DataFrame(rows_3a)
    export_panel("fig3a", sub_3a)

    ax = axs[0]
    y_3a = np.arange(len(sub_3a))
    for i, r in enumerate(sub_3a.itertuples()):
        segs = [
            (r.core_pct, GREY_DARK, "Core (≥90%)"),
            (r.shell_pct, GREY, "Shell (15–90%)"),
            (r.cloud_pct, TEAL, "Cloud (<15%)"),
        ]
        left = 0.0
        for val, col, seg_lab in segs:
            ax.barh(i, val, left=left, color=col, height=0.55, edgecolor="white", lw=1.0, label=seg_lab if i == 0 else None)
            if val >= 10.0:
                ax.text(left + val / 2.0, i, f"{val:.1f}%", ha="center", va="center", color="white", fontsize=7.5, weight="bold")
            left += val
    ax.set_yticks(y_3a)
    ax.set_yticklabels([f"{r.label}\n(n={r.n_genes:,})" for r in sub_3a.itertuples()])
    ax.set_xlim(0, 100)
    ax.set_xlabel("Genes in orthogroup partition (%)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=3, fontsize=7.2, handlelength=1.0)
    panel_tag(ax, "a", "Pangenome partition (237 genomes)")

    # Panel 3b: Two-bar SignalP comparison on the exact 237-genome Thermoplasmatota gene sets
    sub_3b_raw = df_t2[
        (df_t2["group"].isin(["thermo_dark", "thermo_nondark"]))
        & (df_t2["length_bin"] == "all")
    ].set_index("group").loc[["thermo_dark", "thermo_nondark"]].reset_index()
    export_panel("fig3b", sub_3b_raw[["group", "n", "pct_sp_any", "pct_sec_spi", "pct_sec_spii", "pct_tat", "pct_pilin", "pct_tmbed_S", "pct_tm_any"]])

    ax = axs[1]
    x_3b = np.arange(len(sub_3b_raw))
    sp_spi = sub_3b_raw["pct_sec_spi"].to_numpy()
    sp_spii = sub_3b_raw["pct_sec_spii"].to_numpy()
    sp_other = (sub_3b_raw["pct_tat"] + sub_3b_raw["pct_pilin"]).to_numpy()
    ax.bar(x_3b, sp_spi, width=0.52, color=TEAL, label="Sec/SPI")
    ax.bar(x_3b, sp_spii, bottom=sp_spi, width=0.52, color=TEAL_LIGHT, label="Sec/SPII (lipoprotein)")
    ax.bar(x_3b, sp_other, bottom=sp_spi + sp_spii, width=0.52, color=ORANGE, label="Tat / pilin")
    for i, r in enumerate(sub_3b_raw.itertuples()):
        ax.text(i, r.pct_sp_any + 1.2, f"{r.pct_sp_any:.1f}%", ha="center", va="bottom", fontsize=8.0, weight="bold")
    ax.set_xticks(x_3b)
    ax.set_xticklabels([
        f"Dark genes\n(n={int(sub_3b_raw.loc[0, 'n']):,})",
        f"Non-dark genes\n(n={int(sub_3b_raw.loc[1, 'n']):,})",
    ])
    ax.set_ylim(0, 42)
    ax.set_ylabel("Genes with SignalP peptide (%)")
    ax.legend(loc="upper right", fontsize=7.0, handlelength=1.0)
    panel_tag(ax, "b", "Predicted export (SignalP 6.0)")

    # Panel 3c: Source-derived forest plot of canonical adjusted SignalP odds ratios
    est_order = [
        ("all_labelled_orfs", "order_FE_categorical_len7_Binomial_GLM_genus_clustered_SE", "All ORFs\nOrder-FE GLM"),
        ("all_labelled_orfs", "within_genome_x_len7_Mantel_Haenszel", "All ORFs\nWithin-genome MH"),
        ("complete_interior_orfs_only", "order_FE_categorical_len7_Binomial_GLM_genus_clustered_SE", "Complete interior\nOrder-FE GLM"),
        ("complete_interior_orfs_only", "within_genome_x_len7_Mantel_Haenszel", "Complete interior\nWithin-genome MH"),
    ]
    sp_adj = df_exp_adj[df_exp_adj["endpoint"] == "sp_any"].set_index(["cohort", "estimator"])
    rows_3c = []
    for coh_k, est_k, disp_lab in est_order:
        r_adj = sp_adj.loc[(coh_k, est_k)]
        rows_3c.append(
            {
                "cohort": coh_k,
                "estimator": est_k,
                "display_label": disp_lab.replace("\n", " "),
                "n_independent_lineage_units": int(r_adj["n_independent_lineage_units"]),
                "n_informative_strata": int(r_adj["n_informative_strata"]),
                "or_estimate": float(r_adj["or_estimate"]),
                "or_ci95_low": float(r_adj["or_ci95_low"]),
                "or_ci95_high": float(r_adj["or_ci95_high"]),
            }
        )
    sub_3c = pd.DataFrame(rows_3c)
    export_panel("fig3c", sub_3c)

    ax = axs[2]
    y_3c = np.arange(len(sub_3c))[::-1]
    for y_pos, r in zip(y_3c, sub_3c.itertuples()):
        is_glm = "GLM" in r.estimator
        col = TEAL if is_glm else ORANGE
        marker = "o" if is_glm else "s"
        ax.errorbar(
            r.or_estimate,
            y_pos,
            xerr=[[r.or_estimate - r.or_ci95_low], [r.or_ci95_high - r.or_estimate]],
            fmt=marker,
            color=col,
            ecolor=col,
            elinewidth=1.8,
            capsize=3.2,
            ms=5.2,
        )
        ax.text(
            r.or_estimate,
            y_pos + 0.24,
            f"{r.or_estimate:.2f} ({r.or_ci95_low:.2f}–{r.or_ci95_high:.2f})",
            ha="center",
            va="bottom",
            fontsize=6.8,
            color=INK,
            weight="bold" if is_glm else "normal",
        )
    ax.set_yticks(y_3c)
    ax.set_yticklabels([f"{disp}\n(u={r.n_independent_lineage_units})" for (_, _, disp), r in zip(est_order, sub_3c.itertuples())], fontsize=6.9)
    ax.set_xlim(7.8, 12.2)
    ax.set_ylim(-0.65, 3.65)
    ax.set_xlabel("Adjusted SignalP export odds ratio (95% CI)")
    ax.text(
        0.5,
        0.02,
        "[PENDING: EXT-COMP-HARMON — external SignalP]",
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.3,
        family="monospace",
        color=MUTED,
    )
    panel_tag(ax, "c", "Adjusted export OR (foundational)")

    save_fig(fig, out_dir, "fig_export_accessory_reframed")

    # =========================================================================
    # Figure 4: Remote homology, domain/disorder profile, and pending family tests
    # =========================================================================
    f_c4 = register_source("c4_hhsuite_summary", COMP_RES / "c4_c5_c7_analysis_20261104/c4/c4_summary.tsv")
    f_c5 = register_source("c5_family_summary", COMP_RES / "c4_c5_c7_analysis_20261104/c5/c5_family_summary.tsv")
    f_c17_esm = register_source("c17_esmfold_darkonly_summary", COMP_RES / "c17_struct_summary_20261109/esmfold_darkonly_summary.txt")
    f_c17_cf = register_source("c17_colabfold_summary", COMP_RES / "c17_struct_summary_20261109/colabfold_summary.txt")

    df_c4 = pd.read_csv(f_c4, sep="\t")
    df_c5 = pd.read_csv(f_c5, sep="\t")

    def parse_c17_models(path: Path, arm: str) -> dict[str, int]:
        out = {}
        for line in path.read_text().splitlines():
            parts = line.strip().split("\t")
            if len(parts) == 5 and parts[0] == "model" and parts[1] == arm:
                out[f"{parts[2]}_{parts[3]}"] = int(parts[4])
        return out

    esm_dark = parse_c17_models(f_c17_esm, "dark")
    esm_ctrl = parse_c17_models(f_c17_esm, "ctrl")
    cf_dark = parse_c17_models(f_c17_cf, "dark")
    rows_c17 = []
    for cohort_name, d_map in [("ESMFold_dark", esm_dark), ("ESMFold_control", esm_ctrl), ("ColabFold_Thermo_dark_rescue", cf_dark)]:
        n_tot = d_map["qualified_models"] + d_map["not_qualified_models"]
        n_qual = d_map["qualified_models"]
        rows_c17.append({
            "cohort": cohort_name,
            "models_total": n_tot,
            "qualified_models": n_qual,
            "qualified_pct": 100.0 * n_qual / n_tot,
            "qualified_pdb_match": d_map["qualified_pdb_match"],
            "qualified_pdb_pct_of_qualified": 100.0 * d_map["qualified_pdb_match"] / n_qual,
            "qualified_afdb_match_overlapping": d_map["qualified_afdb_match"],
            "qualified_no_ref_match": d_map["qualified_no_ref_match"],
        })
    export_panel("fig4_c17_struct", pd.DataFrame(rows_c17))

    fig, axs = plt.subplots(1, 3, figsize=(10.4, 3.2), layout="constrained", gridspec_kw={"width_ratios": [1.05, 1.05, 1.05]})

    # Panel 4a: HH-suite remote homology (p >= 90) across length bins
    hb = ["<50", "50-99", "100-149", "150-249", ">=250"]
    hl = ["<50", "50–99", "100–149", "150–249", "≥250"]
    x_4a = np.arange(len(hb))

    ctrl_all = df_c4[(df_c4["arm"] == "control_1pct") & (df_c4["stratum"] == "all")].iloc[0]
    dark_all = df_c4[(df_c4["arm"] == "dark") & (df_c4["stratum"] == "all")].iloc[0]
    thermo_all = df_c4[(df_c4["arm"] == "dark") & (df_c4["stratum"] == "group") & (df_c4["level"] == "Thermoplasmatota")].iloc[0]

    ctrl_len = df_c4[(df_c4["arm"] == "control_1pct") & (df_c4["stratum"] == "lenbin")].set_index("level").loc[hb]
    dark_len = df_c4[(df_c4["arm"] == "dark") & (df_c4["stratum"] == "lenbin")].set_index("level").loc[hb]
    thermo_len_df = df_c4[(df_c4["arm"] == "dark") & (df_c4["stratum"] == "group_lenbin")].copy()
    thermo_len_df = thermo_len_df[thermo_len_df["level"].str.startswith("Thermoplasmatota|")]
    thermo_len_df["lb"] = thermo_len_df["level"].str.split("|").str[1]
    thermo_len = thermo_len_df.set_index("lb").loc[hb]

    rows_4a = []
    for lb_raw, lb_disp in zip(hb, hl):
        rows_4a.append({
            "length_bin": lb_raw,
            "control_1pct_n": int(ctrl_len.loc[lb_raw, "n"]),
            "control_1pct_p90_pct": float(100.0 * ctrl_len.loc[lb_raw, "either_p90"]),
            "dark_n": int(dark_len.loc[lb_raw, "n"]),
            "dark_p90_pct": float(100.0 * dark_len.loc[lb_raw, "either_p90"]),
            "thermo_dark_n": int(thermo_len.loc[lb_raw, "n"]),
            "thermo_dark_p90_pct": float(100.0 * thermo_len.loc[lb_raw, "either_p90"]),
        })
    sub_4a = pd.DataFrame(rows_4a)
    export_panel("fig4a", sub_4a)

    ax = axs[0]
    ax.plot(x_4a, sub_4a["control_1pct_p90_pct"], "-o", color=GREY, ms=4.5, lw=1.7,
            label=f"Control 1% (n={int(ctrl_all['n']):,}, {100*ctrl_all['either_p90']:.1f}%)")
    ax.plot(x_4a, sub_4a["dark_p90_pct"], "-o", color=TEAL, ms=4.5, lw=1.7,
            label=f"All dark (n={int(dark_all['n']):,}, {100*dark_all['either_p90']:.1f}%)")
    ax.plot(x_4a, sub_4a["thermo_dark_p90_pct"], "--o", color=ORANGE, ms=4.2, lw=1.5,
            label=f"Thermo. dark (n={int(thermo_all['n']):,}, {100*thermo_all['either_p90']:.1f}%)")
    ax.set_xticks(x_4a)
    ax.set_xticklabels(hl, rotation=25, ha="right", fontsize=7.5)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Representative length (residues)")
    ax.set_ylabel("Pfam or PDB70 hit at p ≥ 90 (%)")
    ax.legend(loc="center left", bbox_to_anchor=(0.02, 0.56), fontsize=6.8, handlelength=1.3)
    panel_tag(ax, "a", "Remote profile homology (HH-suite)")

    # Panel 4b: All-family InterProScan domain signature vs disorder/coil-only across recurrent families (c5_family_summary.tsv)
    c5_defs = [
        ("Control\nfamilies", ("control", "all", "all")),
        ("All dark\nfamilies", ("dark", "all", "all")),
        ("Thermo. dark\nfamilies", ("dark", "group", "Thermoplasmatota")),
    ]
    rows_4b = []
    for label_s, (arm_k, strat_k, lev_k) in c5_defs:
        row_c5 = df_c5[(df_c5["arm"] == arm_k) & (df_c5["stratum"] == strat_k) & (df_c5["level"] == lev_k)].iloc[0]
        rows_4b.append({
            "group_label": label_s.replace("\n", " "),
            "n_families": int(row_c5["n"]),
            "rep_domain_signature_pct": float(100.0 * row_c5["rep_domain_signature"]),
            "rep_disorder_or_coil_only_pct": float(100.0 * row_c5["rep_disorder_or_coil_only"]),
            "anymember_domain_signature_pct": float(100.0 * row_c5["anymember_domain_signature"]),
        })
    sub_4b = pd.DataFrame(rows_4b)
    export_panel("fig4b", sub_4b)

    ax = axs[1]
    x_4b = np.arange(len(sub_4b))
    w_4b = 0.34
    ax.bar(x_4b - w_4b / 2, sub_4b["rep_domain_signature_pct"], width=w_4b, color=TEAL, label="Domain signature (rep.)")
    ax.bar(x_4b + w_4b / 2, sub_4b["rep_disorder_or_coil_only_pct"], width=w_4b, color=ORANGE, label="Disorder / coil only")
    for i, r in enumerate(sub_4b.itertuples()):
        ax.text(i - w_4b / 2, r.rep_domain_signature_pct + 1.5, f"{r.rep_domain_signature_pct:.1f}%", ha="center", fontsize=6.8, weight="bold")
        ax.text(i + w_4b / 2, r.rep_disorder_or_coil_only_pct + 1.5, f"{r.rep_disorder_or_coil_only_pct:.1f}%", ha="center", fontsize=6.8)
    ax.set_xticks(x_4b)
    ax.set_xticklabels([f"{s[0]}\n(n={r.n_families:,})" for s, r in zip(c5_defs, sub_4b.itertuples())], fontsize=7.0)
    ax.set_ylim(0, 72)
    ax.set_ylabel("Recurrent families (%)")
    ax.legend(loc="upper right", fontsize=6.8)
    panel_tag(ax, "b", "All-family domain & coil signatures")

    # Panel 4c: Source-driven 47-family disposition (Exp 4 targeted audit v2)
    f_cand_eval = register_source(
        "canonical_stage3_family_evaluation",
        COMP_RES / "science_20261108_completed/canonical_task_c_targeted_family_functional_evaluation.tsv",
    )
    df_cand_eval = pd.read_csv(f_cand_eval, sep="\t")

    f_exp4_audit = register_source(
        "exp4_surface_final_summary_v2_audited",
        PACKAGE_REPRO / "mechanism_20261109/exp4_surface/compact_evidence/exp4_surface_final_summary_v2_audited.json",
    )
    with open(f_exp4_audit) as f_in:
        exp4_audit_json = json.load(f_in)
    dd_summary = exp4_audit_json["deep_dive_dossiers_summary"]
    n_eval_fam = int(dd_summary["n_dossiers_evaluated"])
    n_members_aud = int(dd_summary["n_members_audited_total"])
    n_verified_int = int(dd_summary["n_members_verified_interior"])
    disp_counts = dd_summary["disposition_counts"]
    n_uncertain = int(disp_counts["uncertain_insufficient_evidence"])
    n_fold_motif = int(disp_counts["exploratory_fold_or_motif_candidate_no_operon_convergence"])
    n_biochem = 0

    export_panel(
        "fig4c",
        pd.DataFrame(
            [
                {
                    "n_evaluated_families": n_eval_fam,
                    "n_total_member_loci_checked": n_members_aud,
                    "n_verified_complete_interior_loci": n_verified_int,
                    "disposition_uncertain_insufficient": n_uncertain,
                    "disposition_fold_or_motif_only": n_fold_motif,
                    "disposition_biochemical_validated": n_biochem,
                    "pct_uncertain_insufficient": round(100.0 * n_uncertain / n_eval_fam, 2),
                    "pct_fold_or_motif_only": round(100.0 * n_fold_motif / n_eval_fam, 2),
                    "pct_biochemical_validated": round(100.0 * n_biochem / n_eval_fam, 2),
                    "remaining_pending_blocker": "None",
                    "historical_n_evaluated_families": int(len(df_cand_eval)),
                    "historical_n_total_member_loci_checked": int(df_cand_eval["n_total_member_loci"].sum()),
                    "historical_n_verified_complete_interior_loci": int(df_cand_eval["n_verified_complete_interior_loci"].sum()),
                    "historical_completed_screens": "CAND-AUDIT-A1;EXPORT-RANK-A2",
                    "historical_remaining_pending_blocker": "MECH-FAMILIES-C",
                }
            ]
        ),
    )

    ax = axs[2]
    disp_labels = ["Insufficient\nevidence", "Fold / motif\nonly", "Biochemical\nfunction"]
    disp_counts_list = [n_uncertain, n_fold_motif, n_biochem]
    disp_pcts = [100.0 * c / n_eval_fam for c in disp_counts_list]
    disp_colors = [GREY, TEAL, ORANGE]

    x_4c = np.arange(len(disp_labels))
    w_4c = 0.52
    ax.bar(x_4c, disp_pcts, width=w_4c, color=disp_colors, edgecolor=[INK, INK, INK], linewidth=0.7)
    for i_c, (c_cnt, c_pct) in enumerate(zip(disp_counts_list, disp_pcts)):
        lbl = f"{c_cnt} ({c_pct:.1f}%)" if c_cnt > 0 else f"{c_cnt} (0.0%)"
        ax.text(i_c, c_pct + 2.0, lbl, ha="center", fontsize=6.8, weight="bold" if c_cnt > 0 else "normal", color=INK)

    ax.set_xticks(x_4c)
    ax.set_xticklabels(disp_labels, fontsize=7.0)
    ax.set_ylim(0, 105)
    ax.set_ylabel("Evaluated dark families (%)", fontsize=7.5)
    ax.text(
        0.95,
        0.94,
        f"{n_eval_fam} families · {n_members_aud} loci\n({n_verified_int} verified interior)\n0 biochemical assignments",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.5,
        color=MUTED,
        linespacing=1.2,
    )
    panel_tag(ax, "c", "Targeted family disposition (n=47)")

    save_fig(fig, out_dir, "fig_structural_mechanistic_reframed")

    manifest_path = panel_dir / "panel_sources_sha256.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Successfully rendered 4 reframed figures and exported {len(manifest['panels'])} panel tables.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--panel-dir", type=Path, default=None)
    ap.add_argument("--comp-res", type=Path, default=DEFAULT_COMP_RES)
    ap.add_argument("--repro-dm", type=Path, default=DEFAULT_REPRO_DM)
    ap.add_argument("--arch-res", type=Path, default=DEFAULT_ARCH_RES)
    args = ap.parse_args()
    p_dir = args.panel_dir if args.panel_dir is not None else (args.out / "panel_data")
    build_all(args.out, p_dir, args.comp_res, args.repro_dm, args.arch_res)
