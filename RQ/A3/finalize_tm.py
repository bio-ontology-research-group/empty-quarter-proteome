#!/usr/bin/env python3
"""Recompute RQ-A3 dual-novel confidence and length sensitivity.

Dual-novel means best TM < 0.5 against both AFDB and the curated PDB; a query
with no hit in a database receives TM=0 for that database. The analysis
universe is the set of ESMFold records with ``status=ok`` and numeric mean
pLDDT, not the union of the two hit tables.

Examples
--------
On IBEX, with ``M=/ibex/scratch/projects/c2014/rob/metagenomics``::

    python3 RQ/A3/finalize_tm.py

With an already concatenated local shard file::

    python3 RQ/A3/finalize_tm.py --plddt /path/to/plddt_all.tsv
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
THRESHOLDS = (0.00, 0.40, 0.50, 0.60, 0.70, 0.80)


def load_tm(path: Path) -> dict[str, float]:
    scores: dict[str, float] = {}
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) < 2 or row[0] in {"rep", "query"}:
                continue
            try:
                scores[row[0]] = float(row[1])
            except ValueError:
                continue
    return scores


def load_plddt(paths: list[Path]) -> tuple[list[tuple[str, int, float]], dict[str, int]]:
    models: dict[str, tuple[int, float]] = {}
    statuses: dict[str, int] = {}
    for path in paths:
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("rep") == "rep":
                    continue
                status = row.get("status", "")
                statuses[status] = statuses.get(status, 0) + 1
                if status != "ok":
                    continue
                try:
                    value = (int(row["len"]), float(row["mean_plddt"]))
                except (KeyError, TypeError, ValueError):
                    continue
                rep = row["rep"]
                if rep in models:
                    raise ValueError(f"duplicate pLDDT representative: {rep}")
                models[rep] = value
    return [(rep, length, plddt) for rep, (length, plddt) in models.items()], statuses


def summarise(
    models: list[tuple[str, int, float]],
    statuses: dict[str, int],
    afdb: dict[str, float],
    pdb: dict[str, float],
) -> tuple[dict[str, str], list[tuple[float, int, int, float]], list[tuple[str, int, int, float]]]:
    rows = [
        (
            rep,
            length,
            plddt,
            afdb.get(rep, 0.0) < 0.5 and pdb.get(rep, 0.0) < 0.5,
        )
        for rep, length, plddt in models
    ]

    threshold_rows = []
    for threshold in THRESHOLDS:
        retained = [row for row in rows if row[2] >= threshold]
        n_dual = sum(row[3] for row in retained)
        rate = 100.0 * n_dual / len(retained)
        threshold_rows.append((threshold, len(retained), n_dual, rate))

    confident = [row for row in rows if row[2] >= 0.60]
    length_specs = (
        ("lt100", "<100", lambda length: length < 100),
        ("100_150", "100–<150", lambda length: 100 <= length < 150),
        ("150_250", "150–<250", lambda length: 150 <= length < 250),
        ("250_400", "250–<400", lambda length: 250 <= length < 400),
        ("ge400", "≥400", lambda length: length >= 400),
    )
    length_rows = []
    for stem, label, predicate in length_specs:
        selected = [row for row in confident if predicate(row[1])]
        n_dual = sum(row[3] for row in selected)
        rate = 100.0 * n_dual / len(selected)
        length_rows.append((stem, label, len(selected), n_dual, rate))

    dual = [row for row in rows if row[3]]
    matched = [row for row in rows if not row[3]]
    summary = {
        "n_plddt_records": str(sum(statuses.values())),
        "n_plddt_ok": str(len(rows)),
        "n_plddt_too_long": str(statuses.get("too_long", 0)),
        "dual_novel_median_plddt": f"{statistics.median(row[2] for row in dual):.3f}",
        "matched_median_plddt": f"{statistics.median(row[2] for row in matched):.3f}",
        "dual_novel_median_length_aa": f"{statistics.median(row[1] for row in dual):.0f}",
        "matched_median_length_aa": f"{statistics.median(row[1] for row in matched):.0f}",
        "plddt_ge_0.60_length_ge100_dual_novel_n": str(
            sum(row[3] for row in confident if row[1] >= 100)
        ),
        "known_fold_matched_false_novelty_control": "NOT_RUN",
    }
    for threshold, denominator, n_dual, rate in threshold_rows:
        stem = f"dual_novel_plddt_ge_{threshold:.2f}"
        summary[f"{stem}_n"] = str(n_dual)
        summary[f"{stem}_denominator"] = str(denominator)
        summary[f"{stem}_pct"] = f"{rate:.1f}"
    for stem, _label, denominator, _n_dual, rate in length_rows:
        summary[f"plddt_ge_0.60_length_{stem}_n"] = str(denominator)
        summary[f"plddt_ge_0.60_length_{stem}_dual_novel_pct"] = f"{rate:.1f}"
    return summary, threshold_rows, length_rows


def write_summary(path: Path, summary: dict[str, str]) -> None:
    ordered = [
        "n_plddt_records",
        "n_plddt_ok",
        "n_plddt_too_long",
    ]
    for threshold in THRESHOLDS:
        stem = f"dual_novel_plddt_ge_{threshold:.2f}"
        ordered.extend((f"{stem}_n", f"{stem}_denominator", f"{stem}_pct"))
    ordered.extend(
        (
            "dual_novel_median_plddt",
            "matched_median_plddt",
            "dual_novel_median_length_aa",
            "matched_median_length_aa",
        )
    )
    for stem in ("lt100", "100_150", "150_250", "250_400", "ge400"):
        ordered.extend(
            (
                f"plddt_ge_0.60_length_{stem}_n",
                f"plddt_ge_0.60_length_{stem}_dual_novel_pct",
            )
        )
    ordered.extend(
        (
            "plddt_ge_0.60_length_ge100_dual_novel_n",
            "known_fold_matched_false_novelty_control",
        )
    )
    with path.open("w") as handle:
        handle.write("metric\tvalue\n")
        handle.write("# RQ-A3 confidence/length sensitivity for candidate novel-fold representatives\n")
        handle.write("# source: best_afdb_tm.tsv + best_pdb_tm.tsv + plddt_shard*.tsv\n")
        handle.write("# dual-novel = best TM < 0.5 against AFDB AND PDB; no hit counts as novel\n")
        for key in ordered:
            handle.write(f"{key}\t{summary[key]}\n")


def draw_figure(
    threshold_rows: list[tuple[float, int, int, float]],
    length_rows: list[tuple[str, str, int, int, float]],
) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.5, 4.2))

    thresholds = [row[0] for row in threshold_rows]
    rates = [row[3] for row in threshold_rows]
    left.plot(thresholds, rates, marker="o", color="#32688e", linewidth=2)
    for threshold, denominator, _n_dual, rate in threshold_rows:
        left.annotate(
            f"n={denominator:,}",
            (threshold, rate),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            fontsize=7,
        )
    left.set_ylim(0, 30)
    left.set_xlabel("minimum mean pLDDT")
    left.set_ylabel("dual-novel representatives (%)")
    left.set_title("(a) Confidence-threshold sensitivity", loc="left", weight="bold", fontsize=10)

    labels = [row[1] for row in length_rows]
    length_rates = [row[4] for row in length_rows]
    bars = right.bar(labels, length_rates, color="#c47a40", edgecolor="black", linewidth=0.4)
    for bar, row in zip(bars, length_rows):
        right.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height(),
            f"{row[4]:.1f}%\nn={row[2]:,}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    right.set_ylim(0, 31)
    right.set_xlabel("protein length (aa), pLDDT≥0.60")
    right.set_ylabel("dual-novel representatives (%)")
    right.set_title("(b) Length sensitivity", loc="left", weight="bold", fontsize=10)

    fig.tight_layout()
    for extension in ("pdf", "png"):
        fig.savefig(HERE / f"fig_fold_validation.{extension}", bbox_inches="tight", dpi=200)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plddt",
        nargs="+",
        type=Path,
        help="one concatenated pLDDT TSV or the individual shard TSVs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    plddt_paths = args.plddt
    if not plddt_paths:
        metagenome_root = os.environ.get(
            "M", "/ibex/scratch/projects/c2014/rob/metagenomics"
        )
        plddt_paths = [
            Path(path)
            for path in sorted(
                glob.glob(f"{metagenome_root}/06_esmfold_novelfold/pdb/plddt_shard*.tsv")
            )
        ]
    if not plddt_paths:
        raise SystemExit("no pLDDT shards found; pass --plddt PATH [PATH ...]")

    afdb = load_tm(HERE / "best_afdb_tm.tsv")
    pdb = load_tm(HERE / "best_pdb_tm.tsv")
    models, statuses = load_plddt(plddt_paths)
    summary, threshold_rows, length_rows = summarise(models, statuses, afdb, pdb)
    write_summary(HERE / "fold_validation_summary.tsv", summary)
    draw_figure(threshold_rows, length_rows)

    for threshold, denominator, n_dual, rate in threshold_rows:
        print(
            f"pLDDT>={threshold:.2f}: n={denominator:,}; "
            f"dual-novel={n_dual:,} ({rate:.1f}%)"
        )
    print("wrote fold_validation_summary.tsv and fig_fold_validation.{pdf,png}")


if __name__ == "__main__":
    main()
