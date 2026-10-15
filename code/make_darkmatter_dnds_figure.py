#!/usr/bin/env python3
"""Build the bounded B2 dN/dS figure without the withdrawn comparator.

The input table is a historical pairwise-NG86 screen.  This figure therefore
shows only the recurrent dark-family distribution and its family-size
robustness.  The vertical omega=1 line is a visual reference, not a calibrated
null for this estimator and filtering pipeline.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


FLOORS = (3, 5, 10, 20, 30)


def load_dark_dnds(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["set"] != "dark" or not row["median_dNdS"]:
                continue
            family_id = row["family_id"]
            if family_id in values:
                raise ValueError(f"duplicate dark family in {path}: {family_id}")
            values[family_id] = float(row["median_dNdS"])
    if not values:
        raise ValueError(f"no evaluable dark families in {path}")
    return values


def load_family_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            family_id = row["family_id"]
            if family_id in sizes:
                raise ValueError(f"duplicate family metadata in {path}: {family_id}")
            sizes[family_id] = int(row["n_members"])
    return sizes


def floor_summary(
    dnds: dict[str, float], sizes: dict[str, int]
) -> list[tuple[int, int, float]]:
    missing = sorted(set(dnds) - set(sizes))
    if missing:
        raise ValueError(
            f"{len(missing)} evaluable dark families lack size metadata; "
            f"first={missing[0]}"
        )
    summary = []
    for floor in FLOORS:
        values = [value for family, value in dnds.items() if sizes[family] >= floor]
        summary.append((floor, len(values), float(np.median(values))))
    return summary


def build_figure(dnds_path: Path, metadata_path: Path, figure_dir: Path) -> None:
    dnds = load_dark_dnds(dnds_path)
    sizes = load_family_sizes(metadata_path)
    summary = floor_summary(dnds, sizes)
    values = np.asarray(list(dnds.values()), dtype=float)

    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))

    shown = values[(values >= 0) & (values <= 2)]
    omitted = len(values) - len(shown)
    axes[0].hist(
        shown,
        bins=np.linspace(0, 2, 60),
        density=True,
        color="#332288",
        alpha=0.78,
        label=(
            f"dark families (total n={len(values):,}; "
            f"all-family median={np.median(values):.3f})"
        ),
    )
    axes[0].axvline(1, color="black", linestyle="--", linewidth=1)
    axes[0].text(
        1.02,
        axes[0].get_ylim()[1] * 0.55,
        "$\\omega=1$ reference\n(not a calibrated null)",
        fontsize=8,
        va="top",
    )
    axes[0].set_xlabel("family median pairwise NG86 dN/dS")
    axes[0].set_ylabel("density among values shown (range: 0--2)")
    axes[0].set_title("Recurrent, alignable dark families")
    axes[0].legend(frameon=False, fontsize=8, loc="upper right")
    if omitted:
        axes[0].text(
            0.98,
            0.03,
            f"{omitted} values >2 not shown",
            transform=axes[0].transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
        )

    floors = [row[0] for row in summary]
    counts = [row[1] for row in summary]
    medians = [row[2] for row in summary]
    axes[1].plot(floors, medians, "o-", color="#cc6677", linewidth=2)
    for floor, count, median in summary:
        first = floor == FLOORS[0]
        axes[1].annotate(
            f"n={count:,}",
            (floor, median),
            xytext=((4 if first else 0), 8),
            textcoords="offset points",
            ha=("left" if first else "center"),
            fontsize=8,
        )
    axes[1].set_xticks(floors)
    axes[1].set_ylim(0.35, 0.53)
    axes[1].set_xlabel("minimum total family members")
    axes[1].set_ylabel("median family dN/dS")
    axes[1].set_title("Family-size sensitivity")

    fig.suptitle("Pairwise NG86 screen: bounded descriptive evidence", fontsize=12)
    fig.tight_layout()
    figure_dir.mkdir(parents=True, exist_ok=True)
    for extension in ("pdf", "png"):
        output = figure_dir / f"fig_darkmatter_dnds.{extension}"
        fig.savefig(output, bbox_inches="tight", dpi=200)
        print(f"wrote {output}")
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--dnds",
        type=Path,
        default=root / "results" / "darkmatter_ecology" / "dnds.tsv",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        default=root / "results" / "darkmatter_ecology" / "darkfam_meta.tsv",
    )
    parser.add_argument("--figure-dir", type=Path, default=root / "figures")
    args = parser.parse_args()
    build_figure(args.dnds, args.metadata, args.figure_dir)


if __name__ == "__main__":
    main()
