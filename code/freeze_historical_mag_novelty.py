#!/usr/bin/env python3
"""Freeze the historical R226 functional-MAG novelty classes used in Figure F5c.

This is a lightweight deposition helper, not a taxonomy analysis.  It reduces the
two GTDB-Tk summary tables to the exact MAG -> deepest-resolved-rank class needed
by the figure and checks complete closure against dark_fraction.tsv.
"""
import argparse
import csv
from pathlib import Path


def novelty_class(classification):
    ranks = {part[:3]: part[3:] for part in classification.split(";")}
    if ranks.get("s__"):
        return "Known species"
    if ranks.get("g__"):
        return "Novel species"
    if ranks.get("f__"):
        return "Novel genus"
    return "Novel family+"


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gtdb",
        action="append",
        default=[],
        help="GTDB-Tk bac120/ar53 summary; may be supplied twice",
    )
    parser.add_argument(
        "--dark-fraction",
        default=root / "results/picrust2_bias/dark_fraction.tsv",
        type=Path,
    )
    parser.add_argument(
        "--output",
        default=root / "reproducibility/historical_mag_novelty_class.tsv",
        type=Path,
    )
    args = parser.parse_args()
    gtdb_paths = [Path(path) for path in args.gtdb] or [
        root / "workbench/results/soil_gtdbtk/gtdbtk.bac120.summary.tsv",
        root / "workbench/results/soil_gtdbtk/gtdbtk.ar53.summary.tsv",
    ]

    taxonomy = {}
    for path in gtdb_paths:
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                mag = row["user_genome"]
                value = novelty_class(row["classification"])
                if mag in taxonomy and taxonomy[mag] != value:
                    raise ValueError(f"conflicting taxonomy for {mag}")
                taxonomy[mag] = value

    with args.dark_fraction.open() as handle:
        mags = {row["mag"] for row in csv.DictReader(handle, delimiter="\t")}
    missing = sorted(mags - set(taxonomy))
    if missing:
        raise ValueError(f"GTDB summaries miss {len(missing)} functional-snapshot MAGs")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["mag", "novelty_class"])
        for mag in sorted(mags):
            writer.writerow([mag, taxonomy[mag]])
    print(f"wrote {args.output}: {len(mags)} MAGs")


if __name__ == "__main__":
    main()
