#!/usr/bin/env python3
"""One Webin-CLI reads manifest per run from wgs_runs.tsv (paired FASTQ via the symlink farm).

  python 04_generate_manifests.py --study PRJEB125153 [--allow-missing-biosample]

Output: WORK/manifests/<run_alias>.manifest and MANIFEST_INDEX.tsv. Decisions: notes/ENA_DEPOSIT.md.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

from ena_common import HERE, WORK

RUNS = HERE / "wgs_runs.tsv"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", required=True, help="child study accession (PRJEB...)")
    ap.add_argument("--library-source", default="METAGENOMIC", choices=["METAGENOMIC", "GENOMIC"])
    ap.add_argument("--allow-missing-biosample", action="store_true")
    args = ap.parse_args()
    if not re.fullmatch(r"PRJEB\d+", args.study):
        sys.exit(f"ERROR: --study must be a PRJEB accession (got {args.study})")

    with open(RUNS) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    missing = [r["run_alias"] for r in rows if not r["biosample_accession"]]
    if missing and not args.allow_missing_biosample:
        sys.exit(f"ERROR: {len(missing)} runs lack a BioSample accession (register the 106 first, or "
                 f"--allow-missing-biosample to skip them). First: {missing[:5]}")
    unknown = [r["run_alias"] for r in rows if r["instrument_model"] == "UNKNOWN"]
    if unknown:
        sys.exit(f"ERROR: instrument model UNKNOWN for {unknown[:5]} ... fix INSTRUMENT_MODEL in 00_build_runs_table.py")

    out = WORK / "manifests"
    out.mkdir(parents=True, exist_ok=True)
    index = []
    n = 0
    for r in rows:
        if not r["biosample_accession"]:
            continue
        r1 = WORK / "fastq" / r["upload_r1"]
        r2 = WORK / "fastq" / r["upload_r2"]
        for p in (r1, r2):
            if not p.exists():
                sys.exit(f"ERROR: {p} missing (run 00_build_runs_table.py --symlinks)")
        desc = f"Shotgun metagenome of soil sample {r['sample_id']}"
        if r["sample_id"] in ("48Dr2",):          # the only public caveat: the deposited files are not raw
            desc += "; fastp-processed reads (raw R1 file lost)"
        lines = [
            f"STUDY\t{args.study}",
            f"SAMPLE\t{r['biosample_accession']}",
            f"NAME\t{r['run_alias']}",
            f"PLATFORM\t{r['platform']}",
            f"INSTRUMENT\t{r['instrument_model']}",
            "LIBRARY_STRATEGY\tWGS",
            f"LIBRARY_SOURCE\t{args.library_source}",
            "LIBRARY_SELECTION\tRANDOM",
            f"LIBRARY_NAME\t{r['run_alias']}",
            f"LIBRARY_CONSTRUCTION_PROTOCOL\t{r['library_kit']}",
            f"DESCRIPTION\t{desc}",
        ]
        if r["insert_size"]:
            lines.append(f"INSERT_SIZE\t{round(float(r['insert_size']))}")
        lines += [f"FASTQ\t{r1}", f"FASTQ\t{r2}"]
        (out / f"{r['run_alias']}.manifest").write_text("\n".join(lines) + "\n")
        index.append((r["run_alias"], r["sample_id"], r["biosample_accession"], r["instrument_model"]))
        n += 1
    with open(out / "MANIFEST_INDEX.tsv", "w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["run_alias", "sample_id", "biosample_accession", "instrument_model"])
        w.writerows(index)
    print(f"{n} manifests written to {out} (study {args.study}); skipped without BioSample: {len(missing)}")
    if missing:
        print("skipped:", " ".join(missing[:20]), "..." if len(missing) > 20 else "")


if __name__ == "__main__":
    main()
