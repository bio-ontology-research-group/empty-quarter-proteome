#!/usr/bin/env python3
"""Build wgs_runs.tsv (one row per run) from the four raw Illumina batches; --symlinks (re)builds the
symlink farm WORK/fastq/<run_alias>_R{1,2}.fastq.gz -> raw file (nothing is renamed on disk).

  python 00_build_runs_table.py [--symlinks]

Outputs: wgs_runs.tsv (tracked), WORK/logs/raw_non_cohort.tsv (excluded tokens). Decisions: notes/ENA_DEPOSIT.md.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from ena_common import HERE, WORK

RAW = Path("/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_reads/raw")
COHORT = Path("/ibex/scratch/projects/c2014/rob/metagenomics/metadata/sample_metadata.tsv")  # 551 canonical
BIOSAMPLE_MAP = HERE / "wgs_biosample_map.tsv"
OUT = HERE / "wgs_runs.tsv"

# Instrument from the FASTQ header serial: A01018 = NovaSeq 6000 (KAUST BCL); LL00134 = NovaSeq X series,
# recorded as "Illumina NovaSeq X" (decision 2026-09-03; X vs X Plus not distinguishable from the data).
INSTRUMENT_MODEL = {"A01018": "Illumina NovaSeq 6000", "LL00134": "Illumina NovaSeq X"}
PLATFORM = "ILLUMINA"

# Library kit per sample from the index-code style in the ORIGINAL file names (meta/md5sum/md5sum.txt lists
# for the renamed 2024 batch; the file names themselves for the later batches). Samples with no index code
# (Nov 2025 batch) are assumed QIAseq FX UDI (decision 2026-09-03) and flagged in the manifest description.
KIT_BY_INDEX = {"QIA-UDI": "QIAseq FX DNA Library UDI Kit (QIAGEN)",
                "D7xx-D5xx": "NEBNext Ultra II DNA Library Prep Kit (NEB)"}
KIT_ASSUMED = "QIAseq FX DNA Library UDI Kit (QIAGEN)"
INDEX_RE = re.compile(r"_([A-Za-z]+-UDI\d+-[A-Za-z]+-UDI\d+|D7\d+-D5\d+)_L00\d_R1_001\.fastq\.gz$")
INSERT_SIZES = WORK / "logs" / "insert_sizes.tsv"     # written by 07_insert_size.sbatch (sample, mean, sd, n)

# batch tag -> (glob for R1 files relative to RAW, regex capturing sample id and lane, merged-lanes flag)
LANE_RE = r"_L00(?P<lane>\d)_R1_001\.fastq\.gz$"
DEMUX_RE = r"^M-25-\d+_(?P<sid>.+?)_[A-Za-z]+-UDI\d+-[A-Za-z]+-UDI\d+"
BATCHES = {
    "A_2024":     ("2024/Lane*/version_01/*_R1_001.fastq.gz",
                   r"^(?P<sid>.+?)" + LANE_RE, False),
    "B_20251127": ("20251127/merged/*_merged_R1.fastq.gz",
                   r"^(?P<sid>.+?)_merged_R1\.fastq\.gz$", True),
    "C_20251230": ("20251230_LL00134_0010_A235YJ5LT4/Lane*/version_01/*_R1_001.fastq.gz",
                   DEMUX_RE + LANE_RE, False),
    "D_20260308": ("20260308_LL00134_0014_A23J2WMLT4/merged_reads/*_merged_R1.fastq.gz",
                   r"^(?P<sid>.+?)_merged_R1\.fastq\.gz$", True),
}
D_LANES = "20260308_LL00134_0014_A23J2WMLT4/Lane*/version_01/*_R1_001.fastq.gz"
# suffix for samples sequenced in more than one batch (decision 2026-09-03: <sample>_wgs_2024 / _2026)
BATCH_TAG = {"A_2024": "2024", "B_20251127": "202511", "C_20251230": "202512", "D_20260308": "2026"}

# Known non-cohort tokens sharing a batch dir (controls, QC reruns). Listed, never deposited.
CASE_COLLISION = {"48PRR1": "plate-5 QC rerun of 48PRr1 (uppercase R), excluded"}

# Samples whose raw file is unusable: deposit the fastp-processed pair instead (decision 2026-09-03).
# 48Dr2: raw R1 (batch A) truncated on disk (gzip EOF at exactly 1,125 MiB); no intact copy exists anywhere
# on IBEX (agent search 2026-09-03); the fastp source directory was deleted after processing.
PROCESSED = Path("/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/processed_reads")
PROCESSED_FALLBACK = {"48Dr2": "raw R1 truncated, no intact copy on IBEX; fastp-processed pair deposited instead"}

COLS = ["sample_id", "run_alias", "batch", "trip", "biosample_accession", "instrument_id", "run_number",
        "flowcell", "lanes", "platform", "instrument_model", "library_kit", "kit_source", "insert_size",
        "insert_size_sd", "r1_path", "r2_path", "r1_bytes", "r2_bytes", "upload_r1", "upload_r2", "note"]


def index_styles() -> dict[str, str]:
    """sample id -> index style ('QIA-UDI' | 'D7xx-D5xx') from original file names in md5 lists and dirs."""
    out: dict[str, str] = {}
    names = [line.split()[-1] for p in RAW.glob("*/Lane*/version_01/meta/md5sum/md5sum.txt") for line in p.read_text().splitlines()]
    names += [p.name for p in RAW.glob("*/Lane*/version_01/*_R1_001.fastq.gz")]
    for n in names:
        m = INDEX_RE.search(n)
        sid = re.match(r"^M-\d+-\d+_(.+?)_", n)
        if m and sid:
            style = "D7xx-D5xx" if m.group(1).startswith("D7") else ("QIA-UDI" if m.group(1).startswith("QIA") else m.group(1).split("-")[0] + "-UDI")
            out.setdefault(sid.group(1), style)
    return out


def insert_sizes() -> dict[str, tuple[str, str]]:
    if not INSERT_SIZES.exists():
        return {}
    with open(INSERT_SIZES) as fh:
        return {r["sample"]: (r["insert_mean"], r["insert_sd"]) for r in csv.DictReader(fh, delimiter="\t") if r.get("insert_mean")}


def first_header(path: Path) -> str:
    with gzip.open(path, "rt") as fh:
        return fh.readline().strip()


def trip_of(sid: str) -> str:
    m = re.match(r"^([FSTV]?)\d", sid)
    return {"": "trip1", "T": "trip2", "F": "trip3", "S": "trip4", "V": "trip5"}[m.group(1)]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--symlinks", action="store_true", help="create WORK/fastq/<upload name> -> raw file")
    args = ap.parse_args()

    with open(COHORT) as fh:
        cohort = [r["sample"] for r in csv.DictReader(fh, delimiter="\t")]
    cohort_set = set(cohort)
    if len(cohort_set) != 551:
        sys.exit(f"ERROR: expected 551 canonical samples, found {len(cohort_set)} in {COHORT}")
    with open(BIOSAMPLE_MAP) as fh:
        biosample = {r["sample_id"]: r["biosample_accession"] for r in csv.DictReader(fh, delimiter="\t")}

    # lane sets for batch D (merged_reads concatenates the per-lane files)
    d_lanes: dict[str, set[str]] = defaultdict(set)
    for f in RAW.glob(D_LANES):
        m = re.match(DEMUX_RE + LANE_RE, f.name)
        if m:
            d_lanes[m.group("sid")].add(m.group("lane"))

    rows: list[dict] = []
    excluded: list[dict] = []
    for batch, (pattern, rx, merged) in BATCHES.items():
        for r1 in sorted(RAW.glob(pattern)):
            m = re.match(rx, r1.name)
            if not m:
                excluded.append({"token": r1.name, "batch": batch, "why": "unparsed filename", "path": str(r1)})
                continue
            sid = m.group("sid")
            if sid not in cohort_set:
                why = CASE_COLLISION.get(sid, "not a canonical soil sample (control / rerun / typo id)")
                excluded.append({"token": sid, "batch": batch, "why": why, "path": str(r1)})
                continue
            note = ""
            if sid in PROCESSED_FALLBACK:
                r1 = PROCESSED / f"{sid}_R1.fastq.gz"
                note = PROCESSED_FALLBACK[sid]
            r2 = Path(str(r1).replace("_R1", "_R2"))
            if not r2.exists():
                sys.exit(f"ERROR: no R2 for {r1}")
            lanes = ";".join(sorted(d_lanes.get(sid, []))) if batch == "D_20260308" else (m.groupdict().get("lane") or "merged")
            hdr = first_header(r1).lstrip("@").split(" ")[0].split(":")
            inst, runno, fc = (hdr[0], hdr[1], hdr[2]) if len(hdr) > 3 else ("", "", "")
            rows.append({"sample_id": sid, "batch": batch, "trip": trip_of(sid), "biosample_accession": biosample.get(sid, ""),
                         "instrument_id": inst, "run_number": runno, "flowcell": fc, "lanes": lanes,
                         "platform": PLATFORM, "instrument_model": INSTRUMENT_MODEL.get(inst, "UNKNOWN"),
                         "r1_path": str(r1), "r2_path": str(r2), "r1_bytes": r1.stat().st_size, "r2_bytes": r2.stat().st_size,
                         "note": note or (CASE_COLLISION.get("48PRR1") if sid == "48PRr1" else "")})

    styles = index_styles()
    isz = insert_sizes()
    for r in rows:
        style = styles.get(r["sample_id"])
        if style in KIT_BY_INDEX:
            r["library_kit"], r["kit_source"] = KIT_BY_INDEX[style], f"index code {style}"
        else:
            r["library_kit"], r["kit_source"] = KIT_ASSUMED, "assumed" + (f" (index {style})" if style else " (no index code)")
        r["insert_size"], r["insert_size_sd"] = isz.get(r["sample_id"], ("", ""))

    # run aliases: <id>_wgs, or <id>_wgs_<batchtag> when the sample has more than one batch
    per_sample = Counter(r["sample_id"] for r in rows)
    for r in rows:
        r["run_alias"] = f"{r['sample_id']}_wgs" + (f"_{BATCH_TAG[r['batch']]}" if per_sample[r["sample_id"]] > 1 else "")
        r["upload_r1"] = f"{r['run_alias']}_R1.fastq.gz"
        r["upload_r2"] = f"{r['run_alias']}_R2.fastq.gz"
        if per_sample[r["sample_id"]] > 1:
            r["note"] = "sequenced in two batches on different instruments; one run per batch (policy: Rob)"
    rows.sort(key=lambda r: (r["sample_id"], r["batch"]))

    missing = sorted(cohort_set - {r["sample_id"] for r in rows})
    if missing:
        sys.exit(f"ERROR: {len(missing)} canonical samples without raw reads: {missing}")
    bad_ratio = [r["run_alias"] for r in rows if not 0.8 < r["r2_bytes"] / r["r1_bytes"] < 1.25]
    unknown_inst = sorted({r["instrument_id"] for r in rows if r["instrument_model"] == "UNKNOWN"})

    with open(OUT, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    WORK.joinpath("logs").mkdir(parents=True, exist_ok=True)
    with open(WORK / "logs" / "raw_non_cohort.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["token", "batch", "why", "path"], delimiter="\t")
        w.writeheader()
        w.writerows(excluded)

    if args.symlinks:
        farm = WORK / "fastq"
        farm.mkdir(parents=True, exist_ok=True)
        wanted = {r["upload_r1"] for r in rows} | {r["upload_r2"] for r in rows}
        for stale in farm.iterdir():                      # drop links from an earlier naming scheme
            if stale.is_symlink() and stale.name not in wanted:
                stale.unlink()
        for r in rows:
            for src, name in ((r["r1_path"], r["upload_r1"]), (r["r2_path"], r["upload_r2"])):
                link = farm / name
                if link.is_symlink() or link.exists():
                    link.unlink()
                os.symlink(src, link)
        print(f"symlink farm: {2 * len(rows)} links in {farm}")

    tab: dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        tab[r["batch"]][r["trip"]] += 1
    print(f"runs: {len(rows)} over {len(per_sample)} samples (551 expected) -> {OUT}")
    print("batch\ttrip1\ttrip2\ttrip3\ttrip4\ttrip5\truns\tTB\tinstrument")
    for b in BATCHES:
        br = [r for r in rows if r["batch"] == b]
        tb = sum(r["r1_bytes"] + r["r2_bytes"] for r in br) / 1e12
        inst = sorted({r["instrument_id"] + "/" + r["flowcell"] for r in br})
        print(b, *[tab[b][t] for t in ("trip1", "trip2", "trip3", "trip4", "trip5")], len(br), f"{tb:.2f}", ",".join(inst), sep="\t")
    print(f"total TB: {sum(r['r1_bytes'] + r['r2_bytes'] for r in rows) / 1e12:.2f}")
    print(f"two-batch samples: {sorted(s for s, n in per_sample.items() if n > 1)}")
    print(f"BioSample known: {sum(1 for r in rows if r['biosample_accession'])} runs; to register: "
          f"{len({r['sample_id'] for r in rows if not r['biosample_accession']})} samples")
    print(f"library kit: {dict(Counter(r['kit_source'].split(' (')[0] + ': ' + r['library_kit'].split(' (')[0] for r in rows))}")
    print(f"insert size measured for {sum(1 for r in rows if r['insert_size'])} runs")
    print(f"excluded non-cohort tokens: {len(excluded)} -> {WORK / 'logs' / 'raw_non_cohort.tsv'}")
    if bad_ratio:
        print(f"WARNING R2/R1 size ratio outside 0.8-1.25: {bad_ratio}")
    if unknown_inst:
        print(f"WARNING unknown instrument ids (add to INSTRUMENT_MODEL): {unknown_inst}")


if __name__ == "__main__":
    main()
