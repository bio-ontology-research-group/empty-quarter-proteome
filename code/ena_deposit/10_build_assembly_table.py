#!/usr/bin/env python3
"""Primary metagenome assemblies (549) under PRJEB125153: table + Webin-CLI genome manifests.

Phases (run in this order; each is idempotent):
  python 10_build_assembly_table.py --prepare      # WORK/assemblies/sample_list.tsv: sample, contig FASTA, BAM
  sbatch --array=1-N%16 11_gzip_assemblies.sbatch  # fasta/<s>_assembly.fa.gz + stats/<s>.tsv (n_contigs, total_bp)
  sbatch --array=1-N%16 12_coverage_bam.sbatch     # coverage/<s>.samtools_coverage.tsv (per-contig mean depth)
  python 10_build_assembly_table.py --study PRJEB125153   # wgs_assemblies.tsv + manifests/<s>_assembly.manifest
  MODE=validate sbatch --array=1-N%16 13_webin_cli_genome.sbatch
  MODE=submit   sbatch --array=1-N%12 13_webin_cli_genome.sbatch   # IRREVERSIBLE, on her word
  python 10_build_assembly_table.py --fold         # ena_analysis (ERZ) from logs/genome_submit_status.tsv

Inputs: the organized data tree EmptyQuarter_Data/processed/sequencing-results/sites/<trip>/site<N>/<sample>/wgs/
(trip and site from the workbench sample_metadata.tsv): assembly/<sample>_contigs.fa = the MEGAHIT assembly.
Coverage: the SemiBin2 self-supervised mapping (BWA-MEM, the sample's own reads on its own assembly), the only
mapping that exists for all 549 assemblies; it lives in soil/wgs_assemblies/binning/<sample>/.
Fields (decisions: notes/ENA_DEPOSIT.md): SAMPLE = soil BioSample (wgs_biosample_map.tsv); RUN_REF = all runs of the
sample (wgs_runs.tsv, comma-separated for the dual-batch samples); COVERAGE = length-weighted mean read depth over
all contigs, sum(depth_i * len_i) / sum(len_i), from samtools coverage on that BAM; PROGRAM / PLATFORM /
MOLECULETYPE fixed; ASSEMBLYNAME <sample>_assembly; no MINGAPLENGTH (contig level, no gaps).
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from ena_common import HERE, WORK

PROCESSED = Path("/ibex/scratch/projects/c2014/EmptyQuarter_Data/processed/sequencing-results/sites")
SAMPLE_META = Path("/ibex/scratch/projects/c2014/rob/metagenomics/metadata/sample_metadata.tsv")
CONTIGS = PROCESSED / "{trip}" / "site{site}" / "{s}" / "wgs" / "assembly" / "{s}_contigs.fa"
BAM = Path("/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/binning") / "{s}" / "{s}_metaG.reads.sorted.bam"
OUT = WORK / "assemblies"
TABLE = HERE / "wgs_assemblies.tsv"
RUNS = HERE / "wgs_runs.tsv"
SAMPLES = HERE / "wgs_biosample_map.tsv"

PROGRAM = "MEGAHIT v1.2.9"
PLATFORM = "ILLUMINA"
MOLECULETYPE = "genomic DNA"
MIN_CONTIG = 500
COLUMNS = ["sample_id", "assembly_name", "biosample_accession", "run_refs", "n_contigs", "total_bp", "coverage",
           "contigs_src", "fasta_gz", "ena_analysis"]


def cohort() -> list[str]:
    with open(SAMPLES) as fh:
        return sorted({r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")})


def run_refs() -> dict[str, list[str]]:
    refs: dict[str, list[str]] = {}
    with open(RUNS) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if not r["ena_run"]:
                sys.exit(f"ERROR: run {r['run_alias']} has no ENA run accession in wgs_runs.tsv")
            refs.setdefault(r["sample_id"], []).append(r["ena_run"])
    return refs


def biosamples() -> dict[str, str]:
    with open(SAMPLES) as fh:
        return {r["sample_id"]: r["biosample_accession"] for r in csv.DictReader(fh, delimiter="\t")}


def sample_paths() -> dict[str, tuple[Path, Path]]:
    """sample -> (contig FASTA in the organized tree, SemiBin self-supervised BAM)."""
    out = {}
    with open(SAMPLE_META) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = r["sample"]
            out[s] = (Path(str(CONTIGS).format(trip=r["trip"], site=r["site"], s=s)), Path(str(BAM).format(s=s)))
    return out


def coverage(s: str) -> str:
    """Length-weighted mean depth from samtools coverage (#rname startpos endpos numreads covbases coverage
    meandepth ...): sum(meandepth_i * len_i) / sum(len_i)."""
    cov = OUT / "coverage" / f"{s}.samtools_coverage.tsv"
    if not (cov.is_file() and cov.stat().st_size):
        sys.exit(f"ERROR: {cov} missing; run 12_coverage_bam.sbatch")
    num = den = 0.0
    with open(cov) as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        i0, i1, idp = head.index("startpos"), head.index("endpos"), head.index("meandepth")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ln = float(f[i1]) - float(f[i0]) + 1
            num += float(f[idp]) * ln
            den += ln
    if den <= 0:
        sys.exit(f"ERROR: zero total length in {cov}")
    return f"{num / den:.2f}"


def stats(s: str) -> tuple[int, int]:
    p = OUT / "stats" / f"{s}.tsv"
    if not p.is_file():
        sys.exit(f"ERROR: {p} missing; run 11_gzip_assemblies.sbatch")
    n, bp = open(p).read().split()
    return int(n), int(bp)


def prepare() -> None:
    """WORK/assemblies/sample_list.tsv: sample, contig FASTA (organized tree), BAM — one line per cohort sample
    that has a non-empty contig file. Every listed sample must have its BAM (no silent fallback)."""
    OUT.mkdir(parents=True, exist_ok=True)
    paths = sample_paths()
    have, skipped = [], []
    for s in cohort():
        if s not in paths:
            sys.exit(f"ERROR: {s} is not in {SAMPLE_META}")
        fa, bam = paths[s]
        if not (fa.is_file() and fa.stat().st_size):
            skipped.append(s)
            continue
        if not bam.is_file():
            sys.exit(f"ERROR: {s}: no self-supervised BAM at {bam}; no coverage source")
        have.append(f"{s}\t{fa}\t{bam}")
    (OUT / "sample_list.tsv").write_text("\n".join(have) + "\n")
    print(f"{len(have)} assemblies listed in sample_list.tsv (skipped, no contig file: {skipped})")


def build(study: str, partial: bool = False) -> None:
    """--partial: skip samples whose coverage table is not there yet (validation can start on the rest);
    the table then holds only the built rows and is completed by a later full run."""
    if not re.fullmatch(r"PRJEB\d+", study):
        sys.exit(f"ERROR: --study must be a PRJEB accession (got {study})")
    refs, bs = run_refs(), biosamples()
    listed = [line.split("\t") for line in (OUT / "sample_list.tsv").read_text().splitlines() if line]
    (OUT / "manifests").mkdir(exist_ok=True)
    old = {}
    if TABLE.is_file():
        with open(TABLE) as fh:
            old = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    rows, skipped = [], []
    for s, fa, _bam in listed:
        cov_file = OUT / "coverage" / f"{s}.samtools_coverage.tsv"
        if partial and not (cov_file.is_file() and cov_file.stat().st_size):
            skipped.append(s)
            continue
        n, bp = stats(s)
        if n < 2:
            sys.exit(f"ERROR: {s} has {n} contigs; a contig-level assembly needs at least 2")
        cov = coverage(s)
        gz = OUT / "fasta" / f"{s}_assembly.fa.gz"
        if not (gz.is_file() and gz.stat().st_size):
            sys.exit(f"ERROR: {gz} missing; run 11_gzip_assemblies.sbatch")
        name = f"{s}_assembly"
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9 _#\-.]*", name) or len(name) > 50:
            sys.exit(f"ERROR: assembly name {name!r} violates the ENA name rule")
        row = dict(sample_id=s, assembly_name=name, biosample_accession=bs[s], run_refs=",".join(refs[s]),
                   n_contigs=n, total_bp=bp, coverage=cov, contigs_src=fa, fasta_gz=str(gz),
                   ena_analysis=old.get(s, {}).get("ena_analysis", ""))
        rows.append(row)
        manifest = [("STUDY", study), ("SAMPLE", row["biosample_accession"]), ("ASSEMBLYNAME", name),
                    ("ASSEMBLY_TYPE", "primary metagenome"), ("COVERAGE", cov), ("PROGRAM", PROGRAM),
                    ("PLATFORM", PLATFORM), ("MOLECULETYPE", MOLECULETYPE), ("RUN_REF", row["run_refs"]),
                    ("DESCRIPTION", f"Primary metagenome assembly of Rub' al-Khali desert soil sample {s}, "
                                    f"contigs >= {MIN_CONTIG} bp"),
                    ("FASTA", str(gz))]
        (OUT / "manifests" / f"{name}.manifest").write_text("".join(f"{k}\t{v}\n" for k, v in manifest))
    with open(TABLE, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} rows -> {TABLE.name}; manifests in {OUT / 'manifests'}"
          + (f"; skipped (no coverage yet): {len(skipped)}" if partial else ""))


def fold() -> None:
    log = WORK / "assemblies" / "logs" / "genome_submit_status.tsv"
    if not log.is_file():
        sys.exit(f"ERROR: {log} missing (nothing submitted yet)")
    erz = {}
    with open(log) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["status"] == "OK" and r["erz"]:
                erz[r["alias"]] = r["erz"]
    with open(TABLE) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows, cols = list(reader), list(reader.fieldnames)
    for r in rows:
        r["ena_analysis"] = erz.get(r["assembly_name"], r["ena_analysis"])
    with open(TABLE, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    done = sum(bool(r["ena_analysis"]) for r in rows)
    print(f"{done}/{len(rows)} assemblies carry an ERZ accession; missing: {[r['sample_id'] for r in rows if not r['ena_analysis']][:10]}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--prepare", action="store_true")
    g.add_argument("--study", help="build table + manifests for this study accession")
    g.add_argument("--fold", action="store_true")
    ap.add_argument("--partial", action="store_true", help="with --study: skip samples whose coverage is not ready")
    a = ap.parse_args()
    if a.prepare:
        prepare()
    elif a.fold:
        fold()
    else:
        build(a.study, a.partial)


if __name__ == "__main__":
    main()
