#!/usr/bin/env python3
"""Webin-CLI genome manifests for the bins (tier binned) and the MAGs (tier mag) from wgs_bins.tsv.

  python3 18_build_bin_manifests.py --tier binned --study PRJEB125153
  python3 18_build_bin_manifests.py --tier mag    --study PRJEB125153
  then: MANIFEST_DIR=.../ena_wgs/bins/manifests_<tier> LOGDIR=.../ena_wgs/bins/logs MODE=validate|submit
        sbatch --array=1-N%k 13_webin_cli_genome.sbatch
  and:  python3 18_build_bin_manifests.py --tier <tier> --fold     # ERZ from logs/genome_submit_status.tsv

Fields (decisions: notes/ENA_DEPOSIT.md, Step 3): STUDY; SAMPLE = the tier's derived BioSample (binned_sample_accession
or mag_sample_accession, registered by 16); ASSEMBLYNAME = bin name (tier mag: <bin>_MAG, since study+sample+name
must be unique per assembly); ASSEMBLY_TYPE "binned metagenome" | "Metagenome-Assembled Genome (MAG)"; COVERAGE =
length-weighted mean read depth over the bin's contigs; PROGRAM / PLATFORM / MOLECULETYPE as the primary assemblies;
RUN_REF = all runs of the sample; FASTA = ena_wgs/bins/fasta/<bin>.fa.gz (17_gzip_bins.sbatch). No DESCRIPTION: the
sample carries the taxon, quality and lineage.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

from ena_common import HERE, WORK

TABLE = HERE / "wgs_bins.tsv"
OUT = WORK / "bins"
PROGRAM = "MEGAHIT v1.2.9"
PLATFORM = "ILLUMINA"
MOLECULETYPE = "genomic DNA"
TYPE = {"binned": "binned metagenome", "mag": "Metagenome-Assembled Genome (MAG)"}
ACC = {"binned": "binned_sample_accession", "mag": "mag_sample_accession"}
ANALYSIS = {"binned": "ena_analysis_binned", "mag": "ena_analysis_mag"}


def rows_and_cols() -> tuple[list[dict], list[str]]:
    with open(TABLE) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return list(r), list(r.fieldnames)


def name_of(row: dict, tier: str) -> str:
    return row["bin"] if tier == "binned" else f"{row['bin']}_MAG"


def build(tier: str, study: str) -> None:
    if not re.fullmatch(r"PRJEB\d+", study):
        sys.exit(f"ERROR: --study must be a PRJEB accession (got {study})")
    rows, _ = rows_and_cols()
    sel = [r for r in rows if r[ACC[tier]] and (tier == "binned" or r["is_representative"] == "1")]
    mdir = OUT / f"manifests_{tier}"
    mdir.mkdir(parents=True, exist_ok=True)
    for r in sel:
        name = name_of(r, tier)
        if len(name) > 50 or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9 _#\-.]*", name):
            sys.exit(f"ERROR: assembly name {name!r} violates the ENA name rule")
        gz = OUT / "fasta" / f"{r['bin']}.fa.gz"
        if not (gz.is_file() and gz.stat().st_size):
            sys.exit(f"ERROR: {gz} missing (17_gzip_bins.sbatch)")
        if not r["coverage"]:
            sys.exit(f"ERROR: {r['bin']} has no coverage")
        manifest = [("STUDY", study), ("SAMPLE", r[ACC[tier]]), ("ASSEMBLYNAME", name), ("ASSEMBLY_TYPE", TYPE[tier]),
                    ("COVERAGE", r["coverage"]), ("PROGRAM", PROGRAM), ("PLATFORM", PLATFORM),
                    ("MOLECULETYPE", MOLECULETYPE), ("RUN_REF", r["run_refs"]), ("FASTA", str(gz))]
        if tier == "mag" and r["n_contigs"] == "1":
            # ENA MAG rule: a contig-level assembly needs >= 2 sequences; a single highly complete sequence is
            # submitted at chromosome level with a chromosome list (name "1", topology linear: circularity not asserted)
            import gzip
            with gzip.open(gz, "rt") as fh:
                contig = next(line for line in fh if line.startswith(">"))[1:].split()[0]
            cl = OUT / "chromlists" / f"{name}.chromlist.gz"
            cl.parent.mkdir(exist_ok=True)
            with gzip.open(cl, "wt") as fh:
                fh.write(f"{contig}\t1\tLinear-Chromosome\n")
            manifest.append(("CHROMOSOME_LIST", str(cl)))
        (mdir / f"{name}.manifest").write_text("".join(f"{k}\t{v}\n" for k, v in manifest))
    skipped = sum(1 for r in rows if not r[ACC[tier]] and (tier == "binned" or r["is_representative"] == "1"))
    print(f"{len(sel)} {tier} manifests -> {mdir}; without a {tier} sample accession (not registered yet): {skipped}")


def fold(tier: str) -> None:
    """ERZ per bin: ENA's live analysis list (logs/ena_live_analyses.tsv: alias webin-genome-<name>, id — the
    authority; a Webin-CLI run can register an analysis and still exit with an FTP error, so the status log alone
    under-counts) completed by the Webin-CLI status log."""
    erz = {}
    live = OUT / "logs" / "ena_live_analyses.tsv"
    if live.is_file():
        for line in live.read_text().splitlines():
            alias, acc = line.split("\t")
            if alias.startswith("webin-genome-"):
                erz[alias[len("webin-genome-"):]] = acc
    log = OUT / "logs" / "genome_submit_status.tsv"
    if log.is_file():
        with open(log) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["status"] == "OK" and r["erz"]:
                    erz.setdefault(r["alias"], r["erz"])
    if not erz:
        sys.exit(f"ERROR: neither {live} nor {log} holds accessions")
    # ENA's "genome assembly has been accessioned" emails (saved verbatim as receipts/emails/ena_email_<date>.txt):
    # pipe table ASSEMBLY_NAME | ASSEMBLY_ACC | STUDY_ID | SAMPLE_ID | CONTIG_ACC | SCAFFOLD_ACC | CHROMOSOME_ACC.
    # For MAGs the ASSEMBLY_ACC is the GCA (the citable accession) and the sequence accessions are the contig range
    # or, for a chromosome-level MAG, the chromosome accession. ERZ rows (primaries, binned) add nothing to the table.
    gca, seqs = {}, {}
    for p in sorted((OUT.parent / "receipts" / "emails").glob("ena_email_*.txt")):
        for line in p.read_text().splitlines():
            f = [x.strip() for x in line.split("|")]
            if len(f) >= 7 and f[1].startswith("GCA_"):
                gca[f[0]] = f[1]
                seqs[f[0]] = f[4] or f[6] or f[5]
    rows, cols = rows_and_cols()
    if tier == "mag":
        for c in ("gca_mag", "sequence_acc_mag"):
            if c not in cols:
                cols.append(c)
    for r in rows:
        r[ANALYSIS[tier]] = erz.get(name_of(r, tier), r[ANALYSIS[tier]])
        if tier == "mag":
            r["gca_mag"] = gca.get(name_of(r, tier), r.get("gca_mag", ""))
            r["sequence_acc_mag"] = seqs.get(name_of(r, tier), r.get("sequence_acc_mag", ""))
    with open(TABLE, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    done = sum(bool(r[ANALYSIS[tier]]) for r in rows)
    print(f"{tier}: {done} bins carry an ERZ accession"
          + (f"; {sum(bool(r['gca_mag']) for r in rows)} carry a GCA from the accession emails" if tier == "mag" else ""))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", required=True, choices=["binned", "mag"])
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--study")
    g.add_argument("--fold", action="store_true")
    a = ap.parse_args()
    fold(a.tier) if a.fold else build(a.tier, a.study)


if __name__ == "__main__":
    main()
