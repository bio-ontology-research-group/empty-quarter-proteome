#!/usr/bin/env python3
"""Bin table for the ENA deposit: the 12,434 HQ+MQ refined bins (binned metagenome tier) incl. the 5,438 dRep
species representatives (MAG tier). Writes code/ena_deposit/wgs_bins.tsv; idempotent, keeps accession columns.

  srun -c 2 --mem 16G -t 01:00:00 python3 14_build_bin_table.py          # bulk parsing: not on a login node

Per bin (decisions: notes/ENA_DEPOSIT.md, Step 3):
  bin, sample, biosample_accession (soil), run_refs            wgs_biosample_map.tsv, wgs_runs.tsv
  is_representative, drep_cluster                              05_refine/_gtdbtk_r232/rep_ids.txt, _drep Cdb.csv
  completeness, contamination, genome_size, n_contigs, n50     05_refine/_checkm2_refined/quality_report.tsv (CheckM2 1.1.0)
  gtdb_lineage                                                 _gtdbtk_r232 (reps) / _gtdbtk_r232_nonreps/run (others)
  has_5S, has_16S, has_23S, n_trna_aa                          mimag_rrna_trna_clean/results/rrna_trna_per_{rep,bin_nonreps}.tsv
  assembly_quality                                             ENA checklist category from the MIMAG rule
  coverage                                                     length-weighted mean depth over the bin's contigs, from the
                                                               sample's samtools coverage table (ena_wgs/assemblies/coverage/)
                                                               and the DAS_Tool contig2bin map
  fasta                                                        05_refine/_refined_bins_all/<bin>.fa
  ncbi_taxid, ncbi_name                                        filled by 15_resolve_bin_taxa.py
  binned_sample_accession, ena_analysis_binned,                filled by the registration / submission steps
  mag_sample_accession, ena_analysis_mag
Missing inputs (GTDB-Tk or rRNA/tRNA not finished) leave the field empty and are counted in the summary.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

from ena_common import HERE, WORK

REFINE = Path("/ibex/scratch/projects/c2014/rob/metagenomics/05_refine")
MIMAG = Path("/ibex/scratch/projects/c2014/rob/metagenomics/mimag_rrna_trna_clean/results")
QUALITY = REFINE / "_checkm2_refined" / "quality_report.tsv"
CDB = REFINE / "_drep" / "out" / "data_tables" / "Cdb.csv"
REP_IDS = REFINE / "_gtdbtk_r232" / "rep_ids.txt"
GTDB_REPS = REFINE / "_gtdbtk_r232"
GTDB_NONREPS = REFINE / "_gtdbtk_r232_nonreps" / "run"
BINS_DIR = REFINE / "_refined_bins_all"
COVERAGE_DIR = WORK / "assemblies" / "coverage"
TABLE = HERE / "wgs_bins.tsv"
SAMPLES = HERE / "wgs_biosample_map.tsv"
RUNS = HERE / "wgs_runs.tsv"

QUALITY_MANY = "Many fragments with little to no review of assembly other than reporting of standard assembly statistics"
QUALITY_RRNA = ("Multiple fragments where gaps span repetitive regions. Presence of the 23S, 16S, and 5S rRNA genes "
                "and at least 18 tRNAs")
QUALITY_SINGLE = "Single contiguous sequence without gaps or ambiguities with a consensus error rate equivalent to Q50 or better"

KEEP = ["ncbi_taxid", "ncbi_name", "ncbi_level", "binned_sample_accession", "ena_analysis_binned",
        "mag_sample_accession", "ena_analysis_mag"]
COLUMNS = ["bin", "sample", "biosample_accession", "run_refs", "is_representative", "drep_cluster", "completeness",
           "contamination", "genome_size", "n_contigs", "n50", "gtdb_lineage", "has_5S", "has_16S", "has_23S",
           "n_trna_aa", "n_trna_total", "assembly_quality", "coverage", "fasta"] + KEEP


def tsv(path: Path, delimiter: str = "\t") -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def pct(x: str) -> str:
    """ENA pattern: integer or up to two decimals, 100 written as '100'."""
    v = round(float(x), 2)
    return "100" if v >= 100 else f"{v:g}"


def gtdb_lineages() -> dict[str, str]:
    out = {}
    for d in (GTDB_REPS, GTDB_NONREPS):
        for name in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
            p = d / name
            if p.is_file():
                for r in tsv(p):
                    out[r["user_genome"]] = r["classification"]
    return out


def rrna_trna() -> dict[str, dict]:
    out = {}
    for name in ("rrna_trna_per_rep.tsv", "rrna_trna_per_bin_nonreps.tsv"):
        p = MIMAG / name
        if p.is_file():
            for r in tsv(p):
                out[r["rep"]] = r
    return out


def bin_coverage(sample: str, contigs_of_bin: dict[str, list[str]]) -> dict[str, str]:
    """Length-weighted mean depth per bin of one sample from its samtools coverage table."""
    p = COVERAGE_DIR / f"{sample}.samtools_coverage.tsv"
    if not (p.is_file() and p.stat().st_size):
        return {}
    depth: dict[str, tuple[float, float]] = {}
    with open(p) as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        ir, i0, i1, idp = head.index("rname"), head.index("startpos"), head.index("endpos"), head.index("meandepth")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            depth[f[ir]] = (float(f[i1]) - float(f[i0]) + 1, float(f[idp]))
    out = {}
    for b, contigs in contigs_of_bin.items():
        num = den = 0.0
        for c in contigs:
            if c not in depth:
                sys.exit(f"ERROR: contig {c} of bin {b} not in {p.name}")
            ln, dp = depth[c]
            num += dp * ln
            den += ln
        out[b] = f"{num / den:.2f}" if den else ""
    return out


def main() -> None:
    quality = {r["Name"]: r for r in tsv(QUALITY)
               if float(r["Completeness"]) >= 50 and float(r["Contamination"]) <= 10}
    reps = set(REP_IDS.read_text().split())
    cluster = {r["genome"][:-3]: r["secondary_cluster"] for r in tsv(CDB, ",")}
    bs = {r["sample_id"]: r["biosample_accession"] for r in tsv(SAMPLES)}
    refs: dict[str, list[str]] = {}
    for r in tsv(RUNS):
        refs.setdefault(r["sample_id"], []).append(r["ena_run"])
    lineage, rt = gtdb_lineages(), rrna_trna()
    old = {r["bin"]: r for r in tsv(TABLE)} if TABLE.is_file() else {}

    bins = sorted(quality)
    by_sample: dict[str, list[str]] = {}
    for b in bins:
        by_sample.setdefault(b.split("__", 1)[0], []).append(b)

    rows, n_cov, n_tax, n_rt = [], 0, 0, 0
    for sample, sample_bins in by_sample.items():
        if sample not in bs:
            sys.exit(f"ERROR: sample {sample} of {sample_bins[0]} has no soil BioSample")
        c2b = REFINE / sample / f"{sample}_DASTool_contig2bin.tsv"
        contigs_of_bin: dict[str, list[str]] = {}
        with open(c2b) as fh:
            for line in fh:
                contig, short = line.rstrip("\n").split("\t")
                contigs_of_bin.setdefault(f"{sample}__{short}", []).append(contig)
        cov = bin_coverage(sample, {b: contigs_of_bin[b] for b in sample_bins})
        for b in sample_bins:
            if b not in contigs_of_bin:
                sys.exit(f"ERROR: {b} not in {c2b.name}")
            q = quality[b]
            fa = BINS_DIR / f"{b}.fa"
            if not fa.is_file():
                sys.exit(f"ERROR: {fa} missing")
            r = rt.get(b)
            comp, cont = float(q["Completeness"]), float(q["Contamination"])
            if r:
                full = comp >= 90 and cont <= 5 and r["has_5S"] == "1" and r["has_16S"] == "1" and r["has_23S"] == "1" \
                    and int(r["n_tRNA_unique_aa"]) >= 18
                aq = QUALITY_RRNA if full else QUALITY_MANY
            else:
                aq = ""
            row = dict(bin=b, sample=sample, biosample_accession=bs[sample], run_refs=",".join(refs[sample]),
                       is_representative=int(b in reps), drep_cluster=cluster.get(b, ""), completeness=pct(comp),
                       contamination=pct(cont), genome_size=q["Genome_Size"], n_contigs=q["Total_Contigs"],
                       n50=q["Contig_N50"], gtdb_lineage=lineage.get(b, ""), has_5S=r["has_5S"] if r else "",
                       has_16S=r["has_16S"] if r else "", has_23S=r["has_23S"] if r else "",
                       n_trna_aa=r["n_tRNA_unique_aa"] if r else "", n_trna_total=r["n_tRNA_total"] if r else "",
                       assembly_quality=aq, coverage=cov.get(b, ""), fasta=str(fa))
            for k in KEEP:
                row[k] = old.get(b, {}).get(k, "")
            n_cov += bool(row["coverage"]); n_tax += bool(row["gtdb_lineage"]); n_rt += bool(r)
            rows.append(row)
    with open(TABLE, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} bins ({sum(r['is_representative'] for r in rows)} representatives) -> {TABLE.name}; "
          f"with coverage {n_cov}, GTDB lineage {n_tax}, rRNA/tRNA {n_rt}")


if __name__ == "__main__":
    main()
