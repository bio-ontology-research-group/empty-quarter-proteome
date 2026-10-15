#!/usr/bin/env python3
"""eco_0 - per-sample sequencing depth from the fastp reports (covariate for eco_3).

Why: the community blocks in eco_3 are computed on TPM-normalised abundances, but a referee will
ask whether a block that looks environmental is really a depth-of-sequencing gradient (deeper
libraries detect more of the rare dark families, which moves a Bray-Curtis profile). This step
writes the depth table once so eco_3 can condition on it.

Depth = summary.after_filtering.total_bases of each `<sample>_report.json`, i.e. bases surviving
fastp, which is the quantity the assembly and the read mapping actually saw. No table of this
exists in the workbench metadata (sample_metadata.tsv carries paths and flags only), so it is
parsed from the reports themselves.

Idempotent: rewrites results/sample_depth.tsv from the reports on every run.

Usage:
  python3 eco_0_depth.py
"""
import json
import os
import sys
import glob

READS = ("/ibex/scratch/projects/c2014/EmptyQuarter_Data/soil/wgs_assemblies/"
         "processed_reads")
B = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter"
RES = os.path.join(B, "darkfam_ecology_clean", "results")
OUT = os.path.join(RES, "sample_depth.tsv")


def parse_reports(reads_dir):
    rows = []
    for path in sorted(glob.glob(os.path.join(reads_dir, "*_report.json"))):
        sample = os.path.basename(path)[: -len("_report.json")]
        with open(path) as fh:
            rep = json.load(fh)
        aft = rep["summary"]["after_filtering"]
        bef = rep["summary"]["before_filtering"]
        rows.append((sample, int(aft["total_bases"]), int(aft["total_reads"]),
                     int(bef["total_bases"])))
    return rows


def main():
    rows = parse_reports(READS)
    if not rows:
        sys.exit(f"[eco_0] no fastp reports under {READS}")
    bad = [s for s, b, r, _ in rows if b <= 0 or r <= 0]
    if bad:
        sys.exit(f"[eco_0] {len(bad)} reports with non-positive depth: {bad[:5]}")
    with open(OUT, "w") as o:
        o.write("sample\ttotal_bases_after_filtering\ttotal_reads_after_filtering\t"
                "total_bases_before_filtering\n")
        for s, b, r, b0 in rows:
            o.write(f"{s}\t{b}\t{r}\t{b0}\n")
    gb = sorted(b / 1e9 for _, b, _, _ in rows)
    n = len(gb)
    print(f"[eco_0] {n} samples; depth after filtering (Gbp): min {gb[0]:.4f} "
          f"q1 {gb[n // 4]:.2f} median {gb[n // 2]:.2f} q3 {gb[3 * n // 4]:.2f} "
          f"max {gb[-1]:.2f}")
    shallow = sorted((b, s, r) for s, b, r, _ in rows if b < 2e9)
    print(f"[eco_0] {len(shallow)} samples below 2 Gbp: "
          + ", ".join(f"{s} ({r:,} reads)" for _, s, r in shallow))
    print(f"[eco_0] written -> {OUT}")


if __name__ == "__main__":
    main()
