#!/usr/bin/env python3
"""
Does the genome-level CoverM detection floor create a SIZE bias, the way the gene-level covered-fraction
rule creates a length bias?

WHY ASK. The gene-level rule (`fam_8_aggregate.sh`, covered fraction >= 0.5) is demonstrably
length-biased — see `presence_rule_length_test.py`, where swapping the rule reverses the
length-occupancy relationship. Sections 2.1-2.4 rest on a DIFFERENT mapping run: `coverm genome` v0.7.0,
`-m relative_abundance`, `-p minimap2-sr`, against `prepared/reps.mmi` (the 5,438 refined species
representatives), profiles at
  rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute/<sample>/abund.tsv
`--min-covered-fraction` was never set in any of the 573 recorded invocations, so CoverM's genome-mode
default of 10% applied. CoverM reports a genome as 0 relative abundance below that floor, so
"detected" in those sections effectively means ">=10% of the genome covered". Same mechanism as the gene
case in principle: a 0.27 Mb genome needs ~27 kb of covered bases to clear the floor, a 13.5 Mb genome
needs ~1.35 Mb. Genome sizes span 49-fold, so the question is empirical, not rhetorical.

THE TEST. For every representative, count the samples in which it is detected (relative abundance > 0),
then look at detection frequency against genome size. If small genomes are detected in systematically
more samples, the bias is real and it propagates into M.8's "fraction of species detected", the
">=15 species of the phylum detected" control, and the compartment/gradient tests that consume
`community_by_sample.tsv`. If detection is flat or rises with size, those sections are clean and M.5
only needs the parameter written down.

Genome sizes come from CheckM2: 05_refine/_checkm2_refined/quality_report.tsv (5,394 of the reps) plus
05_refine/_checkm2_backfill/quality_report.tsv (the remaining 44) — exact-name join, no dedup.

USAGE
  mag_detection_size_test.py <out.tsv>
"""
import sys, os, glob
from statistics import median

M = "/ibex/scratch/projects/c2014/rob/metagenomics"
PROF = f"{M}/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"
QR = [f"{M}/05_refine/_checkm2_refined/quality_report.tsv",
      f"{M}/05_refine/_checkm2_backfill/quality_report.tsv"]
BINS = [(0, 1.5e6, "lt1.5Mb"), (1.5e6, 2.5e6, "1.5-2.5Mb"), (2.5e6, 3.5e6, "2.5-3.5Mb"),
        (3.5e6, 5.0e6, "3.5-5Mb"), (5.0e6, 1e12, "ge5Mb")]

def main():
    outf = sys.argv[1]

    size = {}
    for q in QR:
        if not os.path.exists(q): continue
        with open(q) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ci = {k: i for i, k in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= ci["Genome_Size"]: continue
                nm = f[ci["Name"]]
                if nm not in size:
                    try: size[nm] = int(float(f[ci["Genome_Size"]]))
                    except ValueError: pass
    sys.stderr.write(f"checkm2 genome sizes loaded: {len(size)}\n")

    det, nsamp = {}, 0
    for d in sorted(glob.glob(f"{PROF}/*/abund.tsv")):
        nsamp += 1
        with open(d) as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < 2 or c[0] == "unmapped": continue
                try: v = float(c[1])
                except ValueError: continue
                if c[0] not in det: det[c[0]] = 0
                if v > 0: det[c[0]] += 1
    sys.stderr.write(f"profiles read: {nsamp}   representatives seen: {len(det)}\n")

    matched = {g: n for g, n in det.items() if g in size}
    sys.stderr.write(f"reps with a CheckM2 size: {len(matched)} of {len(det)}\n")

    with open(outf, "w") as out:
        out.write("size_bin\tn_reps\tmedian_genome_size_bp\tmedian_samples_detected\t"
                  "pct_of_samples_detected\tn_never_detected\tn_samples_total\n")
        for lo, hi, lab in BINS:
            sel = [(g, n) for g, n in matched.items() if lo <= size[g] < hi]
            if not sel: continue
            out.write("%s\t%d\t%d\t%d\t%.1f\t%d\t%d\n" % (
                lab, len(sel), int(median([size[g] for g, _ in sel])),
                int(median([n for _, n in sel])),
                100.0 * median([n for _, n in sel]) / nsamp,
                sum(1 for _, n in sel if n == 0), nsamp))
        allv = list(matched.values())
        out.write("ALL\t%d\t%d\t%d\t%.1f\t%d\t%d\n" % (
            len(matched), int(median([size[g] for g in matched])), int(median(allv)),
            100.0 * median(allv) / nsamp, sum(1 for n in allv if n == 0), nsamp))
    sys.stderr.write(f"wrote {outf}\n")

if __name__ == "__main__":
    main()
