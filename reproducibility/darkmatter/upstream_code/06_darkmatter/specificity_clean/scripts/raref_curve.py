#!/usr/bin/env python3
"""
E-value-corrected analytic rarefaction of the §3.1 absence statistic.

THE PROBLEM THIS SOLVES
  A DIAMOND E-value scales with database size: E ~ K*m*n*exp(-lambda*S). So a FIXED --evalue 1e-3 is
  more permissive against a small catalogue than a large one, which flatters exactly the small
  biome-matched catalogues (atacama 6.2M ... vs mgnify90 717.7M). The raw cross-catalogue ordering
  therefore cannot by itself support a "content beats depth" claim. This script asks the size question
  alone: if MGnify90 were smaller, purely by having fewer sequences, how much absence would we see?

THE METHOD (no subsampled databases are ever built)
  Search ONCE against full MGnify90 at --evalue 1 -k 500, keeping every E_full.
  Retain each subject sequence independently with probability p (DB shrinks to fraction p):
    * E-values in the shrunken DB scale as  E_sub = p * E_full
    * a hit is reported at the 1e-3 rule iff  p*E_full <= 1e-3   <=>   E_full <= 1e-3/p
    * so h(p) = #{subjects of this query with E_full <= 1e-3/p}  are the ones that COULD be reported
    * P(query has no reportable hit) = (1-p)^h(p)
  Expected absence at p = mean over families of (1-p)^h(p), computed separately for the 106,042 dark
  families and the 106,042 length/size-matched annotated control.

WHAT IT DOES AND DOES NOT LICENSE
  It isolates the SIZE effect: subsampling is uniform-random, so composition is held constant by
  construction. A real smaller catalogue differs in size AND content; the gap between the observed
  absence for a real catalogue and this curve at its matched p is the part size cannot explain.
  Reading the curve at p = n_target/717738164 gives the size-only expectation for that depth.

BIASES
  NOTE ON DIRECTION (2026-08-17): the -k cap's bias flipped meaning once this curve became the baseline
  for an observed-vs-expected GAP (see results/content_vs_depth.tsv) rather than a bare absence estimate.
  Truncation inflates the EXPECTATION, and therefore inflates the gap. It is negligible — a capped query
  has h >= 500, and (1-p)^500 ~ 1.3% even at the smallest (atacama) anchor — but do not repeat the older
  "conservative for the absence claim" wording, which held only for the bare-absence framing.
  * -k 500 truncates h(p) for families with >500 hits, which can only overstate absence. Those are
    families that plainly have homologs, so the effect on the curve's informative region is small.
  * the search is at --evalue 1, so p < 1e-3 is not representable (threshold 1e-3/p would exceed 1).
    The grid therefore stops at p = 1e-3 (~718k sequences).

USAGE
  raref_curve.py <hits.m8: qseqid evalue> <FAMILIES_FINAL.txt> <CONTROL_FAMILIES.txt> <out.tsv> [name:N ...]
  Optional trailing name:N pairs add anchors at exact catalogue sizes N (added 2026-09-24 for the SMAG
  leg, spec_6c). With none given the grid, and so the output, is identical to rarefaction_curve.tsv.
"""
import sys, math
from array import array
from bisect import bisect_left

FULL_DB = 717738164          # mgnify90 makedb sequence count (June log, reused in spec_1_global)
EVALUE_RULE = 1e-3           # the paper's absence rule

# p grid: a log sweep plus anchors matched to every comparator catalogue we report, so the curve can be
# read directly against the observed rows in specificity_summary.tsv / desert_specificity_summary.tsv.
ANCHORS = [
    ("atacama",             6193359),
    ("desert_autotroph",    7200701),
    ("arid",               30569562),
    ("antarctic_colddesert",49206550),
    ("mojave",             54361006),
    ("gem",               122978191),
    ("deserts_combined",  147531178),
    ("spire_est",         271000000),
    ("mgnify90_full",     FULL_DB),
]

def build_grid(extra=()):
    pts = {}
    for k in range(0, 31):                      # 1e-3 .. 1 log sweep
        p = 10 ** (-3.0 + 3.0 * k / 30.0)
        if p > 1.0: p = 1.0
        pts.setdefault(round(p, 10), "")
    for name, n in list(ANCHORS) + list(extra):
        p = min(1.0, n / FULL_DB)
        key = round(p, 10)
        pts[key] = (pts.get(key, "") + "," + name).lstrip(",")
    grid = sorted(pts.items())                   # ascending p
    return [p for p, _ in grid], [lab for _, lab in grid]

def main():
    m8, darkf, ctrlf, outf = sys.argv[1:5]
    extra = [(a.rsplit(":", 1)[0], int(a.rsplit(":", 1)[1])) for a in sys.argv[5:]]
    ps, labels = build_grid(extra)
    # thresholds 1e-3/p are DESCENDING in p, so reverse to keep an ascending array for bisect
    order = sorted(range(len(ps)), key=lambda j: EVALUE_RULE / ps[j])
    thr = [EVALUE_RULE / ps[j] for j in order]   # ascending
    NP = len(thr)

    idx, sets = {}, []
    for tag, path in (("dark", darkf), ("ctrl", ctrlf)):
        n0 = len(idx)
        with open(path) as fh:
            for line in fh:
                fam = line.split()[0] if line.strip() else None
                if fam and fam not in idx:
                    idx[fam] = len(idx)
        sets.append((tag, n0, len(idx)))
    N = len(idx)
    sys.stderr.write(f"families indexed: {N} ({sets})\nthresholds: {NP}\n")

    # counts[q*NP + j0] = number of subjects whose E_full first qualifies at threshold index j0.
    # A subject with E_full <= thr[j] for all j >= j0 contributes to h at every j >= j0, so h is the
    # prefix sum of counts over j0 <= j. One bisect + one increment per hit row.
    counts = array('i', bytes(4 * N * NP))
    rows = skipped = 0
    with open(m8) as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) < 2: continue
            q = f[0]
            i = idx.get(q)
            if i is None:
                skipped += 1; continue
            try: e = float(f[-1])
            except ValueError: continue
            j0 = bisect_left(thr, e)
            if j0 >= NP: continue                # E_full > 1.0, unreachable at --evalue 1
            counts[i * NP + j0] += 1
            rows += 1
            if rows % 20_000_000 == 0: sys.stderr.write(f"  {rows//1_000_000}M rows\n")
    sys.stderr.write(f"rows used: {rows}  unmatched-query rows: {skipped}\n")

    acc = {tag: [0.0] * NP for tag, _, _ in sets}
    tot = {tag: (hi - lo) for tag, lo, hi in sets}
    for tag, lo, hi in sets:
        a = acc[tag]
        for i in range(lo, hi):
            base = i * NP
            h = 0
            for j in range(NP):
                h += counts[base + j]           # prefix sum == h(p) at this threshold
                if h == 0:
                    a[j] += 1.0
                else:
                    p = ps[order[j]]
                    a[j] += (1.0 - p) ** h if p < 1.0 else 0.0

    with open(outf, 'w') as out:
        out.write("p\teff_db_seqs\tanchor\tdark_pct_absent_expected\tctrl_pct_absent_expected\n")
        for j in range(NP):
            oj = order[j]
            p = ps[oj]
            out.write("%.6g\t%d\t%s\t%.3f\t%.3f\n" % (
                p, round(p * FULL_DB), labels[oj] or "-",
                100.0 * acc["dark"][j] / tot["dark"],
                100.0 * acc["ctrl"][j] / tot["ctrl"]))
    sys.stderr.write(f"wrote {outf}\n")

if __name__ == "__main__":
    main()
