#!/usr/bin/env python3
"""
Does §1.4's "the families are neither rare nor incidental" survive the presence-rule length bias?

BACKGROUND. `fam_presence.tsv` calls a family present in a sample when ANY member gene reaches CoverM
covered fraction >= 0.5 (fam_8_aggregate.sh). A fixed breadth threshold is easier to satisfy on a short
gene, and `presence_rule_length_test.py` DEMONSTRATES the bias (swapping the rule for one whose length
dependence points the other way reverses the length-occupancy relationship). §1.4 quotes three numbers
from this matrix — median 98/551 samples, only 4 single-sample families, 79.4% present in all five
expeditions (94.8% in at least four) — so all three are length-composition dependent.

THE QUESTION THAT ACTUALLY MATTERS is not whether the headline digits shift (they must) but whether the
CLAIM survives when you compare like with like. The claim is that these families are widespread and
recur across expeditions, which is what rules out relic DNA and reagent contamination. If it holds in
every length stratum — including the longest, where the breadth threshold is HARDEST to satisfy and so
the test is most conservative — the claim stands and only the presentation needs fixing.

RESULT (all 551 samples, 5 expeditions; run 2026-08-17):
    stratum    n_fams  med_samples/551  n_single_sample  %all-5-trips  %>=4-trips
    <60aa       24535       201                4            89.7          97.7
    60-99       26099       107                0            81.9          95.6
    100-149     19720        87                0            78.4          94.9
    150-249     21301        71                0            73.6          92.9
    >=250aa     14387        55                0            67.1          90.8
    ALL         106042       98                4            79.4          94.8

CONCLUSION. The claim holds in every stratum. Even the >=250 aa families — the hardest to call present
under a breadth rule — sit at a median of 55 of 551 samples with 67.1% in all five expeditions and
90.8% in at least four. All four single-sample families are under 60 aa; no family at or above 60 aa is
confined to one sample. So the anti-relic-DNA and anti-contaminant inference is not merely intact, it is
STRENGTHENED: the families least favoured by the detection rule still recur across expeditions.
What must change is presentation, not the finding: the headline median 98 is a composition-weighted
number spanning 201 to 55 across length strata, so quote it with that range or quote the stratified
persistence instead.

USAGE
  s14_persistence_by_length.py <basis_results_dir> <out.tsv>
"""
import sys, os
from statistics import median

STRATA = [(0, 60, "lt60aa"), (60, 100, "60-99"), (100, 150, "100-149"),
          (150, 250, "150-249"), (250, 10**9, "ge250aa")]

def main():
    B, outf = sys.argv[1:3]
    with open(os.path.join(B, "fam_presence.tsv")) as fh:
        cols = fh.readline().rstrip("\n").split("\t")[1:]
    trip = {}
    with open(os.path.join(B, "sample_site_trip.tsv")) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) >= 3: trip[f[0]] = f[2]
    trips = sorted({trip[c] for c in cols})
    ti = {t: i for i, t in enumerate(trips)}
    ctrip = [ti[trip[c]] for c in cols]

    L, cur, n = {}, None, 0
    with open(os.path.join(B, "fam_rep_seq.fasta")) as fh:
        for l in fh:
            if l.startswith(">"):
                if cur is not None: L[cur] = n
                cur = l[1:].split()[0]; n = 0
            else:
                n += len(l.strip())
    if cur is not None: L[cur] = n

    rec = []
    with open(os.path.join(B, "fam_presence.tsv")) as fh:
        fh.readline()
        for l in fh:
            f = l.rstrip("\n").split("\t")
            ns, tr = 0, bytearray(len(trips))
            for j, v in enumerate(f[1:]):
                if v == "1":
                    ns += 1; tr[ctrip[j]] = 1
            rec.append((f[0], ns, sum(tr), L.get(f[0])))

    def row(lab, s, tot):
        return "%s\t%d\t%.1f\t%d\t%d\t%.1f\t%.1f\n" % (
            lab, len(s), 100.0 * len(s) / tot, int(median([r[1] for r in s])),
            sum(1 for r in s if r[1] == 1),
            100.0 * sum(1 for r in s if r[2] == len(trips)) / len(s),
            100.0 * sum(1 for r in s if r[2] >= 4) / len(s))

    allr = [r for r in rec if r[3] is not None]
    with open(outf, "w") as out:
        out.write("length_stratum\tn_families\tpct_universe\tmedian_samples_of_551\t"
                  "n_single_sample\tpct_all_5_trips\tpct_ge4_trips\n")
        for lo, hi, lab in STRATA:
            s = [r for r in allr if lo <= r[3] < hi]
            if s: out.write(row(lab, s, len(allr)))
        out.write(row("ALL_as_in_paper", allr, len(allr)))
    sys.stderr.write(f"wrote {outf} ({len(trips)} expeditions, {len(cols)} samples)\n")

if __name__ == "__main__":
    main()
