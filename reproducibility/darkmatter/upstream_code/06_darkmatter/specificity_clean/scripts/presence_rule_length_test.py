#!/usr/bin/env python3
"""
Is the length-vs-occupancy relationship in the dark-family matrix BIOLOGY or the PRESENCE RULE?

THE TEST. Presence in `fam_presence.tsv` is defined in fam_8_aggregate.sh as:
    if cf >= 0.5: pres[i] = 1          # family present in a sample if ANY member CDS has
                                       # CoverM Covered Fraction >= 0.5
A fixed coverage-BREADTH threshold is easier to satisfy on a short feature. If that is what drives the
length-occupancy relationship, then swapping the rule for one whose length dependence points the OTHER
way must flip the relationship. So we recompute occupancy from the raw per-sample CoverM output under
three rules and compare:
    cf50     any member Covered Fraction >= 0.5   (the paper's rule; breadth -> favours SHORT)
    cf90     any member Covered Fraction >= 0.9   (stricter breadth -> should favour SHORT harder)
    tpm_any  family TPM > 0, i.e. any read maps anywhere in the family
             (a longer feature is a BIGGER TARGET for "at least one read" -> favours LONG)

RESULT (62 of 551 samples, every 9th, 39 sites; run 2026-08-17) — median rep length by occupancy bin:
    cf50     146 -> 129 -> 104 ->  66 ->  25 aa   (rises in occupancy => SHORTER)
    cf90     147 -> 117 ->  72 ->  35 ->  23 aa   (same direction, stronger: dose-response)
    tpm_any   76 ->  89 ->  97 -> 114 -> 175 aa   (REVERSED => LONGER)
The direction of the relationship is set by the criterion, not by the data. cf90 being stronger than
cf50 is the dose-response that confirms the breadth threshold is the operative variable.

CONCLUSION. The length-occupancy relationship in `fam_presence.tsv` is an artefact of the breadth
threshold, demonstrated rather than inferred. Neither rule recovers "true" occupancy: breadth rules
favour short features, any-read rules favour long ones. Therefore any claim about how widely a dark
family is distributed must be made WITHIN a length stratum, and no unstratified occupancy or
"% endemic" figure derived from this matrix is safe to quote.
NOT affected by this: `fam_abundance.tsv` (summed member TPM), which is length-normalised by
construction and is a depth measure, not a breadth-thresholded presence call.

USAGE
  presence_rule_length_test.py <basis_results_dir> <out.tsv> [sample_stride]
"""
import sys, os, glob
from statistics import median

def main():
    B = sys.argv[1]
    outf = sys.argv[2]
    stride = int(sys.argv[3]) if len(sys.argv) > 3 else 9

    mfam = {}
    with open(os.path.join(B, "member_fam_final.tsv")) as fh:
        for l in fh:
            m, f = l.rstrip("\n").split("\t")
            mfam[m] = f
    L, cur, n = {}, None, 0
    with open(os.path.join(B, "fam_rep_seq.fasta")) as fh:
        for l in fh:
            if l.startswith(">"):
                if cur is not None: L[cur] = n
                cur = l[1:].split()[0]; n = 0
            else:
                n += len(l.strip())
    if cur is not None: L[cur] = n
    fams = sorted(set(mfam.values())); fi = {f: i for i, f in enumerate(fams)}
    site_of = {}
    with open(os.path.join(B, "sample_site_trip.tsv")) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) >= 2: site_of[f[0]] = f[1]

    files = sorted(glob.glob(os.path.join(B, "coverm", "*.tsv")))[::stride]
    RULES = ["cf50", "cf90", "tpm_any"]
    seen = {r: {} for r in RULES}
    used_sites = set()
    for fp in files:
        s = os.path.basename(fp)[:-4]; site = site_of.get(s)
        if site is None: continue
        used_sites.add(site)
        hit = {r: bytearray(len(fams)) for r in RULES}
        with open(fp) as fh:
            fh.readline()
            for line in fh:
                c = line.rstrip("\n").split("\t")
                f = mfam.get(c[0])
                if f is None: continue
                i = fi[f]; t = float(c[1]); cf = float(c[2])
                if cf >= 0.5: hit["cf50"][i] = 1
                if cf >= 0.9: hit["cf90"][i] = 1
                if t > 0:     hit["tpm_any"][i] = 1
        for r in RULES:
            d = seen[r]
            for i, v in enumerate(hit[r]):
                if v: d.setdefault(i, set()).add(site)
    NS = len(used_sites)
    sys.stderr.write(f"samples={len(files)} sites={NS} families={len(fams)}\n")

    def binq(k):
        q = k / NS
        return "0-5%" if q <= .05 else "5-25%" if q <= .25 else "25-50%" if q <= .5 \
               else "50-90%" if q <= .9 else "90-100%"
    ORD = ["0-5%", "5-25%", "25-50%", "50-90%", "90-100%"]
    with open(outf, "w") as out:
        out.write("presence_rule\toccupancy_bin\tn_families\tmedian_rep_len_aa\tpct_lt100aa\t"
                  "n_samples_used\tn_sites_used\n")
        for r in RULES:
            d = seen[r]
            for b in ORD:
                sel = [i for i, ss in d.items() if binq(len(ss)) == b]
                ln = [L[fams[i]] for i in sel if fams[i] in L]
                if not ln: continue
                out.write("%s\t%s\t%d\t%d\t%.1f\t%d\t%d\n" % (
                    r, b, len(sel), int(median(ln)),
                    100.0 * sum(1 for x in ln if x < 100) / len(ln), len(files), NS))
    sys.stderr.write(f"wrote {outf}\n")

if __name__ == "__main__":
    main()
