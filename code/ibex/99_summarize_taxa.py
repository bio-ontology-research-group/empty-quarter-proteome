#!/usr/bin/env python3
"""Collapse family->member-MAG pairs to a per-family GTDB taxonomy summary.

Usage: 99_summarize_taxa.py fam_mag.tsv gtdbtk.bac120.summary.tsv gtdbtk.ar53.summary.tsv out.tsv

Member MAG names carry a self_/single_ qualifier (e.g. 36Dr1__self_SemiBin_23) and DAS_Tool
'_sub' suffixes that GTDB-Tk user_genome ids (36Dr1__SemiBin_23) lack; norm() reconciles them.
Per family we report taxonomic breadth (n_phyla/genera/species), the dominant phylum/genus
weighted by member count, and the single most common full lineage.
"""
import sys, csv
from collections import defaultdict, Counter

fam_mag, gbac, gar, outp = sys.argv[1:5]


def norm(m):
    m = m.replace("__self_SemiBin", "__SemiBin").replace("__single_SemiBin", "__SemiBin")
    return m[:-4] if m.endswith("_sub") else m


def rank(lin, pre):
    for t in lin.split(";"):
        if t.startswith(pre):
            return t[3:]
    return ""


tax = {}
for g in (gbac, gar):
    with open(g) as fh:
        r = csv.reader(fh, delimiter="\t")
        next(r, None)
        for row in r:
            if len(row) >= 2:
                tax[row[0]] = row[1]
print(f"loaded {len(tax)} GTDB classifications", flush=True)

fam_mags = defaultdict(Counter)
with open(fam_mag) as fh:
    for line in fh:
        fam, mag = line.rstrip("\n").split("\t")
        fam_mags[fam][mag] += 1
print(f"aggregated {len(fam_mags)} families", flush=True)

matched = total = 0
with open(outp, "w") as o:
    o.write("family_id\tn_members\tn_mags\tn_classified_members\tn_phyla\tdom_phylum\t"
            "n_genera\tdom_genus\tn_species\ttop_lineage\ttop_lineage_n\n")
    for fam, mags in fam_mags.items():
        n_members = sum(mags.values())
        phyla, genera, species, lin_c = Counter(), Counter(), Counter(), Counter()
        nclass = 0
        for mag, c in mags.items():
            lin = tax.get(mag) or tax.get(norm(mag))  # raw (refined namespace) then normalized fallback
            total += c
            if not lin:
                continue
            matched += c
            nclass += c
            ph, ge, sp = rank(lin, "p__"), rank(lin, "g__"), rank(lin, "s__")
            if ph:
                phyla[ph] += c
            if ge:
                genera[ge] += c
            if sp:
                species[sp] += c
            lin_c[lin] += c
        domp = phyla.most_common(1)[0][0] if phyla else ""
        domg = genera.most_common(1)[0][0] if genera else ""
        topl = lin_c.most_common(1)[0] if lin_c else ("", 0)
        o.write(f"{fam}\t{n_members}\t{len(mags)}\t{nclass}\t{len(phyla)}\t{domp}\t"
                f"{len(genera)}\t{domg}\t{len(species)}\t{topl[0]}\t{topl[1]}\n")
print(f"member->GTDB match rate: {matched}/{total} ({100*matched/max(total,1):.1f}%)", flush=True)
