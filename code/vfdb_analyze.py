#!/usr/bin/env python3
"""#3 Virulence-factor census from VFDB mmseqs hits — with the dual-use framing built in.

The raw count is huge precisely because VFDB "virulence factors" include generic dual-use
functions (secretion systems, iron/heme uptake, LPS, adhesins) that are ubiquitous in
environmental bacteria. We therefore (a) break hits down by VFDB functional category,
(b) apply identity tiers to show how the signal collapses with stringency, and (c) isolate
genuine pathogenesis categories (exotoxins/toxins). Best hit per query protein only.

In: results/vfdb/vfdb_hits.tsv (query=MAG|MAG::protein, col7=VFDB header w/ category).
Out: results/vfdb/{vf_category_tiers.tsv, vf_toxins.tsv, VFDB.md}
"""
import os, re, csv
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
VF = os.path.join(HERE, "..", "results", "vfdb")
GTDB = [os.path.join(HERE, "..", "workbench", "results", "soil_gtdbtk", x)
        for x in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv")]
HITS = os.path.join(VF, "vfdb_hits.tsv")

CAT_RE = re.compile(r"-\s*([^()\[\]]+?)\s*\(VFC\d+\)\]")     # functional category
VFNAME_RE = re.compile(r"\[([^\[\]]+?)\s*\(VF\d+\)")          # specific VF class
GENE_RE = re.compile(r"\(([^)]+)\)\s+\[")                     # not used; placeholder
TIERS = [40, 50, 70, 90]
TOXIN_KEYWORDS = ("toxin", "Exotoxin", "hemolysin", "haemolysin", "cytolysin",
                  "leukocidin", "Shiga", "anthrax", "diphtheria", "pertussis")


def best_per_query(path):
    """Stream hits, keep best (max bits) per query protein."""
    best = {}
    with open(path) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 7:
                continue
            q = c[0]
            try:
                pid = float(c[2]); bits = float(c[5])
            except ValueError:
                continue
            cur = best.get(q)
            if cur is None or bits > cur[1]:
                best[q] = (pid, bits, c[6])
    return best


def main():
    best = best_per_query(HITS)
    print(f"distinct query proteins with VF hit: {len(best)}")

    # category x tier counts (hits) + MAG sets
    cat_tier_prot = defaultdict(lambda: Counter())          # cat -> tier -> n proteins
    cat_mags = defaultdict(lambda: defaultdict(set))        # cat -> tier -> set(mag)
    tox_rows = []
    all_mag_tier = defaultdict(set)                          # tier -> set(mag) any VF
    for q, (pid, bits, header) in best.items():
        mag = q.split("|")[0]
        m = CAT_RE.search(header)
        cat = m.group(1).strip() if m else "Unclassified"
        vfm = VFNAME_RE.search(header)
        vfname = vfm.group(1).strip() if vfm else ""
        for t in TIERS:
            if pid >= t:
                cat_tier_prot[cat][t] += 1
                cat_mags[cat][t].add(mag)
                all_mag_tier[t].add(mag)
        if any(k.lower() in header.lower() for k in TOXIN_KEYWORDS) and pid >= 50:
            tox_rows.append((mag, round(pid, 1), vfname, cat))

    # category tier table
    with open(os.path.join(VF, "vf_category_tiers.tsv"), "w") as o:
        o.write("vf_category\tprot_id>=40\tprot_id>=50\tprot_id>=70\tprot_id>=90\t"
                "MAGs_id>=40\tMAGs_id>=70\n")
        for cat in sorted(cat_tier_prot, key=lambda c: -cat_tier_prot[c][40]):
            ct = cat_tier_prot[cat]
            o.write(f"{cat}\t{ct[40]}\t{ct[50]}\t{ct[70]}\t{ct[90]}\t"
                    f"{len(cat_mags[cat][40])}\t{len(cat_mags[cat][70])}\n")

    # toxin table (join GTDB phylum)
    tax = {}
    for p in GTDB:
        if os.path.exists(p):
            for r in csv.DictReader(open(p), delimiter="\t"):
                tax[r["user_genome"]] = r["classification"]

    def phy(mag):
        lin = tax.get(mag, "")
        for t in lin.split(";"):
            if t.startswith("p__"):
                return t[3:]
        return ""

    tox_rows.sort(key=lambda r: -r[1])
    with open(os.path.join(VF, "vf_toxins.tsv"), "w") as o:
        o.write("mag\tpident\tvf_class\tcategory\tphylum\n")
        for mag, pid, vfn, cat in tox_rows:
            o.write(f"{mag}\t{pid}\t{vfn}\t{cat}\t{phy(mag)}\n")

    # markdown
    total_mags = 5438
    ml = [f"# Virulence-factor census (VFDB setB) — dual-use framing\n",
          f"mmseqs of {total_mags} dRep-rep proteomes (20.5M proteins) vs VFDB setB (29,804 VF "
          f"proteins). {len(best):,} proteins have a VF hit (id>=40, cov>=70, e<=1e-10). "
          f"**This count is inflated by dual-use genes** — most 'virulence factors' in VFDB are "
          f"secretion/uptake/surface functions ubiquitous in soil bacteria. Identity tiers below "
          f"show how the signal collapses with stringency.\n",
          f"MAGs with >=1 VF hit: id>=40 {len(all_mag_tier[40])}/{total_mags}, "
          f"id>=70 {len(all_mag_tier[70])}, id>=90 {len(all_mag_tier[90])}.\n",
          "## By VFDB functional category (proteins)\n",
          "| category | id>=40 | id>=50 | id>=70 | id>=90 | MAGs(>=40) |",
          "|---|---|---|---|---|---|"]
    for cat in sorted(cat_tier_prot, key=lambda c: -cat_tier_prot[c][40])[:14]:
        ct = cat_tier_prot[cat]
        ml.append(f"| {cat} | {ct[40]} | {ct[50]} | {ct[70]} | {ct[90]} | {len(cat_mags[cat][40])} |")
    ml += [f"\n## Genuine toxins/exotoxins (id>=50): {len(tox_rows)} proteins\n",
           "The pathogenesis-specific slice (exotoxins, hemolysins, etc.) — the part that actually "
           "speaks to 'toxins', not competition/nutrition:\n",
           "| MAG | %id | VF class | phylum |", "|---|---|---|---|"]
    for mag, pid, vfn, cat in tox_rows[:15]:
        ml.append(f"| {mag.split('::')[0]} | {pid} | {vfn[:40]} | {phy(mag)} |")
    open(os.path.join(VF, "VFDB.md"), "w").write("\n".join(ml) + "\n")

    print("MAGs with VF hit: id>=40", len(all_mag_tier[40]), "id>=70", len(all_mag_tier[70]),
          "id>=90", len(all_mag_tier[90]))
    print("top categories (id>=40):")
    for cat in sorted(cat_tier_prot, key=lambda c: -cat_tier_prot[c][40])[:8]:
        ct = cat_tier_prot[cat]
        print(f"  {cat:48s} {ct[40]:7d} -> id>=70 {ct[70]:6d} -> id>=90 {ct[90]:5d}")
    print(f"genuine toxins/exotoxins (id>=50): {len(tox_rows)} proteins")


if __name__ == "__main__":
    main()
