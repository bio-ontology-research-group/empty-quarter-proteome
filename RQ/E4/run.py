#!/usr/bin/env python3
"""RQ-E4: Are dark families taxon-restricted or broadly distributed across the tree?

Uses fam_taxa_summary.tsv (taxonomic breadth per dark protein family).

Run from repo root:  python3 RQ/E4/run.py
Stdlib + matplotlib only.
"""
import csv, os, statistics
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.getcwd()
INP = os.path.join(REPO, "results/darkmatter_ecology/fam_taxa_summary.tsv")
OUT = os.path.join(REPO, "RQ/E4")
os.makedirs(OUT, exist_ok=True)


def gi(v):
    try:
        return int(v)
    except (ValueError, TypeError):
        return 0


def main():
    rows = list(csv.DictReader(open(INP), delimiter="\t"))
    n = len(rows)
    n_phyla = [gi(r["n_phyla"]) for r in rows]
    n_genera = [gi(r["n_genera"]) for r in rows]
    n_species = [gi(r["n_species"]) for r in rows]
    n_mags = [gi(r["n_mags"]) for r in rows]

    def med(x):
        return statistics.median(x)

    print(f"[n] families = {n}")
    print(f"[n_mags]    median={med(n_mags)} mean={statistics.mean(n_mags):.1f} max={max(n_mags)}")
    print(f"[n_phyla]   median={med(n_phyla)} mean={statistics.mean(n_phyla):.2f} max={max(n_phyla)}")
    print(f"[n_genera]  median={med(n_genera)} mean={statistics.mean(n_genera):.2f} max={max(n_genera)}")
    print(f"[n_species] median={med(n_species)} mean={statistics.mean(n_species):.2f} max={max(n_species)}")

    # --- breadth categories ---------------------------------------------
    single_phylum = sum(1 for x in n_phyla if x <= 1)
    multi_phylum = sum(1 for x in n_phyla if x >= 2)
    single_genus = sum(1 for x in n_genera if x == 1)
    le1_genus = sum(1 for x in n_genera if x <= 1)   # incl. unclassified-at-genus
    multi_genus = sum(1 for x in n_genera if x >= 2)
    print(f"\n[breadth] single-phylum (<=1): {single_phylum} ({100*single_phylum/n:.1f}%)")
    print(f"[breadth] >=2 phyla (broad):    {multi_phylum} ({100*multi_phylum/n:.1f}%)")
    print(f"[breadth] single-genus (==1):   {single_genus} ({100*single_genus/n:.1f}%)")
    print(f"[breadth] <=1 genus:            {le1_genus} ({100*le1_genus/n:.1f}%)")
    print(f"[breadth] >=2 genera (multi):   {multi_genus} ({100*multi_genus/n:.1f}%)")

    phyla_hist = Counter(n_phyla)
    genera_hist = Counter(n_genera)
    print(f"\n[n_phyla hist] {dict(sorted(phyla_hist.items()))}")

    # --- carrier phyla ---------------------------------------------------
    dom = Counter(r["dom_phylum"].strip() for r in rows if r["dom_phylum"].strip())
    print("\n[carrier phyla] (dom_phylum)")
    for k, v in dom.most_common(15):
        print(f"  {k:22s} {v:6d}  ({100*v/n:.1f}%)")

    # --- tables ----------------------------------------------------------
    with open(os.path.join(OUT, "breadth_summary.tsv"), "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"n_families\t{n}\n")
        f.write(f"median_n_phyla\t{med(n_phyla)}\n")
        f.write(f"median_n_genera\t{med(n_genera)}\n")
        f.write(f"median_n_species\t{med(n_species)}\n")
        f.write(f"median_n_mags\t{med(n_mags)}\n")
        f.write(f"single_phylum_n\t{single_phylum}\n")
        f.write(f"single_phylum_pct\t{100*single_phylum/n:.2f}\n")
        f.write(f"multi_phylum_n\t{multi_phylum}\n")
        f.write(f"multi_phylum_pct\t{100*multi_phylum/n:.2f}\n")
        f.write(f"single_genus_n\t{single_genus}\n")
        f.write(f"single_genus_pct\t{100*single_genus/n:.2f}\n")
        f.write(f"multi_genus_n\t{multi_genus}\n")
        f.write(f"multi_genus_pct\t{100*multi_genus/n:.2f}\n")

    with open(os.path.join(OUT, "carrier_phyla.tsv"), "w") as f:
        f.write("dom_phylum\tn_families\tpct\n")
        for k, v in dom.most_common():
            f.write(f"{k}\t{v}\t{100*v/n:.2f}\n")

    with open(os.path.join(OUT, "nphyla_hist.tsv"), "w") as f:
        f.write("n_phyla\tn_families\n")
        for k in sorted(phyla_hist):
            f.write(f"{k}\t{phyla_hist[k]}\n")

    # --- figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # n_phyla histogram
    ax = axes[0]
    maxp = max(n_phyla)
    ax.bar(range(1, maxp + 1), [phyla_hist.get(k, 0) for k in range(1, maxp + 1)],
           color="#3b6ea5")
    ax.set_xlabel("Phyla per family (n_phyla)")
    ax.set_ylabel("Dark families")
    ax.set_title(f"Phylum breadth\nmedian={med(n_phyla)}, >=2 phyla = {100*multi_phylum/n:.1f}%")
    ax.set_yscale("log")

    # n_genera histogram (capped)
    ax = axes[1]
    cap = 15
    gh = Counter(min(x, cap) for x in n_genera)
    xs = list(range(0, cap + 1))
    ax.bar(xs, [gh.get(k, 0) for k in xs], color="#4a8c5a")
    ax.set_xlabel(f"Genera per family (n_genera, capped at {cap}+)")
    ax.set_ylabel("Dark families")
    ax.set_title(f"Genus breadth\nmedian={med(n_genera)}, multi-genus(>=2) = {100*multi_genus/n:.1f}%")
    ax.set_yscale("log")

    # carrier phyla
    ax = axes[2]
    top = dom.most_common(12)
    labels = [k for k, _ in top]
    vals = [v for _, v in top]
    ax.barh(range(len(labels)), vals, color="#a5673b")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Dark families (dominant phylum)")
    ax.set_title("Carrier phyla")
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v}", va="center", fontsize=8)

    fig.suptitle("RQ-E4  Taxonomic breadth of dark protein families", fontsize=13, y=1.01)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_taxbreadth.pdf"))
    fig.savefig(os.path.join(OUT, "fig_taxbreadth.png"), dpi=150)
    print("\n[fig] fig_taxbreadth.pdf/png written")

    with open(os.path.join(OUT, "_summary.txt"), "w") as f:
        f.write(f"n={n}\nmedian_n_phyla={med(n_phyla)}\nmedian_n_genera={med(n_genera)}\n"
                f"median_n_species={med(n_species)}\nmedian_n_mags={med(n_mags)}\n"
                f"single_phylum_pct={100*single_phylum/n:.1f}\n"
                f"multi_phylum_pct={100*multi_phylum/n:.1f}\n"
                f"single_genus_pct={100*single_genus/n:.1f}\n"
                f"multi_genus_pct={100*multi_genus/n:.1f}\n")
        for k, v in dom.most_common(10):
            f.write(f"carrier\t{k}\t{v}\n")


if __name__ == "__main__":
    main()
