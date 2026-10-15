#!/usr/bin/env python3
"""Step 1 -- compatible-solute (osmolyte) synthesis census on the 5,438-representative catalogue.

First proper osmolyte census: the historical panel (scripts/82_metabolism_marker_census.sh,
15,469 MAGs, 2026-06-15) covered only otsA/otsB + ectC/ectD. This panel adds the remaining
pathways the manuscript names (glucosylglycerol, mannosylglycerate) plus the alternative
trehalose routes. SYNTHESIS genes only -- no transporters (dual-use, weak osmotic evidence).
Every KO verified against databases/dram/kofam_ko_list.tsv (2026-08-13).

KO table: 04_functional/repko_clean/results/rep_ko_v2.tsv (5,438 rows, 0 unannotated).

Pathway capacity flags (primary = committed synthase present; `_strict` = full route):
  trehalose        otsA K00697 or fused K16055; strict adds otsB K01087 (or K16055 alone)
                   `trehalose_any` also accepts treS K13057 or the treY K06044 + treZ K01236 route
  ectoine          ectC K06720 (committed synthase); ectD K10674 reported separately (hydroxyectoine)
  glucosylglycerol ggpS K03692; strict adds ggpP K05978
  mannosylglycerate mpgS K05947 or single-step mgs K18818; strict = K18818 or (K05947 + mpgP K07026)

Out (results/): marker_prevalence_5438.tsv , pathway_prevalence_5438.tsv , rep_traits.tsv
(rep_traits.tsv has the same layout as uv_survival_clean's -- consumed unchanged by steps 2-3b).
"""
import csv, os, collections

M    = "/ibex/scratch/projects/c2014/rob/metagenomics"
RES  = f"{M}/06_metabolism/osmolyte_clean/results"
REPKO = f"{M}/04_functional/repko_clean/results/rep_ko_v2.tsv"
GTDB  = [f"{M}/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv",
         f"{M}/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"]
os.makedirs(RES, exist_ok=True)

KO_GENE = {
    "K00697": "otsA_trehalose6P_synthase",      "K01087": "otsB_trehalose6P_phosphatase",
    "K16055": "otsAB_fused",                    "K13057": "treS_trehalose_synthase",
    "K06044": "treY_maltooligosyl_synthase",    "K01236": "treZ_trehalohydrolase",
    "K06720": "ectC_ectoine_synthase",          "K10674": "ectD_ectoine_hydroxylase",
    "K03692": "ggpS_glucosylglycerolP_synthase","K05978": "ggpP_glucosylglycerol_phosphatase",
    "K05947": "mpgS_mannosyl3PG_synthase",      "K18818": "mgs_mannosylglycerate_synthase",
    "K07026": "mpgP_mannosyl3PG_phosphatase",
}

def flags(k):
    """pathway capacity flags from one representative's KO set"""
    tre_ots    = bool({"K00697", "K16055"} & k)
    tre_strict = "K16055" in k or {"K00697", "K01087"} <= k
    tre_any    = tre_ots or "K13057" in k or {"K06044", "K01236"} <= k
    gg         = "K03692" in k
    gg_strict  = {"K03692", "K05978"} <= k
    mg         = bool({"K05947", "K18818"} & k)
    mg_strict  = "K18818" in k or {"K05947", "K07026"} <= k
    return collections.OrderedDict([
        ("trehalose", tre_ots), ("trehalose_strict", tre_strict), ("trehalose_any", tre_any),
        ("ectoine", "K06720" in k), ("hydroxyectoine", "K10674" in k),
        ("glucosylglycerol", gg), ("glucosylglycerol_strict", gg_strict),
        ("mannosylglycerate", mg), ("mannosylglycerate_strict", mg_strict),
        ("any_osmolyte", tre_any or "K06720" in k or gg or mg),
    ])

rep_phylum = {}
for p in GTDB:
    for r in csv.DictReader(open(p), delimiter="\t"):
        rep_phylum[r["user_genome"]] = next(
            (t[3:] for t in r["classification"].split(";") if t.startswith("p__")), "")

rep_kos = {}
for line in open(REPKO):
    c = line.rstrip("\n").split("\t")
    if c and c[0]:
        rep_kos[c[0]] = set(c[1].split(",")) & set(KO_GENE) if len(c) > 1 and c[1] else set()
N = len(rep_kos)
print(f"[osmo] representatives: {N} (rep_ko_v2, 0 unannotated)")

marker_n  = collections.Counter()
pathway_n = collections.Counter()
PATHWAYS  = list(flags(set()))
with open(f"{RES}/rep_traits.tsv", "w") as o:
    o.write("rep\tannotated\tphylum\tn_panel_kos\t" + "\t".join(PATHWAYS) + "\n")
    for rep, kos in sorted(rep_kos.items()):
        for ko in kos:
            marker_n[ko] += 1
        f = flags(kos)
        for p, v in f.items():
            pathway_n[p] += int(v)
        o.write(f"{rep}\t1\t{rep_phylum.get(rep,'')}\t{len(kos)}\t"
                + "\t".join(str(int(v)) for v in f.values()) + "\n")

with open(f"{RES}/marker_prevalence_5438.tsv", "w") as o:
    o.write("KO\tgene\tn_species\tpct_of_catalogue\n")
    for ko, n in marker_n.most_common():
        o.write(f"{ko}\t{KO_GENE[ko]}\t{n}\t{n/N*100:.1f}\n")

OLD = {"K00697": 67.8, "K01087": 64.3, "K10674": 17.9, "K06720": 17.0}   # 15,469-MAG basis
print(f"\n{'KO':8s} {'gene':38s} {'n':>6s} {'% 5,438':>8s} {'old % 15,469':>13s}")
for ko, n in marker_n.most_common():
    old = f"{OLD[ko]:.1f}" if ko in OLD else "-"
    print(f"{ko:8s} {KO_GENE[ko]:38s} {n:6d} {n/N*100:7.1f}% {old:>13s}")

with open(f"{RES}/pathway_prevalence_5438.tsv", "w") as o:
    o.write("pathway\tn_species\tpct_of_catalogue\n")
    print(f"\n{'pathway':26s} {'n':>6s} {'% of 5,438':>11s}")
    for p in PATHWAYS:
        o.write(f"{p}\t{pathway_n[p]}\t{pathway_n[p]/N*100:.1f}\n")
        print(f"{p:26s} {pathway_n[p]:6d} {pathway_n[p]/N*100:10.1f}%")
print(f"\n[osmo] wrote {RES}/{{marker,pathway}}_prevalence_5438.tsv + rep_traits.tsv")
