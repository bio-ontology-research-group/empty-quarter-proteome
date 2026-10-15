#!/usr/bin/env python3
r"""RQ-D6 driver: rhodopsin light-energy harvesting & light niche (surface vs deep).

This script (a) documents the IBEX census that produced rhodopsin_hits.tsv, and
(b) processes the pulled summary into rhodopsin_census.tsv + fig_rhodopsin.{pdf,png}.
Run from this directory:  python3 run.py

--------------------------------------------------------------------------------------------
STEP 1 -- locate annotations (IBEX, login node, read-only):
  M=/ibex/scratch/projects/c2014/rob/metagenomics
  # per-MAG functional annotations are EggNOG emapper outputs, one per MAG; paths are derived
  # from $M/06_mobilome/mag_index.tsv (col1=mag_id, col2=bakta .faa). The sibling EggNOG file is:
  #   faa:  .../single_sample/bakta/<bin>/<bin>.faa
  #   egg:  .../single_sample/eggnog/<bin>/<bin>.emapper.annotations   (s/\/bakta\//\/eggnog\//; .faa->.emapper.annotations)
  # emapper columns: 8=Description 9=Preferred_name 12=KEGG_ko 21=PFAMs
  # all 15,469 MAGs have an EggNOG file (0 missing).

STEP 2 -- census (Slurm batch job, NOT login node):
  scripts/ibex/D6_rho_census.sh  ->  sbatch --account=pi-hohndor --partition=batch  (job 47843530)
  Marks a protein as a microbial type-1 rhodopsin if its PFAMs carry Bac_rhodopsin (PF01036) or
  Heliorhodopsin, OR KEGG_ko in K04641-K04644, OR description/name contains "rhodopsin" --
  EXCLUDING ASRT / "...rhodopsin transducer" soluble signalling proteins (not opsins).
  Emits: mag <tab> protein <tab> ko <tab> pname <tab> pfams <tab> description.
  Result on IBEX: 398 candidate proteins -> 341 MAGs.

STEP 3 -- pull summary (no heavy compute pulled, just the small hit table):
  scp/ssh cat $M/06_metabolism/rhodopsin/rhodopsin_hits.tsv  -> RQ/D6/rhodopsin_hits.tsv
  ssh cut -f1 $M/06_mobilome/mag_index.tsv                   -> RQ/D6/all_mags.txt (15,469 ids)

STEP 4 -- local QC + analysis (below):
  Drop 7 eukaryotic GPCR "7-transmembrane receptor (rhodopsin family)" hits (Pfam 7tm_1 /
  Orexin_rec2 / adh_short; eukaryotic contamination in 4 PlantRoot MAGs) -> keep 391 genuine
  type-1 opsins in 337 MAGs. Sub-type from informative KO/name only. Join to compartment via
  MAG->sample (id prefix before "__") -> env_matrix.tsv; Fisher Surface vs Deep.
--------------------------------------------------------------------------------------------
"""
import os, csv
from collections import defaultdict, Counter
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))            # RQ/D6
ROOT = os.path.join(HERE, "..", "..")                        # repo root
ECO = os.path.join(ROOT, "results", "darkmatter_ecology", "env_matrix.tsv")
TOTAL_MAGS = 15469
HITS = os.path.join(HERE, "rhodopsin_hits.tsv")
ALL_MAGS = os.path.join(HERE, "all_mags.txt")

# ---- load census hits; keep only genuine type-1 microbial rhodopsins (Bac_rhodopsin / Helio) ----
mag_hits = defaultdict(list)
n_dropped = 0
for ln in open(HITS):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 6:
        continue
    mag, prot, ko, pname, pfam, desc = c[:6]
    pf = pfam.lower()
    if "bac_rhodopsin" not in pf and "heliorhodopsin" not in pf:
        n_dropped += 1                                       # eukaryotic GPCR rhodopsin-family
        continue
    mag_hits[mag].append((ko, pname, pfam, desc))
print(f"dropped {n_dropped} non-opsin (eukaryotic GPCR rhodopsin-family) proteins")

rho_mags = set(mag_hits)
n_rho = len(rho_mags)
n_prot = sum(len(v) for v in mag_hits.values())
print(f"distinct rhodopsin MAGs: {n_rho}  ({100*n_rho/TOTAL_MAGS:.2f}% of {TOTAL_MAGS})")
print(f"total rhodopsin proteins: {n_prot}")


def subtype(ko, pname, pfam, desc):
    # subtype ONLY from informative KO / explicit name; generic "Bacteriorhodopsin-like protein"
    # is the type-1 superfamily descriptor and does not resolve subtype.
    t = f"{pname} {desc}".lower(); pf = pfam.lower(); kos = ko.replace("ko:", "")
    if "heliorhodopsin" in pf or "helio" in t: return "heliorhodopsin"
    if "xantho" in t: return "xanthorhodopsin"
    if "proteo" in t: return "proteorhodopsin"
    if "K04641" in kos: return "bacteriorhodopsin-type (K04641)"
    if "K04642" in kos or "K04643" in kos: return "sensory-rhodopsin-type (K04642/3)"
    if "K04644" in kos: return "halorhodopsin-type (K04644)"
    return "type-1 rhodopsin (subtype unresolved)"


prot_sub = Counter(); mag_subs = defaultdict(set)
for mag, hits in mag_hits.items():
    for h in hits:
        s = subtype(*h); prot_sub[s] += 1; mag_subs[s].add(mag)
print("\nsub-type (proteins / distinct MAGs):")
for s, n in prot_sub.most_common():
    print(f"  {s}\t{n}\t{len(mag_subs[s])}")

# ---- compartment join (denominator = ALL MAGs in compartment) ----
samp_comp = {r["sample"]: r["compartment"] for r in csv.DictReader(open(ECO), delimiter="\t")}
comp_total = Counter(); comp_rho = Counter()
for mag in (l.strip() for l in open(ALL_MAGS) if l.strip()):
    comp = samp_comp.get(mag.split("__")[0])
    if comp not in ("Surface", "Deep", "PlantRoot"):
        continue
    comp_total[comp] += 1
    if mag in rho_mags:
        comp_rho[comp] += 1
print("\ncompartment  rho/total  pct")
for c in ("Surface", "Deep", "PlantRoot"):
    if comp_total[c]:
        print(f"  {c}\t{comp_rho[c]}/{comp_total[c]}\t{100*comp_rho[c]/comp_total[c]:.2f}%")

sS, tS = comp_rho["Surface"], comp_total["Surface"]
sD, tD = comp_rho["Deep"], comp_total["Deep"]
OR, p = fisher_exact([[sS, tS - sS], [sD, tD - sD]])
pctS, pctD = 100*sS/tS, 100*sD/tD
print(f"\nSurface vs Deep: {pctS:.2f}% vs {pctD:.2f}%  OR={OR:.2f}  Fisher p={p:.3g}")

# ---- write census tsv ----
with open(os.path.join(HERE, "rhodopsin_census.tsv"), "w") as o:
    o.write("# RQ-D6 rhodopsin census (EggNOG of 15,469 MAGs; IBEX job 47843530)\n")
    o.write("metric\tvalue\n")
    o.write(f"total_MAGs\t{TOTAL_MAGS}\n")
    o.write(f"rhodopsin_MAGs\t{n_rho}\n")
    o.write(f"rhodopsin_MAGs_pct\t{100*n_rho/TOTAL_MAGS:.2f}\n")
    o.write(f"rhodopsin_proteins\t{n_prot}\n")
    o.write("#\n# sub-type\tproteins\tdistinct_MAGs\n")
    for s, n in prot_sub.most_common():
        o.write(f"subtype:{s}\t{n}\t{len(mag_subs[s])}\n")
    o.write("#\n# compartment\trho_MAGs\ttotal_MAGs\tpct\n")
    for c in ("Surface", "Deep", "PlantRoot"):
        if comp_total[c]:
            o.write(f"compartment:{c}\t{comp_rho[c]}\t{comp_total[c]}\t{100*comp_rho[c]/comp_total[c]:.2f}\n")
    o.write(f"#\nfisher_surface_vs_deep_OR\t{OR:.3f}\n")
    o.write(f"fisher_surface_vs_deep_p\t{p:.3g}\n")
    o.write(f"pct_surface\t{pctS:.2f}\n")
    o.write(f"pct_deep\t{pctD:.2f}\n")

# ---- figure ----
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4.0))
ax1.bar(["rhodopsin+", "no rhodopsin"], [n_rho, TOTAL_MAGS - n_rho], color=["#d98a00", "#cccccc"])
ax1.set_ylabel("MAGs")
ax1.set_title(f"Rhodopsin prevalence\n{n_rho}/{TOTAL_MAGS} = {100*n_rho/TOTAL_MAGS:.1f}% of MAGs")
for i, v in enumerate([n_rho, TOTAL_MAGS - n_rho]):
    ax1.text(i, v, f"{v}", ha="center", va="bottom", fontsize=9)
comps = ["Surface", "Deep", "PlantRoot"]
rho_pct = [100*comp_rho[c]/comp_total[c] if comp_total[c] else 0 for c in comps]
photo_pct = {"Surface": 22.4, "Deep": 13.1, "PlantRoot": 17.8}   # results/uv_survival/compartment.tsv
x = range(len(comps)); w = 0.38
ax2.bar([i - w/2 for i in x], rho_pct, w, label="rhodopsin", color="#d98a00")
ax2.bar([i + w/2 for i in x], [photo_pct[c] for c in comps], w, label="photolyase (ref)", color="#4477aa")
ax2.set_xticks(list(x)); ax2.set_xticklabels(comps)
ax2.set_ylabel("% of MAGs in compartment")
ax2.set_title(f"Light-niche split\nrho Surface {pctS:.1f}% vs Deep {pctD:.1f}% (Fisher p={p:.2g})")
ax2.legend(fontsize=8)
for i, v in enumerate(rho_pct):
    ax2.text(i - w/2, v, f"{v:.1f}", ha="center", va="bottom", fontsize=7)
plt.tight_layout()
plt.savefig(os.path.join(HERE, "fig_rhodopsin.pdf"))
plt.savefig(os.path.join(HERE, "fig_rhodopsin.png"), dpi=150)
print("\nwrote rhodopsin_census.tsv, fig_rhodopsin.{pdf,png}")
