#!/usr/bin/env python3
"""
RQ-D3: Community atmospheric trace-gas (H2/CO) sink CAPACITY, and is it climate-relevant?

Run from repo root:  python3 RQ/D3/run.py
Inputs:  results/atmospheric/abundance/guild_abundance_by_sample.tsv
         results/atmospheric/abundance/guild_abundance_summary.tsv
         results/atmospheric/{guild_core_mags.txt, hydB_class_summary.tsv, coxL_form_summary.tsv}
         results/darkmatter_ecology/env_matrix.tsv  (sample -> site/compartment)
Outputs: RQ/D3/table_*.tsv, RQ/D3/fig_tracegas_capacity.{pdf,png}
Stdlib + matplotlib only. CAPACITY = genomic potential x abundance, NOT measured flux.
"""
import csv, os, statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ATMO = "results/atmospheric"
OUT = "RQ/D3"
os.makedirs(OUT, exist_ok=True)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def quant(xs, q):
    xs = sorted(xs)
    if not xs:
        return float("nan")
    i = q * (len(xs) - 1)
    lo, hi = int(i), min(int(i) + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (i - lo)


by_sample = read_tsv(os.path.join(ATMO, "abundance/guild_abundance_by_sample.tsv"))

# -------- Capacity index = per-sample, prevalence-weighted community guild fraction.
# Primary (conservative): full community, genus rank, metric "any" (= H2 OR CO high-affinity guild).
# Secondary (full-coverage): full community, phylum rank, metric "any".
def subset(table, rank, metric):
    return [r for r in by_sample
            if r["table"] == table and r["rank"] == rank and r["metric"] == metric]


prim = subset("full", "genus", "any")          # 1227 samples, conservative capacity index
phy = subset("full", "phylum", "any")          # full-coverage upper-ish
co = subset("full", "genus", "CO")
h2 = subset("full", "genus", "H2")
pres = subset("full", "genus", "any_presence")  # presence-based upper bound

prim_vals = [float(r["guild_fraction"]) for r in prim]
cov_vals = [float(r["covered_fraction"]) for r in prim]


def stats_row(label, vals):
    return [label, len(vals), f"{statistics.mean(vals):.4f}", f"{statistics.median(vals):.4f}",
            f"{quant(vals,0.25):.4f}", f"{quant(vals,0.75):.4f}",
            f"{min(vals):.4f}", f"{max(vals):.4f}", f"{statistics.pstdev(vals):.4f}"]


# =============================================== TABLE 1: distribution of the capacity index
with open(os.path.join(OUT, "table_capacity_distribution.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["index", "n_samples", "mean", "median", "q25", "q75", "min", "max", "sd"])
    w.writerow(stats_row("capacity_index (full/genus/any, prevalence-weighted)", prim_vals))
    w.writerow(stats_row("CO_only (full/genus/CO)", [float(r["guild_fraction"]) for r in co]))
    w.writerow(stats_row("H2_only (full/genus/H2)", [float(r["guild_fraction"]) for r in h2]))
    w.writerow(stats_row("presence_upper_bound (full/genus/any_presence)",
                         [float(r["guild_fraction"]) for r in pres]))
    w.writerow(stats_row("phylum_fullcoverage (full/phylum/any)",
                         [float(r["guild_fraction"]) for r in phy]))
    w.writerow(stats_row("covered_fraction (assessable share of community)", cov_vals))

# =============================================== map samples -> compartment / site
# by_sample already carries compartment. env_matrix gives site; sample id = "<run>_<envsample>".
env = {r["sample"]: r for r in read_tsv("results/darkmatter_ecology/env_matrix.tsv")}


def env_sample(sid):
    return sid.split("_", 1)[1] if "_" in sid else sid


# capacity by compartment (full coverage from by_sample)
comp_vals = {}
for r in prim:
    comp_vals.setdefault(r["compartment"], []).append(float(r["guild_fraction"]))

with open(os.path.join(OUT, "table_capacity_by_compartment.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["compartment", "n_samples", "mean", "median", "q25", "q75", "min", "max", "sd"])
    for comp in sorted(comp_vals):
        w.writerow(stats_row(comp, comp_vals[comp]))

# capacity by site (env_matrix mapping; subset that maps)
site_vals = {}
n_mapped = 0
for r in prim:
    es = env_sample(r["sample"])
    if es in env:
        n_mapped += 1
        site_vals.setdefault(env[es]["site"], []).append(float(r["guild_fraction"]))

with open(os.path.join(OUT, "table_capacity_by_site.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["site", "n_samples", "mean_capacity_index"])
    for s in sorted(site_vals, key=lambda x: (len(x), x)):
        w.writerow([s, len(site_vals[s]), f"{statistics.mean(site_vals[s]):.4f}"])

# =============================================== per-sample dump
with open(os.path.join(OUT, "table_capacity_per_sample.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample", "compartment", "site", "capacity_index", "covered_fraction"])
    for r in prim:
        es = env_sample(r["sample"])
        site = env[es]["site"] if es in env else ""
        w.writerow([r["sample"], r["compartment"], site, r["guild_fraction"], r["covered_fraction"]])

# =============================================== guild-genomics context table
hydb = read_tsv_simple = None
with open(os.path.join(ATMO, "hydB_class_summary.tsv")) as fh:
    hyd_rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
with open(os.path.join(ATMO, "coxL_form_summary.tsv")) as fh:
    cox_rows = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
n_core = sum(1 for _ in open(os.path.join(ATMO, "guild_core_mags.txt")) if _.strip())
g1l = next((r for r in hyd_rows if "1l" in r[0]), ["", "0", "0"])
g1h = next((r for r in hyd_rows if "1h" in r[0]), ["", "0", "0"])
formI = next((r for r in cox_rows if r[0] == "form_I"), ["", "0", "0"])

with open(os.path.join(OUT, "table_guild_genomics.tsv"), "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["marker", "gene_hits", "n_MAGs", "note"])
    w.writerow(["hydB [NiFe] Group 1l", g1l[1], g1l[2], "high-affinity H2; desert-dominant clade"])
    w.writerow(["hydB [NiFe] Group 1h", g1h[1], g1h[2], "high-affinity H2; classic Huc"])
    w.writerow(["coxL form-I", formI[1], formI[2], "high-affinity atmospheric CO oxidation"])
    w.writerow(["guild core MAGs (co-encode H2+CO)", "", n_core, "results/atmospheric/guild_core_mags.txt"])

# =============================================== FIGURE
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.6))
med = statistics.median(prim_vals)

# Panel A: distribution across all 1227 samples (histogram)
axA.hist([v * 100 for v in prim_vals], bins=30, color="#3a7d44", edgecolor="black", linewidth=0.4)
axA.axvline(med * 100, color="#b3001b", lw=2, ls="--",
            label=f"median {med*100:.1f}%")
axA.set_xlabel("Trace-gas CAPACITY index per sample\n(% community abundance in high-affinity H₂/CO guild)")
axA.set_ylabel("Samples")
axA.set_title(f"A. Distribution across {len(prim_vals)} samples", fontsize=10.5)
axA.legend(fontsize=9, frameon=False)

# Panel B: by compartment (boxplot)
order = [c for c in ["Surface", "Deep", "PlantRoot"] if c in comp_vals]
data = [[v * 100 for v in comp_vals[c]] for c in order]
bp = axB.boxplot(data, tick_labels=order, showfliers=False, patch_artist=True, widths=0.55)
for patch, c in zip(bp["boxes"], ["#8ecae6", "#bc6c25", "#52b788"]):
    patch.set_facecolor(c)
    patch.set_alpha(0.8)
for med_ln in bp["medians"]:
    med_ln.set_color("black")
for i, c in enumerate(order, 1):
    axB.text(i, statistics.median(comp_vals[c]) * 100 + 0.4,
             f"{statistics.median(comp_vals[c])*100:.1f}%", ha="center", fontsize=9)
axB.set_ylabel("Capacity index (% community abundance)")
axB.set_title("B. By compartment (uniform baseline)", fontsize=10.5)

fig.suptitle("Community atmospheric trace-gas (H₂/CO) sink CAPACITY "
             "(genomic potential × abundance, NOT measured flux)",
             fontsize=11.5, fontweight="bold")
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig_tracegas_capacity.pdf"))
fig.savefig(os.path.join(OUT, "fig_tracegas_capacity.png"), dpi=200)

# =============================================== console
print(f"[D3] capacity index (full/genus/any): n={len(prim_vals)} "
      f"median={med:.4f} ({med*100:.1f}%) range={min(prim_vals)*100:.1f}-{max(prim_vals)*100:.1f}% "
      f"IQR={quant(prim_vals,0.25)*100:.1f}-{quant(prim_vals,0.75)*100:.1f}%")
print(f"[D3] phylum full-coverage median={statistics.median([float(r['guild_fraction']) for r in phy])*100:.1f}% "
      f"presence-upper median={statistics.median([float(r['guild_fraction']) for r in pres])*100:.1f}%")
print(f"[D3] covered_fraction median={statistics.median(cov_vals)*100:.1f}% (assessable share)")
for comp in order:
    print(f"[D3]   {comp}: median {statistics.median(comp_vals[comp])*100:.1f}% (n={len(comp_vals[comp])})")
print(f"[D3] env_matrix mapped {n_mapped}/{len(prim)} samples to a site")
print(f"[D3] genomics MAGs: 1l={g1l[2]} 1h={g1h[2]} coxL-formI={formI[2]} guild-core={n_core}")
print("[D3] wrote tables + fig_tracegas_capacity.{pdf,png}")
