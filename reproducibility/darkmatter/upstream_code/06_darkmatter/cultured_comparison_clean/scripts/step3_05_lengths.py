#!/usr/bin/env python3
"""step3_05 - length distribution of cultured dark proteins vs MAG dark proteins.

Cultured genomes are long-read assemblies, about half of them a single contig, so a cultured dark
protein is essentially never a fragmentation artefact. If the MAG dark set were mostly assembly
noise it would sit far to the left of the cultured dark set; the figure and the table below
measure that gap.

  results/step3_length_summary.tsv
  figures/step3_length_cultured_vs_mag_dark.png
"""
import gzip
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

P = "/ibex/scratch/projects/c2014/rob/metagenomics/06_darkmatter/cultured_comparison_clean"
W = "/ibex/scratch/projects/c2014/rob/metagenomics"
DARK_FAA = f"{W}/06_darkmatter/darkfam_basis_clean/results/dark_reps.faa"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#dcdbd6"
S_MAG = "#2a78d6"      # categorical slot 1
S_CULT = "#eb6834"     # categorical slot 2

# ---- MAG dark rep id, codon count and partial flag, straight off the prodigal headers of
# dark_reps.faa:  ">id # start # end # strand # ID=..;partial=10;.."  ->  id, (end-start+1)//3, 10
# Written once to tmp/ and reused; the exploratory folder had this table but no script for it.
LEN_TSV = f"{P}/tmp/mag_dark_len.tsv.gz"
if not os.path.exists(LEN_TSV):
    os.makedirs(f"{P}/tmp", exist_ok=True)
    n = 0
    with open(DARK_FAA) as fh, gzip.open(LEN_TSV + ".tmp", "wt") as out:
        for line in fh:
            if line[0] != ">":
                continue
            f = line[1:].rstrip("\n").split(" # ")
            aa = (int(f[2]) - int(f[1]) + 1) // 3
            partial = f[4].split("partial=")[1].split(";")[0]
            out.write(f"{f[0]}\t{aa}\t{partial}\n")
            n += 1
    os.rename(LEN_TSV + ".tmp", LEN_TSV)
    print(f"[info] built {os.path.basename(LEN_TSV)}: {n} MAG dark reps")

# ---- MAG dark rep lengths (aa from the prodigal coordinates) + partial flag
mag_len = []
mag_partial = 0
with gzip.open(LEN_TSV, "rt") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        aa = int(f[1])
        if aa < 0:
            continue
        mag_len.append(aa)
        if f[2] != "00":
            mag_partial += 1
mag_len = np.array(mag_len, dtype=np.int32)

# ---- cultured dark protein lengths
with open(f"{P}/results/step3_cultured_counts.tsv") as fh:
    n_cultured_genomes = sum(1 for _ in fh) - 1
dark_ids = set()
with open(f"{P}/tmp/cultured_dark_ids.txt") as fh:
    for line in fh:
        dark_ids.add(line.strip())
cult_len = []
cur = None
n = 0
with open(f"{P}/tmp/cultured_all.faa") as fh:
    for line in fh:
        if line[0] == ">":
            if cur is not None and n:
                cult_len.append(n)
            cur = line[1:].split()[0] in dark_ids
            n = 0
        elif cur:
            n += len(line.strip())
    if cur is not None and n:
        cult_len.append(n)
cult_len = np.array(cult_len, dtype=np.int32)


def stats(a):
    return dict(n=int(a.size), median=float(np.median(a)),
                q1=float(np.percentile(a, 25)), q3=float(np.percentile(a, 75)),
                lt100=float((a < 100).mean()), lt60=float((a < 60).mean()))


sm, sc = stats(mag_len), stats(cult_len)
with open(f"{P}/results/step3_length_summary.tsv", "w") as out:
    out.write("set\tn\tmedian_aa\tq1_aa\tq3_aa\tshare_under_100aa\tshare_under_60aa\t"
              "share_incomplete_gene_call\n")
    out.write(f"mag_dark_reps\t{sm['n']}\t{sm['median']:.0f}\t{sm['q1']:.0f}\t{sm['q3']:.0f}\t"
              f"{sm['lt100']:.4f}\t{sm['lt60']:.4f}\t{mag_partial / max(sm['n'], 1):.4f}\n")
    out.write(f"cultured_dark_proteins\t{sc['n']}\t{sc['median']:.0f}\t{sc['q1']:.0f}\t"
              f"{sc['q3']:.0f}\t{sc['lt100']:.4f}\t{sc['lt60']:.4f}\tNA\n")

# ------------------------------------------------------------------ figure
fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), gridspec_kw={"width_ratios": [2.3, 1]})
fig.patch.set_facecolor(SURFACE)


def ecdf(ax, a, color, label):
    x = np.sort(a)
    y = np.arange(1, x.size + 1) / x.size
    k = max(1, x.size // 20000)
    ax.plot(x[::k], y[::k], color=color, lw=2, label=label, solid_capstyle="round")


ax = axes[0]
ax.set_facecolor(SURFACE)
ecdf(ax, mag_len, S_MAG, "MAG dark reps")
ecdf(ax, cult_len, S_CULT, "Cultured dark proteins")
ax.axvline(100, color=INK2, lw=1, ls=(0, (4, 3)))
ax.text(103, 0.03, "100 aa", color=INK2, fontsize=9)
ax.set_xscale("log")
ax.set_xlim(20, 3000)
ax.set_ylim(0, 1)
ax.set_xlabel("Protein length (aa, log scale)", color=INK2, fontsize=10)
ax.set_ylabel("Cumulative share of proteins", color=INK2, fontsize=10)
ax.set_title("Length distribution of dark proteins", color=INK, fontsize=12, loc="left")
ax.grid(True, color=GRID, lw=0.6)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)
leg = ax.legend(frameon=False, loc="lower right", fontsize=9)
for t in leg.get_texts():
    t.set_color(INK2)
ax.annotate(f"median {sm['median']:.0f} aa", xy=(sm["median"], 0.5),
            xytext=(sm["median"] * 0.30, 0.60), color=INK, fontsize=9,
            arrowprops=dict(arrowstyle="-", color=S_MAG, lw=1.2))
ax.annotate(f"median {sc['median']:.0f} aa", xy=(sc["median"], 0.5),
            xytext=(125, 0.30), color=INK, fontsize=9,
            arrowprops=dict(arrowstyle="-", color=S_CULT, lw=1.2))

ax = axes[1]
ax.set_facecolor(SURFACE)
vals = [sm["lt100"] * 100, sc["lt100"] * 100]
bars = ax.bar([0, 1], vals, width=0.55, color=[S_MAG, S_CULT])
for b in bars:
    b.set_linewidth(2)
    b.set_edgecolor(SURFACE)
for x, v in zip([0, 1], vals):
    ax.text(x, v + 1.5, f"{v:.1f}%", ha="center", color=INK, fontsize=10)
ax.set_xticks([0, 1])
ax.set_xticklabels(["MAG dark\nreps", "Cultured dark\nproteins"], color=INK2, fontsize=9)
ax.set_ylim(0, max(vals) * 1.28)
ax.set_ylabel("Share under 100 aa (%)", color=INK2, fontsize=10)
ax.set_title("Short proteins", color=INK, fontsize=12, loc="left")
ax.grid(True, axis="y", color=GRID, lw=0.6)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=INK2, labelsize=9)

fig.text(0.01, 0.015,
         f"MAG dark reps n={sm['n']:,} ({mag_partial / max(sm['n'], 1) * 100:.1f}% incomplete gene "
         f"calls).  Cultured dark proteins n={sc['n']:,} from {n_cultured_genomes} long-read "
         f"genomes, Bakta gene calls, fragmentation negligible.",
         color=INK2, fontsize=8)
fig.tight_layout(rect=(0, 0.045, 1, 1))
os.makedirs(f"{P}/figures", exist_ok=True)
fig.savefig(f"{P}/figures/step3_length_cultured_vs_mag_dark.png", dpi=200,
            facecolor=SURFACE)
print(f"[info] MAG dark {sm}")
print(f"[info] cultured dark {sc}")
print("DONE_STEP3_05")
