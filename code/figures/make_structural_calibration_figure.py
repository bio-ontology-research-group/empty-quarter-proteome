#!/usr/bin/env python3
"""Section 1.3 figure — consensus design, revision 3 (2026-08-28).

Three-expert blind-panel design, revised twice against their render critiques and a Methods (M.6)
accuracy check. Key decisions on record: plane kept SQUARE (2-1: displacement from the 45-degree
aln-TM = qTM diagonal is the architectural-novelty signal; marginals are repositioned off the
plane's rendered box); panel b titled with the finding; control superposition stays FIRST
(author-approved slate; a right-end move was declined); no em-dashes anywhere in figure text
(author's style rule); violet = matched part of the known structure (panel c only).

Panels:
  a  From 1.5 million dark proteins to 1,031 candidates: constant-width chips (no width encoding),
     clustering rung with the proteins-to-families unit change marked (no false "1.7%"), controls'
     entry annotated with their composition, both gate criteria named, scored-zero footnote.
  b  Novelty is architectural, not local: square TM plane, plain-language axis labels, both
     marginals sqrt-scaled with below-gate tails shaded and annotated, zero-spikes explained.
  c  One covered control, three candidates: fixed A-per-pixel renders (scale.json) with a 20 A bar.
  d  Longer, equally confident: length + pLDDT densities, means marked, rest in neutral grey.

Every plotted number is asserted against summary_05c.json / per-model tables before drawing.
Palette: palette.py. Renders + superposition selection rule: render_structural_supers.py.
"""
import os, csv, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.join(HERE, "..", "..")
FIG = os.path.join(REPO, "figures")
REND = os.path.join(FIG, "renders_structural")
VAL = "/ibex/scratch/projects/c2014/rob/metagenomics/04_structure/00_structure_clean/05_validate"
T1, T2, T3 = 8.0, 7.5, 6.8   # approved 3-tier text grammar (axis/legend . labels . fine notes)
plt.rcParams.update({"font.size": T1, "axes.spines.top": False, "axes.spines.right": False})
from palette import FOCAL, CONTROL, REFKNOWN, GREY, INK, MUTED, GRID

GATE = 0.5

def load(fn):
    out = {}
    for r in csv.DictReader(open(os.path.join(VAL, fn)), delimiter="\t"):
        if r["confidence_qualified"] != "1":
            continue
        out[r["id"]] = (max(float(r["best_pdb_qtm"]), float(r["best_afdb_qtm"])),
                        max(float(r["best_pdb_alntm"]), float(r["best_afdb_alntm"])),
                        int(r["length"]), float(r["mean_plddt"]))
    return out

dark = load("per_rep_structures.tsv")
ctrl = load("per_control_structures.tsv")
S = json.load(open(os.path.join(VAL, "summary_05c.json")))
dx = np.array([v[0] for v in dark.values()]); dy = np.array([v[1] for v in dark.values()])
cx = np.array([v[0] for v in ctrl.values()]); cy = np.array([v[1] for v in ctrl.values()])
dk_q = int((dx < GATE).sum()); dk_a = int((dy < GATE).sum())
ct_q = int((cx < GATE).sum()); ct_a = int((cy < GATE).sum())
plume_d = int(((dx < GATE) & (dy >= GATE)).sum()); plume_c = int(((cx < GATE) & (cy >= GATE)).sum())
assert len(dark) == 2838 and len(ctrl) == 3013
assert dk_q == S["reps_conf_qualified_dual_novel_K"] == 1031 and dk_a == 95
assert ct_q == 101 and ct_a == 68 and plume_d == 936 and plume_c == 33
cand_len = np.array([v[2] for v in dark.values() if v[0] < GATE])
rest_len = np.array([v[2] for v in dark.values() if v[0] >= GATE])
cand_pl = np.array([v[3] for v in dark.values() if v[0] < GATE])
rest_pl = np.array([v[3] for v in dark.values() if v[0] >= GATE])
assert round(cand_len.mean()) == 301 and round(rest_len.mean()) == 215
print("asserted: 2,838/3,013; K 1,031 (plume 936, corner 95); controls 101/68/33; means 301/215")

SUPERS = [("super_control_KF000280", "KF000280", "known-fold control · 333 aa\ncovered end to end · qTM 0.97"),
          ("super_F1Dr1", "F1Dr1::F1Dr1_5662978_6", None),
          ("super_V45Dr1", "V45Dr1::V45Dr1_703031_2", None),
          ("super_22PRr3", "22PRr3::22PRr3_1915643_1", None)]
SUPERS = [(s, lab if lab else
           (lambda v: f"candidate · {v[2]} aa\naln-TM {v[1]:.2f} · qTM {v[0]:.2f}")(dark[pid]))
          for s, pid, lab in SUPERS]
VIEW_A = json.load(open(os.path.join(REND, "scale.json")))["viewport_A"]

fig = plt.figure(figsize=(9.9, 8.4))   # print-true canvas, widened so panel-a arrows breathe
gs = fig.add_gridspec(3, 1, height_ratios=[0.62, 1.55, 1.05],
                      left=0.045, right=0.985, top=0.935, bottom=0.045, hspace=0.34)

# ---- a: illustrative pipeline strip (letters only; explanations live in the caption) ------------------
a = fig.add_subplot(gs[0]); a.axis("off")
a.set_xlim(0.02, 8.92); a.set_ylim(0, 1)   # content extent: strip spans the full figure width
a.text(-0.02, 1.02, "a", transform=a.transAxes, fontsize=12, fontweight="bold", va="top")
def chip(x, y, w, n, col=FOCAL, hollow=False, fs=T2):
    a.add_patch(plt.Rectangle((x - w/2, y - 0.085), w, 0.17, facecolor="white" if hollow else col,
                              edgecolor=col, lw=1.2 if hollow else 0))
    a.text(x, y, n, fontsize=fs, color=col if hollow else "white", ha="center", va="center",
           fontweight="bold")
def arrow(x0, x1, y, lab, col=GREY):
    a.annotate("", (x1, y), (x0, y), arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0))
    a.text((x0 + x1) / 2, y - 0.10, lab, fontsize=T3, color=MUTED, ha="center", va="top")
def icon(fname, x, y, zoom):
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    img = plt.imread(os.path.join(REND, "icons", fname))
    a.add_artist(AnnotationBbox(OffsetImage(img, zoom=zoom), (x, y), frameon=False))

YD = 0.60
chip(0.62, YD, 1.05, "1,466,943")
a.text(0.62, YD + 0.20, "dark representatives", fontsize=T2, color=INK, ha="center")
arrow(1.25, 2.05, YD, "ProstT5 3Di + Foldseek\nvs AFDB v6 and PDB")
chip(2.62, YD, 0.95, "997,646", hollow=True)
a.text(2.62, YD + 0.20, "no match: screen only", fontsize=T2, color=MUTED, ha="center")
arrow(3.22, 4.02, YD, "Foldseek cluster,\nfamilies of 3 or more")
chip(4.62, YD, 0.95, "17,221")
a.text(4.62, YD + 0.20, "fold families", fontsize=T2, color=INK, ha="center")
arrow(5.22, 6.02, YD, "ESMFold;\nconfidence gate")
chip(6.55, YD, 0.85, "2,838")
a.text(6.55, YD + 0.20, "confident models", fontsize=T2, color=INK, ha="center")
arrow(7.10, 7.95, YD, "Foldseek vs AFDB + PDB:\nno covering match")
chip(8.42, YD, 0.85, "1,031")
a.text(8.42, YD + 0.20, "candidate architectures", fontsize=T2, color=INK, ha="center")
# control lane below, joining at the fold stage
YC = 0.16
chip(4.62, YC, 0.95, "3,684", col=CONTROL)
a.text(3.98, YC, "blinded known-fold controls", fontsize=T3, color=CONTROL, ha="right", va="center")
arrow(5.22, 6.02, YC, "same folding,\nsame gate", col=CONTROL)
chip(6.55, YC, 0.85, "3,013", col=CONTROL)
arrow(7.10, 7.95, YC, "same search; false novel:\nthe error floor", col=CONTROL)
chip(8.42, YC, 0.85, "101", col=CONTROL)

mid = gs[1].subgridspec(1, 2, width_ratios=[2.0, 1.15], wspace=0.12)

# ---- b: square hero plane, marginals hung off its rendered box ----------------------------------------
hb = mid[0, 0].subgridspec(2, 2, height_ratios=[0.15, 1.0], width_ratios=[1.0, 0.15],
                           hspace=0.02, wspace=0.02)
b = fig.add_subplot(hb[1, 0]); bt = fig.add_subplot(hb[0, 0]); br = fig.add_subplot(hb[1, 1])
b.fill_between([0, GATE], GATE, 1.02, color=FOCAL, alpha=0.05, zorder=0)
b.fill_between([0, GATE], 0, GATE, color=FOCAL, alpha=0.13, zorder=0)
b.scatter(cx, cy, s=5, color=CONTROL, alpha=0.28, lw=0, zorder=3, rasterized=True)
b.scatter(dx, dy, s=5, color=FOCAL, alpha=0.34, lw=0, zorder=4, rasterized=True)
b.plot([0, 1.02], [0, 1.02], ls=":", lw=0.8, color=GREY, zorder=1)
b.axvline(GATE, ls="--", lw=0.9, color=GREY, zorder=2); b.axhline(GATE, ls="--", lw=0.9, color=GREY, zorder=2)
b.text(0.03, 0.99, "locally familiar,\narchitecturally novel", fontsize=T2, color=INK, va="top")
b.annotate(f"{plume_d} dark models", (0.03, 0.99), xycoords="data",
           xytext=(0, -27), textcoords="offset points", fontsize=T2, color=MUTED, va="top")
b.annotate("novel under both", (0.44, 0.03), xycoords="data", xytext=(0, 12),
           textcoords="offset points", fontsize=T2, color=INK, va="bottom", ha="right")
b.text(0.44, 0.03, f"{dk_a} dark models", fontsize=T2, color=MUTED, va="bottom", ha="right")
b.annotate("covered end to end", (1.0, 0.03), xycoords="data", xytext=(0, 12),
           textcoords="offset points", fontsize=T2, color=INK, va="bottom", ha="right")
b.text(1.0, 0.03, f"{len(dark)-dk_q:,} dark models", fontsize=T2, color=MUTED,
       va="bottom", ha="right")
b.set_xlabel("whole-protein similarity (query-normalised TM)")
b.set_ylabel("local similarity (alignment-normalised TM)")
b.set_xlim(0, 1.02); b.set_ylim(0, 1.02)
b.set_box_aspect(1)

bins = np.linspace(0, 1.02, 42)
centers = (bins[:-1] + bins[1:]) / 2
for ax, vd, vc, horiz in ((bt, dx, cx, False), (br, dy, cy, True)):
    hd = np.sqrt(np.histogram(vd, bins=bins, density=True)[0])
    hc = np.sqrt(np.histogram(vc, bins=bins, density=True)[0])
    tail = centers < GATE
    if not horiz:
        ax.fill_between(centers, hd, color=FOCAL, alpha=0.55, lw=0)
        ax.fill_between(centers, hc, color=CONTROL, alpha=0.45, lw=0)
        ax.fill_between(centers[tail], hd[tail], color=FOCAL, alpha=0.95, lw=0)
        ax.axvline(GATE, ls="--", lw=0.8, color=GREY)
        ax.set_xlim(0, 1.02); ax.set_ylim(0, max(hd.max(), hc.max()) * 1.12)
        ax.annotate("dark 36.3%", (GATE - 0.02, 0.95), xycoords=ax.get_xaxis_transform(),
                    fontsize=T2, color=FOCAL, ha="right", va="top", fontweight="bold")
        ax.annotate("floor 3.35%", (GATE - 0.02, 0.95), xycoords=ax.get_xaxis_transform(),
                    xytext=(0, -11), textcoords="offset points",
                    fontsize=T2, color=CONTROL, ha="right", va="top", fontweight="bold")
        ax.text(1.015, 0.55, "density,\nsqrt scale", fontsize=T3, color=MUTED,
                transform=ax.transAxes, ha="left", va="center", clip_on=False)
    else:
        ax.fill_betweenx(centers, hd, color=FOCAL, alpha=0.55, lw=0)
        ax.fill_betweenx(centers, hc, color=CONTROL, alpha=0.45, lw=0)
        ax.fill_betweenx(centers[tail], hd[tail], color=FOCAL, alpha=0.95, lw=0)
        ax.axhline(GATE, ls="--", lw=0.8, color=GREY)
        ax.set_ylim(0, 1.02); ax.set_xlim(0, max(hd.max(), hc.max()) * 1.12)
        ax.annotate("dark 3.35%", (0.32, GATE - 0.03), xycoords=ax.get_yaxis_transform(),
                    fontsize=T2, color=FOCAL, ha="left", va="top", fontweight="bold",
                    annotation_clip=False)
        ax.annotate("floor 2.26%", (0.32, GATE - 0.03), xycoords=ax.get_yaxis_transform(),
                    xytext=(0, -11), textcoords="offset points", annotation_clip=False,
                    fontsize=T2, color=CONTROL, ha="left", va="top", fontweight="bold")
    ax.axis("off")
bt.set_title("b", loc="left", fontsize=12, fontweight="bold")
fig.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=FOCAL, label="dark fold-family models (2,838)"),
                    plt.Line2D([], [], marker="o", ls="", color=CONTROL, label="blinded known-fold controls (3,013)")],
           loc="upper right", bbox_to_anchor=(0.985, 1.0), frameon=False, fontsize=T1, ncol=2)

# ---- c: superpositions at fixed A-per-pixel (full-width bottom row) -----------------------------------
bot = gs[2].subgridspec(1, 4, wspace=0.05)
imgs = [plt.imread(os.path.join(REND, stem + ".png")) for stem, _ in SUPERS]
margins = []
for img in imgs:
    alpha = img[:, :, 3] if img.shape[2] == 4 else 1 - img[:, :, :3].min(axis=2)
    ys, xs = np.where(alpha > 0.05)
    margins.append(min(ys.min(), xs.min(), img.shape[0] - ys.max(), img.shape[1] - xs.max()))
CROP = max(0, min(margins) - 12)   # same pixels off every side of every image: shared A/px kept
for k, ((stem, lab), img) in enumerate(zip(SUPERS, imgs)):
    ax = fig.add_subplot(bot[0, k])
    img = img[CROP:img.shape[0] - CROP, CROP:img.shape[1] - CROP]
    ax.imshow(img)
    ax.set_xticks([]); ax.set_yticks([])
    anchor = (k == 0)
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color(GRID); sp.set_linewidth(0.7)
    ax.set_xlabel(lab, fontsize=T1, color=INK, labelpad=3)
    if k == 0:
        w = img.shape[1]; h = img.shape[0]
        bar = 20.0 / VIEW_A * (w + 2 * CROP)
        ax.plot([0.05 * w, 0.05 * w + bar], [0.955 * h] * 2, color=INK, lw=1.6)
        ax.text(0.05 * w + bar / 2, 0.935 * h, "20 Å", fontsize=T3, color=INK, ha="center",
                va="bottom")
        ax.text(-0.04, 1.06, "c", transform=ax.transAxes, fontsize=12, fontweight="bold",
                va="top", ha="left")

for _ax in fig.axes:
    if _ax.get_label() == "":
        pass
for k in range(4):
    pass
# square tiles matched to row height
dgs = mid[0, 1].subgridspec(2, 1, height_ratios=[1.0, 0.62], hspace=0.5)
d1 = fig.add_subplot(dgs[0]); d2 = fig.add_subplot(dgs[1])
lb = np.linspace(90, 900, 55)
for vals, col, lab in ((rest_len, GREY, "other confident models"), (cand_len, FOCAL, "candidates")):
    h = np.histogram(vals, bins=lb, density=True)[0]
    d1.fill_between((lb[:-1] + lb[1:]) / 2, h, color=col, alpha=0.55, lw=0, label=lab)
for m, col, ha, off in ((cand_len.mean(), FOCAL, "left", 8), (rest_len.mean(), GREY, "right", -8)):
    d1.axvline(m, color=col, lw=1.1)
    d1.text(m + off, 0.90, f"{m:.0f} aa", transform=d1.get_xaxis_transform(), fontsize=7,
            color=INK, ha=ha)
d1.set_yticks([]); d1.set_xlim(90, 900); d1.tick_params(labelsize=8)
d1.set_xlabel("model length (aa)", fontsize=T1); d1.set_ylabel("share of models", fontsize=T1)
d1.legend(loc="upper right", frameon=False, fontsize=T2)
d1.set_title("d", loc="left", fontsize=12, fontweight="bold")
pb = np.linspace(0.70, 1.0, 32)
for vals, col in ((rest_pl, GREY), (cand_pl, FOCAL)):
    h = np.histogram(vals, bins=pb, density=True)[0]
    d2.fill_between((pb[:-1] + pb[1:]) / 2, h, color=col, alpha=0.55, lw=0)
for m, col, ha, off in ((cand_pl.mean(), FOCAL, "right", -0.004), (rest_pl.mean(), GREY, "left", 0.004)):
    d2.axvline(m, color=col, lw=1.1)
    d2.text(m + off, 0.86, f"{m:.2f}", transform=d2.get_xaxis_transform(), fontsize=T2,
            color=INK, ha=ha)
d2.set_yticks([]); d2.set_xlim(0.70, 1.0); d2.tick_params(labelsize=8)
d2.set_xlabel("model confidence (pLDDT): indistinguishable", fontsize=T1)
d2.set_ylabel("share of models", fontsize=T1)

# hang the marginals off the plane's rendered (square) box
fig.canvas.draw()
pos = b.get_position()
g = 0.006
pt = bt.get_position(); pr = br.get_position()
bt.set_position([pos.x0, pos.y1 + g, pos.width, pt.height])
br.set_position([pos.x1 + g, pos.y0, pr.width, pos.height])

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"fig_structural_calibration.{ext}"), bbox_inches="tight", dpi=200)
print("wrote figures/fig_structural_calibration.{pdf,png} (revision 3)")
