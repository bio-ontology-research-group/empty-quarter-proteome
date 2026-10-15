#!/usr/bin/env python3
"""Rebuild fig_bgc_product + fig_amr_class on the CURRENT 5,438-representative basis
(§2.4 recompute, session 03-bgc-resistome, 2026-08-15).

Reads results/{bgc_product,amr_class,stress_class}_summary_5438.tsv (copied from the
workbench bgc_clean/ and resistome_clean/ result folders; the un-suffixed tables are
the frozen 15,469 historical basis and are not touched). Same style family as
code/make_figures.py, whose sections 1-2 are superseded by this script and pending
the figure-generator harmonisation session.

Run: python3 code/make_bgc_amr_figures.py
"""
import os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', '..', 'results')
FIG = os.path.join(HERE, '..', '..', 'figures')
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'savefig.bbox': 'tight'})
N = 5438

def read(name):
    with open(os.path.join(RES, name)) as f:
        r = csv.reader(f, delimiter='\t'); return next(r), list(r)

def save(fig, stem):
    for ext in ('pdf', 'png'): fig.savefig(os.path.join(FIG, f'{stem}.{ext}'))
    plt.close(fig); print(f'  wrote figures/{stem}.pdf/.png')

# 1. BGC product types by n_MAGs (top 12; terpene-precursor kept separate)
h, rows = read('bgc_product_summary_5438.tsv')
rows = sorted(rows, key=lambda r: -int(r[3]))[:12]
labels = [r[0] for r in rows]; vals = [int(r[3]) for r in rows]
pcts = [float(r[4]) for r in rows]
n = len(labels); fig, ax = plt.subplots(figsize=(7, max(2.5, 0.42*n+1)))
from palette import FOCAL, CONTROL
ax.barh(range(n), vals, color=FOCAL)
ax.set_yticks(range(n)); ax.set_yticklabels(labels); ax.invert_yaxis()
ax.set_xlabel(f'species with product type (of {N:,})')
ax.set_title('Biosynthetic gene clusters by product type', fontweight='bold', loc='left')
for i, (v, p) in enumerate(zip(vals, pcts)):
    ax.text(v, i, f' {v:,} ({p:.0f}%)', va='center', fontsize=9)
ax.margins(x=0.18); save(fig, 'fig_bgc_product')

# 2. Resistome: core-AMR classes (n_MAGs >= 30) + the metal/metalloid STRESS layer
h, amr = read('amr_class_summary_5438.tsv')
amr = [r for r in sorted(amr, key=lambda r: -int(r[2])) if int(r[2]) >= 30]
h, stress = read('stress_class_summary_5438.tsv')
stress = [r for r in sorted(stress, key=lambda r: -int(r[2]))
          if int(r[2]) >= 30 and r[0] != 'NA']
labels = [r[0].title() for r in amr] + [r[0].title() for r in stress]
vals   = [int(r[2]) for r in amr]    + [int(r[2]) for r in stress]
pcts   = [float(r[3]) for r in amr]  + [float(r[3]) for r in stress]
colors = [CONTROL] * len(amr) + [FOCAL] * len(stress)
n = len(labels); fig, ax = plt.subplots(figsize=(7, max(2.5, 0.42*n+1.4)))
ax.barh(range(n), vals, color=colors)
ax.set_yticks(range(n)); ax.set_yticklabels(labels); ax.invert_yaxis()
ax.set_xlabel(f'species with class (of {N:,})')
ax.set_title('Intrinsic resistome by class', fontweight='bold', loc='left')
for i, (v, p) in enumerate(zip(vals, pcts)):
    ax.text(v, i, f' {v:,} ({p:.0f}%)', va='center', fontsize=9)
ax.axhline(len(amr) - 0.5, color='#888888', lw=0.8, ls=':')
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=CONTROL, label='antibiotic (AMR)'),
                   plt.Rectangle((0, 0), 1, 1, color=FOCAL, label='metal/metalloid (STRESS)')],
          frameon=False, loc='lower right', fontsize=9)
ax.margins(x=0.18); save(fig, 'fig_amr_class')
