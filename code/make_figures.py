#!/usr/bin/env python3
"""Generate MAGs-paper figures from ../results/*.tsv into ../figures/. Pure matplotlib + csv
(no pandas). Run: python3 code/make_figures.py"""
import os, csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
FIG = os.path.join(HERE, '..', 'figures')
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'savefig.bbox': 'tight'})

def read(name):
    with open(os.path.join(RES, name)) as f:
        r = csv.reader(f, delimiter='\t'); return next(r), list(r)

def save(fig, stem):
    for ext in ('pdf', 'png'): fig.savefig(os.path.join(FIG, f'{stem}.{ext}'))
    plt.close(fig); print(f'  wrote figures/{stem}.pdf/.png')

def hbar(labels, values, title, xlabel, stem, color='#2c6fbb'):
    n = len(labels); fig, ax = plt.subplots(figsize=(7, max(2.5, 0.42*n+1)))
    ax.barh(range(n), values, color=color)
    ax.set_yticks(range(n)); ax.set_yticklabels(labels); ax.invert_yaxis()
    ax.set_xlabel(xlabel); ax.set_title(title, fontweight='bold', loc='left')
    for i, v in enumerate(values): ax.text(v, i, f' {v:,}', va='center', fontsize=9)
    ax.margins(x=0.13); save(fig, stem)

# 1. AMR class by n_MAGs
h, rows = read('amr_class_summary.tsv')
rows = sorted(rows, key=lambda r: -int(r[2]))[:12]
hbar([r[0].title() for r in rows], [int(r[2]) for r in rows],
     'Intrinsic resistome: AMR classes in EQ MAGs', 'MAGs with class (of 15,469)',
     'fig_amr_class', color='#b4503a')

# 2. BGC product types by n_MAGs (top 12)
h, rows = read('bgc_product_summary.tsv')
rows = sorted(rows, key=lambda r: -int(r[2]))[:12]
hbar([r[0] for r in rows], [int(r[2]) for r in rows],
     'Biosynthetic gene clusters by product type', 'MAGs with BGC type (of 15,469)',
     'fig_bgc_product', color='#7a59a5')

# 3. Life-on-air: coxL vs hydB by phylum (grouped) — top phyla by coxL
def by_phylum(fn):
    h, rows = read(fn); return {r[0]: int(r[1]) for r in rows}
cox = by_phylum('metabolism_coxL_by_phylum.tsv')
hyd = by_phylum('metabolism_hydB_by_phylum.tsv')
phyla = [p for p, _ in sorted(cox.items(), key=lambda x: -x[1])[:8]]
import numpy as np
y = np.arange(len(phyla)); bh = 0.4
fig, ax = plt.subplots(figsize=(7.5, 0.55*len(phyla)+1.2))
ax.barh(y-bh/2, [cox.get(p, 0) for p in phyla], bh, label='coxL (CO oxidation)', color='#2c6fbb')
ax.barh(y+bh/2, [hyd.get(p, 0) for p in phyla], bh, label='hydB (H$_2$ oxidation)', color='#5aa9e6')
ax.set_yticks(y); ax.set_yticklabels(phyla); ax.invert_yaxis()
ax.set_xlabel('MAGs encoding marker'); ax.legend(frameon=False)
ax.set_title('Atmospheric trace-gas metabolism ("life on air") by phylum', fontweight='bold', loc='left')
ax.margins(x=0.13); save(fig, 'fig_lifeonair_by_phylum')

print('mags figures done')
