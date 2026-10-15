#!/usr/bin/env python3
"""Step 6 -- is "salt" answered by gene content or by community turnover?

Step 5 found osmolyte capacity DECLINES weakly with salinity -- opposite to the adaptation
hypothesis. Test the alternative: taxonomic turnover toward salt-in strategists. Per-sample
phylum fractions (of the mapped community) vs the salt axes, partial Spearman controlling
ranked lat/lon, XRF-gated samples. Halobacteriota (haloarchaea, K+ salt-in strategy, no
compatible-solute requirement) is the a-priori candidate.

Out: results/salt_turnover.tsv
"""
import csv, os, collections
import numpy as np
from scipy.stats import rankdata, pearsonr

M="/ibex/scratch/projects/c2014/rob/metagenomics"
REPO="/ibex/project/c2014/alelopezv/repos/empty-quarter-mags"
PROF=f"{M}/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"
RES=os.environ.get("RES_OUT", f"{M}/06_metabolism/osmolyte_clean/results")  # v551 rerun (2026-09-24)
os.makedirs(RES, exist_ok=True)
TOP=["Actinomycetota","Pseudomonadota","Bacteroidota","Chloroflexota","Halobacteriota",
     "Gemmatimonadota","Bacillota","Thermoplasmatota"]
AXES=["salinity_proxy","S","Cl","Na","Ca"]

env={r["sample"]:r for r in csv.DictReader(open(f"{REPO}/results/darkmatter_ecology/env_matrix.tsv"),delimiter="\t")}
phylum={}
for p in (f"{M}/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv",
          f"{M}/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"):
    for r in csv.DictReader(open(p),delimiter="\t"):
        phylum[r["user_genome"]]=next((t[3:] for t in r["classification"].split(";") if t.startswith("p__")),"")

candidates=[s for s in env if env[s].get("has_xrf")=="1" and os.path.isdir(f"{PROF}/{s}")]
frac=collections.defaultdict(dict)
samples=[]
degenerate=[]
for s in candidates:
    byp=collections.Counter(); tot=0.0
    with open(f"{PROF}/{s}/abund.tsv") as fh:
        next(fh)
        for line in fh:
            c=line.rstrip("\n").split("\t")
            if len(c)<2 or c[0]=="unmapped": continue
            try: v=float(c[1])
            except ValueError: continue
            if v>0: byp[phylum.get(c[0],"")]+=v; tot+=v
    if tot<=0:  # 0 reads mapped (e.g. F5Dr1) -> not a real zero, exclude rather than inject a fake 0%
        degenerate.append(s)
        continue
    samples.append(s)
    for p in TOP: frac[p][s]=byp[p]/tot*100
print(f"[turnover] degenerate profiles excluded (0 mapped reads): {degenerate}")

def partial(x,y,lat,lon):
    rx,ry=rankdata(x),rankdata(y)
    C=np.column_stack([np.ones(len(x)),rankdata(lat),rankdata(lon)])
    bx,*_=np.linalg.lstsq(C,rx,rcond=None); by,*_=np.linalg.lstsq(C,ry,rcond=None)
    return pearsonr(rx-C@bx,ry-C@by)

rows=[]
for p in TOP:
    for ax in AXES:
        sub=[s for s in samples if env[s].get(ax) not in ("",None) and env[s].get("lat")]
        x=np.array([frac[p][s] for s in sub]); y=np.array([float(env[s][ax]) for s in sub])
        lat=np.array([float(env[s]["lat"]) for s in sub]); lon=np.array([float(env[s]["lon"]) for s in sub])
        r,pv=partial(x,y,lat,lon)
        rows.append((p,ax,len(sub),r,pv))
with open(f"{RES}/salt_turnover.tsv","w") as o:
    o.write("phylum\taxis\tn\trho_partial\tp_partial\n")
    for p,ax,n,r,pv in rows: o.write(f"{p}\t{ax}\t{n}\t{r:.4f}\t{pv:.4g}\n")
print(f"{'phylum':18s} {'axis':>14s} {'rho_p':>8s} {'p':>10s}")
for p,ax,n,r,pv in sorted(rows,key=lambda t:t[4])[:12]:
    print(f"{p:18s} {ax:>14s} {r:+8.3f} {pv:10.2e}")
print(f"[turnover] wrote {RES}/salt_turnover.tsv  (n={len(samples)} XRF samples)")
