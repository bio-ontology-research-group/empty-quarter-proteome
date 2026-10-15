#!/usr/bin/env python3
"""Step 7 -- compositionality robustness of the Halobacteriota-salt coupling (step 6).

Relative abundances are compositional (fractions sum to 1), so a rising phylum can correlate
with an axis partly because others decline. Two checks, same partial-Spearman design as step 6:
  CLR       : centred log-ratio transform of phylum fractions (pseudocount = half min nonzero)
  detection : Halobacteriota present/absent per sample -- immune to closure entirely.
Out: results/salt_turnover_robustness.tsv
Verified console run 2026-08-13: CLR S +0.371 (1.7e-09), Cl +0.303, salinity +0.168, Ca +0.222;
detection S +0.441 (3.3e-13), Cl +0.369, salinity +0.299, Ca +0.244; detected in 36/248 samples.
"""
import csv, os, collections, math
import numpy as np
from scipy.stats import rankdata, pearsonr

M="/ibex/scratch/projects/c2014/rob/metagenomics"
REPO="/ibex/project/c2014/alelopezv/repos/empty-quarter-mags"
PROF=f"{M}/rebuild_20260729/picrust2_current_20260803/coverm_current/profiles_full_recompute"
RES=os.environ.get("RES_OUT", f"{M}/06_metabolism/osmolyte_clean/results")  # v551 rerun (2026-09-24)
os.makedirs(RES, exist_ok=True)
env={r["sample"]:r for r in csv.DictReader(open(f"{REPO}/results/darkmatter_ecology/env_matrix.tsv"),delimiter="\t")}
phylum={}
for p in (f"{M}/05_refine/_gtdbtk_r232/gtdbtk.bac120.summary.tsv", f"{M}/05_refine/_gtdbtk_r232/gtdbtk.ar53.summary.tsv"):
    for r in csv.DictReader(open(p),delimiter="\t"):
        phylum[r["user_genome"]]=next((t[3:] for t in r["classification"].split(";") if t.startswith("p__")),"")
candidates=[s for s in env if env[s].get("has_xrf")=="1" and os.path.isdir(f"{PROF}/{s}")]
allphy=sorted(set(phylum.values())-{""})
fr={}; det={}
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
    if tot<=0:  # 0 reads mapped (e.g. F5Dr1) -> exclude rather than inject a fake 0
        degenerate.append(s)
        continue
    samples.append(s)
    fr[s]={p:byp[p]/tot for p in allphy}
    det[s]=1 if byp["Halobacteriota"]>0 else 0
print(f"[robustness] degenerate profiles excluded (0 mapped reads): {degenerate}")
nz=[fr[s][p] for s in samples for p in allphy if fr[s][p]>0]; pc=min(nz)/2
clr={}
for s in samples:
    logs={p:math.log(fr[s][p]+pc) for p in allphy}
    gm=sum(logs.values())/len(logs)
    clr[s]=logs["Halobacteriota"]-gm
def partial(x,y,lat,lon):
    rx,ry=rankdata(x),rankdata(y)
    C=np.column_stack([np.ones(len(x)),rankdata(lat),rankdata(lon)])
    bx,*_=np.linalg.lstsq(C,rx,rcond=None); by,*_=np.linalg.lstsq(C,ry,rcond=None)
    return pearsonr(rx-C@bx,ry-C@by)
with open(f"{RES}/salt_turnover_robustness.tsv","w") as o:
    o.write("axis\tvariant\tn\trho_partial\tp_partial\n")
    for ax in ("salinity_proxy","S","Cl","Na","Ca"):
        sub=[s for s in samples if env[s].get(ax) not in ("",None) and env[s].get("lat")]
        y=np.array([float(env[s][ax]) for s in sub]); lat=np.array([float(env[s]["lat"]) for s in sub]); lon=np.array([float(env[s]["lon"]) for s in sub])
        for name,vec in (("clr",np.array([clr[s] for s in sub])),
                         ("detection",np.array([det[s] for s in sub],dtype=float))):
            r,pv=partial(vec,y,lat,lon)
            o.write(f"{ax}\t{name}\t{len(sub)}\t{r:.4f}\t{pv:.4g}\n")
print(f"detected in {sum(det.values())}/{len(samples)}; wrote {RES}/salt_turnover_robustness.tsv")
