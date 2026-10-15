#!/usr/bin/env python3
"""Reviewer-driven sensitivity: permute whole visits of equal replicate count."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');O=W/'06_darkmatter/p4_review_20261013/results/round5_extensions';pins={}
def pin(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 pins[str(p)]=h.hexdigest();return p
meta=pd.read_csv(pin(W/'metadata/sample_metadata.tsv'),sep='\t',dtype=str).set_index('sample')
d=pd.read_csv(pin(W/'06_darkmatter/darkfam_basis_clean/results/fam_presence.tsv'),sep='\t',index_col=0);A=d.to_numpy(dtype=np.float32);m=meta.loc[d.columns].reset_index();trip=m.trip_num.astype(int).to_numpy()-1
visits=[np.asarray(ix) for ix in m.groupby(['site','compartment','trip_num']).indices.values()];v=pd.DataFrame([{'site':m.iloc[ix[0]].site,'compartment':m.iloc[ix[0]].compartment,'n':len(ix),'trip':trip[ix[0]]} for ix in visits]);blocks=[np.asarray(ix) for ix in v.groupby(['site','compartment','n']).indices.values()];vl=v.trip.to_numpy()
def statistic(labels):return float(((A@np.eye(5,dtype=np.float32)[labels])>0).all(axis=1).mean())
obs=statistic(trip);rng=np.random.default_rng(20261017);draw=[]
for _ in range(499):
 new=vl.copy()
 for ix in blocks:new[ix]=rng.permutation(vl[ix])
 labels=trip.copy()
 for i,ix in enumerate(visits):labels[ix]=new[i]
 assert np.array_equal(np.bincount(labels,minlength=5),np.bincount(trip,minlength=5))
 draw.append(statistic(labels))
x=np.asarray(draw);out={'design':'whole site-position-expedition visits permuted within site-position-replicate-count strata','seed':20261017,'permutations':499,'families':len(d),'samples':len(m),'visits':len(visits),'blocks':len(blocks),'exchangeable_visits':sum(len(ix) for ix in blocks if len(ix)>1),'exchangeable_samples':sum(int(v.iloc[ix].n.sum()) for ix in blocks if len(ix)>1),'observed':obs,'null_mean':float(x.mean()),'null_low':float(np.quantile(x,.025)),'null_high':float(np.quantile(x,.975)),'observed_minus_null':obs-float(x.mean()),'p_lower':float((1+(x<=obs).sum())/500),'p_upper':float((1+(x>=obs).sum())/500)}
(O/'temporal_visits_summary.json').write_text(json.dumps(out,indent=2)+'\n');pd.DataFrame({'all_five_frequency':x}).to_csv(O/'temporal_visits_draws.tsv',sep='\t',index=False);(O/'temporal_visits_input_sha256.json').write_text(json.dumps(pins,indent=2)+'\n');v.to_csv(O/'temporal_visits_design.tsv',sep='\t',index=False);print(json.dumps(out,indent=2))
