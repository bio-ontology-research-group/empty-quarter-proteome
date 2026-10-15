#!/usr/bin/env python3
"""Coarsened exact sensitivity matching for family-level descriptive dN/dS.

Exact strata: mean dS (0,.1,.25,.5,1,2), median CDS aa (30,60,100,200,400,inf),
median GC width .05, usable alignment sequence count (3,5,10,20,inf).
Select equal arm counts without replacement within each common stratum using a
fixed seed. The full family set defines the input; dN/dS is never a matching input.
Matching does not control phylogeny or ascertainment during family construction.
"""
from pathlib import Path
import gzip,hashlib,json
import numpy as np
import pandas as pd
from round2_checks import read
root=Path(__file__).resolve().parent;folder=root/'selection_covariates';source=folder/'family_covariates.tsv'
source = source if source.exists() else folder/'family_covariates.tsv.gz'
cov=pd.read_csv(source,sep='\t');d=pd.DataFrame(read('06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv'))
for k in ['n_seq','mean_dS','median_dNdS']:d[k]=pd.to_numeric(d[k])
d=d.merge(cov,on=['family_id','set'],validate='one_to_one');assert len(d)==192847
for label,col,edges in [('ds','mean_dS',[0,.1,.25,.5,1,2.000001]),('aa','median_cds_aa',[0,60,100,200,400,np.inf]),('seq','n_seq',[0,5,10,20,np.inf])]:d[label]=pd.cut(d[col],edges,right=False,labels=False)
d['gc']=np.floor(d['median_gc']/.05).astype(int)
groups=[]
for key,g in d.groupby(['ds','aa','seq','gc'],sort=True):
 a=g.index[g['set']=='dark'].to_numpy();b=g.index[g['set']=='ctrl'].to_numpy()
 if len(a) and len(b):groups.append((a,b))
def match(seed):
 rng=np.random.default_rng(seed);a=[];b=[]
 for x,y in groups:
  n=min(len(x),len(y));a.extend(rng.choice(x,n,replace=False));b.extend(rng.choice(y,n,replace=False))
 return d.loc[a],d.loc[b]
a,b=match(20261013);chosen=pd.concat([a,b]);chosen.to_csv(folder/'matched_families.tsv',sep='\t',index=False)
def summarise(a,b):
 out={'per_arm':len(a),'median_dark':float(a['median_dNdS'].median()),'median_control':float(b['median_dNdS'].median())}
 out['median_difference']=out['median_dark']-out['median_control'];out['balance']={}
 for k in ['mean_dS','median_cds_aa','median_gc','n_seq']:
  scale=np.sqrt((a[k].var()+b[k].var())/2)
  out['balance'][k]={'dark_median':float(a[k].median()),'control_median':float(b[k].median()),'standardised_mean_difference':float((a[k].mean()-b[k].mean())/scale)}
 return out
out={'seed':20261013,'common_strata':len(groups),'unmatched':summarise(d[d['set']=='dark'],d[d['set']=='ctrl']),'matched':summarise(a,b)}
out['unmatched']['arm_counts']={str(k):int(v) for k,v in d['set'].value_counts().items()};del out['unmatched']['per_arm']
deltas=[]
for seed in range(20261013,20261033):
 aa,bb=match(seed);deltas.append(float(aa['median_dNdS'].median()-bb['median_dNdS'].median()))
out['twenty_seed_difference_range']=[min(deltas),max(deltas)];out['covariates_sha256']=hashlib.sha256(gzip.decompress(source.read_bytes()) if source.suffix=='.gz' else source.read_bytes()).hexdigest()
(folder/'matching_summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
for name in ['family_covariates.tsv','matched_families.tsv']:
 p=folder/name
 if p.exists():(folder/(name+'.gz')).write_bytes(gzip.compress(p.read_bytes(),mtime=0))
