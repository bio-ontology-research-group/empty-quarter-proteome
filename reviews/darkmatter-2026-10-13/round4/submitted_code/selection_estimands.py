#!/usr/bin/env python3
"""Stratum-specific descriptive dN/dS, bootstrap intervals and target weighting.

Intervals resample families independently within each arm/stratum; they are
conditional descriptive intervals and do not correct phylogenetic dependence.
"""
import sys,json
from pathlib import Path
import numpy as np,pandas as pd
if len(sys.argv)==1:
 from round2_checks import read
 here=Path(__file__).resolve().parent;out=here/'selection_estimands';out.mkdir(exist_ok=True)
 d=pd.DataFrame(read('06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv'))
 for col in ['mean_dS','median_dNdS','n_seq']:d[col]=pd.to_numeric(d[col])
 c=pd.read_csv(here/'selection_covariates/family_covariates.tsv.gz',sep='\t')
else:
 base=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
 d=pd.read_csv(base/'06_darkmatter/darkfam_basis_clean/results/selection/dnds.tsv',sep='\t');c=pd.read_csv(base/'06_darkmatter/p4_review_20261013/results/selection_covariates/family_covariates.tsv',sep='\t')
d=d.merge(c,on=['family_id','set'],validate='one_to_one')
d['ds']=pd.cut(d.mean_dS,[0,.1,.25,.5,1,2.000001],right=False,labels=False);d['aa']=pd.cut(d.median_cds_aa,[0,60,100,200,400,np.inf],right=False,labels=False);d['gc']=np.floor(d.median_gc/.05);d['seq']=pd.cut(d.n_seq,[0,5,10,20,np.inf],right=False,labels=False)
rng=np.random.default_rng(20261014);records=[]
for k,g in d.groupby('ds',sort=True):
 a=g[g['set']=='dark'].median_dNdS.to_numpy();b=g[g['set']=='ctrl'].median_dNdS.to_numpy();boots=[]
 for _ in range(1000):boots.append(np.median(rng.choice(a,len(a),replace=True))-np.median(rng.choice(b,len(b),replace=True)))
 records.append(dict(stratum=int(k),dark_n=len(a),control_n=len(b),dark_median=float(np.median(a)),control_median=float(np.median(b)),difference=float(np.median(a)-np.median(b)),bootstrap_low=float(np.percentile(boots,2.5)),bootstrap_high=float(np.percentile(boots,97.5))))
strata=[]
for key,g in d.groupby(['ds','aa','gc','seq'],sort=True):
 a=g[g['set']=='dark'].median_dNdS;b=g[g['set']=='ctrl'].median_dNdS
 if len(a) and len(b):strata.append([len(a),len(b),a.mean()-b.mean()])
s=np.array(strata);weighted={label:float(np.average(s[:,2],weights=w)) for label,w in [('dark_distribution',s[:,0]),('control_distribution',s[:,1]),('overlap_distribution',np.minimum(s[:,0],s[:,1]))]}
result={'bootstrap_seed':20261014,'bootstrap_replicates':1000,'bootstrap_scope':'families independent within each arm/stratum; shared ancestry unmodelled','strata':records,'weighted_mean_differences_of_family_medians':weighted,'common_strata':len(s),'common_support_dark':int(s[:,0].sum()),'common_support_control':int(s[:,1].sum())}
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');pd.DataFrame(records).to_csv(out/'strata.tsv',sep='\t',index=False);print(json.dumps(result,indent=2))
