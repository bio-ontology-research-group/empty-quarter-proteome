#!/usr/bin/env python3
"""Exploratory site-aggregated family/environment associations (IBEX/Slurm).

Average abundances and chemistry within site using exactly the samples with each
chemical measurement. Rank-correlate across sites; BH per chemical across families
present in >=25 measured samples and >=5 sites. This removes replicate pseudoreplication, but does not resolve
chemistry measurement provenance, compositionality, or site sampling imbalance.
"""
import sys,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,t
W=Path('/ibex/scratch/projects/c2014/rob/metagenomics');E=Path('/ibex/project/c2014/alelopezv/repos/empty-quarter-mags/results/darkmatter_ecology/env_matrix.tsv');A=W/'06_darkmatter/darkfam_basis_clean/results/fam_abundance.tsv';out=Path(sys.argv[1]);out.mkdir(parents=True,exist_ok=True)
env=pd.read_csv(E,sep='\t').set_index('sample');ab=pd.read_csv(A,sep='\t',index_col=0)
def bh(p):
 order=np.argsort(p);q=np.empty(len(p));q[order]=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(q,1)
summary=[]
for axis in ['S','Ca','Cl','Fe']:
 assert axis in env.columns,env.columns
 e=env.reindex(ab.columns);y=pd.to_numeric(e[axis],errors='coerce');valid=y.notna() & e['site'].notna() & e['has_xrf'].astype(str).str.lower().isin(['true','1','1.0','yes'])
 ev=e[valid];v=y[valid];columns=ab.columns[valid];sites=sorted(ev['site'].unique());x=[];chem=[]
 for site in sites:
  ids=ev.index[ev['site']==site];x.append(ab[ids].mean(axis=1).to_numpy());chem.append(v[ids].median())
 X=np.asarray(x).T;Y=np.asarray(chem);keep=((X>0).sum(axis=1)>=5)&((ab[columns]>0).sum(axis=1).to_numpy()>=25);X=X[keep];ids=ab.index[keep]
 rx=rankdata(X,axis=1);ry=rankdata(Y);rx-=rx.mean(axis=1,keepdims=True);ry-=ry.mean();den=np.sqrt((rx*rx).sum(axis=1)*(ry*ry).sum());corr=np.divide(rx@ry,den,out=np.zeros(len(rx)),where=den>0);corr=np.clip(corr,-.999999,.999999)
 p=2*t.sf(abs(corr)*np.sqrt((len(sites)-2)/(1-corr*corr)),len(sites)-2);q=bh(p)
 table=pd.DataFrame({'family_id':ids,'spearman_r':corr,'p':p,'q':q,'sites_present':(X>0).sum(axis=1)});table.to_csv(out/f'{axis}.tsv.gz',sep='\t',index=False,compression='gzip')
 summary.append(dict(axis=axis,samples=int(valid.sum()),sites=len(sites),families=len(ids),q_lt_005=int((q<.05).sum()),negative=int(((q<.05)&(corr<0)).sum()),positive=int(((q<.05)&(corr>0)).sum()),))
 pd.DataFrame({'site':sites,'median_chemistry':Y,'samples':[int((ev['site']==s).sum()) for s in sites]}).to_csv(out/f'{axis}_site_inputs.tsv',sep='\t',index=False)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(out/'inputs.sha256.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [E,A]},indent=2)+'\n')
print(json.dumps(summary,indent=2))
