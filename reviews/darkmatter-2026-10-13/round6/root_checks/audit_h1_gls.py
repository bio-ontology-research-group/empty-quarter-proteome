"""Independent, small-table audit of H1 v2 using eigen-whitened least squares."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import linalg, stats
R=Path(__file__).resolve().parents[1]
s=R/'snapshots/h1_v2/results'
d=pd.read_csv(s/'denominator_caller_audit.tsv',sep='\t').sort_values('genome')
c=pd.read_csv(s/'phylo_correlation.tsv.gz',sep='\t',index_col=0)
d['log_len']=np.log10(d['length']);d['log_n50']=np.log10(d['N50']);d['frac_rbh_desert']=d['frac_rbh_of_desert_aai']
cols=['completeness','contamination','log_len','log_n50','aai','frac_rbh_desert']
d=d.dropna(subset=cols+['dark_fraction'])
C=c.loc[d.genome,d.genome].to_numpy();X=np.column_stack([np.ones(len(d)),d[cols].to_numpy()]);y=d.dark_fraction.to_numpy()
prof=pd.read_csv(s/'pgls_lambda_profile.tsv',sep='\t');print(prof.columns.tolist())
lam=0.9158 # printed profile MLE; rounded for diagnostic, not a replacement fit
V=lam*C+(1-lam)*np.eye(len(d))+1e-7*np.eye(len(d))
ev,U=linalg.eigh(V);W=(U/np.sqrt(ev)).T
xw,yw=W@X,W@y;b,_,rank,_=linalg.lstsq(xw,yw)
e=yw-xw@b;rss=float(e@e);z=W@np.ones(len(d));mu=float(z@yw/(z@z));en=yw-mu*z;tss=float(en@en)
se=np.sqrt(np.diag(linalg.inv(xw.T@xw))*rss/(len(d)-rank))
o={'snapshot':'h1_v2','n':len(d),'rank':rank,'lambda_rounded_for_audit':lam,'minimum_covariance_eigenvalue':float(ev.min()),'rss':rss,'correct_gls_intercept_tss':tss,'correct_fixed_covariance_r2':1-rss/tss,'dossier_r2':0.99916,'sigma2_ml':rss/len(d),'sigma2_unbiased_fixed_covariance':rss/(len(d)-rank),'variables':['Intercept']+cols,'beta':b.tolist(),'fixed_covariance_unbiased_se':se.tolist(),'audit_only':True}
(R/'root_checks/h1_v2_gls_audit.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))
